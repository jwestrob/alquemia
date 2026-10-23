"""Versioned persistence prototype. No native/scientific/optimizer execution."""
import copy
import hashlib
import json
import os
from pathlib import Path
import threading

POLICY='append_only_component_ledger_v1'


def encode(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def sha(value):
    return hashlib.sha256(encode(value)).hexdigest()


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def directory_sync(path):
    fd=os.open(str(path),os.O_RDONLY|os.O_DIRECTORY)
    try:os.fsync(fd)
    finally:os.close(fd)


def atomic_json(path,value):
    path=Path(path);temporary=path.with_name(path.name+'.tmp')
    with temporary.open('w') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(temporary,path);directory_sync(path.parent)


class JournalError(RuntimeError):pass


def require(ok,message):
    if not ok:raise JournalError(message)


class State:
    def __init__(self):
        self.ledger=None;self.attempts={};self.program_keys=set();self.finished_cells=set()

    def apply(self,event,dry=False):
        kind=event['kind'];data=event['data']
        if kind=='initialize':
            require(self.ledger is None,'Duplicate initialization')
            value=data['ledger']
            require(value['native_program_invocations']==0 and value['vector_keys']==[] and value['requests']==[],'Fresh ledger only')
            if not dry:self.ledger=copy.deepcopy(value)
            return
        require(self.ledger is not None,'Missing initialization')
        ledger=self.ledger
        if kind=='request_start':
            request=data['request'];index=data['index'];key=request['vector_key']
            require(not ledger['stopped'],'Request after stop')
            require(index==len(ledger['requests']) and key not in ledger['vector_keys'],'Duplicate/misordered request')
            require(len(ledger['vector_keys'])<ledger['limits']['new_parameter_vectors'],'Vector cap')
            require(request['status']=='started','Request initial status')
            if not dry:ledger['vector_keys'].append(key);ledger['requests'].append(copy.deepcopy(request))
        elif kind=='attempt_intent':
            attempt=data['attempt'];number=attempt['number'];key=(attempt['cache_key'],attempt['program'])
            require(not ledger['stopped'],'Attempt after stop')
            require(number==len(self.attempts)+1 and number<=ledger['limits']['new_native_programs'],'Attempt sequence/cap')
            require(number not in self.attempts and key not in self.program_keys,'Duplicate native attempt')
            require(attempt['status']=='invoke_boundary_entered','Intent status')
            require(data['request_index']==len(ledger['requests'])-1 and data['request_index']>=0,'Attempt request')
            request=ledger['requests'][data['request_index']]
            require(request['status']=='started' and attempt['cell_id'] in request['cells'],'Attempt cell/request identity')
            if not dry:
                item=copy.deepcopy(attempt);self.attempts[number]=item;self.program_keys.add(key)
                ledger.setdefault('native_attempts',[]).append(item);ledger['native_program_invocations']+=1
        elif kind=='attempt_return':
            attempt=data['attempt'];number=attempt['number']
            require(number in self.attempts,'Return without intent')
            old=self.attempts[number]
            require(old['status']=='invoke_boundary_entered','Duplicate terminal attempt')
            require(all(old[k]==attempt[k] for k in ('number','cache_key','cell_id','program','folder')),'Attempt identity changed')
            require(attempt['status'] in ('returned_success','failed_native_attempt'),'Attempt terminal status')
            require(type(attempt['returncode_observed']) is bool,'Return evidence flag')
            receipts=data['receipt_delta']
            require(attempt['returncode_observed']==any('returncode' in x for x in receipts),'Return evidence mismatch')
            if attempt['status']=='returned_success':require(len(receipts)==1 and receipts[0].get('returncode')==0,'Success receipt')
            if not dry:
                old.clear();old.update(copy.deepcopy(attempt))
                ledger['native_returncodes_observed']+=int(attempt['returncode_observed'])
        elif kind=='request_cell':
            index=data['index'];cid=data['cell_id'];state=data['state']
            require(0<=index<len(ledger['requests']),'Unknown cell request')
            require(cid in ledger['requests'][index]['cells'] and (index,cid) not in self.finished_cells,'Unknown/duplicate cell completion')
            require(state['status'] in ('complete','failed_partial','not_attempted','failed_before_native'),'Cell status')
            if not dry:
                ledger['requests'][index]['cells'][cid]=copy.deepcopy(state);self.finished_cells.add((index,cid))
        elif kind=='request_update':
            index=data['index'];updates=data['updates']
            require(0<=index<len(ledger['requests']),'Unknown request update')
            require(set(updates)<= {'status','objective','error'} and updates.get('status') in ('complete','failed'),'Request update fields')
            require(ledger['requests'][index]['status']=='started','Duplicate terminal request')
            if updates['status']=='complete':
                require(all(x['status']=='complete' for x in ledger['requests'][index]['cells'].values()),'Incomplete request cannot complete')
            if not dry:ledger['requests'][index].update(copy.deepcopy(updates))
        elif kind=='stop':
            if not dry:ledger['stopped']=True;ledger.setdefault('failure',data['error'])
        else:raise JournalError('Unknown journal event '+kind)


class Journal:
    def __init__(self,path,ledger,lock,authority):
        self.path=Path(path);self.lock=lock;self.state=State();self.seq=0;self.previous='0'*64;self.failed=False
        self.fd=os.open(str(self.path),os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_APPEND,0o600)
        try:
            directory_sync(self.path.parent)
            self.append('initialize',{'ledger':ledger,'authority':authority,'policy':POLICY})
        except Exception:
            os.close(self.fd);self.fd=None;raise

    def append(self,kind,data):
        with self.lock:
            require(not self.failed and self.fd is not None,'Journal unavailable after persistence failure')
            event={'kind':kind,'data':copy.deepcopy(data)}
            self.state.apply(event,dry=True)
            body={'sequence':self.seq+1,'previous_sha256':self.previous,'event':event}
            digest=sha(body);raw=encode({**body,'sha256':digest})+b'\n'
            try:
                offset=0
                while offset<len(raw):
                    written=os.write(self.fd,raw[offset:])
                    require(written>0,'Zero journal write');offset+=written
                os.fsync(self.fd)
            except Exception:
                self.failed=True;raise
            self.state.apply(event);self.seq+=1;self.previous=digest
            return {'sequence':self.seq,'sha256':digest}

    def materialize(self,path,expected=None):
        with self.lock:
            require(not self.failed,'Cannot certify materialization after write failure')
            if expected is not None:require(encode(expected)==encode(self.state.ledger),'Journal/live ledger mismatch')
            atomic_json(path,self.state.ledger)
            # Sidecar makes a possibly stale compatibility ledger explicit.
            atomic_json(Path(path).with_name(Path(path).name+'.journal.json'),
                        {'policy':POLICY,'journal_path':str(self.path),'sequence':self.seq,'last_event_sha256':self.previous,'ledger':pin(path)})
            return pin(path)

    def close(self):
        with self.lock:
            if self.fd is not None:os.close(self.fd);self.fd=None


def replay(path):
    """Read-only replay; an incomplete last record is retained, never repaired."""
    raw=Path(path).read_bytes();state=State();previous='0'*64;sequence=0
    chunks=raw.splitlines(keepends=True);partial=b''
    for pos,line in enumerate(chunks):
        if not line.endswith(b'\n'):
            require(pos==len(chunks)-1,'Nonterminal partial record');partial=line;break
        try:record=json.loads(line)
        except Exception as exc:raise JournalError('Invalid committed JSON record') from exc
        require(set(record)=={'sequence','previous_sha256','event','sha256'},'Record fields')
        body={k:v for k,v in record.items() if k!='sha256'}
        require(record['sequence']==sequence+1 and record['previous_sha256']==previous and sha(body)==record['sha256'],'Sequence/checksum chain')
        state.apply(record['event']);sequence+=1;previous=record['sha256']
    require(sequence>0,'No committed initialization')
    unresolved=[number for number,a in state.attempts.items() if a['status']=='invoke_boundary_entered']
    return {'ledger':copy.deepcopy(state.ledger),'sequence':sequence,'last_event_sha256':previous,
            'partial_suffix_bytes':len(partial),'partial_suffix_sha256':hashlib.sha256(partial).hexdigest() if partial else None,
            'unresolved_attempt_numbers':unresolved,'automatic_retry_permitted':False,
            'committed_prefix_only':bool(partial),'budget_charged_attempts':state.ledger['native_program_invocations']}
