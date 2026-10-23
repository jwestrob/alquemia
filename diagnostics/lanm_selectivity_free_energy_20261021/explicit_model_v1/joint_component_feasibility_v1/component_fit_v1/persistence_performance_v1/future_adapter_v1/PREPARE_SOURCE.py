"""Contained future copy: persistence changes only, original files untouched."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ORIGINAL=HERE.parents[1]

def pin(path):return {'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
def replace(text,old,new,count=1):
    if text.count(old)!=count:raise ValueError('Unexpected source shape: '+old[:90])
    return text.replace(old,new)

def main():
    source=(ORIGINAL/'ADAPTER.py').read_text()
    source=replace(source,'import threading\n','import threading\n\n_js=importlib.util.spec_from_file_location("component_ledger_journal",Path(__file__).with_name("JOURNAL.py"))\nJ=importlib.util.module_from_spec(_js);_js.loader.exec_module(J)\n')
    source=replace(source,"    p=Path(p);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\\n');os.replace(tmp,p)","    J.atomic_json(p,x)")
    source=replace(source,"    require(cfg['adapter']==pin(__file__) and cfg['limits']==LIMITS,'New policy source/limits')", "    require(cfg['adapter']==pin(__file__) and cfg['limits']==LIMITS,'New policy source/limits')\n    require(cfg['journal_implementation']==pin(Path(__file__).with_name('JOURNAL.py')) and cfg['persistence_policy']==J.POLICY,'Versioned journal source/policy')")
    source=replace(source,"            save(self.ledger_path,self.ledger)","            self.persistence=J.Journal(self.cache_dir/'EVENTS.jsonl',self.ledger,self.lock,{'configuration_fingerprint':self.batch_fingerprint,'implementation':cfg['adapter'],'journal_implementation':cfg['journal_implementation'],'root_release':cfg.get('root_release')})\n            self.persistence.materialize(self.ledger_path,expected=self.ledger)")
    source=replace(source,"            with self.lock:self.ledger['stopped']=True;self.ledger.setdefault('failure',repr(exc));save(self.ledger_path,self.ledger)","            with self.lock:\n                self.ledger['stopped']=True;self.ledger.setdefault('failure',repr(exc))\n                if hasattr(self,'persistence'):\n                    self.persistence.append('stop',{'error':repr(exc)})\n                    self.persistence.materialize(self.ledger_path)\n                else:save(self.ledger_path,self.ledger)")
    source=replace(source,"                    self.ledger['native_attempts']=copy.deepcopy(self.program_attempts);save(self.ledger_path,self.ledger)","                    self.ledger['native_attempts']=self.program_attempts\n                    self.persistence.append('attempt_intent',{'attempt':attempt,'request_index':len(self.ledger['requests'])-1,'program_inputs':record['program_inputs'][program]})")
    source=replace(source,"                        attempt['returncode_observed']=any('returncode' in r for r in record['receipts'][before:])\n                        self.ledger['native_returncodes_observed']=sum(x.get('returncode_observed',False) for x in self.program_attempts)\n                        self.ledger['native_attempts']=copy.deepcopy(self.program_attempts);save(self.ledger_path,self.ledger)","                        old_observed=attempt.get('returncode_observed',False)\n                        attempt['returncode_observed']=any('returncode' in r for r in record['receipts'][before:])\n                        self.ledger['native_returncodes_observed']+=int(attempt['returncode_observed'])-int(old_observed)\n                        self.ledger['native_attempts']=self.program_attempts\n                        self.persistence.append('attempt_return',{'attempt':attempt,'receipt_delta':record['receipts'][before:],'program_inputs':record['program_inputs'][program]})")
    source=replace(source,"                    self.ledger['requests'].append(request);save(self.ledger_path,self.ledger)","                    self.ledger['requests'].append(request)\n                    request_index=len(self.ledger['requests'])-1\n                    self.persistence.append('request_start',{'index':request_index,'request':request})")
    source=replace(source,"                            with self.lock:request['cells'][cid]=state;save(self.ledger_path,self.ledger)","                            with self.lock:\n                                request['cells'][cid]=state\n                                cell_path=self.cache_dir/'cells'/token/'RESULT.json'\n                                evidence={'result':pin(cell_path)} if cell_path.exists() else {'cache_key':token,'parent_configuration':self.batch['parent_configuration'],'reuse_inventory':self.batch['reuse_inventory'],'immutable_import':token in self.imported_cell_keys}\n                                self.persistence.append('request_cell',{'index':request_index,'cell_id':cid,'state':state,'evidence':evidence})")
    source=replace(source,"                    with self.lock:request.update(status='complete',objective=answer['score']['objective']);save(self.ledger_path,self.ledger)","                    with self.lock:\n                        request.update(status='complete',objective=answer['score']['objective'])\n                        self.persistence.append('request_update',{'index':request_index,'updates':{'status':'complete','objective':answer['score']['objective']}})\n                        self.persistence.materialize(self.ledger_path,expected=self.ledger)")
    source=replace(source,"                    with self.lock:request.update(status='failed',error=repr(exc));save(self.ledger_path,self.ledger)","                    with self.lock:\n                        request.update(status='failed',error=repr(exc))\n                        self.persistence.append('request_update',{'index':request_index,'updates':{'status':'failed','error':repr(exc)}})\n                        self.persistence.materialize(self.ledger_path,expected=self.ledger)")
    source=replace(source,"        def snapshot(self):\n            s=super().snapshot()","        def snapshot(self):\n            with self.lock:self.persistence.materialize(self.ledger_path,expected=self.ledger)\n            s=super().snapshot()")
    for name,text in [('ADAPTER.py',source),('DRIVER.py',(ORIGINAL/'DRIVER.py').read_text())]:
        with (HERE/name).open('x') as f:f.write(text)
    manifest={'status':'prepared_future_persistence_source_only','original_adapter':pin(ORIGINAL/'ADAPTER.py'),'original_driver':pin(ORIGINAL/'DRIVER.py'),
              'implementation':pin(Path(__file__)),'future_adapter':pin(HERE/'ADAPTER.py'),'future_driver':pin(HERE/'DRIVER.py'),'journal':pin(HERE/'JOURNAL.py'),
              'native_calls':0,'scientific_evaluations':0,'fit_execution_released':False,'active_campaign_mutations':0}
    with (HERE/'MANIFEST.json').open('x') as f:json.dump(manifest,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(manifest))
if __name__=='__main__':main()
