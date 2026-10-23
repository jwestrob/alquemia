#!/usr/bin/env python3
"""Recover four saved field parses, then only three untouched serial bypasses."""
import argparse,copy,hashlib,json,os,resource,shutil,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import PREPARE as P
import CHECKS as C
import FIELD_MATRIX as F
LIMITS={'new_native_programs':3,'combined_native_programs':49,'saved_complete_rows':42,'saved_field_reparses':4,'reused_internal_FD_energies':132,'new_internal_FD_energies':0,'native_workers':1,'allocated_CPUs':12,'program_timeout_seconds':120}
ERROR="ValueError('Zero-field export required')"

def fingerprint(cfg):return hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k!='root_release'},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def verify_pins(x):
    if isinstance(x,dict):
        if set(x)=={'path','sha256'}:P.verify(x)
        else:
            for v in x.values():verify_pins(v)
    elif isinstance(x,list):
        for v in x:verify_pins(v)

def parent_evidence(old,cfg,m):
    P.require(old['status']=='radial_drive_numerical_qualification_failed' and old['error']==ERROR and old['global_stop'] is True and old['numerical_checks_passed'] is False,'Preserved exact original failure')
    P.require(old['configuration']==cfg and old['implementation']==cfg['implementation'] and old['root_release']==cfg['root_release'] and old['build']==cfg['build'] and old['executables']==cfg['executables'],'Original execution/source configuration')
    P.require(old['native_program_invocation_attempts']==old['native_programs_returned_zero']==46 and old['native_programs_completed']==42,'Original46/42 accounting')
    P.require(all(old[k]==132 for k in ('internal_FD_energies_observed','internal_FD_energies_requested','internal_FD_energies_completed')),'Original132 FD accounting')
    P.require(old['zero_parity_passed'] is True and not old['cache_bridge_admitted'] and not old['physical_accuracy_qualified'],'No original bridge/physical admission')
    P.require(set(old['rows'])=={j['id'] for j in m['programs']} and len(old['operator_checks'])==6 and len(old['directional_checks'])==4 and old['field_checks']=={},'Original full stage denominator')
    for job in m['programs']:
        rec=old['rows'][job['id']];stage=job['stage'];folder=Path(cfg['output'])/'cells'/job['id']
        if stage=='bypass':
            P.require(rec=={'status':'not_attempted'} and not folder.exists(),'Bypass has never been attempted');continue
        P.require(rec['job']==job and rec['status']==('failed' if stage=='fields' else 'complete'),'Exact original row/job status')
        if stage=='fields':P.require(rec['error']==ERROR and len(rec['exports'])==2 and 'matrix_qualification' not in rec,'Only failed field-header parsing')
        P.require(json.loads((folder/'RESULT.json').read_text())==rec,'Literal original per-cell record')
        verify_pins(rec)
        P.require(all(rec['runtime_input'][k]['sha256']==job['input'][k]['sha256'] for k,_ in P.FILES),'Original literal prepared/runtime bytes')
        receipts=rec['receipts'];args=['system.xyz','Y','Y','0.00001'] if job['fd'] else ['system.xyz','Y','N'] if job['program']=='testgrad' else ['system.xyz','E' if stage=='directions' else 'EM']
        P.require(len(receipts)==1 and receipts[0]['returncode']==0 and receipts[0]['argv']==[cfg['executables'][job['program']]['path']]+args and receipts[0]['cwd']==str(folder) and receipts[0]['output']==rec['actual_output'] and receipts[0]['timeout_seconds']==120,'Exact original native invocation evidence')
        P.require(all(v=='1' for v in receipts[0]['threads'].values()),'Original serial native numerical threads')
    P.require(sum(old['rows'][j['id']]['response_readback']['native_solve_blocks']-1 for j in m['programs'] if j['fd'])==132,'Saved observed FD denominator')

def context(cfg,runtime=False):
    P.require(cfg['implementation']==P.pin(__file__) and cfg['field_adapter']==P.pin(F.__file__) and cfg['declaration']==P.pin(HERE/'DECLARATION.md') and cfg['limits']==LIMITS,'Frozen recovery source/scope')
    inp=P.read(cfg['inputs']);verify_pins(inp)
    P.require(inp['parent_runner']==P.pin(HERE.parent/'RUN.py') and inp['parent_checks']==P.pin(C.__file__) and inp['parent_preparer']==P.pin(P.__file__),'Exact frozen parent code')
    parent=P.module(inp['parent_runner'],'radial_recovery_parent');pc=P.read(inp['parent_configuration']);old=P.read(inp['parent_result']);m,h,ev=parent.context(pc,runtime)
    P.require(P.pin(h['matrix'].__file__)==inp['matrix_helper'],'Exact inherited numerical kernel')
    P.require('parameter (elefield=1439.96455d0)' in P.verify(inp['native_units']).read_text() and 'exfld(i) = exfld(i) / elefield' in P.verify(inp['native_prmkey']).read_text(),'Native field conversion source')
    parent_evidence(old,pc,m)
    diagnosis=P.read(inp['diagnosis']);P.require(diagnosis['terminal_result']==inp['parent_result'] and diagnosis['zero_field_only_parser']==inp['matrix_helper'],'Independent failure localization')
    acct=P.read(inp['parent_sacct']);P.require(acct['job_id']==old['job_id'] and acct['state']=='FAILED' and acct['allocated_CPU_seconds']==288 and acct['parent_only'],'Preserved actual scheduler accounting')
    if runtime:
        review=P.read(cfg['source_review']);P.require(review['status']=='passed_radial_field_recovery_source_review' and all(review[k]==cfg[k] for k in ('implementation','field_adapter','inputs')),'Exact recovery source review')
        release=P.read(cfg['root_release']);P.require(release['status']=='released_radial_field_recovery' and release['configuration_fingerprint']==fingerprint(cfg) and release['limits']==LIMITS and release['output']==cfg['output'] and release['parent_result']==inp['parent_result'],'Exact recovery release')
    return inp,parent,pc,old,m,h,ev

def reparse(rec,job,row,m,h,ev):
    rec.pop('error',None);rec['status']='reparsing_saved_output';hh=dict(h);hh['matrix']=F.scoped(h['matrix'],job)
    C.parse(rec,job,row,m,hh,ev);rec['status']='complete';return rec

def run(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
    inp,parent,pc,old,m,h,ev=context(cfg,True);invoke=parent.invocation(h)
    r={'status':'running','implementation':P.pin(__file__),'configuration':cfg,'inputs':cfg['inputs'],'parent_result':inp['parent_result'],'preserved_parent_error':old['error'],'parent_execution_cost':old['cost'],'parent_scheduler_cost':P.read(inp['parent_sacct']),'manifest':pc['manifest'],'source_manifest':pc['source_manifest'],'build':pc['build'],'executables':pc['executables'],'job_id':os.environ['SLURM_JOB_ID'],'host':os.uname().nodename,'pid':os.getpid(),'rows':copy.deepcopy(old['rows']),'operator_checks':copy.deepcopy(old['operator_checks']),'directional_checks':copy.deepcopy(old['directional_checks']),'field_checks':{},'saved_field_qualification_attempts':0,'saved_fields_qualified':0,'new_native_program_invocation_attempts':0,'new_native_programs_returned_zero':0,'new_native_programs_completed':0,'native_program_invocation_attempts':46,'native_programs_returned_zero':46,'native_programs_completed':42,'internal_FD_energies_observed':132,'internal_FD_energies_requested':132,'internal_FD_energies_completed':132,'new_internal_FD_energies':0,'zero_parity_passed':True,'numerical_checks_passed':False,'cache_bridge_admitted':False,'physical_accuracy_qualified':False,'global_stop':False,'error':None}
    def persist():P.save(out/'RESULT.tmp',r);os.replace(out/'RESULT.tmp',out/'RESULT.json')
    def unchanged42():return all(r['rows'][j['id']]==old['rows'][j['id']] for j in m['programs'] if j['stage'] not in ('fields','bypass'))
    try:
        persist()
        for job in [j for j in m['programs'] if j['stage']=='fields']:
            r['saved_field_qualification_attempts']+=1;persist()
            try:reparse(r['rows'][job['id']],job,m['rows'][job['cell']],m,h,ev)
            except Exception as e:
                r['rows'][job['id']].update(status='failed',recovery_error=repr(e));raise
            r['saved_fields_qualified']+=1;r['native_programs_completed']+=1;persist()
        for g in P.DIMERS:
            r['field_checks'][g]=C.field_checks(r['rows']['active__'+g+'__analyze'],r['rows']['fields__'+g+'__analyze__-1'],r['rows']['fields__'+g+'__analyze__1'],m['rows'][g],m,h)
        P.require(r['saved_fields_qualified']==4 and len(r['field_checks'])==2 and unchanged42(),'All saved field gates before any bypass');persist()
        for job in [j for j in m['programs'] if j['stage']=='bypass']:
            P.require(r['new_native_program_invocation_attempts']<3 and r['rows'][job['id']]=={'status':'not_attempted'},'Only untouched bypass allowlist')
            folder=out/'cells'/job['id'];folder.mkdir(parents=True,exist_ok=False);rec={'status':'not_attempted','job':job,'receipts':[]};r['rows'][job['id']]=rec
            try:
                for k,name in P.FILES:shutil.copyfile(P.verify(job['input'][k]),folder/name)
                rec['runtime_input']={k:P.pin(folder/name) for k,name in P.FILES};rec['status']='started'
                r['new_native_program_invocation_attempts']+=1;r['native_program_invocation_attempts']+=1;persist()
                output=folder/'analyze.out';invoke([pc['executables']['analyze']['path'],'system.xyz','EM'],folder,output,rec['receipts']);rec['actual_output']=P.pin(output)
                r['new_native_programs_returned_zero']+=1;r['native_programs_returned_zero']+=1
                C.parse(rec,job,m['rows'][job['cell']],m,h,ev);rec['status']='complete';r['new_native_programs_completed']+=1;r['native_programs_completed']+=1
            except Exception as e:rec.update(status='failed',error=repr(e));raise
            finally:
                if (folder/'analyze.out').exists():rec['actual_output']=P.pin(folder/'analyze.out')
                P.save(folder/'RESULT.json',rec);persist()
        P.require(unchanged42() and r['saved_fields_qualified']==4 and r['new_native_programs_completed']==3 and r['native_program_invocation_attempts']==r['native_programs_returned_zero']==r['native_programs_completed']==49 and all(v['status']=='complete' for v in r['rows'].values()),'Exact recovered49/132 terminal accounting')
        r.update(status='radial_drive_numerical_qualification_recovered',numerical_checks_passed=True)
    except Exception as e:r.update(status='radial_drive_numerical_qualification_recovery_failed',global_stop=True,error=repr(e))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);r['reused42_rows_exact']=unchanged42();r['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':ch.ru_utime-child.ru_utime,'child_system_CPU_seconds':ch.ru_stime-child.ru_stime,'peak_child_RSS_KiB':ch.ru_maxrss,'allocated_CPUs':12,'native_workers':1,'parent_cost_not_added':True,'per_receipt_child_deltas_not_summed':True};persist()
    return r

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=('preflight','run'));ap.add_argument('--config',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();cfg=json.loads(Path(a.config).read_text());out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    if a.mode=='preflight':
        t,c=time.monotonic(),time.process_time()
        try:
            inp,parent,pc,old,m,h,ev=context(cfg,False)
            # Header-only validation of eight saved exports; no vectors or matrices are solved.
            headers=[]
            for job in [j for j in m['programs'] if j['stage']=='fields']:
                rec=old['rows'][job['id']];row=copy.deepcopy(m['rows'][job['cell']]);row['input']=rec['runtime_input']
                for q in rec['exports']:headers.append(F.header(P.verify(q).read_text().splitlines()[4].split(),row,m['units'],job,h['matrix']))
            r={'status':'passed_radial_recovery_zero_call_preflight','implementation':P.pin(__file__),'configuration':P.pin(a.config),'inputs':cfg['inputs'],'reused_complete_rows':42,'saved_field_outputs':4,'saved_export_headers':headers,'untouched_native_slots':3,'new_native_calls':0,'scientific_matrix_qualifications':0,'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-c}};P.save(out/'RESULT.json',r)
        except Exception as e:P.save(out/'FAILURE.json',{'error':repr(e),'new_native_calls':0});raise
    else:P.require(str(out)==cfg['output'],'Fresh released output path');r=run(cfg,out)
    print(json.dumps({'status':r['status']}));raise SystemExit(1 if 'failed' in r['status'] else 0)
