"""Manufactured lineage and nine-allocation controls; no native or actual collection."""
import copy
import json
from pathlib import Path
import time
import COLLECT as C

def main(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);t,c=time.monotonic(),time.process_time();checks=[]
    def ok(name,value):C.require(value,name);checks.append({'name':name,'passed':True})
    def rejects(name,fn):
        try:fn()
        except (ValueError,KeyError):ok(name,True)
        else:ok(name,False)
    cfg={'manifest':'manifest12','build':'newbuild'};m={'execution_sources':{'RUN.py':'run'},'native_sources':'source12'};sub={'configuration':'newcfg','release':'newrelease','job_id':'101'};row={'inputs':'inputs12','source_inputs':'sourceinputs12'};build={'executables':{'gradient_old':{'path':'never_invoked'}}}
    rec={'status':'stage2_gradient_complete','error':None,'state_id':'central','role':'BYPASS_OLD_ABSENT','implementation':'run','manifest':'manifest12','build':'newbuild','source_manifest':'source12','configuration':'newcfg','release':'newrelease','inputs':'inputs12','source_inputs':'sourceinputs12','certificate':None,'native_attempts':1,'counts':dict(C.COUNTS),'allocation':{'job_id':'101','allocated_CPUs':1,'tasks':1,'nodes':1},'commands':{'native':{'status':'complete','returncode':0,'argv':['never_invoked','system.xyz'],'wall_seconds_limit':3600,'CPU_seconds_limit':3600,'AS_bytes':3*1024**3}}}
    call=lambda x:C.receipt_identity(x,'central',row,cfg,m,sub,build)
    call(rec);ok('new_cutoff12_center_accepted',True)
    for key,value in [('state_id','BYPASS_OLD_ABSENT'),('build','oldbuild'),('inputs','inputs9'),('source_manifest','source9'),('role','ACTIVE'),('certificate',{}),('native_attempts',2),('status','stage2_native_failed')]:
        bad=copy.deepcopy(rec);bad[key]=value;rejects('reject_'+key,lambda bad=bad:call(bad))
    bad=copy.deepcopy(rec);bad['allocation']['job_id']='102';rejects('wrong_job_receipt',lambda:call(bad))
    bad=copy.deepcopy(rec);bad['counts']['energy_only_calls']=1;rejects('extra_native_work',lambda:call(bad))
    bad=copy.deepcopy(rec);bad['commands']['native']['AS_bytes']*=2;rejects('changed_resource_limit',lambda:call(bad))
    subs=[{'state_id':sid,'job_id':str(101+i)} for i,sid in enumerate(C.IDS)]
    lines=[f'{s["job_id"]}|fake|COMPLETED|0:0|2|1|00:00:01' for s in subs]
    lines.extend(f'{s["job_id"]}.batch|batch|COMPLETED|0:0|2|1|00:00:01' for s in subs)
    cost=C.scheduler('\n'.join(lines),subs);ok('nine_parents_only_18_CPU_seconds',cost['allocation_count']==9 and cost['allocated_CPU_seconds']==18 and len(cost['records'])==18)
    rejects('eight_parent_incomplete',lambda:C.scheduler('\n'.join(lines[1:]),subs))
    rejects('duplicate_parent_rejected',lambda:C.scheduler('\n'.join(lines+[lines[0]]),subs))
    rejects('unrelated_parent_rejected',lambda:C.scheduler('\n'.join(lines+['999|fake|COMPLETED|0:0|2|1|00:00:01']),subs))
    rejects('failed_native_allocation',lambda:C.scheduler('\n'.join(lines).replace('101|fake|COMPLETED','101|fake|FAILED'),subs))
    rejects('old_eight_slot_count',lambda:C.scheduler('\n'.join(lines),subs[1:]))
    ok('nine_fixed_states_with_new_center',len(C.IDS)==9 and C.IDS[0]=='central' and len(set(C.IDS))==9)
    r={'status':'passed_cutoff12_collector_manufactured_controls','implementation':C.pin(__file__),'collector':C.pin(C.__file__),'checks':checks,'count':len(checks),'native_calls':0,'actual_collections':0,'scientific_matrix_solves':0,'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-c,'threads':1}}
    C.save(out/'RESULT.json',r);print(json.dumps({'result':C.pin(out/'RESULT.json'),'count':len(checks)}))

if __name__=='__main__':
    import sys
    main(sys.argv[1])
