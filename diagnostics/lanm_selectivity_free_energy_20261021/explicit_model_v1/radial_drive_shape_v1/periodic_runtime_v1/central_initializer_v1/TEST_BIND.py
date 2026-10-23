"""Focused manufactured metadata checks; no native or numerical array operations."""
import ast
import copy
import json
from pathlib import Path
import time
import BIND as B

def run(inputs,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    r={'status':'failed','implementation':B.pin(__file__),'inputs':B.pin(inputs),'checks':[],'error':None,
       'counts':{'native_calls':0,'builds':0,'numerical_arrays':0,'probe_generation':0,'scheduler_calls':0}}
    def check(name,ok):B.require(ok,name);r['checks'].append(name)
    def rejects(name,fn):
        try:fn()
        except (ValueError,KeyError):check(name,True)
        else:raise AssertionError(name)
    try:
        x=json.loads(inputs.read_text());m=B.read(x['central_manifest']);probe=B.read(x['probe_addendum']);d=B.read(x['declaration'])
        for p in x['sources'].values():B.verify(p)
        for name in ('BIND.py','PREPARE.py','TEST_BIND.py'):ast.parse(Path(x['sources'][name]['path']).read_text())
        check('exact source pins and AST parse',True)
        for key in ('central_ready','central_input_admission','central_manifest','runtime_ready','runtime_root_admission','runtime_source_review','build_submission','native_implementation','native_config_template','native_release_template'):B.verify(x[key])
        check('declared existing authorities',True)
        root=B.read(x['central_input_admission']);ready=B.read(x['central_ready'])
        check('input root joins exact manifest/candidate',root['manifest']==x['central_manifest'] and root['candidate']==ready['candidate']==m['model']['candidate']==d['candidate'])
        check('both radial slots and fixed settings retained',set(m['model']['radial'])=={'900','901'} and m['model']['settings']==d['settings'])
        previous=B.read(probe['source_manifest']);old=B.read(probe['historical_result']);admitted=B.read(probe['historical_admission'])
        check('fixed probe provenance joins',admitted['result']==probe['historical_result'] and old['preparation']==probe['source_manifest'])
        B.verify(probe['probe_source']);B.verify(probe['probe_data'])
        check('probe dimension and exact source/data byte pins',probe['dimension']==previous['probe_dimension']==35874 and len(probe['labels'])==10 and previous['probe_source']==probe['probe_source'] and previous['raw_inputs']['probes.dat']==probe['probe_data'])
        fake={'path':'/manufactured/not-a-real-build.json','sha256':'0'*64}
        req={'build':fake,'build_admission':{'path':'/manufactured/not-an-admission.json','sha256':'1'*64}}
        before=copy.deepcopy(m);bound=B.bound_manifest(m,req,probe,x)
        check('immutable bound status with unchanged prepared rows/model',m==before and bound['status']=='bound_periodic_v2_runtime_inputs' and bound['rows']==m['rows'] and bound['model']==m['model'])
        check('both fixed probes bound before initialization',bound['probe_source']==probe['probe_source'] and bound['probe_data']==probe['probe_data'])
        check('downstream result slots remain null',all(bound[k] is None for k in ('runtime_release','initializer','certificate','central_parity_admission','central_scout_admission')))
        bad=copy.deepcopy(m);bad['initializer']=fake
        rejects('old initializer reuse rejected',lambda:B.bound_manifest(bad,req,probe,x))
        bad=copy.deepcopy(m);bad['rows']['central']['role']='BYPASS_CLONE_ABSENT'
        rejects('nonactive central role rejected',lambda:B.bound_manifest(bad,req,probe,x))
        sub=B.read(x['build_submission'])
        build={'status':'periodic_v2_build_complete','error':None,'configuration':sub['configuration'],'release':sub['root_release'],
          'allocation':{'job_id':'20739'},'native_sources':m['native_sources'],'native_programs':0,'attempts':{'configure':1,'archive':1,'compile':8,'link':3},
          'executables':dict.fromkeys(('initialize','parity','gradient_clone','gradient_old'),fake),
          'old_baseline':{'executable':fake,'ticket_schema':'STAGE2_BOUND_INPUT_V1','state_schema':'PERIODIC_PARITY_STATE_V1'}}
        adm={'status':'admitted_periodic_v2_build','result':fake}
        B.validate_build(build,adm,req,sub,m);check('manufactured complete declared build metadata accepted',True)
        for label,change in [('started build',lambda b:b.update(status='started')),('failed build',lambda b:b.update(status='periodic_v2_build_failed',error='failure')),('different job',lambda b:b['allocation'].update(job_id='20740')),('missing role',lambda b:b['executables'].pop('parity'))]:
            bad=copy.deepcopy(build);change(bad);rejects(label+' rejected',lambda:B.validate_build(bad,adm,req,sub,m))
        bad=copy.deepcopy(adm);bad['result']={'path':'/manufactured/wrong','sha256':'2'*64}
        rejects('wrong root result join rejected',lambda:B.validate_build(build,bad,req,sub,m))
        cfg=B.read(x['native_config_template']);cfg.update(status='bound_periodic_v2_native_configuration',operation='initialize',state_id='central',build=fake,build_admission=req['build_admission'])
        check('fixed producer/central initializer without circular phase prerequisites',cfg['implementation']==x['native_implementation'] and cfg['operation']=='initialize' and cfg['state_id']=='central' and all(cfg[k] is None for k in ('certificate','certificate_admission','central_parity_admission','central_scout_review','probe_source','probe_data')))
        check('release remains unissued',B.read(x['native_release_template'])['status']=='unreleased_periodic_v2_native')
        wrapper=Path(x['sources']['RUN.sbatch']['path']).read_text()
        check('one CPU default resources and admitted interpreter',all(s in wrapper for s in ('--cpus-per-task=1','--ntasks=1','--qos=standard_expand','--no-requeue',x['interpreter_path'])) and not any(s in wrapper for s in ('#SBATCH --mem','#SBATCH --time')))
        r['status']='passed_central_initializer_binding_source_controls'
    except BaseException as exc:r['error']=repr(exc)
    r['check_count']=len(r['checks']);r['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
    B.write(out/'RESULT.json',r);return r

if __name__=='__main__':
    import sys
    r=run(Path(sys.argv[1]),Path(sys.argv[2]));print(json.dumps({k:r[k] for k in ('status','check_count','error','cost')}));raise SystemExit(0 if r['error'] is None else 1)
