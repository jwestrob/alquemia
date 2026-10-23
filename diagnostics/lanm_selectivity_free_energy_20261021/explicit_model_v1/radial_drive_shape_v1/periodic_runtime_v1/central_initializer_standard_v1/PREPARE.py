"""Generate the finite source/metadata recovery; no native or scheduler calls."""
import copy
import difflib
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
ROOT_SHA='aca6d8725b3b7d3a03865c5f13a86541f09a6f4ef194460aaacc43e3495ad5af'
QOS_OLD="'Account':'banfield','QOS':'standard_expand','Partition':'standard'"
QOS_NEW="'Account':'banfield','QOS':'standard','Partition':'standard'"
REPLACEMENTS=[
    ("build['release']==submission['root_release']", "build['release']==submission['release']"),
    ("submission['job_id']=='20739','Declared actual build20739'", "submission['job_id']=='20813','Declared actual build20813'"),
    ("    m.update(build=request['build'],build_admission=request['build_admission'],probe_source=probe['probe_source'],probe_data=probe['probe_data'])\n",
     "    m.update(build=request['build'],build_admission=request['build_admission'],probe_source=probe['probe_source'],probe_data=probe['probe_data'])\n"
     "    require(set(inputs['standard_helpers'])=={'slurm1','slurm12'},'Both scheduling helpers fixed before initialization')\n"
     "    for name,p in inputs['standard_helpers'].items():verify(p);m['helpers'][name]=p\n"
     "    m['scheduling_recovery']={'declaration':inputs['scheduling_recovery_declaration'],'source_review':request['binder_source_review'],\n"
     "        'original_helpers':{k:prepared['helpers'][k] for k in ('slurm1','slurm12')},'standard_helpers':inputs['standard_helpers'],\n"
     "        'build_submission':inputs['build_submission']}\n"),
    ("    require(x['implementation']==pin(__file__),'Prepared binder pin')\n",
     "    require(x['implementation']==pin(__file__),'Prepared binder pin')\n"
     "    scheduling=read(x['scheduling_recovery_declaration'])\n"
     "    require(scheduling['status']=='declared_scheduling_and_build_lineage_recovery' and scheduling['parent']==x['declaration']\n"
     "        and scheduling['probe_addendum']==x['probe_addendum'],'Exact scheduling amendment of original scientific declaration')\n"
     "    require(request['build']==x['admitted_build'] and request['build_admission']==x['admitted_build_review'],'Exact admitted replacement20813')\n"
     "    require(set(x['standard_helpers'])=={'slurm1','slurm12'},'Both standard helper pins')\n"
     "    for name,p in x['standard_helpers'].items():require(p==x['sources'][name+'_standard'],'Reviewed helper source binding');verify(p)\n"),
    ("    require(review['status']=='passed_central_initializer_binding_source_review' and review['implementation']==pin(__file__) and review['inputs']==request['inputs'],'Independent binder review')\n",
     "    require(review['status']=='passed_central_initializer_binding_source_review' and review['implementation']==pin(__file__) and review['inputs']==request['inputs'],'Independent binder review')\n"
     "    require(review['scheduling_recovery_declaration']==x['scheduling_recovery_declaration'] and review['standard_helpers']==x['standard_helpers'],'Reviewed scheduling delta')\n")]


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,'Changed pin';return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def write(path,obj):
    with path.open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
def text(path,value):
    with path.open('x') as f:f.write(value)


def main(out):
    tick,cpu=time.monotonic(),time.process_time();out.mkdir(parents=True,exist_ok=False)
    a=json.loads((HERE/'AUTHORITIES.json').read_text());d=read(a['root_declaration'])
    assert a['root_declaration']['sha256']==ROOT_SHA
    oldready=read(a['parent_ready']);old=read(oldready['inputs']);m=read(old['central_manifest'])
    assert old['declaration']==d['parent'] and old['probe_addendum']==d['probe_addendum']
    oldctrl=read(oldready['controls']);assert oldctrl['check_count']==20 and oldctrl['error'] is None
    original=verify(a['parent_binder']).read_text();changed=original
    for before,after in REPLACEMENTS:
        assert changed.count(before)==1,before;changed=changed.replace(before,after,1)
    text(HERE/'BIND.py',changed)
    text(HERE/'BIND.diff',''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),fromfile=a['parent_binder']['path'],tofile=str(HERE/'BIND.py'))))
    one=verify(m['helpers']['slurm1']).read_text()
    assert one.count(QOS_OLD)==1 and one.count("qos='standard_expand'")==1
    one_new=one.replace(QOS_OLD,QOS_NEW,1).replace("qos='standard_expand'","qos='standard'",1)
    text(HERE/'SLURM1_STANDARD.py',one_new)
    text(HERE/'SLURM1_STANDARD.diff',''.join(difflib.unified_diff(one.splitlines(True),one_new.splitlines(True),fromfile=m['helpers']['slurm1']['path'],tofile=str(HERE/'SLURM1_STANDARD.py'))))
    twelve=verify(m['helpers']['slurm12']).read_text();twelve_new=verify(a['standard12_helper']).read_text()
    assert twelve.count(QOS_OLD)==1 and twelve_new==twelve.replace(QOS_OLD,QOS_NEW,1)
    standard_ready=read(a['standard12_ready']);standard_review=read(a['standard12_review'])
    assert standard_ready['allocation_parser']==a['standard12_helper'] and standard_review['status']=='passed_periodic_v2_build_source_review'
    assert d['build_source_review']==a['standard12_review']
    for target,authority in [('RUN.sbatch',oldready['wrapper']),('CERTIFICATE.sbatch',a['old_certificate_wrapper'])]:
        wrapper=verify(authority).read_text();assert wrapper.count('#SBATCH --qos=standard_expand\n')==1
        text(HERE/target,wrapper.replace('#SBATCH --qos=standard_expand\n','#SBATCH --qos=standard\n',1))

    build=read(a['build']);admission=read(a['build_admission']);submission=read(a['build_submission'])
    assert build['status']=='periodic_v2_build_complete' and build['error'] is None
    assert admission['status']=='admitted_periodic_v2_build' and admission['result']==a['build']
    assert str(build['allocation']['job_id'])==submission['job_id']=='20813'
    assert build['configuration']==submission['configuration'] and build['release']==submission['release']
    assert build['native_sources']==m['native_sources'] and build['native_programs']==0
    assert set(build['executables'])=={'initialize','parity','gradient_clone','gradient_old'}
    x=copy.deepcopy(old);x['implementation']=pin(HERE/'BIND.py');x['build_submission']=a['build_submission']
    x['scheduling_recovery_declaration']=a['root_declaration']
    x['admitted_build']=a['build'];x['admitted_build_review']=a['build_admission']
    x['standard_helpers']={'slurm1':pin(HERE/'SLURM1_STANDARD.py'),'slurm12':a['standard12_helper']}
    x['sources']={n:pin(HERE/n) for n in ['BIND.py','PREPARE.py','TEST_POLICY.py','RUN.sbatch','CERTIFICATE.sbatch']}
    x['sources'].update({k+'_standard':v for k,v in x['standard_helpers'].items()})
    x['scheduling_recovery_evidence']={'old_ready':a['parent_ready'],'old_controls_reused':oldready['controls'],
        'standard12_ready':a['standard12_ready'],'standard12_review':a['standard12_review'],'policy_snapshot':a['policy_snapshot']}
    for p in x['sources'].values():verify(p)
    write(out/'INPUTS.json',x)
    req=read(oldready['request_template']);req.update(implementation=x['implementation'],inputs=pin(out/'INPUTS.json'),build=a['build'],build_admission=a['build_admission'])
    assert req['status']=='unbound_single_central_initializer_binding' and req['binder_source_review'] is None and req['native_output'] is None
    write(HERE/'REQUEST_TEMPLATE.json',req)
    result={'status':'prepared_standard_central_initializer_source_metadata','implementation':pin(__file__),'declaration':a['root_declaration'],
        'inputs':pin(out/'INPUTS.json'),'binder':pin(HERE/'BIND.py'),'standard_helpers':x['standard_helpers'],
        'actual_build':a['build'],'actual_build_admission':a['build_admission'],'old_controls_reused':oldready['controls'],
        'model_and_rows_unchanged':True,'shared_manifest_not_yet_bound':True,
        'counts':{'native_calls':0,'compiler_calls':0,'solver_calls':0,'numerical_arrays':0,'submissions':0,'scheduler_queries_or_mutations':0},
        'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    write(out/'RESULT.json',result);print(json.dumps({'result':pin(out/'RESULT.json'),'inputs':result['inputs']}))


if __name__=='__main__':
    import sys
    main(Path(sys.argv[1]).resolve())
