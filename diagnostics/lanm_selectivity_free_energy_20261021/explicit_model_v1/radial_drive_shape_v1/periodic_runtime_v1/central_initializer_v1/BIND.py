"""Immutable metadata binding for one central initializer; invokes no program."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(path),'sha256':h}
def verify(p):
    require(isinstance(p,dict) and set(p)=={'path','sha256'} and pin(p['path'])==p,'Missing/changed immutable pin')
    return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def write(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')

def validate_build(build,admission,request,submission,manifest):
    require(build['status']=='periodic_v2_build_complete' and build.get('error') is None,'Actual successful build required')
    require(admission['status']=='admitted_periodic_v2_build' and admission['result']==request['build'],'Exact root build admission')
    require(build['configuration']==submission['configuration'] and build['release']==submission['root_release'],'Declared build configuration/release')
    require(str(build['allocation']['job_id'])==submission['job_id']=='20739','Declared actual build20739')
    require(build['native_sources']==manifest['native_sources'],'Same runtime native source')
    require(build['native_programs']==0 and build['attempts']=={'configure':1,'archive':1,'compile':8,'link':3},'Actual build-only finite denominator')
    require(set(build['executables'])=={'initialize','parity','gradient_clone','gradient_old'},'Complete four-role build bundle')
    old=build['old_baseline'];require(old['executable']==build['executables']['gradient_old'],'Matching old baseline executable')
    require(old['ticket_schema']=='STAGE2_BOUND_INPUT_V1' and old['state_schema']=='PERIODIC_PARITY_STATE_V1','Old baseline lineage remains V1')

def bound_manifest(prepared,request,probe,inputs):
    require(prepared['status']=='prepared_periodic_v2_central_inputs_not_executable','Admitted prepared manifest')
    require(prepared['rows']['central']['role']=='ACTIVE','Central ACTIVE row')
    require(all(prepared.get(k) is None for k in ('build','build_admission','runtime_release','initializer','certificate','central_parity_admission','central_scout_admission','probe_source','probe_data')),'No prior state/result chain')
    m=copy.deepcopy(prepared);m['status']='bound_periodic_v2_runtime_inputs'
    m.update(build=request['build'],build_admission=request['build_admission'],probe_source=probe['probe_source'],probe_data=probe['probe_data'])
    m['binding']={'implementation':pin(__file__),'declaration':inputs['declaration'],'prepared_manifest':inputs['central_manifest'],
       'input_admission':inputs['central_input_admission'],'probe_addendum':inputs['probe_addendum']}
    return m

def context(request_path):
    request=json.loads(request_path.read_text());require(request['status']=='requested_single_central_initializer_binding','Unbound request')
    require(request['implementation']==pin(__file__),'Exact binder')
    x=read(request['inputs']);require(x['status']=='prepared_central_initializer_binding_inputs','Prepared source inputs')
    require(x['implementation']==pin(__file__),'Prepared binder pin')
    for p in x['sources'].values():verify(p)
    declaration=read(x['declaration']);probe=read(x['probe_addendum'])
    require(probe['parent_declaration']==x['declaration'] and probe['status']=='declared_fixed_probe_metadata_binding','Root probe addendum')
    review=read(request['binder_source_review'])
    require(review['status']=='passed_central_initializer_binding_source_review' and review['implementation']==pin(__file__) and review['inputs']==request['inputs'],'Independent binder review')
    prepared=read(x['central_manifest']);central=read(x['central_input_admission']);ready=read(x['central_ready'])
    require(central['status']=='admitted_periodic_v2_central_inputs' and central['manifest']==x['central_manifest'] and central['ready']==x['central_ready'],'Root central input admission')
    independent=read(central['independent_review'])
    require(independent['status']=='passed_radial_periodic_central_binding_source_review' and independent['manifest']==x['central_manifest'],'Independent central input admission')
    require(ready['candidate']==central['candidate']==prepared['model']['candidate']==declaration['candidate'],'Exact declared terminal candidate')
    require(prepared['model']['settings']==declaration['settings'],'Fixed cutoff12 settings')
    runtime=read(x['runtime_ready']);runtime_root=read(x['runtime_root_admission']);runtime_review=read(x['runtime_source_review'])
    require(runtime_root['status']=='admitted_periodic_v2_runtime_source' and runtime_root['ready']==x['runtime_ready'] and runtime_root['independent_review']==x['runtime_source_review'],'Root runtime source admission')
    require(runtime_review['status']=='passed_periodic_v2_runtime_source_review','Independent runtime source admission')
    require(runtime_review['execution_sources']==runtime['execution_sources']==prepared['execution_sources'],'Same execution sources')
    require(runtime_review['native_sources']==runtime['native_sources']==prepared['native_sources'],'Same native sources')
    require(runtime['implementations']['native']==x['native_implementation'],'Exact native producer')
    # NATIVE.source_context/native_build retain their existing raw pin/runtime
    # checks. This binder joins their admitted metadata without re-admitting it.
    env=read(prepared['environment']);require(str(Path(env['prefix'])/'bin/python')==x['interpreter_path'],'Admitted interpreter path')
    verify(x['interpreter'])
    previous=read(probe['source_manifest']);history=read(probe['historical_admission']);old_result=read(probe['historical_result'])
    require(history['result']==probe['historical_result'] and old_result['preparation']==probe['source_manifest'],'Admitted fixed-probe provenance')
    require(previous['probe_dimension']==probe['dimension']==35874 and previous['labels']==probe['labels'] and len(probe['labels'])==10,'Fixed ten-probe dimension/order')
    require(previous['probe_source']==probe['probe_source'] and previous['raw_inputs']['probes.dat']==probe['probe_data'],'Exact prior source/data only')
    verify(probe['probe_source']);verify(probe['probe_data'])
    submission=read(x['build_submission']);build=read(request['build']);admission=read(request['build_admission'])
    validate_build(build,admission,request,submission,prepared)
    require(all(request.get(k) is None for k in ('initializer','certificate','central_parity_admission','central_scout_review')),'Fresh initializer without downstream results')
    require(Path(request['native_output']).is_absolute() and not Path(request['native_output']).exists(),'Fresh absolute native output required')
    return request,x,prepared,probe

def run(request_path,out):
    tick=time.monotonic();cpu=time.process_time()
    request,x,prepared,probe=context(request_path)
    out.mkdir(parents=True,exist_ok=False)
    m=bound_manifest(prepared,request,probe,x);write(out/'MANIFEST.json',m)
    cfg=read(x['native_config_template']);cfg.update(status='bound_periodic_v2_native_configuration',manifest=pin(out/'MANIFEST.json'),
        build=request['build'],build_admission=request['build_admission'],state_id='central',operation='initialize')
    require(all(cfg[k] is None for k in ('central_parity_admission','central_scout_review','certificate','certificate_admission','probe_source','probe_data')),'No bootstrap circularity or probe execution')
    write(out/'CONFIG_RUN.json',cfg)
    release=read(x['native_release_template']);release.update(configuration=pin(out/'CONFIG_RUN.json'),output=str(Path(request['native_output']).resolve()),
        source_review=x['runtime_source_review'],state_id='central',operation='initialize',binding_request=pin(request_path),binder_source_review=request['binder_source_review'])
    require(release['status']=='unreleased_periodic_v2_native','Binder cannot release a native call')
    write(out/'ROOT_RELEASE_TEMPLATE.json',release)
    result={'status':'bound_single_central_initializer_not_released','implementation':pin(__file__),'inputs':request['inputs'],'request':pin(request_path),
        'manifest':pin(out/'MANIFEST.json'),'configuration':pin(out/'CONFIG_RUN.json'),'release_template':pin(out/'ROOT_RELEASE_TEMPLATE.json'),
        'native_implementation':x['native_implementation'],'candidate':prepared['model']['candidate'],'build':request['build'],'build_admission':request['build_admission'],
        'prospective_argv':[x['interpreter_path'],x['native_implementation']['path'],'--config',str(out/'CONFIG_RUN.json'),'--release','<ROOT_RELEASE>','--out',release['output']],
        'scope':'One direct initializer only. Full V2 body/h/H admission remains CERTIFY; later phase artifacts belong in configurations, never in this immutable manifest.',
        'counts':{'native_calls':0,'builds':0,'probe_generation':0,'numerical_replay':0,'submissions':0},
        'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    write(out/'RESULT.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--request',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    print(json.dumps(run(a.request.resolve(),a.out.resolve())))
