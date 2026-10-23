"""Literal Stage2 runtime keys and legacy central certificate; zero numerical calls."""
import argparse
import ast
from COMMON import *

HERE=Path(__file__).resolve().parent

def run(out):
    tick=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False)
    pn=HERE.parent.parent;ps=pn.parent/'periodic_stability_definition_v1'
    inputs=json.loads((HERE.parent/'SOURCE_INPUTS_v2.json').read_text())
    original=read(inputs['displacement_manifest'])
    stage1cfg=json.loads((pn/'slurm_stage1_v1/CONFIG.json').read_text())
    wstage2=Path(inputs['source_output']).parent
    wpn=wstage2.parent
    stage1root=pin(wpn/'slurm_stage1_v1/ROOT_ACTUAL_REVIEW.json');actual=read(stage1root)
    require(stage1root['sha256']=='3bb90d1b6b83e28bdb692ed624858d474dff2135c758ff923ecdeb5756b5bba9','frozen Stage1 actual')
    require(actual['status']=='admitted_finite_periodic_operator_parity' and actual['candidate']==inputs['candidate'],'old La parity admitted')
    bypass_pin=pin(wstage2/'bypass_preparation_v1/MANIFEST.json');bypass=read(bypass_pin)
    native_source_review=pin(HERE.parent/'ROOT_NATIVE_SOURCE_REVIEW.json');nr=read(native_source_review)
    require(nr['source_manifest']==pin(Path(inputs['source_output'])/'MANIFEST.json'),'native source review')
    ascfg=json.loads((ps/'candidate_sparse_assembly_v1/CONFIG.json').read_text())
    facdir=ps/'sparse_factorization_design_v1/static_factor_worker_nominal_v1/slurm_execution_v1'
    fcfg=json.loads((facdir/'CONFIG.json').read_text())
    helpers={'slurm1':pin(pn/'slurm_stage1_v1/RUN.py'),'slurm12':pin(facdir/'RUN.py'),
        'snapshot':pin(ps/'pme_export_v1/SNAPSHOT.py'),'pme_checks':pin(ps/'pme_export_v1/full_box_initialization_v1/EXECUTE_v2.py'),
        'assembly':pin(ps/'sparse_assembly_v1/ASSEMBLE.py'),'response_environment':pin(ps/'candidate_sparse_assembly_v1/RUN.py'),
        'backend':fcfg['backend']}
    rows={'central':{'role':'ACTIVE','source_inputs':original['central_prepared_input']}}
    for row in original['rows']:rows[row['id']]={'role':'ACTIVE','source_inputs':row['input']}
    for row in bypass['rows']:
        rows[row['role']]={'role':row['role'],'source_inputs':{k:row['inputs'][v] for k,v in [('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm')]}}
    require(len(rows)==36 and len(original['rows'])==32,'fixed36/32 denominator')
    for state_id,row in rows.items():
        target=out/'inputs'/state_id;target.mkdir(parents=True)
        source=row['source_inputs'];key=verify(source['key']).read_bytes()
        require(b'PME-NORM-EXPORT' not in key.upper(),'observer flag not already present')
        (target/'system.key').write_bytes(key+b'\nPME-NORM-EXPORT 1\n')
        verify(source['xyz']);verify(source['parameters'])
        row['inputs']={'xyz':source['xyz'],'key':pin(target/'system.key'),'parameters':source['parameters']}
    collector=pin(HERE.parent/'collector_v1/COLLECT.py')
    require(collector['sha256']=='66c6259d428511b6b399d30aae4dfb9c65ef470135a1324812dc5c0dbc5218a8','frozen collector parser')
    execution={name:pin(HERE/name) for name in ('COMMON.py','BUILD.py','NATIVE.py','BIND.py','RUN.sbatch')}
    for name in ('COMMON.py','BUILD.py','NATIVE.py','BIND.py'):ast.parse((HERE/name).read_text())
    manifest={'status':'prepared_stage2_execution_inputs_no_numerical_calls','producer':pin(__file__),
        'candidate':inputs['candidate'],'composition':inputs['composition'],'native_sources':nr['source_manifest'],
        'native_source_review':native_source_review,'stage1_actual_admission':stage1root,'stage1_link_recipe':inputs['link_recipe'],
        'central_static_input':inputs['central_static_input'],'numerical_contract':inputs['numerical_contract'],
        'displaced_manifest':inputs['displacement_manifest'],'directions':inputs['directions'],'bypass_manifest':bypass_pin,
        'environment':stage1cfg['environment'],'native_runtime':stage1cfg['native_runtime'],'loader_command':stage1cfg['loader_command'],
        'physics_metadata':ascfg['metadata'],'nominal_contract':fcfg['nominal_contract'],'helpers':helpers,'collector':collector,
        'execution_sources':execution,'rows':rows,'counts':{'active_gradients':33,'bypass_gradients':3,'changed_initializers':32,'changed_assemblies':32,'changed_factors':32,'native_programs':68},
        'resources':{'build_or_native_CPUs':1,'default_memory_MiB':4000,'native_AS_bytes':3*1024**3,'native_wall_seconds':3600,'factor_CPUs':12,'factor_default_memory_MiB':48000,'factor_AS_bytes':32*1024**3,'factor_wall_seconds':900,'factor_CPU_seconds':10800,'maximum_own_factor_jobs':16,'current_queue_review_required':True},
        'source_preparation_actual':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'native_calls':0,'assemblies':0,'factors':0,'solves':0}}
    save(out/'MANIFEST.json',manifest);mp=pin(out/'MANIFEST.json')
    oldstatic=read(actual['static_root_admission']);operator=read(actual['static_input']);factor=read(oldstatic['worker_result'])
    require(oldstatic['status']=='admitted_nominal_saved_La_static_numerical_test' and oldstatic['sufficient_test_passed'] is True,'actual central factor admission')
    require(factor['status']=='passed_nominal_static_sufficient_numerical_test' and factor['sufficient_numerical_test_passed'] is True and factor['input_manifest']==actual['static_input'],'actual central factor input/result')
    center={'status':'stage2_state_nominal_certificate_passed','implementation':pin(__file__),'state_id':'central','candidate':inputs['candidate'],
        'manifest':mp,'source_manifest':nr['source_manifest'],'source_inputs':rows['central']['source_inputs'],'inputs':rows['central']['inputs'],
        'initializer':None,'assembly':operator['operator_preparation_result'],'factor':oldstatic['worker_result'],
        'root_admission':stage1root,'factor_root_admission':actual['static_root_admission'],'expected_state':actual['expected_central_state'],
        'pme_export':actual['initialized_PME'],'nominal':factor['operator'],'sufficient_numerical_test_passed':True,
        'e_app':None,'application_error_bound_established':False,'reused_central_authorities_only':True,
        'cost':{'new_native_calls':0,'new_assemblies':0,'new_factors':0,'new_solves':0},'commands':{}}
    verify(center['expected_state']);verify(center['pme_export']);save(out/'CENTRAL_CERTIFICATE.json',center)
    result={'status':'stage2_execution_preparation_complete','manifest':mp,'central_certificate':pin(out/'CENTRAL_CERTIFICATE.json'),
        'checks':['36 exact input rows','observer flag only appended to literal keys','all32 frozen geometries reused','native source root review joined','actual Stage1 and20053 central authorities reused','all execution Python AST parsed'],
        'actual':manifest['source_preparation_actual']}
    save(out/'RESULT.json',result);print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.out)
