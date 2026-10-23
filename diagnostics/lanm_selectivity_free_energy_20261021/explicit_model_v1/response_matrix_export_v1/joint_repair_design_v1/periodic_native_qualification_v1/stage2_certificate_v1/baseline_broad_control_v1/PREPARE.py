"""Metadata-only preparation of the prospectively declared eight baseline controls."""
import argparse
import ast
import copy
import difflib
from BINDING import *

def prepare(out):
    tick=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False)
    authorities=read(pin(HERE/'SOURCE_INPUTS.json'));parent=read(authorities['parent_manifest'])
    original=verify(parent['execution_sources']['NATIVE.py']).read_text()
    actual=(HERE/'RUN.py').read_text()
    require(actual==original.replace('from COMMON import *','from BINDING import *'),'Only native import changes')
    require(authorities['root_declaration']['sha256']==DECLARATION_SHA,'Frozen prospective declaration')
    m=copy.deepcopy(parent)
    m.update(status='prepared_eight_response_disabled_baseline_inputs',producer=pin(__file__),
             parent_manifest=authorities['parent_manifest'],build=authorities['build'],
             root_declaration=authorities['root_declaration'],active_full_result=authorities['active_full_result'],
             central_analytic_reuse=authorities['central_analytic_reuse'],central_scout_review=authorities['central_scout_review'],
             parent_execution_sources=parent['execution_sources'],
             execution_sources={f:pin(HERE/f) for f in ('RUN.py','BINDING.py','PREPARE.py')},
             maximum_native_programs=8,rows={},counts={'native_programs_maximum':8,'new_initializers':0,'new_factors':0,
             'central_gradient_reused':1,'directions':2,'half_steps':2,'signs':2},
             resources={'cpus':1,'AS_bytes':3*1024**3,'wall_seconds':3600,'CPU_seconds':3600,
             'Slurm':'original standard_expand one node/task, default memory/time, no-requeue'},
             source_preparation_actual=None)
    for sid in IDS:
        old=parent['rows'][sid];folder=out/'inputs'/sid;folder.mkdir(parents=True)
        (folder/'system.prepared.key').write_bytes(disabled_key(verify(old['source_inputs']['key']).read_bytes()))
        (folder/'system.key').write_bytes(disabled_key(verify(old['inputs']['key']).read_bytes()))
        source=dict(old['source_inputs']);source['key']=pin(folder/'system.prepared.key')
        inputs=dict(old['inputs']);inputs['key']=pin(folder/'system.key')
        m['rows'][sid]={'role':'BYPASS_OLD_ABSENT','source_inputs':source,'inputs':inputs,
                       'parent_active_inputs':old['inputs'],
                       'output':str(Path(authorities['task_output_parent'])/sid/'execution')}
    validate_manifest(m)
    require(len({r['output'] for r in m['rows'].values()})==8,'Eight disjoint immutable outputs')
    save(out/'MANIFEST.json',m);mp=pin(out/'MANIFEST.json');configs={}
    for sid in IDS:
        folder=out/'bindings'/sid;folder.mkdir(parents=True)
        cfg={'status':'bound_stage2_baseline_control_configuration','manifest':mp,'implementation':m['execution_sources']['RUN.py'],
             'state_id':sid,'operation':'gradient','build':m['build'],'central_scout_review':m['central_scout_review'],
             'certificate':None,'certificate_admission':None}
        save(folder/'CONFIG.json',cfg);configs[sid]=pin(folder/'CONFIG.json')
        save(folder/'ROOT_RELEASE_UNBOUND.json',{'status':None,'implementation':cfg['implementation'],
             'configuration':configs[sid],'output':m['rows'][sid]['output'],'source_review':None,
             'operation':'gradient','state_id':sid,'maximum_native_programs':1})
    (out/'RUN.diff').write_text(''.join(difflib.unified_diff(original.splitlines(keepends=True),actual.splitlines(keepends=True),
                                                       fromfile=str(PARENT/'NATIVE.py'),tofile=str(HERE/'RUN.py'))))
    result={'status':'prepared_eight_response_disabled_baseline_controls','implementation':pin(__file__),
            'source_inputs':pin(HERE/'SOURCE_INPUTS.json'),'manifest':mp,'configurations':configs,
            'original_run_functions_AST_identical':all(ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False)
               for a,b in zip([x for x in ast.parse(original).body if isinstance(x,ast.FunctionDef)],
                              [x for x in ast.parse(actual).body if isinstance(x,ast.FunctionDef)])),
            'original_runtime_helpers':parent['execution_sources']['COMMON.py'],'original_parser':parent['collector'],
            'original_wrapper':parent['execution_sources']['RUN.sbatch'],'build':m['build'],
            'slot_count':8,'new_native_calls':0,'new_initializations':0,'new_factors':0,'submissions':0,
            'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    require(result['original_run_functions_AST_identical'],'Native run AST identity')
    save(out/'RESULT.json',result);print(json.dumps({'result':pin(out/'RESULT.json'),'manifest':mp,'slots':8,'native_calls':0}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);prepare(parser.parse_args().out)
