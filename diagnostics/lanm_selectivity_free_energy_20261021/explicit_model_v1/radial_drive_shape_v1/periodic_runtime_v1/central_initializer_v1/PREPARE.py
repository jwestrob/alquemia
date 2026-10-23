"""Source and existing-pin inventory only; does not read a live build result."""
import json
from pathlib import Path
import BIND as B

HERE=Path(__file__).resolve().parent
def run(out):
    out.mkdir(parents=True,exist_ok=False)
    decl=B.read(B.pin(HERE/'ROOT_DECLARATION.json'));p=decl['parents']
    ready=B.read(p['runtime_ready']);m=B.read(p['central_manifest']);env=B.read(m['environment'])
    interpreter=str(Path(env['prefix'])/'bin/python')
    x={'status':'prepared_central_initializer_binding_inputs','implementation':B.pin(HERE/'BIND.py'),
       'declaration':B.pin(HERE/'ROOT_DECLARATION.json'),'probe_addendum':B.pin(HERE/'PROBE_BINDING_ADDENDUM.json'),
       'central_ready':p['central_ready'],'central_input_admission':p['central_input_admission'],'central_manifest':p['central_manifest'],
       'runtime_ready':p['runtime_ready'],'runtime_root_admission':p['runtime_root_source_admission'],'runtime_source_review':p['runtime_independent_review'],
       'build_submission':p['build_submission'],'native_implementation':ready['implementations']['native'],
       'native_config_template':B.pin(HERE.parent/'NATIVE_CONFIG_TEMPLATE.json'),'native_release_template':B.pin(HERE.parent/'NATIVE_RELEASE_TEMPLATE.json'),
       'interpreter_path':interpreter,'interpreter':B.pin(interpreter),
       'sources':{n:B.pin(HERE/n) for n in ['BIND.py','PREPARE.py','TEST_BIND.py','RUN.sbatch']},
       'build':None,'build_admission':None,'runtime_release':None,
       'counts':{'native_calls':0,'builds':0,'numerical_arrays':0,'submissions':0}}
    B.write(out/'INPUTS.json',x)
    req={'status':'unbound_single_central_initializer_binding','implementation':x['implementation'],'inputs':B.pin(out/'INPUTS.json'),
       'binder_source_review':None,'build':None,'build_admission':None,'native_output':None,
       'initializer':None,'certificate':None,'central_parity_admission':None,'central_scout_review':None}
    B.write(HERE/'REQUEST_TEMPLATE.json',req)
    print(json.dumps({'inputs':B.pin(out/'INPUTS.json'),'request':B.pin(HERE/'REQUEST_TEMPLATE.json')}))

if __name__=='__main__':
    import sys
    run(Path(sys.argv[1]).resolve())
