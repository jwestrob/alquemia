"""Prepare only literal keys, a one-constant harness delta and future bindings."""
import argparse
import ast
import difflib
from BINDING import *

def prepare(out):
    tick=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False)
    a=read(pin(HERE/'SOURCE_INPUTS.json'));parent=read(a['parent_manifest']);oldbuild=read(a['parent_build'])
    decl=read(a['root_declaration']);base=read(a['baseline_manifest']);central=read(decl['central_baseline_result'])
    original=verify(parent['execution_sources']['NATIVE.py']).read_text();actual=(HERE/'RUN.py').read_text()
    require(actual==original.replace('from COMMON import *','from BINDING import *'),'Only runtime import substitution')
    native=read(parent['native_sources']);sourceout=Path(a['native_source_output']);sourceout.mkdir(parents=True,exist_ok=False)
    gate=verify(native['files']['GATE.f90']).read_bytes();require(gate.count(b'ewaldcut/=9d0')==1,'One literal cutoff guard')
    changed=gate.replace(b'ewaldcut/=9d0',b'ewaldcut/=12d0');(sourceout/'GATE.f90').write_bytes(changed)
    (sourceout/'GATE.diff').write_text(''.join(difflib.unified_diff(gate.decode().splitlines(keepends=True),changed.decode().splitlines(keepends=True),fromfile=native['files']['GATE.f90']['path'],tofile=str(sourceout/'GATE.f90'))))
    nm={'status':'prepared_cutoff12_harness_only','producer':pin(__file__),'parent_manifest':parent['native_sources'],
        'files':{'GATE.f90':pin(sourceout/'GATE.f90'),**{n:native['files'][n] for n in ('GRADIENT.f90','STATE_WRITER.f90')}},
        'patch':pin(sourceout/'GATE.diff'),'only_change':'ewaldcut/=9d0 to ewaldcut/=12d0','native_kernel_changes':0}
    save(sourceout/'MANIFEST.json',nm)
    inherited=('helpers','collector','environment','native_runtime','loader_command','candidate','composition','stage1_link_recipe','directions','displaced_manifest')
    m={k:parent[k] for k in inherited}
    m.update(status='prepared_nine_cutoff12_baseline_inputs',producer=pin(__file__),source_inputs=pin(HERE/'SOURCE_INPUTS.json'),
        root_declaration=a['root_declaration'],baseline_manifest=a['baseline_manifest'],parent_manifest=a['parent_manifest'],parent_build=a['parent_build'],
        baseline_failed_collection=a['baseline_failed_collection'],native_sources=pin(sourceout/'MANIFEST.json'),
        execution_sources={n:pin(HERE/n) for n in ('BUILD.py','BINDING.py','RUN.py','PREPARE.py','BIND.py')},
        original_wrapper=parent['execution_sources']['RUN.sbatch'],reused_objects={},maximum_native_programs=9,rows={},build_output=a['build_output'],
        counts={'compile_attempts_maximum':1,'link_attempts_maximum':1,'native_gradient_programs_maximum':9,'new_initializers':0,'new_factors':0,'new_central_gradient':1,'displaced_gradients':8},
        resources={'CPUs_per_job':1,'gradient_AS_bytes':3*1024**3,'gradient_wall_seconds':3600,'gradient_CPU_seconds':3600,
                   'compile_link_phase_wall_seconds':180,'compile_link_phase_CPU_seconds':180,'Slurm':'one node/task, banfield/standard/standard_expand, default memory/time, no-requeue'},
        interpretation='Baseline cutoff diagnostic only; candidate metadata retained for lineage, local response absent in all nine keys.')
    for name in ('GRADIENT.f90','STATE_WRITER.f90'):
        argv=oldbuild['commands']['compile_'+name.replace('.','_')]['argv'];m['reused_objects'][name]=pin(argv[argv.index('-o')+1])
    for sid in IDS:
        old=central if sid=='central' else base['rows'][sid];folder=out/'inputs'/sid;folder.mkdir(parents=True)
        (folder/'system.prepared.key').write_bytes(cutoff_key(verify(old['source_inputs']['key']).read_bytes()))
        (folder/'system.key').write_bytes(cutoff_key(verify(old['inputs']['key']).read_bytes()))
        source={**old['source_inputs'],'key':pin(folder/'system.prepared.key')};inputs={**old['inputs'],'key':pin(folder/'system.key')}
        m['rows'][sid]={'role':'BYPASS_OLD_ABSENT','source_inputs':source,'inputs':inputs,'cutoff9_inputs':old['inputs'],
                       'output':str(Path(a['task_output_parent'])/sid/'execution')}
    validate_manifest(m);save(out/'MANIFEST.json',m);mp=pin(out/'MANIFEST.json')
    cfg={'status':'bound_cutoff12_build_configuration','manifest':mp,'implementation':m['execution_sources']['BUILD.py']}
    save(out/'CONFIG_BUILD.json',cfg)
    save(out/'ROOT_RELEASE_BUILD_UNBOUND.json',{'status':None,'implementation':cfg['implementation'],'configuration':pin(out/'CONFIG_BUILD.json'),'output':m['build_output'],'source_review':None})
    save(out/'NATIVE_BINDING_TEMPLATE.json',{'status':None,'manifest':mp,'implementation':m['execution_sources']['RUN.py'],'build':None,'build_root_review':None,'state_ids':list(IDS),'native_release':None})
    (out/'RUN.diff').write_text(''.join(difflib.unified_diff(original.splitlines(keepends=True),actual.splitlines(keepends=True),fromfile=parent['execution_sources']['NATIVE.py']['path'],tofile=str(HERE/'RUN.py'))))
    result={'status':'prepared_nine_cutoff12_baseline_no_native','implementation':pin(__file__),'manifest':mp,'source_inputs':pin(HERE/'SOURCE_INPUTS.json'),
        'native_sources':m['native_sources'],'build_configuration':pin(out/'CONFIG_BUILD.json'),'slot_count':9,
        'original_run_functions_AST_identical':[ast.dump(x,include_attributes=False) for x in ast.parse(original).body if isinstance(x,ast.FunctionDef)]==[ast.dump(x,include_attributes=False) for x in ast.parse(actual).body if isinstance(x,ast.FunctionDef)],
        'new_native_calls':0,'new_initializations':0,'new_factors':0,'compiles':0,'links':0,'submissions':0,
        'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    require(result['original_run_functions_AST_identical'],'Native run AST identity');save(out/'RESULT.json',result)
    print(json.dumps({'result':pin(out/'RESULT.json'),'manifest':mp,'build_configuration':result['build_configuration'],'native_calls':0}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);prepare(p.parse_args().out)
