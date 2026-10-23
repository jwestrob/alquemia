"""Prepare exact immutable import metadata and unbound runtime configuration."""
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('physical_fit_adapter',HERE/'ADAPTER.py')
A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)


def prepare(out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    inputs=json.loads((HERE/'DESIGN_INPUTS.json').read_text());oldcfg=A.read(inputs['old_config']);oldinv=A.read(oldcfg['reuse_inventory'])
    newroot=out.parent;declaration=A.pin(HERE/'ROOT_DECLARATION.json');execution=A.pin(HERE/'ROOT_EXECUTION_CONTEXT.json')
    full_root=Path(inputs['pilot']['path']).parent.parent/'ROOT_ACTUAL_REVIEW.json'
    authorities={'original_inventory':oldcfg['reuse_inventory'],'component_fit':inputs['completed_fit'],'component_review':inputs['completed_fit_review'],'full_pilot_root':A.pin(full_root)}
    registry=A.R.build(authorities,oldcfg['parent_configuration'],oldcfg['parent_evaluator'])
    pilot=A.read(A.R.read(authorities['full_pilot_root'])['result'])
    inv={'status':'prepared_immutable_exact_reuse_index','declaration':declaration,'execution_context':execution,'parent_configuration':oldcfg['parent_configuration'],
         'authorities':authorities,'registry':registry,'origin_key':oldinv['origin_key'],'derivative_keys':oldinv['derivative_keys'],
         'pilot_key':oldinv['pilot_key'],'full_pilot_key':pilot['candidate']['key'],
         'counts':{'seed_vectors':27,'origin':1,'initial_probes':24,'additional_pilots':2,'available_exact_vectors':len(registry)},
         'full_pilot_physical_failure_preserved':True,'mutable_ledger_imported':False,'old_vector_bodies_read_for_index':0,'new_native_calls':0}
    A.save(out/'REUSE_INVENTORY.json',inv)
    deps=json.loads((HERE/'CORE_DEPENDENCIES.json').read_text())
    cfg={k:oldcfg[k] for k in ('parent_evaluator','parent_configuration','decomposition','interpreter','hostname','limits','reserved_outcomes_authorized')}
    cfg.update(implementation=A.pin(HERE/'DRIVER.py'),adapter=A.pin(HERE/'ADAPTER.py'),core=A.pin(HERE/'CORE.py'),reuse_implementation=A.pin(HERE/'REUSE.py'),
        declaration=declaration,execution_context=execution,reuse_inventory=A.pin(out/'REUSE_INVENTORY.json'),journal_implementation=A.pin(HERE/'JOURNAL.py'),
        persistence_policy=A.J.POLICY,solver_environment=inputs['conic_environment'],saved_initial_proposal=deps['saved_initial_proposal'],
        output_directory=str(newroot/'run_v1/execution'),root_release=None,source_review=None,preparation=None)
    A.save(out/'CONFIG_PREPARE.json',cfg)
    result={'status':'prepared_exact_reuse_metadata','implementation':A.pin(__file__),'design_inputs':A.pin(HERE/'DESIGN_INPUTS.json'),'registry':A.pin(out/'REUSE_INVENTORY.json'),
            'configuration':A.pin(out/'CONFIG_PREPARE.json'),'counts':inv['counts'],'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1},'native_calls':0,'solver_calls':0}
    A.save(out/'RESULT.json',result);print(json.dumps({'result':A.pin(out/'RESULT.json'),'configuration':result['configuration'],'counts':inv['counts']}))


if __name__=='__main__':
    import sys
    prepare(sys.argv[1])
