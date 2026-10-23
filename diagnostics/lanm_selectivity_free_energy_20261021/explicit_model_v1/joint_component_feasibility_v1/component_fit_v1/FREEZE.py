#!/usr/bin/env python3
"""Bind already completed source preparation; never execute a fit."""
import importlib.util
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('freeze_policy',HERE/'ADAPTER.py');A=importlib.util.module_from_spec(s);s.loader.exec_module(A)
W=HERE.parents[4]/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1'
cfg=json.loads((W/'inputs_v2/CONFIG_PREPARE.json').read_text());prep=A.pin(W/'preparation_v2/RESULT.json');p=A.read(prep)
A.require(p['status']=='prepared_component_fit_zero_native' and p['configuration_fingerprint']==A.fingerprint(cfg) and p['native_programs']==0 and p['real_proposal_calls']==0,'Final exact zero-call preparation')
for k in ('implementation','adapter','core'):A.verify(cfg[k])
cc=A.read(A.pin(HERE/'CORE_READY.json'));A.require(cc['implementation']==cfg['core'] and cc['source_frozen'],'Final frozen numerical core')
controls=A.pin(W/'execution_controls_v2/RESULT.json');t=A.read(controls);A.require(t['adapter']==cfg['adapter'] and t['status']=='passed_component_execution_controls' and t['native_programs']==0,'Final execution controls')
cfg['preparation']=prep
A.require(not (HERE/'READY.json').exists(),'Freeze once')
A.save(HERE/'CONFIG_RUN_TEMPLATE.json',cfg)
keys=('implementation','adapter','core','declaration','execution_context','parent_configuration','reuse_inventory','source_review','preparation','limits','output_directory')
release={k:cfg[k] for k in keys};release.update(status='released_bounded_component_fit',configuration_fingerprint=A.fingerprint(cfg),cache_directory=cfg['output_directory']+'/cache',reserved_outcomes_authorized=False)
A.save(HERE/'ROOT_RELEASE_TEMPLATE.json',release)
config_run=W/'run_v1/CONFIG_RUN.json';release_run=W/'run_v1/ROOT_RELEASE.json'
ready={'status':'component_fit_source_prepared_pending_root_release','implementation':cfg['implementation'],'adapter':cfg['adapter'],'core':cfg['core'],'configuration_fingerprint':A.fingerprint(cfg),
       'declaration':cfg['declaration'],'execution_context':cfg['execution_context'],'reuse_inventory':cfg['reuse_inventory'],'preparation':prep,'preparation_configuration':A.pin(W/'inputs_v2/CONFIG_PREPARE.json'),
       'execution_controls':controls,'core_ready':A.pin(HERE/'CORE_READY.json'),'core_controls':cc['controls'],'configuration_template':A.pin(HERE/'CONFIG_RUN_TEMPLATE.json'),'release_template':A.pin(HERE/'ROOT_RELEASE_TEMPLATE.json'),
       'inventory_preparer':A.pin(HERE/'PREPARE_INPUTS.py'),'controls_source':A.pin(HERE/'CONTROLS.py'),'implementation_notes':A.pin(HERE/'IMPLEMENTATION.md'),'freeze_source':A.pin(__file__),
       'runtime_binding':{'source_review_status':'passed_component_fit_source_review','source_review_fields':['implementation','adapter','core','configuration_fingerprint'],
                          'configuration_path':str(config_run),'release_path':str(release_run),'instructions':'Bind source_review in config and release. Write release, pin it as config.root_release, then write config. Preserve all scientific fields and fingerprint. Root checks live aggregate login resources before release.'},
       'argv':[cfg['interpreter']['path'],cfg['implementation']['path'],'run','--config',str(config_run),'--out',cfg['output_directory']],
       'environment':{k:'1' for k in (*A.THREADS,'BLIS_NUM_THREADS','PYTHONDONTWRITEBYTECODE')},'limits':A.LIMITS,
       'reuse_counts':p['import_counts'],'initial_metrics':p['initial_metrics'],'native_programs_executed':0,'actual_optimizer_calls':0,'no_automatic_retry_or_extension':True,
       'remaining_gate':'Independent exact source/preparation review and explicit root execution release. No actual optimizer or native call has run.'}
A.save(HERE/'READY.json',ready);print(json.dumps({'ready':A.pin(HERE/'READY.json'),'configuration_template':ready['configuration_template'],'release_template':ready['release_template'],'preparation':prep}))
