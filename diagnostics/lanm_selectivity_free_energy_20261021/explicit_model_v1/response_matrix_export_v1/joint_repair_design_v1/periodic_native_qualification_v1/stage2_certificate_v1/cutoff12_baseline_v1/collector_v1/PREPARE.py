"""Metadata-only unbound templates and source freeze; no actual endpoint reads."""
import argparse
import ast
from pathlib import Path
import COLLECT as C

HERE=Path(__file__).resolve().parent

def templates():
    ready=C.read(C.pin(HERE.parent/'READY.json'));m=C.read(ready['manifest']);decl=C.read(ready['root_declaration'])
    C.require(ready['manifest']['sha256']=='1ef94367e66a3dd564515852bac4282bc4db06b3a09fad7d3231ac0ccf66aab9' and tuple(m['rows'])==C.IDS,'Exact frozen source panel')
    base=HERE.parent.parent/'baseline_broad_control_v1/collector_v1/COLLECT.py'
    inp={'manifest':ready['manifest'],'root_declaration':ready['root_declaration'],'parent_collector':m['collector'],'component_helper':C.pin(base),'cutoff9_collection':m['baseline_failed_collection'],'cutoff9_actual_review':decl['evidence'],'directions':m['directions'],'native_source_ready':C.pin(HERE.parent/'READY.json')}
    C.save(HERE/'INPUTS.json',inp)
    cfg={'status':'bound_cutoff12_baseline_collection','implementation':C.pin(C.__file__),'inputs':C.pin(HERE/'INPUTS.json'),'manifest':ready['manifest'],'build':None,'build_root_review':None,'submissions':None,'scheduler_accounting':None,'build_scheduler_accounting':None,'results':{sid:None for sid in C.IDS}}
    C.save(HERE/'CONFIG_TEMPLATE.json',cfg)
    C.save(HERE/'ROOT_RELEASE_TEMPLATE.json',{'status':'released_cutoff12_baseline_collection','implementation':C.pin(C.__file__),'configuration':None,'source_review':None,'output':None})
    print('Unbound templates prepared; no native/collection calls.')

def freeze(controls):
    inp=C.read(C.pin(HERE/'INPUTS.json'));control_pin=C.pin(controls);r=C.read(control_pin)
    C.require(r['status']=='passed_cutoff12_collector_manufactured_controls' and r['collector']==C.pin(C.__file__) and r['implementation']==C.pin(HERE/'TEST_POLICY.py') and r['count']==len(r['checks']) and all(x['passed'] for x in r['checks']),'Exact passing manufactured controls')
    C.require(r['native_calls']==r['actual_collections']==r['scientific_matrix_solves']==0,'No scientific execution')
    for name in ('COLLECT.py','TEST_POLICY.py','PREPARE.py'):ast.parse((HERE/name).read_text())
    cfg=C.read(C.pin(HERE/'CONFIG_TEMPLATE.json'));C.require(cfg['implementation']==C.pin(C.__file__) and cfg['inputs']==C.pin(HERE/'INPUTS.json') and cfg['manifest']==inp['manifest'],'Template freeze joins')
    for item in inp.values():C.verify(item)
    ready={'status':'ready_cutoff12_baseline_collector_source_only','implementation':C.pin(C.__file__),'inputs':C.pin(HERE/'INPUTS.json'),'specification':C.pin(HERE/'SPEC.txt'),'preparer':C.pin(__file__),'controls_source':C.pin(HERE/'TEST_POLICY.py'),'controls':control_pin,'configuration_template':C.pin(HERE/'CONFIG_TEMPLATE.json'),'release_template':C.pin(HERE/'ROOT_RELEASE_TEMPLATE.json'),'native_source_ready':inp['native_source_ready'],'manifest':inp['manifest'],'required_source_review':{'status':'passed_cutoff12_baseline_collector_source_review','implementation':C.pin(C.__file__),'inputs':C.pin(HERE/'INPUTS.json')},'counts':{'new_gradient_results_required':9,'new_central_results_required':1,'original_derivative_checks':6,'directions':2,'steps_per_direction':2},'bindings_pending':['build','build_root_review','nine_terminal_results','nine_ordered_submissions','scheduler_accounting','build_scheduler_accounting','source_review','root_release'],'run_argv':['/home/jwestrob/.pyenv/versions/3.12.0/bin/python3',str(HERE/'COLLECT.py'),'ROOT_BOUND_CONFIG','ROOT_BOUND_RELEASE','FRESH_COLLECTION_OUTPUT'],'environment':{k:'1' for k in C.THREADS},'actual_cutoff12_outcomes_read':0,'native_calls':0,'actual_collections':0,'scope':'Source package only; root owns collection release. No active-response or physical qualification.'}
    C.save(HERE/'READY.json',ready);print(C.pin(HERE/'READY.json'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('templates','freeze'));p.add_argument('--controls',type=Path);a=p.parse_args()
    if a.mode=='templates':templates()
    else:freeze(a.controls)
