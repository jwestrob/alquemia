"""Bound template and closed-gate checks only; no fit/new reference outcomes."""
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('fresh_reference',HERE/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
start,cpu=time.monotonic(),time.process_time()
base=m.W/'explicit_model_v1/water_coordination_radial_validation_v1/reference'
out=base/'saved_controls_v1';out.mkdir(parents=True,exist_ok=False)
prep=base/'preparation_v1';manifest=json.loads((prep/'MANIFEST.json').read_text());cfg=json.loads((HERE/'CONFIG_SOURCE.json').read_text());checks=[]
def check(ok,name):m.require(ok,name);checks.append(name)
panel=json.loads(m.verify(manifest['native_panel']).read_text())
check(manifest['implementation']==m.pin(HERE/'RUN.py') and manifest['configuration']==m.pin(HERE/'CONFIG_SOURCE.json'),'Prepared source/config frozen')
check(set(manifest['input_templates'])==set(m.CELLS) and len(panel['cells'])==4,'Exact four templates and native cells')
for c in panel['cells']:
    template=m.verify(manifest['input_templates'][c['cell_id']]['input']).read_text()
    central=json.loads(m.verify(cfg['central_references'][c['metal']]).read_text())
    check(template.replace('@FLOOR_0.75_GRANTED_MIB_DIV_4@','3000')==m.verify(central['input']).read_text(),'Exact central deck '+c['cell_id'])
    check(m.verify(manifest['input_templates'][c['cell_id']]['geometry']).read_bytes()==m.verify(c['geometry']['xyz']).read_bytes(),'Exact radial geometry bytes '+c['cell_id'])
    check(c['explicit_electrons']==98 and c['multiplicity']==1 and c['charge']==3,'Effective state '+c['cell_id'])
for metal,a in manifest['anchor_reference'].items():
    check(a['energy_stages']['COSX_grid_settings_signature'] and a['energy_stages']['analytic_gradient_stage'] and len(a['energy_stages']['final_exchange_change_hartree'])==1,'Saved central native stage metadata '+metal)
opened=[];oldread=m.read
def traced(p):opened.append(p['path']);return oldread(p)
m.read=traced
try:
    m.gate(HERE/'CONFIG_SOURCE.json',prep/'MANIFEST.json',prep/'ROOT_RELEASE.template.json')
except ValueError as exc:
    check(str(exc)=='Exact root release required','Unreleased template rejected before candidate binding')
else:raise AssertionError('Unreleased execution was admitted')
check(not any('/paired_water_coordination_joint_v1/run_v1/' in p for p in opened),'No fit trajectory/result read')
check(all(manifest[k] is None for k in ('candidate_result','selected_review','root_release')),'Prepared candidate/review/release remain null')
receipt={'status':'saved_reference_controls_passed','implementation':m.pin(HERE/'RUN.py'),'test_source':m.pin(__file__),'configuration':m.pin(HERE/'CONFIG_SOURCE.json'),'prepared_manifest':m.pin(prep/'MANIFEST.json'),'checks':checks,'check_count':len(checks),'native_calls':0,'candidate_or_new_validation_outcomes_read':False,'local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
p=m.save(out/'RESULT.json',receipt);print(json.dumps({'result':p,'checks':len(checks),'wall':receipt['wall_seconds'],'CPU':receipt['CPU_seconds']}))
