#!/usr/bin/env python3
"""Single mock positive-path regression for the collector status/denominator fix."""
import argparse,importlib.util,json,time
from pathlib import Path
s=importlib.util.spec_from_file_location('collector_local_cc',Path(__file__).with_name('REFERENCE.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();tick,cpu=time.monotonic(),time.process_time();p=m.load(a.panel);a.out.mkdir(parents=True,exist_ok=False)
def fake(path,c,p,panel):
    r=p['controls'][c['canonical_id']]
    return dict(status='complete',energy_hartree=r['energy_hartree'],RHF_energies_hartree=[r['RHF_energy_hartree']],aux_C={'synthetic_same':True},allocated_cpu_seconds=0.)
m.admit_result=fake
inv={c['cell_id']:'synthetic-not-an-actual-result' for c in p['cells'] if c['stage']=='controls'};m.b.write(a.out/'INVENTORY.json',inv)
r=m.collect(a.panel,a.out/'INVENTORY.json',a.out/'synthetic_result','controls')
checks={'all18_rows_complete':len(r['rows'])==r['complete']==r['declared']==18 and all(x['status']=='complete' for x in r['rows'].values()),'35_controls_pass':len(r['checks'])==35 and all(r['checks'].values()),'no_molecular_calls':r['new_native_calls']==0,'qualified_mock_status':r['status']=='local_controls_passed'}
m.b.check(all(checks.values()),'Collector mock regression');m.b.write(a.out/'RESULT.json',dict(status='passed_saved_collector_control',implementation=m.b.pin(m.__file__),test=m.b.pin(__file__),panel=m.b.pin(a.panel),checks=checks,new_native_calls=0,synthetic_result_only=True,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,local_threads=1))
print(json.dumps(checks))
