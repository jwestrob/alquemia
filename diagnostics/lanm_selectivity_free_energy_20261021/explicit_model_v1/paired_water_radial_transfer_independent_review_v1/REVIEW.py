#!/usr/bin/env python3
"""Freeze the narrow conditional-adapter review; do not read radial outcomes."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[4]
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_radial_transfer_v1'
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
OUT=W/'paired_water_radial_transfer_independent_review_v1'
checks=[]
def pin(path):
    path=Path(path).resolve();return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
def check(ok,name):
    checks.append({'check':name,'passed':bool(ok)})
    if not ok:raise ValueError(name)

tick=time.monotonic();OUT.mkdir(parents=True,exist_ok=False)
check(pin(D/'RUN.py')['sha256']=='67a9052b29914cb82b13e669dd9190fc7f3429c95cd88e48af7b6df72dad7d6b','Final reviewed source')
pins=json.loads((D/'PINS.json').read_text())
for name,p in pins['pins'].items():check(pin(p['path'])==p,'Protected pin '+name)
controls_path=W/'paired_water_radial_transfer_v1/saved_controls_v1/RESULT.json'
controls=json.loads(controls_path.read_text())
check(pin(controls_path)['sha256']=='8cbd4fa5ddf42ab7d9094eb267a6db5aee59fca44a0c4c39fcd8465d610075ab','Frozen narrow-control receipt')
check(controls['run_source']==pin(D/'RUN.py') and controls['implementation']==pin(D/'TEST_BINDING.py'),'Controls bind reviewed sources')
check(controls['status']=='passed' and len(controls['records'])==31 and all(x['passed'] for x in controls['records']),'31 narrow controls all passed')
check(controls['new_native_calls']==0 and not controls['candidate_bound'] and not controls['radial_execution_outcome_opened'] and not controls['reserved_validation_opened'],'Controls remain saved-only and unbound')
sp=importlib.util.spec_from_file_location('conditional_transfer_review',D/'RUN.py');r=importlib.util.module_from_spec(sp);sp.loader.exec_module(r)
check(all(getattr(r,k) is getattr(r.old,k) for k in ('loaded','gradient','program_admission','replay_native')),'Native/parser helpers unchanged')
data=json.loads(Path(pins['pins']['prepared32']['path']).read_text())
check(set(r.WATER)|set(r.d.HYDRATED)|set(r.d.DRY)==set(data['rows']) and len(data['rows'])==32,'Complete32 denominator')
check(len(r.WATER)==8 and len(r.d.HYDRATED)==18 and len(r.d.DRY)==6 and len(r.d.CENTERS)==4,'Changed memberships')
check(len(r.d.CENTERS)+2*len(r.d.HYDRATED)==40 and len(r.d.FD)==10 and sum(6*data['rows'][g]['geometry']['atoms'] for g in r.d.FD)==1380,'Independent finite inventory calculation')
ready=json.loads((D/'READY.json').read_text())
old=json.loads(Path(pins['pins']['previous_saved_controls']['path']).read_text())
check(old['status']=='passed' and len(old['records'])==49 and old['new_native_calls']==0,'Original49 parser/arithmetic controls reused without rerun')
source=(D/'RUN.py').read_text()
check("for gid in d.RADIAL:cell(gid)" in source and source.index("'Eight water admission failed'")<source.index('pool.map(cell,d.HYDRATED)'),'Eight-water admission precedes donor stage')
check("result['native_program_invocations']<40" in source and "result['internal_FD_energies_requested']<=1380" in source,'Native ceilings enforced')
check("verify(reuse['cartesian_check_authority'])" in source and "g==saved_fd and g['fd_passed']" in source,'Separate saved Cartesian authority and gate retained')
check("explicit_diagnostic_candidate_decision" in source and "all(bound['candidate_flags'][k]" in source,'Development-failing candidate retains explicit root decision')
result={'status':'passed_conditional_source_review','review_source':pin(__file__),'implementation':pin(D/'RUN.py'),'plan':pin(D/'PLAN.md'),'pin_inventory':pin(D/'PINS.json'),'ready':pin(D/'READY.json'),'saved_controls':pin(controls_path),'previous_controls':pins['pins']['previous_saved_controls'],'checks':checks,'passed_checks':len(checks),'resolved_finding':'Whole-family equality replaced by separate checks of original names/bounds/epsilons, radial start/spec and original bound origin.','reviewed_authority_branches':['new complete radial native cache','exact17578 starting-candidate baseline','identical previous donor input with correct analytic/Cartesian evidence'],'manual_source_findings':['Selected radial objective uses pinned radial training; original four-center score replay is separate.','Separate prior Cartesian evidence is reused without replacing selected analytic model authority.','Only represented-metal water-O pair lines change; original DATA/SCORE/native helpers remain pinned.','First native/Cartesian failure stops new calls; locked result persistence preserves sibling receipts and32slots.','Original centers, all eight water development cells,18 consumed donors and6 unchanged dry controls remain separately labeled.','Bounds and fixed epsilons remain unchanged; no final-point rank is inferred from an earlier Jacobian.'],'no_pending_blocker':True,'candidate_bound':False,'radial_execution_outcome_opened':False,'reserved_validation_opened':False,'new_native_calls':0,'new_submissions':0,'local_threads':1,'elapsed_seconds':time.monotonic()-tick,'command':sys.argv,'remaining_gate':'Actual candidate binding/raw replay and root release are still required; this source review authorizes no native run.'}
(OUT/'REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'checks':len(checks),'result':pin(OUT/'REVIEW.json'),'elapsed_seconds':result['elapsed_seconds']}))
