#!/usr/bin/env python3
"""Finite saved component differences only; no raw numerical reparse/model/solver call."""
import hashlib,json,math,time
from pathlib import Path
from decimal import Decimal
ROOT=Path(__file__).resolve().parents[5]
# Resolve project root without relying on a mutable working directory.
while not (ROOT/'diagnostics').is_dir() or not (ROOT/'workspaces').is_dir():ROOT=ROOT.parent
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021';W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021';REL=Path('explicit_model_v1/radial_drive_shape_v1/terminal_trimer_regression_v1');OUT=Path(__file__).resolve().parent

def pin(p):return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(q):
 p=Path(q['path']);assert pin(p)==q;return json.loads(p.read_text())
def val(x):return Decimal(str(x))
def main():
 tick,cpu=time.monotonic(),time.process_time();s=json.loads((D/REL/'SOURCE_INPUTS.json').read_text());trpin=pin(W/REL/'run_v1/execution/RESULT.json');tr=read(trpin);inv=read(tr['inventory']);sel=read(tr['candidate']);feca=read(s['feca_result']);fsel=read(feca['selected'])
 assert tr['status']=='radial14_consumed_trimer_regression_complete' and tr['actual_counts']['admitted_native_induction_blocks']==12
 current=dict(sel['cells'],**tr['rows']);old=dict(fsel['cells'],**feca['rows']);definitions=dict({k:v['definition'] for k,v in inv['reused_selected'].items()},**inv['rows']);rows={};raw=[]
 for cid,cycle in tr['cycles'].items():
  ids={k:c for k,c in cycle['model_terms'].items() if k in definitions};names=set().union(*(set(current[k]['moments']['components']) for k in ids));deltas={}
  for label in ('original_pre_response_baseline','feca12_physical_gate_terminal'):
   baseline={k:(definitions[k]['baseline'] if label=='original_pre_response_baseline' else old[k]['moments']) for k in ids}
   delta={name:sum((val(coef)*(val(current[k]['moments']['components'].get(name,{}).get('energy_kcal_mol',0))-val(baseline[k]['components'].get(name,{}).get('energy_kcal_mol',0))) for k,coef in ids.items()),Decimal(0)) for name in sorted(names)}
   total=val(tr['comparison']['cycles'][cid]['models']['frozen_radial14_terminal5d992']['E3_kcal_mol'])-val(tr['comparison']['cycles'][cid]['models'][label]['E3_kcal_mol'])
   closure=total-sum(delta.values());assert abs(closure)<Decimal('0.000001') # Report original printed closure; this is an arithmetic sanity guard, not an accuracy gate.
   deltas[label]={'component_differences_kcal_mol':{k:str(v) for k,v in delta.items()},'total_E3_difference_kcal_mol':str(total),'total_minus_component_sum_kcal_mol':str(closure),'nonpolar_component_sum_kcal_mol':str(sum(v for k,v in delta.items() if k!='Polarization'))}
  rows[cid]={'model_terms':cycle['model_terms'],'changed_terms':ids,'unchanged_static_terms':[k for k in cycle['model_terms'] if k not in ids],'difference_from':deltas}
  for k in ids:raw.append(current[k]['moments']['analyze_output'])
 evidence={k:s[k] for k in ('terminal_candidate','feca_result','parent_evaluator','native_source_manifest','original_trimer_comparison','original_trimer_native_admission')};evidence.update(terminal_trimer_result=trpin,scientific_inventory=tr['inventory'],feca_selected=feca['selected'])
 result={'status':'saved_component_attribution_not_new_native_admission','implementation':pin(Path(__file__).resolve()),'evidence':evidence,'rows':rows,'current_raw_output_pins':list({x['sha256']:x for x in raw}.values()),'scope':'Only saved printed component differences; static22 cancel unchanged. No recomputation of energy, induced moments, solver or raw admission. Polarization includes all native/custom response and its relaxation.','new_native_calls':0,'new_solver_calls':0,'new_reference_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
 (OUT/'SAVED_COMPONENTS.json').write_text(json.dumps(result,indent=2)+'\n')
 for cid,row in rows.items():print(cid,[(k,v['component_differences_kcal_mol']['Polarization'],v['nonpolar_component_sum_kcal_mol']) for k,v in row['difference_from'].items()])
if __name__=='__main__':main()
