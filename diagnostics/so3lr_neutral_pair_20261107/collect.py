"""Saved-only analysis of declared controls. No model, fit or molecular calls."""
import os, re, subprocess
from util import *
import numpy as np
assert os.environ.get('SLURM_JOB_ID')
report=Path(os.environ.get('TAURI_REPORT_DIR',str(D)))
report.mkdir(exist_ok=True,parents=True)

local=read(D/'LOCAL_RESULTS.json')
for name in ('predictions','frames'):
 assert pin(local[name]['path'])==local[name]
assert local['model_calls']==128 and abs(local['FD']['residual'])<.003
local_summary=dict(
 preparations=3,model_calls=local['model_calls'],
 maximum_spectator_energy_change_kcal_mol=max(abs(r['water_effect_kcal_mol']) for r in local['rows']),
 maximum_force_difference_component_change_kcal_mol_A=max(r['max_force_difference_change_kcal_mol_A'] for r in local['rows']),
 maximum_trans_metal_force_difference_vector_change_kcal_mol_A=max(r['trans_metal_force_difference_change_kcal_mol_A'] for r in local['rows']),
 directional_derivative_check=local['FD'],passes_energy_criterion=local['status']=='pass')
assert local_summary['maximum_spectator_energy_change_kcal_mol']==local['max_water_effect_kcal_mol']

prep=read(W/'preparation_v1/RESULT.json')
periodic=[]
raw_cache={}
def arrays(row):
 p=row['raw_arrays'];assert pin(p['path'])==p
 a=dict(np.load(p['path']))
 assert a['forces_eV_A'].dtype==np.float64 and np.isfinite(a['forces_eV_A']).all()
 assert a['positions_A'].shape==a['forces_eV_A'].shape==(row['natoms'],3)
 assert abs(a['partial_charges'].sum())<1e-6 and row['state_head_verified']
 return a
for tag in ('periodic_cpu_v2','periodic_cpu_controls_v1','periodic_gpu_v1'):
 p=W/tag/'RESULT.json'
 complete=p.exists()
 if tag=='periodic_gpu_v1' and not complete:continue
 if complete:
  result=read(p)
 else:
  p=W/tag/'PARTIAL.json'
  if not p.exists():continue
  result=read(p)
  result['job']=read(W/tag/'START.json')['job']
  sel=W/tag/'REMOTE_SELECTIONS.json'
  result['remote_selections']=read(sel) if sel.exists() else []
 by={row['id']:arrays(row) for row in result['rows']}
 raw_cache[tag]=by
 for start in range(3):
  a,b=f's{start}_A',f's{start}_B'
  if a not in by or b not in by:continue
  A,B=by[a],by[b]
  assert np.array_equal(A['positions_A'],B['positions_A'])
  assert sorted(A['Z'])==sorted(B['Z'])
  assert np.flatnonzero(A['Z']!=B['Z']).tolist()==[0,43]
  gap=float(B['energy_eV']-A['energy_eV'])*EV_KCAL
  df=(B['forces_eV_A']-A['forces_eV_A'])*EV_KCAL
  row=dict(backend=tag,start=start,job=result['job'],job_result_complete=complete,
     fixed_configuration_B_minus_A_kcal_mol=gap,
     free_energy_estimate=False,translation=None,remote_waters=[])
  ta,tb=f's{start}_translated_A',f's{start}_translated_B'
  if ta in by and tb in by:
   tg=float(by[tb]['energy_eV']-by[ta]['energy_eV'])*EV_KCAL
   row['translation']=dict(swap_energy_change_kcal_mol=tg-gap,
      maximum_endpoint_energy_change_kcal_mol=max(abs(float(by[t]['energy_eV']-by[s]['energy_eV']))*EV_KCAL for t,s in [(ta,a),(tb,b)]),
      maximum_force_difference_component_change_kcal_mol_A=float(abs((by[tb]['forces_eV_A']-by[ta]['forces_eV_A'])*EV_KCAL-df).max()))
  sel=next((s for s in result['remote_selections'] if s['start']==start),None)
  if sel:
   for j,dist in zip(sel['water_ids'],sel['solute_clearances_A']):
    ra,rb=f's{start}_remove{j}_A',f's{start}_remove{j}_B'
    if ra not in by or rb not in by:continue
    inds=[i for i in range(len(A['Z'])) if i not in (88+3*j,89+3*j,90+3*j)]
    assert np.array_equal(by[ra]['Z'],A['Z'][inds]) and np.array_equal(by[rb]['Z'],B['Z'][inds])
    assert np.array_equal(by[ra]['positions_A'],A['positions_A'][inds])
    rg=float(by[rb]['energy_eV']-by[ra]['energy_eV'])*EV_KCAL
    delta=(by[rb]['forces_eV_A']-by[ra]['forces_eV_A'])*EV_KCAL-df[inds]
    row['remote_waters'].append(dict(water_id=j,solute_heavy_atom_clearance_A=dist,
      swap_energy_change_kcal_mol=rg-gap,
      maximum_force_difference_component_change_kcal_mol_A=float(abs(delta).max()),
      metal_force_difference_vector_changes_kcal_mol_A=[float(np.linalg.norm(delta[k])) for k in (0,43)]))
  row['primary_robustness_criterion_kcal_mol']=.5
  row['complete_controls']=row['translation'] is not None and len(row['remote_waters'])==3
  if row['complete_controls']:
   row['maximum_swap_change_kcal_mol']=max(abs(row['translation']['swap_energy_change_kcal_mol']),*(abs(r['swap_energy_change_kcal_mol']) for r in row['remote_waters']))
   row['passes_energy_criterion']=row['maximum_swap_change_kcal_mol']<=.5
  periodic.append(row)

agreement=[]
if 'periodic_gpu_v1' in raw_cache:
 cpu=raw_cache['periodic_cpu_v2'];gpu=raw_cache['periodic_gpu_v1']
 for cid in ('s0_A','s0_B'):
  if cid in gpu:
   ce=float(abs(cpu[cid]['energy_eV']-gpu[cid]['energy_eV']))
   cf=float(abs(cpu[cid]['forces_eV_A']-gpu[cid]['forces_eV_A']).max())
   agreement.append(dict(id=cid,energy_difference_eV=ce,max_force_difference_eV_A=cf,
     passed=ce<=1e-5 and cf<=1e-5))
summary=dict(local=local_summary,periodic=periodic,
 cpu_GPU_agreement=agreement,analysis_job=os.environ['SLURM_JOB_ID'],molecular_calls=0,
 interpretation='Robustness and implementation checks only. Unrelaxed endpoint gaps are not selectivity predictions. Remote-water removal can have a physical electrostatic response.')
save(report/'CONTROL_SUMMARY.json',summary)
original=read(W/'saved_checks_52991/CONTROL_SUMMARY.json')
assert original['local']==summary['local']
for old in original['periodic']:
 assert old==next(r for r in periodic if (r['backend'],r['start'])==(old['backend'],old['start']))
save(report/'REPRODUCTION_CHECKS.json',dict(passed=True,
 original=pin(W/'saved_checks_52991/CONTROL_SUMMARY.json'),
 exact_local_summary=True,exact_CPU_periodic_summaries=len(original['periodic']),
 molecular_calls=0,synthetic_fixtures=0))
jobs=['52753','52797','52830','52877','52948','52844',os.environ['SLURM_JOB_ID']]
s=subprocess.check_output(['sacct','-n','-P','-j',','.join(jobs),
 '--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS'],text=True)
(report/'ACCOUNTING.txt').write_text(s)
save(report/'ACCOUNTING.json',dict(rows=[dict(zip(['job','state','elapsed_seconds','allocated_cpus','measured_CPU','max_RSS'],line.split('|')[:6])) for line in s.splitlines() if line.split('|')[0] in jobs],
 note='The current collector may still be RUNNING in this snapshot. GPU allocation time is zero while pending; refresh after completion.'))
print('Saved-only analysis complete:',local_summary,'periodic controls',[(r['backend'],r.get('maximum_swap_change_kcal_mol')) for r in periodic])
