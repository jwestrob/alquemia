"""Summarize real failed/qualified gates and make a saved-data diagnostic figure."""
import os,argparse,re
from util import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ase.io import read as ase_read
assert os.environ.get('SLURM_JOB_ID')
ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
out=args.output_dir;out.mkdir(exist_ok=False,parents=True)
for p in [Path(__file__),D/'util.py',D/'CHEMICAL_INTEGRITY_AUDIT.md',D/'PREPARATION_COMPARISON.md']:(out/p.name).write_bytes(p.read_bytes())
controls=read(W/'saved_checks_65911/CONTROL_SUMMARY.json')
audit=read(W/'integrity_audit_66010/RESULT.json')
scout=read(W/'md_scout_analysis_65972/RESULT.json')
capped=read(W/'preparation_comparison_ASE_v1/RESULT.json')
bridge=read(W/'preparation_comparison_v1/BRIDGE.json')
full=[r for r in controls['periodic'] if r['backend']=='periodic_gpu_v2']
hydrated=next(r for r in controls['periodic'] if r['backend']=='periodic_hydrated_gpu_v1')
assert len(full)==3 and all(r['passes_energy_criterion'] for r in full) and hydrated['passes_energy_criterion']
assert all(r['passed'] for r in controls['cpu_GPU_agreement'])
assert all(audit[k] for k in ['atom_order_verified','initial_coordinates_verified','checkpoint_HDF_coordinates_verified'])
assert len(audit['altered_bonds'])==3 and {r['ligand_id'] for r in audit['altered_bonds']}=={600}
assert all(r['first_bad_phase']=='relaxation' and r['first_bad_step']==20 for r in audit['altered_bonds'])
assert capped['status']=='chemical_integrity_failure' and capped['rows'][-1]['step']==37
bad=capped['rows'][-1]['ligands'][0]['bad_bonds'];assert len(bad)==1 and (bad[0]['atom_a'],bad[0]['atom_b'])==(6,7)
assert max(r['max_atom_step_A'] for r in capped['rows'])<=.05000001 and not bridge['passed']
cap_dir=W/'preparation_comparison_ASE_v1';cap_frames=ase_read(cap_dir/'optimizer.traj',index=':')
assert len(cap_frames)==len(capped['rows'])
prep=read(W/'hydrated_preparation_v1/RESULT.json')['rows'][0]
initial=ase_read(prep['A']['path']);graphs={r['ligand_id']:r for r in read(W/'chemistry.json')['rows']}
previous=initial.positions
for frame,row in zip(cap_frames,capped['rows']):
 assert np.array_equal(frame.numbers,initial.numbers)
 dx=frame.positions-previous;dx-=42*np.round(dx/42)
 assert abs(float(np.linalg.norm(dx,axis=1).max())-row['max_atom_step_A'])<1e-10
 for fragment in prep['fragments']:
  lid=int(fragment['ligand_id']);ids=np.array(fragment['ligand_indices']);pairs=ids[np.array(graphs[lid]['bonds'])]
  dx=frame.positions[pairs[:,0]]-frame.positions[pairs[:,1]];dx-=42*np.round(dx/42)
  dd=np.linalg.norm(dx,axis=1);saved=next(r for r in row['ligands'] if r['ligand_id']==lid)
  assert abs(float(dd.max())-saved['bond_max_A'])<1e-10 and abs(float(dd.min())-saved['bond_min_A'])<1e-10
 previous=frame.positions
md=[r for r in scout['rows'] if r['phase']=='MD'];assert len(md)==200 and max(r['time_ps'] for r in md)==.25
native_log=W/'md_scout_v1/native.log';text=native_log.read_text()
avg=float(re.findall(r'Average time per step:\s+(\S+) seconds',text)[-1])
warm=read(W/'periodic_gpu_v2/RESULT.json')['warm_seconds'];assert len(warm)==5
fig,axes=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
colors={600:'#b54431',599:'#247a78'}
labels={600:'trans CDTA (La)',599:'cis CDTA (Dy)'}
for lid in [600,599]:
 rows=[r for r in audit['rows'] if r['phase']=='relaxation' and r['ligand_id']==lid]
 axes[0].plot([r['step'] for r in rows],[r['bond_max_A'] for r in rows],color=colors[lid],label=labels[lid])
 axes[1].plot([r['step'] for r in capped['rows']],[next(l['bond_max_A'] for l in r['ligands'] if l['ligand_id']==lid) for r in capped['rows']],color=colors[lid])
 rows=[r for r in md if int(r['ligand_id'])==lid]
 axes[2].plot([r['time_ps'] for r in rows],[r['bond_max_A'] for r in rows],color=colors[lid])
for ax in axes:
 ax.axhline(2,color='#666666',ls='--',lw=1);ax.grid(alpha=.2);ax.set_ylabel('Longest original covalent bond (Å)')
axes[0].set(title='Native FIRE preparation',xlabel='Optimization step');axes[0].legend(frameon=False,fontsize=8)
axes[1].set(title='ASE LBFGS, 0.05 Å step cap',xlabel='Optimization step')
axes[2].set(title='MD after the failed FIRE preparation',xlabel='MD time (ps)')
fig.savefig(out/'CHEMICAL_INTEGRITY.png',dpi=170);fig.savefig(out/'CHEMICAL_INTEGRITY.pdf');plt.close(fig)
payload=dict(decision='do_not_advance_current_Tauri_setup',classification='validation_insufficient_for_intact_ligand_selectivity',
 reason='The target trans ligand fails the declared covalent-integrity test under both attempted preparation paths; the MD trajectory begins with altered chemistry.',
 global_potential_rejection_claim=False,selectivity_free_energy_kcal_mol=None,
 original_periodic_max_swap_changes_kcal_mol=[r['maximum_swap_change_kcal_mol'] for r in full],
 hydrated_max_swap_change_kcal_mol=hydrated['maximum_swap_change_kcal_mol'],
 CPU_GPU_agreement=controls['cpu_GPU_agreement'],
 native_ASE_force_bridge=dict(max_difference_eV_A=bridge['max_force_difference_eV_A'],tolerance_eV_A=1e-5,passed=False),
 native_integrity=audit['altered_bonds'],capped_integrity=dict(step=37,bad_bonds=bad,maximum_atom_step_A=max(r['max_atom_step_A'] for r in capped['rows'])),
 capped_raw_trajectory_verified=True,capped_raw_trajectory_frames=len(cap_frames),
 completed_MD_ps=.25,independent_MD_starts=1,planned_starts_not_run=5,
 ASE_full_box_warm_median_seconds=float(np.median(warm)),ASE_full_box_warm_range_seconds=[float(min(warm)),float(max(warm))],
 native_MD_reported_average_seconds_per_step=avg,native_MD_steps=1000,preparation_steps=200,
 MD_used_for_affinity=False,qualified_Dy_reference_available=False,experimental_common_medium_admitted=False,
 intact_ensemble_sampling_demonstrated=False,production_baseline_changed=False,
 further_MD_released=False,new_DFT_calls=0,model_calls=0,analysis_job=os.environ['SLURM_JOB_ID'],
 sources=[pin(W/p) for p in ['saved_checks_65911/CONTROL_SUMMARY.json','integrity_audit_66010/RESULT.json','md_scout_analysis_65972/RESULT.json','preparation_comparison_ASE_v1/RESULT.json','preparation_comparison_v1/BRIDGE.json']])
save(out/'FINAL_RESULTS.json',payload)
print('Saved-data final review: stop current setup; original failed gates preserved. Figure and numerical summary written.',flush=True)
