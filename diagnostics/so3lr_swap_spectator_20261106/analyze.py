"""Analyze saved real-frame predictions; no additional model evaluations."""
import os,sys,csv,argparse
from pathlib import Path
D=Path(__file__).resolve().parent;ROOT=D.parents[1]
OLD=ROOT/'diagnostics/so3lr_tauri_screen_20261106'
sys.path.insert(0,str(OLD))
from common import read,save,pin,EV_KCAL
import numpy as np
assert os.environ.get('SLURM_JOB_ID')
ap=argparse.ArgumentParser();ap.add_argument('--output-dir',default=str(D));args=ap.parse_args()
OUT=Path(args.output_dir);OUT.mkdir(parents=True,exist_ok=True)
W=ROOT/'workspaces/so3lr_swap_spectator_20261106';attempt=W/'run_v1'
result=read(attempt/'RESULT.json');frames={r['id']:r for r in read(attempt/'FRAMES.json')['rows']}
entries={(r['model'],r['cutoff_A'],r['id']):r for r in result['rows']}
assert len(entries)==62
components=['nn_energy','zbl_repulsion','electrostatic_energy','dispersion_energy']
for row in result['rows']:
 assert pin(row['raw_arrays']['path'])==row['raw_arrays']

def raw(model,cut,cid):
 with np.load(entries[model,cut,cid]['raw_arrays']['path']) as a:return {k:a[k].copy() for k in a.files}

def table(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

changes=[];overview=[];energy_rows=[]
for model in ('l','m'):
 for cut in (1000.,None):
  A=entries[model,cut,'A'];B=entries[model,cut,'B'];D0=(B['energy_eV']-A['energy_eV'])*EV_KCAL
  qa=raw(model,cut,'A');qb=raw(model,cut,'B');qw=raw(model,cut,'free_spectator_water')['partial_charges'][:3]
  bare_force=qb['forces'][:54]-qa['forces'][:54]
  base_xyz=np.array(frames['A']['positions_A'])
  these=[]
  for distance in (100.,200.,400.):
   for orientation in ('original','reversed'):
    ida=f'A_water_{distance:g}_{orientation}';idb=f'B_water_{distance:g}_{orientation}'
    a=entries[model,cut,ida];b=entries[model,cut,idb];ar=raw(model,cut,ida);br=raw(model,cut,idb)
    Dw=(b['energy_eV']-a['energy_eV'])*EV_KCAL;delta=Dw-D0
    wx=np.array(frames[ida]['positions_A'])[54:]
    distances=np.linalg.norm(base_xyz[:,None,:]-wx[None,:,:],axis=-1)
    ordinary=14.399645351950548*EV_KCAL*float(((qb['partial_charges'][:54]-qa['partial_charges'][:54])[:,None]*qw[None,:]/distances).sum())
    fchange=(br['forces'][:54]-ar['forces'][:54]-bare_force)*EV_KCAL
    parts={k:(b['component_sums_eV'][k]-a['component_sums_eV'][k]-B['component_sums_eV'][k]+A['component_sums_eV'][k])*EV_KCAL for k in components}
    assert abs(sum(parts.values())-delta)<1e-7
    row=dict(model=model,cutoff_A=cut,distance_A=distance,orientation=orientation,
      bare_swap_kcal_mol=D0,with_water_swap_kcal_mol=Dw,water_change_kcal_mol=delta,
      frozen_charge_Coulomb_change_kcal_mol=ordinary,
      max_swap_force_change_kcal_mol_A=float(abs(fchange).max()),
      acetate_metal_swap_force_change_kcal_mol_A=np.linalg.norm(fchange[0]).item(),
      aquo_metal_swap_force_change_kcal_mol_A=np.linalg.norm(fchange[29]).item(),
      spectator_charge_A=a['charge_on_spectator'],spectator_charge_B=b['charge_on_spectator'],
      acetate_charge_shift_A=a['charge_on_acetate_complex']-A['charge_on_acetate_complex'],
      acetate_charge_shift_B=b['charge_on_acetate_complex']-B['charge_on_acetate_complex'],
      **{k+'_change_kcal_mol':v for k,v in parts.items()})
    changes.append(row);these.append(row)
  far=[r for r in these if r['distance_A']>=200]
  peak=max(abs(r['water_change_kcal_mol']) for r in far)
  physical=max(abs(r['frozen_charge_Coulomb_change_kcal_mol']) for r in far)
  status=('pass' if peak<=.5 else 'fail') if physical<=.05 else 'inconclusive_remote_interaction'
  overview.append(dict(model=model,cutoff_A=cut,bare_swap_kcal_mol=D0,
    maximum_far_water_effect_kcal_mol=peak,maximum_frozen_Coulomb_effect_kcal_mol=physical,
    criterion_kcal_mol=.5,status=status,
    distance_change_200_to_400_kcal_mol={o:next(r['water_change_kcal_mol'] for r in these if r['distance_A']==400 and r['orientation']==o)-next(r['water_change_kcal_mol'] for r in these if r['distance_A']==200 and r['orientation']==o) for o in ('original','reversed')}))
for r in result['rows']:
 energy_rows.append({k:r[k] for k in ['model','cutoff_A','id','natoms','charge','multiplicity','energy_eV','evaluation_seconds']})
table(OUT/'swap_changes.tsv',changes);table(OUT/'energies.tsv',energy_rows)
primary=[r for r in overview if r['model']=='l'];assert len(primary)==2
robust=primary[0]['status']==primary[1]['status'] and primary[0]['status'] in ('pass','fail')
decision=primary[0]['status'] if robust else 'inconclusive_cutoff_or_remote_interaction'
out=dict(decision=decision,scope='One distant neutral water, charged isolated acetate/aquo swap prototype; no chemical accuracy or neutral periodic-box validation.',
 primary_model='so3lr-tauri-l-2606',sensitivity_model='so3lr-tauri-m-2606',
 overview=overview,criterion_kcal_mol=.5,rows=changes,
 implementation_checks=dict(all_frames_state_coordinates_head_verified=True,complete_LR_pairs=True,
   float64_finite=True,charge_conservation=True,energy_component_closure=True,
   saved_LaW8_replays=[r for r in result['checks'] if r.get('test')=='saved_real_LaW8_replay']),
 native_predictions=pin(attempt/'RESULT.json'),frame_manifest=pin(attempt/'FRAMES.json'),plan=pin(D/'PLAN.md'),
 job=os.environ['SLURM_JOB_ID'],model_calls=62,DFT_calls=0,MD_steps=0)
save(OUT/'RESULTS.json',out)
for r in overview:print('SUMMARY',r,flush=True)
print('DECISION',decision,'for the declared prototype only',flush=True)
