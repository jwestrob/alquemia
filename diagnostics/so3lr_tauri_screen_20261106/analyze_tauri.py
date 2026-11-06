"""Saved-array analysis of the frozen real-fixture panel; no model calls."""
import os,csv,argparse
import numpy as np
from common import *
assert os.environ.get('SLURM_JOB_ID')
ap=argparse.ArgumentParser();ap.add_argument('--result',default=str(W/'panel_cpu_v1/RESULT.json'));args=ap.parse_args()
source=read(args.result);panel=read(W/'PANEL.json')['rows'];by={r['id']:r for r in panel}
startup=read(Path(args.result).parent/'START.json')
assert startup['panel']==pin(W/'PANEL.json') and startup['environment']==pin(D/'environment-lock.txt')
pred={(r['model'],r['id']):r for r in source['rows']}

def write_tsv(name,rows):
 with (RESULTS_DIR/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def contrast(energy):
 return dict(C_W=(energy['aquo9_central_Dy']-energy['aquo8_central_Dy'])-(energy['aquo9_central_La']-energy['aquo8_central_La']),
             C_Ac=(energy['hydrated_acetate_Dy']-energy['aquo8_central_Dy'])-(energy['hydrated_acetate_La']-energy['aquo8_central_La']))

def readout(r,forces):
 """Forces in kcal/mol/A. Positive radial force moves the donor outward."""
 x=np.asarray(r['positions_A']);f=np.asarray(forces);out={}
 donors=[i for i,z in enumerate(r['Z']) if z==8]
 for i in donors:
  u=x[i]-x[0];u/=np.linalg.norm(u)
  out[f'O{i}_distance_A']=float(np.linalg.norm(x[i]-x[0]))
  out[f'O{i}_radial_force_kcal_mol_A']=float(f[i]@u)
  out[f'O{i}_metal_toward_donor_kcal_mol_A']=float(f[0]@u)
 if r['kind'] in ('aquo','hydrated_acetate'):
  first=1 if r['kind']=='aquo' else 8
  assert all(r['Z'][i:i+3]==[8,1,1] for i in range(first,len(x),3))
  out['shell_dE_dR_kcal_mol_A']=float(sum(-f[i:i+3].sum(0)@((x[i]-x[0])/np.linalg.norm(x[i]-x[0])) for i in range(first,len(x),3)))
 if r['kind'] in ('dry_acetate','hydrated_acetate'):
  center=by['dry_acetate_'+('Dy' if 66 in r['Z'] else 'La')]
  fr=center['historical_reference']['mechanical_modes']['analytic']['frame']
  ii=fr['ligand_indices'];rad=np.asarray(fr['radial_unit']);axis=np.asarray(fr['rotation_unit']);cent=np.asarray(fr['centroid_A'])
  out['ligand_radial_dE_dq_kcal_mol_A']=-float(f[ii].sum(0)@rad)
  tangent=np.cross(axis,x[ii]-cent)/fr['rotation_length_A']
  out['ligand_angular_dE_dq_kcal_mol_A']=-float((f[ii]*tangent).sum())
 return out

energies=[];local=[];deltas=[];historical=[];contrasts=[];sensitivities=[];dry_curves=[]
for model in ('l','m'):
 for r in panel:
  p=pred[model,r['id']];assert p['charge']==r['charge'] and p['multiplicity']==r['multiplicity']
  assert pin(p['raw_arrays']['path'])==p['raw_arrays']
  with np.load(p['raw_arrays']['path']) as raw:
   assert raw['energy'].dtype==np.float64 and raw['forces'].dtype==np.float64
   assert float(raw['energy'][0])==p['energy_eV']
   assert np.array_equal(raw['forces'][:len(r['Z'])],p['forces_eV_A'])
  energies.append(dict(model=model,id=r['id'],charge=r['charge'],multiplicity=r['multiplicity'],energy_eV=p['energy_eV'],
                       energy_kcal_mol=p['energy_eV']*EV_KCAL,cutoff_A=1000,source_npz_sha256=p['raw_arrays']['sha256']))
  ro=readout(r,np.asarray(p['forces_eV_A'])*EV_KCAL)
  for key,value in ro.items():local.append(dict(model=model,id=r['id'],observable=key,value=value))
  h=r.get('historical_reference') or {};hf=h.get('forces_kcal_mol_A',h.get('forces'))
  if hf is not None:
   hr=readout(r,hf)
   for key,value in ro.items():
    if not key.endswith('distance_A'):
     historical.append(dict(model=model,id=r['id'],observable=key,Tauri=value,historical=hr[key],difference=value-hr[key],
                            status='descriptive_only_unqualified_PBE0_D4_reference',Dy_state='effective_singlet_ECP55' if 66 in r['Z'] else 'La_singlet_ECP46'))
  if 'La' in r['id']:
   other=by[r['id'].replace('La','Dy')]
   assert r['positions_A']==other['positions_A'] and r['Z'][1:]==other['Z'][1:] and r['charge']==other['charge']
   dp=pred[model,other['id']];df=(np.asarray(dp['forces_eV_A'])-np.asarray(p['forces_eV_A']))*EV_KCAL
   dr=readout(r,df)
   deltas.append(dict(model=model,La_id=r['id'],Dy_id=other['id'],coordinates_atom_order_state_matched_except_metal=True,
                      charge=r['charge'],multiplicities=[1,6],forces_Dy_minus_La_kcal_mol_A=df.tolist(),
                      observables={k:v for k,v in dr.items() if not k.endswith('distance_A')}))
 # Same six total-energy states, differenced before any rounding.
 states=['aquo8_central_La','aquo8_central_Dy','aquo9_central_La','aquo9_central_Dy','hydrated_acetate_La','hydrated_acetate_Dy']
 for cut in (1000.,100.,None):
  suff='' if cut==1000 else '_cutoff_'+str(cut)
  en={cid:pred[model,cid+suff]['energy_eV']*EV_KCAL for cid in states}
  c=contrast(en);contrasts.append(dict(model=model,LR_cutoff_A=cut,contrasts_kcal_mol=c,
     six_significant_digit_export_kcal_mol=contrast({cid:float(format(pred[model,cid+suff]['energy_eV'],'.6g'))*EV_KCAL for cid in states})))
 for met in ('La','Dy'):
  for direction in ('radial','rotation'):
   for step in (.005,.01):
    lab='005' if step==.005 else '010'
    a=pred[model,f'pbe0__{met}__{direction}__minus_{lab}__complex']['energy_eV']
    b=pred[model,f'pbe0__{met}__{direction}__plus_{lab}__complex']['energy_eV']
    ctr=by['dry_acetate_'+met];ro=readout(ctr,np.asarray(pred[model,ctr['id']]['forces_eV_A'])*EV_KCAL)
    key='ligand_radial_dE_dq_kcal_mol_A' if direction=='radial' else 'ligand_angular_dE_dq_kcal_mol_A'
    dry_curves.append(dict(model=model,metal=met,direction=direction,step_A=step,central_finite_difference_kcal_mol_A=(b-a)*EV_KCAL/(2*step),
                           center_analytic_kcal_mol_A=ro[key],historical_reference_status='Dy_effective_singlet_not_physical_sextet; all_old_roots_unqualified'))
 base=pred[model,'aquo8_central_La'];bf=np.asarray(base['forces_eV_A']);we=pred[model,'spectator_water_free']['energy_eV']
 for distance in (50.,100.):
  p=pred[model,'spectator_'+str(distance)];df=(np.asarray(p['forces_eV_A'])[:len(bf)]-bf)*EV_KCAL
  sensitivities.append(dict(model=model,test='real_water_spectator',distance_A=distance,
      nonadditive_energy_kcal_mol=(p['energy_eV']-base['energy_eV']-we)*EV_KCAL,
      max_base_atom_force_change_kcal_mol_A=float(abs(df).max()),metal_force_change_kcal_mol_A=df[0].tolist()))
 a=pred[model,'translated_dry_fragment'];b=pred[model,'translated_dry_unconstrained']
 sensitivities.append(dict(model=model,test='separated_real_metal_acetate_charge_conditioning',separation_translation_A=50,
      explicit_fragment_minus_plain_energy_kcal_mol=(a['energy_eV']-b['energy_eV'])*EV_KCAL,
      max_force_difference_kcal_mol_A=float(abs(np.asarray(a['forces_eV_A'])-np.asarray(b['forces_eV_A'])).max()*EV_KCAL)))

checks=[x for x in source['checks'] if 'test' in x]
assert len(energies)==48 and len(deltas)==24
assert all(x.get('passed',False) for x in checks),checks
write_tsv('energies.tsv',energies);write_tsv('force_readouts.tsv',local);write_tsv('historical_force_comparison.tsv',historical);write_tsv('dry_curves.tsv',dry_curves)
delta_path=W/'DELTA_FORCES.json' if RESULTS_DIR==D else RESULTS_DIR/'DELTA_FORCES.json'
save(delta_path,dict(rows=deltas,definitions='Dy minus La; exact common coordinates and charge, physical multiplicities 6 and 1'))
summary=dict(source=pin(args.result),panel=pin(W/'PANEL.json'),job=os.environ['SLURM_JOB_ID'],
    primary_counts=dict(configurations=24,models=2,predictions=48,identical_coordinate_Dy_La_pairs_per_model=12),
    implementation_checks=checks,balanced_contrasts=contrasts,conditioning_sensitivity=sensitivities,
    throughput=[x for x in source['timings'] if x['stage']=='warm_throughput'],
    model_load_and_compile_calls=[x for x in source['timings'] if x['stage']=='load' or (x['stage']=='evaluate' and x['wall_seconds']>1)],
    archived_delta_forces=pin(delta_path),all_implementation_checks_pass=True,
    scientific_status=('Closed reference attempt; see REFERENCE_COMPARISON.json. No MD released.' if (RESULTS_DIR/'REFERENCE_RESULTS.json').exists()
       else 'Await physical-spin references; old PBE0-D4 data descriptive only. No MD released.'))
save(RESULTS_DIR/'TAURI_RESULTS.json',summary)
for c in contrasts:print(c['model'],c['LR_cutoff_A'],c['contrasts_kcal_mol'])
for x in summary['throughput']:print('warm CPU',x['model'],x['mean_seconds'],'seconds/evaluation')
print('all checks pass',len(checks))

refpath=RESULTS_DIR/'REFERENCE_RESULTS.json'
if refpath.exists():
 refs=read(refpath);rr={r['id']:r for r in refs['rows']};fr=[];balanced=[];primary_wrong=[]
 for model in ('l','m'):
  for cid,r in rr.items():
   if not r.get('engrad_complete'):continue
   pr=readout(by[cid],np.asarray(pred[model,cid]['forces_eV_A'])*EV_KCAL)
   qr=readout(by[cid],r['forces_kcal_mol_A'])
   for key in pr:
    if key.endswith('distance_A'):continue
    opposite=pr[key]*qr[key]<0 and abs(qr[key])>1
    radial=key.endswith('_radial_force_kcal_mol_A') or key=='shell_dE_dR_kcal_mol_A'
    stop=bool(model=='l' and r['qualified'] and radial and opposite)
    row=dict(model=model,id=cid,observable=key,Tauri=pr[key],native_reference=qr[key],difference=pr[key]-qr[key],
       reference_qualified=r['qualified'],provisional_primary_force_gate_failed=stop)
    fr.append(row)
    if stop:primary_wrong.append(row)
  if all(r.get('energy_Eh') is not None for r in rr.values()):
   rc=contrast({cid:r['energy_Eh']*EH_KCAL for cid,r in rr.items()})
   tc=next(c for c in contrasts if c['model']==model and c['LR_cutoff_A']==1000.)['contrasts_kcal_mol']
   for name,ids in [('C_W',['aquo8_central_La','aquo8_central_Dy','aquo9_central_La','aquo9_central_Dy']),
                    ('C_Ac',['aquo8_central_La','aquo8_central_Dy','hydrated_acetate_La','hydrated_acetate_Dy'])]:
    qualified=all(rr[cid]['qualified'] for cid in ids)
    error=tc[name]-rc[name]
    balanced.append(dict(model=model,contrast=name,Tauri_kcal_mol=tc[name],native_reference_kcal_mol=rc[name],
      difference_kcal_mol=error,reference_qualified=qualified,provisional_budget_kcal_mol=1.,
      status=('reject' if abs(error)>1 else 'within_provisional_budget') if qualified else 'unresolved_reference'))
  else:
   tc=next(c for c in contrasts if c['model']==model and c['LR_cutoff_A']==1000.)['contrasts_kcal_mol']
   for name in ('C_W','C_Ac'):
    balanced.append(dict(model=model,contrast=name,Tauri_kcal_mol=tc[name],native_reference_kcal_mol=None,
      difference_kcal_mol=None,reference_qualified=False,provisional_budget_kcal_mol=1.,status='unresolved_reference'))
 primary=[r for r in balanced if r['model']=='l']
 if primary_wrong or any(r['status']=='reject' for r in primary):decision='reject_potential'
 elif len(primary)==2 and all(r['reference_qualified'] for r in primary):decision='screen_pass_only_no_solvent_release_without_experimental_admission'
 else:decision='unresolved_reference'
 write_tsv('reference_force_comparison.tsv',fr)
 save(RESULTS_DIR/'REFERENCE_COMPARISON.json',dict(reference=pin(refpath),Tauri=pin(RESULTS_DIR/'TAURI_RESULTS.json'),
     balanced_contrasts=balanced,qualified_primary_wrong_force_responses=primary_wrong,screen_decision=decision,
     experimental_admission=read(D/'EXPERIMENTAL_ADMISSION.json')['admitted_pair'],MD_released=False,
     interpretation='Unqualified numerical differences are descriptive and cannot establish potential failure or success.'))
 print('reference comparison',decision,balanced)
