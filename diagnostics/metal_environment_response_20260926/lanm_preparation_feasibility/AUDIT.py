"""Read-only audit/proposed source membership for three consumed LanM structures."""
import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
manifest=ROOT/'workspaces/lanm_global_occupancy_20260923/prepared_v1/manifest.json';m=read_json(manifest)
repairs={
 'Hans_8DQ2':'workspaces/benchmark_set_20260915/prepared/hans_lanm_v1/EF3/amide_v3/repair_manifest.json',
 'Hans_8FNR':'workspaces/lanm_series_followup_20260923/dy_transfer_v1/sources/EF3/amide_v3/repair_manifest.json',
 'Mex_8FNS':'workspaces/lanm_series_followup_20260923/sources_v2/Mex_EF3/amide_v3/repair_manifest.json'}
def stats(v):return dict(count=len(v),min=float(min(v)),median=float(np.median(v)),max=float(max(v)))
results=[]
for s in m['sources']:
 data=read_json(verify(s['atoms']));audit=read_json(verify(s['protein_audit']));pm=read_json(verify(s['protonation_receipt']))
 verify(s['crystal']);verify(s['protonated_source']);a=data['atoms'];ids={v['id']:i for i,v in enumerate(a)}
 x=np.array([v['xyz_A'] for v in a]);old=np.array([v['source_xyz_A'] for v in a]);ef3=np.array(s['sites']['EF3']['xyz_A'])
 lo,hi=(83,94) if s['protein']=='Hans' else (84,95)
 graph=defaultdict(set);groups=defaultdict(list);donorbonds=[]
 for b in data['bonds']:
  i,j=b['indices'];graph[i].add(j);graph[j].add(i)
  kind='-'.join(sorted([a[i]['element'],a[j]['element']]))
  before=np.linalg.norm(old[i]-old[j]);after=np.linalg.norm(x[i]-x[j]);eq=b['forcefield_equilibrium_A']
  groups[kind].append((before,after,after-eq))
  if a[i]['element']!='H' and a[j]['element']!='H' and any(lo<=a[k]['source']['resid']<=hi for k in (i,j)):
   donorbonds.append(dict(ids=[a[k]['id'] for k in (i,j)],distance_A=float(after),equilibrium_A=eq))
 bondstats={k:dict(source_A=stats([v[0] for v in vv]),prepared_A=stats([v[1] for v in vv]),prepared_residual_A=stats([v[2] for v in vv])) for k,vv in groups.items()}
 rp=ROOT/repairs[s['source_id']];r=read_json(rp)
 oldcore=[]
 for v in r['atom_graph']['source_to_qm']:
  if v['kind']!='source':continue
  q=v['source'];identifier=f"{q['chain']}/{q['resnum']}/{q['insertion_code']}/{q['atom']}"
  if identifier not in ids:raise ValueError('core source absent '+identifier)
  oldcore.append(ids[identifier])
 seed=[i for i in oldcore if a[i]['element'] in ('O','N')]
 loop={i for i,v in enumerate(a) if v['kind']=='protein_source' and lo<=v['source']['resid']<=hi}
 # Preserve full bounding peptide units; only flanking CA-C and N-CA sigma cuts.
 flanking=[f'A/{lo-1}//C',f'A/{lo-1}//O',f'A/{hi+1}//N',f'A/{hi+1}//H']
 if any(v not in ids for v in flanking):raise ValueError('unsupported flanking peptide unit')
 selected=loop|{ids[v] for v in flanking}
 # All real waters already in direct core, plus O within3.4A of complete core N/O.
 support=[]
 for i,v in enumerate(a):
  if v['element'] not in ('O','N') or i in oldcore:continue
  dd=np.linalg.norm(x[seed]-x[i],axis=1);jj=int(np.argmin(dd))
  if dd[jj]<=3.4:
   support.append(dict(id=v['id'],kind=v['kind'],distance_A=float(dd[jj]),core_contact=a[seed[jj]]['id'],inside_loop=i in loop))
 water_res={a[i]['source']['resid'] for i in oldcore if a[i]['kind']=='retained_crystal_water'}
 water_res.update(a[ids[v['id']]]['source']['resid'] for v in support if v['kind']=='retained_crystal_water')
 for i,v in enumerate(a):
  if v['kind']=='retained_crystal_water' and v['source']['resid'] in water_res:selected.add(i)
 outside=[v for v in support if v['kind']=='protein_source' and ids[v['id']] not in selected]
 extra_caps=[];extra_amides=[]
 for contact in outside:
  oi=ids[contact['id']]
  if a[oi]['source']['name']!='O':raise ValueError('unsupported nonlocal Hbond support group')
  carbons=[j for j in graph[oi] if a[j]['element']=='C' and a[j]['source']['name']=='C']
  if len(carbons)!=1:raise ValueError('ambiguous carbonyl connectivity')
  ci=carbons[0];nn=[j for j in graph[ci] if a[j]['element']=='N']
  if len(nn)!=1:raise ValueError('missing actual amide nitrogen')
  ni=nn[0];hh=[j for j in graph[ni] if a[j]['element']=='H']
  selected.update([oi,ci,ni,*hh])
  capids=[]
  for retained in (ci,ni):
   ca=[j for j in graph[retained] if a[j]['source']['name']=='CA']
   if len(ca)!=1:raise ValueError('unsupported nonlocal amide sigma boundary')
   capids.append(dict(retained=a[retained]['id'],omitted=a[ca[0]]['id']))
  extra_caps.extend(capids);extra_amides.append([a[j]['id'] for j in [ci,oi,ni,*hh]])
 # Distance-selected support contacts are not asserted hydrogen bonds.
 nearest=[]
 heavy=[i for i,v in enumerate(a) if v['element']!='H']
 for ni,i in enumerate(heavy):
  excluded={i}|graph[i]|set(j for k in graph[i] for j in graph[k])
  for j in heavy[ni+1:]:
   if j in excluded:continue
   d=np.linalg.norm(x[i]-x[j])
   if d<1.8:nearest.append(dict(ids=[a[i]['id'],a[j]['id']],distance_A=float(d)))
 spectator=[dict(site=k,source_element=v['source']['element'],distance_from_EF3_A=float(np.linalg.norm(np.asarray(v['xyz_A'])-ef3)),source=v['source'],representation='not_yet_assigned') for k,v in s['sites'].items() if k!='EF3']
 # One deterministic intact-sidechain perturbation proposal, no coordinates generated.
 target=85 if s['protein']=='Hans' else 86
 target_ids=[v['id'] for v in a if v['kind']=='protein_source' and v['source']['resid']==target]
 rot_ids=[v for v in target_ids if v.split('/')[-1] in ('CG','OD1','OD2','HB2','HB3')]
 results.append(dict(source_id=s['source_id'],pins={k:s[k] for k in ('atoms','protein_audit','crystal','protonated_source','protonation_receipt','water_inventory')},old_EF3_graph=record(rp),protonation={k:pm[k] for k in ('ph','add_missing_residues','repaired_missing_atom_count','repaired_missing_terminal_atom_count')},protein_residues=s['protein_residues'],protein_atoms=s['protein_atoms'],water_count=s['water_count'],protein_charge=s['protein_charge'],H_bond_statistics=bondstats,prepared_heavy_max_displacement_A=float(max(np.linalg.norm(x[i]-old[i])for i in heavy)),EF3_heavy_bonds=donorbonds,nonbonded_heavy_contacts_below_1p8_A=sorted(nearest,key=lambda q:q['distance_A']),direct_donors=r['coordination']['typed_donors'],spectators=spectator,existing_core_H_coordinate_max_change_A=float(max(np.linalg.norm(x[i]-old[i]) for i in oldcore if a[i]['element']=='H')),support_contacts=support,proposal=dict(status='source_membership_only_not_prepared_QM_input',rule='complete EF3 12-residue loop; flanking peptide units; complete direct-core waters and waters within3.4A of complete core N/O',loop_residues=[lo,hi],source_atoms=[dict(physical_index=i,id=a[i]['id']) for i in sorted(selected)],caps=[dict(retained=f'A/{lo-1}//C',omitted=f'A/{lo-1}//CA'),dict(retained=f'A/{hi+1}//N',omitted=f'A/{hi+1}//CA')]+extra_caps,nonlocal_amide_units=extra_amides,water_residues=sorted(water_res),real_nonmetal_atoms=len(selected),proposed_QM_atoms_with_metal_and_caps=len(selected)+3+len(extra_caps),outside_loop_protein_support=outside,formal_region_charge_proposal=(-1 if s['protein']=='Hans' else 0),charge_multiplicity='formal bookkeeping proposal; verify electrons after actual preparation; La physical1 Dy physical6 only if supported explicit state',perturbation=dict(type='intact carboxylate torsion',axis=[f'A/{target}//CA',f'A/{target}//CB'],rotate_source_atoms=rot_ids,right_hand_angle_degrees=2.,status='proposal only; validate actual full covalent graph and contacts before any scoring'))))
result=dict(status='read_only_LanM_source_and_region_feasibility',manifest=record(manifest),sources=results,energy_force_calls=0,new_preparations=0,reserved_outcomes_accessed=False)
write_new(Path(__file__).with_name('RESULT_v2.json'),result)
for r in results:
 print(r['source_id'], 'Hdelta',r['existing_core_H_coordinate_max_change_A'],'clashes',r['nonbonded_heavy_contacts_below_1p8_A'],'proposed atoms',r['proposal']['proposed_QM_atoms_with_metal_and_caps'],'water',r['proposal']['water_residues'],'other support',r['proposal']['outside_loop_protein_support'])
