"""Read existing coordinates, native parameters and scout forces; no energy calls."""
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
w=ROOT/'workspaces/metal_environment_response_20260926'
l=read_json(w/'hybrid_preparation_v1/ledger_v1/LEDGER.json')
a=read_json(verify(l['artifacts']['particles.json']));terms=read_json(verify(l['artifacts']['bonded_terms.json']))
r=read_json(w/'component_checks_v1/results/Ca_A/receipt.json');x=np.load(verify(r['task']['coordinates']));f=np.load(verify(r['forces']))['forces_kcal_mol_A']
ids={v['id']:i for i,v in enumerate(a)}
def summary(v):return dict(n=len(v),minimum=float(np.min(v)),median=float(np.median(v)),maximum=float(np.max(v)))
bonds=[];angles=[]
for t in terms:
 ix=t['atoms']
 if t['force']=='HarmonicBondForce':
  i,j=ix;dist=float(np.linalg.norm(x[i]-x[j]));r0=float(t['parameters'][0].split()[0])*10
  bonds.append(dict(ids=[a[k]['id'] for k in ix],elements='-'.join(sorted(a[k]['element'] for k in ix)),distance_A=dist,equilibrium_A=r0,residual_A=dist-r0,source_index=t['index']))
 elif t['force']=='HarmonicAngleForce':
  i,j,k=ix;u=x[i]-x[j];v=x[k]-x[j];angle=np.degrees(np.arccos(np.clip(u@v/(np.linalg.norm(u)*np.linalg.norm(v)),-1,1)));eq=np.degrees(float(t['parameters'][0].split()[0]))
  angles.append(dict(ids=[a[k]['id'] for k in ix],angle_deg=float(angle),equilibrium_deg=float(eq),residual_deg=float(angle-eq),includes_H=any(a[k]['element']=='H' for k in ix)))
bond_stats={kind:dict(actual_A=summary([v['distance_A'] for v in bonds if v['elements']==kind]),residual_A=summary([v['residual_A'] for v in bonds if v['elements']==kind])) for kind in sorted({v['elements'] for v in bonds})}
angle_stats={str(h):summary([abs(v['residual_deg']) for v in angles if v['includes_H']==h]) for h in (True,False)}
top_forces={}
for gi,name in enumerate(('retained_bonded','MM_LJ','MM_Coulomb','QM_MM_LJ')):
 norm=np.linalg.norm(f[gi],axis=1);ii=np.argsort(norm)[-15:][::-1]
 top_forces[name]=[dict(id=a[i]['id'],magnitude_kcal_mol_A=float(norm[i])) for i in ii]
original=ROOT/'workspaces/mxaf_qm/1H4I.pdb';protonated=ROOT/'workspaces/mxaf_qm/1H4I_protonated.pdb'
def pdb(path):
 result={}
 for s in path.read_text().splitlines():
  if s.startswith(('ATOM  ','HETATM')) and s[21]=='A':
   key=f"A/{int(s[22:26])}/{s[26]}/{s[12:16].strip()}"
   result[key]=dict(xyz_A=[float(s[k:k+8]) for k in (30,38,46)],element=s[76:78].strip(),resname=s[17:20].strip())
 return result
orig=pdb(original);prot=pdb(protonated)
common=[k for k in orig.keys()&prot.keys() if orig[k]['element']!='H']
added=[k for k in prot if k not in orig and prot[k]['element']!='H']
examples=[]
exceptions={tuple(sorted(e['atoms'])):e for e in read_json(verify(l['artifacts']['exceptions.json']))}
for one,two in [('A/595/ /CE','A/596/ /O'),('A/599/ /CB','A/599/ /O'),('A/595/ /CG','A/596/ /CA')]:
 i,j=ids[one],ids[two];examples.append(dict(ids=[one,two],distance_A=float(np.linalg.norm(x[i]-x[j])),source_PDB_distance_A=float(np.linalg.norm(np.asarray(prot[one]['xyz_A'])-prot[two]['xyz_A'])),native_exception=exceptions.get(tuple(sorted((i,j)))),distance_to_QM_metal_A=[float(np.linalg.norm(x[k]-x[ids['metal']])) for k in (i,j)]))
result=dict(status='read_only_geometry_and_saved_force_diagnosis',inputs=dict(ledger=record(w/'hybrid_preparation_v1/ledger_v1/LEDGER.json'),scout=record(w/'component_checks_v1/results/Ca_A/receipt.json'),original_PDB=record(original),protonated_PDB=record(protonated)),bond_statistics=bond_stats,angle_absolute_residual_degrees=angle_stats,top_bond_residuals=sorted(bonds,key=lambda v:abs(v['residual_A']),reverse=True)[:20],top_angle_residuals=sorted(angles,key=lambda v:abs(v['residual_deg']),reverse=True)[:20],saved_top_forces=top_forces,clash_examples=examples,source_comparison=dict(common_heavy_atom_count=len(common),max_common_heavy_displacement_A=max(float(np.linalg.norm(np.asarray(orig[k]['xyz_A'])-prot[k]['xyz_A'])) for k in common),added_heavy_atoms=added,missing_residue_remarks=[s for s in original.read_text().splitlines() if s.startswith('REMARK 465') and (' A   59' in s)]),energy_force_calls=0)
write_new(Path(__file__).with_name('SOURCE_DIAGNOSIS.json'),result)
print(json.dumps({k:v for k,v in result.items() if k in ('bond_statistics','angle_absolute_residual_degrees','clash_examples','source_comparison')},indent=2))
