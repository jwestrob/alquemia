"""Reuse normalized 52-atom 1F6S with exact archived physical field, no energy calls."""
import argparse,copy,json,math,sys
from pathlib import Path
import numpy as np

def main():
 p=argparse.ArgumentParser();p.add_argument('--repository',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.repository.resolve();sys.path.insert(0,str(root/'scripts'))
 from affordable_common import read_json,record,verify,write_new
 from affordable_response import cap_jacobians
 from mace_responsive_field import pointcharge_text
 from metal_environment_force_checks import xyz_data,pc_data
 mp=root/'workspaces/mace_omol_20260917/responsive_quantum_v2/manifest.json';m=read_json(mp);tasks={t['metal']:t for t in m['tasks'] if t['case_id']=='ALPHA_1F6S'};s=read_json(verify(tasks['Ca']['source_state']));case=read_json(verify(s['source_case']));phys={v['id']:v for v in s['physical_atoms']};prep=read_json(verify(s['normalized_preparation']));assert prep['physical_atoms']==s['physical_atoms'];verify(s['forcefield_background']['source']);verify(s['forcefield_background']['forcefield'])
 for z,t in tasks.items():assert verify(t['pointcharges']).read_text()==pointcharge_text(t['source_task']);assert t['charge']==(-1 if z=='Ca' else 0) and t['multiplicity']==1
 assert tasks['La']['pointcharges']['sha256']==tasks['Ca']['pointcharges']['sha256']
 _,els,qm=xyz_data(verify(tasks['La']['xyz']));_,ce,cq=xyz_data(verify(tasks['Ca']['xyz']));assert np.array_equal(qm,cq) and els[1:]==ce[1:] and len(qm)==52
 def pid(r):return f"{r['chain']}/{r['resnum']}/{r['insertion_code']}/{r['atom']}"
 graph=case['source_graph'];mapping=[{'qm_index':0,'kind':'metal','id':'metal','element':'La','xyz_A':qm[0].tolist()}];sourceids=set();capchecks=[]
 for v in graph['source_to_qm']:
  i=v['qm_index'];r={'qm_index':i,'element':els[i],'xyz_A':qm[i].tolist()}
  if v['kind']=='source':
   sid=pid(v['source']);actual=phys[sid];assert actual['element']==els[i] and np.max(abs(qm[i]-actual['xyz_A']))<1e-9;sourceids.add(sid);r.update(kind='source',id=sid,source_id=sid,original_source=v['source'],normalization='archived source-normalized coordinates, not original raw H')
  else:
   assert v['kind']=='sigma_link_H';left,right=pid(v['retained']),pid(v['omitted']);x,y=[np.array(phys[k]['xyz_A']) for k in (left,right)];l=v['length_A'];expected=x+l*(y-x)/np.linalg.norm(y-x);err=float(abs(expected-qm[i]).max());assert err<1e-9;ja,jb=cap_jacobians(x,y,l);r.update(kind='cap',id='cap/'+str(i),retained_source_id=left,omitted_source_id=right,retained_xyz_A=x.tolist(),omitted_xyz_A=y.tolist(),length_A=l,jacobian_retained=ja.tolist(),jacobian_omitted=jb.tolist(),source_graph_record=v);capchecks.append(err)
  mapping.append(r)
 mapping.sort(key=lambda v:v['qm_index']);assert [v['qm_index'] for v in mapping]==list(range(52));assert np.max(abs(qm[0]-phys['metal']['xyz_A']))<1e-9
 # Retained peptide chemistry must include the actual source N/H, not carbonyl CH2O.
 amides=[('A/79//C','A/79//O','A/80//N','A/80//H'),('A/84//C','A/84//O','A/85//N','A/85//H')]
 for group in amides:assert set(group)<=sourceids
 background_bonds={frozenset(b) for b in s['forcefield_background']['bonds']}
 for carbon,oxygen,nitrogen,hydrogen in amides:assert all(frozenset(b) in background_bonds for b in [(carbon,oxygen),(carbon,nitrogen),(nitrogen,hydrogen)])
 for resid in (211,212):assert {f'A/{resid}//'+n for n in ('O','H1','H2')}<=sourceids
 weights=read_json(verify(tasks['Ca']['source_task']['weights']));_,q,xyz=pc_data(verify(tasks['Ca']['pointcharges']));assert len(q)==1880 and abs(float(q.sum())+4)<1e-10
 env=[]
 for i,(sid,j,charge) in enumerate(zip(weights['physical_ids'],weights['physical_indices'],weights['weights_e'])):
  v=s['physical_atoms'][j];assert v['id']==sid and np.array_equal(xyz[i],np.array(v['xyz_A'])) and q[i]==charge;assert sid not in set(s['projection_support_ids']);env.append(dict(v,charge_e=charge,physical_index=j))
 assert set(weights['physical_ids'])|set(s['projection_support_ids'])==set(phys)
 byid={v['id']:v for v in env};heavy=qm[[i for i,e in enumerate(els) if e!='H']];candidates=[]
 for v in env:
  spec={'SER':('OG','HG','CB'),'THR':('OG1','HG1','CB'),'TYR':('OH','HH','CZ')}.get(v.get('resname'))
  if not spec or not v['id'].endswith('/'+spec[0]):continue
  prefix=v['id'].rsplit('/',1)[0]+'/';ids=[prefix+n for n in spec]
  if all(sid in byid for sid in ids):candidates.append({'distance_to_core_heavy_A':float(np.linalg.norm(heavy-v['xyz_A'],axis=1).min()),'oxygen_id':ids[0],'hydrogen_id':ids[1],'axis_carbon_id':ids[2]})
 candidates.sort(key=lambda v:(v['distance_to_core_heavy_A'],v['oxygen_id']));sel=candidates[0];assert sel['oxygen_id']=='A/86//OG1';o,h,c=[np.array(byid[sel[k]]['xyz_A']) for k in ('oxygen_id','hydrogen_id','axis_carbon_id')];axis=(o-c)/np.linalg.norm(o-c);r=h-o;t=math.radians(10);hp=o+r*math.cos(t)+np.cross(axis,r)*math.sin(t)+axis*np.dot(axis,r)*(1-math.cos(t));b=copy.deepcopy(env);idx=next(i for i,v in enumerate(b) if v['id']==sel['hydrogen_id']);b[idx]['xyz_A']=hp.tolist();assert abs(np.linalg.norm(hp-o)-np.linalg.norm(h-o))<1e-12
 nonbonded=np.array([v['xyz_A'] for v in env if v['element']!='H' and v['id']!=sel['oxygen_id']]+heavy.tolist());mind=float(np.linalg.norm(nonbonded-hp,axis=1).min());assert mind>=1.2
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'core_mapping.json',mapping);write_new(out/'boundary_mapping.json',{'source_state':tasks['Ca']['source_state'],'field_weights':tasks['Ca']['source_task']['weights'],'boundary_ledger':s['boundary_ledger'],'projection_support_ids':s['projection_support_ids'],'policy':'existing normalized full-boundary graph redistribution, no new charges or PQQ boundary substitution'})
 envs={}
 for label,atoms in [('A',env),('B',b)]:
  write_new(out/f'environment_{label}.json',atoms);pc=out/f'environment_{label}.pc';pc.write_text(str(len(atoms))+'\n'+''.join(' '.join(format(v,'.17g') for v in (x['charge_e'],*x['xyz_A']))+'\n' for x in atoms));envs[label]={'pointcharges':record(pc),'atoms':record(out/f'environment_{label}.json'),'count':len(atoms),'charge_e':sum(x['charge_e'] for x in atoms)}
 assert envs['A']['pointcharges']['sha256']==tasks['Ca']['pointcharges']['sha256']
 result={'protocol_id':'nikasha_1F6S_normalized52_field_response_preparation_v1','source_state':tasks['Ca']['source_state'],'source_manifest':record(mp),'source_case':s['source_case'],'source':s['forcefield_background']['source'],'source_normalized_preparation':s['normalized_preparation'],'core_mapping':record(out/'core_mapping.json'),'boundary_mapping':record(out/'boundary_mapping.json'),'endpoints':{z:{'xyz':t['xyz'],'charge':t['charge'],'multiplicity':1} for z,t in tasks.items()},'environments':envs,'assembly':s['assembly'],'microstate':s['microstate'],'explicit_waters':s['explicit_waters'],'evidence':s['evidence'],'evidence_use':'consumed_method_development_not_blind','perturbation':dict(sel,pointcharge_index=idx,axis_unit=axis.tolist(),pivot_xyz_A=o.tolist(),angle_degrees=10.,hydrogen_A_A=h.tolist(),hydrogen_B_A=hp.tolist(),OH_length_A=float(np.linalg.norm(r)),displacement_A=float(np.linalg.norm(hp-h)),minimum_nonbonded_H_heavy_after_A=mind,selection_rule='nearest complete exterior Ser/Thr/Tyr hydroxyl to QM heavy atom; identity tiebreak; positive10degree carbon-to-oxygen rotation',candidate_ranking=candidates),'mapping_checks':{'all_52_atoms_accounted':True,'actual_source_peptide_amides':amides,'both_3atom_waters_in_QM':True,'cap_count':len(capchecks),'max_cap_coordinate_error_A':max(capchecks),'A_pointcharges_exact_archived_bytes':True},'plan':record(Path(__file__).with_name('PLAN.md')),'implementation':record(__file__),'molecular_evaluations':0,'limitations':['Distinct52atom normalized NMA representation, not original40atom repair.','General peptide boundary redistribution differs from PQQ sidechain N/C scheme.','Thr86OH is4.64A from core; remote weak perturbation is not direct-affinity validation.','Existing source/protein-water H normalization is retained; neither hydration equilibrium nor full hybrid is modeled.'],'numerical_qualification_status':'not_established','full_hybrid_status':'unsupported_complete_cross_and_boundary_terms'}
 write_new(out/'INPUTS.json',result);print(json.dumps({'inputs':str(out/'INPUTS.json'),'mapping':result['mapping_checks'],'perturbation':{k:v for k,v in result['perturbation'].items() if k!='candidate_ranking'}},indent=2))
if __name__=='__main__':main()
