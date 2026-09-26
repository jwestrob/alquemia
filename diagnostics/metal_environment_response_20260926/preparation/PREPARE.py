"""Prepare real fixed-core field-response inputs; no molecular evaluations."""
from pathlib import Path
import argparse, copy, hashlib, json, math, sys
import numpy as np

def rec(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def load(p):return json.loads(Path(p).read_text())
def verify(r):
 assert rec(r['path'])==r,r['path'];return Path(r['path'])
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def rotate(v,k,t):return v*math.cos(t)+np.cross(k,v)*math.sin(t)+k*np.dot(k,v)*(1-math.cos(t))
def main():
 a=argparse.ArgumentParser();a.add_argument('--repository',type=Path,required=True);a.add_argument('--output',type=Path,required=True);a=a.parse_args();root=a.repository.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 sys.path.insert(0,str(root/'scripts'));from affordable_response import cap_jacobians
 mpath=root/'workspaces/density_embedding_20260916/prepared_v1/manifest.json';m=load(mpath);tasks={t['metal']:t for t in m['tasks'] if t['case']=='1h4i_qm36'}
 states={z:load(verify(t['source_state'])) for z,t in tasks.items()};s=states['La'];c=states['Ca'];assert s['environment_atoms']==c['environment_atoms'];assert load(verify(s['boundary_mapping']))==load(verify(c['boundary_mapping']));assert s['core_total_charge_e']-c['core_total_charge_e']==1
 boundary=load(verify(s['boundary_mapping']));source=verify(s['source']);verify(boundary['forcefield']);assert len(s['core_atoms'])==54 and len(s['environment_atoms'])==9087
 assert all(x['xyz_A']==y['xyz_A'] and (i==0 or x['element']==y['element']) for i,(x,y) in enumerate(zip(s['core_atoms'],c['core_atoms'])))
 env=s['environment_atoms'];ids={x['id']:i for i,x in enumerate(env)};core=np.array([x['xyz_A'] for x in s['core_atoms']]);physical={x['id']:x for x in s['physical_atoms']};source_rows=[]
 for line in source.read_text().splitlines():
  if line.startswith(('ATOM  ','HETATM')):
   source_rows.append({'id':f'{line[21]}/{int(line[22:26])}/{line[26]}/{line[12:16].strip()}','resname':line[17:20].strip(),'xyz_A':[float(line[30:38]),float(line[38:46]),float(line[46:54])],'element':line[76:78].strip()})
 # Selection uses structural chemistry only: closest complete exterior hydroxyl to a QM heavy atom.
 candidates=[]
 for x in source_rows:
  if x['resname'] not in ('THR','SER','TYR'):continue
  name=x['id'].split('/')[-1];spec={'THR':('OG1','HG1','CB'),'SER':('OG','HG','CB'),'TYR':('OH','HH','CZ')}[x['resname']]
  if name!=spec[0]:continue
  prefix=x['id'].rsplit('/',1)[0]+'/';keys=[prefix+n for n in spec]
  if not all(k in ids for k in keys):continue
  heavy=core[[i for i,x in enumerate(s['core_atoms']) if x['element']!='H']];dist=float(np.linalg.norm(heavy-env[ids[keys[0]]]['xyz_A'],axis=1).min());candidates.append({'distance_to_core_heavy_A':dist,'resname':x['resname'],'oxygen_id':keys[0],'hydrogen_id':keys[1],'axis_carbon_id':keys[2]})
 candidates.sort(key=lambda x:(x['distance_to_core_heavy_A'],x['oxygen_id']));sel=candidates[0];assert sel['oxygen_id']=='A/159/ /OG1'
 oi,hi,ci=[ids[sel[k]] for k in ('oxygen_id','hydrogen_id','axis_carbon_id')];o,h,carbon=[np.array(env[i]['xyz_A']) for i in (oi,hi,ci)];axis=(o-carbon)/np.linalg.norm(o-carbon);hp=o+rotate(h-o,axis,math.radians(10));b=copy.deepcopy(env);b[hi]['xyz_A']=hp.tolist()
 # Only one H rotates; source bond length and angle preserved (no hydrogen re-preparation).
 assert abs(np.linalg.norm(hp-o)-np.linalg.norm(h-o))<1e-12
 assert abs(np.dot(hp-o,axis)-np.dot(h-o,axis))<1e-12
 assert [x['charge_e'] for x in b]==[x['charge_e'] for x in env]
 heavy_env=np.array([x['xyz_A'] for i,x in enumerate(env) if x['element']!='H' and i!=oi]);nearest_before=float(np.linalg.norm(heavy_env-h,axis=1).min());nearest_after=float(np.linalg.norm(heavy_env-hp,axis=1).min());assert nearest_after>=1.2,'new nonbonded H-heavy distance below declared 1.2 A'
 # Explicit source mapping; cap chain rule uses physical CB/CA, no independent cap motion.
 maps=[]
 for i,x in enumerate(s['core_atoms']):
  row={'qm_index':i,'element':x['element'],'kind':x['kind'],'id':x['id'],'xyz_A':x['xyz_A']}
  if x['kind']=='protein_source':row.update(source_id=x['id'],source_index=x['source_index'])
  elif x['kind']=='cap':
   frag=x['id'].split('/',1)[1];ledger=next(l for l in boundary['ledgers'] if l['fragment_id']==frag);retained=next(t for t in s['core_atoms'] if t.get('source_index') in ledger['source_qm_indices'] and t['id'].endswith('/CB'));prefix=retained['id'].rsplit('/',1)[0]+'/';pa=np.array(physical[prefix+'CB']['xyz_A']);pb=np.array(physical[prefix+'CA']['xyz_A']);ja,jb=cap_jacobians(pa,pb,1.09);expected=pa+1.09*(pb-pa)/np.linalg.norm(pb-pa)
   row.update(retained_source_id=prefix+'CB',omitted_source_id=prefix+'CA',retained_xyz_A=pa.tolist(),omitted_xyz_A=pb.tolist(),length_A=1.09,jacobian_retained=ja.tolist(),jacobian_omitted=jb.tolist(),serialized_origin_roundoff_A=float(np.linalg.norm(expected-x['xyz_A'])),mapping='r_link = r_CB + 1.09 (r_CA-r_CB)/norm(r_CA-r_CB)')
   assert row['serialized_origin_roundoff_A']<1e-6
  else:
   matches=[r for r in source_rows if r['resname'] in ('PQQ','CA','LA') and np.linalg.norm(np.array(r['xyz_A'])-x['xyz_A'])<.002 and (r['element']==x['element'] or x['kind']=='metal')]
   row['source_matches']=matches;row['mapping_status']='source_coordinate_match' if len(matches)==1 else 'carver_local_hydrogen' if x['element']=='H' and not matches else 'explicitly_unresolved'
  maps.append(row)
 write(out/'core_mapping.json',maps);write(out/'boundary_mapping.json',boundary);write(out/'environment_A.json',env);write(out/'environment_B.json',b)
 env_records={}
 for label,atoms in [('A',env),('B',b)]:
  pc=out/f'environment_{label}.pc';pc.write_text(str(len(atoms))+'\n'+''.join(' '.join(format(v,'.17g') for v in (x['charge_e'],*x['xyz_A']))+'\n' for x in atoms));env_records[label]={'pointcharges':rec(pc),'atoms':rec(out/f'environment_{label}.json'),'charge_e':sum(x['charge_e'] for x in atoms),'count':len(atoms)}
 assert rec(out/'environment_A.pc')['sha256']==tasks['La']['pointcharges']['sha256']
 perturb={**sel,'pointcharge_index':hi,'axis_carbon_index':ci,'oxygen_index':oi,'axis_xyz_A':carbon.tolist(),'axis_unit':axis.tolist(),'pivot_xyz_A':o.tolist(),'angle_degrees':10.,'hydrogen_A_A':h.tolist(),'hydrogen_B_A':hp.tolist(),'hydrogen_displacement_A':float(np.linalg.norm(hp-h)),'OH_length_A':float(np.linalg.norm(h-o)),'nearest_nonbonded_environment_heavy_before_A':nearest_before,'nearest_nonbonded_environment_heavy_after_A':nearest_after,'selection_rule':'nearest complete exterior Thr/Ser/Tyr hydroxyl oxygen to any QM heavy atom; lexicographic identity tie-break; +10 degrees right-hand around carbon-to-oxygen axis','candidate_ranking':candidates,'interpretation':'modeled orientation perturbation, not equilibrium sample or experimental observation'}
 result={'protocol_id':'nikasha_fixed_core_electronic_embedding_response_preparation_v1','source_state':tasks['La']['source_state'],'source_states':{z:t['source_state'] for z,t in tasks.items()},'source':s['source'],'source_manifest':rec(mpath),'assembly':s['assembly'],'microstate':s['microstate'],'explicit_waters':s['explicit_waters'],'endpoints':{z:{'xyz':t['xyz'],'charge':t['charge'],'multiplicity':1,'archived_embedded_output_path':t['output_path'],'archived_isolated_task':t['vacuum_task']} for z,t in tasks.items()},'environments':env_records,'core_mapping':rec(out/'core_mapping.json'),'boundary_mapping':rec(out/'boundary_mapping.json'),'source_boundary_mapping':s['boundary_mapping'],'perturbation':perturb,'orca':m['orca'],'execution_policy':m['execution_policy'],'implementation':rec(__file__),'molecular_evaluations':0,'limitations':['Electronic embedding component only: full hybrid energy unavailable without cross metal/PQQ repulsion/dispersion parameters and defined covalent boundary mechanics.','Source H coordinates are inherited unchanged, including 1.18537 A Thr159 OH; no source H normalization or chemistry repair here.','Selected catalytic protein chain A plus its coordinate-matched PQQ/metal (source heterogen chain B); other chains excluded from the MM protein environment.','Archived core_xyz precision is 1e-6 A; physical cap reconstruction has sub-microangstrom rounding mismatch, preserved and recorded.']}
 for z in tasks:verify(result['endpoints'][z]['xyz'])
 write(out/'INPUTS.json',result);print(json.dumps({'inputs':str(out/'INPUTS.json'),'perturbation':{k:v for k,v in perturb.items() if k!='candidate_ranking'},'core_map_statuses':[r.get('mapping_status') for r in maps if 'mapping_status' in r]},indent=2))
if __name__=='__main__':main()
