"""Build a field pair from existing, matched 4MAE core and parent charges only."""
from pathlib import Path
import argparse,copy,json,math,re,sys
import numpy as np

def main():
 p=argparse.ArgumentParser();p.add_argument('--repository',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.repository.resolve();sys.path.insert(0,str(root/'scripts'))
 from affordable_common import read_json,record,verify,write_new
 from affordable_response import cap_jacobians
 ip=root/'workspaces/scaffold_environment_20260922/inventory_v2/INVENTORY.json';inventory=read_json(ip);row=next(r for r in inventory['rows'] if r['case_id']=='4MAE');source=verify(row['source']);carve=read_json(verify(row['source_carve']));parent=read_json(verify(row['artifacts']['parent_atoms.json']));verify(row['artifacts']['parent_system.xml']);verify(inventory['forcefield'])
 assert read_json(verify(row['source_protonation']))['output']==row['source'];assert read_json(verify(row['terminal_repair']))['source']==row['source']
 byid={v['source_id']:v for v in parent};corepath=verify(carve['outputs']['La_xyz']);lines=corepath.read_text().splitlines();qxyz=np.array([[float(t) for t in line.split()[1:]] for line in lines[2:]]);assert len(qxyz)==80
 mapping=[{'qm_index':0,'id':'metal','kind':'metal','element':'La','xyz_A':qxyz[0].tolist(),'source_selector':carve['selected_site']}];offset=1;removed=set();deltas={};ledger=[]
 for frag in carve['qm_fragments']:
  name=frag['id'];n=frag['atom_count'];atoms=frag['atom_records'];assert len(atoms)==n
  if frag['kind']=='fixed_core_pqq':
   for j,atom in enumerate(atoms):
    assert np.max(abs(qxyz[offset+j]-atom['xyz_A']))<1e-6
    mapping.append({'qm_index':offset+j,'id':'pqq/'+atom['name'],'kind':'pqq','element':atom['element'],'xyz_A':qxyz[offset+j].tolist(),'origin':atom['origin'],'source_selector':name+'/'+atom['name'] if atom['origin']=='source_heavy_atom' else None})
  else:
   match=re.fullmatch(r'([^:]+):([A-Z]+)(\d+)',name);chain,resname,num=match.groups();prefix=f'{chain}/{num}//';localsource=[]
   for j,atom in enumerate(atoms):
    pos=qxyz[offset+j];assert np.max(abs(pos-atom['xyz_A']))<1e-6
    m={'qm_index':offset+j,'id':prefix+atom['name'],'element':atom['element'],'xyz_A':pos.tolist()}
    if atom['origin']=='generated_Cbeta_link_cap':
     ra=np.array(byid[prefix+'CB']['xyz_A']);rb=np.array(byid[prefix+'CA']['xyz_A']);ja,jb=cap_jacobians(ra,rb,1.09);err=float(np.linalg.norm(ra+1.09*(rb-ra)/np.linalg.norm(rb-ra)-pos));assert err<1e-6
     m.update(kind='cap',id='cap/'+name,retained_source_id=prefix+'CB',omitted_source_id=prefix+'CA',retained_xyz_A=ra.tolist(),omitted_xyz_A=rb.tolist(),length_A=1.09,jacobian_retained=ja.tolist(),jacobian_omitted=jb.tolist(),serialized_origin_roundoff_A=err)
    else:
     v=byid[m['id']];assert v['element']==atom['element'] and np.max(abs(np.array(v['xyz_A'])-pos))<1e-12;localsource.append(m['id']);m.update(kind='protein_source',source_id=m['id'],parent_index=v['parent_index'])
    mapping.append(m)
   local=set(localsource)|{prefix+'CA'};assert not removed&local;removed|=local
   res=[v for v in parent if v['source_id'].startswith(prefix)];q0=sum(v['charge_e'] for v in res);retained=[v for v in res if v['source_id'] not in local];change=q0-frag['formal_charge']-sum(v['charge_e'] for v in retained);recipients=[prefix+'N',prefix+'C']
   for sid in recipients:deltas[sid]=change/2
   ledger.append({'fragment_id':name,'formal_charge':frag['formal_charge'],'removed_source_ids':sorted(local),'recipients':recipients,'each_increment_e':change/2,'original_residue_charge_e':q0})
  offset+=n
 assert offset==80
 env=[]
 for v in parent:
  if v['source_id'] in removed:continue
  env.append({'id':v['source_id'],'element':v['element'],'resname':v['resname'],'parent_index':v['parent_index'],'xyz_A':v['xyz_A'],'charge_e':v['charge_e']+deltas.get(v['source_id'],0),'source_or_archived_completion':v['source_or_archived_completion']})
 expected=sum(v['charge_e'] for v in parent)-sum(l['formal_charge'] for l in ledger);assert abs(sum(v['charge_e'] for v in env)-expected)<1e-10
 byenv={v['id']:v for v in env};heavy=qxyz[[i for i,v in enumerate(mapping) if v['element']!='H']];candidates=[]
 for v in env:
  spec={'THR':('OG1','HG1','CB'),'SER':('OG','HG','CB'),'TYR':('OH','HH','CZ')}.get(v['resname'])
  if not spec or not v['id'].endswith('/'+spec[0]):continue
  prefix=v['id'].rsplit('/',1)[0]+'/';ids=[prefix+s for s in spec]
  if all(s in byenv for s in ids):candidates.append({'distance_to_core_heavy_A':float(np.linalg.norm(heavy-v['xyz_A'],axis=1).min()),'oxygen_id':ids[0],'hydrogen_id':ids[1],'axis_carbon_id':ids[2]})
 candidates.sort(key=lambda x:(x['distance_to_core_heavy_A'],x['oxygen_id']));selected=candidates[0];o,h,c=[np.array(byenv[selected[k]]['xyz_A']) for k in ('oxygen_id','hydrogen_id','axis_carbon_id')];axis=(o-c)/np.linalg.norm(o-c);r=h-o;t=math.radians(10);hp=o+r*math.cos(t)+np.cross(axis,r)*math.sin(t)+axis*np.dot(axis,r)*(1-math.cos(t));b=copy.deepcopy(env);idx=next(i for i,v in enumerate(b) if v['id']==selected['hydrogen_id']);b[idx]['xyz_A']=hp.tolist()
 nonbonded=np.array([v['xyz_A'] for v in env if v['element']!='H' and v['id']!=selected['oxygen_id']]+heavy.tolist());mind=float(np.linalg.norm(nonbonded-hp,axis=1).min());assert mind>=1.2;assert abs(np.linalg.norm(hp-o)-np.linalg.norm(h-o))<1e-12
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'core_mapping.json',mapping);write_new(out/'boundary_mapping.json',{'source':row['source'],'forcefield':inventory['forcefield'],'parent_inventory':record(ip),'parent_atoms':row['artifacts']['parent_atoms.json'],'parent_system':row['artifacts']['parent_system.xml'],'terminal_repair':row['terminal_repair'],'ledgers':ledger,'ff_chain_charge_e':sum(v['charge_e'] for v in parent),'environment_charge_e':expected,'policy':'original sidechain charge removal plus CA and equal N/C residue-local redistribution'})
 environments={}
 for label,atoms in [('A',env),('B',b)]:
  write_new(out/f'environment_{label}.json',atoms);pc=out/f'environment_{label}.pc';pc.write_text(str(len(atoms))+'\n'+''.join(' '.join(format(q,'.17g') for q in (v['charge_e'],*v['xyz_A']))+'\n' for v in atoms));environments[label]={'pointcharges':record(pc),'atoms':record(out/f'environment_{label}.json'),'count':len(atoms),'charge_e':sum(v['charge_e'] for v in atoms)}
 endpoints={}
 for metal,charge in [('Ca',-3),('La',-2)]:
  xyz=verify(carve['outputs'][metal+'_xyz']);coords=np.array([[float(v) for v in l.split()[1:]] for l in xyz.read_text().splitlines()[2:]]);assert np.array_equal(coords,qxyz);endpoints[metal]={'xyz':record(xyz),'charge':charge,'multiplicity':1}
 result={'protocol_id':'nikasha_4MAE_canonical_dry_field_response_preparation_v1','source':row['source'],'source_carve':row['source_carve'],'source_protonation':row['source_protonation'],'source_parent_inventory':record(ip),'source_state':record(out/'boundary_mapping.json'),'core_mapping':record(out/'core_mapping.json'),'boundary_mapping':record(out/'boundary_mapping.json'),'endpoints':endpoints,'environments':environments,'perturbation':dict(selected,pointcharge_index=idx,axis_unit=axis.tolist(),pivot_xyz_A=o.tolist(),hydrogen_A_A=h.tolist(),hydrogen_B_A=hp.tolist(),angle_degrees=10.,OH_length_A=float(np.linalg.norm(r)),displacement_A=float(np.linalg.norm(hp-h)),minimum_nonbonded_H_heavy_after_A=mind,candidate_ranking=candidates),'assembly':'archived catalytic chain A with selected PQQ/metal; conditional dry normalization','explicit_waters':[],'excluded_source_chemistry':['all source waters','directly coordinating crystallization adduct15P603; archived exclusion without replacement'],'microstate':'canonical PQQ3minus, frozen source protein protonation','implementation':record(__file__),'plan':record(Path(__file__).with_name('PLAN.md')),'full_hybrid_status':'unsupported_cross_parameters_and_boundary_reference','molecular_evaluations':0}
 write_new(out/'INPUTS.json',result);print(json.dumps({'inputs':str(out/'INPUTS.json'),'core_atoms':80,'environment_atoms':len(env),'environment_charge_e':expected,'perturbation':{k:v for k,v in result['perturbation'].items() if k!='candidate_ranking'}},indent=2))
if __name__=='__main__':main()
