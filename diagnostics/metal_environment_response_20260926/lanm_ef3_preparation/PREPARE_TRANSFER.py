"""Transferred EF3 exact-source local electronic and finite-field preparation, no Context."""
from pathlib import Path
import argparse,copy,json,math,sys
import numpy as np

def main():
 p=argparse.ArgumentParser();p.add_argument('--source-id',choices=['Hans_8FNR','Mex_8FNS'],required=True);p.add_argument('--repository',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.repository.resolve();sys.path.insert(0,str(root/'scripts'))
 from affordable_common import read_json,record,verify,write_new
 from affordable_response import cap_jacobians
 import openmm as mm
 from openmm import app,unit
 auditpath=root/'diagnostics/metal_environment_response_20260926/lanm_preparation_feasibility/RESULT_v2.json';audit=read_json(auditpath);r=next(v for v in audit['sources'] if v['source_id']==a.source_id);source=next(v for v in read_json(verify(audit['manifest']))['sources'] if v['source_id']==a.source_id);data=read_json(verify(r['pins']['atoms']));atoms=data['atoms'];byid={v['id']:v for v in atoms};nparent=r['protein_atoms']+3*r['water_count'];qparent=r['protein_charge'];qcore=r['proposal']['formal_region_charge_proposal'];nqm=r['proposal']['proposed_QM_atoms_with_metal_and_caps'];assert len(atoms)==nparent and len(byid)==nparent
 ffbase=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/site-packages/openmm/app/data');ffpaths=[ffbase/'amber19/protein.ff19SB.xml',ffbase/'tip3p.xml'];pdb=app.PDBFile(str(verify(r['pins']['protonated_source'])));model=app.Modeller(pdb.topology,pdb.positions)
 def aid(v):return f'{v.residue.chain.id}/{v.residue.id}/{v.residue.insertionCode.strip()}/{v.name}'
 model.delete([v for v in model.topology.atoms() if aid(v) not in byid]);topatoms=list(model.topology.atoms());assert len(topatoms)==nparent and {aid(v) for v in topatoms}==set(byid)
 ff=app.ForceField(*(str(p) for p in ffpaths));system=ff.createSystem(model.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False);nb=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce));charges={aid(v):nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i,v in enumerate(topatoms)}
 assert abs(sum(charges.values())-qparent)<1e-10;assert not any(v.residue.name in ('LA','DY','ND','NA') for v in topatoms)
 bonds={frozenset((aid(u),aid(v))) for u,v in model.topology.bonds()};adj={sid:set() for sid in byid}
 for pair in bonds:
  u,v=tuple(pair);adj[u].add(v);adj[v].add(u)
 exported=[dict(v,charge_e=charges[v['id']]) for v in atoms]
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);ex=out/'export_v1';ex.mkdir();write_new(ex/'atoms_charges.json',exported);write_new(ex/'bonds.json',[sorted(v) for v in sorted(bonds,key=lambda b:tuple(sorted(b)))]);(ex/'protein_water_system.xml').write_text(mm.XmlSerializer.serialize(system))
 write_new(ex/'EXPORT.json',{'source':r['pins']['protonated_source'],'normalized_coordinates':r['pins']['atoms'],'forcefields':[record(p) for p in ffpaths],'atom_count':nparent,'charge_e':sum(charges.values()),'openmm_version':mm.__version__,'atoms':record(ex/'atoms_charges.json'),'bonds':record(ex/'bonds.json'),'system':record(ex/'protein_water_system.xml'),'energy_force_calls':0,'hydrogen_additions':0,'coordinate_changes':0,'source_IDs_exact':True})
 selected=[v['id'] for v in r['proposal']['source_atoms']];qmset=set(selected);assert len(selected)==r['proposal']['real_nonmetal_atoms'] and len(qmset)==len(selected)
 for item in r['proposal']['source_atoms']:assert atoms[item['physical_index']]['id']==item['id']
 cap_specs=r['proposal']['caps'];omitted={v['omitted'] for v in cap_specs};assert len(omitted)==4
 assert {tuple(sorted((u,v))) for u in qmset for v in adj[u] if v not in qmset}=={tuple(sorted((v['retained'],v['omitted']))) for v in cap_specs}
 removed=qmset|omitted;newcharges=charges.copy();ledger=[]
 for cut in cap_specs:
  retained,omit=cut['retained'],cut['omitted'];prefix=omit.rsplit('/',1)[0]+'/';resids={sid for sid in byid if sid.startswith(prefix)};local=resids&removed;remain=resids-removed;recipients=sorted(sid for sid in adj[omit] if sid in remain and byid[sid]['element']!='H');assert len(recipients)==2
  wanted=sum(charges[sid] for sid in resids);before=sum(charges[sid] for sid in remain);increment=(wanted-before)/len(recipients)
  for sid in recipients:newcharges[sid]+=increment
  before_dip=sum((charges[sid]*np.array(byid[sid]['xyz_A']) for sid in resids),np.zeros(3));after_dip=sum((newcharges[sid]*np.array(byid[sid]['xyz_A']) for sid in remain),np.zeros(3))
  assert abs(sum(newcharges[sid] for sid in remain)-wanted)<1e-12
  ledger.append({'partial_residue':prefix,'formal_QM_fragment_charge':0,'original_residue_charge_e':wanted,'removed_source_ids':sorted(local),'retained_MM_charge_before_e':before,'target_retained_MM_charge_e':wanted,'recipients':recipients,'each_increment_e':increment,'dipole_before_eA':before_dip.tolist(),'retained_MM_dipole_after_eA':after_dip.tolist(),'interpretation':'change of MM charge representation; QM contribution not included, no dipole conservation asserted'})
 # Complete loop formal charge is independently checked from template charges.
 loopids=[sid for sid in byid if r['proposal']['loop_residues'][0]<=int(sid.split('/')[1])<=r['proposal']['loop_residues'][1] and byid[sid]['kind']=='protein_source'];assert set(loopids)<=qmset;assert abs(sum(charges[sid] for sid in loopids)-(qcore-3))<1e-10
 sourcewaters={v['id'].rsplit('/',1)[0] for v in atoms if v['kind']=='retained_crystal_water'};selectedwaters={sid.rsplit('/',1)[0] for sid in qmset if byid[sid]['kind']=='retained_crystal_water'}
 assert len(selectedwaters)==({'Hans_8FNR':3,'Mex_8FNS':7}[a.source_id])
 for residue in sourcewaters:assert len([sid for sid in byid if sid.rsplit('/',1)[0]==residue])==3
 spectators=[]
 for v in r['spectators']:
  z=v['source_element'];charge={'La':3.,'Dy':3.,'Nd':3.,'Na':1.}[z];src=v['source'];spectators.append({'id':f"{src['chain']}/{src['resid']}/{src['insertion_code']}/{src['name']}",'element':z,'xyz_A':src['xyz_A'],'charge_e':charge,'site':v['site'],'kind':'spectator_formal_monopole','representation':'frozen formal charge, no FF ion parameters/LJ or polarization'})
 env=[dict(v,charge_e=newcharges[v['id']]) for v in atoms if v['id'] not in removed]+spectators
 totalcharge=qparent+3+sum(v['charge_e'] for v in spectators);qfield=totalcharge-qcore;assert abs(sum(v['charge_e'] for v in env)-qfield)<1e-10
 # Retrieve target from actual source site inventory; never use an occupancy hypothesis XYZ.
 target=source['sites']['EF3'] if isinstance(source['sites'],dict) else next(v for v in source['sites'] if v['site']=='EF3')
 if 'source' in target:target=target['source']
 physicalA={sid:np.array(v['xyz_A']) for sid,v in byid.items()};physicalB={sid:pos.copy() for sid,pos in physicalA.items()};pert=r['proposal']['perturbation'];ca,cb=[physicalA[sid] for sid in pert['axis']];axis=(cb-ca)/np.linalg.norm(cb-ca);angle=math.radians(pert['right_hand_angle_degrees'])
 expected_moving=set(pert['rotate_source_atoms']);assert expected_moving=={pert['axis'][0].rsplit('/',1)[0]+'/'+name for name in ['HB2','HB3','CG','OD1','OD2']};assert expected_moving<=qmset
 for sid in expected_moving:
  v=physicalA[sid]-cb;physicalB[sid]=cb+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+axis*np.dot(axis,v)*(1-np.cos(angle))
 for pair in bonds:
  u,v=tuple(pair)
  if u in expected_moving or v in expected_moving:assert abs(np.linalg.norm(physicalA[u]-physicalA[v])-np.linalg.norm(physicalB[u]-physicalB[v]))<1e-12
 pc=out/'environment.pc';pc.write_text(str(len(env))+'\n'+''.join(' '.join(format(q,'.17g') for q in (v['charge_e'],*v['xyz_A']))+'\n' for v in env));write_new(out/'environment_atoms.json',env);write_new(out/'boundary_mapping.json',{'policy':'residue-local heavy-neighbor charge shift v1','ledger':ledger,'all_selected_QM_source_ids':selected,'removed_MM1_CA_ids':sorted(omitted),'retained_physical_CA_coordinates':{sid:byid[sid]['xyz_A'] for sid in omitted},'spectators':spectators,'MM_field_charge_e':sum(v['charge_e'] for v in env),'parent_protein_water_charge_e':sum(charges.values()),'whole_physical_source_charge_e':totalcharge,'QM_charge_e':qcore,'full_hybrid':False})
 configurations={};maps={};Z={'H':1,'C':6,'N':7,'O':8,'S':16,'La':57,'Dy':66}
 for label,physical in [('A',physicalA),('B',physicalB)]:
  mapping=[{'qm_index':0,'id':'target_EF3','kind':'metal','xyz_A':target['xyz_A'],'source_site':target}]
  for sid in selected:mapping.append({'qm_index':len(mapping),'id':sid,'kind':'source','element':byid[sid]['element'],'xyz_A':physical[sid].tolist(),'source_id':sid,'parent_index':byid[sid]['index']})
  for cut in cap_specs:
   u,v=cut['retained'],cut['omitted'];x,y=physical[u],physical[v];length=1.09 if byid[u]['element']=='C' else 1.01;ja,jb=cap_jacobians(x,y,length);mapping.append({'qm_index':len(mapping),'id':'cap/'+u,'kind':'cap','element':'H','xyz_A':(x+length*(y-x)/np.linalg.norm(y-x)).tolist(),'retained_source_id':u,'omitted_source_id':v,'retained_xyz_A':x.tolist(),'omitted_xyz_A':y.tolist(),'length_A':length,'jacobian_retained':ja.tolist(),'jacobian_omitted':jb.tolist()})
  assert len(mapping)==nqm;write_new(out/f'core_mapping_{label}.json',mapping);maps[label]=record(out/f'core_mapping_{label}.json');endpoints={}
  for metal,mult in [('La',1),('Dy',6)]:
   elements=[metal]+[v['element'] for v in mapping[1:]];electrons=sum(Z[e] for e in elements)-qcore;assert (electrons-(mult-1))%2==0;xyz=out/f'{metal}_{label}.xyz';xyz.write_text(f'{nqm}\n{a.source_id} EF3 {metal}(III) charge={qcore} mult={mult} {label}; actual source spectators\n'+''.join(el+' '+' '.join(format(p,'.17g') for p in v['xyz_A'])+'\n' for el,v in zip(elements,mapping)));endpoints[metal]={'xyz':record(xyz),'charge':qcore,'multiplicity':mult,'all_electron_count':electrons,'oxidation_state':3,'physical_spin_2S':mult-1,'electronic_state':'LaIII closed shell' if metal=='La' else 'DyIII high-spin sextet hypothesis; no spin-orbit treatment qualified'}
  configurations[label]={'endpoints':endpoints,'pointcharges':record(pc),'core_mapping':maps[label]}
 result={'protocol_id':f'nikasha_{a.source_id}_EF3_normalized{nqm}_finite_field_preparation_v1','source_id':a.source_id,'source_state':record(auditpath),'source_pins':r['pins'],'source_manifest':audit['manifest'],'export':record(ex/'EXPORT.json'),'configurations':configurations,'core_mapping':maps,'boundary_mapping':record(out/'boundary_mapping.json'),'physical_source':r['pins']['atoms'],'environment_atoms':record(out/'environment_atoms.json'),'perturbation':dict(pert,axis_unit=axis.tolist(),pivot_xyz_A=cb.tolist(),changed_source_ids=sorted(expected_moving)),'assembly':source['assembly_limitation'],'spectator_occupancy':spectators,'waters':{'total_source':r['water_count'],'QM_selected_residues':sorted(selectedwaters),'inventory':r['pins']['water_inventory']},'charge_accounting':{'protein_water_e':qparent,'all_site_ions_e':totalcharge-qparent,'total_e':totalcharge,'QM_e':qcore,'MM_e':qfield},'implementation':record(__file__),'energy_force_calls':0,'full_hybrid_status':'unavailable_no_cross_LJ_or_bonded_total','numerical_state_qualification':'pending_native_reference','reserved_outcomes_accessed':False}
 write_new(out/'INPUTS.json',result);print(json.dumps({'inputs':str(out/'INPUTS.json'),'field_atoms':len(env),'boundary':ledger,'states':{k:{z:{f:t[f] for f in ['charge','multiplicity','all_electron_count']} for z,t in v['endpoints'].items()} for k,v in configurations.items()}},indent=2))
if __name__=='__main__':main()
