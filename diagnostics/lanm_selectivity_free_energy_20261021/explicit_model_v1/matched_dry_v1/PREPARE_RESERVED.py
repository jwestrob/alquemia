#!/usr/bin/env python3
"""Reserve24 fixed validation cells; geometry-only and no result reads."""
import argparse,copy,hashlib,importlib.util,json,math,os,shutil,time
from pathlib import Path
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();b=a.root.resolve();out=a.output.resolve();tick=time.monotonic();out.mkdir(parents=True,exist_ok=False);(out/'xyz').mkdir()
 assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
 source=b/'diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/REFERENCE.py';assert pin(source)['sha256']=='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129';sp=importlib.util.spec_from_file_location('reference_saved_geometry',source);ref=importlib.util.module_from_spec(sp);sp.loader.exec_module(ref)
 mixedp=b/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/prepared_v2/PANEL.json';assert pin(mixedp)['sha256']=='1191e3a00e3c02aa861faf83aae22eb250c94064e7ac780bf2fbebb00eb6c21a';mixed=json.loads(mixedp.read_text())
 proteinp=b/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/peptide_fragments_v3/PANEL.json';assert pin(proteinp)['sha256']=='44bc0eef9e0686c8d1a52c4b3184254f2ff7401d9d5f89aed2f0f6d7eac56b22';protein=json.loads(proteinp.read_text());assert protein['status']=='prepared_no_scoring_release' and protein['energy_calls']==protein['force_calls']==0
 parentp=b/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/bound_v1/NATIVE_PANEL.json';parent=ref.load_panel(parentp);rows=[];gates=[]
 def axis_from(radial):
  u=ref.unit(radial);axis=ref.cross(u,[0,0,1])
  if ref.norm(axis)<.1:axis=ref.cross(u,[1,0,0])
  return ref.unit(axis)
 def rotate(x,center,axis,phi):
  v=ref.sub(x,center);c=math.cos(phi);s=math.sin(phi);cross=ref.cross(axis,v);projection=ref.dot(axis,v)
  return [center[k]+v[k]*c+cross[k]*s+axis[k]*projection*(1-c) for k in range(3)]
 def contacts(g,groups):
  values=[]
  for k,left in enumerate(groups):
   for right in groups[k+1:]:
    for i in left:
     for j in right:
      dist=ref.norm(ref.sub(g['positions_A'][i],g['positions_A'][j]));limit=1.2 if 'H' in (g['symbols'][i],g['symbols'][j]) else 1.8;values.append((dist,limit,i,j))
  return {'passed':all(d>=l for d,l,i,j in values),'minimum_margin_A':min(d-l for d,l,i,j in values),'closest':list(min(values))}
 def emit(g,old,transforms,contact_groups):
  oldpos=old['positions_A'];changed={i for t in transforms for i in t['indices']};xp=ref.xyz_write(out/'xyz'/(g['geometry_id']+'.xyz'),g['symbols'],g['positions_A'],g['geometry_id']+' reserved; unopened validation');ss,xx=ref.xyz_read(xp['path']);g['positions_A']=xx;g['xyz']=xp
  assert ss==g['symbols'] and all(xx[i]==oldpos[i] for i in range(g['atoms']) if i not in changed)
  maxerr=0.
  for t in transforms:
   ids=t['indices'];maxerr=max(maxerr,max((abs(ref.norm(ref.sub(xx[i],xx[j]))-ref.norm(ref.sub(oldpos[i],oldpos[j]))) for k,i in enumerate(ids) for j in ids[k+1:]),default=0))
  assert maxerr<1e-9
  check=contacts(g,contact_groups);g['validation_status']='reserved' if check['passed'] else 'unavailable_contact_gate';g['parent_geometry']=old['xyz'];g['rotation_construction']=transforms;rows.append(g);gates.append({'id':g['geometry_id'],'max_internal_distance_error_A':maxerr,**check})
 for label,angle in [('minus',-15.),('plus',15.)]:
  for metal in ('La','Dy'):
   old=next(g for g in mixed['reused_controls'] if g['metal']==metal);ref.verify(old['xyz']);g=copy.deepcopy(old);g.update(geometry_id=f'validation_aquo8_orientation_{label}_{metal}',case_id='validation_aquo8_orientation',distortion=label);trans=[];groups=[]
   for j in range(8):
    ids=[1+3*j+k for k in range(3)];o=ids[0];center=old['positions_A'][o];axis=axis_from(ref.sub(center,old['positions_A'][0]));phi=math.radians(angle)
    for i in ids[1:]:g['positions_A'][i]=rotate(old['positions_A'][i],center,axis,phi)
    trans.append({'indices':ids,'moved_indices':ids[1:],'center':center,'axis':axis,'degrees':angle});groups.append(ids)
   # Oxygen coordinates are included in rigidity checks but remain unchanged.
   emit(g,old,trans,groups);assert all(g['positions_A'][1+3*j]==old['positions_A'][1+3*j] for j in range(8))
 for family in ('acetate_mono','acetate_bi','acetamide'):
  for label,angle in [('minus',-15.),('plus',15.)]:
   for metal in ('La','Dy'):
    old=next(g for g in mixed['geometries'] if g['donor_family']==family and g['distortion']=='central' and g['metal']==metal);ref.verify(old['xyz']);g=copy.deepcopy(old);g.update(geometry_id=f'validation_hydrated_{family}_orientation_{label}_{metal}',case_id=f'validation_hydrated_{family}_orientation',distortion=label)
    ds=g['donor_indices'];center=[sum(old['positions_A'][i][k] for i in ds)/len(ds) for k in range(3)];axis=axis_from(ref.sub(center,old['positions_A'][0]));lig=g['ligand_indices']
    for i in lig:g['positions_A'][i]=rotate(old['positions_A'][i],center,axis,math.radians(angle))
    waters=[i for i in range(g['atoms']) if i!=0 and i not in lig];emit(g,old,[{'indices':lig,'moved_indices':lig,'center':center,'axis':axis,'degrees':angle}],[lig,waters])
 for old in protein['geometries']:
  ref.verify(old['xyz']);g=copy.deepcopy(old);g['parent_geometry']=old['xyz'];g.update(geometry_id='validation_'+old['geometry_id'],case_id='validation_'+old['case_id'],distortion='central',water_count=0,ligand_indices=list(range(1,old['atoms'])),validation_status='reserved')
  dest=out/'xyz'/(g['geometry_id']+'.xyz');shutil.copyfile(ref.verify(old['xyz']),dest);g['xyz']=pin(dest);assert g['xyz']['sha256']==old['xyz']['sha256'];rows.append(g);gates.append({'id':g['geometry_id'],'passed':True,'cap_v3':pin(proteinp)})
 assert len(rows)==24 and len({g['geometry_id'] for g in rows})==24
 for g in rows:
  peer=[p for p in rows if p['case_id']==g['case_id'] and p['distortion']==g['distortion'] and p['metal']!=g['metal']];assert len(peer)==1 and peer[0]['positions_A']==g['positions_A']
 specpin=pin(Path(__file__).with_name('RESERVED_VALIDATION_SPEC.md'));result={'status':'reserved_complete' if all(x['passed'] for x in gates) else 'reserved_with_contact_exclusions','implementation':pin(__file__),'specification':specpin,'source_mixed_panel':pin(mixedp),'source_protein_panel':pin(proteinp),'geometries':rows,'checks':gates,'declared_cells':24,'new_molecular_calls':0,'local_threads':1,'candidate_frozen':False,'outcomes_opened':False,'scoring_release':False,'elapsed_seconds':time.monotonic()-tick};rp=ref.write(out/'RESERVED_PANEL.json',result)
 native={k:copy.deepcopy(parent[k]) for k in ('protocol','adapter','published_coordinate_source','archived_recipe_source','orca','assets')};cells=[]
 for g in rows:
  core=46 if g['metal']=='La' else 55;ne=sum(ref.Z[s] for s in g['symbols'])-g['charge']-core;assert ne%2==0
  cells.append({'cell_id':'ecp_pbe0__'+g['geometry_id'],'recipe':'ecp_pbe0','geometry_id':g['geometry_id'],'charge':g['charge'],'multiplicity':1,'physical_multiplicity':g['physical_multiplicity'],'metal':g['metal'],'ecp_core':core,'explicit_electrons':ne,'geometry':g,'pilot':False})
 native.update(specification=specpin,source_inventory=rp,geometries=rows,free_fragments=[],cells=cells,pilot_cell_ids=[],scoring_release=False,candidate_frozen=False,interpretation='Reserved validation: no execution/readout before candidate freeze; contact-excluded slots remain unavailable')
 np=ref.write(out/'NATIVE_PANEL.json',native);vp=ref.write(out/'NATIVE_VALIDATION.json',ref.validate(np['path']));ready={'status':result['status'],'reserved_panel':rp,'native_panel':np,'native_saved_validation':vp,'declared_cells':24,'eligible_cells':sum(x['passed'] for x in gates),'new_molecular_calls':0,'scoring_release':False};ref.write(out/'READY.json',ready);print(json.dumps(ready))
if __name__=='__main__':main()
