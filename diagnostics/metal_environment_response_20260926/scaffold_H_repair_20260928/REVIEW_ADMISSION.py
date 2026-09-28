"""Separate chemical review of recovered real geometry, no optimizer rerun."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
base=ROOT/'workspaces/metal_environment_response_20260926/scaffold_H_repair_v1';rows=[]
for source in ['Hans8DQ2','Hans8FNR','Mex8FNS']:
 p=base/source;r=read_json(p/'RECOVERED.json');c=read_json(verify(r['inputs']));ex=read_json(verify(c['export']));a={v['id']:v for v in read_json(verify(ex['atoms']))};excluded=[];unresolved=[]
 for inv in r['inversions']:
  h=[sid for sid in inv['neighbors'] if a[sid]['element']=='H'];center=a[inv['id']]
  if center['element']=='C' and len(h)>=2:excluded.append(dict(inv,review='nonstereogenic carbon with chemically equivalent protium neighbors; label-volume change retained'))
  else:unresolved.append(inv)
 checks={k:v for k,v in r['checks'].items() if k not in ('optimizer_success','no_inversions')};checks['no_possible_stereocenter_volume_inversion']=not unresolved;checks['saved_force_stationarity']=r['max_freeH_gradient_kcal_mol_A']<=.1
 row={'source_id':source,'admitted':all(checks.values()),'status':'chemically_reviewed_geometry_admitted' if all(checks.values()) else 'not_admitted','recovered_record':record(p/'RECOVERED.json'),'inputs':r['inputs'],'parent_export':c['export'],'saved_coordinates':r['saved_coordinates'],'saved_forces':r['saved_forces'],'free_H_ids':r['free_H_ids'],'checks':checks,'nonstereogenic_label_volume_changes':excluded,'unresolved_inversions':unresolved,'optimizer_status':None,'original_admission_status':r['status'],'policy':'Separate root-authorized chemistry review; original optimizer-status null/generic-volume failures retained. Saved gradient substitutes for unavailable convergence receipt as geometric stationarity evidence only. No molecular score used.','implementation':record(__file__),'new_energy_force_calls':0};write_new(p/'REVIEWED_ADMISSION.json',row);rows.append({'source_id':source,'admitted':row['admitted'],'admission':record(p/'REVIEWED_ADMISSION.json')})
write_new(base/'REVIEWED_ADMISSION.json',{'sources':rows,'new_energy_force_calls':0});print(rows)
