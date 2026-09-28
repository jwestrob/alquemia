"""Locate native H-angle defects in actual QM/exterior and original source."""
from pathlib import Path
import sys,json
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
p=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2';ex=read_json(p/'export_v1/EXPORT.json');atoms=read_json(verify(ex['atoms']));source=read_json(p/'INPUTS.json');qset=set(read_json(verify(source['boundary_mapping']))['all_selected_QM_source_ids']);s=mm.XmlSerializer.deserialize(verify(ex['system']).read_text());f=next(v for v in s.getForces() if isinstance(v,mm.HarmonicAngleForce));rows=[]
def angle(points):
 u=points[0]-points[1];v=points[2]-points[1];return float(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1)))
for i in range(f.getNumAngles()):
 a,b,c,eq,k=f.getAngleParameters(i);ids=[int(a),int(b),int(c)];aa=[atoms[j] for j in ids];h=any(v['element']=='H' for v in aa);n=sum(v['id'] in qset for v in aa);theta=angle(np.array([v['xyz_A'] for v in aa]));old=angle(np.array([v['source_xyz_A'] for v in aa]));target=eq.value_in_unit(unit.radian);e=.5*k.value_in_unit(unit.kilojoule_per_mole/unit.radian**2)*(theta-target)**2/4.184
 rows.append({'ids':[v['id'] for v in aa],'H_containing':h,'region':'QM' if n==3 else 'MM' if n==0 else 'cross','current_deg':np.degrees(theta),'original_deg':np.degrees(old),'target_deg':np.degrees(target),'energy_kcal_mol':e})
by={reg:{'H_angle_energy_kcal_mol':sum(r['energy_kcal_mol'] for r in rows if r['region']==reg and r['H_containing']),'worst_H_angles':sorted([r for r in rows if r['region']==reg and r['H_containing']],key=lambda r:r['energy_kcal_mol'],reverse=True)[:8]} for reg in ('QM','cross','MM')}
r={'inputs':record(p/'INPUTS.json'),'implementation':record(__file__),'by_region':by,'Asp85_H_angles':[r for r in rows if r['H_containing'] and all(x.startswith('A/85/') for x in r['ids'])],'protonation_receipt':source['source_pins']['protonation_receipt'],'protein_H_source_angle_change_max_deg':max(abs(r['current_deg']-r['original_deg']) for r in rows if r['H_containing'] and not any('/33' in x and False for x in r['ids'])),'molecular_context_calls':0,'modified_coordinates':False};write_new(ROOT/'workspaces/metal_environment_response_20260926/scaffold_decomposition_v1/H_ORIGIN_AUDIT.json',r);print(json.dumps(by['QM'],indent=2));print('Asp85',json.dumps(r['Asp85_H_angles'],indent=2))
