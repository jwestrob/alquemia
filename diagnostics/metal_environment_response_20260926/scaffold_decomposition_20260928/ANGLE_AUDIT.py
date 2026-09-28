"""Source native angle strain from exact XML; no Context, no new geometry."""
import json,sys
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import unit
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
p=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2';ex=read_json(p/'export_v1/EXPORT.json');atoms=read_json(verify(ex['atoms']));s=mm.XmlSerializer.deserialize(verify(ex['system']).read_text());f=next(v for v in s.getForces() if isinstance(v,mm.HarmonicAngleForce));x=np.array([a['xyz_A'] for a in atoms]);rows=[]
for i in range(f.getNumAngles()):
 a,b,c,t,k=f.getAngleParameters(i);u=x[int(a)]-x[int(b)];v=x[int(c)]-x[int(b)];theta=np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1));eq=t.value_in_unit(unit.radian);energy=.5*k.value_in_unit(unit.kilojoule_per_mole/unit.radian**2)*(theta-eq)**2/4.184;rows.append({'indices':[int(a),int(b),int(c)],'source_ids':[atoms[int(j)]['id'] for j in (a,b,c)],'elements':[atoms[int(j)]['element'] for j in (a,b,c)],'angle_degrees':float(np.degrees(theta)),'equilibrium_degrees':float(np.degrees(eq)),'energy_kcal_mol':energy})
write_new(ROOT/'workspaces/metal_environment_response_20260926/scaffold_decomposition_v1/ANGLE_AUDIT.json',{'implementation':record(__file__),'export':record(p/'export_v1/EXPORT.json'),'total_kcal_mol':sum(r['energy_kcal_mol'] for r in rows),'hydrogen_containing_kcal_mol':sum(r['energy_kcal_mol'] for r in rows if 'H' in r['elements']),'all_heavy_kcal_mol':sum(r['energy_kcal_mol'] for r in rows if 'H' not in r['elements']),'largest':sorted(rows,key=lambda r:r['energy_kcal_mol'],reverse=True)[:20],'context_evaluations':0,'new_geometries':0})
print(json.dumps(sorted(rows,key=lambda r:r['energy_kcal_mol'],reverse=True)[:5],indent=2))
