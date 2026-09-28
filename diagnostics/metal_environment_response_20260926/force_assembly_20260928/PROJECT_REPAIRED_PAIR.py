"""Project matched admitted La/Dy origin gradients onto the declared Glu91 tangent."""
import argparse, sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/'scripts'))
from affordable_common import read_json,record,verify,write_new,InvalidArtifact
p=argparse.ArgumentParser()
for k in ('la','dy','motion-inputs','output'):p.add_argument('--'+k,required=True)
a=p.parse_args();c=read_json(a.motion_inputs);t=read_json(verify(c['tangents']));verify(t['origin_inputs'])
rows={};physical=None
for metal,path in [('La',a.la),('Dy',a.dy)]:
 r=read_json(path)
 if r['inputs']['sha256']!=t['origin_inputs']['sha256'] or r['metal']!=metal or r['configuration']!='A':raise InvalidArtifact('Origin/state mismatch')
 ids=[v['id'] for v in r['physical_atoms']]
 if len(ids)!=len(set(ids)) or set(ids)!=set(t['physical_source_ids']):raise InvalidArtifact('Physical inventory mismatch')
 ix={s:i for i,s in enumerate(t['physical_source_ids'])};order=[ix[s] for s in ids]
 u=np.array(t['states']['origin']['physical_unit_tangent'])[order];q=np.array(t['states']['origin']['raw_physical_A_per_radian'])[order]
 coords={v['id']:v['xyz_A'] for v in r['physical_atoms']}
 if physical is not None and coords!=physical:raise InvalidArtifact('Paired physical coordinates differ')
 physical=coords
 if not np.isclose(np.linalg.norm(u),1.,atol=1e-12):raise InvalidArtifact('Nonunit tangent')
 rows[metal]={}
 for name,g in [('electronic',r['electronic']['gradient_kcal_mol_A']),('classical',r['classical_gradient_kcal_mol_A']),('total',r['total_gradient_kcal_mol_A'])]:
  g=np.array(g)
  if g.shape!=u.shape or not np.isfinite(g).all():raise InvalidArtifact('Unavailable gradient')
  rows[metal][name]={'normalized_gradient_kcal_mol_A':float((g*u).sum()),'gradient_kcal_mol_radian':float((g*q).sum())}
difference={k:{unit:rows['Dy'][k][unit]-rows['La'][k][unit] for unit in rows['La'][k]} for k in rows['La']}
write_new(a.output,dict(implementation=record(__file__),assemblies={'La':record(a.la),'Dy':record(a.dy)},motion_inputs=record(a.motion_inputs),rows=rows,Dy_minus_La=difference,force_sign='force is negative energy gradient',full_hybrid_derivative_qualified=False,solvent_consistent=False,affinity=None))
