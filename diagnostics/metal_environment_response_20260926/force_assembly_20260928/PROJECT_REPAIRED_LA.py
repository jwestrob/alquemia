"""Project an admitted repaired La gradient onto the already-prepared Glu91 motion."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,InvalidArtifact
p=argparse.ArgumentParser();p.add_argument('--assembly',required=True);p.add_argument('--motion-inputs',required=True);p.add_argument('--old-projection',required=True);p.add_argument('--output',required=True);a=p.parse_args()
r=read_json(a.assembly);c=read_json(a.motion_inputs);t=read_json(verify(c['tangents']));verify(t['origin_inputs'])
if r['inputs']['sha256']!=t['origin_inputs']['sha256'] or r['metal']!='La' or r['configuration']!='A':raise InvalidArtifact('Different physical origin/state')
ids=[v['id'] for v in r['physical_atoms']]
if len(ids)!=len(set(ids)) or set(ids)!=set(t['physical_source_ids']):raise InvalidArtifact('Different physical inventory')
ix={s:i for i,s in enumerate(t['physical_source_ids'])};order=[ix[s] for s in ids];u=np.array(t['states']['origin']['physical_unit_tangent'])[order];q=np.array(t['states']['origin']['raw_physical_A_per_radian'])[order]
if not np.isclose(np.linalg.norm(u),1.,atol=1e-12):raise InvalidArtifact('Nonunit tangent')
components={}
for name,g in [('electronic',r['electronic']['gradient_kcal_mol_A']),('classical',r['classical_gradient_kcal_mol_A']),('total',r['total_gradient_kcal_mol_A'])]:
 g=np.array(g)
 if g.shape!=u.shape or not np.isfinite(g).all():raise InvalidArtifact('Missing/invalid gradient')
 components[name]=dict(normalized_gradient_kcal_mol_A=float((g*u).sum()),gradient_kcal_mol_radian=float((g*q).sum()))
old=read_json(a.old_projection)['rows']['A']['modes']['Glu91_CB_CG']['components']['total']
write_new(a.output,dict(assembly=record(a.assembly),motion_inputs=record(a.motion_inputs),old_projection=record(a.old_projection),implementation=record(__file__),components=components,old_La_normalized_gradient_kcal_mol_A=old['gradient_load_kcal_mol_A']['La'],old_La_gradient_kcal_mol_radian=old['raw_generalized_gradient']['La'],Dy_minus_La_gradient=None,affinity=None,full_hybrid_derivative_qualified=False,solvent_consistent=False,interpretation='One-metal response to preparation repair, not differential selectivity or a stationary minimum'))
