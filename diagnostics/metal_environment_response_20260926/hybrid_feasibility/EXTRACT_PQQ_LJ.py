"""Validate real typed graph and export only GAFF2 cross-LJ parameters."""
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new
base=ROOT/'workspaces/metal_environment_response_20260926/hybrid_preparation_v1/pqq_types_v1'
graph=read_json(base/'source_graph.json');text=(base/'typed.mol2').read_text();atoms=[s.split() for s in text.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0].splitlines() if s.strip()]
bonds=[s.split() for s in text.split('@<TRIPOS>BOND\n')[1].split('@<TRIPOS>')[0].splitlines() if s.strip()]
if len(atoms)!=27 or len(bonds)!=29:raise ValueError('typed inventory differs')
ids={a['name']:i+1 for i,a in enumerate(graph['atoms'])}
expected={(min(ids[a],ids[b]),max(ids[a],ids[b]),{'SING':1,'DOUB':2,'TRIP':3}[o]) for a,b,o,flag in graph['bonds']}
observed={(min(int(r[1]),int(r[2])),max(int(r[1]),int(r[2])),int(r[3])) for r in bonds}
if observed!=expected:raise ValueError('typing changed pinned bond-order graph')
ff=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/dat/leap/parm/gaff2.dat');section=ff.read_text().split('MOD4      RE',1)[1].split('\n\n',1)[0];lj={}
for line in section.splitlines():
 v=line.split()
 if len(v)>=3:lj[v[0]]=(float(v[1]),float(v[2]))
rows=[]
for i,(typed,source) in enumerate(zip(atoms,graph['atoms'])):
 expected_serialized=[float(f'{float(f"{v:.4f}"):.3f}') for v in source['xyz_A']]
 if int(typed[0])!=i+1 or not np.allclose([float(v) for v in typed[2:5]],expected_serialized,atol=1e-12,rtol=0):raise ValueError('typed atom ordering/known serialization differs')
 atype=typed[5]
 if atype not in lj:raise ValueError('LJ type absent in pinned GAFF2 table')
 r,e=lj[atype];rows.append(dict(source,gaff2_type=atype,Rmin_half_A=r,epsilon_kcal_mol=e,partial_charge_e=None))
write_new(base/'PQQ_CROSS_LJ.json',dict(status='standard_gaff2_atom_types_and_LJ_only',atoms=rows,bonds=graph['bonds'],formal_charge=-3,microstate=graph['state'],gaff2_table=record(ff),typed_intermediate=record(base/'typed.mol2'),graph=record(base/'source_graph.json'),implementation=record(__file__),partial_charge_model=None,energy_model=None,pair_mixing_rule=None,metal_parameters=None,limitations=['Typing-only mol2 has zero placeholder charge columns written by antechamber; these are not physical partial charges and are not exported or used.','Antechamber AC intermediary serializes coordinates at 1e-3A precision; returned artifact restores exact pinned source coordinates by verified atom order and graph.','No whole PQQ MM bonded energy or charge model is provided.','No cross-pair energies, mixing rule or metal LJ family chosen.']))
print('27 exact-source atoms,29 preserved CCD-derived bonds; types',sorted({a['gaff2_type'] for a in rows}))
