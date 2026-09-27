"""Source-exact PQQ3- graph and GAFF2 atom types only; no charge calculation."""
import json,sys,subprocess,os,time
from pathlib import Path
import numpy as np,gemmi
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,write_new,verify
base=ROOT/'workspaces/metal_environment_response_20260926/hybrid_preparation_v1';out=base/'pqq_types_v1';out.mkdir(exist_ok=False)
ccd=base/'PQQ_CCD.cif';block=gemmi.cif.read(str(ccd)).sole_block()
atomtable=block.find('_chem_comp_atom.',['atom_id','type_symbol','charge','pdbx_aromatic_flag'])
ccd_atoms={str(r[0]):dict(element=str(r[1]),formal_charge=int(r[2]),aromatic=str(r[3])) for r in atomtable}
bonds=[(str(r[0]),str(r[1]),str(r[2]),str(r[3])) for r in block.find('_chem_comp_bond.',['atom_id_1','atom_id_2','value_order','pdbx_aromatic_flag'])]
mapping=read_json(ROOT/'workspaces/metal_environment_response_20260926/preparation/scout_v3/core_mapping.json');pqq=[a for a in mapping if a['kind']=='pqq'];heavy={a['source_matches'][0]['id'].split('/')[-1]:a for a in pqq if a['element']!='H'}
if set(heavy)!={k for k,a in ccd_atoms.items() if a['element']!='H'}:raise ValueError('CCD/source heavy atom names differ')
# Frozen current state is carboxyl-deprotonated oxidized PQQ3-: remove CCD acid
# protons only, assign formal -1 at their singly bonded oxygen; no redox change.
removed={'HOB2','HOB7','HOB9'};negative={'O2B','O7B','O9B'}
for h,o in zip(sorted(removed),sorted(negative)):
 if not any({a,b}=={h,o} and order=='SING' for a,b,order,_ in bonds):raise ValueError('CCD acid attachment differs')
rows=[]
for a in pqq:
 if a['element']!='H':name=a['source_matches'][0]['id'].split('/')[-1]
 else:
  parent=min(heavy,key=lambda n:np.linalg.norm(np.asarray(a['xyz_A'])-heavy[n]['xyz_A']))
  options=[b if x==parent else x for x,b,_,_ in bonds if parent in (x,b) and ccd_atoms[b if x==parent else x]['element']=='H' and (b if x==parent else x) not in removed]
  if len(options)!=1:raise ValueError('source H has ambiguous CCD mapping')
  if np.linalg.norm(np.asarray(a['xyz_A'])-heavy[parent]['xyz_A'])>1.25:raise ValueError('source H not bonded at declared parent')
  name=options[0]
 rows.append(dict(name=name,qm_index=a['qm_index'],element=a['element'],xyz_A=a['xyz_A'],formal_charge=-1 if name in negative else 0,ccd_aromatic_flag=ccd_atoms[name]['aromatic']))
if len(rows)!=27 or len({r['name'] for r in rows})!=27 or sum(r['formal_charge'] for r in rows)!=-3:raise ValueError('PQQ state differs')
ids={r['name']:i+1 for i,r in enumerate(rows)};keep=[(a,b,order,arom) for a,b,order,arom in bonds if a in ids and b in ids]
lines=['PQQ_oxidized_3minus_source_exact','Nikasha graph from CCD, no charge fit','',f'{len(rows):3d}{len(keep):3d}  0  0  0  0            999 V2000']
for a in rows:
 x,y,z=a['xyz_A'];lines.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a["element"]:<3s} 0  0  0  0  0  0  0  0  0  0  0  0')
for a,b,o,_ in keep:lines.append(f'{ids[a]:3d}{ids[b]:3d}{dict(SING=1,DOUB=2,TRIP=3)[o]:3d}  0  0  0  0')
lines.append('M  CHG  3'+''.join(f'{ids[n]:4d}{-1:4d}' for n in sorted(negative)));lines+=['M  END','$$$$']
(out/'graph.mol').write_text('\n'.join(lines)+'\n')
write_new(out/'source_graph.json',dict(atoms=rows,bonds=keep,ccd=record(ccd),ccd_removed_acid_protons=sorted(removed),state='oxidized_PQQ_3minus',formal_charge=-3,charge_model=None,typing_coordinate_precision_note='MDL intermediary 1e-4A only; final cross-LJ artifact uses exact source coordinates, no energy evaluates intermediary'))
exe=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/antechamber')
cmd=[str(exe),'-i','graph.mol','-fi','mdl','-o','typed.mol2','-fo','mol2','-at','gaff2','-j','1','-nc','-3','-m','1','-seq','n','-an','n','-pf','n','-s','2']
env=os.environ.copy();env['AMBERHOME']=str(exe.parent.parent)
start=time.monotonic();r=subprocess.run(cmd,cwd=out,env=env,capture_output=True,text=True)
(out/'antechamber.stdout').write_text(r.stdout);(out/'antechamber.stderr').write_text(r.stderr)
write_new(out/'EXECUTION.json',dict(command=cmd,executable=record(exe),returncode=r.returncode,elapsed_seconds=time.monotonic()-start,charge_calculation_requested=False,molecular_energy_force_calls=0))
if r.returncode:raise RuntimeError('atom typing failed; inspect retained attempt')
print(out)
