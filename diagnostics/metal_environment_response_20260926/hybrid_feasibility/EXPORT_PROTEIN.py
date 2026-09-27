"""Re-export the System already created by unchanged source_protein; no energies."""
import json,sys,time,functools
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
import affordable_state as source
from affordable_common import read_json,verify,record,write_new
from openmm import XmlSerializer,NonbondedForce,unit
c=read_json(ROOT/'workspaces/metal_environment_response_20260926/preparation/scout_v3/INPUTS.json')
out=ROOT/'workspaces/metal_environment_response_20260926/hybrid_preparation_v1/protein_v2'
out.mkdir(exist_ok=False)
start=time.monotonic();captured=[];original=source.app.ForceField.createSystem
# Observe and retain the return value; do not alter arguments, topology or System.
@functools.wraps(original)
def capture(self,*args,**kwargs):
 system=original(self,*args,**kwargs);captured.append(system);return system
source.app.ForceField.createSystem=capture
try:atoms,pos,charges,physical,excluded=source.source_protein(verify(c['source']))
finally:source.app.ForceField.createSystem=original
assert len(captured)==1
system=captured[0];nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
state=read_json(verify(c['source_state']));old={a['source_index']:a for a in state['physical_atoms'] if a.get('kind')=='protein_source'}
# Depending on source schema, protein physical atoms retain source_index rather than kind.
if not old:old={a['source_index']:a for a in state['physical_atoms'] if 'source_index' in a}
if len(old)!=len(atoms):raise ValueError('source protein atom inventory differs')
maxdr=0.
for i,a in enumerate(physical):
 if a['id']!=old[i]['id']:raise ValueError('source order/name differs')
 dr=np.max(abs(np.asarray(a['xyz_A'])-old[i]['xyz_A']));maxdr=max(maxdr,float(dr))
 if dr>1e-12:raise ValueError('source coordinate differs')
boundary=read_json(verify(c['boundary_mapping']));adjusted=charges.copy();removed=set()
for l in boundary['ledgers']:
 removed.update(l['removed_source_indices'])
 for i in l['recipients']:adjusted[i]+=l['each_increment_e']
env=read_json(verify(c['environments']['A']['atoms']));expect=[i for i in range(len(atoms)) if i not in removed]
if len(expect)!=len(env):raise ValueError('environment complement differs')
maxdq=0.
for i,e in zip(expect,env):
 if i!=e['source_index'] or physical[i]['id']!=e['id']:raise ValueError('environment mapping differs')
 maxdq=max(maxdq,abs(float(adjusted[i])-e['charge_e']))
 if abs(adjusted[i]-e['charge_e'])>1e-12:raise ValueError('source-derived environment charge differs')
qm={a['source_index'] for a in read_json(verify(c['core_mapping'])) if a.get('kind')=='protein_source'}
rows=[]
for i,a in enumerate(physical):
 q,s,e=nb.getParticleParameters(i)
 rows.append(dict(a,charge_e=q.value_in_unit(unit.elementary_charge),sigma_nm=s.value_in_unit(unit.nanometer),epsilon_kJ_mol=e.value_in_unit(unit.kilojoule_per_mole),region='QM_protein' if i in qm else 'MM_protein'))
exceptions=[]
for i in range(nb.getNumExceptions()):
 a,b,q,s,e=nb.getExceptionParameters(i);exceptions.append(dict(atoms=[a,b],chargeprod_e2=q.value_in_unit(unit.elementary_charge**2),sigma_nm=s.value_in_unit(unit.nanometer),epsilon_kJ_mol=e.value_in_unit(unit.kilojoule_per_mole)))
def support(ids):return 'QM' if all(i in qm for i in ids) else ('MM' if all(i not in qm for i in ids) else 'cross')
terms=[]
for f in system.getForces():
 name=f.__class__.__name__;groups=[]
 if name=='HarmonicBondForce':groups=[([int(a),int(b)],[str(k) for k in (r,k)]) for a,b,r,k in (f.getBondParameters(i) for i in range(f.getNumBonds()))]
 elif name=='HarmonicAngleForce':groups=[([int(a),int(b),int(c)],[str(v) for v in (theta,k)]) for a,b,c,theta,k in (f.getAngleParameters(i) for i in range(f.getNumAngles()))]
 elif name=='PeriodicTorsionForce':groups=[([int(a),int(b),int(c),int(d)],[str(v) for v in (period,phase,k)]) for a,b,c,d,period,phase,k in (f.getTorsionParameters(i) for i in range(f.getNumTorsions()))]
 elif name=='CMAPTorsionForce':groups=[([int(v) for v in vals[1:]],{'map_index':int(vals[0])}) for vals in (f.getTorsionParameters(i) for i in range(f.getNumTorsions()))]
 elif name not in ('NonbondedForce','CMMotionRemover'):raise ValueError('unrecognized force '+name)
 for i,(ids,params) in enumerate(groups):terms.append(dict(force=name,index=i,atoms=ids,region=support(ids),parameters=params))
(out/'system.xml').write_text(XmlSerializer.serialize(system))
write_new(out/'atoms.json',rows);write_new(out/'exceptions.json',exceptions);write_new(out/'bonded_supports.json',terms)
report=dict(status='exact_source_protein_parameters_exported',inputs=record(ROOT/'workspaces/metal_environment_response_20260926/preparation/scout_v3/INPUTS.json'),source=record(verify(c['source'])),implementation=record(__file__),source_protein_implementation=record(Path(source.__file__)),forcefield=record(source.FF),atom_count=len(atoms),charge_e=float(charges.sum()),coordinate_max_difference_A=maxdr,environment_charge_max_difference_e=maxdq,qm_protein_atom_count=len(qm),cross_bonded_counts={n:sum(t['force']==n and t['region']=='cross' for t in terms) for n in sorted({t['force'] for t in terms})},exception_count=len(exceptions),artifact_pins={p.name:record(p) for p in sorted(out.iterdir())},molecular_energy_force_calls=0,wall_seconds=time.monotonic()-start,limitation='Raw protein parent parameters; not an adopted additive boundary Hamiltonian. No PQQ/metal terms included.')
write_new(out/'RESULT.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('artifact_pins',)},indent=2))
