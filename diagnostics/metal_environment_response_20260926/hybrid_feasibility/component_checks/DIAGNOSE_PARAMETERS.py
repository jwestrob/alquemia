"""Compare pinned serialized force parameters; never evaluate a potential."""
import sys
from pathlib import Path
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new
w=ROOT/'workspaces/metal_environment_response_20260926'
l=read_json(w/'hybrid_preparation_v1/ledger_v1/LEDGER.json');r=read_json(w/'component_checks_v1/results/Ca_A/receipt.json')
p=ET.parse(verify(l['native_system']));c=ET.parse(verify(r['model']))
parent={f.attrib['type']:f for f in p.find('Forces')};new={f.attrib['name']:f for f in c.find('Forces')}
terms=read_json(verify(l['artifacts']['bonded_terms.json']));results={}
for name,tag in [('HarmonicBondForce','Bonds'),('HarmonicAngleForce','Angles'),('PeriodicTorsionForce','Torsions'),('CMAPTorsionForce','Torsions')]:
 ix=sorted(t['index'] for t in terms if t['force']==name)
 expected=[parent[name].find(tag)[i].attrib for i in ix];actual=[t.attrib for t in new['retained_'+name].find(tag)]
 if expected!=actual:raise ValueError('serialized native parameter/index differs '+name)
 results[name]=dict(count=len(ix),all_native_term_indices_and_coefficients_exact=True)
 if name=='CMAPTorsionForce':
  oldmaps=ET.tostring(parent[name].find('Maps'));newmaps=ET.tostring(new['retained_'+name].find('Maps'))
  if oldmaps!=newmaps:raise ValueError('CMAP table serialization differs')
  results[name]['all_CMAP_tables_exact']=True
particles=read_json(verify(l['artifacts']['particles.json']))
oldnb=parent['NonbondedForce'];newnb=new['MM_LJ'];count=0
for i,a in enumerate(particles[:9113]):
 if a['region']=='MM':
  expected=oldnb.find('Particles')[i].attrib;actual=newnb.find('Particles')[i].attrib
  if any(expected[k]!=actual[k] for k in ('sig','eps')):raise ValueError('MM LJ per-particle differs')
  count+=1
exceptions=read_json(verify(l['artifacts']['exceptions.json']));ecount=0
for i,e in enumerate(exceptions):
 if e['region']=='MM':
  expected=oldnb.find('Exceptions')[i].attrib;actual=newnb.find('Exceptions')[i].attrib
  if any(expected[k]!=actual[k] for k in ('p1','p2','sig','eps')):raise ValueError('MM LJ exception differs')
  ecount+=1
result=dict(status='native_parameters_exact_for_retained_terms',inputs=dict(native_system=l['native_system'],scout_model=r['model']),bonded=results,MM_LJ_particle_parameters_exact=count,MM_LJ_exception_parameters_exact=ecount,energy_force_calls=0)
write_new(Path(__file__).with_name('PARAMETER_MAPPING_AUDIT.json'),result);print(result)
