"""Read final printed populations; no electronic calculations or affinity inference."""
import argparse,json,re,hashlib
from pathlib import Path

def pin(p):
 p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def read_population(s,title,n):
 section=s.rsplit(title,1)[1];rows=[]
 for line in section.splitlines():
  m=re.match(r'^\s*(\d+)\s+([A-Z][a-z]?)\s*:\s*([-+\d.]+)\s*$',line)
  if m:rows.append((int(m[1]),m[2],float(m[3])))
  elif rows:break
 if [r[0] for r in rows]!=list(range(n)):raise ValueError('Incomplete/duplicate population atoms')
 return rows

def audit(path):
 m=json.loads(Path(path).read_text());out=[]
 for t in m['tasks']:
  p=Path(t['output_path']);s=p.read_text()
  if 'ORCA TERMINATED NORMALLY' not in s or 'SCF CONVERGED AFTER' not in s:raise ValueError('Unfinished electronic output')
  mp=Path(t['core_mapping']['path'])
  if pin(mp)['sha256']!=t['core_mapping']['sha256']:raise ValueError('Mapping hash differs')
  atoms=sorted(json.loads(mp.read_text()),key=lambda x:x['qm_index']);n=len(atoms)
  orbital_section=s.rsplit('ORBITAL ENERGIES',1)[1].split('MULLIKEN POPULATION ANALYSIS',1)[0]
  orbitals=[tuple(map(float,x)) for x in re.findall(r'^\s*(\d+)\s+(\d+\.\d+)\s+([-\d.]+)\s+([-\d.]+)\s*$',orbital_section,re.M)]
  occupied=[x for x in orbitals if x[1]>0];virtual=[x for x in orbitals if x[1]==0]
  if not occupied or not virtual:raise ValueError('Missing frontier orbital data')
  populations={}
  for scheme in ['MULLIKEN','LOEWDIN']:
   rows=read_population(s,scheme+' ATOMIC CHARGES',n)
   closure=sum(x[2] for x in rows)-t['charge']
   if abs(closure)>n*0.00000051:raise ValueError('Charge closure exceeds printed rounding bound')
   groups={}
   for a,(_,el,q) in zip(atoms,rows):
    expected=t['metal'] if a['kind']=='metal' else a['element']
    if el!=expected:raise ValueError('Population/mapping element mismatch')
    group='metal' if a['kind']=='metal' else ('caps' if a['kind']=='cap' else '/'.join(a['source_id'].split('/')[:3]))
    groups[group]=groups.get(group,0)+q
   populations[scheme]=dict(groups_e=groups,sum_e=sum(x[2] for x in rows),charge_closure_e=closure,atom_charges_e=[x[2] for x in rows])
  homo=max(occupied,key=lambda x:x[2]);lumo=min(virtual,key=lambda x:x[2])
  out.append(dict(task=t['task_id'],output=pin(p),mapping=pin(mp),atom_count=n,printed_electron_occupancy=sum(x[1] for x in occupied),homo_Eh=homo[2],lumo_Eh=lumo[2],gap_Eh=lumo[2]-homo[2],populations=populations))
 return dict(manifest=pin(path),implementation=pin(__file__),rows=out,interpretation='Descriptive final-state populations only; not oxidation states, SCF stability, affinity, or evidence about repaired inputs')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--output',required=True);a=p.parse_args();r=audit(a.manifest)
 with Path(a.output).open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 for v in r['rows']:print(v['task'],'gap_Eh',v['gap_Eh'],'electron_occupancy',v['printed_electron_occupancy'],'metal_charges',{k:x['groups_e']['metal'] for k,x in v['populations'].items()})
