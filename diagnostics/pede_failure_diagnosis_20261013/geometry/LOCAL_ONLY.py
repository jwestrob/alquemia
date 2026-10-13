#!/usr/bin/env python3
"""Complete only v1's unavailable local geometry metrics using exact physical maps."""
import argparse,ast,hashlib,importlib.util,json,os,resource,socket,sys,time,traceback
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def pin(p):
 p=Path(p).resolve();h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return {'path':str(p),'sha256':h.hexdigest()}
def read(p):return json.loads(Path(p).read_text())
def verify(p):
 if pin(p['path'])!=p:raise ValueError('changed '+p['path'])
 return Path(p['path'])
def write(p,x):
 with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def freeze(out):
 out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);v1=ROOT/'workspaces/pede_failure_diagnosis_20261013/geometry_v1';m=read(v1/'manifest.json');nat=read(ROOT/'workspaces/ca_coverage_recovery_20261012/two_ca_scoring_c5q88_v1/NATIVE_COLLECTION.json');ns={c['source_case_id']:c for c in nat['cases']}
 proposals={c['case_id']:ns[c['case_id']]['proposals'] for c in m['cases'] if c['family']=='new_2Ca'};deps={};todo=['mace_site_kinematics']
 while todo:
  n=todo.pop();p=ROOT/'scripts'/(n+'.py')
  if n in deps or not p.exists():continue
  deps[n]=pin(p)
  for a in ast.walk(ast.parse(p.read_text())):
   if isinstance(a,ast.Import):todo.extend(b.name.split('.')[0] for b in a.names)
   elif isinstance(a,ast.ImportFrom) and a.module:todo.append(a.module.split('.')[0])
 cfg={'driver':pin(__file__),'previous_manifest':pin(v1/'manifest.json'),'previous_result':pin(v1/'execution/RESULT.json'),'previous_receipt':pin(v1/'execution/EXECUTION.json'),'helper':m['driver'],'plan':pin(HERE/'LOCAL_RECOVERY.md'),'code':deps,'proposals':proposals,'output':str(out/'execution'),'scope':'only19 unavailable saved local geometries; completed raw/confidence/alignment copied unchanged','numerical_threads':1,'new_energy_calls':0}
 for ps in proposals.values():
  for p in ps.values():verify(p)
 write(out/'manifest.json',cfg);print(json.dumps(pin(out/'manifest.json')))
def execute(path):
 if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURM_JOB_NODELIST'):raise RuntimeError('real Slurm required before imports')
 for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
 cfg=read(path);m=read(verify(cfg['previous_manifest']));old=read(verify(cfg['previous_result']))
 for p in [cfg['driver'],cfg['helper'],cfg['plan'],cfg['previous_receipt'],*cfg['code'].values(),*m['inputs'],*m['code'].values(),*m['input_collections'],*sum([list(v.values()) for v in cfg['proposals'].values()],[])]:verify(p)
 out=Path(cfg['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();receipt={'manifest':pin(path),'job_id':os.environ['SLURM_JOB_ID'],'hostname':socket.gethostname(),'status':'running','new_energy_calls':0,'new_folds':0,'numerical_threads':1};write(out/'START.json',receipt)
 try:
  import numpy as np,gemmi
  sys.path.insert(0,str(ROOT/'scripts'));from mace_site_kinematics import Kinematics
  import pqq_microstates as pqq,coordination_policy as policy
  spec=importlib.util.spec_from_file_location('frozen_geometry_metrics',cfg['helper']['path']);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b);b.np=np;b.gemmi=gemmi;b.pqq=pqq;b.policy=policy;reader=b.xyz_read
  rows=[]
  for c,r0 in zip(m['cases'],old['rows']):
   assert c['case_id']==r0['case_id'];r=dict(r0);r['prior_local_status']=r0['status']
   if not c['mapping']:rows.append(r);continue
   try:
    data=read(c['mapping']['path'])['context'];kin=Kinematics(data);cache={};replay={};meta=data['source_atom_metadata']
    for name,rec in c['geometries'].items():
     symbols,xyz=reader(rec);q=np.zeros(len(data['modes'])) if name=='origin' else np.asarray(read(cfg['proposals'][c['case_id']][name.removeprefix('adaptive_')]['path'])['proposal']['full_q'])
     physical,rendered,_,_=kin.evaluate(q);replay[name]={'saved_XYZ_atoms':len(xyz),'physical_atoms':len(physical),'max_saved_vs_exact_kinematic_replay_A':float(np.max(np.abs(xyz-rendered))),'physical_mode_ids':[v['id'] for v in data['modes']],'saved_full_q':q.tolist()}
     physical[kin.metal]=xyz[0]
     psymbols=[v['element'] if v else symbols[0] for v in meta]
     for link in data['core_links']:
      i,kind,*rest=link
      if kind=='source':
       j=rest[0]
       if symbols[i]!=psymbols[j]:raise ValueError('source link element mismatch')
       physical[j]=xyz[i]
     cache[rec['path']]=(psymbols,physical)
    st=gemmi.read_structure(c['raw']['path']);rawmap={};schema=None
    for ch in st[0]:
     for res in ch:
      if (ch.name,res.seqid.num,res.seqid.icode.strip())==b.selkey(c['pqq'],'')[:3]:schema=pqq.identify_naming_schema(res)
      for a in res:rawmap[(ch.name,res.seqid.num,res.seqid.icode.strip(),a.name.strip())]={'xyz':np.array([a.pos.x,a.pos.y,a.pos.z]),'element':a.element.name,'resname':res.name}
    if schema is None:raise ValueError('PQQ schema missing')
    b.xyz_read=lambda rec:cache[rec['path']]
    r['local']=b.local_case(c,schema,rawmap);r['local_map_replay']=replay;r['status']='complete';r.pop('reason',None);r.pop('traceback',None)
   except Exception as e:r.update(status='partial_or_unavailable',reason=str(e),traceback=traceback.format_exc())
   finally:b.xyz_read=reader
   rows.append(r)
  result={'manifest':pin(path),'terminal':True,'previous_raw_result':cfg['previous_result'],'canonical_control':old['canonical_control'],'rows':rows,'new_energies':0,'new_folds':0,'interpretation':old['interpretation']}
  write(out/'RESULT.json',result);receipt.update(status='complete' if all(r['status']=='complete' for r in rows) else 'partial',result=pin(out/'RESULT.json'))
 except Exception as e:receipt.update(status='failed',reason=str(e),traceback=traceback.format_exc());raise
 finally:
  u=resource.getrusage(resource.RUSAGE_SELF);receipt.update(wall_seconds=time.monotonic()-start,CPU_seconds=u.ru_utime+u.ru_stime,peak_RSS_KiB=u.ru_maxrss);write(out/'EXECUTION.json',receipt)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('operation',choices=['freeze','execute']);p.add_argument('path');a=p.parse_args();freeze(a.path) if a.operation=='freeze' else execute(a.path)
