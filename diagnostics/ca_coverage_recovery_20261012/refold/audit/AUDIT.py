#!/usr/bin/env python3
"""Finite generated-source geometry audit. Molecular imports only inside Slurm."""
from pathlib import Path
import argparse, ast, hashlib, json, os, resource, socket, sys, time, traceback
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def pin(path):
 p=Path(path).resolve();h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return {'path':str(p),'sha256':h.hexdigest()}
def read(p):return json.loads(Path(p).read_text())
def verify(p):
 if pin(p['path'])!=p:raise ValueError('changed pin: '+p['path'])
 return Path(p['path'])
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def freeze(receipt,workspace):
 r=read(receipt)
 assert r['status']=='execution_complete_not_structurally_validated' and len(r['models'])==5
 assert [x['sample'] for x in r['models']]==list(range(5))
 out=Path(workspace).resolve();out.mkdir(parents=True,exist_ok=False)
 pending=['normalize_af3_cif','detect_pqq','pqq_microstates','carve_with_pqq'];deps={}
 while pending:
  name=pending.pop();p=ROOT/'scripts'/(name+'.py')
  if name in deps or not p.is_file():continue
  deps[name]=pin(p)
  for n in ast.walk(ast.parse(p.read_text())):
   if isinstance(n,ast.Import):pending.extend(a.name.split('.')[0] for a in n.names)
   elif isinstance(n,ast.ImportFrom) and n.module:pending.append(n.module.split('.')[0])
 for model in r['models']:verify(model['model'])
 m={'protocol':'two_Ca_generated_geometry_audit_v1','execution':pin(receipt),'fold_manifest':r['manifest'],
    'target':r['target'],'models':r['models'],'driver':pin(__file__),'plan':pin(HERE/'PLAN.md'),'sbatch':pin(HERE/'RUN.sbatch'),
    'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'),'code':deps,
    'PQQ_association_A':4.0,'association_comparison':'strictly_less_than','typed_contact_A':3.1,
    'output':str(out/'execution'),'new_hydrogens':False,'new_scores':0,'all_samples_retained':True}
 write(out/'manifest.json',m);print(json.dumps({'manifest':pin(out/'manifest.json'),'executed':False}))
def loc(c,r):return {'chain':c.name,'resnum':r.seqid.num,'icode':r.seqid.icode.strip(),'resname':r.name}
def atoms(st):
 return [(c.name,r.seqid.num,r.seqid.icode.strip(),a.name,a.element.name,str(a.altloc),a.occ,a.pos.x,a.pos.y,a.pos.z)
         for c in st[0] for r in c for a in r]
def one(item,out,gemmi,normalizer,detect,base):
 out.mkdir(parents=True,exist_ok=False)
 result={**item,'status':'unavailable','scope':'source geometry only; no protonation/context/scoring admission','new_hydrogens':False,'energy_calls':0}
 try:
  source=verify(item['model']);st=gemmi.read_structure(str(source))
  if len(st)!=1:raise ValueError('exactly one model required')
  before=atoms(st);rawpqq=detect.inspect_pqq_candidates(st[0])
  result['raw_PQQ_candidates']=[x.to_dict() for x in rawpqq]
  if len(rawpqq)!=1:raise ValueError('exactly one schema-valid PQQ required')
  normalized_pdb=out/'normalizer_intermediate.pdb'
  result['normalizer']=normalizer.normalize(source,normalized_pdb)
  norm=gemmi.read_structure(str(normalized_pdb));matches=detect.inspect_pqq_candidates(norm[0])
  if len(matches)!=1:raise ValueError('normalizer did not retain unique PQQ')
  p=rawpqq[0];expected=(p.locator.chain,p.locator.resnum,p.locator.icode)
  target=matches[0].locator
  if (target.chain,target.resnum,target.icode)!=expected or target.resname!='PQQ':raise ValueError('PQQ normalizer locator changed')
  # Transfer the validated name only; audit and output retain full raw coordinates.
  result['residue_name_alias']={'before':p.locator.to_dict(),'after':target.to_dict()}
  for c in st[0]:
   for r in c:
    if (c.name,r.seqid.num,r.seqid.icode.strip())==expected:r.name='PQQ'
  full=out/'full_source_PQQ_normalized.cif';st.make_mmcif_document().write_file(str(full))
  model=gemmi.read_structure(str(full))
  if atoms(model)!=before:raise ValueError('full-source atom identity/element/occupancy/coordinate changed')
  result['full_source']=pin(full);result['all_source_atoms_preserved_exactly']=True
  result['source_atom_count']=len(before);result['source_heavy_atom_count']=sum(a[4] not in ('H','D') for a in before)
  result['source_hydrogen_count']=sum(a[4] in ('H','D') for a in before)
  pqq=detect.inspect_pqq_candidates(model[0])[0];result['PQQ']={**pqq.to_dict(),'validated_microstate_metadata':pqq.prepared.metadata(),
      'source_atom_count':len(pqq.residue),'source_heavy_atom_count':sum(a.element.name not in ('H','D') for a in pqq.residue),
      'microstate_metadata_is_not_source_protonation':True}
  sites=base.find_metal_sites(model[0]);result['metal_sites']=[]
  if len(sites)!=2 or any(s.element!='Ca' for s in sites):raise ValueError('exactly two native Ca sites required')
  for site in sites:
   residue=next(r for c in model[0] for r in c if (c.name,r.seqid.num,r.seqid.icode.strip(),r.name)==(site.chain,site.resnum,site.icode,site.resname))
   if len(residue)!=1:raise ValueError('Ca source residue must remain monatomic')
   nearby=base._prepare_nearby_pqq(model[0],site.atom.pos,base.DEFAULT_PQQ_MICROSTATE)
   qualified,included,sulfur,untyped=base._classify_contacts(model[0],site,nearby)
   nearest=min(site.atom.pos.dist(a.pos) for a in pqq.residue if a.element.name not in ('H','D'))
   result['metal_sites'].append({'selector':{'chain':site.chain,'resnum':site.resnum,'icode':site.icode,'resname':site.resname,'atom':site.atom_name,'element':site.element},
       'label':site.label(),'xyz_A':list(site.atom.pos),'nearest_PQQ_heavy_A':nearest,'PQQ_associated':nearest<4.0,
       'PQQ_donor_distances_A':{a.name:site.atom.pos.dist(a.pos) for a in pqq.residue if a.name in pqq.prepared.direct_donors},
       'typed_CN':len(qualified),'typed_contacts':[c.metadata() for c in qualified],
       'untyped_ON_3p1A':[c.metadata() for c in untyped],'sulfur_contacts':[c.metadata() for c in sulfur]})
  result['Ca_Ca_distance_A']=sites[0].atom.pos.dist(sites[1].atom.pos)
  active=[s for s in result['metal_sites'] if s['PQQ_associated']]
  result['PQQ_associated_Ca_count']=len(active)
  if len(active)!=1:raise ValueError('unique Ca strictly within4.0A of PQQ not established')
  result['active_Ca']=active[0]['selector'];result['spectator_Ca']=next(s['selector'] for s in result['metal_sites'] if not s['PQQ_associated'])
  result['status']='unique_PQQ_Ca_source_geometry';result['spectator_context_admission']='not_evaluated'
 except Exception as e:result.update(reason=str(e),exception=traceback.format_exc())
 write(out/'AUDIT.json',result);return result

def execute(manifest):
 if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURM_JOB_NODELIST'):raise RuntimeError('real Slurm allocation required')
 for k in THREADS:os.environ[k]='1'
 m=read(manifest)
 for k in ('driver','plan','sbatch','python','execution','fold_manifest'):verify(m[k])
 for p in m['code'].values():verify(p)
 out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 receipt={'manifest':pin(manifest),'job_id':os.environ['SLURM_JOB_ID'],'hostname':socket.gethostname(),'status':'running','new_scores':0,'new_hydrogens':False,'numerical_threads':1}
 write(out/'START.json',receipt)
 try:
  sys.path.insert(0,str(ROOT/'scripts'))
  import gemmi, normalize_af3_cif as normalizer, detect_pqq as detect, carve_with_pqq as base
  rows=[one(i,out/('sample-'+str(i['sample'])),gemmi,normalizer,detect,base) for i in m['models']]
  summary={'manifest':pin(manifest),'terminal':True,'rows':rows,'denominator':5,'unique_PQQ_Ca_sources':sum(r['status']=='unique_PQQ_Ca_source_geometry' for r in rows),
           'scope':'source geometry only; no protonation/context/scoring admission','new_scores':0,'new_hydrogens':False,'all_samples_retained':True}
  write(out/'FINAL_AUDIT.json',summary);receipt.update(status='complete',result=pin(out/'FINAL_AUDIT.json'))
 except Exception as e:receipt.update(status='failed',reason=str(e),exception=traceback.format_exc());raise
 finally:
  usage=resource.getrusage(resource.RUSAGE_SELF);receipt.update(wall_seconds=time.monotonic()-started,CPU_seconds=usage.ru_utime+usage.ru_stime,peak_RSS_KiB=usage.ru_maxrss)
  write(out/'EXECUTION.json',receipt)
if __name__=='__main__':
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='command',required=True)
 a=s.add_parser('freeze');a.add_argument('--receipt',required=True);a.add_argument('--workspace',required=True)
 a=s.add_parser('execute');a.add_argument('--manifest',required=True)
 a=p.parse_args();freeze(a.receipt,a.workspace) if a.command=='freeze' else execute(a.manifest)
