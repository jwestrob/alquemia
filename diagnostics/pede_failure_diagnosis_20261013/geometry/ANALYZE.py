#!/usr/bin/env python3
"""Finite saved PedE geometry/confidence inspection. Real Slurm required for analysis."""
from pathlib import Path
import argparse, ast, hashlib, json, os, resource, socket, sys, time, traceback
from types import SimpleNamespace
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def pin(p):
 p=Path(p).resolve();h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
 if pin(p['path'])!=p:raise ValueError('pin changed: '+p['path'])
 return Path(p['path'])
def read(p):return json.loads(Path(p).read_text())
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def confidence_pins(model):
 p=Path(model['path']);stem=p.name.removesuffix('model.cif')
 return {name:pin(p.parent/(stem+name+'.json')) if (p.parent/(stem+name+'.json')).exists() else None for name in ('confidences','summary_confidences')}
def freeze(workspace):
 out=Path(workspace).resolve();out.mkdir(parents=True,exist_ok=False)
 rec=ROOT/'workspaces/ca_coverage_recovery_20261012';oldpath=ROOT/'workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json';envpath=ROOT/'workspaces/motion_envelope_20260923/run_v2/PREPARATION.json'
 nativepath=rec/'two_ca_scoring_c5q88_v1/NATIVE_COLLECTION.json';solpath=rec/'two_ca_solvent_c5q88_v1/FINAL_COLLECTION.json';auditpath=rec/'refold/audit_q88jh5_v1/execution/FINAL_AUDIT.json';preppath=rec/'spectator_preparation_c5q88_v1/execution/FINAL_PREPARATION.json'
 native={r['source_case_id']:r for r in read(nativepath)['cases']};scores={r['source_case_id']:r for r in read(solpath)['rows']};preps={r['case_id']:r for r in read(preppath)['rows']};cases=[]
 for a in read(auditpath)['rows']:
  cid=a['new_source_id'];n=native[cid];p=preps[cid];s=scores[cid]
  cases.append({'case_id':cid,'family':'new_2Ca','sample':a['sample'],'canonical':False,'raw':a['model'],'confidence':confidence_pins(a['model']),
      'metal':a['active_Ca'],'spectator':a['spectator_Ca'],'roles':p['source']['roles'],'pqq':p['source']['pqq'],
      'existing_admission':p['status'],'mapping':p['maps']['Ca'],'core':p['original_core']['parent'],
      'geometries':{c:n['matrix']['Ca'][c]['xyz'] for c in ('origin','adaptive_Ca','adaptive_La')},
      'selected':{z:s['pool']['operational'][z]['candidate'] for z in ('Ca','La')},'decisions':s.get('decisions',{}),'saved_score_pool':s['pool']})
 env={r['case_id']:r for r in read(envpath)['cases'] if r['case_id'].startswith('q88jh5-') and '__conditioned_La__' in r['case_id']}
 for r in read(oldpath)['cases']:
  if not r['case_id'].startswith('q88jh5-') or '__conditioned_La__' not in r['case_id']:continue
  s=r['source'];cid=r['case_id'];p=env.get(cid);raw=s['source_structure'];i=int(cid.rsplit('-',1)[-1])
  cases.append({'case_id':cid,'family':'saved_La_control','sample':i,'canonical':bool(s.get('canonical_coordinate_match')),'raw':raw,'confidence':confidence_pins(raw),
     'metal':s['raw_source_metal'],'spectator':None,'roles':s['roles'],'pqq':s['pqq'],'existing_admission':r['status'],'existing_reason':r.get('reason'),
     'mapping':p['maps']['Ca'] if p else None,'core':p['original_core']['parent'] if p else None,
     'geometries':{'origin':p['representations']['context']['endpoints']['Ca']['xyz']} if p else {},'selected':{},'decisions':{}})
 assert len(cases)==10 and sum(c['canonical'] for c in cases)==1
 deps={};pending=['carve_with_pqq','pqq_microstates','detect_pqq','coordination_policy']
 while pending:
  name=pending.pop();path=ROOT/'scripts'/(name+'.py')
  if name in deps or not path.exists():continue
  deps[name]=pin(path)
  for n in ast.walk(ast.parse(path.read_text())):
   if isinstance(n,ast.Import):pending.extend(a.name.split('.')[0] for a in n.names)
   elif isinstance(n,ast.ImportFrom) and n.module:pending.append(n.module.split('.')[0])
 inputs={}
 for c in cases:
  for r in [c['raw'],*c['confidence'].values(),c['mapping'],c['core'],*c['geometries'].values()]:
   if r:verify(r);inputs[r['path']]=r
 m={'protocol':'PedE_saved_geometry_confidence_v1','driver':pin(__file__),'plan':pin(HERE/'PLAN.md'),'python':pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'),
    'input_collections':[pin(p) for p in (nativepath,solpath,auditpath,preppath,oldpath,envpath)],'inputs':list(inputs.values()),'code':deps,'cases':sorted(cases,key=lambda c:(c['family'],c['sample'])),
    'output':str(out/'execution'),'new_energy_calls':0,'new_folds':0,'new_thresholds':False,'numerical_threads':1}
 write(out/'manifest.json',m);print(json.dumps({'manifest':pin(out/'manifest.json'),'cases':len(cases),'executed':False}))
def key(meta):return (meta['chain'],int(meta['resnum']),meta.get('icode',meta.get('insertion_code','')),meta['atom'])
def selkey(sel,atom):return (sel['chain'],int(sel['resnum']),sel.get('icode',''),atom)
def stat(v):
 a=np.asarray(v,float)
 return {'n':int(a.size),'mean':float(a.mean()),'min':float(a.min()),'max':float(a.max())} if a.size else {'n':0,'unavailable':True}
def dihedral(a,b,c,d):
 b0=-(b-a);b1=c-b;b1=b1/np.linalg.norm(b1);b2=d-c;v=b0-np.dot(b0,b1)*b1;w=b2-np.dot(b2,b1)*b1
 return float(np.degrees(np.arctan2(np.dot(np.cross(b1,v),w),np.dot(v,w))))
def metrics(atommap,metal_xyz,roles,pqqsel,schema,protonated=frozenset()):
 donors={};chi={};groups={};qualified=[]
 for role,sel in {**roles,'pqq_cofactor':pqqsel}.items():
  members={k[3]:v for k,v in atommap.items() if k[:3]==selkey(sel,'')[:3]};groups[role]=members
  d={}
  for name,v in members.items():
   dtype=('cofactor_'+v['element'] if name in schema.direct_donors else None) if role=='pqq_cofactor' else policy.donor_type(sel['resname'],name,v['element'])
   if selkey(sel,name) in protonated and (dtype=='sidechain_N' or (sel.get('resname') in ('ASP','GLU') and dtype=='sidechain_O') or (dtype=='backbone_O' and name in ('O','OXT'))):dtype=None
   if v['element'] in ('O','N'):
    dist=float(np.linalg.norm(v['xyz']-metal_xyz));d[name]={'distance_A':dist,'typed_donor':dtype,'qualifies_3p1A':bool(dtype and dist<=3.1)}
    if dtype and dist<=3.1:qualified.append({'role':role,'atom':name,'distance_A':dist,'type':dtype})
  donors[role]={'atoms':d,'typed_denticity_3p1A':sum(v['qualifies_3p1A'] for v in d.values()),'sorted_O_distances_A':sorted(v['distance_A'] for n,v in d.items() if members[n]['element']=='O')}
  definitions={'chi1':('N','CA','CB','CG'),'chi2':('CA','CB','CG','CD'),'chi3':('CB','CG','CD','OE1')}
  if sel.get('resname') in ('ASP','ASN'):definitions['chi2']=('CA','CB','CG','OD1');definitions.pop('chi3')
  if sel.get('resname')=='SER':definitions={'chi1':('N','CA','CB','OG')}
  if sel.get('resname')=='ARG':definitions={'chi1':('N','CA','CB','CG'),'chi2':('CA','CB','CG','CD'),'chi3':('CB','CG','CD','NE')}
  chi[role]={n:dihedral(*(members[a]['xyz'] for a in names)) for n,names in definitions.items() if all(a in members for a in names)}
 p=groups['pqq_cofactor'];bl=[{'atoms':[a,b],'distance_A':float(np.linalg.norm(p[a]['xyz']-p[b]['xyz']))} for a,b in schema.bonds]
 carboxyl={x[0] for x in schema.carboxyls.values()};ring=sorted(n for n,v in p.items() if v['element'] in ('C','N') and n not in carboxyl)
 xyz=np.asarray([p[n]['xyz'] for n in ring]);u,s,vh=np.linalg.svd(xyz-xyz.mean(axis=0),full_matrices=False);deviation=(xyz-xyz.mean(axis=0))@vh[-1]
 return {'donors':donors,'typed_declared_role_CN_3p1A':len(qualified),'typed_contacts':qualified,'chi_degrees':chi,
   'PQQ':{'schema':schema.schema_id,'bond_lengths':bl,'existing_connectivity_gate_pass':all(.9<=b['distance_A']<=1.85 for b in bl),
          'ring_plane_atoms':ring,'ring_plane_RMS_A':float(np.sqrt(np.mean(deviation**2))),'ring_plane_max_A':float(np.max(np.abs(deviation)))}}
def raw_case(c):
 st=gemmi.read_structure(str(verify(c['raw'])));assert len(st)==1
 rows=[];amap={};metal=None;metalatom=None;pqqres=None
 for ch in st[0]:
  for r in ch:
   for a in r:
    k=(ch.name,r.seqid.num,r.seqid.icode.strip(),a.name.strip());v={'xyz':np.array([a.pos.x,a.pos.y,a.pos.z]),'element':a.element.name,'resname':r.name,'plddt':float(a.b_iso)}
    if k in amap:raise ValueError('ambiguous source atom '+str(k))
    amap[k]=v;rows.append((k,v))
    if k==selkey(c['metal'],c['metal']['atom']):metal=v['xyz'];metalatom=a
   if (ch.name,r.seqid.num,r.seqid.icode.strip())==selkey(c['pqq'],'')[:3]:pqqres=r
 if metal is None or pqqres is None:raise ValueError('declared source site/PQQ missing')
 schema=pqq.identify_naming_schema(pqqres);out=metrics(amap,metal,c['roles'],c['pqq'],schema)
 pqqres.name='PQQ'
 ms=c['metal'];site=carver.MetalSite(ms['chain'],ms['resnum'],ms.get('icode',''),ms['resname'],ms['atom'],metalatom.element.name,metalatom)
 ps=c['pqq'];pk=carver.ResidueKey(ps['chain'],ps['resnum'],ps.get('icode',''),'PQQ')
 contacts=carver._classify_contacts(st[0],site,{pk:(pqqres,SimpleNamespace(direct_donors=schema.direct_donors))})
 out['full_source_existing_typed_contacts']={name:[v.metadata() for v in vs] for name,vs in zip(('qualified_3p1A','included_3p3A','sulfur_diagnostic','untyped_ON_3p1A'),contacts)}
 out['full_source_contact_scope']='raw source heavy geometry; selected PQQ name aliased in memory only, oxidized baseline donor schema; no source protonation inferred'
 out['source_atom_count']=len(rows);out['source_metal']=c['metal'];out['existing_admission']=c['existing_admission'];out['existing_reason']=c.get('existing_reason')
 confidence={};cf=c['confidence'].get('confidences');sc=c['confidence'].get('summary_confidences')
 groups={'active_metal':[selkey(c['metal'],c['metal']['atom'])],'PQQ':[k for k in amap if k[:3]==selkey(c['pqq'],'')[:3]],**{role:[k for k in amap if k[:3]==selkey(sel,'')[:3]] for role,sel in c['roles'].items()}}
 if c.get('spectator'):groups['spectator']=[selkey(c['spectator'],c['spectator']['atom'])]
 groups['declared_typed_donor_atoms']=[selkey(({**c['roles'],'pqq_cofactor':c['pqq']})[role],name) for role,v in out['donors'].items() for name,d in v['atoms'].items() if d['typed_donor']]
 confidence['CIF_atom_plddt']={name:stat([amap[k]['plddt'] for k in keys]) for name,keys in groups.items()}
 if sc:confidence['summary']=read(verify(sc))
 if cf:
  x=read(verify(cf));chains=[k[0] for k,v in rows]
  if x['atom_chain_ids']!=chains:raise ValueError('confidence atom-chain ordering differs from exact model')
  if len(x['atom_plddts'])!=len(rows):raise ValueError('confidence atom count differs')
  confidence['atom_plddt_CIF_max_difference']=float(np.max(np.abs(np.asarray(x['atom_plddts'])-np.array([v['plddt'] for k,v in rows]))))
  confidence['atom_plddt']={name:stat([x['atom_plddts'][i] for i,(k,v) in enumerate(rows) if k in keys]) for name,keys in groups.items()}
  tokens=[(ch,int(r)) for ch,r in zip(x['token_chain_ids'],x['token_res_ids'])];tgroups={name:[i for i,k in enumerate(tokens) if k in {(a[0],a[1]) for a in keys}] for name,keys in groups.items()}
  pae=np.asarray(x['pae']);contacts=np.asarray(x['contact_probs']);confidence['directed_PAE_A']={};confidence['contact_probability']={}
  for name,indices in tgroups.items():
   if name=='active_metal':continue
   ii=tgroups['active_metal'];confidence['directed_PAE_A']['metal_to_'+name]=stat(pae[np.ix_(ii,indices)].ravel());confidence['directed_PAE_A'][name+'_to_metal']=stat(pae[np.ix_(indices,ii)].ravel());confidence['contact_probability']['metal_to_'+name]=stat(contacts[np.ix_(ii,indices)].ravel())
 else:confidence['detail_status']='confidence JSON unavailable; raw CIF atom values retained'
 out['confidence']=confidence
 return out,amap,metal,schema

def xyz_read(rec):
 lines=verify(rec).read_text().splitlines();n=int(lines[0]);rows=[x.split() for x in lines[2:] if x.strip()];assert len(rows)==n
 return [r[0] for r in rows],np.asarray([[float(v) for v in r[1:4]] for r in rows])
def local_case(c,schema,rawmap):
 if not c['mapping']:return {'status':'unavailable_no_existing_prepared_map','existing_reason':c.get('existing_reason')}
 mapping=read(verify(c['mapping']))['context'];meta=mapping['source_atom_metadata'];out={};orig=None;orig_metrics=None
 for name,pinrec in c['geometries'].items():
  symbols,xyz=xyz_read(pinrec);assert len(xyz)==len(meta)
  if name=='origin':orig=xyz
  if orig is None:raise ValueError('origin must precede candidates')
  amap={key(m):{'xyz':xyz[i],'element':symbols[i],'resname':m['resname']} for i,m in enumerate(meta) if m}
  protonated=set()
  for i,j in mapping['bonds']:
   if symbols[i] in ('H','D') and meta[j]:protonated.add(key(meta[j]))
   if symbols[j] in ('H','D') and meta[i]:protonated.add(key(meta[i]))
  m=metrics(amap,xyz[0],c['roles'],c['pqq'],schema,protonated)
  anchored=dict(amap);borrowed=[]
  for sel in c['roles'].values():
   for atom in ('N','CA','C','O'):
    k=selkey(sel,atom)
    if k not in anchored and k in rawmap:anchored[k]=rawmap[k];borrowed.append(list(k))
  m['chi_with_fixed_source_backbone_degrees']=metrics(anchored,xyz[0],c['roles'],c['pqq'],schema,protonated)['chi_degrees'];m['fixed_source_backbone_reference_atoms']=borrowed
  matches=[(i,key(v)) for i,v in enumerate(meta) if v and symbols[i] not in ('H','D') and key(v) in rawmap]
  m['raw_source_mapped_heavy_displacements']={'n':len(matches),'max_A':max((float(np.linalg.norm(xyz[i]-rawmap[k]['xyz'])) for i,k in matches),default=None),'by_role':{role:stat([float(np.linalg.norm(xyz[i]-rawmap[k]['xyz'])) for i,k in matches if k[:3]==selkey(sel,'')[:3]]) for role,sel in {**c['roles'],'PQQ':c['pqq']}.items()}}
  delta=xyz-orig;role_motion={}
  for role,sel in {**c['roles'],'PQQ':c['pqq']}.items():
   inds=[i for i,v in enumerate(meta) if v and key(v)[:3]==selkey(sel,'')[:3] and symbols[i] not in ('H','D')]
   role_motion[role]={'mapped_heavy_atoms':len(inds),'RMS_displacement_A':float(np.sqrt(np.mean(np.sum(delta[inds]**2,axis=1)))) if inds else None,'max_displacement_A':float(np.max(np.linalg.norm(delta[inds],axis=1))) if inds else None}
  bonds=[{'indices':[i,j],'before_A':float(np.linalg.norm(orig[i]-orig[j])),'after_A':float(np.linalg.norm(xyz[i]-xyz[j]))} for i,j in mapping['bonds']]
  m.update(xyz=pinrec,selected_by=[z for z,v in c['selected'].items() if v==name],metal_displacement_A=float(np.linalg.norm(delta[0])),role_motion=role_motion,
           bond_max_change_A=max((abs(v['after_A']-v['before_A']) for v in bonds),default=0),mapped_bond_lengths=bonds)
  if name=='origin':orig_metrics=m
  m['chi_change_degrees']={role:{k:float((value-orig_metrics['chi_degrees'][role][k]+180)%360-180) for k,value in d.items() if k in orig_metrics['chi_degrees'].get(role,{})} for role,d in m['chi_degrees'].items()}
  m['chi_change_with_fixed_source_backbone_degrees']={role:{k:float((value-orig_metrics['chi_with_fixed_source_backbone_degrees'][role][k]+180)%360-180) for k,value in d.items() if k in orig_metrics['chi_with_fixed_source_backbone_degrees'].get(role,{})} for role,d in m['chi_with_fixed_source_backbone_degrees'].items()}
  out[name]=m
 return out

def aligned_compare(c,amap,metal,ref_amap,ref_metal,schema):
 keys=[selkey(sel,a) for sel in c['roles'].values() for a in ('N','CA','C','O')]
 common=[k for k in keys if k in amap and k in ref_amap];assert len(common)==20
 a=np.asarray([amap[k]['xyz'] for k in common]);b=np.asarray([ref_amap[k]['xyz'] for k in common]);ac=a.mean(axis=0);bc=b.mean(axis=0)
 u,s,vh=np.linalg.svd((a-ac).T@(b-bc));d=np.eye(3);d[-1,-1]=np.linalg.det(u@vh);rot=u@d@vh
 transform=lambda x:(x-ac)@rot+bc
 names=sorted(schema.expected_elements);pp=c['pqq'];aP={n:transform(amap[selkey(pp,n)]['xyz']) for n in names};bP={n:ref_amap[selkey(pp,n)]['xyz'] for n in names}
 sq={n:float(np.sum((aP[n]-bP[n])**2)) for n in names};sym=sum(sq.values())
 for carbon,left,right in schema.carboxyls.values():
  direct=sq[left]+sq[right];swap=float(np.sum((aP[left]-bP[right])**2)+np.sum((aP[right]-bP[left])**2));sym+=min(direct,swap)-direct
 return {'anchor_keys':[list(k) for k in common],'anchor_RMSD_A':float(np.sqrt(np.mean(np.sum((transform(a)-b)**2,axis=1)))),
         'metal_displacement_after_anchor_fit_A':float(np.linalg.norm(transform(metal)-ref_metal)),
         'PQQ_named_atom_RMSD_A':float(np.sqrt(sum(sq.values())/len(names))),'PQQ_symmetry_corrected_RMSD_A':float(np.sqrt(sym/len(names)))}

def execute(manifest):
 if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURM_JOB_NODELIST'):raise RuntimeError('real Slurm allocation required before geometry imports')
 for key in THREADS:os.environ[key]='1'
 m=read(manifest)
 for p in [m['driver'],m['plan'],m['python'],*m['inputs'],*m['input_collections'],*m['code'].values()]:verify(p)
 out=Path(m['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();receipt={'manifest':pin(manifest),'job_id':os.environ['SLURM_JOB_ID'],'hostname':socket.gethostname(),'status':'running','new_energy_calls':0,'new_folds':0,'numerical_threads':1};write(out/'START.json',receipt)
 try:
  global np,gemmi,pqq,policy,carver
  import numpy as np,gemmi
  sys.path.insert(0,str(ROOT/'scripts'));import pqq_microstates as pqq,coordination_policy as policy,carve_with_pqq as carver
  results=[];raw={}
  for c in m['cases']:
   r={'case_id':c['case_id'],'family':c['family'],'sample':c['sample'],'canonical':c['canonical'],'selected':c['selected'],'decisions':c['decisions'],'raw_source':c['raw']}
   try:
    r['raw'],amap,metal,schema=raw_case(c);raw[c['case_id']]=(amap,metal,schema);r['local']=local_case(c,schema,amap);r['status']='complete'
   except Exception as e:r.update(status='partial_or_unavailable',reason=str(e),traceback=traceback.format_exc())
   results.append(r)
  ref=next(c for c in m['cases'] if c['canonical']);ra,rm,rs=raw[ref['case_id']]
  for c,r in zip(m['cases'],results):
   if c['case_id'] in raw:
    aa,mm,ss=raw[c['case_id']];r['canonical_alignment']=aligned_compare(c,aa,mm,ra,rm,ss)
   write(out/(c['family']+'_sample'+str(c['sample'])+'.json'),r)
  result={'manifest':pin(manifest),'terminal':True,'canonical_control':ref['case_id'],'rows':results,'new_energies':0,'new_folds':0,'interpretation':'descriptive saved geometry/confidence; no unique structural/model causal attribution or new quality threshold'}
  write(out/'RESULT.json',result);receipt.update(status='complete',result=pin(out/'RESULT.json'))
 except Exception as e:receipt.update(status='failed',reason=str(e),traceback=traceback.format_exc());raise
 finally:
  u=resource.getrusage(resource.RUSAGE_SELF);receipt.update(wall_seconds=time.monotonic()-start,CPU_seconds=u.ru_utime+u.ru_stime,peak_RSS_KiB=u.ru_maxrss);write(out/'EXECUTION.json',receipt)
if __name__=='__main__':
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);a=s.add_parser('freeze');a.add_argument('--workspace',required=True);a=s.add_parser('execute');a.add_argument('--manifest',required=True);a=p.parse_args();freeze(a.workspace) if a.op=='freeze' else execute(a.manifest)
