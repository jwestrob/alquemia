#!/usr/bin/env python3
"""Separate same-policy q0 restarts for the two NODE_FAIL-interrupted MxaF La searches."""
import argparse,copy,fcntl,hashlib,importlib.util,json,os,socket,subprocess,sys,time,traceback
from pathlib import Path

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
IDS=[f'p12293-pqq-la_model__conditioned_Ca__seed-1_sample-{i}__envelope_e88ef960dad9c35e' for i in (1,4)]
ORDER=['origin','adaptive_Ca','adaptive_La']
THREAD_KEYS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')
PROTOCOL='MxaF_Ca_conditioned_same_policy_La_q0_restart_v1'
def pin(path):
 p=Path(path).resolve();h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
 if pin(p['path'])!=p:raise ValueError('changed artifact: '+p['path'])
 return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);temp=p.with_name(p.name+'.tmp.'+str(os.getpid()))
 with temp.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 try:os.link(temp,p)
 finally:temp.unlink()
def reuse_check(cell,model):
 for k in ('MACE','xyz','forces'):verify(cell[k])
 r=read(cell['MACE']);req=read(r['request'])
 assert cell['status']==r['status']=='complete' and r['model']==model
 assert req['xyz']==cell['xyz'] and req['charge']==cell['charge'] and req['multiplicity']==cell['multiplicity']
 assert r['energy_eV']==cell['energy_eV'] and r['forces']==cell['forces']
 return cell

def prepare(primary,collection,output):
 m=read(primary);old=read(collection);assert old['manifest']==pin(primary) and old['terminal'] is True
 assert m['implementation']['sha256']=='e3070a408777d8754506076440414f974ed01fc1207122bad812b88bad6718d0';verify(m['implementation'])
 out=output.resolve();out.mkdir(parents=True,exist_ok=False);records=[]
 for cid in IDS:
  source=next(c for c in old['cases'] if c['case_id']==cid)
  sourcefile=Path(primary).parent/'sources'/cid/'SOURCE_RESULT.json';assert read(sourcefile)==source
  assert source['native_pool_status']=='unavailable';pm=read(source['proposal_manifest'])
  assert pm['batch_manifest']==pin(primary) and pm['settings']==m['settings'] and pm['model']==m['model']
  la=next(t for t in pm['tasks'] if t['metal']=='La');ca=read(source['proposals']['Ca']);interrupted=read(source['proposals']['La'])
  assert ca['status']=='proposal_available' and ca['manifest']==source['proposal_manifest']
  assert interrupted['status']=='unavailable' and interrupted['reason']=='interrupted_started_search_preserved_no_q0_restart'
  assert interrupted['manifest']==source['proposal_manifest'] and interrupted['task_id']==la['task_id']
  for z in ('Ca','La'):
   for q in ('origin','adaptive_Ca'):reuse_check(source['matrix'][z][q],m['model'])
  for t in pm['tasks']:
   for k in ('xyz','mapping','source_preparation'):verify(t[k])
   verify(t['origin_reuse']['point']['MACE']);verify(t['origin_reuse']['point']['forces'])
  records.append({'case_id':cid,'primary_source':pin(sourcefile),'primary_proposal_manifest':source['proposal_manifest'],
    'reused_Ca_proposal':source['proposals']['Ca'],'interrupted_La_search':source['proposals']['La']})
 wm=copy.deepcopy(m);wm['cases']=[dict(next(c for c in m['cases'] if c['case_id']==cid),lane=0) for cid in IDS]
 bases=[t for t in m['tasks'][:m['origin_task_count']] if t['case_id'] in IDS]
 fixtures=[next(t for t in m['tasks'] if t['task_id']==f['task_id']) for f in m['fixtures']]
 wm.update(protocol_id=PROTOCOL+'_transport',cases=wm['cases'],tasks=bases+fixtures,origin_task_count=4,lanes=1,maximum_active_numerical_threads=7,
           parent_manifest=pin(primary),transport_only_changes=['source subset','single lane','seven total numerical threads'],
           finite_scope={'sources':2,'new_La_optimizer_starts':2,'new_Ca_optimizer_starts':0,'qualification_calls':2,'maximum_new_cross_calls':2,'new_DFT_calls':0,'new_GFN2_calls':0})
 wp=out/'worker_manifest.json';write(wp,wm)
 for r in records:
  pm=copy.deepcopy(read(r['primary_proposal_manifest']));pm['batch_manifest']=pin(wp)
  pm['technical_recovery']={'source_proposal_manifest':r['primary_proposal_manifest'],'interrupted_La_search':r['interrupted_La_search'],'restart':'same optimizer from original q0; no optimizer-state continuation'}
  p=out/'sources'/r['case_id']/'proposal_manifest.json';write(p,pm);r['restart_proposal_manifest']=pin(p)
 plan={'protocol_id':PROTOCOL,'implementation':pin(__file__),'primary_manifest':pin(primary),'primary_collection':pin(collection),'worker_manifest':pin(wp),'cases':records,
       'settings':m['settings'],'pool_settings':m['pool_settings'],'model':m['model'],'gates':m['gates'],
       'resources':{'CPUs':7,'model_workers':1,'threads_per_model':7,'parent_threads':1,'parent_and_model_work_sequential':True},
       'scope':wm['finite_scope'],'authorization':'Jacob explicitly authorized separate technical redo of previously unavailable cases, including same-policy q0 restart',
       'old_attempts_preserved':True,'calibration_changed':False,'production_changed':False,'partial_old_trial_reuse':'No optimizer history or partial La trials reused; repeated trial calls are fresh measured work.'}
 write(out/'manifest.json',plan)
 write(out/'PREPARATION.json',{'manifest':pin(out/'manifest.json'),'sources':2,'reused_native_cells':8,'reused_Ca_proposals':2,'archived_q0_reuse':2,'prepared_metadata_only':True,'new_molecular_calls':0,'scientific_preflight':'unrun; required under Slurm before model load'})
 print(json.dumps({'manifest':str(out/'manifest.json'),'sources':2,'new_molecular_calls':0}))

def load(plan):
 p=read(plan);assert p['protocol_id']==PROTOCOL and p['implementation']==pin(__file__)
 old=read(p['primary_manifest']);m=read(p['worker_manifest']);assert m['parent_manifest']==p['primary_manifest']
 assert p['settings']==m['settings']==old['settings'] and p['model']==m['model']==old['model'] and m['gates']==old['gates']
 for key in ('model','settings','pool_settings','optimizer_software','optimizer_python','model_python','calculator_source','source_modules','source_snapshot','optimizer_executable','model_executable','fixtures','implementation'):
  assert m[key]==old[key], 'transport changed scientific/runtime field: '+key
 expected=[t for t in old['tasks'][:old['origin_task_count']] if t['case_id'] in IDS]+[next(t for t in old['tasks'] if t['task_id']==f['task_id']) for f in old['fixtures']]
 assert m['tasks']==expected and [r['case_id'] for r in p['cases']]==IDS
 assert [c['case_id'] for c in m['cases']]==IDS and all(c['lane']==0 for c in m['cases']) and len(m['tasks'])==6 and m['origin_task_count']==4
 for key in ('implementation','calculator_source','optimizer_executable','model_executable'):verify(m[key])
 verify(m['model']['checkpoint'])
 for v in m['source_modules'].values():verify(v)
 for f in m['fixtures']:
  for key in ('xyz','forces','native_receipt'):verify(f[key])
 for r in p['cases']:
  source=read(r['primary_source']);prior=read(r['primary_proposal_manifest']);fresh=read(r['restart_proposal_manifest'])
  assert fresh['tasks']==prior['tasks'] and fresh['source']==prior['source'] and fresh['selector']==prior['selector']
  assert fresh['batch_manifest']==p['worker_manifest'] and fresh['settings']==m['settings']
  for z in ('Ca','La'):
   for q in ('origin','adaptive_Ca'):reuse_check(source['matrix'][z][q],m['model'])
  read(r['reused_Ca_proposal']);read(r['interrupted_La_search'])
 return p,m

def allocation():
 job=os.environ.get('SLURM_JOB_ID');assert job and job.isdigit(),'scientific operations require actual Slurm'
 raw=subprocess.run(['scontrol','show','job','-o',job],text=True,capture_output=True,check=True).stdout
 fields=dict(t.split('=',1) for t in raw.split() if '=' in t)
 assert fields['JobState']=='RUNNING' and int(fields['NumCPUs'])==7 and int(fields['NumNodes'])==1
 assert int(os.environ['SLURM_CPUS_ON_NODE'])==7 and len(os.sched_getaffinity(0))>=7
 assert all(os.environ.get(k)=='1' for k in THREAD_KEYS)
 return {'job_id':job,'host':socket.gethostname(),'affinity':sorted(os.sched_getaffinity(0)),'scheduler_record':raw,'allocated_CPUs':7}

def backend(m):
 spec=importlib.util.spec_from_file_location('mxaf_frozen_cpu',verify(m['implementation']));b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b);b.modules(m);return b

def scientific_preflight(p,m,b):
 import numpy as np,scipy
 import scipy.optimize._slsqplib as kernel
 import scipy.optimize._slsqp_py as wrapper
 from mace_site_kinematics import Kinematics
 from accommodation_torsion_profiles import geometry_check
 import motion_envelope_run as builder
 import slsqp_precision
 assert scipy.__version__==m['optimizer_software']['version']
 assert pin(kernel.__file__)==m['optimizer_software']['kernel'] and pin(wrapper.__file__)==m['optimizer_software']['wrapper']
 assert slsqp_precision.engine()[0].completion.SETTINGS==m['settings']
 checks=[]
 for r in p['cases']:
  source=read(r['primary_source']);pm=read(r['restart_proposal_manifest'])
  tasks,selector=builder.build_tasks(copy.deepcopy(source['source']['origin_row']),m['model'])
  assert tasks==pm['tasks'] and selector==pm['selector'],'paired force selector/task changed'
  for t in tasks:
   kin=Kinematics(read(t['mapping'])['context']);q=np.zeros(len(kin.modes));a=b.xyz(verify(t['xyz']))
   residual=float(np.max(np.abs(kin.evaluate(q)[1]-np.asarray([v[1:] for v in a]))))
   assert residual<=1e-12 and geometry_check(kin,q,[v[0] for v in a])['pass']
   assert all(x==0 for x in t['origin_reuse']['point']['active_q_radian'])
   checks.append({'task_id':t['task_id'],'q0_coordinate_max_error_A':residual,'selector_unchanged':True})
 return {'status':'pass','checks':checks,'new_molecular_calls':0}

def recover_source(p,m,b,record,transport):
 from nikasha_pool import same_geometry
 from mace_hybrid import write_xyz
 import slsqp_precision
 pp=verify(record['restart_proposal_manifest']);dest=pp.parent/'SOURCE_RESULT.json'
 if dest.exists():
  result=read(dest);assert result['technical_recovery']['manifest']==pin(Path(p['_path']));return result
 old=read(record['primary_source']);result=copy.deepcopy(old);pm=read(pp);task=next(t for t in pm['tasks'] if t['metal']=='La')
 rp=b.proposal(pp,task,transport,slsqp_precision.engine()[0].completion);proposal=read(rp)
 result.update(batch_manifest=p['worker_manifest'],primary_batch_manifest=old['batch_manifest'],technical_recovery={'protocol_id':PROTOCOL,'manifest':pin(Path(p['_path'])),'primary_source':record['primary_source'],'primary_La_interruption':record['interrupted_La_search'],'original_selection_unchanged':True},
               proposal_manifest=record['restart_proposal_manifest'],proposal_manifests={'Ca':record['primary_proposal_manifest'],'La':record['restart_proposal_manifest']})
 result['proposals']['La']=rp
 if proposal['status']=='proposal_available':
  point=proposal['proposal'];name='adaptive_La';candidate={'status':'available','xyz':point['coordinate'],'origin':'La','proposal_receipt':rp,'boundary_flag':proposal['boundary_flag']}
  result['candidates'][name]=candidate;native=read(point['MACE']);result['matrix']['La'][name]=b.native_cell(task,native,point['MACE'],point['coordinate'],name)
  representative=next((q for q in ORDER[:-1] if same_geometry(b.xyz(verify(candidate['xyz'])),b.xyz(verify(result['candidates'][q]['xyz'])))),name)
  result['aliases'][name]={'representative':representative,'coordinate_tolerance_A':1e-12}
  if representative!=name:
   result['matrix']['Ca'][name]={**result['matrix']['Ca'][representative],'alias_of':representative,'candidate_origin':name}
  else:
   t=next(t for t in m['tasks'] if t['case_id']==record['case_id'] and t['metal']=='Ca');cd=pp.parent/'cross_pool'/'adaptive_La__Ca';cd.mkdir(parents=True,exist_ok=True);xp=cd/'context.xyz'
   rows=b.xyz(verify(candidate['xyz']));rows[0]=('Ca',*rows[0][1:])
   if not xp.exists():write_xyz(xp,rows)
   assert same_geometry(b.xyz(xp),rows) and b.xyz(xp)[0][0]=='Ca','existing cross coordinates differ'
   ct={'task_id':record['case_id']+'__adaptive_La__Ca','case_id':record['case_id'],'metal':'Ca','xyz':pin(xp),'charge':t['charge'],'multiplicity':t['multiplicity']}
   cp=cd/'manifest.json';cm={'kind':'cross_pool','batch_manifest':p['worker_manifest'],'model':m['model'],'tasks':[ct],'candidate_source':candidate}
   b.write_or_check(cp,cm);req={'manifest':pin(cp),'task_id':ct['task_id'],'xyz':ct['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity']}
   raw,receipt=b.cell(cd,req,transport,m['model']);result['matrix']['Ca'][name]=b.native_cell(t,raw,receipt,ct['xyz'],name)
 else:
  result['candidates']['adaptive_La']={'status':'unavailable','reason':proposal.get('reason',proposal['status']),'proposal_receipt':rp}
  for z in ('Ca','La'):result['matrix'][z]['adaptive_La']={'status':'unavailable','reason':proposal.get('reason',proposal['status'])}
 complete=all(result['matrix'][z][q]['status']=='complete' for z in ('Ca','La') for q in ORDER)
 result['status']='native_matrix_complete' if complete else 'unavailable';result['native_pool_status']='available' if complete else 'unavailable'
 write(dest,result);return result

def collect(plan,attempt,error=None):
 p,m=load(plan);base=Path(plan).parent;cases=[]
 for rec in p['cases']:
  path=Path(rec['restart_proposal_manifest']['path']).parent/'SOURCE_RESULT.json'
  if path.exists():cases.append(read(path))
 paths=list(base.glob('sources/**/mace_result.json'))+list(base.glob('attempts/*/lane0/qualification/*/mace_result.json'))
 receipts=[read(x) for x in paths]
 return {'protocol_id':PROTOCOL,'manifest':pin(plan),'worker_manifest':p['worker_manifest'],'primary_collection':p['primary_collection'],'terminal':len(cases)==2,'denominator':2,'cases':cases,
         'pending_source_ids':[r['case_id'] for r in p['cases'] if r['case_id'] not in {c['case_id'] for c in cases}],
         'native_matrix_complete':sum(c['native_pool_status']=='available' for c in cases),'attempt':str(attempt),'error':error,
         'costs':{'new_recorded_model_calls':len(receipts),'new_recorded_forwards':sum(r['actual_model_forwards'] for r in receipts),'new_recorded_wall_seconds_sum':sum(r['wall_seconds'] for r in receipts),'new_recorded_process_cpu_seconds':sum(r['process_cpu_seconds'] for r in receipts),
                  'reused_native_cells':8,'prior_costs_in_primary_collection':p['primary_collection'],'started_calls_without_terminal_receipt':[str(x) for x in base.glob('**/CPU_CALL_STARTED.json') if not (x.parent/'mace_result.json').exists()]},
         'production_changed':False,'calibration_changed':False,'new_GFN2_calls':0,'new_DFT_calls':0}

def execute(plan):
 placement=allocation();p,m=load(plan);base=Path(plan).parent
 with (base/'RUN.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (base/'FINAL_COLLECTION.json').exists():raise RuntimeError('recovery already collected; inspect existing results')
  attempt=base/'attempts'/'001';attempt.mkdir(parents=True,exist_ok=False);(attempt/'lane0').mkdir()
  write(attempt/'PLACEMENT.json',placement)
  transport=None;error=None;started=time.monotonic();b=backend(m);old_handlers=b.handlers()
  try:
   write(attempt/'PREFLIGHT.json',scientific_preflight(p,m,b))
   transport=b.WarmCPU(verify(p['worker_manifest']),m,attempt,0);checks=[]
   import numpy as np
   for f in m['fixtures']:
    t=next(t for t in m['tasks'] if t['task_id']==f['task_id']);req={'manifest':p['worker_manifest'],'task_id':t['task_id'],'xyz':f['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity']}
    r,rp=b.cell(attempt/'lane0'/'qualification'/f['id'],req,transport,m['model']);check={'fixture':f['id'],'receipt':rp,'pass':False}
    if r['status']=='complete':
     de=(r['energy_eV']-f['energy_eV'])*b.EV_TO_KCAL;df=float(np.max(np.abs(np.load(verify(r['forces']),allow_pickle=False)-np.load(verify(f['forces']),allow_pickle=False))))
     check.update(energy_error_kcal_mol=de,force_max_error_eV_A=df,**{'pass':abs(de)<=m['gates']['energy_kcal_mol'] and df<=m['gates']['force_max_eV_A']})
    checks.append(check)
   write(attempt/'QUALIFICATION.json',{'checks':checks,'pass':all(x['pass'] for x in checks)})
   if not all(x['pass'] for x in checks):raise RuntimeError('original CPU qualification failed; no restart performed')
   p['_path']=str(plan)
   for rec in p['cases']:recover_source(p,m,b,rec,transport)
  except BaseException as exc:error=type(exc).__name__+': '+str(exc);traceback.print_exc()
  finally:
   if transport:
    try:
     if error:transport.abort()
     else:transport.close()
    except BaseException as exc:error=(error+'; ' if error else '')+str(exc);transport.abort()
   import signal
   for s,h in old_handlers.items():signal.signal(s,h)
   result=collect(plan,attempt,error);elapsed=time.monotonic()-started
   result['execution_receipt']={'placement':placement,'elapsed_seconds':elapsed,'allocated_CPU_seconds':elapsed*7,'worker_summary':pin(attempt/'lane0/CPU_WORKER_SUMMARY.json') if (attempt/'lane0/CPU_WORKER_SUMMARY.json').exists() else None}
   write(base/'FINAL_COLLECTION.json',result)
  if error:raise RuntimeError(error)
  print(json.dumps({'terminal':result['terminal'],'native_matrix_complete':result['native_matrix_complete'],'denominator':2}))

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='op',required=True)
 a=sub.add_parser('prepare');a.add_argument('--primary',type=Path,required=True);a.add_argument('--collection',type=Path,required=True);a.add_argument('--output',type=Path,required=True)
 a=sub.add_parser('execute');a.add_argument('--manifest',type=Path,required=True)
 args=parser.parse_args()
 if args.op=='prepare':prepare(args.primary.resolve(),args.collection.resolve(),args.output)
 else:execute(args.manifest.resolve())
