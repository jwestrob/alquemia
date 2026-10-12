"""Reconstruct the completed finite scout from actual pinned receipts only."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[k]='1'
import hashlib,json,statistics
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
RUN=ROOT/'workspaces/orca_parameter_audit_20261012/cpu_proposal_scout_v1'
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,p
    return Path(p['path'])
def read(p):return json.loads(Path(p).read_text())
final=read(RUN/'FINAL.json');m=read(verify(final['manifest']));worker=read(RUN/'CPU_WORKER_SUMMARY.json')
assert final['status']=='complete' and final['all_scientific_gates_pass'] and len(final['proposals'])==2
assert final['new_DFT_calls']==final['new_solvent_calls']==0
rows=[];checks=[]
for p in final['model_receipts']:
    r=read(verify(p));req=read(verify(r['request']));verify(req['xyz']);verify(r['forces'])
    assert r['manifest']==final['manifest']==req['manifest']
    assert r['status']=='complete' and r['model']==m['model'] and r['actual_model_forwards']==1
    t=next(t for t in m['tasks'] if t['task_id']==r['task_id'])
    assert (r['state']['charge'],r['state']['spin_multiplicity'],r['state']['atoms'])==(t['charge'],t['multiplicity'],202)
    assert r['state']['head']=='omol' and r['state']['graph_count']==1 and r['parameter_versions_unchanged']
    rows.append({'receipt':p,'wall_seconds':r['wall_seconds'],'CPU_seconds':r['process_cpu_seconds'],
                 'peak_RSS_KiB':r['peak_host_RSS_KiB'],'forwards':r['actual_model_forwards']})
for f,c in zip(m['fixtures'],final['replays']):
    r=read(verify(c['receipt']));a=np.load(verify(r['forces']),allow_pickle=False);b=np.load(verify(f['forces']),allow_pickle=False)
    assert a.shape==b.shape==(202,3)
    de=(r['energy_eV']-f['energy_eV'])*23.06054783061903;df=float(np.max(np.abs(a-b)))
    assert de==c['energy_error_kcal_mol'] and df==c['force_max_error_eV_A']
    checks.append({'id':f['id'],'energy_error_kcal_mol':de,'force_max_error_eV_A':df,
                   'pass':abs(de)<=m['gates']['energy_kcal_mol'] and df<=m['gates']['force_max_eV_A']})
assert len(rows)==worker['attempted_requests']==worker['actual_model_forwards']==final['actual_requests']
assert not worker['failed_requests'] and all(c['pass'] for c in checks)
proposals=[]
for r in final['proposals']:
    p=read(verify(r['receipt']));assert p['status']=='proposal_available' and r['pass']
    calls=[x for x in p['requests'] if not x.get('reused')]
    proposals.append({'metal':p['metal'],'receipt':r['receipt'],'fresh_force_calls':len(calls),
       'optimizer':p['optimizer'],'wall_seconds':p['wall_seconds'],'energy_error_kcal_mol':r['energy_error_kcal_mol'],
       'maximum_coordinate_difference_A':r['maximum_coordinate_difference_A'],'boundary_matches':r['boundary_matches'],
       'boundary_flag':r['boundary_flag'],'physical_feasible':r['final_geometry']['physical_feasible']})
summary={'status':'passed_scoped_CPU_force_and_two_proposal_qualification','implementation':pin(__file__),
 'manifest':final['manifest'],'final_receipt':pin(RUN/'FINAL.json'),'worker_summary':pin(RUN/'CPU_WORKER_SUMMARY.json'),
 'phase1_receipt':pin(RUN/'PHASE1.json'),'reconstructed_force_checks':checks,'proposals':proposals,
 'cost':{'actual_new_forward_calls':len(rows),'failed_calls':worker['failed_requests'],'model_load_seconds':worker['model_load_seconds'],
   'wall_seconds':final['wall_seconds'],'worker_CPU_seconds':worker['process_cpu_seconds'],
   'parent_CPU_seconds':final['parent_process_cpu_seconds'],'total_CPU_seconds':worker['process_cpu_seconds']+final['parent_process_cpu_seconds'],
   'peak_worker_RSS_KiB':worker['peak_host_RSS_KiB'],'peak_parent_RSS_KiB':final['parent_peak_RSS_KiB'],
   'per_call_seconds':{'min':min(r['wall_seconds'] for r in rows),'median':statistics.median(r['wall_seconds'] for r in rows),
      'max':max(r['wall_seconds'] for r in rows),'sum':sum(r['wall_seconds'] for r in rows)},
   'threads':8,'GPU_seconds':0},
 'call_receipts':rows,'interpretation':'Qualifies this CPU runtime on this202atom paired fixture and preserved optimizer method; no full107 campaign, broad throughput prediction or biological validation.',
 'new_calls_in_collection':0}
with (OUT/'SCOUT_RESULT.json').open('x') as f:json.dump(summary,f,indent=2);f.write('\n')
print(json.dumps({k:summary[k] for k in ['status','proposals','cost']},indent=2))
