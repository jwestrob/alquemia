#!/usr/bin/env python3
"""Read saved old/rebased maps on Slurm; no geometry edits or model calls."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/cap_diagnostic'
SOURCE=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_v1'
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def require(v,s):
    if not v:raise ValueError(s)
def verify(p):
    require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')

def prepare(output):
    c=read(SOURCE/'CONFIG.json');old=read(c['old_config']);worker=read(old['worker_manifest']);actual=read(SOURCE/'preparation/PREPARATION.json')
    rows=[]
    for s in c['sources']:
        r=s['row'];sid=r['source_case_id'];status,=[x for x in actual['rows'] if x['source_case_id']==sid]
        rows.append({'source_case_id':sid,'prior_status':status['status'],'prior_reason':status.get('reason'),
            'old_mapping':r['maps']['Ca'],'rebased_mapping':pin(SOURCE/'preparation'/sid/'mapping.json'),
            'before_xyz':r['matrix']['Ca']['amide_flip']['original_xyz'],'after_xyz':r['matrix']['Ca']['amide_flip']['xyz'],
            'discrete_admission':r['geometry_admission']})
    require(len(rows)==10,'Expected all ten real sources')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('DIAGNOSE.py','RUN.sbatch','PLAN.md'):shutil.copyfile(HERE/name,out/name)
    m={'protocol':'saved_Asn_cap_rebase_roundoff_diagnostic_v1','implementation':pin(out/'DIAGNOSE.py'),
       'wrapper':pin(out/'RUN.sbatch'),'plan':pin(out/'PLAN.md'),'source_config':pin(SOURCE/'CONFIG.json'),
       'source_preparation':pin(SOURCE/'preparation/PREPARATION.json'),'source_execution':pin(SOURCE/'preparation/EXECUTION.json'),
       'python':c['python'],'worker_manifest':old['worker_manifest'],'source_snapshot':worker['source_snapshot'],
       'source_modules':{k:worker['source_modules'][k] for k in ('mace_site_kinematics.py','affordable_response.py','affordable_common.py')},
       'rows':rows,'existing_q0_replay_tolerance_A':1e-12,'new_energy_calls':0,'output':str(out/'execution')}
    write(out/'CONFIG.json',m);print(json.dumps(pin(out/'CONFIG.json')))

def analyze(row,np,Kinematics):
    old=read(row['old_mapping'])['context'];new=read(row['rebased_mapping'])['context'];ad=read(row['discrete_admission']);caps=ad['cap_indices']
    old_lines=verify(row['before_xyz']).read_text().splitlines();new_lines=verify(row['after_xyz']).read_text().splitlines()
    before=np.asarray([[float(v) for v in x.split()[1:]] for x in old_lines[2:]]);after=np.asarray([[float(v) for v in x.split()[1:]] for x in new_lines[2:]])
    require(before.shape==after.shape==(204,3) and np.isfinite(before).all() and np.isfinite(after).all(),'Invalid saved coordinate inventory')
    ko,kn=Kinematics(old),Kinematics(new);require(ko.modes==kn.modes and len(caps)==14,'Changed modes or cap inventory')
    po,co,jo,jco=ko.evaluate(np.zeros(len(ko.modes)));pn,cn,jn,jcn=kn.evaluate(np.zeros(len(kn.modes)))
    for a in (po,co,jo,jco,pn,cn,jn,jcn):require(np.isfinite(a).all(),'Nonfinite kinematic output')
    caprows=[]
    for i in caps:
        lo,=[x for x in old['core_links'] if x[0]==i and x[1]=='cap'];ln,=[x for x in new['core_links'] if x[0]==i and x[1]=='cap'];anchors=lo[2:4]
        caprows.append({'qm_index':i,'old_link':lo,'new_link':ln,'links_exact':lo==ln,
            'raw_before_line':old_lines[i+2],'raw_after_line':new_lines[i+2],'raw_lines_exact':old_lines[i+2]==new_lines[i+2],
            'raw_before_A':before[i].tolist(),'raw_after_A':after[i].tolist(),'raw_coordinates_exact':bool(np.array_equal(before[i],after[i])),
            'old_core_input_A':old['core_positions_A'][i],'new_core_input_A':new['core_positions_A'][i],
            'core_input_exact':old['core_positions_A'][i]==new['core_positions_A'][i],
            'old_q0_A':co[i].tolist(),'new_q0_A':cn[i].tolist(),'q0_exact_between_maps':bool(np.array_equal(co[i],cn[i])),
            'old_q0_minus_raw_A':(co[i]-before[i]).tolist(),'new_q0_minus_raw_A':(cn[i]-after[i]).tolist(),
            'new_minus_old_q0_A':(cn[i]-co[i]).tolist(),'old_anchor_positions_A':ko.positions[anchors].tolist(),
            'new_anchor_positions_A':kn.positions[anchors].tolist(),'anchors_exact':bool(np.array_equal(ko.positions[anchors],kn.positions[anchors])),
            'cap_jacobian_max_change':float(np.max(np.abs(jcn[:,i]-jco[:,i])))})
    return {**row,'status':'diagnosed','cap_count':len(caps),'caps':caprows,
        'raw_caps_exact':bool(np.array_equal(before[caps],after[caps])),
        'old_q0_caps_equal_raw':bool(np.array_equal(co[caps],before[caps])),
        'new_q0_caps_equal_raw':bool(np.array_equal(cn[caps],after[caps])),
        'old_new_q0_caps_exact':bool(np.array_equal(co[caps],cn[caps])),
        'old_q0_cap_max_residual_A':float(np.max(np.abs(co[caps]-before[caps]))),
        'new_q0_cap_max_residual_A':float(np.max(np.abs(cn[caps]-after[caps]))),
        'old_new_q0_cap_max_change_A':float(np.max(np.abs(cn[caps]-co[caps]))),
        'old_full_q0_max_residual_A':float(np.max(np.abs(co-before))),
        'new_full_q0_max_residual_A':float(np.max(np.abs(cn-after))),
        'cap_links_exact':all(v['links_exact'] for v in caprows),'cap_anchors_exact':all(v['anchors_exact'] for v in caprows),
        'raw_input_caps_exact':all(v['core_input_exact'] for v in caprows),
        'frozen_exact_predicate_reproduces_status':bool(np.array_equal(cn[caps],before[caps]))==(row['prior_status']=='prepared')}

def execute(path,sha):
    job=os.environ.get('SLURM_JOB_ID');require(job and job.isdigit(),'Real Slurm required before numerical imports')
    raw=subprocess.check_output(['scontrol','show','job','-o',job],text=True);fields=dict(x.split('=',1) for x in raw.split() if '=' in x)
    require(fields['JobState']=='RUNNING' and fields['NumCPUs']=='1' and fields['NumNodes']=='1','One running allocated CPU required')
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread required')
    mp={'path':str(Path(path).resolve()),'sha256':sha};m=read(mp);require(m['implementation']==pin(__file__) and m['python']==pin(sys.executable),'Implementation/interpreter mismatch')
    for k in ('wrapper','plan','source_config','source_preparation','source_execution','worker_manifest'):verify(m[k])
    for p in m['source_modules'].values():verify(p)
    out=Path(m['output']);out.mkdir(exist_ok=False);start=time.monotonic();rows=[]
    ex={'status':'failed','config':mp,'job_id':job,'host':socket.gethostname(),'scontrol':raw,'new_energy_calls':0}
    try:
        sys.path.insert(0,m['source_snapshot']);import numpy as np
        from mace_site_kinematics import Kinematics
        require(pin(sys.modules[Kinematics.__module__].__file__)==m['source_modules']['mace_site_kinematics.py'],'Wrong imported kinematics')
        for r in m['rows']:
            try:rows.append(analyze(r,np,Kinematics))
            except Exception as e:rows.append({**r,'status':'unavailable','reason':str(e),'traceback':traceback.format_exc()})
        write(out/'RESULT.json',{'config':mp,'declared_sources':10,'rows':rows,'new_energy_calls':0,
            'scope':'Saved q0/cap algebra only; no repaired inputs, changed admission or model calls'})
        ex['status']='complete' if all(r['status']=='diagnosed' for r in rows) else 'partial';ex['result']=pin(out/'RESULT.json')
    except BaseException:ex['traceback']=traceback.format_exc();raise
    finally:ex['wall_seconds']=time.monotonic()-start;write(out/'EXECUTION.json',ex)

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);q=s.add_parser('prepare');q.add_argument('--output',required=True)
    q=s.add_parser('execute');q.add_argument('--config',required=True);q.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.action=='prepare':prepare(a.output)
    else:execute(a.config,a.config_sha)
