#!/usr/bin/env python3
"""Admit exact eight consumed water rows and17578 candidate baseline, zero calls."""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3];WS=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
TRANSFER=WS/'explicit_model_v1/paired_water_transfer_v1/run_v1/execution/RESULT.json'
TRANSFER_SHA='702d97ababc232336e3b5405a1326cec4cafeccc9570ce9f514076cf7908e986'
CELLS=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central') for m in ('La','Dy'))


def require(ok,why):
    if not ok:raise ValueError(why)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed evidence '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def prepare(out):
    tick=time.monotonic();require(pin(TRANSFER)['sha256']==TRANSFER_SHA,'Wrong transfer result');transfer=json.loads(TRANSFER.read_text())
    require(transfer['status']=='consumed_transfer_complete' and transfer['stage_A_passed'] and not transfer['global_stop'],'Transfer unavailable')
    bound=read(transfer['bound']);data=read(bound['data']);cfg=data['config'];candidate=read(bound['candidate_result'])
    require(candidate['status']=='bounded_joint_refinement_complete' and bound['parameters']==candidate['selected']['parameters'],'Candidate identity differs')
    helper=module(bound['data_source'],'radial_transfer_data');native=module(cfg['native_helper'],'radial_native_parser')
    body=module(cfg['body_helper'],'radial_rigid_body');reference=module(cfg['reference_helper'],'radial_reference')
    objective=module(pin(HERE/'OBJECTIVE.py'),'radial_objective');rows={};checks=[]
    for gid in CELLS:
        row=copy.deepcopy(data['rows'][gid]);actual=transfer['rows'][gid];require(actual['status']=='complete','Candidate row unavailable')
        for name in ('xyz','key','parameters'):verify(row['input'][name]);verify(actual['input']['files'][name])
        expected,render=helper.render(cfg,row,bound['parameters']);require(verify(actual['input']['files']['key']).read_text()==expected,'Candidate key differs')
        require(actual['input']['files']['xyz']['sha256']==row['input']['xyz']['sha256'] and
                actual['input']['files']['parameters']['sha256']==row['input']['parameters']['sha256'],'Candidate geometry/parameters differ')
        for receipt in actual['receipts']:require(receipt['returncode']==0,'Failed candidate native receipt');verify(receipt['output'])
        if '_central_' in gid:require(actual['cartesian_check']['fd_passed'],'Center Cartesian check failed');verify(actual['cartesian_check']['output'])
        model=actual['model'];parsed=native.parse(verify(model['native']['analyze_output']),verify(model['native']['gradient_output']),dict(row['input'],derivative_test=False))
        require(parsed==model['native'] and parsed['total_energy_kcal_mol']==model['energy'] and parsed['force_kcal_mol_A']==model['forces'],'Candidate raw E/F replay differs')
        rr=read(row['reference_result']);panel=reference.load_panel(verify(rr['panel']));cell=next(c for c in panel['cells'] if c['cell_id']==rr['cell_id'])
        replay=reference.collect_one(panel,cell,Path(rr['output']['path']).parent,rr['returncode'])
        require(replay['status']=='complete' and all(rr[k]==v for k,v in replay.items()),'Reference raw replay differs')
        require(rr['forces_kcal_mol_A']==row['reference']['forces'] and rr['energy_hartree']*reference.HA==row['reference']['energy'],'Reference arrays/units differ')
        require(body.bodies(row)==row['bodies'],'Bodies changed')
        row.update(original_baseline=row['baseline'],original_baseline_result=row['baseline_result'],baseline=copy.deepcopy(model),
            baseline_result=pin(TRANSFER),baseline_input=actual['input']['files'],baseline_parameters=bound['parameters'],baseline_parse_mode='analytic_only')
        rows[gid]=row;checks.append({'cell':gid,'candidate_raw_replay':True,'reference_raw_replay':True,'source_key_render_exact':True})
    result={'status':'radial_training_complete','implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'objective':pin(HERE/'OBJECTIVE.py'),
        'transfer_result':pin(TRANSFER),'transfer_bound':transfer['bound'],'transfer_data':bound['data'],'candidate_result':bound['candidate_result'],
        'source_bindings':{k:cfg[k] for k in ('native_helper','body_helper','resolver','reference_helper','executables')},
        'baseline_parameters':bound['parameters'],'family':candidate['family'],'rows':rows,'complete':8,'new_native_calls':0,'local_threads':1,
        'validation_outcomes_opened':False,'physical_accuracy_qualified':False,'checks':checks}
    model={k:r['baseline'] for k,r in rows.items()};score=objective.evaluate(result,model,body)
    target={k:r['reference'] for k,r in rows.items()};selfscore=objective.evaluate(result,target,body);require(selfscore['objective']==0 and selfscore['development_checks_passed'],'Self-target control failed')
    # A common E shift within a metal cancels every balanced/work target.
    shifted=copy.deepcopy(target)
    for k in CELLS:shifted[k]['energy']+=1 if k.endswith('_La') else 2
    require(objective.evaluate(result,shifted,body)['objective']<1e-18,'Balanced energy-offset cancellation failed')
    internal=copy.deepcopy(target);row=rows['water_central_La'];ids=next(x['indices'] for x in row['bodies'] if x['category']=='water');i,j=ids[:2]
    vector=[row['geometry']['positions_A'][j][k]-row['geometry']['positions_A'][i][k] for k in range(3)]
    for k in range(3):internal['water_central_La']['forces'][i][k]+=vector[k];internal['water_central_La']['forces'][j][k]-=vector[k]
    require(objective.evaluate(result,internal,body)['objective']<1e-20,'Internal force leaked into rigid target')
    result.update(baseline_score=score,saved_controls_passed=True,saved_controls=['eight candidate/reference raw replays','self_target_zero','constant_metal_energy_offset_cancellation','internal_central_pair_rigid_cancellation'],
        elapsed_seconds=time.monotonic()-tick)
    out.mkdir(parents=True,exist_ok=False);(out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json'),'objective':score['objective'],'native_calls':0}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();prepare(a.out)
