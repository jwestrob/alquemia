#!/usr/bin/env python3
"""Narrow saved interface controls; no radial outcome binding or native calls."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('radial_transfer_interface',HERE/'RUN.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r)


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);checks=[]
    def check(ok,label):
        assert ok,label;checks.append({'case':label,'passed':True})
    def rejects(call,label):
        try:call()
        except (ValueError,AssertionError):check(True,label)
        else:raise AssertionError(label)
    data=r.loaded(r.verify(r.PINS['pins']['prepared32']));cfg=data['config']
    training=r.read(r.PINS['pins']['radial_training'])
    radial=HERE.parent/'paired_water_radial_repair_v1';pcfg=json.loads((radial/'CONFIG_PREPARE.json').read_text())
    prep_path=Path(r.PINS['pins']['radial_training']['path']).parents[1]/'preparation_v1/RESULT.json'
    prep=json.loads(prep_path.read_text());ecfg=r.read(pcfg['evaluator_config'])
    check(prep['status']=='prepared' and prep['new_native_program_invocations']==0 and prep['implementation']==r.PINS['pins']['radial_source'],'Use only zero-call radial preparation, not execution outcome')
    controls=r.read(r.PINS['pins']['previous_saved_controls'])
    check(controls['status']=='passed' and len(controls['records'])==49 and controls['new_native_calls']==0,'Reuse49 original saved parser/arithmetic controls')
    check(all(getattr(r,n)is getattr(r.old,n)for n in ('loaded','gradient','program_admission','replay_native')),'Original native helpers are identical objects')
    check(len(r.WATER)==8 and len(r.d.CELLS)==32 and len(r.d.FD)==10,'Eight water/32 full/10FD denominator')
    check(len(r.d.CENTERS)+2*len(r.d.HYDRATED)==40 and sum(6*data['rows'][k]['geometry']['atoms']for k in r.d.FD)==1380,'40 program and1380 internal FD ceilings')
    check(set(r.MEMBERSHIPS['radial_water_development_8'])==set(r.WATER) and set(r.MEMBERSHIPS['hydrated_donor_transfer_18'])==set(r.d.HYDRATED),'Distinct development and transfer memberships')
    r.family_gate(prep['family'],training,ecfg);check(True,'Start differs from original bound origins without widening')
    bad=copy.deepcopy(prep['family']);bad['parameters'][0]['upper']+=.01
    rejects(lambda:r.family_gate(bad,training,ecfg),'Reject widened bounds')
    bad=copy.deepcopy(prep['family']);bad['parameters'][0]['initial']=training['family']['parameters'][0]['initial']
    rejects(lambda:r.family_gate(bad,training,ecfg),'Reject old bound origin substituted for actual start')
    previous=r.read(r.PINS['pins']['previous_transfer_result']);native=r.module(cfg['native_helper'],'binding_test_native');body=r.module(cfg['body_helper'],'binding_test_body')
    objective=r.module(r.PINS['pins']['radial_objective'],'binding_test_objective')
    check(objective.evaluate(training,prep['initial']['model'],body)==prep['initial']['score'],'Eight-water prepared score replays exactly')
    for gid in r.WATER:
        value=prep['initial']['model'][gid];inp=previous['rows'][gid]['input']['files']
        check(r.selected_authority(gid,value,data['rows'][gid],inp,ecfg,prep_path,previous)=='exact17578_candidate_baseline','Exact17578 authority branch '+gid)
        reuse=r.previous_reuse(gid,inp,previous,data['rows'][gid],native)
        check(reuse is not None and (gid not in r.d.FD or reuse['cartesian_check']['fd_passed']),'Exact previous physical input/FD reuse '+gid)
    gid='water_central_La';value=copy.deepcopy(prep['initial']['model'][gid]);inp=previous['rows'][gid]['input']['files']
    check(value['native']['finite_difference']is None and previous['rows'][gid]['cartesian_check']['fd_passed'],'Analytic model and prior full FD remain separate')
    bad=copy.deepcopy(value);bad['provenance']['native_authority']=data['rows'][gid]['baseline_result']
    rejects(lambda:r.selected_authority(gid,bad,data['rows'][gid],inp,ecfg,prep_path,previous),'Reject original SI authority as17578 candidate baseline')
    badinput=copy.deepcopy(inp);badinput['key']['sha256']='0'*64
    check(r.previous_reuse(gid,badinput,previous,data['rows'][gid],native)is None,'Different key cannot reuse previous completed chemistry')
    # Synthetic authority metadata exercises the new-cache interface without reading the radial execution.
    candidate=out/'synthetic_candidate/RESULT.json';directory=candidate.parent/'cache/cells/synthetic';directory.mkdir(parents=True)
    value['provenance']['reused_baseline']=False
    ar={'status':'complete','cell_id':gid,'cache_key':value['provenance']['cache_key'],'input':inp,
        'cache_identity':{k:ecfg[k]for k in ('implementation','native_helper','resolver','executables','training_result')},
        **{k:value[k]for k in ('energy','forces','native')},
        'receipts':[{'returncode':0,'output':value['native'][k]}for k in ('analyze_output','gradient_output')],
        'synthetic_interface_control_not_native_execution':True}
    ar['cache_identity'].update(cell_id=gid,effective_values=value['provenance']['effective_parameters'])
    r.save(directory/'RESULT.json',ar);value['provenance']['native_authority']=r.pin(directory/'RESULT.json')
    check(r.selected_authority(gid,value,data['rows'][gid],inp,ecfg,candidate,previous)=='new_radial_native_cache','Synthetic new-cache authority branch')
    ar['cache_identity']['implementation']=cfg['evaluator_helper'];r.save(directory/'RESULT.json',ar);value['provenance']['native_authority']=r.pin(directory/'RESULT.json')
    rejects(lambda:r.selected_authority(gid,value,data['rows'][gid],inp,ecfg,candidate,previous),'Reject relabeled old evaluator cache')
    result={'status':'passed','implementation':r.pin(__file__),'run_source':r.pin(HERE/'RUN.py'),'plan':r.PINS['pins']['plan'],
        'prepared_data':r.PINS['pins']['prepared32'],'radial_zero_call_preparation':r.pin(prep_path),
        'previous_saved_controls':r.PINS['pins']['previous_saved_controls'],'records':checks,
        'new_native_calls':0,'radial_execution_outcome_opened':False,'candidate_bound':False,'reserved_validation_opened':False,
        'elapsed_seconds':time.monotonic()-start}
    r.save(out/'RESULT.json',result);print(json.dumps({'status':'passed','checks':len(checks),'result':r.pin(out/'RESULT.json')}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args();main(a.out)
