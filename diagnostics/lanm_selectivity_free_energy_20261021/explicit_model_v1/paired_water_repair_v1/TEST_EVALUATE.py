#!/usr/bin/env python3
"""Saved-only evaluator admission and failure controls. Never runs a native tool."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('paired_evaluator',HERE/'EVALUATE.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def main():
    tick=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    cfg=m.default_config();m.save(a.out/'CONFIG.json',cfg);e=m.Evaluator(cfg,a.out/'cache');records=[];initial=list(m.INITIAL)
    # Any accidental native invocation is an immediate test failure.
    def forbidden(*args,**kwargs):raise AssertionError('Saved validation attempted a native invocation')
    e.native.invoke=forbidden
    with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(e.evaluate,[initial]*8))
    for model in results:
        for gid,value in model.items():
            assert value['energy']==e.data['rows'][gid]['baseline']['energy'] and value['forces']==e.data['rows'][gid]['baseline']['forces']
            assert value['provenance']['reused_baseline'] and value['provenance']['requested_vector']==initial
    snap=e.snapshot();assert snap['distinct_parameter_vectors']==1 and snap['native_program_invocations']==0 and len(snap['requests'])==8
    assert all(set(x['cells'])==set(m.CELLS) and x['status']=='complete' for x in snap['requests'])
    records.append({'case':'four_worker_baseline_replay_and_cache','requests':8,'native_calls':0,'passed':True})
    for gid in m.CELLS:
        basekey,baserec=e.render(initial,gid);row=e.data['rows'][gid];assert basekey==m.verify(row['input']['key']).read_text()
        records.append({'case':'original_native_baseline_'+gid,'evidence':e.baseline_review[gid],'passed':True})
        own=0 if gid.endswith('_La') else 2;other=2-own
        v=initial[:];v[own]*=1.01;v[own+1]+=.01;key,record=e.render(v,gid)
        assert key.startswith(basekey.rstrip()+'\n\n') and len(key[len(basekey.rstrip()):].strip().splitlines())==3
        assert record['nonwater_parameters_unchanged'] and record['resolved_pair']['epsilon_kcal_mol']==m.EPSILON[row['input']['metal']]
        a1=e.cache_identity(v,gid)[0];v2=v[:];v2[other]*=.99;v2[other+1]+=.02;assert e.cache_identity(v2,gid)[0]==a1
        v3=v[:];v3[own]=math.nextafter(v3[own],math.inf);assert e.cache_identity(v3,gid)[0]!=a1
        assert e.cache_identity(initial,gid)[0]!=a1
        records.append({'case':'active_pair_only_render_and_exact_effective_cache_'+gid,'passed':True})
        for program,keyname in (('analyze','analyze_output'),('testgrad','gradient_output')):
            e._program_admission(m.verify(row['baseline']['native'][keyname]),program)
    # Exact four source rows only; no reserved or donor fixture is accessible.
    try:e.render(initial,'hydrated_acetamide_central_La')
    except ValueError:records.append({'case':'outside_panel_rejected','passed':True})
    else:raise AssertionError('Unbounded cell admitted')
    for bad in ([*initial,float('nan')],[float('nan'),*initial[1:]],[1.,*initial[1:]],initial[:3]):
        try:e.checked_values(bad)
        except ValueError:pass
        else:raise AssertionError('Invalid parameter admitted')
    records.append({'case':'invalid_shape_nonfinite_bounds_rejected','passed':True})
    native_path=Path(e.data['rows']['water_central_La']['baseline']['native']['analyze_output']['path']);text=native_path.read_text()
    for name,modified in (('polarization_failure',text+'\nInduced dipoles not converged\n'),('nonfinite_marker',text+'\nNaN\n'),('missing_energy',text.replace('Total Potential Energy','Missing Energy')),('missing_components',text.replace('Energy Component Breakdown :','Missing breakdown'))):
        f=a.out/(name+'.synthetic.txt');f.write_text(modified)
        try:e._program_admission(f,'analyze')
        except ValueError:records.append({'case':name,'passed':True,'synthetic':True})
        else:raise AssertionError(name)
    # Gate failures are tested without forging an allocation or creating a native record.
    v=initial[:];v[0]*=1.01
    try:e.evaluate(v)
    except ValueError:pass
    else:raise AssertionError('Unreleased native calculation allowed')
    snap=e.snapshot();assert snap['native_program_invocations']==0 and snap['stopped'];assert set(snap['requests'][-1]['cells'])==set(m.CELLS)
    assert not list((a.out/'cache').glob('cells/*/RESULT.json'))
    records.append({'case':'no_release_or_allocation_stops_new_native_calls_and_retains_denominator','passed':True})
    e2=m.Evaluator(cfg,a.out/'cache');assert e2.stop_event.is_set() and e2.snapshot()['native_program_invocations']==0
    records.append({'case':'persistent_failure_not_automatically_retried','passed':True})
    m.save(a.out/'SNAPSHOT.json',snap)
    result={'status':'passed','source':m.pin(__file__),'implementation':m.pin(HERE/'EVALUATE.py'),'config':m.pin(a.out/'CONFIG.json'),'config_sha256':m.config_fingerprint(cfg),
        'records':records,'family':e.family,'new_native_program_invocations':0,'saved_baseline_native_rows_replayed':4,'maximum_local_worker_threads':4,
        'numerical_environment':{k:os.environ.get(k) for k in m.THREADS},'command':sys.argv,'elapsed_seconds':time.monotonic()-tick,'reserved_validation_outcomes_opened':False}
    m.save(a.out/'RESULT.json',result);print(json.dumps({'status':'passed','controls':len(records),'result':m.pin(a.out/'RESULT.json'),'native_calls':0}))

if __name__=='__main__':main()
