#!/usr/bin/env python3
"""Saved-only comparison using the original exact four-cell/twelve-check scorer."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import time
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('radial_binding',HERE/'BIND.py');B=importlib.util.module_from_spec(s);s.loader.exec_module(B)

def collect(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();B.require(all(os.environ.get(k)=='1' for k in B.THREADS),'One saved collection thread')
    B.require(cfg['implementation']==B.pin(__file__),'Frozen collector source')
    b,ctx=B.bound_context(cfg['bound']);m,inp=ctx[:2]
    ref=B.module(m['sources']['REFERENCE'],'radial_ref_collect');model=B.module(m['sources']['MODEL'],'radial_model_collect')
    old=B.module(inp['old_collection'],'radial_old_compare');score=B.module(inp['scorer'],'radial_fixed_score');body=ctx[9].h['body_helper']
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    rp=ref.collect(cfg['bound']['path'],m['reference_manifest']['path'],B.verify(cfg['reference_release']),B.verify(cfg['reference_result_inventory']),out/'REFERENCE_ADMISSION.json');refs=B.read(rp)
    B.require(refs['binding']['bound']==cfg['bound'] and refs['binding']['selected']==b['selected'] and refs['binding']['parameters']==b['parameters'] and refs['binding']['root_acceptance']==b['root_acceptance'],'Reference/model selected candidate join')
    raw=None;mp=None
    try:
        mp=B.pin(cfg['model_result_path']);raw=B.read(mp);admitted=model.admit_result(cfg['bound'],raw)
    except Exception as exc:admitted=old.unavailable_model(repr(exc),b,mp)
    ap=B.save(out/'MODEL_ADMISSION.json',admitted)
    scored=old.compare(b,refs,admitted,score,body)
    status='validation_unavailable' if scored['available_checks']<12 else ('validation_passed' if scored['all_validation_checks_passed'] else 'validation_failed')
    cost={k:raw.get(k) for k in ('job_id','native_program_invocations','new_FD_energies','elapsed_seconds','parent_CPU_seconds','aggregate_child_user_CPU_seconds','aggregate_child_system_CPU_seconds','per_program_child_CPU_may_overlap_not_summed')} if isinstance(raw,dict) else {}
    r=dict(status=status,implementation=B.pin(__file__),configuration=cfg,bound=cfg['bound'],candidate_result=b['candidate_result'],selected_review=b['candidate_review'],selected=b['selected'],parameters=b['parameters'],root_acceptance=b['root_acceptance'],
           declared_cells=4,score=scored,reference_admission=refs,reference_admission_file=rp,model_admission=admitted,model_admission_file=ap,raw_model_result=mp,model_actual_cost=cost,
           new_native_calls=0,local_threads=1,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,physical_accuracy_qualified=False,automatic_fitting_release=False,production_released=False)
    B.save(out/'RESULT.json',r);(out/'REPORT.md').write_text(old.report(r));return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    r=collect(json.loads(Path(a.config).read_text()),a.out);print(json.dumps({'status':r['status'],'new_native_calls':0}))
