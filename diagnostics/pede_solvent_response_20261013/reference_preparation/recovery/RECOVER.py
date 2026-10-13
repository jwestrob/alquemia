#!/usr/bin/env python3
"""One-CPU, saved-artifact parser regression and collection; no executable calls."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import traceback
from COLLECT_V2 import actual_trah_threshold,gradient_criterion,check,collect,original,pin,write,TOLG,TRAH

def rejected(text):
    try: actual_trah_threshold(text)
    except ValueError: return True
    return False

def main():
    p=argparse.ArgumentParser();p.add_argument('config',type=Path);p.add_argument('--sha',required=True);a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Actual Slurm required')
    check({'path':str(a.config.resolve()),'sha256':a.sha});c=json.loads(a.config.read_text())
    for record in c['code'].values(): check(record)
    if c['code']['RECOVER.py']!=pin(__file__): raise ValueError('Wrong recovery executable')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False)
    status={'status':'started','config':pin(a.config),'code':c['code'],'job_id':os.environ['SLURM_JOB_ID'],
            'started_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'new_molecular_calls':0}
    try:
        for record in c['pinned_inputs']:check(record)
        mp=check(c['manifest']);m=json.loads(mp.read_text());module=original(m)
        pf=json.loads(check(c['preflight']).read_text());raw=check(c['output_stdout']).read_text();old=check(c['older_TRAH_fixture']).read_text()
        threshold=actual_trah_threshold(raw);old_threshold=actual_trah_threshold(old)
        macro=module.trah_macro(raw);old_macro=module.trah_macro(old)
        tests={
            'current_summary_omits_orbital_gradient':'Last Orbital Gradient' not in raw,
            'current_printed_threshold':threshold['threshold']==1e-5,
            'current_actual_macro':bool(macro) and module.numeric(macro[-1][2])==5.663e-7,
            'old_summary_present':'Last Orbital Gradient' in old,
            'old_printed_threshold':old_threshold['threshold']==1e-5,
            'old_actual_macro':bool(old_macro) and module.numeric(old_macro[-1][2])==5.528534e-7,
            'missing_threshold_rejected':rejected(re.sub(TOLG,'',re.sub(TRAH,'',raw,flags=re.M),flags=re.M)),
            'conflicting_threshold_rejected':rejected(re.sub(TRAH,' Converg. threshold (grad. norm) .... 1.000e-06',raw,flags=re.M)),
            'above_threshold_rejected':not gradient_criterion(raw,threshold['threshold']*2)['passed'],
            'negative_norm_rejected':not gradient_criterion(raw,-1e-7)['passed']}
        baseline=module.collect(m,mp,pf);result=collect(module,m,mp,pf)
        before=next(r for r in baseline['cells'] if r['task_id']=='sample1_La_water')
        after=next(r for r in result['cells'] if r['task_id']=='sample1_La_water')
        tests['actual_v1_bug_reproduced']=before['checks']['final_TRAH_gradient'] is False and before['status']=='failed_collection'
        tests['all_other_checks_preserved']=all(before['checks'][k]==v for k,v in after['checks'].items() if k!='final_TRAH_gradient')
        tests['all_other_cell_fields_preserved']=all(after[k]==v for k,v in before.items() if k not in ('checks','status'))
        tests['current_recovered']=after['status']=='complete_diagnostic'
        tests['remaining7_not_completed']=all(r['status']=='not_completed' for r in result['cells'] if r['task_id']!='sample1_La_water')
        evidence={'tests':tests,'all_passed':all(tests.values()),'current_threshold':threshold,'older_threshold':old_threshold,
                  'negative_fixture_scope':'in-memory parser-only text modifications, never scientific inputs or saved output'}
        write(out/'REGRESSION.json',evidence);write(out/'COLLECTION.json',result)
        if not all(tests.values()): raise RuntimeError('Real-artifact regression failed')
        status['status']='complete';status['regression']=pin(out/'REGRESSION.json');status['collection']=pin(out/'COLLECTION.json')
    except BaseException as e:
        status['status']='failed';status['error']=f'{type(e).__name__}: {e}';status['traceback']=traceback.format_exc();raise
    finally:
        status['finished_utc']=dt.datetime.now(dt.timezone.utc).isoformat();write(out/'EXECUTION.json',status)

if __name__=='__main__':main()
