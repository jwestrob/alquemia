#!/usr/bin/env python3
"""Manufactured proposal-count control only; no campaign/module/solver reads."""
import argparse, ast, hashlib, json, time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def pin(p):p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def require(ok,msg):
    if not ok:raise ValueError(msg)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,type=Path);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    t,c=time.monotonic(),time.process_time();checks=[];source=HERE/'REVIEW_v2.py'
    try:
        tree=ast.parse(source.read_text());old=ast.parse((HERE/'REVIEW.py').read_text());fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='check_proposal_counts')
        ns={'require':require};exec(compile(ast.Module(body=[fn],type_ignores=[]),str(source),'exec'),ns);f=ns['check_proposal_counts']
        def expect(args,reported,name):
            value=f(*args,reported);checks.append(name);return value
        def reject(args,reported,name):
            try:f(*args,reported)
            except ValueError:checks.append(name);return
            raise AssertionError(name)
        def report(logical,real,imported):return {'logical_proposal_calls':logical,'real_proposal_calls':real,'imported_proposal_calls':imported}
        expect((1,0,1,0),report(1,0,1),'first_import_only')
        expect((8,7,1,0),report(8,7,1),'full_eight_step_solver_cap')
        value=expect((2,0,1,1),report(2,0,1),'later_explicit_zero_call_failure_counts_as_logical')
        require(value['zero_call_failed_proposals']==1,'Explicit failure retained')
        expect((3,1,1,1),report(3,1,1),'mixed_real_imported_and_failed_counts')
        reject((2,0,1,0),report(2,0,1),'missing_logical_record_category_rejected')
        reject((2,0,1,1),report(1,0,1),'false_reported_logical_count_rejected')
        reject((8,8,0,0),report(8,8,0),'eight_new_solves_rejected')
        reject((2,0,2,0),report(2,0,2),'duplicate_import_rejected')
        reject((9,7,1,1),report(9,7,1),'nine_logical_records_rejected')
        oldf={x.name:x for x in old.body if isinstance(x,ast.FunctionDef)};newf={x.name:x for x in tree.body if isinstance(x,ast.FunctionDef)}
        unchanged=set(oldf)-{'run'}
        require(all(ast.dump(oldf[k])==ast.dump(newf[k]) for k in unchanged),'Only run counting and new helper changed');checks.append('all_existing_other_function_ASTs_unchanged')
        status='passed_proposal_count_source_controls';error=None
    except Exception as exc:status='failed_proposal_count_source_controls';error=repr(exc)
    result={'status':status,'source':pin(source),'previous_source':pin(HERE/'REVIEW.py'),'controls':pin(__file__),'checks':checks,'error':error,'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-c,'threads':1},'new_native_calls':0,'new_optimizer_calls':0,'scientific_vector_reads':0}
    p=args.out/'RESULT.json';p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':status,'checks':len(checks),'error':error,'result':pin(p)}));return int(error is not None)
if __name__=='__main__':raise SystemExit(main())
