#!/usr/bin/env python3
"""Artificial source controls only. No scientific vectors, proposal or native APIs."""
import argparse, ast, copy, hashlib, importlib.util, json, os, tempfile, time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def mod(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    t,c=time.monotonic(),time.process_time();r=mod(HERE/'REVIEW.py','artificial_reviewer');j=mod(HERE.parent/'JOURNAL.py','artificial_journal');checks=[]
    def ok(value,name):
        if not value:raise AssertionError(name)
        checks.append(name)
    def rejects(fun,name):
        try:fun()
        except (ValueError,j.JournalError):checks.append(name);return
        raise AssertionError(name)
    try:
        ok(all(os.environ.get(k)=='1' for k in r.THREADS),'one_thread')
        ok(r.actual_choice(1.,[{'q':.5,'all25':False},None,None,None])==(None,1.),'physical_failure_never_selected')
        ok(r.actual_choice(1.,[{'q':.5,'all25':True,'old_ceiling':False},None,None,None])==(0,.5),'old_ceiling_descriptive')
        ok(r.actual_choice(1.,[{'q':1.-.5e-10,'all25':True},None,None,None])==(None,1.),'current_wins_epsilon_tie')
        ok(r.actual_choice(1.,[{'q':.5,'all25':True},{'q':.5-.5e-10,'all25':True},None,None])==(0,.5),'earlier_fraction_wins_tie')
        ok(r.actual_choice(1.,[{'q':.6,'all25':True},{'q':.5,'all25':True},None,None])==(1,.5),'strict_later_improvement')
        ledger={'native_program_invocations':0,'native_returncodes_observed':0,'vector_keys':[],'requests':[],'stopped':False,'limits':{'new_parameter_vectors':1,'new_native_programs':1}}
        request={'vector_key':'vector','status':'started','cells':{'cell':{'status':'not_attempted'}}}
        attempt={'number':1,'cache_key':'cache','cell_id':'cell','program':'analyze','folder':'artificial','status':'invoke_boundary_entered'}
        completed=dict(attempt,status='returned_success',returncode_observed=True)
        events=[('initialize',{'ledger':ledger,'authority':{},'policy':j.POLICY}),('request_start',{'index':0,'request':request}),
                ('attempt_intent',{'attempt':attempt,'request_index':0,'program_inputs':{}}),
                ('attempt_return',{'attempt':completed,'receipt_delta':[{'returncode':0}],'program_inputs':{}}),
                ('request_cell',{'index':0,'cell_id':'cell','state':{'status':'complete'},'evidence':{}}),
                ('request_update',{'index':0,'updates':{'status':'complete','objective':0.}})]
        def write(path,items,suffix=b''):
            previous='0'*64;lines=[]
            for n,(kind,data) in enumerate(items,1):
                body={'sequence':n,'previous_sha256':previous,'event':{'kind':kind,'data':data}};previous=j.sha(body);lines.append(j.encode(dict(body,sha256=previous))+b'\n')
            path.write_bytes(b''.join(lines)+suffix)
        with tempfile.TemporaryDirectory(prefix='physical_review_artificial_') as td:
            p=Path(td)/'events';write(p,events);actual=r.stream_journal(p,j);old=j.replay(p)
            ok(actual[0]==old['ledger'] and actual[1]['sequence']==6 and actual[1]['last_event_sha256']==old['last_event_sha256'],'one_pass_lossless_journal')
            ok(actual[1]['attempt_status_counts']=={'returned_success':1} and len(actual[3])==len(actual[4])==1,'intent_return_counts')
            write(p,events[:3],b'{partial');q=r.stream_journal(p,j)
            ok(q[1]['partial_suffix_bytes']==8 and q[1]['unresolved_attempt_numbers']==[1] and q[0]['native_program_invocations']==1,'crash_intent_budget_and_partial_preserved')
            write(p,events+[('attempt_return',events[3][1])]);rejects(lambda:r.stream_journal(p,j),'duplicate_return_rejected')
            write(p,events[:3]+[('attempt_intent',events[2][1])]);rejects(lambda:r.stream_journal(p,j),'duplicate_attempt_rejected')
            other=dict(attempt,number=2,cache_key='other')
            write(p,events[:3]+[('attempt_intent',{'attempt':other,'request_index':0,'program_inputs':{}})]);rejects(lambda:r.stream_journal(p,j),'attempt_cap_rejected')
            write(p,events);raw=p.read_bytes();p.write_bytes(raw.replace(b'"returncode":0',b'"returncode":1'));rejects(lambda:r.stream_journal(p,j),'committed_checksum_corruption_rejected')
        oldsrc=HERE.parent.parent/'independent_actual_review_v1/REVIEW.py'
        def function(path,name):return next(x for x in ast.parse(path.read_text()).body if isinstance(x,ast.FunctionDef) and x.name==name)
        ok(ast.dump(function(HERE/'REVIEW.py','independent25'))==ast.dump(function(oldsrc,'independent25')),'original25_independent_arithmetic_unchanged')
        tree=ast.parse((HERE/'REVIEW.py').read_text());calls=[x for x in ast.walk(tree) if isinstance(x,ast.Call)]
        ok(not any(isinstance(x.func,ast.Attribute) and (x.func.attr in ('propose','evaluate','invoke','DefaultSolver') or (x.func.attr=='solve' and ast.unparse(x.func)!='np.linalg.solve')) for x in calls),'no_solver_evaluator_native_call_sites')
        result={'status':'passed_artificial_terminal_review_controls','checks':checks,'error':None}
    except Exception as exc:result={'status':'failed_artificial_terminal_review_controls','checks':checks,'error':repr(exc)}
    result.update(source=r.pin(HERE/'REVIEW.py'),controls=r.pin(__file__),journal=r.pin(HERE.parent/'JOURNAL.py'),cost={'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-c,'threads':1},new_native_calls=0,new_optimizer_calls=0,scientific_vector_reads=0)
    r.save(args.out/'RESULT.json',result);print(json.dumps({'status':result['status'],'checks':len(checks),'error':result['error'],'result':r.pin(args.out/'RESULT.json')}));return int(result['error'] is not None)
if __name__=='__main__':raise SystemExit(main())
