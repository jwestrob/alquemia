#!/usr/bin/env python3
"""Finite fraction routing/source delta/final pin controls; no scientific inputs."""
import argparse,ast,hashlib,json,os,time,types
from pathlib import Path
import RUN_v2 as R

def main(out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);checks=[]
    def check(n,ok):R.require(ok,n);checks.append(n)
    def reject(n,fn):
        try:fn()
        except ValueError:check(n,True)
        else:raise AssertionError(n)
    ns={'FRACTIONS':(1.,.5,.25,.125)}
    exec("def certificate(p,x,f,ok):\n if f not in FRACTIONS: raise ValueError('fraction')\n p['seen'].append((f,FRACTIONS,list(x),ok))\n return {'f':f}\n",ns)
    core=types.SimpleNamespace(certify_fraction=ns['certificate']);point=[2.,3.]
    def qualified(c,a,p,x,f,ok):
        raw=c.certify_fraction(p,x,f,ok);aux=c.certify_fraction(p,[.125,.1875],1.,ok);return {'raw':raw,'aux':aux}
    p={'seen':[]};r=R.fixed_fraction_certificate(core,None,p,point,.0625,False,qualified)
    check('raw_private_allowlist_exact1_16',p['seen'][0]==(.0625,(.0625,),point,False))
    check('internal_full_auxiliary_uses_unchanged_original',p['seen'][1]==(1.,(1.,.5,.25,.125),[.125,.1875],False))
    check('original_function_globals_and_identity_unchanged',core.certify_fraction is ns['certificate'] and ns['FRACTIONS']==(1.,.5,.25,.125))
    for f in (1.,.5,.25,.125,.03125):reject('external_fraction_rejected_'+str(f),lambda f=f:R.fixed_fraction_certificate(core,None,{'seen':[]},point,f,False,qualified))
    reject('missing_auxiliary_call_rejected',lambda:R.fixed_fraction_certificate(core,None,{'seen':[]},point,.0625,False,lambda c,a,p,x,f,o:c.certify_fraction(p,x,f,o)))
    def extra(c,a,p,x,f,o):qualified(c,a,p,x,f,o);return c.certify_fraction(p,x,1.,o)
    reject('third_helper_call_rejected',lambda:R.fixed_fraction_certificate(core,None,{'seen':[]},point,.0625,False,extra))
    def wrong(c,a,p,x,f,o):c.certify_fraction(p,x,f,o);return c.certify_fraction(p,x,.5,o)
    reject('nonfull_internal_call_rejected',lambda:R.fixed_fraction_certificate(core,None,{'seen':[]},point,.0625,False,wrong))
    src=Path(R.__file__).read_text();old=(R.HERE/'RUN.py').read_text();tree=ast.parse(src);otree=ast.parse(old)
    newfn={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)};oldfn={n.name:n for n in otree.body if isinstance(n,ast.FunctionDef)}
    check('all_unchanged_functions_AST_exact',all(ast.dump(newfn[k],include_attributes=False)==ast.dump(v,include_attributes=False) for k,v in oldfn.items() if k not in ('metadata','context')))
    helper=ast.get_source_segment(src,newfn['fixed_fraction_certificate'])
    inverse=src.replace('import types\n','').replace("pin(HERE/'RUN_v2.sbatch')","pin(HERE/'RUN.sbatch')").replace(helper+'\n\n\n','')
    replacement="scalar=types.SimpleNamespace(certificate=lambda c,a,p,x,f,ok:fixed_fraction_certificate(c,a,p,x,f,ok,solver.certificate))\n    policy=C.Policy(core,attr,M,scalar,fixed,inv,read(cfg['targets']),{'settings':cfg['settings']},None,cfg['settings'])"
    original="policy=C.Policy(core,attr,M,solver,fixed,inv,read(cfg['targets']),{'settings':cfg['settings']},None,cfg['settings'])"
    check('inverse_exact_raw_source_recovery',inverse.replace(replacement,original)==old)
    # Execute exact nested final persist code with a mock that changes ledger bytes.
    persist=next(n for n in ast.walk(newfn['run']) if isinstance(n,ast.FunctionDef) and n.name=='persist')
    ledger=out/'LEDGER.json';ledger.write_text('before\n');journal=out/'EVENTS.jsonl';journal.write_text('event\n')
    ev=types.SimpleNamespace(cache_dir=out,ledger_path=ledger,ledger=dict(vector_keys=['one'],native_program_invocations=56,native_returncodes_observed=56))
    def materialize(out,result,ev,terminal):ledger.write_text('after authoritative snapshot\n');R.save(out/'RESULT.json',result)
    result={};scope=dict(vars(R),out=out,result=result,ev=ev,P=types.SimpleNamespace(persist_result=materialize));exec(compile(ast.fix_missing_locations(ast.Module(body=[persist],type_ignores=[])),__file__,'exec'),scope);scope['persist'](True)
    check('final_ledger_pin_after_materialization',result['journal_reconciliation']['ledger']==R.pin(ledger) and result['journal_reconciliation']['final_authoritative_snapshot'])
    check('final_saved_result_has_authoritative_pin',json.loads((out/'RESULT.json').read_text())['journal_reconciliation']==result['journal_reconciliation'])
    def failed(out,result,ev,terminal):result['terminal_snapshot_error']='mock mismatch';result['error']='mock mismatch'
    bad={'complete_native_admission':True,'witness_scientific_passed':True,'journal_reconciliation':{}};scope.update(result=bad,P=types.SimpleNamespace(persist_result=failed));scope['persist'](True)
    check('unreconciled_terminal_cannot_be_admitted',bad['complete_native_admission'] is False and bad['witness_scientific_passed'] is False and 'journal_reconciliation' not in bad and bad['status']=='failed_boundary_fraction_witness')
    res=dict(status='passed_fixed_fraction_recovery_controls',implementation=R.pin(__file__),run=R.pin(R.__file__),preserved_run=R.pin(R.HERE/'RUN.py'),checks=checks,checks_passed=len(checks),native_programs=0,solver_calls=0,scientific_vectors_read=0,threads=1,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu));R.save(out/'RESULT.json',res);print(json.dumps(dict(status=res['status'],checks=len(checks))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);main(p.parse_args().out)
