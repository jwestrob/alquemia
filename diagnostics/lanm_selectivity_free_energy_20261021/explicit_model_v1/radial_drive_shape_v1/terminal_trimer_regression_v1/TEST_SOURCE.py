#!/usr/bin/env python3
"""Small manufactured accounting/failure controls; no native/scientific replay."""
import argparse,copy,json,math,os,tempfile,time,types
from pathlib import Path
import RUN as r

def main(out):
    start,cpu=time.monotonic(),time.process_time();out.mkdir(parents=True,exist_ok=False);checks=[]
    def check(name,ok):r.require(ok,name);checks.append({'name':name,'passed':True})
    def rejects(name,fun):
        try:fun()
        except (ValueError,RuntimeError):check(name,True)
        else:raise AssertionError(name)
    s=r.read(r.pin(Path(__file__).parent/'SOURCE_INPUTS.json'))
    inv=r.read(r.read(r.pin(Path(__file__).parent/'SCIENTIFIC_READY.json'))['inventory'])
    r.validate_partition(inv);check('exact_14_22_6_inventory_partition',True)
    bad=copy.deepcopy(inv);bad['new_analyze_order']=list(reversed(r.IDS));rejects('reject_reordered_or_extra_calls',lambda:r.validate_partition(bad))
    bad=copy.deepcopy(inv);bad['rows'][r.IDS[0]]['needs_export']=True;rejects('reject_unreleased_matrix_branch',lambda:r.validate_partition(bad))
    bad=copy.deepcopy(inv);bad['historical_six_row_import_authorized']=True;rejects('reject_old_six_row_cache_import',lambda:r.validate_partition(bad))
    scalar=r.extract(s['parent_adapter'],('cycles','comparisons'),{'require':r.require,'math':math})
    terms={};vals={};previous={}
    for n in (8,9):
        for metal in ('La','Dy'):
            for pair in ('near','far'):
                cid=f'aquo{n}_central_{metal}__{pair}';offset=1 if metal=='La' else 3
                terms[cid]={'model_terms':{'trimer_'+cid:1,'fragment':-1},'reference_electronic_E3_kcal_mol':0.,'model_E3_kcal_mol':-1.,'reference_mu3_Debye':[0.,0.,0.],'model_mu3_Debye':[0.,0.,0.],'reference_print_bound_Debye_per_component':.001,'model_print_bound_Debye_per_component':.001}
                vals['trimer_'+cid]={'status':'complete','energy':offset+2.,'mu_Debye':[offset+1.,-offset-2.,offset+3.]}
                previous[cid]={'candidate_E3_kcal_mol':.5,'candidate_mu3_Debye':[.1,.2,.3]}
    vals['fragment']={'status':'complete','energy':2.,'mu_Debye':[1.,-2.,3.]}
    cycles=scalar['cycles'](terms,vals);comparison=r.relabel(scalar['comparisons'](cycles,previous))
    check('eight_signed_manufactured_cycles',len(cycles)==8 and cycles['aquo8_central_Dy__near']['candidate_mu3_Debye']==[3.,-3.,3.])
    pair=comparison['paired_Dy_minus_La']['aquo8__near']['models']['frozen_radial14_terminal5d992']
    check('four_signed_Dy_minus_La_contrasts',len(comparison['paired_Dy_minus_La'])==4 and pair['Dy_minus_La_E3_kcal_mol']==2. and pair['Dy_minus_La_mu3_Debye']==[2.,-2.,2.])
    check('explicit_feca12_radial14_baseline_labels',set(comparison['cycles']['aquo8_central_La__near']['models'])=={'feca12_physical_gate_terminal','frozen_radial14_terminal5d992','original_pre_response_baseline'})
    missing=copy.deepcopy(vals);missing['trimer_aquo8_central_La__near']['status']='failed';c=scalar['cycles'](terms,missing)
    check('partial_failure_not_a_complete_cycle',c['aquo8_central_La__near']['status']=='unavailable' and 'candidate_E3_kcal_mol' not in c['aquo8_central_La__near'])
    with tempfile.TemporaryDirectory() as td:
        tmp=Path(td);(tmp/'source.xyz').write_text('manufactured xyz\n');(tmp/'source.prm').write_text('manufactured prm\n')
        source={'xyz':r.pin(tmp/'source.xyz'),'parameters':r.pin(tmp/'source.prm')}
        def exercise(name,fail_invoke=0,fail_parse=0,initial=0):
            dest=tmp/name;dest.mkdir();calls=[];snapshots=[]
            result={'rows':{k:{'status':'not_attempted'} for k in r.IDS},'native_program_invocations':initial,'executor':{'mode':'manufactured'}}
            def persist():snapshots.append(copy.deepcopy(result))
            def invoke(argv,cwd,output,receipts):
                calls.append(argv);r.require(snapshots[-1]['native_program_invocations']==initial+len(calls),'Intent is saved before native entry')
                output.write_text('manufactured only\n');receipts.append({'argv':argv,'returncode':1 if len(calls)==fail_invoke else 0,'output':r.pin(output)})
                if len(calls)==fail_invoke:raise RuntimeError('manufactured native failure')
            def parse(record,row):
                if len(calls)==fail_parse:raise RuntimeError('manufactured parser failure')
                return {'energy':1.,'mu_Debye':[0.,0.,0.],'radial_response_readback':{'analyze':{'native_solve_blocks':2}}}
            proxy=types.SimpleNamespace(render=lambda *a,**k:('mock radial key\n',{'radial_drive':.1}),_parse=parse,h={'native_helper':types.SimpleNamespace(invoke=invoke)})
            ctx=types.SimpleNamespace(inv={'rows':{cid:{'source_input':source} for cid in r.IDS}},proxy=proxy,values=[0.]*14,s={'analyze':{'path':'MOCK_NEVER_EXECUTED'}})
            error=None
            try:r.execute_six(ctx,result,dest,persist)
            except (ValueError,RuntimeError) as exc:error=str(exc)
            return result,calls,error,snapshots
        result,calls,error,snaps=exercise('success');check('exact_six_serial_native_slots',len(calls)==6 and error is None and all(x['status']=='complete' for x in result['rows'].values()))
        check('twelve_induction_blocks_distinct_from_zero_optimizer',sum(x['radial_response_readback']['analyze']['native_solve_blocks'] for x in result['rows'].values())==12)
        result,calls,error,snaps=exercise('invoke_failure',fail_invoke=2);check('stop_after_second_native_failure_no_retry',len(calls)==2 and error is not None and result['rows'][r.IDS[1]]['status']=='failed' and all(result['rows'][k]['status']=='not_attempted' for k in r.IDS[2:]))
        check('failed_native_receipt_preserved',result['rows'][r.IDS[1]]['receipts'][0]['returncode']==1)
        result,calls,error,snaps=exercise('parse_failure',fail_parse=2);check('successful_native_output_survives_parser_failure',len(calls)==2 and result['rows'][r.IDS[1]]['status']=='failed' and result['rows'][r.IDS[1]]['receipts'][0]['returncode']==0)
        result,calls,error,snaps=exercise('budget',initial=6);check('hard_six_cap_no_seventh_invoke',not calls and error is not None)
        cfg={'implementation':{'path':'unbound','sha256':'0'},'candidate':None,'inventory':None,'root_release':None}
        bad=r.run(cfg,tmp/'startup_failure',preflight=True);check('startup_failure_has_zero_call_terminal_receipt',bad['status']=='radial14_consumed_trimer_regression_failed' and bad['native_program_invocations']==0 and (tmp/'startup_failure/RESULT.json').is_file())
    r.save(out/'RESULT.json',{'status':'passed_radial14_terminal_trimer_source_controls','implementation':r.pin(__file__),'runner':r.pin(r.__file__),'checks':checks,'checks_passed':len(checks),'actual_native_calls':0,'actual_optimizer_calls':0,'actual_scientific_replays':0,'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}})
    print(json.dumps({'checks_passed':len(checks),'actual_native_calls':0}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();main(Path(a.out))
