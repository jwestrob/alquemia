#!/usr/bin/env python3
"""Real-input guard integration and owned-stop/failure-receipt checks; no chemistry."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json, shutil, signal, threading
import guarded_native_gfn2_dispatch as dispatch
import orca_parameter_runtime_qualification as b

def main():
    out=b.W/'guarded_dispatch_validation_v2';out.mkdir(exist_ok=False)
    def prepare(jid,name=None):
        source=b.W/'results'/jid;r=b.read(source/'receipt.json');d=out/(name or jid);d.mkdir()
        for f in ['endpoint.inp','core.xyz']:shutil.copyfile(source/f,d/f)
        b.write(d/'expected.json',r['job']['expected'])
        return SimpleNamespace(input=str(d/'endpoint.inp'),xyz=str(d/'core.xyz'),
            parameters=str(b.ROOT/'params/native_gfn2_parameters_v1/corrected_la.xtb.json'),
            policy=str(b.ROOT/'params/native_gfn2_parameters_v1/la_shell_three_arrays_v1.json'),
            approval=str(b.ROOT/'diagnostics/orca_parameter_audit_20261012/parameter_guard/runtime_registry_ca_la_corrected_20261012.json'),
            runtime_binding=str(b.W/'RUNTIME_BINDING.json'),observed_runtime=str(b.W/'OBSERVED_RUNTIME_FINAL.json'),
            protocol_id='guarded-dispatch-final-real-input-validation-v2',
            medium='vacuum' if r['job']['medium']=='vacuum' else 'ALPB(Water)',
            expected=str(d/'expected.json'),source=[str(source/'receipt.json')],
            charge=r['job']['expected']['charge'],multiplicity=r['job']['expected']['multiplicity'],
            timeout_seconds=240,execute=False)
    records=[]
    with patch.object(b,'upstream',side_effect=AssertionError('Molecular executor forbidden in validation')):
        for jid in ['Ca_vacuum_corrected_file','Ca_water_corrected_file','La_vacuum_corrected','La_water_corrected']:
            args=prepare(jid);receipt=dispatch.dispatch(args)
            assert receipt['numerical_calls']==0 and receipt['scientific_result_accepted'] is None
            assert receipt['guard_preflight']['parameter_preflight_passed']
            records.append({'case':jid,'passed':True,'receipt':b.pin(Path(args.input).parent/'dispatch_receipt.json')})
        for kind in ['wrong_parameter_file','wrong_expected_state']:
            args=prepare('Ca_vacuum_corrected_file',kind);args.execute=True
            if kind=='wrong_parameter_file':args.parameters=str(b.W/'inputs/La_vacuum_original/control.xtb.json')
            else:
                e=b.read(args.expected);e['charge']+=1;Path(args.expected).write_text(json.dumps(e))
            try:dispatch.dispatch(args)
            except ValueError as error:records.append({'case':kind,'passed':True,'rejected':True,'error':str(error)})
            else:raise AssertionError('invalid binding passed')
    args=prepare('Ca_vacuum_corrected_file','simulated_parser_failure_only');args.execute=True
    event=threading.Event();sig_before={s:signal.getsignal(s) for s in [signal.SIGINT,signal.SIGTERM,signal.SIGHUP]}
    def executor_stub(*a,**kw):
        handler=signal.getsignal(signal.SIGTERM);assert callable(handler)
        handler(signal.SIGTERM,None);assert event.is_set()
        return {'exit_code':0,'stop_reason':None,'error':None,'simulation_only':True,'actual_numerical_calls':0}
    module=SimpleNamespace(STOP=event,bounded_process=executor_stub,
        parse_output=lambda *a,**k:(_ for _ in ()).throw(ValueError('intentional simulated parser failure')))
    with patch.object(b,'upstream',return_value=module),patch.object(dispatch.guard,'postflight',return_value={'parameter_postflight_passed':True,'simulation_only':True}):
        try:dispatch.dispatch(args)
        except ValueError as error:assert str(error)=='intentional simulated parser failure'
        else:raise AssertionError('parser fault did not surface')
    assert all(signal.getsignal(s)==h for s,h in sig_before.items())
    p=Path(args.input).parent/'dispatch_receipt.json';r=b.read(p)
    assert r['scientific_result_accepted'] is False and r['interrupt_signals']==['SIGTERM']
    assert 'intentional simulated parser failure' in r['error'] and r['artifacts']
    records.append({'case':'simulated_parser_failure_only','passed':True,'owned_stop_flag_set':event.is_set(),
        'original_signal_handlers_restored':True,'failure_receipt_preserved':b.pin(p),
        'simulation_only':'Executor and postflight mocked solely to exercise parser-failure cleanup; no scientific value or molecular process created.'})
    result={'implementation':b.pin(dispatch.__file__),'validation_implementation':b.pin(__file__),
        'checks':records,'all_passed':all(r['passed'] for r in records),'actual_numerical_calls':0,
        'supersedes':'Earlier v1 dry-run/rejection proofs predate reviewed dispatcher cleanup repair; retained as historical evidence only.'}
    b.write(b.REPORT/'GUARDED_DISPATCH_VALIDATION_V2.json',result)
    print(json.dumps({'checks':len(records),'all_passed':result['all_passed'],'actual_numerical_calls':0,'implementation_sha256':result['implementation']['sha256']}))

if __name__=='__main__':main()
