#!/usr/bin/env python3
"""Stage/preflight a fresh native-GFN2 endpoint; --execute opts into one bounded call.

Every path requires the shared parameter guard. There is no unguarded/default
fallback. Scientific SCF/state/population admission remains a separate check.
"""
from pathlib import Path
import argparse, ctypes, json, os, re, signal, traceback
import native_gfn2_parameters as guard
import orca_parameter_runtime_qualification as runtime


def dispatch(args):
    inp=Path(args.input).resolve();directory=inp.parent;xyz=Path(args.xyz).resolve()
    if inp.name!='endpoint.inp':raise ValueError('bounded parser requires endpoint.inp')
    if xyz.parent!=directory:raise ValueError('XYZ must be copied beside the fresh input')
    if any(directory.glob('*.gbw')) or any((directory/n).exists() for n in ['endpoint.out','dispatch_receipt.json']):
        raise ValueError('fresh endpoint directory required; restart/output files exist')
    text=inp.read_text()
    ranks=re.findall(r'^\s*nprocs\s+(\d+)\s*$',text,re.I|re.M)
    if ranks!=['1'] or re.search(r'\bPAL\d+\b',text,re.I):raise ValueError('exact serial nprocs1 required')
    if not re.search(r'\bNoAutostart\b',text,re.I):raise ValueError('explicit NoAutostart required')
    expected=runtime.read(args.expected)
    required=['charge','multiplicity','valence_electrons','atoms','basis_functions']
    if not all(k in expected for k in required):raise ValueError('complete independent expected-state record required')
    if expected['charge']!=args.charge or expected['multiplicity']!=args.multiplicity:
        raise ValueError('expected state differs from declared binding')
    sources=[*args.source,args.expected,__file__,runtime.__file__,str(runtime.W/'source/run_native_controls.py')]
    binding=guard.build_binding(input_path=inp,xyz=xyz,parameters=args.parameters,
        runtime_binding=args.runtime_binding,observed_runtime=args.observed_runtime,source_inputs=sources,
        protocol_id=args.protocol_id,charge=args.charge,multiplicity=args.multiplicity,
        medium=args.medium,parameter_filename='control.xtb.json')
    if any(binding['environment'].get(k)!='1' for k in runtime.THREAD_KEYS):
        raise ValueError('runtime approval does not enforce single-thread children')
    saved={k:os.environ.get(k) for k in binding['environment']}
    handlers={}
    receipt={'schema':'guarded_native_gfn2_dispatch.v1','created_utc':runtime.now(),
             'execute_requested':args.execute,'dispatch_implementation':guard.pin(__file__),
             'scientific_result_accepted':None,'numerical_calls':0,'expected':expected}
    try:
        for k,v in binding['environment'].items():
            if v is None:os.environ.pop(k,None)
            else:os.environ[k]=v
        staged=guard.stage(binding,args.parameters,args.policy,approval=guard.pin(args.approval),mode='research')
        receipt['guard_stage']=staged
        manifest=staged['manifest']['path'];before=guard.preflight(manifest)
        receipt['guard_preflight']=before
        if args.execute:
            libc=ctypes.CDLL(None,use_errno=True)
            if libc.prctl(36,1,0,0,0)!=0:raise OSError('subreaper unavailable')
            module=runtime.upstream();observed={};max_threads=0
            module.STOP.clear()
            def interrupted(sig,frame):
                receipt.setdefault('interrupt_signals',[]).append(signal.Signals(sig).name)
                module.STOP.set()
            for sig in [signal.SIGINT,signal.SIGTERM,signal.SIGHUP]:
                handlers[sig]=signal.getsignal(sig)
                signal.signal(sig,interrupted)
            def monitor(pgid):
                nonlocal max_threads
                rss,members=module.group_snapshot(pgid)
                for member in members:
                    p=Path('/proc')/str(member['pid'])/'status'
                    try:
                        n=int(re.search(r'^Threads:\s+(\d+)',p.read_text(),re.M)[1])
                        max_threads=max(max_threads,n);observed[member['pid']]=max(n,observed.get(member['pid'],0))
                    except FileNotFoundError:continue
                    if n>1:raise RuntimeError('owned process exceeded one numerical thread')
                return rss,members
            env=os.environ.copy();env['ORCA_TMPDIR']=str(directory)
            for key in ['OMP_PLACES','OMP_PROC_BIND']:env.pop(key,None)
            process=module.bounded_process([binding['executable']['path'],inp.name,'--bind-to none'],
                directory,env,args.timeout_seconds,6442450944,monitor)
            receipt.update(numerical_calls=1,process=process,max_process_threads=max_threads,observed_threads=observed)
            try:
                receipt['guard_postflight']=guard.postflight(manifest,directory/'endpoint.out',process['exit_code'])
            except guard.GuardError as error:
                receipt['guard_postflight']={'parameter_postflight_passed':False,'error':str(error)}
            output=module.parse_output(directory/'endpoint.out',expected);prop=runtime.properties(directory)
            q=prop.get('charges_e')
            checks={'normal':prop['normal_termination'],'converged':prop['converged'],
                    'state':all(prop.get(k)==expected[k] for k in ['charge','multiplicity','valence_electrons','atoms']),
                    'all_charges':q is not None and len(q)==expected['atoms'] and prop['charge_indices']==list(range(expected['atoms'])),
                    'charge_sum':bool(q) and abs(sum(q)-expected['charge'])<=1e-5,
                    'energy_agreement':prop.get('energy_hartree') is not None and output.get('energy_hartree') is not None and abs(prop['energy_hartree']-output['energy_hartree'])<=1e-9}
            receipt.update(output=output,property=prop,property_checks=checks)
            receipt['scientific_result_accepted']=bool(receipt['guard_postflight']['parameter_postflight_passed'] and
                process['exit_code']==0 and not process['stop_reason'] and not process['error'] and output['passed'] and all(checks.values()))
        return receipt
    except BaseException as error:
        receipt.update(scientific_result_accepted=False,error=repr(error),traceback=traceback.format_exc())
        raise
    finally:
        for sig,handler in handlers.items():signal.signal(sig,handler)
        receipt['finished_utc']=runtime.now()
        receipt['artifacts']={p.name:guard.pin(p) for p in sorted(directory.iterdir())
                              if p.is_file() and p.name!='dispatch_receipt.json'}
        runtime.write(directory/'dispatch_receipt.json',receipt)
        for k,v in saved.items():
            if v is None:os.environ.pop(k,None)
            else:os.environ[k]=v


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['input','xyz','parameters','policy','approval','runtime-binding','observed-runtime','protocol-id','medium','expected']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--source',action='append',required=True)
    p.add_argument('--charge',type=int,required=True);p.add_argument('--multiplicity',type=int,required=True)
    p.add_argument('--timeout-seconds',type=int,default=240,choices=range(1,241),metavar='1..240')
    p.add_argument('--execute',action='store_true')
    args=p.parse_args();r=dispatch(args)
    print(json.dumps({'receipt':str(Path(args.input).parent/'dispatch_receipt.json'),
                      'guard_preflight':r['guard_preflight']['parameter_preflight_passed'],
                      'numerical_calls':r['numerical_calls'],'scientific_result_accepted':r['scientific_result_accepted']}))
    if r['execute_requested'] and not r['scientific_result_accepted']:raise SystemExit(2)


if __name__=='__main__':main()
