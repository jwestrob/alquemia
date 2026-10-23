#!/usr/bin/env python3
"""Bind only a terminal fit to the unchanged saved report; never score or execute it."""
import argparse,copy,hashlib,importlib.util,json,shlex
from pathlib import Path
HERE=Path(__file__).resolve().parent
REPORT=HERE.parent/'REPORT.py';TEMPLATE=HERE.parent/'CONFIG_UNBOUND.json'
REPORT_SHA='d1989966c4d14dfc433b5b50d21f603dd9702b85b1aa3a8c8b5b603e7282a5e1'
TEMPLATE_SHA='ea1c8455231c3b59e3a793c67f4d5c51b882feba06c640e9e04f52c258499c07'
PYTHON='/home/jwestrob/.pyenv/versions/3.12.0/bin/python'

def require(ok,message):
    if not ok:raise ValueError(message)

def pin(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def terminal_artifacts(fit):
    """Inspect terminal metadata first; never open vector artifacts in this function."""
    require(fit.get('stage')=='run' and fit.get('status') in ('bounded_joint_response_fit_complete','failed_or_unavailable') and isinstance(fit.get('cost'),dict),'Actual terminal run required; started, preparation and nonterminal results rejected')
    initial=None;selected=None
    if fit.get('initial_key'):
        row=fit['vectors'][fit['initial_key']];require(row['status']=='complete','Initial artifact is not complete');initial=row['artifact']
    if fit['status']=='bounded_joint_response_fit_complete':
        require(fit.get('error') is None and fit.get('candidate_frozen') is True,'Completed fit lacks frozen selection')
        row=fit['vectors'][fit['selected_key']];selected=fit['selected']
        require(row['status']=='complete' and selected==row['artifact'],'Selected artifact differs from terminal selection')
        require(initial is not None,'Completed fit lacks initial artifact')
    # Failed terminal fits deliberately bind no selected candidate, even if a current/previous point exists.
    return initial,selected


def bind(fit_path,out,report_out,regression=None,scheduler=None):
    require(pin(REPORT)['sha256']==REPORT_SHA and pin(TEMPLATE)['sha256']==TEMPLATE_SHA,'Frozen report/template changed')
    spec=importlib.util.spec_from_file_location('terminal_report_binding',REPORT);r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
    out=Path(out).resolve();report_out=Path(report_out).resolve()
    require(not out.exists(),'Binding output must be new');require(not report_out.exists(),'Report output must be absent')
    require(report_out!=out and report_out not in out.parents,'Report output cannot be binding directory or its ancestor')
    cfg=json.loads(TEMPLATE.read_text());require(all(cfg[k] is None for k in ('fit_result','initial_artifact','selected_artifact','regression_result','scheduler_accounting')),'Template is already bound')
    fit_pin=pin(fit_path);fit=r.read(fit_pin);initial,selected=terminal_artifacts(fit)
    require(fit['implementation']==cfg['driver'],'Terminal fit driver differs')
    dc=r.read(fit['configuration']);ec=r.read(dc['evaluator_config'])
    require(fit['limits']==dc['limits']==dict(outer_steps=24,new_vectors=672,native_programs=24000,workers=12),'Frozen fit limits differ')
    require(dc['evaluator']==cfg['evaluator'] and dc['parameter_domain']==cfg['domain'],'Driver source/domain differs')
    require(all(ec[k]==cfg[k] for k in ('objective','old_objective','body_helper','targets','development_data')),'Evaluator source/target joins differ')
    require(fit['scientific_fingerprint']==r.digest({k:v for k,v in dc.items() if k not in ('root_release','preparation_result')}),'Scientific fingerprint differs')
    for q in (initial,selected):
        if q is not None:r.verify(q)  # bytes only; the existing report owns saved objective replay.
    cfg.update(fit_result=fit_pin,initial_artifact=initial,selected_artifact=selected,
        regression_result=None if regression is None else pin(regression),scheduler_accounting=None if scheduler is None else pin(scheduler))
    if cfg['regression_result']:
        reg=r.read(cfg['regression_result']);require(reg['fit_result']==fit_pin and reg['selected']==selected,'Optional regression candidate/fit differs')
    if cfg['scheduler_accounting']:r.read(cfg['scheduler_accounting'])
    for q in cfg.values():
        if isinstance(q,dict) and set(q)=={'path','sha256'}:r.verify(q)
    require(pin(fit_path)==fit_pin,'Terminal result changed during binding')
    out.mkdir(parents=True,exist_ok=False);config=out/'CONFIG_BOUND.json';r.save(config,cfg)
    env={k:'1' for k in (*r.THREADS,'BLIS_NUM_THREADS')};env['PYTHONDONTWRITEBYTECODE']='1'
    argv=[PYTHON,str(REPORT),'--config',str(config),'--out',str(report_out)]
    shell='env '+' '.join(shlex.quote(k+'='+v) for k,v in env.items())+' '+shlex.join(argv)
    receipt=dict(status='terminal_report_bound_not_executed',implementation=pin(__file__),report_implementation=pin(REPORT),configuration_template=pin(TEMPLATE),configuration=pin(config),fit_result=fit_pin,
        initial_artifact=initial,selected_artifact=selected,fit_status=fit['status'],report_output=str(report_out),thread_environment=env,argv=argv,shell_command=shell,new_native_calls=0,report_executed=False,scoring_executed=False,vector_artifact_contents_read=False)
    r.save(out/'BINDING.json',receipt);(out/'RUN_REPORT.sh').write_text('#!/bin/bash\nset -euo pipefail\nexec '+shell+'\n')
    print(json.dumps(dict(binding=pin(out/'BINDING.json'),configuration=pin(config),argv=argv,thread_environment=env),indent=2));return receipt


def main():
    p=argparse.ArgumentParser();p.add_argument('--fit-result',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--report-out',type=Path,required=True)
    p.add_argument('--regression-result',type=Path);p.add_argument('--scheduler-accounting',type=Path);a=p.parse_args()
    bind(a.fit_result,a.out,a.report_out,a.regression_result,a.scheduler_accounting)
if __name__=='__main__':main()
