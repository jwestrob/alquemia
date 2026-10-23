#!/usr/bin/env python3
"""Six consumed postselection trimers; reuse the frozen evaluator's native API."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
IDS=tuple(f'aquo{n}_central_{metal}__{pair}__metal_two_water'
          for n in (8,9) for metal in ('La','Dy') for pair in (('near','far') if n==8 else ('near',)))

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):
    p=Path(p);t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');os.replace(t,p)
def module(p):
    s=importlib.util.spec_from_file_location('consumed_regression_evaluator',verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def cycles(terms,values):
    result={}
    for cid,target in terms.items():
        coeff=target['model_terms'];missing=[k for k in coeff if k not in values or values[k].get('status','complete')!='complete']
        row=dict(status='unavailable' if missing else 'complete',model_terms=coeff,missing=missing,
                 reference_electronic_E3_kcal_mol=target['reference_electronic_E3_kcal_mol'],baseline_E3_kcal_mol=target['model_E3_kcal_mol'],
                 reference_mu3_Debye=target['reference_mu3_Debye'],baseline_mu3_Debye=target['model_mu3_Debye'],
                 original_reference_print_bound_Debye_per_component=target['reference_print_bound_Debye_per_component'],
                 original_model_print_bound_Debye_per_component=target['model_print_bound_Debye_per_component'])
        if not missing:
            energy=math.fsum(c*values[k]['energy'] for k,c in coeff.items());mu=[math.fsum(c*values[k]['mu_Debye'][j] for k,c in coeff.items()) for j in range(3)]
            require(math.isfinite(energy) and all(math.isfinite(v) for v in mu),'Nonfinite consumed cycle')
            row.update(candidate_E3_kcal_mol=energy,energy_error_kcal_mol=energy-target['reference_electronic_E3_kcal_mol'],candidate_mu3_Debye=mu,
                       dipole_error_Debye=[v-t for v,t in zip(mu,target['reference_mu3_Debye'])])
        result[cid]=row
    return result

def prepare(out):
    inputs=read(pin(HERE/'SOURCE_INPUTS.json'));ev=module(inputs['evaluator']);ready=read(inputs['evaluator_ready']);review=read(inputs['evaluator_review'])
    require(ready['implementation']==review['implementation']==inputs['evaluator'] and review['source_review_passed'],'Frozen evaluator source review')
    cfg=read(ready['configuration_template']);em=read(cfg['manifest']);src=em['sources'];domain=read(src['domain']);inv=read(src['fixture_inventory'])
    target=read(src['target_inventory']);rm=read(src['baseline_manifest']);mp=read(src['moment_preparation']);full=read(src['baseline_result']);actual=read(src['baseline_review'])
    require(actual['status']=='passed_baseline_matrix_inventory_review' and actual['result']==src['baseline_result'] and actual['manifest']==src['baseline_manifest'],'Admitted54 baseline lineage')
    require(rm['inventory']==mp['inventory']==src['fixture_inventory'] and mp['targets']==src['target_inventory'],'Literal inventory/frame lineage')
    terms=target['groups']['trimer_energy_and_dipoles']['cycles'];require(len(terms)==8,'Eight consumed regression cycles')
    all_trimers={next(k for k,c in t['model_terms'].items() if c==1 and k.endswith('__metal_two_water')) for t in terms.values()}
    require(all_trimers-set(em['rows'])==set(IDS),'Exact six excluded trimers')
    rows={};baseline={}
    for cid in IDS:
        row=copy.deepcopy(rm['rows'][cid]);support=mp['rows'][cid];old=full['rows'][cid]
        require(old['status']=='complete' and row['source_input']==support['input'] and row['types']==support['types'] and row['positions_A']==support['positions_A'],'Literal admitted trimer/frame join')
        row.update(metal=inv['rows'][cid]['metal'],needs_export=False,needs_gradient=False,permanent_lab_Debye=support['permanent_lab_Debye'],energy_replay_kcal_mol=None)
        row['fixed_POLPAIR']=domain['fixed_native_POLPAIR'][row['metal']]
        for p in row['source_input'].values():verify(p)
        moments=old['moments']
        for k in ('analyze_output','induced_output'):verify(moments[k])
        baseline[cid]=dict(status='complete',inherited_native=True,producer_result=src['baseline_result'],producer_review=src['baseline_review'],
                           producer_cell_id=cid,energy=moments['energy_kcal_mol'],mu_Debye=[support['permanent_lab_Debye'][j]+moments['induced_sum_Debye'][j] for j in range(3)])
        rows[cid]=row
    m=dict(status='prepared_consumed_regression_unbound',implementation=pin(__file__),specification=pin(HERE/'SPEC.md'),inputs=inputs,evaluator_template=ready['configuration_template'],evaluator_manifest=cfg['manifest'],
           sources=src,rows=rows,baseline_cells=baseline,cycles=terms,units=em['units'],initial_vector=em['initial_vector'],maximum_native_calls=6,workers=1,new_native_calls=0)
    save(out/'MANIFEST.json',m)
    save(out/'CONFIG_RUN_TEMPLATE.json',dict(implementation=pin(__file__),specification=m['specification'],manifest=pin(out/'MANIFEST.json'),evaluator_config=None,fit_result=None,root_release=None))
    return m

def context(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Regression source/specification')
    m=read(cfg['manifest']);require(m['implementation']==cfg['implementation'] and m['specification']==cfg['specification'] and set(m['rows'])==set(IDS),'Frozen six-row preparation')
    ev=module(m['inputs']['evaluator']);ec=read(cfg['evaluator_config']);em=read(ec['manifest']);template=read(m['evaluator_template'])
    require(ec['implementation']==m['inputs']['evaluator'] and ec['manifest']==m['evaluator_manifest'] and em['sources']==m['sources'],'Evaluator input lineage')
    for key in ev.SOURCES:require(ec[key]==template[key],'Changed evaluator source '+key);verify(ec[key])
    fr=read(cfg['fit_result']);fc=read(fr['configuration']);selected=read(fr['selected'])
    require(fr['status']=='bounded_joint_response_fit_complete' and fr['candidate_frozen'] is True and fr['error'] is None,'Terminal frozen candidate required')
    require(fc['evaluator']==m['inputs']['evaluator'] and fc['evaluator_config']==cfg['evaluator_config'],'Actual fit/evaluator configuration join')
    require(fr['selected']==fr['vectors'][fr['selected_key']]['artifact'] and selected['status']=='complete','Exact selected artifact')
    values=ev.checked(selected['parameters'],read(ec['domain']))
    require(set(selected['cells'])==set(em['rows']) and all(r['status']=='complete' for r in selected['cells'].values()),'Complete frozen46 candidate')
    require(selected['score']['objective']==fr['selected_objective'] and selected['score']['original_EF_checks_passed']==fr['original_EF_checks_passed'],'Selected objective/gate receipt')
    q=read(ec['qualification_review']);build=read(ec['build_receipt'])
    require(q['status']=='passed_local_response_native_qualification' and all(q[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Completed native numerical qualification')
    require(all(q[k]==ec[k] for k in ('source_manifest','build_receipt','executables')) and build['status']=='complete' and build['manifest']==ec['source_manifest'],'Exact qualified native build')
    for k,p in ec['executables'].items():require(build['executables'][k]==p,'Build executable join');verify(p)
    for p in (ec['source_manifest'],q['implementation'],q['result']):verify(p)
    h={k:ev.module(ec[k],'regression_'+k) for k in ev.HELPERS};domain=read(ec['domain'])
    proxy=SimpleNamespace(h=h,m=em,config=ec)
    proxy.render=lambda vals,cid,exports=None:ev.render(vals,em['rows'][cid],domain,h['resolver'],exports)
    needed={cid for t in m['cycles'].values() for cid in t['model_terms']} & set(selected['cells'])
    for cid in sorted(needed):
        r=selected['cells'][cid];row=em['rows'][cid]
        if r['inherited_native']:
            require(ev.effective(values,row['metal'],domain)==ev.effective(m['initial_vector'],row['metal'],domain),'Inherited selected cell effective parameters')
            require(all(r[k]==v for k,v in em['baseline_cells'][cid].items()),'Inherited selected baseline equality')
        else:
            directory=Path(r['receipts'][0]['cwd']).parent;ev.Evaluator._saved_inputs(proxy,r,row,values,directory)
            parsed=ev.Evaluator._parse(proxy,r,row);require(all(r[k]==v for k,v in parsed.items()),'Selected dependency raw replay')
    return m,ev,ec,em,fr,selected,values,domain,proxy

def run(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();m,ev,ec,em,fr,selected,values,domain,proxy=context(cfg)
    release=read(cfg['root_release'])
    require(release['status']=='released_consumed_trimer_regression' and all(release[k]==cfg[k] for k in ('implementation','manifest','fit_result','evaluator_config')) and release['selected']==fr['selected'],'Separate exact frozen-candidate release')
    require(release['maximum_native_calls']==6 and release['workers']==1 and release['reserved_outcomes_authorized'] is False,'Finite consumed-only release')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and os.environ.get('SLURM_NTASKS','1')=='1' and os.environ.get('SLURM_CPUS_PER_TASK')=='1' and all(os.environ.get(k)=='1' for k in THREADS),'One serial single-thread Slurm CPU')
    result=dict(status='running',implementation=pin(__file__),configuration=cfg,manifest=cfg['manifest'],fit_result=cfg['fit_result'],selected=fr['selected'],parameters=values,
                rows={cid:dict(status='not_attempted') for cid in IDS},cycles={},native_program_invocations=0,original_EF_checks_passed=fr['original_EF_checks_passed'],
                original_EF=selected['score']['old_EF'],error=None,job_id=os.environ['SLURM_JOB_ID'],fit=False,physical_accuracy_qualified=False,reserved_outcomes_opened=False)
    vals={cid:dict(status='complete',energy=r['energy_kcal_mol'],mu_Debye=r['mu_Debye']) for cid,r in em['static'].items()};vals.update(selected['cells'])
    def persist():
        result['cycles']=cycles(m['cycles'],dict(vals,**result['rows']));save(out/'RESULT.json',result)
    try:
        persist()
        for cid in IDS:
            row=m['rows'][cid]
            if ev.effective(values,row['metal'],domain)==ev.effective(m['initial_vector'],row['metal'],domain):
                result['rows'][cid]=copy.deepcopy(m['baseline_cells'][cid]);persist();continue
            folder=out/'cells'/cid;folder.mkdir(parents=True,exist_ok=False);key,meta=ev.render(values,row,domain,proxy.h['resolver'],exports=False)
            for k,name in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['source_input'][k]),folder/name)
            (folder/'system.key').write_text(key);inp={k:pin(folder/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            rec=dict(status='started',cell_id=cid,inherited_native=False,render=meta,program_inputs={'analyze':inp},receipts=[]);result['rows'][cid]=rec
            require(result['native_program_invocations']<6,'Six-call hard ceiling');result['native_program_invocations']+=1;persist()
            try:
                proxy.h['native_helper'].invoke([ec['executables']['analyze']['path'],'system.xyz','EM'],folder,folder/'analyze.out',rec['receipts'])
                rec.update(ev.Evaluator._parse(proxy,rec,row));rec['status']='complete'
            except Exception as exc:rec.update(status='failed',error=repr(exc));raise
            finally:save(folder/'RESULT.json',rec);persist()
        require(all(r['status']=='complete' for r in result['rows'].values()) and all(r['status']=='complete' for r in result['cycles'].values()),'Complete six/eight denominator')
        result['status']='consumed_trimer_regression_complete'
    except Exception as exc:result.update(status='consumed_trimer_regression_failed',error=repr(exc))
    result['cost']=dict(wall_seconds=time.monotonic()-tick,parent_CPU_seconds=time.process_time()-cpu,allocated_CPUs=1);persist();return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));p.add_argument('--config');p.add_argument('--out',required=True);a=p.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
    r=prepare(out) if a.action=='prepare' else run(json.loads(Path(a.config).read_text()),out)
    print(json.dumps({'status':r['status'],'output':str(out)}))
    raise SystemExit(0 if r['status'] in ('prepared_consumed_regression_unbound','consumed_trimer_regression_complete') else 1)
