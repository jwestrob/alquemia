#!/usr/bin/env python3
"""Finite terminal saved audit: two witnesses, no solver/native/parser rerun."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace

CHECKS=[]
PINS={}
DATA={}
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')

def require(ok,msg):
    CHECKS.append({'name':msg,'passed':bool(ok)})
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();s=p.stat();key=(str(p),s.st_size,s.st_mtime_ns)
    if key not in PINS:PINS[key]={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    return PINS[key]
def verify(p):require(pin(p['path'])==p,'Immutable '+p['path']);return Path(p['path'])
def read(p):
    key=(p['path'],p['sha256'])
    if key not in DATA:DATA[key]=json.loads(verify(p).read_text())
    return DATA[key]
def token(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def forbidden(*args,**kwargs):raise ValueError('Saved audit forbids solver/native calls')
def near(a,b,msg):
    require(math.isfinite(float(a)) and math.isfinite(float(b)) and abs(a-b)<=1e-10*max(1.,abs(a),abs(b)),msg+' inherited arithmetic')
def compare(a,b,msg):
    if isinstance(a,dict):
        require(isinstance(b,dict) and set(a)==set(b),msg+' keys')
        for k in a:compare(a[k],b[k],msg+'/'+k)
    elif isinstance(a,list):
        require(isinstance(b,list) and len(a)==len(b),msg+' length')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,msg+'/'+str(i))
    elif isinstance(a,bool) or a is None or isinstance(a,str):require(a==b,msg+' exact')
    elif isinstance(a,(int,float)):near(a,b,msg)
    else:raise TypeError(msg)

def independent25(vector,data):
    """Reconstruct original body projections/cycles, preserving failed gates."""
    model=vector['model'];old=vector['score']['old_EF'];gates={}
    for kind in ('La','Dy','paired'):
        for geo in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central'):
            la,dy=data[geo+'_La'],data[geo+'_Dy'];xyz=la['geometry']['positions_A'];residual=[]
            for i in range(len(xyz)):
                residual.append([(model[geo+'_Dy']['forces'][i][j]-model[geo+'_La']['forces'][i][j])-(dy['reference']['forces'][i][j]-la['reference']['forces'][i][j]) if kind=='paired' else model[geo+'_'+kind]['forces'][i][j]-data[geo+'_'+kind]['reference']['forces'][i][j] for j in range(3)])
            raw=[]
            for category in ('ion','water'):
                for body in (b for b in la['bodies'] if b['category']==category):
                    ids=body['indices'];center=[math.fsum(xyz[i][j] for i in ids)/len(ids) for j in range(3)]
                    raw.extend(math.fsum(residual[i][j] for i in ids) for j in range(3))
                    if category=='water':
                        for a,b in ((1,2),(2,0),(0,1)):raw.append(math.fsum((xyz[i][a]-center[a])*residual[i][b]-(xyz[i][b]-center[b])*residual[i][a] for i in ids))
            nwater=(len(xyz)-1)//3;weights=[1/6]*3+[1/(12*nwater)]*(6*nwater)
            rms=math.sqrt(math.fsum(w*x*x for w,x in zip(weights,raw)));maximum=max(abs(x) for x in raw);saved=old['per_fixture'][kind][geo]
            compare(raw,saved['raw_body_components_kcal_mol_A'],'Cartesian-to-body '+kind+'/'+geo)
            near(rms,saved['rigid_RMS_kcal_mol_A'],'Body RMS');near(maximum,saved['maximum_body_component_error_kcal_mol_A'],'Body max')
            passed=saved['rigid_RMS_kcal_mol_A']<=2 and saved['maximum_body_component_error_kcal_mol_A']<=5
            require(saved['passed']==passed,'Unchanged exact force gate')
            gates['force/'+kind+'/'+geo]={'rigid_RMS_kcal_mol_A':saved['rigid_RMS_kcal_mol_A'],'maximum_body_component_error_kcal_mol_A':saved['maximum_body_component_error_kcal_mol_A'],'passed':passed}
    def E(x,g,m):return x[g+'_'+m]['energy']
    refs={k:r['reference'] for k,r in data.items()};errors={}
    for g in ('water_minus','water_central','water_plus','aquo9_central'):
        errors['balanced_exchange/'+g]=(E(model,g,'Dy')-E(model,g,'La')-E(model,'aquo8_central','Dy')+E(model,'aquo8_central','La'))-(E(refs,g,'Dy')-E(refs,g,'La')-E(refs,'aquo8_central','Dy')+E(refs,'aquo8_central','La'))
    for m in ('La','Dy'):
        for sign in ('minus','plus'):errors['individual_radial_work/'+m+'_'+sign]=(E(model,'water_'+sign,m)-E(model,'water_central',m))-(E(refs,'water_'+sign,m)-E(refs,'water_central',m))
    for sign in ('minus','plus'):errors['paired_radial_work/'+sign]=errors['individual_radial_work/Dy_'+sign]-errors['individual_radial_work/La_'+sign]
    for name,error in errors.items():
        group,key=name.split('/');saved=old['energy_controls'][group][key];near(error,saved['error_kcal_mol'],'Independent energy cycle')
        passed=abs(saved['error_kcal_mol'])<=1;require(saved['passed']==passed,'Unchanged exact energy gate')
        gates['energy/'+name]={'error_kcal_mol':saved['error_kcal_mol'],'passed':passed}
    require(len(gates)==25 and all(v['passed'] for v in gates.values())==vector['score']['original_EF_checks_passed'],'All25 flags including failures')
    return gates

def run(c,out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    report={'status':'started','reviewer_implementation':pin(__file__),'review_configuration':c,'new_native_calls':0,'new_solver_calls':0,'new_reference_calls':0,'selected_candidate':None,'model_admitted':False,'error':None}
    try:
        require(c['implementation']==pin(__file__) and str(out)==c['output_directory'],'Exact reviewer/output')
        require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
        r=read(c['result']);require(r['status']!='started','Terminal result only')
        cfg=read(r['configuration']);ready=read(c['runner_ready']);inv=read(cfg['inventory']);sources=read(inv['inputs'])['sources']
        require(r['implementation']==cfg['implementation']==ready['implementation'],'Reviewed producer')
        fp=token({k:v for k,v in cfg.items() if k not in ('root_release','source_review','preparation')})
        require(r['configuration_fingerprint']==fp==ready['configuration_fingerprint'],'Frozen runtime fingerprint')
        require(Path(c['result']['path']).parent==Path(cfg['output_directory']),'Exact result directory')
        runner=module(cfg['implementation'],'saved_two_witness_runner');runner.check_release(cfg)
        require(r['limits']==cfg['limits']==runner.LIMITS and r['solver_calls']==r['reference_calls']==r['heldout_reads']==0,'Finite scope no solver/reference/heldout')
        require(r['selected_candidate'] is None and r['job_id'] is None,'No selected candidate; truthful login')
        plan=read(cfg['plan']);affine=read(plan['affine_result']);center=read(sources['selected'])
        report.update(implementation=r['implementation'],result=c['result'],configuration=r['configuration'],root_release=cfg['root_release'],source_review=cfg['source_review'],center=sources['selected'],affine_result=plan['affine_result'])
        old=module(c['journal_reviewer'],'saved_witness_journal_reviewer');journal=module(cfg['journal_implementation'],'saved_witness_journal')
        successful=r['status']=='completed_two_direction_native_validation' and r['error'] is None
        if 'evaluator' not in r:
            require(not successful and r['native_programs']==0,'Failed-before-evaluator zero attempts')
            report.update(status='reviewed_failed_two_direction_before_evaluator',producer_status=r['status'],producer_error=r['error'],producer_cost=r['cost']);return report
        ledger=read(r['evaluator']['ledger']);lp=Path(r['evaluator']['ledger']['path'])
        replay,summary,initial,intents,returns,cellpins=old.stream_journal(lp.with_name('EVENTS.jsonl'),journal)
        require(replay==ledger,'One journal replay equals full ledger')
        side=read(pin(lp.with_name(lp.name+'.journal.json')))
        require(side['ledger']==r['evaluator']['ledger'] and side['sequence']==summary['sequence'] and side['last_event_sha256']==summary['last_event_sha256'],'Committed ledger sidecar')
        require(initial['authority']=={'configuration_fingerprint':fp,'implementation':cfg['implementation'],'root_release':cfg['root_release']},'Journal exact authority')
        attempts=ledger.get('native_attempts',[])
        require(len(attempts)==ledger['native_program_invocations']==r['native_programs']<=112,'Actual attempt ceiling112')
        require(set(intents)==set(range(1,len(attempts)+1)) and len(ledger['requests'])==len(ledger['vector_keys'])<=2,'Finite intent/request counts')
        require(ledger['native_returncodes_observed']==sum(bool(a.get('returncode_observed')) for a in attempts),'Actual observed return counts')
        for num,returned in returns.items():
            require(num in intents,'Return has durable intent')
            for key in ('number','cache_key','cell_id','program','folder'):require(returned['attempt'][key]==intents[num]['attempt'][key],'Intent/return '+key)
            require(returned['program_inputs']==intents[num]['program_inputs'],'Inputs at invocation/return')
            for p in returned['program_inputs'].values():verify(p)
            for receipt in returned['receipt_delta']:verify(receipt['output'])
        report.update(journal=summary,actual_native_programs=r['native_programs'],actual_attempts=attempts,producer_cost=r['cost'])
        if not successful:
            report.update(status='reviewed_failed_or_incomplete_two_direction_validation',producer_status=r['status'],producer_error=r['error'],remaining_rows=r['rows']);return report
        require(not ledger['stopped'] and not summary['partial_suffix_bytes'] and not summary['unresolved_attempt_numbers'],'No hidden failure')
        require(r['native_programs']==112==inv['potential_new_programs'] and len(returns)==112 and len(r['vectors'])==2,'Completed2/112')
        require(all(a['status']=='returned_success' and a['returncode_observed'] for a in attempts),'All attempted programs returned')
        require(ledger['vector_keys']==[x['key'] for x in inv['rows']] and list(r['rows'])==['B_half','D_half'],'Declared fixed ordering')
        require(all(q['status']=='complete' and len(q['cells'])==46 and all(v['status']=='complete' for v in q['cells'].values()) for q in ledger['requests']),'Both complete46-cell batches')
        parent=read(cfg['parent_configuration']);base=module(cfg['parent_evaluator'],'saved_witness_evaluator')
        require(parent['implementation']==cfg['parent_evaluator'],'Unchanged evaluator')
        for p in parent['executables'].values():verify(p)
        build=read(parent['build_receipt']);qual=read(parent['qualification_review'])
        require(build['status']=='complete' and build['manifest']==parent['source_manifest'] and all(build['executables'][k]==p for k,p in parent['executables'].items()),'Unchanged native build')
        require(qual['status']=='passed_local_response_native_qualification' and all(qual[k]==parent[k] for k in ('source_manifest','build_receipt','executables')) and all(qual[k] for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Original numerical qualification')
        verify(qual['implementation']);verify(qual['result'])
        h={k:module(parent[k],'saved_witness_'+k) for k in base.HELPERS};h['native_helper'].invoke=forbidden
        manifest=read(parent['manifest']);domain=read(parent['domain']);data=read(parent['development_data']);targets=read(parent['targets'])
        ctx=SimpleNamespace(h=h,m=manifest,config=parent,data=data,targets=targets,target_inventory=read(parent['target_inventory']))
        def render(values,cid,exports=None):return base.render(values,manifest['rows'][cid],domain,h['resolver'],exports)
        ctx.render=render
        core=module(cfg['core'],'saved_witness_core');core.clarabel.DefaultSolver=forbidden
        declaration=read(cfg['campaign_declaration']);origin=read(declaration['start']['artifact'])
        fixed=core.make_fixed(origin,manifest,read(sources['decomposition']),domain,declaration)
        compare(r['center_metrics'],core.assess(center,fixed),'Admitted final center metrics')
        cells={};findings={}
        for row in inv['rows']:
            meta=r['rows'][row['id']];require(meta==r['vectors'][row['key']] and meta['status']=='complete' and meta['parameters']==row['parameters'],'Fixed complete row')
            v=read(meta['artifact']);require(v['status']=='complete' and v['parameters']==row['parameters'] and token(v['parameters'])==row['key'],'Exact native point')
            require(set(v['cells'])==set(manifest['rows']) and len(v['cells'])==46,'Full46 denominator')
            assembled=base.Evaluator._assemble(ctx,v['cells'])
            compare({k:v[k] for k in ('model','response','score')},{k:assembled[k] for k in ('model','response','score')},'Stored response and531 score')
            for cid,rec in v['cells'].items():
                key=rec['cache_key'];native_row=manifest['rows'][cid];text,render_meta=render(v['parameters'],cid)
                identity={'cell_id':cid,'effective_values':render_meta['physical_values_a_b_c_d_R_epsilon'],'xyz':native_row['source_input']['xyz'],'parameters':native_row['source_input']['parameters'],'key_sha256':hashlib.sha256(text.encode()).hexdigest(),'implementation':parent['implementation'],'source_manifest':parent['source_manifest'],'executables':parent['executables'],'helpers':{k:parent[k] for k in base.HELPERS}}
                rp=rec['request_provenance'];require(rec['status']=='complete' and key==token(identity)==rp['cache_key'] and rec['cache_identity']==identity==rp['cache_identity'] and rp['vector']==v['parameters'] and rp['render']==render_meta,'Exact cell/cache/render')
                require(key not in cells and not rec['inherited_native'] and rec['root_release']==cfg['root_release'] and rec['job_id'] is None,'Unique new cell; no hidden inherited cost')
                body=read(pin(lp.parent/'cells'/key/'RESULT.json'));require(body=={k:x for k,x in rec.items() if k!='request_provenance'},'Immutable final cell body')
                base.Evaluator._saved_inputs(ctx,rec,native_row,v['parameters'],lp.parent/'cells'/key)
                cells[key]=rec
            prediction=next(x for x in affine['cases'][row['case']]['predictions'] if x['fraction']==row['fraction'])
            require(prediction['parameters']==v['parameters'] and prediction['key']==row['key'],'Exact affine physical vector')
            measured=runner.summarize(v,prediction,core,fixed,targets);compare(meta['comparison'],measured,'Saved actual-vs-affine comparison')
            gates=independent25(v,data['rows']);compare(gates,measured['metrics']['original25'],'Independent original25')
            require(measured['model_admitted'] is False and measured['selected_candidate'] is None,'Neither witness promoted')
            center_scores=r['center_metrics']['original_block_scores'];scores=measured['metrics']['original_block_scores']
            findings[row['id']]={'key':row['key'],'vector':meta['artifact'],'actual':measured,'independent_original25':gates,'physical_pass_count':sum(x['passed'] for x in gates.values()),'physical_failures':{k:x for k,x in gates.items() if not x['passed']},'component_q_decreases_from_center':measured['metrics']['q']<r['center_metrics']['q'],'electrical_norm_nonincrease':{k:scores[k]<=center_scores[k] for k in ('pair_mu','C_mu','chi')},'affine_solver_status':affine['cases'][row['case']]['solver_status'],'affine_solver_success':affine['cases'][row['case']]['solver_success'],'diagnostic_only':True}
        require(len(cells)==92,'Two full distinct46-cell executions')
        observed=[]
        for num,returned in returns.items():
            attempt=returned['attempt'];rec=cells[attempt['cache_key']];receipts=returned['receipt_delta']
            require(len(receipts)==1 and receipts[0]['returncode']==0 and receipts[0] in rec['receipts'],'One successful child receipt')
            require(returned['program_inputs']==rec['program_inputs'][attempt['program']],'Exact program inputs')
            require(receipts[0]['executor']['mode']=='login' and receipts[0]['executor']['job_id'] is None,'Truthful child allocation')
            observed.append((attempt['cache_key'],attempt['program']))
        require(len(set(observed))==sum(len(x['receipts']) for x in cells.values())==112,'112 distinct native child receipts')
        report.update(status='passed_two_direction_native_validation_actual_review',rows=findings,accounting={'new_vectors':2,'new_cells':92,'native_invocation_attempts':112,'child_receipts':112,'analyze':sum(p=='analyze' for _,p in observed),'testgrad':sum(p=='testgrad' for _,p in observed),'historical_vectors_read':2,'historical_scope':'admitted final center and original fixed-normalization origin only; no1090-body re-admission','native_parser_repeated':False},measured_total_CPU_seconds=math.fsum(r['cost'][k] for k in ('controller_CPU_seconds','child_user_CPU_seconds','child_system_CPU_seconds')),interpretation='Saved input/receipt/score/diagnostic audit passes independently of physical witness success. No candidate selected, new fit, physical-accuracy or global-capacity conclusion.')
        return report
    except Exception as exc:report.update(status='failed_two_direction_actual_review',error=repr(exc));return report
    finally:
        report['checks']=CHECKS;report['review_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
        report['scope_limits']=['Raw native admission is reused from the unchanged qualified evaluator; this review checks literal inputs/output hashes/receipts and recomputes saved response/objective arithmetic, without rerunning molecular parsers or matrix solves.','The inherited1e-10 arithmetic comparison does not alter physical gates, electrical acceptance or native admission thresholds.','B and D are fixed consumed-development diagnostics. Physical failures and AlmostSolved/local-affine limitations remain explicit; neither is automatically promoted.','References retain common-basis/electronic/D4/SCF/state limitations; no heldout or new reference read.']
        (out/'REVIEW.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);args=p.parse_args()
    result=run(json.loads(Path(args.config).read_text()),args.out)
    print(json.dumps({'status':result['status'],'review':pin(Path(args.out)/'REVIEW.json'),'error':result['error']}))
    raise SystemExit(bool(result['error']))
