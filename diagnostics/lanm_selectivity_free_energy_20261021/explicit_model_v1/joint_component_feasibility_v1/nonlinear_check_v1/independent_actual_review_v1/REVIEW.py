#!/usr/bin/env python3
"""Saved56-output admission with inherited parsers and independent physical gates."""
import copy
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import tempfile
import time

BASE=Path(__file__).resolve().parents[6]
REL=Path('explicit_model_v1/joint_component_feasibility_v1/nonlinear_check_v1')
D=BASE/'diagnostics/lanm_selectivity_free_energy_20261021'/REL
W=BASE/'workspaces/lanm_selectivity_free_energy_20261021'/REL

def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def read(p):
    if isinstance(p,dict):assert pin(p['path'])==p;p=p['path']
    return json.loads(Path(p).read_text())

def load(p):
    assert pin(p['path'])==p
    spec=importlib.util.spec_from_file_location('reviewed_fixed_half',p['path']);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def normalized(v):return json.loads(json.dumps(v))

def near(a,b):
    assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=1e-10*max(1,abs(a),abs(b))
    return abs(a-b)

def main():
    wall,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    resultpin=pin(W/'run_v1/execution/RESULT.json');r=read(resultpin);cfg=read(r['configuration']);m=load(r['implementation'])
    assert r['status']=='complete_one_fixed_half_native_pilot' and r['error'] is None and r['native_programs']==56 and r['job_id'] is None
    base,parent,cert,origin=m.context(cfg)
    assert r['configuration_fingerprint']==m.fingerprint(cfg)==r['snapshot']['batch_fingerprint']
    assert r['candidate']==cfg['candidate'] and r['candidate']['parameters']==cert['physical_half_parameters']
    release=read(cfg['root_release']);review=read(cfg['source_review']);prep=read(cfg['preparation'])
    assert release['status']=='released_one_fixed_half_native_pilot' and review['status']=='passed_fixed_half_native_source_review'
    for value in (release,review,prep):assert value['configuration_fingerprint']==r['configuration_fingerprint']
    for name in ('parent_configuration','half_review','origin','declaration','source_review','preparation','candidate','limits','executor'):assert release[name]==cfg[name]
    assert release['implementation']==cfg['implementation']==r['implementation']==review['implementation']
    assert release['output_directory']==cfg['output_directory']==str(W/'run_v1/execution') and release['cache_directory']==str(W/'run_v1/execution/cache')
    build=read(parent['build_receipt']);qualification=read(parent['qualification_review'])
    assert build['status']=='complete' and build['manifest']==parent['source_manifest']
    for k in ('analyze','testgrad'):assert pin(parent['executables'][k]['path'])==parent['executables'][k]==build['executables'][k]
    for k in ('source_manifest','build_receipt','executables'):assert qualification[k]==parent[k]
    assert qualification['status']=='passed_local_response_native_qualification' and all(qualification[k] for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed'))
    for p in (parent['source_manifest'],qualification['implementation'],qualification['result']):assert pin(p['path'])==p
    ledger=read(r['snapshot']['ledger']);assert ledger['stopped'] is False and ledger['vector_keys']==[cfg['candidate']['key']] and ledger['native_program_invocations']==ledger['native_returncodes_observed']==56
    assert len(ledger['requests'])==1 and ledger['requests'][0]['status']=='complete' and len(ledger['requests'][0]['cells'])==46
    assert all(v['status']=='complete' for v in ledger['requests'][0]['cells'].values())
    attempts=ledger['native_attempts'];assert len(attempts)==56 and [v['number'] for v in attempts]==list(range(1,57))
    assert all(v['status']=='returned_success' and v['receipt_count']==1 and v['returncode_observed'] for v in attempts)
    assert len({(v['cell_id'],v['program'],v['folder']) for v in attempts})==56
    rows={};events=[];program_counts={'analyze':0,'testgrad':0};maxind=0.;minratio=1.;maxeig=0.;maxchi=0.;matrix_count=0;max_gate_difference=0.
    with tempfile.TemporaryDirectory(prefix='fixed_half_saved_review_') as temp:
        ev=m.evaluator_class(base)(cfg,Path(temp)/'cache','prepare')
        def forbidden(*a,**kw):raise AssertionError('Review cannot invoke native programs')
        ev.h['native_helper'].invoke=forbidden
        cells=r['answer']['cells'];assert set(cells)==set(ev.m['rows']) and len(cells)==46
        for cid,record in cells.items():
            row=ev.m['rows'][cid];assert record['status']=='complete' and record['inherited_native'] is False and record['first_requested_vector']==cfg['candidate']['parameters']
            key,identity,_,meta=ev.cache_identity(cfg['candidate']['parameters'],cid)
            assert (key,identity,meta)==(record['cache_key'],record['cache_identity'],record['render'])
            directory=W/'run_v1/execution/cache/cells'/key;disk=read(directory/'RESULT.json')
            assert disk=={k:v for k,v in record.items() if k!='request_provenance'}
            assert record['executor']==r['executor'] and record['job_id'] is None and record['root_release']==cfg['root_release']
            ev._saved_inputs(record,row,cfg['candidate']['parameters'],directory)
            parsed=ev._parse(record,row)
            for name in ('energy','moments','mu_Debye','local_response_readback'):
                assert normalized(parsed[name])==record[name],(cid,name)
            if row['needs_gradient']:assert normalized(parsed['model'])==record['model']
            if row['needs_export']:
                old,new=record['export_qualification'],parsed['export_qualification']
                assert old['passed'] is True and new['passed'] is True and all(v['passed'] for v in old['checks'].values()) and all(v['passed'] for v in new['checks'].values())
                assert old['two_solve_exact']==new['two_solve_exact']
                for a,b in zip(old['exports'],new['exports']):
                    for field in ('source','sequence','sites','vectors','A','checks','restored','peek'):assert a[field]==normalized(b[field]),(cid,field)
                    matrix_count+=1
                delta=abs(record['chi_a0_cubed']-parsed['chi_a0_cubed']);bound=ev.h['scout_helper'].arithmetic(3*row['atoms'],record['chi_a0_cubed'],parsed['chi_a0_cubed'])
                assert delta<=bound;maxchi=max(maxchi,delta)
            for spectrum in parsed['local_response_readback']['spectra']:
                minratio=min(minratio,spectrum[0]/spectrum[1]);maxeig=max(maxeig,spectrum[3])
            assert len(record['receipts'])==1+int(row['needs_gradient'])
            for program,receipt in zip(('analyze','testgrad'),record['receipts']):
                assert receipt['executor']==r['executor'] and receipt['timeout_seconds']==300 and receipt['returncode']==0
                assert pin(receipt['output']['path'])==receipt['output'] and all(v=='1' for v in receipt['threads'].values())
                assert any(a['cell_id']==cid and a['program']==program and a['folder']==receipt['cwd'] for a in attempts)
                nativeind=ev.h['induction_helper'].induction(receipt['output']['path']);maxind=max(maxind,nativeind['maximum_printed_RMS_residual_Debye'])
                start=datetime.datetime.fromisoformat(receipt['started_utc']);end=datetime.datetime.fromisoformat(receipt['ended_utc']);assert end>=start and receipt['elapsed_seconds']>=0
                events.extend([(start,1),(end,-1)]);program_counts[program]+=1
            rows[cid]=dict(result=pin(directory/'RESULT.json'),outputs=[v['output'] for v in record['receipts']],raw_replay_passed=True)
        assembled=ev._assemble(copy.deepcopy(cells))
        assert normalized(assembled['model'])==r['answer']['model'] and normalized(assembled['response'])==r['answer']['response']
        assert normalized(assembled['score'])==r['answer']['score']
        # Independent physical force projection, paired differences and category weights.
        model=assembled['model'];data=ev.data['rows'];old=assembled['score']['old_EF'];gates={}
        for kind in ('La','Dy','paired'):
            for geo in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central'):
                la,dy=data[geo+'_La'],data[geo+'_Dy'];n=len(la['geometry']['positions_A'])
                residual=[]
                for i in range(n):
                    v=[]
                    for j in range(3):
                        if kind=='paired':v.append((model[geo+'_Dy']['forces'][i][j]-model[geo+'_La']['forces'][i][j])-(dy['reference']['forces'][i][j]-la['reference']['forces'][i][j]))
                        else:v.append(model[geo+'_'+kind]['forces'][i][j]-data[geo+'_'+kind]['reference']['forces'][i][j])
                    residual.append(v)
                raw=[];xyz=la['geometry']['positions_A']
                for category in ('ion','water'):
                    for body in (b for b in la['bodies'] if b['category']==category):
                        ids=body['indices'];center=[math.fsum(xyz[i][j] for i in ids)/len(ids) for j in range(3)]
                        raw.extend(math.fsum(residual[i][j] for i in ids) for j in range(3))
                        if category=='water':
                            for a,b in ((1,2),(2,0),(0,1)):raw.append(math.fsum((xyz[i][a]-center[a])*residual[i][b]-(xyz[i][b]-center[b])*residual[i][a] for i in ids))
                waters=(n-1)//3;weights=[1/6]*3+[1/(12*waters)]*(6*waters);rms=math.sqrt(math.fsum(w*x*x for w,x in zip(weights,raw)));maximum=max(abs(x) for x in raw)
                stored=old['per_fixture'][kind][geo]
                for a,b in zip(raw,stored['raw_body_components_kcal_mol_A']):max_gate_difference=max(max_gate_difference,near(a,b))
                near(rms,stored['rigid_RMS_kcal_mol_A']);near(maximum,stored['maximum_body_component_error_kcal_mol_A'])
                assert rms<=2 and maximum<=5 and stored['passed'] is True
                gates['force/'+kind+'/'+geo]=dict(RMS=rms,maximum=maximum,passed=True)
        def E(values,geo,metal):return values[geo+'_'+metal]['energy']
        refs={k:v['reference'] for k,v in data.items()};errors={}
        for geo in ('water_minus','water_central','water_plus','aquo9_central'):
            errors['balanced_exchange/'+geo]=(E(model,geo,'Dy')-E(model,geo,'La')-E(model,'aquo8_central','Dy')+E(model,'aquo8_central','La'))-(E(refs,geo,'Dy')-E(refs,geo,'La')-E(refs,'aquo8_central','Dy')+E(refs,'aquo8_central','La'))
        for metal in ('La','Dy'):
            for sign in ('minus','plus'):errors['individual_radial_work/'+metal+'_'+sign]=(E(model,'water_'+sign,metal)-E(model,'water_central',metal))-(E(refs,'water_'+sign,metal)-E(refs,'water_central',metal))
        for sign in ('minus','plus'):errors['paired_radial_work/'+sign]=errors['individual_radial_work/Dy_'+sign]-errors['individual_radial_work/La_'+sign]
        for name,error in errors.items():
            group,key=name.split('/');stored=old['energy_controls'][group][key];near(error,stored['error_kcal_mol']);assert abs(error)<=1 and stored['passed'] is True
            gates['energy/'+name]=dict(error=error,passed=True)
        assert len(gates)==25 and assembled['score']['original_EF_checks_passed'] is True and r['original25_pass'] is True
    running=maximum=0
    for _,change in sorted(events,key=lambda item:(item[0],item[1])):running+=change;maximum=max(maximum,running)
    assert running==0 and maximum<=12 and program_counts=={'analyze':46,'testgrad':10} and matrix_count==8
    cost=r['cost'];assert cost['scheduler_allocated_CPU_seconds']==0 and cost['native_worker_limit']==12 and cost['aggregate_login_budget']==16
    assert all(cost[k]>=0 for k in ('wall_seconds','controller_CPU_seconds','native_child_user_CPU_seconds','native_child_system_CPU_seconds'))
    out=W/'independent_actual_review_v1';out.mkdir(parents=True,exist_ok=False)
    report=dict(status='passed_fixed_half_native_actual_review',implementation=pin(__file__),result=resultpin,configuration=r['configuration'],root_release=cfg['root_release'],producer=r['implementation'],parent_evaluator=cfg['parent_evaluator'],native_source=parent['source_manifest'],build_receipt=parent['build_receipt'],executables=parent['executables'],qualification_review=parent['qualification_review'],candidate=r['candidate'],ledger=r['snapshot']['ledger'],cells=rows,counts=dict(cells=46,programs=program_counts,total_receipts=56,matrix_snapshots=matrix_count,physical_gates=25,maximum_receipt_overlap=maximum),
        numeric_summary=dict(maximum_induction_residual_Debye=maxind,minimum_native_normalized_spectral_ratio=minratio,maximum_native_eigensystem_residual=maxeig,maximum_replayed_chi_difference=maxchi,maximum_independent_body_component_difference=max_gate_difference),original25=gates,original25_pass=True,native_actual_cost=cost,
        native_CPU_seconds_not_double_counted=cost['native_child_user_CPU_seconds']+cost['native_child_system_CPU_seconds'],total_controller_plus_native_CPU_seconds=cost['controller_CPU_seconds']+cost['native_child_user_CPU_seconds']+cost['native_child_system_CPU_seconds'],
        scope=['Inherited unchanged _saved_inputs/_parse admission on all46 outputs, plus unchanged pure assembly; authorship disclosed.','Independent candidate/release/build/receipt/ledger/slot checks and Cartesian-body projection of all15 force gates plus all10 energy gates.','Raw export fields compared exactly; both old/new numerical qualification flags required. Derived LAPACK objects are not compared by exact equality.','Component/prediction interpretation reviewed separately by Singer; no duplicate whole-component audit.'],limits=['Consumed-development native numerical/25-gate admission only; no fresh validation or production promotion.','Original failed SLSQP remains failed; this separately declared fixed-half pilot is distinct.','Current native forces are analytic; prior qualification supplies derivative consistency, no new FD suite run.'],review_cost=dict(wall_seconds=time.monotonic()-wall,CPU_seconds=time.process_time()-cpu,threads=1,native_calls=0,optimizer_calls=0))
    p=out/'REVIEW.json';p.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(review=pin(p),counts=report['counts'],numeric_summary=report['numeric_summary'],review_cost=report['review_cost']),indent=2))

if __name__=='__main__':main()
