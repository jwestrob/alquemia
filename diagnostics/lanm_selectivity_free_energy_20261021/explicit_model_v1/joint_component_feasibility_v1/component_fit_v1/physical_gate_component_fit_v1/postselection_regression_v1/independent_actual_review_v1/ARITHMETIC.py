"""Direct saved-token audit of eight consumed cycles; no admission suite or native calls."""
import argparse,hashlib,json,math,os,re,sys,time,importlib.util,datetime
from types import SimpleNamespace
from decimal import Decimal,getcontext
from pathlib import Path
getcontext().prec=50
HERE=Path(__file__).resolve().parent;ROOT=next(p for p in HERE.parents if (p/'diagnostics').is_dir() and (p/'workspaces').is_dir())
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/postselection_regression_v1'
PINS={}
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):
    assert pin(p['path'])==p,p;PINS[p['path']]=p['sha256'];return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def txt(p):return verify(p).read_text()
def D(v):return Decimal(str(v).replace('D','E'))
def norm(v):return math.sqrt(math.fsum(float(x)**2 for x in v))
def arithmetic(*v):return 128*24*sys.float_info.epsilon*max([1.]+[abs(float(x)) for x in v])

def main(args,out):
    tick,cpu=time.monotonic(),time.process_time()
    for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:assert os.environ.get(key)=='1'
    binding=read(pin(args.config));assert binding['implementation']==pin(__file__)
    ready=read(binding['producer_ready']);rp=binding['result'];assert rp is not None
    assert binding['output_directory']==str(out)
    r=read(rp);assert r['status']=='consumed_trimer_regression_complete' and 0<=r['native_program_invocations']<=6 and r['error'] is None
    verify(r['implementation']);m=read(r['manifest']);ec=read(r['configuration']['evaluator_config']);em=read(ec['manifest']);selected=read(r['selected'])
    assert r['parameters']==selected['parameters'] and r['fit'] is False and r['reserved_outcomes_opened'] is False
    assert len(r['rows'])==6 and len(r['cycles'])==8 and all(x['status']=='complete' for x in r['rows'].values())
    verify(ec['executables']['analyze'])
    cfg=read(binding['producer_configuration']);assert cfg==r['configuration']
    assert cfg['implementation']==ready['implementation'] and cfg['manifest']==ready['manifest']
    assert cfg['specification']==ready['specification']==m['specification'];verify(cfg['specification'])
    assert set(r['rows'])==set(m['rows']) and set(r['cycles'])==set(m['cycles'])
    assert set(selected['cells'])==set(em['rows']) and len(selected['cells'])==46 and all(c['status']=='complete' for c in selected['cells'].values())
    release=read(cfg['root_release']);source_review=read(cfg['source_review'])
    assert release['status']=='released_physical_gate_postselection_regression' and release['maximum_native_calls']==6 and release['workers']==1 and release['reserved_outcomes_authorized'] is False
    for key in ('implementation','manifest','fit_result','fit_review','evaluator_config','source_review','hostname','output_directory'):assert release[key]==cfg[key]
    assert release['selected']==r['selected'] and r['implementation']==cfg['implementation']==m['implementation']==source_review['implementation']
    assert source_review['status']=='passed_physical_gate_postselection_regression_source_review' and source_review['manifest']==r['manifest']
    fit=read(r['fit_result']);fitreview=read(r['fit_review']);fitcfg=read(fit['configuration'])
    assert fit['status']=='bounded_physical_gate_component_fit_complete' and fit['mode']=='run' and fit['candidate_frozen'] and fit['error'] is None and fit['selected']==r['selected']
    assert fitreview['status']=='passed_physical_gate_component_fit_actual_review' and fitreview['selected_native_admission_passed'] is True
    assert fitreview['result']==r['fit_result'] and fitreview['selected']==r['selected'] and fitreview['configuration']==fit['configuration'] and fitreview['implementation']==fit['implementation']
    assert fitcfg['parent_configuration']==cfg['evaluator_config'] and fitcfg['parent_evaluator']==ec['implementation'] and fitcfg['declaration']==m['campaign_declaration']
    assert r['selected_metrics']==fit['selected_metrics'] and r['parameters']==fit['vectors'][fit['selected_key']]['parameters']
    assert r['job_id'] is None and r['executor']['mode']=='login' and r['executor']['job_id'] is None and r['executor']['native_workers']==1
    assert r['cost']['scheduler_allocation'] is None and r['cost']['scheduler_allocated_CPU_seconds']==0
    assert release['execution_mode']=='login' and release['aggregate_login_limit']==16 and r['executor']['hostname']==cfg['hostname'] and r['executor']['run_directory']==cfg['output_directory']
    assert r['cost']['inherited_costs_not_added'] is True and r['cost']['native_workers']==r['cost']['numerical_threads']==1
    histories=[]
    for authority in m['inputs']['historical_regressions']:
        saved=read(authority['result']);review=read(authority['review'])
        assert saved['status']=='consumed_trimer_regression_complete' and saved['error'] is None
        assert review['status']=='passed_consumed_trimer_arithmetic_review' and review['result']==authority['result']
        assert set(saved['rows'])==set(m['rows']) and len(saved['cycles'])==8 and all(x['status']=='complete' for x in saved['rows'].values())
        histories.append((authority,saved))
    assert len(histories)==2 and histories[0][0]['result']==m['previous_regression'] and histories[0][0]['review']==m['previous_regression_review']
    previous=histories[0][1];previous_review=read(m['previous_regression_review']);joint=histories[1][1]
    assert r['previous_selected_regression']==m['previous_regression'] and r['previous_selected_cycles']==previous['cycles']
    assert fit['reserved_outcomes_opened'] is False and fit['automatic_extension'] is False
    assert r['original_EF_checks_passed']==selected['score']['original_EF_checks_passed']==fit['selected_metrics']['all25_original_physical_gates'] is True
    campaign=read(m['inputs']['campaign_ready'])
    for key in ('implementation','adapter','core'):
        assert fitcfg[key]==campaign[key];verify(fitcfg[key])
    assert fit['selected_key']==fit['current_key'] and fit['selected']==fit['vectors'][fit['selected_key']]['artifact']
    assert hashlib.sha256(json.dumps(selected['parameters'],sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()==fit['selected_key']
    
    # Qualified parser/render reuse is confined to the six trimer rows; never invoke native.
    def module(ref,name):
        spec=importlib.util.spec_from_file_location(name,verify(ref));obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj
    def forbidden(*args,**kwargs):raise AssertionError('Saved review forbids molecular calls')
    ev=module(ec['implementation'],'review_scientific_evaluator');h={k:module(ec[k],'regression_review_'+k) for k in ev.HELPERS}
    h['native_helper'].invoke=forbidden
    domain=read(ec['domain']);proxy=SimpleNamespace(h=h,m=em,config=ec)
    qualification=read(ec['qualification_review']);build=read(ec['build_receipt'])
    assert qualification['status']=='passed_local_response_native_qualification' and build['status']=='complete'
    for key in ('source_manifest','build_receipt','executables'):assert qualification[key]==ec[key]
    assert build['manifest']==ec['source_manifest'] and build['executables']==ec['executables']
    for key in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed'):assert qualification[key] is True
    last_end=None;native_replays={};row_kinds={};new_calls=0
    for cid,cell in r['rows'].items():
        row=m['rows'][cid];metal=row['metal'];current=ev.effective(r['parameters'],metal,domain)
        match=next(((a,s) for a,s in histories if current==ev.effective(s['parameters'],metal,domain)),None)
        if current==ev.effective(m['initial_vector'],metal,domain):
            assert cell==m['baseline_cells'][cid] and cell['inherited_native']
            row_kinds[cid]='inherited_original_baseline'
            native_replays[cid]=dict(status='reused_admitted_exact_baseline',producer_result=cell['producer_result'],producer_review=cell['producer_review'])
            verify(cell['producer_result']);verify(cell['producer_review']);continue
        if match is not None:
            authority,saved=match;expected=dict(saved['rows'][cid],reused_native=True,reuse_authority=authority['result'],reuse_review=authority['review'])
            assert cell==expected and cell['status']=='complete'
            row_kinds[cid]='inherited_historical_regression'
            if cell['inherited_native']:
                native_replays[cid]=dict(status='reused_admitted_historical_baseline',authority=authority);continue
        else:
            assert not cell['inherited_native'] and not cell.get('reused_native',False)
            row_kinds[cid]='new_native';new_calls+=1
        inp=cell['program_inputs']['analyze'];assert len(cell['receipts'])==1;receipt=cell['receipts'][0]
        assert receipt['returncode']==0 and receipt['argv']==[ec['executables']['analyze']['path'],'system.xyz','EM']
        rendered,meta=ev.render(r['parameters'],row,domain,h['resolver'],exports=False)
        assert txt(inp['key'])==rendered and cell['render']==meta
        for key in ('xyz','parameters'):assert verify(inp[key]).read_bytes()==verify(row['source_input'][key]).read_bytes()
        assert receipt['timeout_seconds']==300 and all(v=='1' for v in receipt['threads'].values())
        if row_kinds[cid]=='new_native':
            assert receipt['executor']==cell['executor']==r['executor']
            directory=Path(receipt['cwd']);assert directory==Path(cfg['output_directory'])/'cells'/cid
            assert json.loads((directory/'RESULT.json').read_text())==cell
            start=datetime.datetime.fromisoformat(receipt['started_utc']);end=datetime.datetime.fromisoformat(receipt['ended_utc'])
            assert end>=start and (last_end is None or start>=last_end);last_end=end
        parsed=ev.Evaluator._parse(proxy,cell,row);assert all(cell[k]==v for k,v in parsed.items())
        native_replays[cid]=dict(status='passed_saved_native_input_output_replay',kind=row_kinds[cid],input=inp,output=receipt['output'])
    assert new_calls==r['native_program_invocations'] and new_calls+sum(v!='new_native' for v in row_kinds.values())==6
    targets=read(m['sources']['target_inventory']);old=read(targets['artifacts']['trimer_comparison'])
    artifacts={'model_moments':read(targets['artifacts']['model_moments']),'trimer_model':read(old['inputs']['model_result'])}
    baseline=read(m['sources']['baseline_result']);mp=read(m['sources']['moment_preparation'])
    needed={k for c in m['cycles'].values() for k in c['model_terms']};values={};terms={};max_term_mu=0.
    assert len(needed)==42
    
    def raw_energy(mom):
        text=txt(mom['analyze_output']);v=re.findall(r'^\s*Total Potential Energy\s*:\s*([-+\d.EeDd]+)\s+Kcal/mole',text,re.M)
        assert len(v)==1 and float(D(v[0]))==mom['energy_kcal_mol'];return D(v[0]),text
    
    for cid in sorted(needed):
        if cid in r['rows'] or cid in selected['cells']:
            cell=r['rows'].get(cid,selected['cells'].get(cid)); assert cell['status']=='complete'
            row=m['rows'].get(cid,em['rows'].get(cid));perm=[D(x) for x in row['permanent_lab_Debye']]
            if cell['inherited_native']:
                basecell=m['baseline_cells'].get(cid,em['baseline_cells'].get(cid));assert all(cell[k]==v for k,v in basecell.items())
                producer=read(cell['producer_result']);review=read(cell['producer_review'])
                assert review['status']=='passed_baseline_matrix_inventory_review' and review['result']==cell['producer_result']
                raw=producer['rows'][cell['producer_cell_id']];assert raw['status']=='complete';mom=raw['moments']
                energy,text=raw_energy(mom);lines=[x.split() for x in txt(mom['induced_output']).splitlines() if x.strip()]
                n=int(lines[0][0]);assert n==row['atoms'] and len(lines)==n+1 and [int(x[0]) for x in lines[1:]]==list(range(1,n+1))
                mu=[perm[j]+sum(D(x[j+2]) for x in lines[1:]) for j in range(3)]
                authority='admitted_inherited_original_baseline'
            else:
                mom=cell['moments'];energy,text=raw_energy(mom)
                n=row['atoms'];dual=[x.split()[2:] for x in text.splitlines() if x.startswith('LOCAL_RESPONSE_INDUCED')][-n:]
                assert len(dual)==n and all(len(x)==6 for x in dual)
                mu=[perm[j]+D(em['units']['debye_per_eA'])*sum((D(x[j])+D(x[j+3]))/2 for x in dual) for j in range(3)]
                authority=row_kinds[cid] if cid in r['rows'] else 'frozen_selected_artifact'
            stored=cell['mu_Debye'];assert float(energy)==cell['energy'];source=mom['analyze_output']
        else:
            cell=em['static'][cid];src=cell['source']
            if 'artifact' in src:
                mom=artifacts[src['artifact']]
                for k in src['json_path']:mom=mom[k]
                perm=[D(x) for x in cell['permanent_lab_Debye']]
            else:
                assert src['result']==m['sources']['baseline_result'];mom=baseline['rows'][src['cell']]['moments']
                perm=[D(x) for x in mp['rows'][cid]['permanent_lab_Debye']]
            energy,text=raw_energy(mom);lines=[x.split() for x in txt(mom['induced_output']).splitlines() if x.strip()]
            n=int(lines[0][0]);assert len(lines)==n+1 and [int(x[0]) for x in lines[1:]]==list(range(1,n+1))
            mu=[perm[j]+sum(D(x[j+2]) for x in lines[1:]) for j in range(3)]
            assert float(energy)==cell['energy_kcal_mol'];stored=cell['mu_Debye'];authority='unchanged_saved_fragment';source=mom['analyze_output']
        discrepancy=max(abs(float(mu[j])-stored[j]) for j in range(3));max_term_mu=max(max_term_mu,discrepancy)
        assert discrepancy<=arithmetic(*mu,*stored)
        values[cid]=(energy,mu);terms[cid]=dict(authority=authority,raw_output=source,energy_decimal=str(energy),mu_Debye_decimal=[str(v) for v in mu],maximum_mu_arithmetic_difference=discrepancy)
    cycles={};maxE=maxM=0.
    for cid,c in m['cycles'].items():
        a=r['cycles'][cid];o=old['cycles'][cid]
        assert a['status']=='complete' and not a['missing'] and a['model_terms']==c['model_terms']==o['model_terms'] and len(c['model_terms'])==7 and sorted(c['model_terms'].values())==[-1,-1,-1,1,1,1,1]
        for key in ['reference_electronic_E3_kcal_mol','reference_mu3_Debye']:assert a[key]==c[key]==o[key]
        assert a['baseline_E3_kcal_mol']==c['model_E3_kcal_mol']==o['model_E3_kcal_mol'] and a['baseline_mu3_Debye']==c['model_mu3_Debye']==o['model_mu3_Debye']
        E=sum(D(co)*values[k][0] for k,co in c['model_terms'].items());M=[sum(D(co)*values[k][1][j] for k,co in c['model_terms'].items()) for j in range(3)]
        ed=abs(float(E)-a['candidate_E3_kcal_mol']);md=max(abs(float(M[j])-a['candidate_mu3_Debye'][j]) for j in range(3));maxE=max(maxE,ed);maxM=max(maxM,md)
        assert ed<=arithmetic(E,a['candidate_E3_kcal_mol']) and md<=arithmetic(*M,*a['candidate_mu3_Debye'])
        ref=D(c['reference_electronic_E3_kcal_mol']);bm=[D(x) for x in c['model_mu3_Debye']];qm=[D(x) for x in c['reference_mu3_Debye']]
        err=E-ref;be=D(c['model_E3_kcal_mol'])-ref;me=[M[j]-qm[j] for j in range(3)];bme=[bm[j]-qm[j] for j in range(3)]
        assert abs(float(err)-a['energy_error_kcal_mol'])<=arithmetic(err) and max(abs(float(me[j])-a['dipole_error_Debye'][j]) for j in range(3))<=arithmetic(*me)
        cycles[cid]=dict(candidate_E3_kcal_mol=float(E),candidate_E3_decimal=str(E),baseline_E3_kcal_mol=float(c['model_E3_kcal_mol']),reference_E3_kcal_mol=float(ref),candidate_energy_error=float(err),baseline_energy_error=float(be),absolute_energy_error_change=float(abs(err)-abs(be)),candidate_mu3_Debye=[float(x) for x in M],baseline_mu3_Debye=[float(x) for x in bm],reference_mu3_Debye=[float(x) for x in qm],candidate_mu_error_norm_Debye=norm(me),baseline_mu_error_norm_Debye=norm(bme),mu_error_norm_change_Debye=norm(me)-norm(bme),original_reference_print_bound_Debye_per_component=a['original_reference_print_bound_Debye_per_component'],original_model_print_bound_Debye_per_component=a['original_model_print_bound_Debye_per_component'],coefficients=c['model_terms'])
        prior=previous['cycles'][cid]
        for key in ('model_terms','reference_electronic_E3_kcal_mol','reference_mu3_Debye','baseline_E3_kcal_mol','baseline_mu3_Debye'):assert prior[key]==a[key]
        pe=D(prior['candidate_E3_kcal_mol'])-ref;pm=[D(x)-qm[j] for j,x in enumerate(prior['candidate_mu3_Debye'])]
        cycles[cid].update(previous_selected_E3_kcal_mol=prior['candidate_E3_kcal_mol'],previous_selected_mu3_Debye=prior['candidate_mu3_Debye'],previous_selected_energy_error=float(pe),previous_selected_mu_error_norm_Debye=norm(pm),absolute_energy_error_change_from_previous=float(abs(err)-abs(pe)),mu_error_norm_change_from_previous_Debye=norm(me)-norm(pm),candidate_mu_error_vector_Debye=[float(x) for x in me],previous_selected_mu_error_vector_Debye=[float(x) for x in pm],baseline_mu_error_vector_Debye=[float(x) for x in bme])
        jp=joint['cycles'][cid]
        for key in ('model_terms','reference_electronic_E3_kcal_mol','reference_mu3_Debye','baseline_E3_kcal_mol','baseline_mu3_Debye'):assert jp[key]==a[key]
        for key in ('original_reference_print_bound_Debye_per_component','original_model_print_bound_Debye_per_component'):assert a[key]==prior[key]==jp[key]
        assert a['original_reference_print_bound_Debye_per_component']==c['reference_print_bound_Debye_per_component'] and a['original_model_print_bound_Debye_per_component']==c['model_print_bound_Debye_per_component']
        je=D(jp['candidate_E3_kcal_mol'])-ref;jm=[D(x)-qm[j] for j,x in enumerate(jp['candidate_mu3_Debye'])]
        cycles[cid].update(prior_joint_fit_E3_kcal_mol=jp['candidate_E3_kcal_mol'],prior_joint_fit_mu3_Debye=jp['candidate_mu3_Debye'],prior_joint_fit_energy_error=float(je),prior_joint_fit_mu_error_norm_Debye=norm(jm),absolute_energy_error_change_from_prior_joint_fit=float(abs(err)-abs(je)),mu_error_norm_change_from_prior_joint_fit_Debye=norm(me)-norm(jm))
    paired={}
    for n in [8,9]:
        for context in ['near','far']:
            la=cycles[f'aquo{n}_central_La__{context}'];dy=cycles[f'aquo{n}_central_Dy__{context}'];p={}
            for kind in ['candidate','baseline','reference','previous_selected','prior_joint_fit']:
                p[kind+'_E3_Dy_minus_La']=float(D(dy[kind+'_E3_kcal_mol'])-D(la[kind+'_E3_kcal_mol']))
                p[kind+'_mu3_Dy_minus_La_Debye']=[float(D(y)-D(x)) for x,y in zip(la[kind+'_mu3_Debye'],dy[kind+'_mu3_Debye'])]
            for kind in ['candidate','baseline','previous_selected','prior_joint_fit']:
                p[kind+'_energy_error']=p[kind+'_E3_Dy_minus_La']-p['reference_E3_Dy_minus_La']
                p[kind+'_mu_error_norm_Debye']=norm([y-x for x,y in zip(p['reference_mu3_Dy_minus_La_Debye'],p[kind+'_mu3_Dy_minus_La_Debye'])])
            p['absolute_energy_error_change']=abs(p['candidate_energy_error'])-abs(p['baseline_energy_error']);p['mu_error_norm_change_Debye']=p['candidate_mu_error_norm_Debye']-p['baseline_mu_error_norm_Debye']
            p['absolute_energy_error_change_from_previous']=abs(p['candidate_energy_error'])-abs(p['previous_selected_energy_error'])
            p['mu_error_norm_change_from_previous_Debye']=p['candidate_mu_error_norm_Debye']-p['previous_selected_mu_error_norm_Debye']
            p['candidate_mu_error_vector_Debye']=[y-x for x,y in zip(p['reference_mu3_Dy_minus_La_Debye'],p['candidate_mu3_Dy_minus_La_Debye'])]
            p['absolute_energy_error_change_from_prior_joint_fit']=abs(p['candidate_energy_error'])-abs(p['prior_joint_fit_energy_error'])
            p['mu_error_norm_change_from_prior_joint_fit_Debye']=p['candidate_mu_error_norm_Debye']-p['prior_joint_fit_mu_error_norm_Debye']
            paired[f'W{n}_{context}']=p
    # Check the producer's additional comparison table without introducing a score gate.
    comp=r['comparison'];assert comp['cycle_count']==8 and comp['paired_count']==4 and comp['scientific_threshold_added'] is False
    labels={'original_pre_response_baseline':'baseline','old_component_terminal':'previous_selected','new_terminal':'candidate'}
    for cid,v in cycles.items():
        report=comp['cycles'][cid];assert report['status']=='complete' and set(report['models'])==set(labels)
        for label,kind in labels.items():
            got=report['models'][label]
            for key,target in [('E3_kcal_mol',v[kind+'_E3_kcal_mol']),('energy_error_kcal_mol',v[kind+'_energy_error']),('dipole_error_norm_Debye',v[kind+'_mu_error_norm_Debye'])]:assert abs(got[key]-target)<=arithmetic(got[key],target)
            assert max(abs(a-b) for a,b in zip(got['mu3_Debye'],v[kind+'_mu3_Debye']))<=arithmetic(*got['mu3_Debye'],*v[kind+'_mu3_Debye'])
    for cid,v in paired.items():
        n,context=cid[1:].split('_');report=comp['paired_Dy_minus_La'][f'aquo{n}__{context}'];assert report['status']=='complete' and set(report['models'])==set(labels)
        for label,kind in labels.items():
            got=report['models'][label]
            for key,target in [('Dy_minus_La_E3_kcal_mol',v[kind+'_E3_Dy_minus_La']),('energy_error_kcal_mol',v[kind+'_energy_error']),('dipole_error_norm_Debye',v[kind+'_mu_error_norm_Debye'])]:assert abs(got[key]-target)<=arithmetic(got[key],target)
            assert max(abs(a-b) for a,b in zip(got['Dy_minus_La_mu3_Debye'],v[kind+'_mu3_Dy_minus_La_Debye']))<=arithmetic(*got['Dy_minus_La_mu3_Debye'],*v[kind+'_mu3_Dy_minus_La_Debye'])
    result=dict(status='passed_consumed_trimer_arithmetic_review',implementation=pin(__file__),result=rp,producer=r['implementation'],manifest=r['manifest'],selected=r['selected'],cycles=cycles,paired_Dy_minus_La=paired,raw_terms=terms,term_count=len(terms),cycle_count=8,paired_count=4,source_pins=[dict(path=k,sha256=v) for k,v in PINS.items()],maximum_energy_sum_replay_difference_kcal_mol=maxE,maximum_mu_sum_replay_difference_Debye=maxM,maximum_term_mu_arithmetic_difference_Debye=max_term_mu,energy_absolute_errors_worsened=sum(c['absolute_energy_error_change']>0 for c in cycles.values()),dipole_vector_norm_errors_improved=sum(c['mu_error_norm_change_Debye']<0 for c in cycles.values()),paired_energy_absolute_errors_worsened=sum(c['absolute_energy_error_change']>0 for c in paired.values()),new_native_calls=0,original_EF_checks_passed=r['original_EF_checks_passed'],scientific_threshold_added=False,physical_accuracy_qualified=False,limits=['Arithmetic review reuses admitted baseline/reference and permanent-frame values; it does not repeat native/QM admission or the full fit audit.','Eight consumed E3 cycles are outside the objective; their dependencies overlap objective cells, so this is not held-out validation.','No E3 acceptance threshold was declared; worsening is descriptive evidence for a progression decision.','PBE0 and frozen electronic-state limitations persist; no unique missing interaction or parameter mechanism is established.','Static fragments retain native six-decimal induced moments and the original printed-precision envelopes.'],cost=dict(saved_review_wall_seconds=time.monotonic()-tick,saved_review_CPU_seconds=time.process_time()-cpu,threads=1,producer_application=r['cost'],scheduler_allocation=None,scheduler_allocated_CPU_seconds=0))
    
    result.update(configuration=binding['producer_configuration'],root_release=cfg['root_release'],source_review=cfg['source_review'],fit_result=r['fit_result'],fit_review=r['fit_review'],previous_selected_regression=m['previous_regression'],previous_selected_review=m['previous_regression_review'],native_input_output_replay=native_replays,new_cell_native_admission_passed=True,reviewer_parent=pin(ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/postselection_regression_v1/independent_actual_review_v1/ARITHMETIC.py'),row_kinds=row_kinds,new_native_programs_reviewed=new_calls,inherited_trimer_rows=6-new_calls,prior_joint_fit_regression=histories[1][0]['result'],prior_joint_fit_review=histories[1][0]['review'],comparison_counts={
        'cycle_energy_abs_improved_from_previous':sum(v['absolute_energy_error_change_from_previous']<0 for v in cycles.values()),
        'cycle_mu_norm_improved_from_previous':sum(v['mu_error_norm_change_from_previous_Debye']<0 for v in cycles.values()),
        'paired_energy_abs_improved_from_previous':sum(v['absolute_energy_error_change_from_previous']<0 for v in paired.values()),
        'paired_mu_norm_improved_from_previous':sum(v['mu_error_norm_change_from_previous_Debye']<0 for v in paired.values()),
        'cycle_energy_abs_worse_than_original_baseline':sum(v['absolute_energy_error_change']>0 for v in cycles.values()),
        'cycle_mu_norm_worse_than_original_baseline':sum(v['mu_error_norm_change_Debye']>0 for v in cycles.values()),
        'paired_energy_abs_worse_than_original_baseline':sum(v['absolute_energy_error_change']>0 for v in paired.values()),
        'paired_mu_norm_worse_than_original_baseline':sum(v['mu_error_norm_change_Debye']>0 for v in paired.values())})
    result['limits'][0]='Zero to six new cells and exact historical effective-metal matches replay rendering/input bytes and the qualified native parser; original baseline reuse binds its admitted authority. All42 cycle terms replay raw tokens. Selected/reference/static admissions are reused without another full fit or quantum audit.'
    lines=['Consumed development regression: eight E3/mu3 cycles; not holdout validation.',
           'Signed energy errors are model minus reference (kcal/mol); dipole errors below are Euclidean vector-error norms (Debye).',
           'Order: original baseline / prior selected joint fit / old component terminal / new physical-gate terminal.',
           '', 'Cycle | signed E3 error: baseline / joint / old component / new | mu3 error norm: baseline / joint / old component / new']
    for cid,v in cycles.items():
        lines.append(f"{cid} | {v['baseline_energy_error']:.8f} / {v['prior_joint_fit_energy_error']:.8f} / {v['previous_selected_energy_error']:.8f} / {v['candidate_energy_error']:.8f} | {v['baseline_mu_error_norm_Debye']:.8f} / {v['prior_joint_fit_mu_error_norm_Debye']:.8f} / {v['previous_selected_mu_error_norm_Debye']:.8f} / {v['candidate_mu_error_norm_Debye']:.8f}")
    lines.extend(['','Paired Dy-minus-La | signed E3 error: baseline / joint / old component / new | mu3 error norm: baseline / joint / old component / new'])
    for cid,v in paired.items():
        lines.append(f"{cid} | {v['baseline_energy_error']:.8f} / {v['prior_joint_fit_energy_error']:.8f} / {v['previous_selected_energy_error']:.8f} / {v['candidate_energy_error']:.8f} | {v['baseline_mu_error_norm_Debye']:.8f} / {v['prior_joint_fit_mu_error_norm_Debye']:.8f} / {v['previous_selected_mu_error_norm_Debye']:.8f} / {v['candidate_mu_error_norm_Debye']:.8f}")
    lines.extend(['','Comparison counts: '+json.dumps(result['comparison_counts'],sort_keys=True),
                  'No physical E3/mu3 acceptance threshold was introduced. These cycles were consumed earlier and share objective dependencies.',
                  'Original PBE0/frozen-electronic-state and print-precision limits remain; no unique mechanism or physical accuracy follows.',
                  'Actual regression cost: '+json.dumps(r['cost'],sort_keys=True),
                  'Saved review cost: '+json.dumps(result['cost'],sort_keys=True)])
    (out/'REPORT.txt').write_text('\n'.join(lines)+'\n')
    result['report']=pin(out/'REPORT.txt')
    (out/'REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(review=pin(out/'REVIEW.json'),counts={k:result[k] for k in ['term_count','energy_absolute_errors_worsened','dipole_vector_norm_errors_improved','paired_energy_absolute_errors_worsened']},paired=paired,maximum_replay_differences=[maxE,maxM],cost=result['cost'])))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);args=ap.parse_args()
    binding=json.loads(args.config.read_text());out=Path(binding['output_directory']).resolve();out.mkdir(parents=True,exist_ok=False)
    try:main(args,out)
    except Exception as exc:
        (out/'FAILED_REVIEW.json').write_text(json.dumps(dict(status='failed_consumed_trimer_arithmetic_review',implementation=pin(__file__),binding=pin(args.config),error=dict(type=type(exc).__name__,message=str(exc)),new_native_calls=0),indent=2)+'\n')
        raise
