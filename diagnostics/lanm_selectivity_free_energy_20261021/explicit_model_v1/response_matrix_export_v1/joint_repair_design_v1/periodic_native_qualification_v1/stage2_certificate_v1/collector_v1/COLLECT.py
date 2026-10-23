"""Finite saved Stage2 collection. No native, FFT, factor or solver call."""
import hashlib
import json
import math
import re
import sys
import time
from decimal import Decimal, localcontext
from pathlib import Path

N=11958
DIM=3*N
STEPS=(Decimal('0.0001'),Decimal('0.00005'))
GATE=Decimal('0.001')
EPS=Decimal.from_float(sys.float_info.epsilon)
BYPASSES=('BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO')
DIRECTION_IDS=('metal_0','metal_1','metal_2','water_rotation','whole_box','strain_0','strain_1','strain_2')


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as stream:sha=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),sha256=sha)


def verify(item):
    require(isinstance(item,dict) and set(item)=={'path','sha256'} and pin(item['path'])==item,'Missing/changed pin')
    return Path(item['path'])


def read(item):return json.loads(verify(item).read_text())


def number(token):
    require(re.fullmatch(r'[+-]?\d+\.\d+(?:[EeDd][+-]?\d+)?',token) is not None,'Malformed native decimal')
    value=Decimal(token.replace('D','E').replace('d','E'))
    require(value.is_finite() and math.isfinite(float(value)),'Nonfinite native value')
    return value


def rounding(token):
    value=number(token)
    return Decimal(5).scaleb(value.as_tuple().exponent-1)


def arithmetic(dim,*values):return Decimal(128)*dim*EPS*max(Decimal(1),*(abs(x) for x in values))


def verdict(label,error,bound,**extra):
    require(error.is_finite() and bound.is_finite() and error>=0 and bound>=0,'Finite numerical check')
    return dict(label=label,error=float(error),bound=float(bound),passed=error<=bound,**extra)


def parse_gradient(path,n=N):
    """Single strict read of full output; no implicit displaced/last-solve choice."""
    with Path(path).open() as stream,localcontext() as ctx:
        ctx.prec=70
        def literal(value):require(stream.readline().strip()==value,'Native line: '+value)
        def row(prefix,width):
            fields=stream.readline().split();head=prefix.split()
            require(fields[:len(head)]==head and len(fields)==len(head)+width,'Native record: '+prefix)
            return fields[len(head):]
        literal('STAGE2_ONE_ANALYTIC_GRADIENT_V1')
        literal('UNITS E=kcal/mol G=positive_dE/dr_kcal/mol/angstrom VIRIAL=kcal/mol')
        energy=row('TOTAL_ENERGY',1)[0];ep=row('POLARIZATION_ENERGY',1)[0]
        electric,dielec=row('CONTRACTION_CONSTANTS electric dielec',2)
        require(number(dielec)>0 and number(electric)>0,'Positive contraction constants')
        number(energy);number(ep);g=[];vir=[]
        for i in range(1,n+1):
            fields=row('TOTAL_GRADIENT',4);require(fields[0]==str(i),'Gradient atom order')
            for token in fields[1:]:number(token)
            g.extend(fields[1:])
        for i in range(1,4):
            fields=row('TOTAL_VIRIAL',4);require(fields[0]==str(i),'Virial row order')
            for token in fields[1:]:number(token)
            vir.extend(fields[1:])
        role=row('ROLE',1)[0];state_id=row('STATE_ID',1)[0];binding=row('INPUT_BINDING',1)[0]
        require(role in ('ACTIVE',)+BYPASSES,'Declared gradient role')
        require(re.fullmatch('[a-f0-9]{64}',binding) is not None,'Input binding digest')
        result=dict(role=role,state_id=state_id,input_binding=binding,energy_token=energy,polarization_token=ep,
                    electric_token=electric,dielec_token=dielec,gradient_tokens=g,virial_tokens=vir,
                    E_G_virial_token_sha256=hashlib.sha256(('\n'.join([energy]+g+vir)+'\n').encode()).hexdigest(),
                    counters=None,PCG_pre_peek=None,shared_h=None,polarization_contraction=None)
        if role=='BYPASS_OLD_ABSENT':
            literal('INTERNAL_COUNTERS unavailable_unchanged_old_archive')
        else:
            counts=row('COUNTS induced_solve direct mutual',3)
            require(all(x.isdecimal() for x in counts),'Nonnegative counters')
            solve,direct,mutual=map(int,counts)
            require(solve==direct==1 and mutual>=1,'One direct/solve; actual mutual applications')
            literal('COUNTS extra_diagnostic_operator_applications=0 internal_FD=0')
            p=row('PCG_PRE_PEEK',5);require(p[0].isdecimal(),'PCG iteration count')
            iterations=int(p[0]);eps,epsd,epsp,old=map(number,p[1:])
            require(3<=iterations<100 and 0<=eps<Decimal('1e-8') and eps<=old and epsd>=0 and epsp>=0,
                    'Original native induction termination')
            result['counters']=dict(induced_solves=solve,direct_initializations=direct,mutual_applications=mutual,
                                    extra_diagnostic_operator_applications=0,internal_FD=0)
            result['PCG_pre_peek']=dict(iterations=iterations,RMS_Debye=float(eps),squared_channel_residuals=[float(epsd),float(epsp)],
                                       previous_RMS_Debye=float(old),tokens=p,post_peek_residual_not_recomputed=True)
            if role=='ACTIVE':
                reported=number(row('SHARED_H_READBACK_MAX_ERROR',1)[0]);require(reported>=0,'Nonnegative shared-h readback')
                hook_error=Decimal(0);scale=Decimal(1);contraction=Decimal(0);dual_max=Decimal(0)
                for i in range(1,n+1):
                    records=[]
                    for prefix,width in (('RAW_DIRECT',7),('CORRECTED_RHS',7),('FINAL_POST_PEEK',8),('CORRECTION',5)):
                        tokens=row(prefix,width);require(tokens[0]==str(i),'Dual/correction site order')
                        records.append(list(map(number,tokens[1:])))
                    raw,rhs,final,correction=records;alpha=final[0];u=final[1:4];v=final[4:];h=correction[1:]
                    require(alpha>=Decimal('1e-8'),'Positive actual intrinsic alpha')
                    scale=max(scale,*(abs(x) for x in raw+rhs+h))
                    hook_error=max(hook_error,*(abs(rhs[j]-raw[j]-h[j%3]) for j in range(6)))
                    dual_max=max(dual_max,*(abs(u[j]-v[j]) for j in range(3)))
                    contraction+=sum((u[j]*rhs[3+j] for j in range(3)),Decimal(0))
                bound=arithmetic(3*n,scale)
                result['shared_h']=verdict('same h in both direct RHS',max(hook_error,reported),bound,
                                            reconstructed_max_error=float(hook_error),reported_max_error=float(reported))
                prediction=-number(electric)/(2*number(dielec))*contraction
                result['polarization_contraction']=verdict('stationary dual polarization energy',abs(number(ep)-prediction),
                    Decimal('1e-8')+arithmetic(3*n,number(ep),prediction),native_kcal_mol=float(number(ep)),
                    contraction_kcal_mol=float(prediction),maximum_u_v_difference=float(dual_max),
                    convention='-electric/(2*dielec) sum_i u_i dot corrected_Ep_i; final post-peek u/v')
        literal('END gradient_calls=1 energy_only_calls=0 internal_FD=0')
        require(not stream.read().strip(),'Trailing native output')
        return result


def projection(central,direction,tangent):
    with localcontext() as ctx:
        ctx.prec=70
        if direction['kind']=='strain':
            coefficients=[Decimal.from_float(float(x))/Decimal.from_float(float(direction['L0_A']))
                          for row in direction['strain_matrix'] for x in row]
            tokens=central['virial_tokens']
        else:
            coefficients=[Decimal.from_float(float(x)) for x in tangent];tokens=central['gradient_tokens']
        require(len(tokens)==len(coefficients),'Projection dimension')
        value=sum((c*number(t) for c,t in zip(coefficients,tokens)),Decimal(0))
        bound=sum((abs(c)*rounding(t) for c,t in zip(coefficients,tokens)),Decimal(0))
        return value,bound


def slope_checks(direction,analytic,analytic_print,energies):
    with localcontext() as ctx:
        ctx.prec=70;steps=[]
        for step in STEPS:
            minus,plus=energies[step]
            fd=(number(plus)-number(minus))/(2*step)
            ep=(rounding(plus)+rounding(minus))/(2*step)
            allowance=ep+analytic_print
            steps.append(dict(half_step_A=float(step),energy_minus_token=minus,energy_plus_token=plus,
                FD_kcal_mol_A=float(fd),analytic_kcal_mol_A=float(analytic),energy_print_bound=float(ep),
                analytic_print_bound=float(analytic_print),printing_bound=float(allowance),
                check=verdict(direction+' FD '+str(step),abs(fd-analytic),GATE+allowance),_fd=fd,_printing=allowance))
        two=verdict(direction+' two-step',abs(steps[0]['_fd']-steps[1]['_fd']),GATE+sum((s['_printing'] for s in steps),Decimal(0)))
        for s in steps:del s['_fd'];del s['_printing']
        return dict(steps=steps,two_step=two,passed=all(s['check']['passed'] for s in steps) and two['passed'])


def load_tangent(item):
    import numpy as np
    a=np.load(verify(item),allow_pickle=False)
    require(a.shape==(N,3) and a.dtype==np.float64 and np.isfinite(a).all(),'Frozen tangent shape/type')
    return a.reshape(-1).tolist()


def expected_states(displaced,bypass,mode):
    require(mode in ('central_scout','full'),'Declared collection mode')
    expected={'central':dict(role='ACTIVE',inputs=displaced['central_prepared_input'])}
    if mode=='full':
        require(len(displaced['rows'])==32 and len({r['id'] for r in displaced['rows']})==32,'All32 unique prepared displacements')
        for row in displaced['rows']:expected[row['id']]=dict(role='ACTIVE',inputs=row['input'],direction=row['direction'],
                                                           step=Decimal(str(row['half_step_A'])),sign=row['sign'])
    require(tuple(r['role'] for r in bypass['rows'])==BYPASSES,'Exact ordered three bypasses')
    for row in bypass['rows']:
        expected[row['role']]=dict(role=row['role'],inputs={k:row['inputs'][v] for k,v in
                                 [('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm')]})
    return expected


def same_inputs(actual,expected):
    require(set(actual)==set(expected)=={'xyz','key','parameters'},'Exact three input files')
    for name,item in actual.items():
        verify(item);verify(expected[name]);require(item['sha256']==expected[name]['sha256'],'Exact prepared input '+name)


def runtime_inputs(actual,expected):
    require(set(actual)==set(expected)=={'xyz','key','parameters'},'Exact runtime input set')
    for name in ('xyz','parameters'):
        verify(actual[name]);verify(expected[name]);require(actual[name]['sha256']==expected[name]['sha256'],'Unchanged runtime '+name)
    require(verify(actual['key']).read_bytes()==verify(expected['key']).read_bytes()+b'\nPME-NORM-EXPORT 1\n',
            'Only declared read-only PME export flag appended')


def summarize_panel(parsed,expected,directions,mode):
    """Pure arithmetic over admitted raw records; missing cells never pass."""
    result=dict(bypasses=None,central_polarization_contraction=None,shared_h_checks={},directions={},FD_checks_completed=0)
    if all(k in parsed for k in BYPASSES):
        digests={k:parsed[k]['E_G_virial_token_sha256'] for k in BYPASSES}
        result['bypasses']=dict(token_digests=digests,passed=len(set(digests.values()))==1,
            compared='literal total-energy, full Cartesian gradient and full virial tokens; roles/counters excluded')
    if 'central' in parsed:result['central_polarization_contraction']=parsed['central']['polarization_contraction']
    result['shared_h_checks']={k:v['shared_h'] for k,v in parsed.items() if v['role']=='ACTIVE'}
    if mode=='full':
        require(tuple(x['id'] for x in directions['directions'])==DIRECTION_IDS,'Exact eight fixed directions')
        for direction in directions['directions']:
            name=direction['id'];needed={}
            for step in STEPS:
                keys={sign:[k for k,v in expected.items() if v.get('direction')==name and v.get('step')==step and v.get('sign')==sign]
                      for sign in (-1,1)}
                require(all(len(v)==1 for v in keys.values()),'Exact paired step slots')
                needed[step]=(keys[-1][0],keys[1][0])
            missing=[k for pair in needed.values() for k in pair if k not in parsed]
            if 'central' not in parsed:missing.append('central')
            if missing:result['directions'][name]=dict(status='unavailable',missing=missing,passed=False);continue
            tangent=load_tangent(direction['coefficients_dr_dq'])
            analytic,printing=projection(parsed['central'],direction,tangent)
            energies={s:tuple(parsed[k]['energy_token'] for k in keys) for s,keys in needed.items()}
            result['directions'][name]=dict(status='complete',**slope_checks(name,analytic,printing,energies))
            result['FD_checks_completed']+=3
    required=4 if mode=='central_scout' else 36
    result['complete']=len(parsed)==required and all(k in parsed for k in expected)
    result['passed']=bool(result['complete'] and result['bypasses'] and result['bypasses']['passed'] and
        result['central_polarization_contraction'] and result['central_polarization_contraction']['passed'] and
        all(v is not None and v['passed'] for v in result['shared_h_checks'].values()) and
        (mode=='central_scout' or result['FD_checks_completed']==24 and all(v['passed'] for v in result['directions'].values())))
    return result


def admit_certificate(item,state_id,expected,cfg,authorities):
    certificate=read(item)
    producer=cfg['central_certificate_source'] if state_id=='central' else cfg['certificate_source']
    require(certificate['implementation']==producer,'Certificate producer source pin');verify(producer)
    require(certificate['status']=='stage2_state_nominal_certificate_passed' and
            certificate['sufficient_numerical_test_passed'] is True,'Passed state-specific nominal certificate')
    require(certificate['state_id']==state_id and certificate['candidate']==authorities['candidate'] and
            certificate['source_manifest']==authorities['source_manifest'],'Certificate exact state/source/candidate')
    runtime_inputs(certificate['inputs'],expected['inputs']);verify(certificate['manifest'])
    verify(certificate['expected_state']);verify(certificate['pme_export'])
    nominal=certificate['nominal']
    require(nominal['e_app'] is None and nominal['application_error_bound_established'] is False,'Nominal-only certificate scope')
    require(all(isinstance(nominal[k],(float,int)) and math.isfinite(nominal[k]) and nominal[k]>=0
                for k in ('U_B','U_R','delta_odd','U_nominal','tau')) and nominal['tau']>0,'Finite nominal certificate values')
    for name in ('assembly','factor'):verify(certificate[name])
    if state_id=='central':
        require(certificate['initializer'] is None and certificate['root_admission']==cfg['stage1_actual_admission'] and
                certificate['factor_root_admission']==authorities['central_factor_admission'],'Exact inherited central certificate')
        verify(certificate['factor_root_admission']);initializer=None
    else:
        require(certificate['build']==cfg['build'],'Displaced certificate build')
        initializer=read(certificate['initializer'])
        require(initializer['status']=='stage2_initialization_complete' and initializer['state_id']==state_id,
                'Complete own displaced initialization')
        require(initializer['state_output']==certificate['expected_state'] and initializer['pme_export']==certificate['pme_export'],
                'Initializer actual state and PME certificate join')
        runtime_inputs(initializer['inputs'],expected['inputs'])
    return certificate,initializer


def collect(config_path,result):
    cfg=json.loads(Path(config_path).read_text());a=read(cfg['authorities'])
    require(cfg['status']=='bound_stage2_saved_collection' and cfg['implementation']==pin(__file__),'Bound collector/source')
    for item in a.values():verify(item)
    for name in ('execution_source','build','root_release','execution_manifest'):verify(cfg[name])
    stage1=read(cfg['stage1_actual_admission'])
    require(stage1['status']=='admitted_finite_periodic_operator_parity' and stage1['sampled_operator_parity_passed'] is True and
            stage1['candidate']==a['candidate'],'Admitted actual Stage1')
    displaced=read(a['displaced_manifest']);directions=read(a['directions']);bypass=read(a['bypass_manifest'])
    require(displaced['contract']==a['numerical_contract'] and displaced['directions']==a['directions'] and
            displaced['candidate']==a['candidate'],'Unchanged contract/directions/candidate')
    expected=expected_states(displaced,bypass,cfg['mode']);ids=list(expected)
    require(set(cfg['gradient_results'])==set(ids),'Exact declared gradient denominator')
    displaced_ids=[r['id'] for r in displaced['rows']] if cfg['mode']=='full' else []
    require(set(cfg['certificates'])==set(displaced_ids),'Exact declared certificate denominator')
    result.update(mode=cfg['mode'],configuration=pin(config_path),authorities=cfg['authorities'],candidate=a['candidate'],
        root_release=cfg['root_release'],execution_source=cfg['execution_source'],build=cfg['build'],
        stage1_actual_admission=cfg['stage1_actual_admission'],rows={},certificates={},failures=[],counts={},
        central_certificate=cfg['central_certificate'],scheduler_accounting=cfg.get('scheduler_accounting'))
    parsed={};certificates={};initializers={};native_counts=[]
    for state_id in ['central']+displaced_ids:
        item=cfg['central_certificate'] if state_id=='central' else cfg['certificates'][state_id]
        if item is None:
            result['certificates'][state_id]=dict(status='unavailable',result=None);continue
        try:
            saved=read(item)
            result['certificates'][state_id]=dict(status='pending_admission',result=item,producer_status=saved.get('status'),
                producer_error=saved.get('error'),cost=saved.get('cost'),commands=saved.get('commands'),
                initializer=saved.get('initializer'),assembly=saved.get('assembly'),factor=saved.get('factor'))
            certificate,initializer=admit_certificate(item,state_id,expected[state_id],cfg,a)
            require(certificate['manifest']==cfg['execution_manifest'],'Exact prepared execution inventory')
            if state_id=='central':require(certificate['expected_state']==stage1['expected_central_state'] and
                certificate['pme_export']==stage1['initialized_PME'],'Actual Stage1 central state reused exactly')
            certificates[state_id]=certificate
            if initializer is not None:initializers[state_id]=initializer
            result['certificates'][state_id]=dict(status='admitted_saved_certificate',result=item,
                expected_state=certificate['expected_state'],nominal=certificate['nominal'],
                initializer=certificate['initializer'],assembly=certificate['assembly'],factor=certificate['factor'],
                cost=certificate.get('cost'),commands=certificate.get('commands'))
        except Exception as error:
            result['certificates'].setdefault(state_id,dict(result=item)).update(status='failed_admission',error=repr(error))
            result['failures'].append(dict(stage='certificate',state_id=state_id,error=repr(error)))
    for state_id in ids:
        item=cfg['gradient_results'][state_id]
        if item is None:
            result['rows'][state_id]=dict(status='unavailable',result=None);continue
        row=dict(result=item,status='failed_admission');result['rows'][state_id]=row
        try:
            native=read(item);row.update(producer_status=native.get('status'),producer_error=native.get('error'),
                recorded_counts=native.get('counts'),cost=native.get('cost'),commands=native.get('commands'))
            counts=native.get('counts')
            if isinstance(counts,dict) and type(counts.get('native_programs')) is int:native_counts.append(counts['native_programs'])
            require(native['status']=='stage2_gradient_complete','Native gradient incomplete/failed')
            require(native['state_id']==state_id and native['role']==expected[state_id]['role'],'Native declared slot/role')
            require(native['source_manifest']==a['source_manifest'] and native['build']==cfg['build'],'Native source/build pins')
            for name in ('configuration','release'):verify(native[name])
            require(native['implementation']==cfg['execution_source'],'Native execution adapter')
            same_inputs(native['source_inputs'],expected[state_id]['inputs'])
            runtime_inputs(native['inputs'],expected[state_id]['inputs'])
            require(counts==dict(native_programs=1,gradient_calls=1,energy_only_calls=0,extra_operator_applications=0),
                    'Single native gradient program; no extra diagnostics')
            command=native['commands']['native'];require(command['returncode']==0,'Native process successful')
            raw=verify(native['gradient_output']);value=parse_gradient(raw)
            require(value['state_id']==state_id and value['role']==native['role'] and value['input_binding']==native['input_binding'],
                    'Raw identity/binding matches released record')
            if native['role']=='ACTIVE':
                require(state_id in certificates,'Own state certificate admitted')
                certificate=certificates[state_id]
                expected_pin=cfg['central_certificate'] if state_id=='central' else cfg['certificates'][state_id]
                require(native['certificate']==expected_pin and native['expected_state']==certificate['expected_state'],
                        'Actual solve uses own admitted expected state')
                verify(native['actual_state'])
                require(native['actual_state']['sha256']==certificate['expected_state']['sha256'],'Exact native post-direct state gate')
            else:require(native['certificate'] is None,'Bypass has no active certificate')
            parsed[state_id]=value
            row.update(status='admitted_saved_gradient',gradient_output=native['gradient_output'],raw_identity=dict(
                role=value['role'],state_id=value['state_id'],input_binding=value['input_binding']),
                energy_token=value['energy_token'],polarization_token=value['polarization_token'],
                E_G_virial_token_sha256=value['E_G_virial_token_sha256'],counters=value['counters'],
                PCG_pre_peek=value['PCG_pre_peek'],shared_h=value['shared_h'])
        except Exception as error:
            row['error']=repr(error);result['failures'].append(dict(stage='gradient',state_id=state_id,error=repr(error)))
    result['numerical']=summarize_panel(parsed,expected,directions,cfg['mode'])
    init_counts=[v.get('counts',{}).get('native_programs') for v in initializers.values()]
    require(all(type(x) is int and x==1 for x in init_counts),'One native initializer per admitted displaced certificate')
    result['counts']=dict(expected_gradient_programs=len(ids),admitted_gradients=len(parsed),
        recorded_gradient_native_programs=sum(native_counts),gradient_receipts_with_count=len(native_counts),
        expected_new_initializers=len(displaced_ids),admitted_new_initializers=len(initializers),
        recorded_initializer_native_programs=sum(init_counts),expected_new_certificates=len(displaced_ids),
        admitted_new_certificates=sum(k!='central' for k in certificates),central_certificate_reused='central' in certificates,
        known_total_new_native_programs=sum(native_counts)+sum(init_counts),FD_checks_completed=result['numerical']['FD_checks_completed'],
        measured_clone_induced_solves=sum(v['counters']['induced_solves'] for v in parsed.values() if v['counters']),
        measured_clone_direct_initializations=sum(v['counters']['direct_initializations'] for v in parsed.values() if v['counters']),
        measured_clone_mutual_applications=sum(v['counters']['mutual_applications'] for v in parsed.values() if v['counters']),
        old_library_internal_counters='unavailable; one gradient/source induction path, not fabricated zero',
        additional_collector_native_calls=0,additional_collector_operator_applications=0,additional_collector_factorizations=0)
    result['initializer_receipts']={k:dict(counts=v.get('counts'),cost=v.get('cost'),commands=v.get('commands')) for k,v in initializers.items()}
    if cfg.get('scheduler_accounting') is not None:result['saved_scheduler_accounting']=read(cfg['scheduler_accounting'])
    complete=len(parsed)==len(ids) and len(certificates)==1+len(displaced_ids) and not result['failures']
    passed=complete and result['numerical']['passed']
    result.update(complete=complete,numerical_checks_passed=passed,
        status=('stage2_central_bypass_scout_passed' if cfg['mode']=='central_scout' else 'stage2_finite_derivative_qualification_passed') if passed
        else ('stage2_saved_collection_incomplete' if not complete else 'stage2_finite_numerical_qualification_failed'))
    return result


def main(config,out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result=dict(status='stage2_saved_collection_failed',implementation=pin(__file__),configuration=pin(config),
                new_native_calls=0,new_FFT_replays=0,new_factorizations=0,e_app=None,
                physical_accuracy_qualified=False,trajectory_stability_qualified=False,
                limitations=['Finite numerical consistency at exact consumed old e41f La fixture only.',
                             'Only eight declared projections; full Cartesian gradient saved but not individually qualified.',
                             'Pre-peek residuals retained; no extra post-peek residual/operator application.',
                             'Certificate/replay authorities reused, no repeated sparse assembly or factorization.',
                             'Per-task and scheduler cost records retained separately; overlapping accounting not summed.'])
    try:collect(config,result)
    except Exception as error:result.update(status='stage2_saved_collection_failed',error=repr(error))
    result['collector_cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1)
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (out/'REPORT.txt').write_text(result['status']+'\n'+json.dumps({k:result.get(k) for k in
        ('mode','counts','numerical','failures','limitations','collector_cost')},indent=2)+'\n')
    print(result['status']);return 0 if result.get('numerical_checks_passed') is True else 1


if __name__=='__main__':sys.exit(main(*sys.argv[1:]))
