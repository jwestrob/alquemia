"""Metadata-only binder for exactly the original32 changed-state slots."""
import argparse
from contextlib import redirect_stdout
import io
from pathlib import Path
import sys
from types import SimpleNamespace
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from COMMON import *


def ids(m):
    result=[k for k,v in m['rows'].items() if k!='central' and v['role']=='ACTIVE']
    require(len(result)==32 and len(set(result))==32,'exact32 changed-state inventory');return result


def base(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['status']=='bound_stage2_metadata_batch','exact bound metadata source')
    sources=read(cfg['sources']);m=read(sources['manifest']);slots=ids(m)
    require(set(cfg['upstream_results'])==set(cfg['existing_attempts'])==set(slots),'complete32-slot maps, with nulls explicit')
    require(cfg['phase'] in ('initializers','certificates','gradients'),'fixed three stages')
    for p in sources['execution_sources'].values():verify(p)
    require(sources['execution_sources']==m['execution_sources'],'unchanged execution sources')
    cs=read(sources['certificate_source_inputs'])
    require(cs['execution_manifest']==sources['manifest'],'certificate producer recovered manifest')
    for p in cs['certificate_sources'].values():verify(p)
    b=read(cfg['build']);sc=read(cfg['central_scout_review'])
    require(b['status']=='stage2_build_complete' and b['manifest']==sources['manifest'] and b['source_manifest']==m['native_sources'],'exact successful build')
    require(sc['status']=='admitted_stage2_central_and_bypass_scout' and sc['manifest']==sources['manifest'] and sc['build']==cfg['build'],'actual central scout identity')
    require(sc['all_three_bypass_tokens_equal'] is True and sc['central_dual_energy_passed'] is True,'actual central scout gates')
    return sources,m,slots,cs,b


def common_result(item,state,cfg,sources,m,implementation,status):
    r=read(item);row=m['rows'][state]
    require(r['status']==status and r.get('error') is None and r['implementation']==implementation,'complete exact producer')
    require(r['state_id']==state and r['manifest']==sources['manifest'] and r['build']==cfg['build'],'same state/build/manifest')
    require(r['inputs']==row['inputs'] and r['source_inputs']==row['source_inputs'] and r['source_manifest']==m['native_sources'],'literal state/source inputs')
    c=read(r['configuration']);release=read(r['release'])
    require(c['manifest']==sources['manifest'] and c['state_id']==state and c['build']==cfg['build'] and c['central_scout_review']==cfg['central_scout_review'],'actual configured state and scout')
    require(c['implementation']==implementation and release['implementation']==implementation and release['configuration']==r['configuration'],'actual configuration/release source')
    require(Path(release['output']).resolve()==verify(item).parent and release['state_id']==state,'actual released output/state')
    return r,c,release


def admit_initializer(item,state,cfg,sources,m):
    r,c,release=common_result(item,state,cfg,sources,m,m['execution_sources']['NATIVE.py'],'stage2_initialization_complete')
    require(c['status']=='bound_stage2_native_configuration' and c['operation']=='initialize' and release['status']=='released_stage2_native' and release['operation']=='initialize' and release['maximum_native_programs']==1,'one released initializer')
    # Existing finite initializer admission, including saved state/export pins.
    initialization({'initializer':item,'state_id':state,'build':cfg['build'],'manifest':sources['manifest']},m)
    require(r['native_attempts']==1 and set(r['commands'])=={'native'},'one actual native attempt')
    cmd=r['commands']['native'];require(cmd['status']=='complete' and cmd['returncode']==0,'native return code')
    verify(cmd['stdout']);verify(cmd['stderr']);verify(r['initialization_receipt'])
    require(r['PME_checks']['operator_norms']['q_negative_count']==0 and r['PME_checks']['operator_norms']['q_min']>=0,'existing native PME admission')
    return r


def admit_certificate(item,state,cfg,sources,m,cs):
    producer=cs['certificate_sources']['CERTIFY.py']
    r,c,release=common_result(item,state,cfg,sources,m,producer,'stage2_state_nominal_certificate_passed')
    require(c['status']=='bound_stage2_certificate_configuration' and c['certificate_sources']==cs['certificate_sources'],'reviewed exact certificate sources')
    require(release['status']=='released_stage2_certificate' and release['maximum_assemblies']==release['maximum_numeric_factors']==1 and release['maximum_solve_probes']==10,'original finite certificate release')
    require(r['sufficient_numerical_test_passed'] is True and r['native_calls']==0 and r['e_app'] is None and r['application_error_bound_established'] is False,'nominal-only completed certificate')
    require(r['candidate']==m['candidate'] and c['initializer']==r['initializer'],'candidate/initializer join')
    initializer=admit_initializer(r['initializer'],state,cfg,sources,m)
    require(r['expected_state']==initializer['state_output'] and r['pme_export']==initializer['pme_export'],'own initialized state')
    require(set(r['commands'])=={'assembly','factor'},'two exact phase attempts')
    for phase in ('assembly','factor'):
        cmd=r['commands'][phase];require(cmd['status']=='complete' and cmd['returncode']==0,'completed phase process');verify(cmd['stdout']);verify(cmd['stderr'])
    a=read(r['assembly']);f=read(r['factor'])
    for row in (a,f):
        require(row['implementation']==producer and row['configuration']==r['configuration'] and row['release']==r['release'],'phase producer/config/release')
        require(row['state_id']==state and row['manifest']==sources['manifest'] and row['build']==cfg['build'] and row['initializer']==r['initializer'],'phase state lineage')
        require(row['inputs']==r['inputs'] and row['source_inputs']==r['source_inputs'] and row['source_manifest']==m['native_sources'],'phase literal inputs')
        require(row['native_calls']==0 and row.get('error') is None,'no hidden phase native calls/error')
    require(a['status']=='stage2_assembly_complete' and a['assembly_attempts']==1 and a['symbolic_attempts']==a['numeric_attempts']==a['solve_attempts']==0,'exact assembly counts')
    require(a['expected_state']==r['expected_state'] and a['pme_export']==r['pme_export'] and len(a['state_checks'])==4 and all(x['passed'] is True for x in a['state_checks']),'existing state correspondence flags')
    require(f['status']=='passed_nominal_static_sufficient_numerical_test' and f['assembly']==r['assembly'] and f['assembly_attempts']==0 and f['symbolic_attempts']==f['numeric_attempts']==1 and f['solve_attempts']==10,'exact factor/ten-probe counts')
    require(f['sufficient_numerical_test_passed'] is True and f['backend_exception'] is None and not f['backend_warnings'] and not f['audit_warnings'],'original factor admission')
    audit=f['audit'];require(audit['passed'] is True and audit['finite_positive_diagonal'] is True and audit['original_factor_flags_preserved'] is True and len(audit['probes'])==10 and all(x['passed'] is True for x in audit['probes']),'original sampled residual flags')
    require(r['nominal']==f['operator'] and all(r['nominal'][k]==a[k] for k in ('U_B','U_R','delta_odd')),'declared operator joins; no numerical recomputation')
    require(r['nominal']['e_app'] is None and r['nominal']['application_error_bound_established'] is False,'unchanged nominal interpretation')
    return r


def root_acceptance(item,proposal_pin,proposal):
    root=read(item)
    require(root['status']=='accepted_stage2_certificate_batch_admission' and root['proposal']==proposal_pin,'root accepts exact proposed batch')
    for key in ('manifest','build','central_scout_review','accepted_certificates'):require(root[key]==proposal[key],'exact root batch '+key)
    return root


def run(config,out):
    tick=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False)
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config),'rows':{},'error':None,'native_calls':0,'assembly_calls':0,'factor_calls':0,'submissions':0}
    save(out/'RESULT.json',result)
    try:
        cfg=json.loads(Path(config).read_text());sources,m,slots,cs,build=base(cfg)
        result.update(phase=cfg['phase'],manifest=sources['manifest'],build=cfg['build'],central_scout_review=cfg['central_scout_review'])
        accepted={}
        for state in slots:
            row={'upstream_result':cfg['upstream_results'][state],'existing_attempt':cfg['existing_attempts'][state],'configuration':None,'release_template':None,'error':None};result['rows'][state]=row
            try:
                existing=row['existing_attempt']
                if existing is not None:
                    saved=read(existing);row.update(producer_status=saved.get('status'),producer_error=saved.get('error'),producer_cost=saved.get('cost'),producer_commands=saved.get('commands'))
                    require(saved.get('state_id')==state and saved.get('manifest')==sources['manifest'],'existing attempt exact slot identity')
                    row['status']='existing_attempt_retained_no_new_binding';continue
                if cfg['phase']=='initializers':
                    require(row['upstream_result'] is None,'initializer phase has no upstream numerical result');row['status']='eligible_for_binding'
                elif row['upstream_result'] is None:row['status']='missing_upstream_result';continue
                else:
                    saved=read(row['upstream_result']);row.update(producer_status=saved.get('status'),producer_error=saved.get('error'),producer_cost=saved.get('cost'),producer_commands=saved.get('commands'))
                    if cfg['phase']=='certificates':admit_initializer(row['upstream_result'],state,cfg,sources,m)
                    else:admit_certificate(row['upstream_result'],state,cfg,sources,m,cs)
                    row['status']='eligible_for_binding';accepted[state]=row['upstream_result']
            except Exception as exc:row.update(status='rejected_upstream_preserved',error=repr(exc))
        root=None;proposal_pin=None
        if cfg['phase']=='gradients':
            proposal={'status':'proposed_stage2_certificate_batch_admission','manifest':sources['manifest'],'build':cfg['build'],'central_scout_review':cfg['central_scout_review'],'accepted_certificates':accepted,'binder':pin(__file__),'configuration':pin(config),'original_slot_count':32,'missing_or_rejected_states':[s for s in slots if s not in accepted]}
            # To avoid a circular acceptance/config hash, an explicit earlier proposal
            # is reused only when its scientific identity and accepted map match.
            if cfg.get('root_batch_admission') is not None:
                actual=read(cfg['root_batch_admission']);proposal_pin=actual['proposal'];previous=read(proposal_pin)
                for k in ('status','manifest','build','central_scout_review','accepted_certificates','binder','original_slot_count','missing_or_rejected_states'):require(previous[k]==proposal[k],'unchanged root-reviewed proposal '+k)
                root=root_acceptance(cfg['root_batch_admission'],proposal_pin,previous)
            save(out/'PROPOSED_ADMISSION.json',proposal);result['proposed_admission']=pin(out/'PROPOSED_ADMISSION.json');result['root_batch_admission']=cfg.get('root_batch_admission')
        native_binder=load(m['execution_sources']['BIND.py'],'stage2_saved_native_binder')
        if cfg['phase']=='certificates':
            sys.path.insert(0,str(Path(cs['certificate_sources']['BIND.py']['path']).parent));cert_binder=load(cs['certificate_sources']['BIND.py'],'stage2_saved_certificate_binder')
        for state,row in result['rows'].items():
            if row['status']!='eligible_for_binding':continue
            if cfg['phase']=='gradients' and root is None:row['status']='metadata_admitted_awaiting_root_batch_acceptance';continue
            target=out/'bindings'/state
            try:
                if cfg['phase']=='certificates':
                    args=SimpleNamespace(sources=verify(sources['certificate_source_inputs']),state_id=state,build=verify(cfg['build']),scout_review=verify(cfg['central_scout_review']),initializer=verify(row['upstream_result']),out=target)
                    with redirect_stdout(io.StringIO()):cert_binder.run(args)
                else:
                    admission=None
                    if cfg['phase']=='gradients':
                        receipt={'status':'accepted_stage2_exact_state_certificate','certificate':row['upstream_result'],'state_id':state,'root_batch_admission':cfg['root_batch_admission'],'root_proposal':proposal_pin,'manifest':sources['manifest'],'build':cfg['build']}
                        admissions=out/'state_admissions';admissions.mkdir(exist_ok=True);admission=admissions/(state+'.json');save(admission,receipt);row['state_admission']=pin(admission)
                    args=SimpleNamespace(manifest=verify(sources['manifest']),kind='native',state_id=state,operation='initialize' if cfg['phase']=='initializers' else 'gradient',build=verify(cfg['build']),scout_review=verify(cfg['central_scout_review']),certificate=None if cfg['phase']=='initializers' else verify(row['upstream_result']),certificate_admission=admission,out=target)
                    with redirect_stdout(io.StringIO()):native_binder.run(args)
                row.update(status='configuration_bound_not_released',configuration=pin(target/'CONFIG.json'),release_template=pin(target/'ROOT_RELEASE_UNBOUND.json'))
            except Exception as exc:row.update(status='binding_failed_preserved',error=repr(exc))
            save(out/'RESULT.json',result)
        counts={status:sum(row['status']==status for row in result['rows'].values()) for status in sorted({row['status'] for row in result['rows'].values()})}
        result['counts']={'original_slots':32,'statuses':counts,'bound_configurations':sum(x['configuration'] is not None for x in result['rows'].values())}
        result['status']='stage2_metadata_batch_prepared' if not any(x.get('error') for x in result['rows'].values()) else 'stage2_metadata_batch_partial_rejections_preserved'
    except Exception as exc:result.update(status='stage2_metadata_batch_failed',error=repr(exc))
    result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
    save(out/'RESULT.json',result);print(json.dumps({'result':pin(out/'RESULT.json'),'status':result['status'],'counts':result.get('counts')}));return 1 if result['status']=='stage2_metadata_batch_failed' else 0

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();raise SystemExit(run(a.config.resolve(),a.out.resolve()))
