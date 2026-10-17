#!/usr/bin/env python3
"""Versioned saved-parent replay: archived whitening T is authoritative.

The two loader bodies below are copied from pinned parent response/pool inputs;
only the archived-metric admission replaces regenerated float-list equality.
No historical source, physical model, evaluator, search or guard is modified.
"""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    p=Path(record['path']);require(pin(p)['sha256']==record['sha256'],'Changed pin: '+str(p));return p


def read(record):return json.loads(verify(record).read_text())


def module(name,record):
    spec=importlib.util.spec_from_file_location(name,verify(record))
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj


def validate_recovery(record):
    evidence=read(record);failed=read(evidence['failed_preflight']);diagnostic=read(evidence['diagnostic_result'])
    verify(evidence['response_implementation']);verify(evidence['failed_implementation'])
    require(failed['slurm_job_id']=='12954' and failed['status']=='execution_failure'
            and failed['failure']=="ValueError('Exact ten-mode terminal lineage differs')"
            and failed['new_composite_evaluations']==0,'Expected precall failure missing')
    require(diagnostic['status']=='complete' and diagnostic['slurm_job_id']=='12956'
            and diagnostic['failed_preflight']==evidence['failed_preflight']
            and diagnostic['failed_frame_source']==evidence['response_implementation']
            and diagnostic['new_molecular_calls']==0 and len(diagnostic['rows'])==4,
            'Actual roundoff diagnostic lineage differs')
    require({(r['source_id'],r['metal']) for r in diagnostic['rows']}==
            {(s,m) for s in ('7CCO','7CCN') for m in ('La','Dy')},'Diagnostic state denominator differs')
    for r in diagnostic['rows']:
        require(r['checks']['coordinate_scaling_T_exact'] is False and
                all(v for k,v in r['checks'].items() if k!='coordinate_scaling_T_exact') and
                r['existing_whitening_tolerance']==1e-10 and
                all(v<=1e-10 for v in r['original_whitening_maximum_residual'].values()),
                'Original failure is not isolated to admitted whitening roundoff')
    return evidence


def archived_metric(start,search,preflight,sid,metal):
    """Retain the exact saved transform and validate its original physical map."""
    row=next(r for r in preflight['rows'] if (r['source_id'],r['metal'],r['arm'])==(sid,metal,'candidate'))
    require(preflight['status']=='passed' and row['pass'] and
            search['active_indices']==start['active']==row['active_indices'] and
            search['coordinate_scaling_T']==row['coordinate_scaling_T'],
            'Archived metric identity/active membership differs')
    saved=np.asarray(search['coordinate_scaling_T'],dtype=float)
    regenerated=np.asarray(start['metric']['coordinate_scaling_T'],dtype=float)
    n=len(start['active'])
    require(saved.shape==regenerated.shape==(n,n) and n==10 and
            np.isfinite(saved).all() and np.isfinite(regenerated).all(),'Invalid metric arrays')
    # This is a float64 representation bound, not a physical/whitening tolerance.
    bound=64*np.finfo(np.float64).eps*max(1.,float(np.abs(saved).max()),float(np.abs(regenerated).max()))
    difference=float(np.abs(saved-regenerated).max())
    require(difference<=bound,'Regenerated metric differs beyond float64 roundoff')
    _,_,jac0,_=start['source']['kin'].evaluate(np.zeros(64))
    matrix=jac0[start['active']].reshape(n,-1).T
    a=matrix@saved;b=matrix@regenerated
    residual={name:float(np.abs(t.T@t-np.eye(n)).max()) for name,t in [('archived',a),('regenerated',b)]}
    require(all(v<=1e-10 for v in residual.values()),'Existing original-source whitening gate failed')
    # Componentwise dot-product error scale includes all ten local operands.
    physical_bound=64*np.finfo(np.float64).eps*np.maximum(1.,np.maximum(np.abs(matrix)@np.abs(saved),np.abs(matrix)@np.abs(regenerated)))
    require(np.all(np.abs(a-b)<=physical_bound),'Archived/regenerated physical basis differs beyond roundoff')
    current=start['jac'][start['active']].reshape(n,-1).T
    ca=current@saved;cb=current@regenerated
    current_bound=64*np.finfo(np.float64).eps*np.maximum(1.,np.maximum(np.abs(current)@np.abs(saved),np.abs(current)@np.abs(regenerated)))
    require(np.all(np.abs(ca-cb)<=current_bound),'Current physical basis differs beyond roundoff')
    value=dict(start);value['transform']=saved.copy();value['metric']=dict(start['metric'],coordinate_scaling_T=saved.tolist())
    return value,{'archived_transform_authoritative':True,'archived_basis_preflight_row':row,
        'maximum_T_difference':difference,'T_float64_roundoff_bound':float(bound),
        'original_whitening_maximum_residual':residual,'existing_whitening_tolerance':1e-10,
        'original_physical_basis_maximum_difference':float(np.abs(a-b).max()),
        'current_physical_basis_maximum_difference':float(np.abs(ca-cb).max()),
        'physical_basis_roundoff_checks_passed':True,'regenerated_transform_used_for_science':False}


def matched_inputs(response,cfg,out,recovery):
    STATES, START_IDS = response.STATES, response.START_IDS
    BASIS_SHA, MATCHED_SHA = response.BASIS_SHA, response.MATCHED_SHA
    archive_saved = response.archive_saved; relative = response.relative
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    require(pin(response.__file__)['sha256']==cfg['implementation']['sha256']==recovery['response_implementation']['sha256'],'Response adapter differs'); verify(cfg['plan'])
    matched=read(cfg['matched_pool_result']); mpcfg=read(matched['config'])
    verify(matched['implementation']);verify(mpcfg['implementation'])
    require(matched['status']=='complete' and matched['slurm_job_id']=='12720'
            and matched['implementation']['sha256']==mpcfg['implementation']['sha256']==MATCHED_SHA
            and len(matched['rows'])==36 and all(r['status']=='complete' for r in matched['rows']), 'Complete frozen matched pool required')
    mp=module('matched_response_pool_parent',matched['implementation']); d=out/'matched_parent_admission'; d.mkdir()
    h,cells,_,scfg,build,binary,old=mp.inputs(mpcfg,d)
    basis_parent=read(old['basis_scout_result']); bcfg=read(basis_parent['config'])
    require(basis_parent['implementation']['sha256']==BASIS_SHA and basis_parent['implementation']==bcfg['implementation'],'Frozen basis helper differs')
    basis=module('matched_response_basis',basis_parent['implementation'])
    (_,range_module,pool,base,kinworld,settings,old_starts,_,_,_,lbfgs,old_cfg)=basis.inputs(verify(basis_parent['config']))
    global np,FALLBACK
    import numpy as np
    FALLBACK=basis.FALLBACK
    response.np=np;response.FALLBACK=FALLBACK
    archived_basis_preflight=read(basis_parent['preflight_result'])
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=scfg['library_directory'])
    for k in ('LD_PRELOAD','PYTHONHOME','PYTHONPATH'):env.pop(k,None)
    for lib in build['linked_libraries']:verify(lib)
    by={(r['candidate_id'],r['metal']):r for r in matched['rows']}
    cb={(r['candidate_id'],r['metal']):c for c,r in cells}
    world={'config':cfg,'matched_pool':matched,'old_pool':old,'h':h,'binary':binary,'build':build,'env':env,
        'basis':basis,'base':base,'range':range_module,'pool':pool,'kinworld':kinworld,'sources':kinworld['sources'],
        'settings':settings,'mcfg':kinworld['mcfg'],'native_helper':kinworld['fd'],'lbfgs':lbfgs,'starts':{},
        'base_cells':{(sid,m):cb[('geometry_00' if sid=='7CCO' else 'geometry_01',m)] for sid,m in STATES},
        'matched_cells':cb,'basis_implementation':basis_parent['implementation'],'matched_implementation':matched['implementation']}
    checks=[]; saved=out/'matched_saved_starts';saved.mkdir()
    from scaffold_accommodation import Geometry, SETTINGS as GEOMETRY_SETTINGS
    for sid,metal in STATES:
        source=world['sources'][sid];old_start=old_starts[(sid,metal,'candidate')];cid=START_IDS[(sid,metal)]
        slot=next(r for r in old['searches'] if r['source_id']==sid and r['metal']==metal and r['arm']=='candidate')
        search_old=read(slot['search']); row=by[(cid,metal)]; c=cb[(cid,metal)]
        old_start,metric_replay=archived_metric(old_start,search_old,archived_basis_preflight,sid,metal)
        require(slot['status']=='admitted' and slot['candidate_id']==cid and slot['point']==search_old['candidate']
                and row['MACE_evidence']['evidence_record']==slot['point'] and search_old['candidate_admitted']
                and search_old['final_geometry']['pass'] and search_old['source_rebased'] is False
                and search_old['active_indices']==old_start['active'] and search_old['coordinate_scaling_T']==old_start['metric']['coordinate_scaling_T'],
                'Exact ten-mode terminal lineage differs')
        q=np.asarray(search_old['full_q'],dtype=float);x=np.asarray(c['x'],dtype=float)
        require(q.shape==(64,) and np.isfinite(q).all(),'Invalid old terminal coordinates')
        active=old_start['active'];require(len(active)==10 and np.count_nonzero(q[[i for i in range(64) if i not in active]])==0,'Mode inventory differs')
        replay,_,jac,_=source['kin'].evaluate(q)
        require(np.allclose(replay,x,rtol=0,atol=1e-12),'Terminal source-map replay differs')
        guard=range_module.range_guard(base,source,active,Geometry,GEOMETRY_SETTINGS,metal)
        origin_guard=guard.check(np.zeros(64),source['x'],source['x']);start_guard=guard.check(q,x,x)
        require(origin_guard['pass'] and start_guard['pass'],'Original/warm geometry admission failed')
        origin_id='geometry_00' if sid=='7CCO' else 'geometry_01'
        origin,op,_=archive_saved(world,source,metal,cb[(origin_id,metal)],by[(origin_id,metal)],saved/(sid+'_'+metal+'_origin'))
        point,pp,gradient=archive_saved(world,source,metal,c,row,saved/(sid+'_'+metal+'_start'))
        start=dict(old_start,search=search_old,search_pin=slot['search'],point=point,point_pin=pp,
            origin=origin,origin_pin=op,x=x,q=q,jac=jac,gradient=gradient,guard=guard,origin_guard=origin_guard,
            start_guard=start_guard,cache={origin['cache_key']:op,point['cache_key']:pp},records=[op,pp],
            old_objective_terminal=slot['point'])
        world['starts'][(sid,metal)]=start
        checks.append({'source_id':sid,'metal':metal,'pass':True,'candidate_id':cid,'previous_search':slot['search'],
            'old_objective_terminal':slot['point'],'matched_start':pp,'matched_origin':op,'initialization_molecular_calls':0,
            'old_termination':search_old['termination'],'active_indices':active,'active_names':search_old['active_names'],
            'coordinate_scaling_T':start['transform'].tolist(),'origin_guard':origin_guard,'start_guard':start_guard,
            'matched_work_from_origin':relative(point['components'],origin['components']),
            'matched_force_diagnostic':basis.projected_diagnostic(jac,gradient,active,source['kin'].modes),
            'old_composite_gradient_reused':False,'fresh_history':True,'archived_metric_replay':metric_replay})
    world['preflight']={'status':'passed','declared_states':4,'molecular_calls':0,'model_loaded':False,'rows':checks}
    return world


def pool_inputs(frozen_pool,cfg,out,recovery_pin):
    recovery=validate_recovery(recovery_pin)
    SOURCES, METALS=frozen_pool.SOURCES,frozen_pool.METALS
    require(pin(frozen_pool.__file__)['sha256'] == cfg['implementation']['sha256'], 'Pool implementation differs')
    require(pin(cfg['plan']['path']) == cfg['plan'], 'Plan differs')
    helper = cfg['response_helper']
    require(pin(helper['path'])['sha256'] == helper['sha256'], 'Response helper differs')
    spec = importlib.util.spec_from_file_location('matched_response_pool_helper', helper['path'])
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    parents = [read(cfg[key]) for key in ('scout_result', 'remaining_result')]
    configs = [read(parent['config']) for parent in parents]
    for parent, pcfg, stage in zip(parents, configs, ('scout', 'remaining')):
        require(parent['status'] in ('complete', 'unavailable') and parent['stage'] == stage,
                'Execution failure cannot supply terminal candidates')
        require(parent['implementation']['sha256'] == pcfg['implementation']['sha256'] == helper['sha256']
                and parent['plan']['sha256'] == pcfg['plan']['sha256'] == cfg['plan']['sha256']
                and parent['matched_pool_result'] == pcfg['matched_pool_result'] == cfg['matched_pool_result'],
                'Parent implementation/plan/model lineage differs')
        for record in (parent['implementation'], pcfg['implementation'], parent['plan'], pcfg['plan']):
            base.verify(record)
        pre = read(parent['preflight_result'])
        require(parent['protocol_id'] == base.PROTOCOL and pre == parent['preflight']
                and pre['status'] == 'passed' and pre['declared_states'] == 4
                and pre['molecular_calls'] == 0 and pre['model_loaded'] is False
                and len(pre['rows']) == 4 and all(r['pass'] for r in pre['rows'])
                and {(r['source_id'], r['metal']) for r in pre['rows']}
                == {(s, m) for s in SOURCES for m in METALS}, 'Parent four-start preflight differs')
        require((parent['status'] == 'complete')
                == all(r['status'] == 'complete' and r['candidate_admitted'] for r in parent['rows']),
                'Parent/row availability differs')
    require(parents[0]['settings'] == parents[1]['settings'], 'Parent search settings differ')
    require(parents[1]['scout_result'] == configs[1]['scout_result'] == cfg['scout_result'], 'Remaining scout gate differs')
    slots = [dict(row, parent=cfg[key]) for key, parent in zip(('scout_result', 'remaining_result'), parents)
             for row in parent['rows']]
    require(len(parents[0]['rows']) == 1 and (slots[0]['source_id'], slots[0]['metal']) == ('7CCO', 'La')
            and len(slots) == 4 and {(r['source_id'], r['metal']) for r in slots}
            == {(s, m) for s in SOURCES for m in METALS}, 'Four-search denominator differs')
    directory = out / 'source_admission'
    directory.mkdir()
    world = matched_inputs(base, configs[0], directory,recovery)
    require(parents[0]['settings'] == world['settings']
            and parents[0]['fallback_rule'] == parents[1]['fallback_rule'] == base.FALLBACK,
            'Frozen search settings/fallback differ')
    prior = read(cfg['matched_pool_result'])
    require(world['matched_pool'] == prior and prior['status'] == 'complete'
            and prior['slurm_job_id'] == '12720' and prior['declared_geometries'] == 18
            and prior['declared_cells'] == len(prior['rows']) == 36
            and all(r['status'] == 'complete' for r in prior['rows']), 'Frozen18/36 pool differs')
    inventory = read(read(prior['config'])['inventory'])
    candidates = copy.deepcopy(inventory['candidates'])
    require(len(candidates) == 18 and [c['candidate_id'] for c in candidates]
            == ['geometry_%02d' % i for i in range(18)], 'Old candidate order differs')
    cells = {(r['candidate_id'], r['metal']): copy.deepcopy(r) for r in prior['rows']}
    require(len(cells) == 36 and set(cells) == {(c['candidate_id'], m) for c in candidates for m in METALS},
            'Old paired matrix differs')
    return base, world, prior, candidates, cells, slots, parents

