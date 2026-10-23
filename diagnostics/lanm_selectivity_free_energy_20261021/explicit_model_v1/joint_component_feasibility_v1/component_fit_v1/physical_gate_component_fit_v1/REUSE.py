"""Immutable exact-vector metadata index; no numerical or native replay."""
import copy
import hashlib
import json
from pathlib import Path


def require(ok,message):
    if not ok:raise ValueError(message)


def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def read(authority):
    require(pin(authority['path'])==authority,'Immutable authority '+authority['path'])
    return json.loads(Path(authority['path']).read_text())


def build(authorities,parent_configuration,parent_evaluator):
    """Read terminal indexes/reviews once; leave historical vector bodies unopened."""
    old=read(authorities['original_inventory']);root=read(old['root_fit_review'])
    audit=read(root['independent_review']);fit=read(root['fit_result'])
    require(root['status']=='admitted_completed_frozen_fit_provenance_and_numerical_review' and audit['status']=='passed_terminal_joint_fit_review_continuation','Original fit aggregate admission')
    require(audit['fit_result']==root['fit_result']==old['fit_result'] and fit['status']=='bounded_joint_response_fit_complete' and fit['error'] is None,'Complete original fit')
    original_cfg=read(fit['configuration'])
    require(original_cfg['evaluator_config']==parent_configuration and original_cfg['evaluator']==parent_evaluator,'Original scientific evaluator')
    completed=read(authorities['component_fit']);review=read(authorities['component_review'])
    require(completed['status']=='bounded_component_fit_complete' and completed['error'] is None and completed['candidate_frozen'] is True,'Terminal component fit')
    require(review['status']=='passed_component_fit_actual_review' and review['result']==authorities['component_fit'] and review['configuration']==completed['configuration'] and review['implementation']==completed['implementation'] and review['selected']==completed['selected'] and review['selected_native_admission_passed'] is True,'Completed component admission')
    component_cfg=read(completed['configuration'])
    require(component_cfg['parent_configuration']==parent_configuration and component_cfg['parent_evaluator']==parent_evaluator,'Component scientific evaluator')
    index={}
    def add(key,parameters,artifact,kind,source):
        require(digest(parameters)==key and len(parameters)==12,'Exact twelve-vector metadata')
        row={'parameters':parameters,'artifact':artifact,'kind':kind,'source':source}
        if key in index:require(index[key]['parameters']==parameters,'Duplicate vector must be exact')
        else:index[key]=row
    for result,source in ((fit,old['fit_result']),(completed,authorities['component_fit'])):
        for key,row in result['vectors'].items():
            require(row['status']=='complete','No failed historical vector imported')
            add(key,row['parameters'],row['artifact'],'vector',source)
    half_root=read(old['pilot_root_review']);half=read(half_root['result'])
    require(half_root['status']=='admitted_one_fixed_half_native_consumed_development_pilot','Half-pilot actual admission')
    for p in (half_root['native_independent_review'],half_root['comparison_independent_review']):read(p)
    full_root=read(authorities['full_pilot_root']);full=read(full_root['result']);full_review=read(full_root['independent_actual_review'])
    require(full_root['status']=='reviewed_complete_native_tradeoff_physical_gate_failure_preserved' and full['status']=='complete_one_ceiling_tradeoff_native_pilot' and full['error'] is None,'Complete failed-physical full pilot')
    require(full_review['status']=='passed_ceiling_tradeoff_native_actual_review' and full_review['result']==full_root['result'] and full_review['original25_pass'] is False and full_root['physical_gates_passed']==24,'Native admission never erases physical failure')
    for result,source in ((half,half_root['result']),(full,full_root['result'])):
        config=read(result['configuration'])
        require(config['parent_configuration']==parent_configuration and config['parent_evaluator']==parent_evaluator,'Pilot scientific evaluator')
        add(result['candidate']['key'],result['candidate']['parameters'],source,'inline_answer',source)
    # Deterministic source precedence: original fit, component fit, then pilots.
    # JSON key ordering is not used to choose which independently executed cell wins.
    return index


def load(entry):
    value=read(entry['artifact'])
    if entry['kind']=='inline_answer':
        require(value['candidate']['parameters']==entry['parameters'],'Inline exact vector')
        value=copy.deepcopy(value['answer']);value.update(status='complete',parameters=entry['parameters'][:])
    else:require(entry['kind']=='vector','Known authority representation')
    require(value['status']=='complete' and value['parameters']==entry['parameters'],'Complete exact saved vector body')
    require(set(('model','response','score','cells'))<=set(value),'Complete scientific vector contract')
    return value
