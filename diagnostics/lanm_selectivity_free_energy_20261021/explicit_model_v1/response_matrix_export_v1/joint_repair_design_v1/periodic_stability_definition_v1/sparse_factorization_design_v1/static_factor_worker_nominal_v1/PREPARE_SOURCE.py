"""Source-only, literal delta from the preserved static worker; no numeric imports."""
from pathlib import Path
import difflib
import hashlib

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'static_factor_worker_v1/RUN.py'
EXPECTED = '0462442495751f77c6b0d98f1844d297669422c2cd9ac4f3a2a9355bc813b3fe'
assert hashlib.sha256(OLD.read_bytes()).hexdigest() == EXPECTED
old = OLD.read_text()
new = old

def replace(a, b):
    global new
    assert new.count(a) == 1, a
    new = new.replace(a, b)

replace('"""One null-bound static factor worker; no assembly, native engine or optimizer."""',
        '"""One nominal saved-operator factor attempt; no assembly or native engine."""')
replace("def context(config,release_path,out):", '''NOMINAL_SHA = '2aa3ac434f3667302fb5cf70d3adda1444e942c6f839e052efef82305ec2e10d'


def nominal_parameters(inputs, review):
    # e_app is unknown, not zero. This worker implements only the nominal contract.
    for item in (inputs, review):
        assert item['e_app'] is None and item['application_error_bound_established'] is False
    for key in ('U_B','U_R','delta_odd'):
        assert inputs[key] == review[key]
        assert type(inputs[key]) in (int,float) and math.isfinite(inputs[key])
    ub,ur,odd=(float(inputs[key]) for key in ('U_B','U_R','delta_odd'))
    upper=ub+ur
    assert ur>=0 and odd>=0 and math.isfinite(upper) and upper>0
    tau=1e-6*upper+odd
    assert math.isfinite(tau) and tau>0
    return dict(U_B=ub,U_R=ur,delta_odd=odd,U_nominal=upper,tau=tau,
                e_app=None,application_error_bound_established=False,
                scope='nominal_saved_operator_only')


def nominal_provenance(cfg, inputs, review):
    assert cfg['nominal_contract']['sha256']==NOMINAL_SHA
    contract=read(cfg['nominal_contract'])['recommended_contract']
    assert contract['tau_nominal']=='1e-6*(U_B+U_R)+delta_odd'
    assert contract['e_app'] is None and contract['application_error_bound_established'] is False
    for key in ('nominal_contract','pme_result','pme_root_review','pme_bound_review'):
        assert inputs[key]==review[key]==cfg[key]
    assert inputs['metal']==review['metal']==cfg['metal']=='La'
    root=read(cfg['pme_root_review']);actual=read(cfg['pme_result']);independent=read(cfg['pme_bound_review'])
    assert root['status']=='admitted_full_box_operator_initialization_nominal_data'
    assert root['execution']==independent['result']==cfg['pme_result']
    assert root['independent_actual_review']==cfg['pme_bound_review']
    assert actual['status']=='full_box_operator_initialization_complete' and actual['error'] is None
    assert root['composition']==actual['composition']==cfg['composition_manifest']
    assert root['e_app'] is None and root['application_error_bound_established'] is False
    for key,source in (('result','campaign_result'),('configuration','campaign_configuration'),
                       ('selected','candidate'),('actual_review','terminal_fit_review')):
        assert actual['candidate'][key]==cfg[source]
    row=actual['cells']['La'];rr=root['cells']['La'];ir=independent['cells']['La']
    assert row['status']=='complete' and row['site_checks']['passed'] is True
    assert row['export']==rr['export']==ir['export']==cfg['pme_export']
    assert row['inputs']==rr['inputs'] and row['inputs']['system.xyz']['sha256']==inputs['geometry']['sha256']
    norms=row['operator_norms']
    assert norms==rr['operator_norms']==ir['operator_norms']
    assert norms['arithmetic']=='ordinary_double_not_outward'
    assert norms['q_negative_count']==0 and norms['q_min']>=0
    for key in ('q_min','q_max','q_abs_max','q_conjugate_difference_max','U_R_nominal','evenness_perturbation_nominal'):
        assert math.isfinite(norms[key]) and norms[key]>=0
    assert inputs['U_R']==norms['U_R_nominal'] and inputs['delta_odd']==norms['evenness_perturbation_nominal']
    for key in ('operator_producer','operator_preparation_result'):
        verify(inputs[key]);assert review[key]==inputs[key]
    nominal_parameters(inputs,review)


def context(config,release_path,out):''')
replace("assert all(review[k] is True for k in ('terminal_candidate_reviewed','initialized_PME_PSD_and_bound_reviewed',",
        "nominal_provenance(cfg,inputs,review)\n    assert all(review[k] is True for k in ('terminal_candidate_reviewed','initialized_PME_nominal_data_and_bound_reviewed',")
replace("assert source_review['status']=='passed_single_static_factor_worker_source_review'",
        "assert source_review['status']=='passed_nominal_static_factor_worker_source_review'")
replace("input_manifest=cfg['input_manifest'],input_review=cfg['input_review'],limits=limits,",
        "input_manifest=cfg['input_manifest'],input_review=cfg['input_review'],limits=limits,\n                nominal_contract=cfg['nominal_contract'],nominal_only=True,e_app=None,application_error_bound_established=False,")
replace("ub,ur=float(inputs['U_B']),float(inputs['U_R']);tau=1e-6*(ub+ur)",
        "nominal=nominal_parameters(inputs,read(cfg['input_review']));tau=nominal['tau']")
replace("result['operator']=dict(shape=[N,N],nnz=NNZ,U_B=ub,U_R=ur,tau=tau,",
        "result['operator']=dict(shape=[N,N],nnz=NNZ,**nominal,")
replace("interpretation='Positive shifted real part is sufficient only; failure does not establish full-PME instability.')",
        "interpretation='Nominal saved-operator sufficient numerical test only; no floating application bound, trajectory or energy/force admission. Failure does not establish full-PME instability.')")
replace("result['status']='passed_static_sufficient_numerical_test' if passed else 'inconclusive_static_arithmetic_or_boundary'",
        "result['status']='passed_nominal_static_sufficient_numerical_test' if passed else 'inconclusive_static_arithmetic_or_boundary'")
for name,text in [('RUN.py',new),('NOMINAL.patch',''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(OLD),tofile=str(HERE/'RUN.py'))))]:
    with (HERE/name).open('x') as handle:handle.write(text)
print(hashlib.sha256(new.encode()).hexdigest())
