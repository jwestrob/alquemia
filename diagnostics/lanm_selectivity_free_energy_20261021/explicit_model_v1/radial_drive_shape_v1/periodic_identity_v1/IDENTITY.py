"""Scalar V2 identity policy. No native, matrix, FFT, solve or factor routines.

This module does not authorize execution. Future wrappers must additionally bind
release, parser and binary authority and retain ordinary native admission.
"""
import hashlib
import json
import math
import numbers
from pathlib import Path

SCHEMA='PERIODIC_PARITY_STATE_V2'
SETTINGS={'ewaldcut_A':12.,'vdwcut_A':9.,'vdwtaper_A':8.1,
          'aeewald_inverse_A':.4,'apewald_inverse_A':.4,'grid':[64,64,64],
          'order':5,'boundary':'TINFOIL','poleps':1e-8,'politer':100}

def require(condition,message):
    if not condition:raise ValueError(message)

def finite(x):
    require(isinstance(x,numbers.Real) and not isinstance(x,bool) and math.isfinite(x),'Finite scalar required')
    return float(x)

def coefficients(value):
    require(isinstance(value,dict) and set(value)=={'900','901'},'Both coefficient slots required')
    for row in value.values():
        require(isinstance(row,list) and len(row)==4,'Four original coefficients')
        for x in row:finite(x)
    return value

def radial(value):
    require(isinstance(value,dict) and set(value)=={'900','901'},'Both radial slots required')
    for x in value.values():require(abs(finite(x))<=1.,'Radial coefficient domain')
    return value

def cell(value):
    require(isinstance(value,list) and len(value)==3,'Orthogonal cell shape')
    for x in value:require(finite(x)>=44.,'Orthogonal cell lower bound')
    return value

def settings(value, fixed=True):
    require(isinstance(value,dict) and set(value)==set(SETTINGS),'Complete settings identity required')
    for k in ('ewaldcut_A','vdwcut_A','vdwtaper_A','aeewald_inverse_A','apewald_inverse_A','poleps'):
        require(finite(value[k])>0,'Positive setting '+k)
    require(value['boundary']=='TINFOIL','Boundary')
    require(type(value['order']) is int and value['order']==5,'Spline order')
    require(type(value['politer']) is int and value['politer']>0,'Iteration count')
    require(isinstance(value['grid'],list) and len(value['grid'])==3 and all(type(x) is int and x>0 for x in value['grid']),'Grid')
    if fixed:require(value==SETTINGS,'Frozen cutoff12 settings mismatch')
    return value

def active_model(coefs, e):
    coefficients(coefs);radial(e)
    return any(x!=0 for row in coefs.values() for x in row) or any(x!=0 for x in e.values())

def validate_role(role,coefs,e,lractive,lrperiodic):
    require(type(lractive) is bool and type(lrperiodic) is bool,'Native activation flags')
    active=active_model(coefs,e)
    if role in ('ACTIVE','INITIALIZE'):
        require(active and lractive and lrperiodic,'Active model/flags required')
    elif role in ('BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO'):
        require(not active and not lractive and not lrperiodic,'Active model cannot bypass')
    else:raise ValueError('Unknown role')
    return True

def validate_settings_joins(**stages):
    require(set(stages)=={'key','gate','serialized','parser','assembler','certificate'},'All six settings joins required')
    for name,value in stages.items():
        try:settings(value)
        except ValueError as exc:raise ValueError(name+': '+str(exc)) from exc
    return True

def pin_shape(p):
    require(isinstance(p,dict) and set(p)=={'path','sha256'},'Exact pin fields')
    require(isinstance(p['path'],str) and Path(p['path']).is_absolute(),'Absolute pin path')
    require(isinstance(p['sha256'],str) and len(p['sha256'])==64 and all(c in '0123456789abcdef' for c in p['sha256']),'SHA256 pin')
    return p

def verify_pin(p):
    pin_shape(p)
    h=hashlib.sha256()
    with Path(p['path']).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    require(h.hexdigest()==p['sha256'],'Immutable file pin mismatch')
    return Path(p['path'])

def validate_identity(observed,expected):
    keys={'schema','state_id','cell','settings','coefficients','radial','inputs','source_manifest','build','state_output','pme_export'}
    for value in (observed,expected):
        require(isinstance(value,dict) and set(value)==keys,'Complete V2 identity required')
        require(value['schema']==SCHEMA,'V2 state schema required')
        require(isinstance(value['state_id'],str) and bool(value['state_id']),'State id')
        cell(value['cell']);settings(value['settings']);coefficients(value['coefficients']);radial(value['radial'])
        require(set(value['inputs'])=={'xyz','key','parameters'},'Input authority')
        for p in list(value['inputs'].values())+[value[k] for k in ('source_manifest','build','state_output','pme_export')]:pin_shape(p)
    require(observed==expected,'Full state/model identity differs')
    return True

def nominal(value):
    require(value['e_app'] is None and value['application_error_bound_established'] is False,'Nominal-only scope')
    ub,ur,odd=[finite(value[k]) for k in ('U_B','U_R','delta_odd')]
    require(min(ub,ur,odd)>=0,'Nonnegative nominal bounds')
    require(math.isfinite(ub+ur) and value['U_nominal']==ub+ur,'Nominal bound sum')
    require(value['tau']==1e-6*(ub+ur)+odd,'Unchanged nominal shift')
    return True

def _validate_chain_records(expected, records, pins):
    """Pure metadata helper; public bound entry below verifies every raw file pin."""
    statuses={'initializer':'stage2_initialization_complete','assembly':'stage2_assembly_complete',
              'factor':'passed_nominal_static_sufficient_numerical_test',
              'certificate':'stage2_state_nominal_certificate_passed',
              'root_admission':'admitted_periodic_identity_v2_certificate'}
    require(set(records)==set(statuses) and set(pins)==set(statuses),'Fresh full certificate chain required')
    for label,status in statuses.items():
        pin_shape(pins[label]);r=records[label]
        require(r['status']==status,'Terminal status '+label)
        validate_identity(r['identity'],expected)
    init,asm,fac,cert,root=[records[k] for k in statuses]
    require(init['state_output']==expected['state_output'] and init['pme_export']==expected['pme_export'],'Raw initializer joins')
    require(asm['initializer']==pins['initializer'],'Assembly initializer join')
    require(asm['cutoff_A']==12. and asm['ewald_inverse_A']==.4,'Actual assembler arguments')
    require(fac['assembly']==pins['assembly'],'Factor assembly join')
    for k in ('initializer','assembly','factor'):require(cert[k]==pins[k],'Certificate '+k+' join')
    require(cert['expected_state']==expected['state_output'] and cert['pme_export']==expected['pme_export'],'Certificate raw state joins')
    require(cert['sufficient_numerical_test_passed'] is True,'Actual nominal sufficient test')
    nominal(fac['nominal']);nominal(cert['nominal'])
    require(fac['nominal']==cert['nominal'],'Factor/certificate nominal join')
    require(root['certificate']==pins['certificate'] and root['initializer']==pins['initializer'],'Root chain joins')
    # Same branch for central and every displacement. No historical implicit bridge.
    require(root['fresh_same_state_chain'] is True,'No implicit central/legacy reuse')
    return {'status':'validated_periodic_identity_v2_metadata','state_id':expected['state_id'],
            'native_admission_performed':False,'physical_qualification':False}

def validate_bound_chain(expected,pins):
    records={name:json.loads(verify_pin(p).read_text()) for name,p in pins.items()}
    result=_validate_chain_records(expected,records,pins)
    for p in list(expected['inputs'].values())+[expected[k] for k in ('source_manifest','build','state_output','pme_export')]:verify_pin(p)
    return result

def validate_legacy_authority(authority,v2_header):
    """Explicit diagnostic authority, never creates an expected-state artifact."""
    require(authority['legacy_schema']=='PERIODIC_PARITY_STATE_V1','Explicit V1 diagnostic')
    require(authority['zero_radial_proven'] is True,'Legacy radial absence/zero proof required')
    require(set(authority['proof_pins'])=={'state','key','composition','source','admission'},'Complete legacy proof pins')
    for p in authority['proof_pins'].values():verify_pin(p)
    settings(authority['settings'],fixed=False);cell(authority['cell'])
    require(authority['settings']==v2_header['settings'] and authority['cell']==v2_header['cell'],'Legacy comparison requires SAME settings/cell')
    require(radial(v2_header['radial'])=={'900':0.,'901':0.},'Exact zero-only legacy comparison')
    return True
