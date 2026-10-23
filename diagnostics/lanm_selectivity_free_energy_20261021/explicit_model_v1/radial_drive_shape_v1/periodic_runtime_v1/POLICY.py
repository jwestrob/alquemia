"""Pure campaign/geometry admission policy; callers verify bound raw artifacts."""
import IDENTITY as I

BYPASSES=('BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO')

def model(m):
    x=m['model'];I.coefficients(x['coefficients']);I.radial(x['radial']);I.settings(x['settings'])
    for k in ('candidate','composition'):I.pin_shape(x[k])
    I.pin_shape(m['campaign']);I.pin_shape(m['native_sources'])
    I.require(m['state_schema']==I.SCHEMA,'V2 manifest schema')
    return x

def campaign_identity(m,build):
    x=model(m);I.pin_shape(build)
    return {'campaign':m['campaign'],'candidate':x['candidate'],'composition':x['composition'],
            'coefficients':x['coefficients'],'radial':x['radial'],'settings':x['settings'],
            'source_manifest':m['native_sources'],'build':build,'schema':I.SCHEMA}

def full_identity(m,build,state_id,state_output,pme_export):
    x=model(m);row=m['rows'][state_id];I.cell(row['cell'])
    result={'schema':I.SCHEMA,'state_id':state_id,'cell':row['cell'],'settings':x['settings'],
            'coefficients':x['coefficients'],'radial':x['radial'],'inputs':row['inputs'],
            'source_manifest':m['native_sources'],'build':build,'state_output':state_output,'pme_export':pme_export}
    I.validate_identity(result,result)
    return result

def phase_needs(operation,state_id,role):
    I.require(role in ('ACTIVE',)+BYPASSES,'Declared row role')
    I.require(operation in ('initialize','certificate','parity','gradient'),'Declared phase')
    if operation in ('initialize','certificate','parity'):I.require(role=='ACTIVE','Active-only initialization/certificate/parity')
    if operation=='parity':I.require(state_id=='central','One central parity program')
    changed=state_id!='central' and role=='ACTIVE'
    return {'central_parity':operation=='gradient' or changed,
            'central_scout':changed,'own_certificate':operation in ('parity','gradient') and role=='ACTIVE'}

def parity_scope(admission,m,build):
    I.require(admission['status']=='admitted_periodic_v2_central_operator_parity','New central parity admission')
    I.require(admission['campaign_identity']==campaign_identity(m,build),'Same campaign/model/build/settings parity')
    identity=admission['central_identity']
    I.require(identity['state_id']=='central','Central parity state')
    expected=full_identity(m,build,'central',identity['state_output'],identity['pme_export'])
    I.validate_identity(identity,expected)
    I.require(admission['sampled_operator_parity_passed'] is True,'Actual sampled parity')
    I.pin_shape(admission['result']);I.pin_shape(admission['certificate_admission'])
    return True

def scout_scope(admission,m,build,parity_pin):
    I.require(admission['status']=='admitted_periodic_v2_central_and_bypass_scout','New central scout')
    I.require(admission['campaign_identity']==campaign_identity(m,build),'Same campaign scout')
    I.require(admission['central_parity_admission']==parity_pin,'Exact parity root join')
    I.require(admission['all_three_bypass_tokens_equal'] is True and admission['central_dual_energy_passed'] is True,'Actual scout checks')
    return True

def ticket(operation,role,state_id,binding):
    I.require(len(binding)==64 and all(x in '0123456789abcdef' for x in binding),'Configuration digest')
    I.require(isinstance(state_id,str) and state_id and '\n' not in state_id and '\r' not in state_id,'Literal state id')
    phase_needs(operation,state_id,role)
    actual='INITIALIZE' if operation=='initialize' else role
    if role=='BYPASS_OLD_ABSENT':
        I.require(operation=='gradient','Old baseline gradient only')
        return '\n'.join(['STAGE2_BOUND_INPUT_V1',actual,state_id,binding,'END',''])
    return '\n'.join(['STAGE2_BOUND_INPUT_V2',I.SCHEMA,actual,state_id,binding,'END',''])
