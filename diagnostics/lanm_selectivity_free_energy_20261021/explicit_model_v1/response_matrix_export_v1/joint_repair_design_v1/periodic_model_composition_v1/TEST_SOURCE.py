#!/usr/bin/env python3
"""Small synthetic contracts only; never calls prepare or opens a fit outcome."""
import copy
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
def main():
    tick,cpu=time.monotonic(),time.process_time();s=importlib.util.spec_from_file_location('composition_under_test',HERE/'COMPOSE.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    inputs=json.loads((HERE/'INPUTS.json').read_text());domain=m.read(inputs['domain']);checks=[]
    def ck(ok,name):
        if not ok:raise AssertionError(name)
        checks.append(name)
    def reject(fn,name):
        try:fn()
        except (ValueError,KeyError,TypeError):checks.append(name);return
        raise AssertionError('Accepted '+name)
    v=[p['initial'] for p in domain['parameters']];v[:4]=[.1,-.2,.5,-.25];v[4:8]=[-.1,.2,-.5,.25]
    la=m.physical(v,'La',domain);dy=m.physical(v,'Dy',domain)
    ck(la[:4]==[.1,-.2,domain['S']/16,3*domain['S']/16] and dy[:4]==[-.1,.2,-domain['S']/16,-3*domain['S']/16],'Independent box-to-physical coefficient mapping')
    ck(la[-1]==domain['fixed_native_POLPAIR']['La'] and dy[-1]==domain['fixed_native_POLPAIR']['Dy'],'Parent fixed damping preserved')
    for bad in (float('nan'),float('inf'),True,2.):
        q=v[:];q[0]=bad;reject(lambda:m.physical(q,'La',domain),'Reject bad domain value '+str(bad))
    d=copy.deepcopy(domain);d['parameters'][0]['name']='wrong';reject(lambda:m.physical(v,'La',d),'Reject parameter order ambiguity')
    prm='\n'.join(f'atom {typ} {cls} X "fixture" 1 1.0 0' for typ,cls in m.ROLES.items())
    ck(m.namespace(prm)['reserved_metal_namespace_clear'],'Tiny role map accepted')
    for extra in ('atom 900 900 La "collision" 57 138.91 0','atom 999 901 X "class collision" 1 1 0','atom 402 73 O "duplicate" 8 16 2'):
        reject(lambda:m.namespace(prm+'\n'+extra),'Reject type/class/duplicate collision '+extra.split()[1])
    reject(lambda:m.namespace(prm.replace('atom 402 73','atom 402 74')),'Reject O class/H class confusion')
    base='parameters amoebabio09.prm\n# retained comment\natom 900 900 La "mass" 57 138.91 0\natom 901 901 Dy "mass" 66 138.91 0\npolar-eps 1e-8\nvdwpr 73 900 3.5 .2\nvdwpr 73 901 3.5 .2\npolpair 402 900 .28\npolpair 402 901 .28\n'
    key,overlay=m.render(base,v,domain)
    preserved=lambda text:[x for x in text.splitlines() if not x.lower().startswith(('vdwpr ','polpair ','polar-local-response '))]
    ck(preserved(key)==preserved(base),'Every nonpair line including masses unchanged')
    lines=overlay.splitlines();ck(len(lines)==6 and [x for x in lines if x.startswith('POLAR-LOCAL-RESPONSE')]==[f'POLAR-LOCAL-RESPONSE {mt} 402 '+' '.join(f'{x:.17g}' for x in z[:4]) for mt,z in ((900,la),(901,dy))],'Exactly six translated overlay records')
    ck(all(' 73 ' in x for x in lines if x.startswith('vdwpr')) and all(' 402 ' in x for x in lines if x.startswith('polpair')),'VDW classes versus POLPAIR types')
    for bad in (base.replace('vdwpr 73 900 3.5 .2\n',''),base+'vdwpr 900 73 3.5 .2\n',base+'vdwpr 76 900 3 .2\n',base+'POLAR-DONOR-RESPONSE 900 403 .1 .1\n'):
        reject(lambda:m.render(bad,v,domain),'Reject missing/duplicate/salt/donor override')
    def p(name):return {'path':'/synthetic/'+name,'sha256':name}
    sources={k:p(k) for k in ('campaign_driver','campaign_adapter','campaign_core','vacuum_evaluator','parent_configuration')}
    config={'implementation':sources['campaign_driver'],'adapter':sources['campaign_adapter'],'core':sources['campaign_core'],'parent_evaluator':sources['vacuum_evaluator'],'parent_configuration':sources['parent_configuration'],'root_release':p('root')}
    binding={k:p(k) for k in ('result','configuration','selected','actual_review')};token=m.digest(v);vector={'status':'complete','parameters':v,'cells':{str(i):{'status':'complete'} for i in range(46)}}
    result={'status':'bounded_component_fit_complete','mode':'run','error':None,'candidate_frozen':True,'implementation':sources['campaign_driver'],'configuration':binding['configuration'],'configuration_fingerprint':m.campaign_fingerprint(config),'selected_key':token,'selected':binding['selected'],'vectors':{token:{'status':'complete','artifact':binding['selected'],'parameters':v}}}
    review={'status':'passed_component_fit_actual_review','selected_native_admission_passed':True,'result':binding['result'],'configuration':binding['configuration'],'selected':binding['selected'],'implementation':sources['campaign_driver'],'selected_key':token}
    ck(m.candidate_contract(result,config,review,vector,binding,sources)==v,'Synthetic complete terminal/selected46 contract')
    for field,bad in [('status','started'),('status','failed_component_fit'),('candidate_frozen',False),('error','failed'),('selected_key','wrong'),('configuration_fingerprint','wrong')]:
        r=copy.deepcopy(result);r[field]=bad;reject(lambda:m.candidate_contract(r,config,review,vector,binding,sources),'Reject terminal/index/hash '+field+str(bad))
    r=copy.deepcopy(review);r['selected_native_admission_passed']=False;reject(lambda:m.candidate_contract(result,config,r,vector,binding,sources),'Reject unadmitted selected candidate')
    r=copy.deepcopy(review);r['configuration']=p('wrong');reject(lambda:m.candidate_contract(result,config,r,vector,binding,sources),'Reject review/config mismatch')
    q=copy.deepcopy(vector);q['cells']['0']['status']='failed';reject(lambda:m.candidate_contract(result,config,review,q,binding,sources),'Reject incomplete selected vector')
    reject(lambda:m.validate_engine_binding(None,p('manifest')),'Reject unbound engine')
    # Pure synthetic capability checks use temporary in-memory pin readers,
    # not any candidate/native artifact or an engine invocation.
    eb={k:p(k) for k in ('source_manifest','executable','review','result')};composition=p('composition')
    er={'status':'passed_periodic_BIO_response_static_composition_review','source_manifest':eb['source_manifest'],'executable':eb['executable'],'result':eb['result'],'composition':composition,'capabilities':{k:True for k in m.CAPABILITIES},'periodic_settings':{'pinned':True},'image_convention':'declared','stability_convention':'declared'}
    nr={'composition':composition,'numerical_checks_passed':True};original_verify,original_read=m.verify,m.read
    try:
        m.verify=lambda x:Path(x['path']);m.read=lambda x:er if x==eb['review'] else nr
        ck(m.validate_engine_binding(eb,composition)['supported_bound_composition'],'Synthetic complete capability contract only')
        for capability in ('response_900_901_402','periodic_PME','Na405_Cl414_self_consistent_polarization','native_metal_salt_mixing'):
            er['capabilities'][capability]=False;reject(lambda:m.validate_engine_binding(eb,composition),'Reject unsupported '+capability);er['capabilities'][capability]=True
    finally:m.verify,m.read=original_verify,original_read
    dest=HERE/'source_controls_v1';dest.mkdir(exist_ok=False);result={'status':'passed_synthetic_source_controls','implementation':m.pin(HERE/'COMPOSE.py'),'test_source':m.pin(__file__),'checks':checks,'count':len(checks),'actual_candidate_compositions':0,'fit_outcomes_read':False,'native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}};(dest/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'count':len(checks),'cost':result['cost'],'result':m.pin(dest/'RESULT.json')}))

if __name__=='__main__':main()
