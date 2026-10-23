#!/usr/bin/env python3
"""Bounded Stage1 source/metadata review; no imported producer or numerical code."""
import ast
import datetime
import hashlib
import json
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
CHECKS=[]
VERIFIED={}

def require(value,label):
    if not value: raise ValueError(label)
    CHECKS.append(label)

def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as handle: digest=hashlib.file_digest(handle,'sha256').hexdigest()
    return {'path':str(path),'sha256':digest}

def verify(value):
    if value['path'] not in VERIFIED:
        require(pin(value['path'])==value,'pin '+value['path'])
        VERIFIED[value['path']]=value
    else: require(VERIFIED[value['path']]==value,'consistent repeated pin '+value['path'])
    return Path(value['path'])

def read(value): return json.loads(verify(value).read_text())

def pins_in(value):
    if isinstance(value,dict):
        if set(value)=={'path','sha256'}: verify(value)
        else:
            for child in value.values(): pins_in(child)
    elif isinstance(value,list):
        for child in value: pins_in(child)

def main():
    tick,cpu=time.monotonic(),time.process_time()
    ready_pin=pin(BASE/'STAGE1_READY.json');ready=read(ready_pin)
    require(ready_pin['sha256']=='8566dec062015e130c58fc4b17a6949e52a6edad1b53c324f391bbe842974dff','exact frozen READY')
    manifest=read(ready['manifest']);prep=read(ready['preparation']);controls=read(ready['controls']);link=read(ready['link_recipe'])
    for obj in (ready,manifest,prep,link): pins_in(obj)
    sources=manifest['sources'];auth=read(manifest['authorities']);static=read(ready['static_input'])
    root=read(auth['static_root_admission']);composition=read(prep['composition']);build=read(auth['parent_build'])
    require(root['job_id']==20053 and root['input_manifest']==ready['static_input'] and root['sufficient_test_passed'] is True,'job20053 actual admission')
    require(prep['candidate']==ready['candidate']==static['candidate']==composition['candidate']['selected'],'exact e41f candidate chain')
    require(composition['selected_parameter_key']=='e41fcd3ddff3b8882350b05298864465223ab8362290b3e13363e34cb71d1800','old candidate only')
    require(static['pme_export']['sha256']=='6689dd39e78801284ac1a6c77870e45d0ccbb345def10baf36926a987cd61a47' and static['pme_result']['sha256']=='1d9338035b6aae842f99f6563f43614bc4a760cb13134a529b94045ce8f8fb9e','exact admitted PME authority')
    require(build['status']=='complete' and build['manifest']==auth['parent_source'] and build['library']==ready['library']==link['library'],'admitted source/build/archive identity')
    provenance=read(sources['SETUP_PROVENANCE.json']);original=verify(provenance['parent_induce']).read_text()
    start=original.index('      subroutine periodicresponsesetup (handled)\n');stop=original.index('      subroutine periodicresponsebackendguard\n',start)
    old=original[start:stop];new=verify(ready['setup']).read_text()
    expected=old.replace('subroutine periodicresponsesetup (handled)','subroutine parityresponsesetup (handled)',1).replace('c     This source has no admission of a real periodic execution.\n      call periodicresponsebackendguard\n','c     Standalone explicit-operator initialization only.\nc     Normal setup and both fatal guards remain unchanged.\n')
    require(new==expected and verify(sources['SETUP_AUTHORITY.f']).read_text()==old,'exact two-edit isolated setup clone')
    for name in ('periodicresponsebackendguard','periodicresponsestability'):
        body=original.split('      subroutine '+name+'\n',1)[1].split('      end\n',1)[0]
        require('      call fatal\n' in body,'original fatal guard '+name)
    harness=verify(ready['implementation']).read_text();calls=re.findall(r'^\s*call\s+(\w+)',harness,re.M|re.I)
    require(not set(map(str.lower,calls))&{'induce','gradient','energy','analysis','epolar','epolar1','epolar3','localresponsesetup'},'no solve/energy/gradient entry')
    require(calls.count('dfield0c')==1 and calls.count('ufield0')==1 and calls.count('reciprocal_apply')==1,'one direct initialization and two bounded action call sites')
    require('do mode=1,3' in harness and 'do col=1,10' in harness and 'other=mod(col,10)+1' in harness,'30 full plus30 reciprocal reference application loop')
    require(ready['execution']=={'native_programs':1,'complete_operator_applications':30,'separate_reciprocal_applications':30,'total_reciprocal_applications':60,'direct_initializations':1,'induced_solves':0,'energies':0,'gradients':0,'factorizations':0},'exact program denominator')
    for name in ('PARITY.f90','PARITY_SETUP.f'):
        text=verify(sources[name]).read_text();width=132 if name.endswith('f90') else 72
        require(all(len(line)<=width for line in text.splitlines() if not(name.endswith('.f') and line[:1].lower()=='c')),'compiler source line width '+name)
    for name in ('PREPARE.py','REPLAY.py','DERIVE_SETUP.py','CONTROLS.py'): ast.parse(verify(sources[name]).read_text())
    setting_tree=ast.parse(verify(prep['settings']).read_text())
    settings=next(ast.literal_eval(node.value) for node in setting_tree.body if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SETTINGS' for t in node.targets))
    original_key=verify(prep['original_inputs']['key']).read_text();actual_key=verify(prep['raw_inputs']['system.key']).read_text()
    require(actual_key==original_key+settings,'literal composed key plus exact admitted settings')
    require(prep['raw_inputs']['system.xyz']['sha256']==static['geometry']['sha256']==prep['original_inputs']['xyz']['sha256'],'same literal La coordinates')
    require(prep['raw_inputs']['amoebabio09.prm']['sha256']==prep['original_inputs']['parameters']['sha256'],'BIO parameter bytes unchanged')
    xyz=verify(prep['raw_inputs']['system.xyz']).read_text().splitlines();state=verify(prep['raw_inputs']['expected_state.dat']).read_text().splitlines()
    require(int(xyz[0].split()[0])==11958 and len(xyz)==11960 and state[0]=='PERIODIC_PARITY_EXPECTED_V1','literal site denominator')
    require(state[1].split()==['11958',format(static['metal_stiffness']['normalized_beta'],'.17g')],'expected H tied to job20053')
    alpha={}
    for text in (verify(prep['raw_inputs']['amoebabio09.prm']).read_text(),actual_key):
        for line in text.splitlines():
            fields=line.split()
            if fields and fields[0].lower()=='polarize': alpha[int(fields[1])]=float(fields[2])
    require(len(state)==11963 and state[-1]=='END','exact expected-state inventory')
    for xline,sline in zip(xyz[2:],state[4:-1]):
        x,s=xline.split(),sline.split();typ=int(x[5])
        if not(int(s[0])==int(x[0]) and int(s[1])==typ and list(map(float,s[2:5]))==list(map(float,x[2:5])) and float(s[5])==alpha[typ]): raise ValueError('literal site join '+x[0])
    CHECKS.append('all11958 atom/type/coordinate/alpha literal expectations')
    require([list(map(float,state[i].split())) for i in (2,3)]==[prep['coefficients'][m] for m in ('900','901')],'both original coefficient records')
    require(controls['status']=='passed_stage1_source_and_manufactured_controls' and controls['count']==32 and all(x['passed'] for x in controls['checks']),'32 existing controls read without rerun')
    require(link['compile_and_link_commands']==1 and link['source_rebuilds_required']==0 and link['source']==[ready['implementation'],ready['setup']],'single external link, no archive rebuild')
    require(all(v=='1' for v in link['runtime_environment'].values()),'one-thread prescribed runtime')
    require(link['standard_expand_compatibility']['node_runtime_admitted'] is False and ready['root_build_release'] is None and ready['root_execution_release'] is None,'runtime/build/execution remain unadmitted')
    report=pin(HERE/'REPORT.md')
    result={'status':'passed_stage1_periodic_operator_source_review','implementation':ready['implementation'],'setup':ready['setup'],'replay':ready['replay'],'source_manifest':ready['manifest'],'stage1_ready':ready_pin,'preparation':ready['preparation'],'controls':ready['controls'],'numerical_contract':ready['numerical_contract'],'candidate':ready['candidate'],'static_input':ready['static_input'],'library':ready['library'],'link_recipe':ready['link_recipe'],'reviewer_implementation':pin(__file__),'report':report,'checks':CHECKS,'verified_pin_count':len(VERIFIED),'literal_sites_joined':11958,'normal_guards_byte_preserved':True,'source_only':True,'controls_rerun':False,'builds':0,'native_calls':0,'operator_applications':0,'factorizations':0,'new_fit_outcomes_read':False,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'authorship':'Planck independently reviewed Singer harness/setup/replay; Planck authored the candidate sparse assembly and earlier static scheduling adapter. Mathematical contract is Arendt-authored.','limits':['Stage1 finite operator source readiness only; no compiled/native result or Stage2 release.','Raw direct RHS plus h readback; common direct h addition during a solve is not numerically tested here.','Root launcher must bind exact pins, admitted static state, executable and actual runtime before the one native program.','No uniform e_app, physical accuracy, new candidate, Dy, force or trajectory admission.'],'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    result['status']='passed_stage1_complete_periodic_operator_source_review'
    result['ready']=ready_pin
    with (HERE/'REVIEW.json').open('x') as handle: handle.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    frozen={'status':result['status'],'review':pin(HERE/'REVIEW.json'),'report':report,'implementation':ready['implementation'],'setup':ready['setup'],'replay':ready['replay'],'source_manifest':ready['manifest'],'stage1_ready':ready_pin,'source_only':True,'builds':0,'native_calls':0}
    with (HERE/'READY.json').open('x') as handle: handle.write(json.dumps(frozen,indent=2)+'\n')
    print(json.dumps({'review':frozen['review'],'ready':pin(HERE/'READY.json'),'checks':len(CHECKS),'pins':len(VERIFIED),'cost':result['cost']}))

if __name__=='__main__': main()
