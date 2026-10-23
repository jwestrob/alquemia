#!/usr/bin/env python3
"""Only literal central input copies and exact absent/zero response-key deltas."""
import difflib
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
def pin(p):
    p=Path(p).resolve()
    with p.open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(p),'sha256':digest}
def checked(ref):
    if pin(ref['path'])!=ref: raise ValueError('input changed '+ref['path'])
    return Path(ref['path'])
def write(p,value):
    with p.open('x') as f:f.write(value)
def main():
    wall,cpu=time.monotonic(),time.process_time()
    inputs=json.loads((HERE/'SOURCE_INPUTS_v2.json').read_text())
    manifest=json.loads(checked(inputs['displacement_manifest']).read_text())
    center=manifest['central_prepared_input']
    key=checked(center['key']).read_text()
    lines=key.splitlines(True)
    selected=[i for i,x in enumerate(lines) if x.split() and x.split()[0].upper()=='POLAR-LOCAL-RESPONSE']
    if len(selected)!=2 or [lines[i].split()[1:3] for i in selected]!=[['900','402'],['901','402']]:
        raise ValueError('exact response record inventory')
    absent=''.join(x for i,x in enumerate(lines) if i not in selected)
    zero=''.join((' '.join(x.split()[:3])+' 0 0 0 0\n') if i in selected else x for i,x in enumerate(lines))
    output=Path(inputs['source_output']).parent/'bypass_preparation_v1'
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for role,content in [('BYPASS_OLD_ABSENT',absent),('BYPASS_CLONE_ABSENT',absent),('BYPASS_CLONE_ZERO',zero)]:
        out=output/role;out.mkdir()
        for name,field in [('system.xyz','xyz'),('amoebabio09.prm','parameters')]:
            with (out/name).open('xb') as f:f.write(checked(center[field]).read_bytes())
        write(out/'system.key',content)
        write(out/'key.diff',''.join(difflib.unified_diff(key.splitlines(True),content.splitlines(True),fromfile='central/system.prepared.key',tofile=role+'/system.key')))
        rows.append({'role':role,'inputs':{name:pin(out/name) for name in ('system.xyz','system.key','amoebabio09.prm')},'key_diff':pin(out/'key.diff'),'native_output':None,'certificate_required':False,'all_runtime_and_literal_fixture_pins_required':True})
    if rows[0]['inputs']['system.key']['sha256']!=rows[1]['inputs']['system.key']['sha256']:raise ValueError('absent input parity')
    result={'status':'prepared_exact_three_central_bypass_inputs_not_executable','producer':pin(__file__),'source_inputs':pin(HERE/'SOURCE_INPUTS_v2.json'),'central_input':center,'rows':rows,'normal_archive_baseline':inputs['library'],'no_active_certificate_required':True,'build':None,'wrapper':None,'root_release':None,'checks':['exact two 900/901 O402 records','two absent keys identical','zero changes only eight response coefficients','three exact geometry/parameter copies'],'actual':{'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1,'native_calls':0}}
    write(output/'MANIFEST.json',json.dumps(result,indent=2)+'\n')
    print(json.dumps({'manifest':pin(output/'MANIFEST.json'),'actual':result['actual']}))
if __name__=='__main__':main()
