#!/usr/bin/env python3
"""Read-only review of the exact prepared 121-program input inventory."""
import collections
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np

ROOT = Path.cwd()
SUB = Path('lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/donor_extension_qualification_v1')
D, W = ROOT/'diagnostics'/SUB, ROOT/'workspaces'/SUB
MP = W/'preparation_v1/MANIFEST.json'
seen, checks = {}, {}

def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def check(label, value):
    checks[label]=bool(value)
    if not value: raise ValueError(label)

def pins(x):
    if isinstance(x,dict):
        if 'path' in x and 'sha256' in x:
            p=x['path']
            if p not in seen: seen[p]=pin(p)
            check('pin:'+p,seen[p]['sha256']==x['sha256'])
        else:
            for v in x.values(): pins(v)
    elif isinstance(x,list):
        for v in x:pins(v)

def read(p):return json.loads(Path(p['path']).read_text())
def xyz(p):
    lines=Path(p['path']).read_text().splitlines()
    n=int(lines[0].split()[0]);a=[s.split() for s in lines[1:]]
    check('xyz_order:'+p['path'],len(a)==n and [int(s[0]) for s in a]==list(range(1,n+1)))
    return a,np.array([[float(x) for x in s[2:5]] for s in a])
def meta(a):return [(s[1],int(s[5]),list(map(int,s[6:]))) for s in a]
def rot(axis,angle):
    n=np.array(axis,dtype=float);n/=np.linalg.norm(n)
    k=np.array([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
    return np.eye(3)+math.sin(angle)*k+(1-math.cos(angle))*(k@k)
def near(a,b,tol=2e-14):return np.max(np.abs(np.asarray(a)-np.asarray(b)))<=tol

def main():
    wall,cpu=time.monotonic(),time.process_time()
    check('frozen_manifest',pin(MP)['sha256']=='12bc8b408c01eafd34284ae91137ae09eaea60759b2fd5dbd3f5f285b7781b54')
    m=json.loads(MP.read_text());pins(m)
    rows=m['rows'];ps=m['programs'];dm=read(m['sources']['donor_manifest']);dr=read(m['sources']['donor_result'])
    wm=read(m['sources']['water_manifest']);wr=read(m['sources']['water_result'])
    expected_profiles={'absent':None,'zero':[0.,0.],'drive':[.1,0.],'stiffness':[0.,.1],'combined':[.1,.1]}
    check('profiles',m['profiles']==expected_profiles)
    check('counts',len(rows)==16 and len(ps)==121 and len({p['id'] for p in ps})==121)
    check('stage_counts',dict(collections.Counter(p['stage'] for p in ps))=={'A':24,'B':19,'C':18,'D':32,'E':16,'F':4,'G':8}==m['stage_counts'])
    check('program_counts',dict(collections.Counter(p['program'] for p in ps))=={'analyze':66,'testgrad':55}==m['program_counts'])
    check('caps',m['maximum_native_programs']==121 and m['maximum_internal_FD_energies']==576 and m['maximum_snapshots']==20 and m['maximum_workers']==12)
    check('unbound_execution',all(m[k] is None for k in ('native_source_manifest','build_receipt','executables','root_release','runner_implementation')) and m['new_native_calls']==0)
    centers={metal:rows[metal+'_central'] for metal in ('La','Dy')}
    for metal,r in centers.items():
        c=dm['central'][metal];a=read(r['baseline']['authority'])
        check('selected17621:'+metal,r['input']==c['source']['effective_input'] and r['baseline']['authority']==c['source']['authority'] and r['baseline']['energy']==a['energy'] and r['baseline']['forces']==a['forces'])
        check('selected_key:'+metal,r['input']['key']['sha256']=={'La':'ec7e7644150e9382aa57c77f3341ee8d7b904b797b89828365e3b64ff15f9b08','Dy':'d1b823cbdc6ff210ea2b1ed3b21f8c24e5ba148740016b1be30a89104bd799ca'}[metal])
    geometry_errors={}
    for cid,r in rows.items():
        a,x=xyz(r['input']['xyz']);check('row_metadata:'+cid,meta(a)==list(zip(r['names'],r['types'],r['bonds_1based'])) and x.tolist()==r['positions_A'] and len(a)==r['atoms'])
        if r['kind']=='donor_rotation':
            old=dm['rows'][cid];central=centers[r['metal']]; frame=r['rotation_frame'];cx=np.array(central['positions_A']);origin=np.array(frame['centroid_A'])
            pred=cx.copy();pred[1:]=(cx[1:]-origin)@rot(frame['rotation_unit'],r['q_A']/frame['rotation_length_A']).T+origin
            geometry_errors[cid]=float(np.max(np.abs(pred-x)))
            check('literal_endpoint:'+cid,r['input']==old['input'] and r['q_A']==old['q_A'] and r['baseline']['energy']==dr['rows'][cid]['parsed']['native']['total_energy_kcal_mol'] and r['baseline']['forces'] is None)
            check('rotation_endpoint:'+cid,near(pred,x,1.1e-12))
        if r['kind']=='unaffected_active_water':
            old=wr['rows']['combined__'+cid]
            check('active_water_authority:'+cid,r['input']==old['runtime_input'] and r['baseline']['forces']==old['value']['forces'] and r['baseline']['energy']==old['value']['energy'])
            lines=[s.split() for s in Path(r['input']['key']['path']).read_text().splitlines() if s.strip()]
            coeff=[s for s in lines if s[0].upper()=='POLAR-LOCAL-RESPONSE']
            check('actual_nonzero_water:'+cid,len(coeff)==2 and {int(s[1]) for s in coeff}=={500,501} and all(int(s[2])==36 and list(map(float,s[3:]))==[.05,-.04,.08,-.06] for s in coeff))
        if r['kind']=='unaffected':
            old=wm['rows'][cid]
            check('unaffected_authority:'+cid,r['input']==old['baseline_gradient_input'] and r['baseline']['forces']==old['baseline_gradient']['forces'] and r['baseline']['energy']==old['baseline_gradient']['energy'])
        if r['kind']=='unaffected_literal_deletion':
            pa,px=xyz(centers['La']['input']['xyz'])
            check('acetate7_literal',a[0][2:5]==pa[1][2:5] and all(s[1:6]==t[1:6] and list(map(int,s[6:]))==[int(z)-1 for z in t[6:]] for s,t in zip(a,pa[1:])) and len(a)==7 and np.array_equal(x,px[1:]) and r['baseline'] is None)
    u=m['units'];factor=627.509474063*1439.96455/(332.0637133*.529177210544)
    direction=np.array([.7562812564509392,-.15594410634600103,-.6353897204369424])
    check('field_units',u=={'bohr_A':.529177210544,'hartree_kcal_mol':627.509474063,'electric':332.0637133,'elefield':1439.96455,'debye_per_eA':4.803206802,'MV_cm_per_field_au':factor} and m['field_direction']==direction.tolist() and m['origin_A']==[0,0,0])
    domains=set();counts=collections.Counter();semantic=[];max_transform_error=0.
    for p in ps:
        cid=p['cell'];r=rows[cid];a,x=xyz(p['input']['xyz']);srca,sx=xyz(r['input']['xyz']);tag=p['id'];f=p['field_au'];fd=p['fd_total_separation_A'];tr=p['transform']
        check('source_join:'+tag,p['source_input']==r['input'] and p['input']['parameters']['sha256']==r['input']['parameters']['sha256'] and meta(a)==meta(srca))
        base=Path(r['input']['key']['path']).read_text();text=Path(p['input']['key']['path']).read_text();overlay=p['key_overlay']
        check('exact_key_overlay:'+tag,text==base+('\n'+'\n'.join(overlay)+'\n' if overlay else ''))
        donor=[s.split() for s in overlay if s.startswith('POLAR-DONOR-RESPONSE')];coeff=expected_profiles[p['profile']]
        check('donor_profile:'+tag,(not donor if coeff is None else len(donor)==2 and {int(s[1]) for s in donor}=={500,501} and all(int(s[2])==403 and list(map(float,s[3:]))==coeff for s in donor)))
        allowed={'POLAR-PRINT','SAVE-UINDUCE','POLAR-DONOR-RESPONSE','EXTERNAL-FIELD','RESPONSE-MATRIX-EXPORT'}
        check('only_declared_overlays:'+tag,all(s.split()[0] in allowed for s in overlay))
        fields=[list(map(float,s.split()[1:])) for s in overlay if s.startswith('EXTERNAL-FIELD')]
        check('field:'+tag,not fields if f is None else f in (-.001,-.0005,.0005,.001) and len(fields)==1 and near(fields[0],factor*f*direction,1e-15) and p['stage']=='D')
        if tr is None:
            check('literal_xyz:'+tag,p['input']['xyz']['sha256']==r['input']['xyz']['sha256'])
        else:
            if tr['kind']=='translation':pred=sx+[.37,-.21,.16];check('translation_metadata:'+tag,tr['vector_A']==[.37,-.21,.16])
            elif tr['kind']=='rotation':
                rotation=rot([1,2,-1],.23);pred=sx@rotation.T;check('rotation_metadata:'+tag,near(rotation,tr['matrix']) and tr['origin_A']==[0,0,0])
            else:
                s=np.array([[1,.2,-.3],[.2,-.7,.4],[-.3,.4,.5]]);raw=(sx-sx[0])@s.T;length=np.linalg.norm(raw);tangent=raw/length;pred=sx+tr['q_A']*tangent
                check('strain_metadata:'+tag,tr['S']==s.tolist() and abs(length-tr['normalization_L_A'])<2e-14 and near(tangent,tr['tangent']) and tr['q_A'] in (-1e-4,-5e-5,5e-5,1e-4))
            err=float(np.max(np.abs(x-pred)));max_transform_error=max(max_transform_error,err);check('transformed_xyz:'+tag,err<2e-14)
        expected_arg=['system.xyz','EM' if p['exports'] or f is not None or p['stage']=='E' else 'E'] if p['program']=='analyze' else ['system.xyz','Y','Y',format(fd,'.8f')] if fd is not None else ['system.xyz','Y','N']
        check('argv:'+tag,p['argv_tail']==expected_arg)
        check('FD_ledger:'+tag,(fd in (1e-4,5e-5) and p['fd_each_side_A']==fd/2 and p['internal_FD_energies']==48 and len(a)==8 and p['stage']=='C') if fd is not None else p['internal_FD_energies']==0 and p['fd_each_side_A'] is None)
        export='RESPONSE-MATRIX-EXPORT' in overlay
        check('exports:'+tag,p['exports']==(2 if export else 0) and (not export or p['program']=='analyze' and f is None and p['stage'] in ('A','C')))
        parent=Path(p['input']['xyz']['path']).parent;domains.add(parent)
        check('distinct_input_files:'+tag,all(Path(z['path']).parent==parent for z in p['input'].values()))
        counts['FD']+=p['internal_FD_energies'];counts['snapshots']+=p['exports'];counts['old']+=p['binary_role']=='original_testgrad'
        semantic.append((p['stage'],cid,p['profile'],p['program'],f,fd,None if tr is None else tr['kind'],None if tr is None else tr.get('q_A'),p['binary_role']))
    expected=[]
    def add(stage,cid,profile,program,field=None,fd=None,transform=None,q=None,role=None):expected.append((stage,cid,profile,program,field,fd,transform,q,role or program))
    donor=[k for k,r in rows.items() if r['kind'].startswith('donor')]
    for cid in donor:
        for profile in ('absent','zero'):
            add('A',cid,profile,'analyze')
            if cid.endswith('_central'):add('A',cid,profile,'testgrad')
    for cid in set(rows)-set(donor):
        for profile in ('absent','zero','combined'):add('B',cid,profile,'testgrad')
    add('B','acetate7','absent','testgrad',role='original_testgrad')
    for metal in ('La','Dy'):
        cid=metal+'_central'
        for profile in ('drive','stiffness','combined'):
            for fd in (1e-4,5e-5):add('C',cid,profile,'testgrad',fd=fd)
            add('C',cid,profile,'analyze')
            for f in (-.001,-.0005,.0005,.001):
                add('D',cid,profile,'analyze',field=f)
                if profile=='combined':add('D',cid,profile,'testgrad',field=f)
        for tr in ('translation','rotation'):add('F',cid,'combined','testgrad',transform=tr)
        for q in (-1e-4,-5e-5,5e-5,1e-4):add('G',cid,'combined','analyze',transform='strain',q=q)
    for cid in donor:
        if not cid.endswith('_central'):
            for program in ('analyze','testgrad'):add('E',cid,'combined',program)
    check('complete_semantic_inventory',collections.Counter(expected)==collections.Counter(semantic))
    check('distinct_domains_and_total_counts',len(domains)==121 and counts=={'FD':576,'snapshots':20,'old':1})
    out=W/'independent_inventory_review_v1';out.mkdir(exist_ok=False)
    result={'status':'passed_donor_qualification_prepared_inventory_review','implementation':pin(__file__),'manifest':pin(MP),'plan':m['plan'],'preparation_implementation':m['implementation'],'checks_passed':len(checks),'unique_file_pins_verified':len(seen),'checks':checks,'maximum_rigid_or_strain_coordinate_error_A':max_transform_error,'rotation_endpoint_coordinate_errors_A':geometry_errors,'counts':{'base_fixtures':16,'programs':121,'analyze':66,'testgrad':55,'internal_FD_energies':576,'snapshots':20,'original_binary_programs':1},'blockers':[],'scope':'Actual frozen input inventory only; runner and compiled numerical behavior are separately reviewed. No author suite repeated.','new_native_calls':0,'new_quantum_calls':0,'builds':0,'author_suite_reruns':0,'cost':{'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    with (out/'REVIEW.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(pin(out/'REVIEW.json')));print(json.dumps({k:result[k] for k in ('status','checks_passed','unique_file_pins_verified','counts','cost')}))

if __name__=='__main__':main()
