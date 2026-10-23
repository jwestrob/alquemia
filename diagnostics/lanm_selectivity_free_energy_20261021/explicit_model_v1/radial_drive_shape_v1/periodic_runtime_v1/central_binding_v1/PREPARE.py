#!/usr/bin/env python3
"""One immutable central input binding; no native or scientific-array work."""
import argparse
import ast
from collections import Counter
import copy
import difflib
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import socket
import sys
import time

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
def require(ok,message):
    if not ok:raise ValueError(message)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(set(p)=={'path','sha256'} and pin(p['path'])==p,'Changed input '+str(p));return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def load(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def tag(line):
    words=line.split('#',1)[0].split();return words[0].lower() if words else ''
def strip_response(text):return ''.join(line for line in text.splitlines(keepends=True) if tag(line) not in ('polar-local-response','polar-radial-drive'))

def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(config.read_text())
    require(cfg['implementation']==pin(__file__) and cfg['inputs']==pin(HERE/'INPUTS.json') and cfg['specification']==pin(HERE/'SPEC.txt'),'Exact prepared source')
    require(all(os.environ.get(k)=='1' for k in THREADS),'One-thread environment')
    require(Path(cfg['output']).resolve()==out.resolve() and not out.exists(),'New immutable output')
    out.mkdir(parents=True);checks=[]
    result={'status':'central_binding_started','implementation':pin(__file__),'configuration':pin(config),'error':None,'checks':checks,'counts':dict.fromkeys(('native','build','FFT','assembly','factor','solver','submissions'),0),'build':None,'runtime_release':None}
    save(out/'RESULT.json',result)
    def check(ok,label):require(ok,label);checks.append(label)
    try:
        s=read(cfg['inputs'])
        for p in s.values():verify(p)
        decl=read(s['root_declaration']);check(decl['status']=='declared_periodic_central_input_binding_only','Root declared input-only scope')
        for k,p in decl['parents'].items():check(s[k]==p,'Root parent '+k)
        admission=read(s['fit_root_admission']);candidate=read(s['candidate']);producer=read(admission['producer']);fitcfg=read(s['fit_configuration'])
        check(admission['status']=='admitted_bounded_radial14_development_fit','Terminal fit admission')
        check(admission['selected_artifact']==producer['selected']==s['candidate'],'Exact admitted terminal artifact')
        values=candidate['parameters'];key=digest(values)
        check(key==decl['candidate_key']==producer['selected_key']==admission['selected_key'],'Exact 20605 terminal key; no pending witness')
        check(producer['status']=='bounded_radial14_iterative_fit_complete' and producer['error'] is None and producer['candidate_frozen'] is True,'Complete frozen producer')
        check(producer['configuration']==s['fit_configuration'] and fitcfg['parent_evaluator']==s['vacuum_evaluator'],'Fit configuration/helper joins')
        check(candidate['status']=='complete' and len(values)==14 and len(candidate['cells'])==46 and all(x['status']=='complete' for x in candidate['cells'].values()),'Admitted complete46; no raw replay')
        check(all(type(x) is float and math.isfinite(x) for x in values),'Finite literal14 vector')
        domain=read(s['domain']);parent=read(fitcfg['parent_configuration']);check(parent['domain']==s['domain'],'Exact original12 domain')
        compose=load(s['composition_helper'],'central_compose')
        # Execute only the exact pure native mapping functions, with no NumPy import.
        tree=ast.parse(verify(s['vacuum_evaluator']).read_text());names={'checked','effective'}
        nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
        check({n.name for n in nodes}==names,'Pinned native pure helper extraction')
        env={'require':require,'math':math};exec(compile(ast.Module(body=nodes,type_ignores=[]),s['vacuum_evaluator']['path'],'exec'),env)
        effective={m:compose.physical(values[:12],m,domain) for m in ('La','Dy')}
        check(all(effective[m][:6]==env['effective'](values[:12],m,domain) for m in effective),'Composition/native physical mappings exact')
        baseline=read(s['baseline_manifest']);base_row=baseline['rows']['central'];base_inputs=base_row['inputs']
        for p in base_inputs.values():verify(p)
        mapping=read(s['BIO_mapping']);mapped=next(x for x in mapping['aqueous_endpoints'] if x['id']=='La_aqueous_BIO')
        check(base_inputs['xyz']['sha256']==mapped['files']['xyz']['sha256'],'Exact original La geometry/topology/cell bytes')
        check(base_inputs['parameters']['sha256']==mapped['files']['parameters']['sha256']==mapping['BIO_parameters_preserved_byte_for_byte']['sha256'],'Exact unchanged BIO library bytes')
        check(mapping['common_metal_mass_Da']==138.91 and mapped['formal_charge']==0,'Existing mass convention and neutral composition')
        base=verify(base_inputs['key']).read_text();prm=verify(base_inputs['parameters']).read_text()
        check(strip_response(base)==base,'Old absent baseline has no response records')
        namespace=compose.namespace(prm);check(namespace['reserved_metal_namespace_clear'],'Reserved900/901 namespace')
        active,overlay=compose.render(base,values[:12],domain)
        coefficients={str(t):effective[m][:4] for m,t in [('La',900),('Dy',901)]}
        radial={'900':values[12],'901':values[13]}
        radial_text=''.join(f'POLAR-RADIAL-DRIVE {t} 402 {radial[t]:.17g}\n' for t in ('900','901'))
        active+=radial_text;overlay+=radial_text
        def outside_pairs(text):return ''.join(l for l in strip_response(text).splitlines(keepends=True) if tag(l) not in ('vdwpr','polpair'))
        check(outside_pairs(active)==outside_pairs(base),'Only declared pair/response key records changed')
        old_pair=[l for l in base.splitlines() if tag(l)=='polpair'];new_pair=[l for l in active.splitlines() if tag(l)=='polpair']
        check(old_pair==new_pair,'POLPAIR bytes unchanged')
        check(sum(tag(l)=='vdwpr' for l in active.splitlines())==2,'Exactly two candidate metal-water pair overrides')
        absent=strip_response(active)
        zero=absent+''.join(f'POLAR-LOCAL-RESPONSE {t} 402 0 0 0 0\n' for t in (900,901))+''.join(f'POLAR-RADIAL-DRIVE {t} 402 0\n' for t in (900,901))
        rows={};roles={'central':'ACTIVE','BYPASS_OLD_ABSENT':'BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT':'BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO':'BYPASS_CLONE_ZERO'}
        cell=list(map(float,verify(base_inputs['xyz']).read_text().splitlines()[1].split()))
        check(cell==[49.324]*3+[90.]*3,'Exact saved cell49.324A orthogonal')
        for sid,role in roles.items():
            folder=out/'inputs'/sid;folder.mkdir(parents=True)
            for k,n in [('xyz','system.xyz'),('parameters','amoebabio09.prm')]:shutil.copyfile(verify(base_inputs[k]),folder/n)
            text=active if role=='ACTIVE' else zero if role=='BYPASS_CLONE_ZERO' else absent
            (folder/'system.key').write_text(text)
            inputs={k:pin(folder/n) for k,n in [('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm')]}
            # These prepared literal files are runtime inputs; the prior baseline is separately pinned.
            row=read(s['row_template']);row.update(state_id=sid,role=role,cell=cell[:3],source_inputs=copy.deepcopy(inputs),inputs=inputs,baseline_source_inputs=base_inputs)
            rows[sid]=row
        (out/'active_overlay.key').write_text(overlay)
        (out/'BASELINE_TO_ACTIVE.diff').write_text(''.join(difflib.unified_diff(base.splitlines(keepends=True),active.splitlines(keepends=True),fromfile='cutoff12_old_absent',tofile='terminal20605_active')))
        composition={'status':'prepared_only_not_executable','implementation':pin(__file__),'candidate':s['candidate'],'candidate_key':key,'fit_root_admission':s['fit_root_admission'],'fit_configuration':s['fit_configuration'],'source_inputs':cfg['inputs'],'baseline_manifest':s['baseline_manifest'],'prior_composition':s['old_composition'],'mapping':s['BIO_mapping'],'parameters':values,'effective_a_b_c_d_R_epsilon_POLPAIR':effective,'coefficients':coefficients,'radial':radial,'roles':rows,'overlay':pin(out/'active_overlay.key'),'delta':pin(out/'BASELINE_TO_ACTIVE.diff'),'common_metal_mass_Da':138.91,'Dy_original_mass_Da':162.5,'mass_convention_requires_dynamics_decision':True,'BIO_terms_unchanged':True,'geometry_changed':False,'build':None,'runtime_release':None,'physical_accuracy_qualified':False}
        save(out/'COMPOSITION.json',composition)
        oldmeta=read(s['old_physics_metadata'])
        keep=('type_parameters','mutual_scales','pair_overrides','metal_site_index','sites','water_count','all_alpha_positive')
        meta={k:oldmeta[k] for k in keep};meta.update(status='prepared_radial14_literal_physics_metadata_only',candidate=s['candidate'],composition=pin(out/'COMPOSITION.json'),geometry=rows['central']['inputs']['xyz'],coefficients=coefficients,radial=radial,cell_A=cell[:3],cutoff_A=12.,ewald_inverse_A=.4,initialized_PME=None,H_evaluated=False,old_fixed_metadata_authority=s['old_physics_metadata'])
        save(out/'PHYSICS_METADATA.json',meta)
        m=read(s['manifest_template']);m.update(status='prepared_periodic_v2_central_inputs_not_executable',campaign=s['root_declaration'],rows=rows,physics_metadata=pin(out/'PHYSICS_METADATA.json'))
        m['model'].update(candidate=s['candidate'],composition=pin(out/'COMPOSITION.json'),coefficients=coefficients,radial=radial)
        for k in ('environment','native_runtime','loader_command'):m[k]=baseline[k]
        m.update(build=None,build_admission=None,runtime_release=None,initializer=None,certificate=None,central_parity_admission=None,central_scout_admission=None)
        # Pinned current metadata validators only: no native initialization or arrays.
        sys.path.insert(0,str(verify(s['runtime_METADATA']).parent))
        md=load(s['runtime_METADATA'],'central_metadata');native=load(s['runtime_NATIVE'],'central_native_metadata')
        md.policy.model(m) if hasattr(md,'policy') else None
        sites,observed_cell,observed_meta,observed_coeff,observed_radial,pairs=md.physical_metadata(rows['central'],m)
        check(observed_coeff==coefficients and observed_radial==radial,'Both explicit coefficient/radial slots read back')
        counts=Counter(x['type'] for x in sites)
        check(dict(counts)=={900:1,402:3980,403:7960,405:7,414:10},'Exact11958 site types and neutral1La+7Na-10Cl composition')
        for sid in roles:
            if sid!='central':native.bypass_inputs(rows[sid],m);check(True,'Exact bypass delta '+sid)
        check(md.key_settings(active)==decl['settings']==m['model']['settings'],'Exact declared cutoff12 settings')
        for sid in roles:check(rows[sid]['inputs']['xyz']['sha256']==base_inputs['xyz']['sha256'] and rows[sid]['inputs']['parameters']['sha256']==base_inputs['parameters']['sha256'],'Unchanged xyz/prm '+sid)
        check(not any(k.startswith('numpy') or k.startswith('scipy') for k in sys.modules),'No numerical array library loaded')
        save(out/'MANIFEST.json',m)
        result.update(status='prepared_periodic_central_binding_only',manifest=pin(out/'MANIFEST.json'),composition=pin(out/'COMPOSITION.json'),physics_metadata=pin(out/'PHYSICS_METADATA.json'),candidate=s['candidate'],candidate_key=key,coefficients=coefficients,radial=radial,effective=effective,roles=list(roles),type_counts=dict(counts),settings=m['model']['settings'],unexpected_input_mismatches=[],fresh_V2_chain_required=True,raw_RHS_channels_equated=False,scientific_admission=False)
    except BaseException as exc:
        result.update(status='failed_periodic_central_binding',error=repr(exc));raise
    finally:
        result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'hostname':socket.gethostname(),'PID':os.getpid(),'job_id':None}
        save(out/'RESULT.json',result)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    r=prepare(a.config,a.out);print(json.dumps({'status':r['status'],'checks':len(r['checks']),'manifest':r.get('manifest'),'cost':r['cost']}))
