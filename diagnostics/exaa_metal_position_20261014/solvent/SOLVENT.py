#!/usr/bin/env python3
"""Eight fixed transplants: existing component-v2 SCC and native frozen ALPB."""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/exaa_metal_position_20261014/solvent'
BASE=ROOT/'workspaces/exaa_metal_position_20261014'
PACK=ROOT/'workspaces/exaa_source_packing_20261014/packing_score_v1/CONFIG.json'
PROTOCOL='ExaA_reciprocal_metal_position_frozen_ALPB_v1'
ORDER=('origin','adaptive_La');METALS=('Ca','La');SCOPES=('original','transplanted')
EV=23.06054783061903;HA=627.509474;TERMS=('Gelec','Ghb','Gsasa','Gshift')

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,data):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(output):
    """Metadata-only pinning; all scientific verification remains inside Slurm."""
    sp=BASE/'scoring_v1/CONFIG.json';rp=sp.parent/'execution/FINAL_COLLECTION.json'
    score=read(sp);result=read(rp);old=read(PACK);cells=[]
    if result['config']!=pin(sp) or not result['terminal'] or result['complete_new_cells']!=8:raise ValueError('Eight completed MACE cells required')
    for t in score['cells']:
        row,=[r for r in score['rows'] if r['sample']==t['sample']]
        actual,=[r for r in result['rows'] if r['sample']==t['sample']]
        original=row['original']['matrix'][t['metal']][t['candidate']]
        new=actual['matrix'][t['metal']][t['candidate']]['transplanted']
        if new['status']!='complete' or new['xyz']!=t['xyz']:raise ValueError('MACE/preparation join failed')
        cells.append({**t,'valence_electrons':original['GFN2_valence_electrons'],
            **{k:row[k] for k in ('maps','state_signature','source_preparation')},
            'MACE':new['receipt'],'MACE_eV':new['energy_eV'],'original':original})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'SOLVENT.py');shutil.copyfile(HERE/'RUN.sbatch',out/'RUN.sbatch')
    c={k:old[k] for k in ('runtime','settings','factorial','component','native_manifest','native_build','native_binary','native_adapter','MACE_model')}
    c.update(protocol_id=PROTOCOL,implementation=pin(out/'SOLVENT.py'),wrapper=pin(out/'RUN.sbatch'),
        plan=pin(HERE.parent/'SOLVENT_PLAN.md'),matrix=pin(ROOT/'diagnostics/exaa_source_packing_20261014/MATRIX.py'),
        score_config=pin(sp),score_collection=pin(rp),score_adapter=score['implementation'],
        input_audit=pin(HERE.parent/'solvent_input_audit/EVIDENCE.json'),
        corrected_ORCA_parameter_provenance=pin(ROOT/'params/native_gfn2_parameters_v1/corrected_la.xtb.json'),
        ORCA_correction_policy=pin(ROOT/'params/native_gfn2_parameters_v1/la_shell_three_arrays_v1.json'),
        parameter_role='SCC consumes published standalone parameter file in runtime; native ALPB uses qualified internal factory; corrected ORCA JSON is audit provenance, not an executed input',
        cells=cells,rows=score['rows'],CPUs=7,maximum_SCC_calls=32,maximum_native_component_calls=16,
        new_MACE_calls=0,new_DFT_calls=0,new_optimization=False,reused_original_cells=8)
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def mace(cell,c):
    r=read(cell['MACE']);q=read(r['request']);verify(cell['xyz']);verify(r['forces'])
    if r['status']!='complete' or r['model']!=c['MACE_model'] or r['energy_eV']!=cell['MACE_eV']:raise ValueError('MACE result/model differs')
    if q['xyz']['sha256']!=cell['xyz']['sha256'] or (q['charge'],q['multiplicity'])!=(cell['charge'],cell['multiplicity']):raise ValueError('MACE exact state/geometry differs')
    if (r['state']['atoms'],r['state']['charge'],r['state']['spin_multiplicity'],r['state']['head'])!=(192,cell['charge'],1,'omol'):raise ValueError('MACE batch state differs')
    return {'status':'complete','receipt':cell['MACE'],'energy_eV':r['energy_eV']}

def reusable(c,cell,n,r):
    old=cell['original'];p=old['native_frozen_ALPB_result'];a=read(p);t=a['input'];nm=read(c['native_manifest'])
    if a['status']!='available' or a['manifest']!=c['native_manifest']:raise ValueError('Original native component protocol unavailable')
    for k in ('xyz','charge','multiplicity','atom_count'):
        if t[k]!=old[k]:raise ValueError('Original component/source join differs: '+k)
    for k in ('state_signature','maps','source_preparation'):
        if t[k]!=cell[k]:raise ValueError('Original physical state differs: '+k)
    if t['scientific_input_key']['runtime']!=c['runtime']:raise ValueError('Original SCC runtime differs')
    if not any(u['MACE']==old['MACE'] and u['MACE_eV']==old['MACE_eV'] for u in t['uses']):raise ValueError('Original MACE/native join missing')
    archived=n.task_preflight(nm,t)
    if n.checks(nm,a['water'],archived)!=a['checks'] or not all(v['passed'] for v in a['checks']):raise ValueError('Original native water diagonal failed')
    if not all(v['informative'] for v in a['vacuum_precision']):raise ValueError('Original vacuum charge precision failed')
    for medium in ('vacuum','water'):
        rr=a[medium];stored=read(rr['result'])
        if stored!={k:v for k,v in rr.items() if k!='result'}:raise ValueError('Original native replay record differs')
        for k in ('xyz','charges','execution','cavity','born_matrix_with_HB'):verify(rr[k])
        ex=read(rr['execution']);verify(ex['stdout']);verify(ex['stderr'])
        if rr['charge_source']!=medium or rr['xyz']['sha256']!=t['xyz']['sha256'] or rr['charges']['sha256']!=t['endpoints'][medium]['charges']['sha256']:raise ValueError('Original charge/coordinate identity differs')
    if any(a['water'][k]['sha256']!=a['vacuum'][k]['sha256'] for k in ('cavity','born_matrix_with_HB')):raise ValueError('Original cavity mismatch')
    if a['G_qvac_hartree']!=a['vacuum']['G_Eh'] or a['G_qvac_hartree']!=old['native_frozen_ALPB_hartree']:raise ValueError('Original frozen-charge scalar differs')
    endpoints=t['endpoints'];v,w=[endpoints[q]['task'] for q in ('vacuum','water')]
    vp=verify(endpoints['vacuum']['point_manifest']);wp=verify(endpoints['water']['point_manifest'])
    if vp!=wp:raise ValueError('Original endpoint manifest differs')
    component=r.component(vp.parent,v,w)
    if not component['accepted']:raise ValueError('Original component-v2 qualification failed')
    absolute={q:read(endpoints[q]['point'])['accepted'] for q in ('vacuum','water')}
    return {'status':'available','MACE':old['MACE'],'MACE_eV':old['MACE_eV'],'xyz':old['xyz'],
      'native_result':p,'G_qvac_hartree':a['G_qvac_hartree'],'components_hartree':a['vacuum']['components_Eh'],
      'component_v2':component,'historical_absolute_point_accepted':absolute}

def preflight(c,matrix,f,r,n):
    for k in ('plan','wrapper','score_config','score_collection','score_adapter','input_audit','corrected_ORCA_parameter_provenance','ORCA_correction_policy'):verify(c[k])
    if read(c['ORCA_correction_policy'])['candidate_sha256']!=c['corrected_ORCA_parameter_provenance']['sha256']:raise ValueError('Corrected full parameter provenance differs')
    sc=read(c['score_config']);sr=read(c['score_collection'])
    if sr['config']!=c['score_config'] or not sr['terminal'] or sr['complete_new_cells']!=8 or len(c['cells'])!=8:raise ValueError('Finite completed MACE scope differs')
    scoring=load('_metal_position_scoring_preflight',c['score_adapter']);b=load('_metal_position_cpu_metadata',read(sc['worker_manifest'])['implementation'])
    old_pf=scoring.preflight(sc,b)
    runtime=read(c['runtime']);nm=read(c['native_manifest']);build=read(c['native_build'])
    for k in ('xtb','xtb_parameter_file','tblite_python'):verify(runtime[k])
    for p in runtime['xtb_dynamic']['libraries']+runtime['tblite_distribution']['distribution_files']:verify(p)
    for d in runtime['tblite_dynamic'].values():
        for p in d['libraries']:verify(p)
    for p in nm['implementation'].values():verify(p)
    for p in nm['build']['pins']:verify(p)
    for p in build['linked_libraries']:verify(p)
    verify(build['binary'])
    if build['manifest']!=c['native_manifest'] or build['binary']!=c['native_binary'] or nm['MACE_model']!=c['MACE_model']:raise ValueError('Qualified native build/model differs')
    original={};mr={}
    for cell in c['cells']:
        t,=[t for t in sc['cells'] if t['task_id']==cell['task_id']]
        if any(cell[k]!=t[k] for k in t):raise ValueError('New cell differs from completed MACE input')
        symbols,coords=f.coordinates(cell['xyz'])
        if sum(f.VALENCE[z] for z in symbols)-cell['charge']!=cell['valence_electrons'] or len(symbols)!=192:raise ValueError('Electronic inventory differs')
        mr[cell['task_id']]=mace(cell,c)
        original[cell['task_id']]=reusable(c,cell,n,r)
    return {'status':'passed','original_components':original,'reused_MACE':mr,'saved_geometry_checks':old_pf['rows'],
       'gate':'existing component-v2 plus native water diagonal and informative vacuum precision; historical absolute gate retained separately','new_molecular_calls':0}

def energy_cell(c,cell,s):
    if not s or s['status']!='available':return {'status':'unavailable','reason':s.get('reason') if s else 'not executed'}
    nr=read(s['result']);t=nr['input']
    if t['xyz']!=cell['xyz'] or nr['G_qvac_hartree']!=s['G_qvac_hartree']:raise ValueError('New native result identity differs')
    return {'status':'available','MACE':cell['MACE'],'MACE_eV':cell['MACE_eV'],'xyz':cell['xyz'],
      'native_result':s['result'],'G_qvac_hartree':s['G_qvac_hartree'],'components_hartree':nr['vacuum']['components_Eh'],
      'component_v2':s.get('component_parity'),
      'historical_absolute_point_accepted':{q:read(t['endpoints'][q]['point'])['accepted'] for q in ('vacuum','water')}}

def work(p,o):
    if p['status']!='available' or o['status']!='available':return None
    dm=(p['MACE_eV']-o['MACE_eV'])*EV;ds=(p['G_qvac_hartree']-o['G_qvac_hartree'])*HA
    return {'MACE':dm,'frozen_solvent':ds,'total':dm+ds,
            **{k:(p['components_hartree'][k]-o['components_hartree'][k])*HA for k in TERMS}}
def difference(a,b):return {k:a[k]-b[k] for k in a} if a is not None and b is not None else None

def summarize(c,pf,sr):
    new={r['task_id']:r for r in sr};rows=[]
    for source in c['rows']:
        matrices={s:{z:{} for z in METALS} for s in SCOPES};works={s:{} for s in SCOPES}
        for cell in (t for t in c['cells'] if t['sample']==source['sample']):
            z,q,tid=cell['metal'],cell['candidate'],cell['task_id']
            matrices['original'][z][q]=pf['original_components'][tid]
            matrices['transplanted'][z][q]=energy_cell(c,cell,new.get(tid))
        for s in SCOPES:
            works[s]={z:work(matrices[s][z]['adaptive_La'],matrices[s][z]['origin']) for z in METALS}
        d={s:difference(works[s]['Ca'],works[s]['La']) for s in SCOPES}
        pools={};order=[(s,q) for s in SCOPES for q in ORDER]
        for name,keys in (('original_two',order[:2]),('common_four',order)):
            if not all(matrices[s][z][q]['status']=='available' for s,q in keys for z in METALS):pools[name]={'status':'unavailable'};continue
            chosen={}
            for z in METALS:
                options=[(s,q,work(matrices[s][z][q],matrices['original'][z]['origin'])) for s,q in keys]
                s,q,w=min(options,key=lambda x:x[2]['total'])
                chosen[z]={'placement':s,'candidate':q,'xyz':matrices[s][z][q]['xyz'],'work_from_original_origin_kcal_mol':w}
            pools[name]={'status':'available','selected':chosen,'differential_work_kcal_mol':difference(chosen['Ca']['work_from_original_origin_kcal_mol'],chosen['La']['work_from_original_origin_kcal_mol'])}
        pools['four_minus_two_differential_work_kcal_mol']=difference(pools['common_four'].get('differential_work_kcal_mol'),pools['original_two'].get('differential_work_kcal_mol'))
        rows.append({'sample':source['sample'],'role':source['role'],'source_case_id':source['source_case_id'],'matrix':matrices,
          'donor_work_kcal_mol':works,'D_kcal_mol':d,'delta_D_kcal_mol':difference(d['transplanted'],d['original']),
          'secondary_mathematical_pools':pools})
    lookup={r['sample']:r for r in rows}
    return {'protocol_id':PROTOCOL,'terminal':True,'declared_new_cells':8,'reused_original_cells':8,
       'available_new_cells':sum(r['status']=='available' for r in sr),'rows':rows,
       'target_minus_control_D_kcal_mol':{s:difference(lookup[3]['D_kcal_mol'][s],lookup[0]['D_kcal_mol'][s]) for s in SCOPES},
       'pool_tie_order':[list(k) for k in order],'units':{'stored_MACE':'eV','stored_solvent':'hartree','works':'kcal/mol'},
       'conversions':{'eV':EV,'hartree':HA},'classification':None,'scope':'paired donor response and four-geometry descriptive minima; no historical full-pool equivalence, bands or new fit'}

def execute(path,action):
    if not os.environ.get('SLURM_JOB_ID'):raise ValueError('Actual Slurm required before scientific imports')
    c=read(path)
    if c['implementation']!=pin(__file__):raise ValueError('Wrong frozen adapter')
    matrix=load('_existing_fixed_pool_solvent',c['matrix']);place=matrix.placement()
    base=Path(path).parent/('execution' if action=='execute' else 'preflight_'+os.environ['SLURM_JOB_ID'])
    start=time.monotonic();pf=None;sr=[];receipt={'status':'failed','config':pin(path),'placement':place,'action':action}
    with (Path(path).parent/'SOLVENT.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);base.mkdir(exist_ok=False)
        try:
            sys.path.insert(0,str(ROOT/'scripts'))
            f=load('_metal_position_factorial',c['factorial']);r=load('_metal_position_component_v2',c['component']);n=load('_metal_position_native_replay',c['native_adapter'])
            pf=preflight(c,matrix,f,r,n);write(base/'PREFLIGHT.json',pf)
            if action=='execute':sr=matrix.solvent(c,base,f,r,n,pf['reused_MACE'])
            receipt['status']='complete'
        except BaseException as e:
            receipt.update(reason=str(e),traceback=traceback.format_exc());raise
        finally:
            if pf and action=='execute':
                # Raw endpoint/partial replay evidence is retained on an interrupted stage.
                if not sr:
                    for cell in c['cells']:
                        p=base/'solvent/replay/tasks'/cell['task_id']/'RESULT.json'
                        if p.exists():
                            row=read(p);sr.append({k:row[k] for k in ('task_id','status','reason','G_qvac_hartree')}|{'result':pin(p)})
                write(base/'FINAL_COLLECTION.json',{'config':pin(path),**summarize(c,pf,sr)})
            elapsed=time.monotonic()-start
            receipt.update(elapsed_seconds=elapsed,allocated_CPU_seconds=7*elapsed,
              SCC_started=sum(1 for _ in base.glob('solvent/points/*/*/START.json')),
              SCC_terminal=sum(1 for _ in base.glob('solvent/points/*/*/EXECUTION.json')),
              native_component_calls=sum(1 for _ in base.glob('solvent/replay/tasks/*/*/native_EXECUTION.json')),
              new_MACE_calls=0,new_DFT_calls=0)
            write(base/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--output',type=Path,required=True)
    for name in ('preflight','execute'):
        a=sub.add_parser(name);a.add_argument('--config',type=Path,required=True);a.add_argument('--config-sha',required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.output)
    else:execute(verify({'path':str(a.config.resolve()),'sha256':a.config_sha}),a.action)
