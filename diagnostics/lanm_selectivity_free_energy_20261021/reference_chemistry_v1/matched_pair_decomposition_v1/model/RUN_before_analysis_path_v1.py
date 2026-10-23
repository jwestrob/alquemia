#!/usr/bin/env python3
"""Literal 34-fragment, energy-only native model adapter."""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import threading
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
GIDS=tuple(f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy'))
IDS=tuple(f'pbe0__aquo{n}__{m}__water_{i:02d}' for n in (8,9) for m in ('La','Dy') for i in range(1,n+1))
POINT=[3.919394849221829,.17662704152555941,3.640166322240097,.15457541153280463,.15863887317162664,.202479823314524]
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
NUMBER=r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def fingerprint(cfg):return hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k not in ('root_release','preparation_result')},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def parse_analyze(path):
    """Analyze subset of pinned RUN_PANEL.parse; no invented gradient record."""
    text=Path(path).read_text()
    require('Tinker' in text and 'Version 26.2' in text,'Native version missing')
    require(not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b',text),'Native failure/nonfinite marker')
    es=re.findall(r'Total Potential Energy\s*:\s*('+NUMBER+')',text)
    require(len(es)==1 and text.count('Energy Component Breakdown :')==1,'Energy/components denominator differs')
    components={}
    for line in text.split('Energy Component Breakdown :',1)[1].splitlines():
        m=re.fullmatch(r'\s*([A-Za-z][A-Za-z /()\-]+?)\s+('+NUMBER+r')\s+(\d+)\s*',line)
        if m:
            name=m[1].strip();require(name not in components,'Duplicate component');components[name]={'energy_kcal_mol':float(m[2].replace('D','E').replace('d','e')),'interactions':int(m[3])}
    energy=float(es[0].replace('D','E').replace('d','e'));interm=re.findall(r'Intermolecular Energy\s*:\s*('+NUMBER+')',text)
    require(len(interm)<=1 and {'Van der Waals','Atomic Multipoles','Polarization'}<=set(components),'Required native components absent')
    require(all(math.isfinite(x) for x in [energy]+[v['energy_kcal_mol'] for v in components.values()]),'Nonfinite energy')
    remainder=energy-sum(v['energy_kcal_mol'] for v in components.values());require(abs(remainder)<=1e-6,'Printed component closure failed')
    inter=float(interm[0].replace('D','E').replace('d','e')) if interm else None
    require(inter is None or math.isfinite(inter),'Nonfinite intermolecular energy')
    return {'total_energy_kcal_mol':energy,'components':components,'native_intermolecular_energy_kcal_mol':inter,
        'printed_component_sum_remainder_kcal_mol':remainder,'analyze_output':pin(path),
        'polarization':{'method':'mutual induced dipoles','declared_polar_eps':1e-8,'no_native_nonconvergence_warning':True,'printed_iteration_residual':None}}

def typed_rows(text):
    lines=text.splitlines();n=int(lines[0].split()[0]);rows=[s.split() for s in lines[1:] if s.strip()]
    require(len(rows)==n and [int(r[0]) for r in rows]==list(range(1,n+1)),'Typed XYZ indices/count differ')
    require(all(len(r)>=6 for r in rows),'Typed XYZ fields absent')
    return rows

def extract(cell,parent,parent_xyz):
    gid=cell['parent_geometry_id'];n=int(gid.split('_')[0][4:]);metal=gid.rsplit('_',1)[1];i=cell['water_index'];idx=[0,3*i-2,3*i-1,3*i]
    require(1<=i<=n and cell['parent_atom_indices']==idx and cell['cell_id']==f'pbe0__aquo{n}__{metal}__water_{i:02d}','Literal fragment identity differs')
    geom=cell['geometry'];symbols=[metal,'O','H','H'];positions=[parent['positions_A'][j] for j in idx]
    require(geom['symbols']==symbols and geom['positions_A']==positions and geom['atoms']==4 and geom['metal']==metal and geom['charge']==3,'Shared geometry differs from literal parent')
    ordinary=verify(geom['xyz']).read_text().splitlines();actual=[s.split() for s in ordinary[2:] if s.strip()]
    require(int(ordinary[0])==4 and len(actual)==4 and [r[0] for r in actual]==symbols and [[float(x) for x in r[1:4]] for r in actual]==positions,'Ordinary shared XYZ differs')
    old=typed_rows(parent_xyz);require(len(old)==1+3*n,'Parent typed count differs');mapping={a:j+1 for j,a in enumerate(idx)};rows=[]
    for j,a in enumerate(idx):
        r=old[a];require([float(x) for x in r[2:5]]==positions[j],'Parent typed coordinates differ')
        require(all(int(b)-1 in mapping for b in r[6:]),'Extracted water has an outside bond')
        rows.append([str(j+1),r[1],*r[2:6],*[str(mapping[int(b)-1]) for b in r[6:]]])
    types=[int(r[5]) for r in rows];bonds=[[int(b)-1 for b in r[6:]] for r in rows]
    require(types==[500 if metal=='La' else 501,36,37,37] and bonds==[[],[2,3],[1],[1]],'Extracted typing/graph differs')
    return '4  literal '+cell['cell_id']+'\n'+'\n'.join(' '.join(r) for r in rows)+'\n',types,bonds

def bare_ion_proof(key,prm,metal,resolver,sources):
    p=resolver.parse_parameters(prm,key);mt=500 if metal=='La' else 501
    require(p['atom'][mt]['valence']==0 and p['metal_multipoles'][mt]==[3.]+[0.]*9 and p['polarize'][mt]['groups']==[],'Ion identity/isolated topology differs')
    for text in (key,prm):
        for line in text.splitlines():
            t=line.split('#',1)[0].strip().lower().split()
            if t:require(t[0] not in ('ewald','a-axis','b-axis','c-axis','alpha','beta','gamma','octahedron','dodecahedron','solvate','solvation','external-field','electric-field','reaction-field') and not t[0].startswith(('restrain','position-restraint')),'Nonvacuum/external term prevents zero-ion convention')
    require('polarization            MUTUAL' in prm,'Mutual polarization context differs')
    return {'status':'analytically_zero','energy_kcal_mol':0.0,'metal':metal,'atoms':1,'bonds':0,
        'distinct_pairs':0,'periodic':False,'external_field':False,'solvation':False,
        'reason':'energy.f initializes every term to zero; no bonded entries; distinct-pair vdW/multipole/polarization sums are empty; no field means no induced dipole. No constant elemental self energy in this Hamiltonian.',
        'source_pins':sources,'missing_data_substitution':False}

def default_config(panel):
    root=read(pin(HERE.parent/'ROOT_INPUTS.json'));data=read(root['consumed_input_data']);ecfg=json.loads((ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1/EVALUATOR_CONFIG.json').read_text())
    src=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/tinker_source/source'
    return {'protocol':'matched_literal_pair_model_v1','implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),
        'root_inputs':pin(HERE.parent/'ROOT_INPUTS.json'),'declaration':root['declaration'],'data':root['consumed_input_data'],'panel':pin(panel),
        'candidate_result':data['provenance']['candidate_result'],'selected_review':data['provenance']['selected_review'],
        'native_helper':ecfg['native_helper'],'resolver':ecfg['resolver'],'analyze':ecfg['executables']['analyze'],
        'physics_sources':{n:pin(src/n) for n in ('energy.f','ehal.f','empole.f','epolar.f','induce.f')},
        'maximum_native_programs':34,'workers':4,'root_release':None,'preparation_result':None}

def checked_config(cfg,threads=False):
    require(cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Model source/spec differs')
    require(cfg['protocol']=='matched_literal_pair_model_v1' and cfg['maximum_native_programs']==34 and cfg['workers']==4,'Finite model scope differs')
    if threads:require(all(os.environ.get(k)=='1' for k in THREADS),'Thread controls differ')
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
    for v in cfg['physics_sources'].values():verify(v)
    root=read(cfg['root_inputs']);require(root['declaration']==cfg['declaration'] and root['consumed_input_data']==cfg['data'],'Root input join differs')
    require(cfg['candidate_result']['sha256']=='5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832' and cfg['selected_review']['sha256']=='fd7ea77167c833c3990d7624ff136042ce03f3ccf9dc4169ba7cff2aa51bb6a2','Frozen selection differs')

def prepare(cfg,out):
    checked_config(cfg,True);data=read(cfg['data']);panel=read(cfg['panel']);review=read(cfg['selected_review']);resolver=module(cfg['resolver'],'pair_model_resolver')
    require(panel['root_inputs']==cfg['root_inputs'],'Shared panel lineage differs')
    require(data['baseline_parameters']==POINT and data['provenance']['candidate_result']==cfg['candidate_result'] and data['provenance']['selected_review']==cfg['selected_review'],'Compact selected lineage differs')
    require(review['status']=='selected_ten_raw_review_passed' and review['result']==cfg['candidate_result'] and review['all_selected_raw_replays_exact'] and review['exact_selected_score_replayed'] and review['all_development_controls_passed'] is False and review['error'] is None,'Raw pass/science failure admission differs')
    require(len(panel['cells'])==36 and sum(c['role']=='bare_ion' for c in panel['cells'])==2,'Shared 36-cell denominator differs')
    cells=[c for c in panel['cells'] if c['role']=='dimer'];require([c['cell_id'] for c in cells]==list(IDS),'Exact 34 literal dimers required')
    clusters={};sources={};ions={}
    for gid in GIDS:
        m=data['baseline_model'][gid];inp=m['provenance']['effective_input'];require(m['provenance']['requested_vector']==POINT,'Selected native vector differs')
        for p in inp.values():verify(p)
        verify(m['provenance']['native_authority']);raw=parse_analyze(verify(m['native']['analyze_output']))
        require(raw['total_energy_kcal_mol']==m['energy'] and raw['components']==m['native']['components'],'Saved cluster energy replay differs')
        clusters[gid]={'energy_kcal_mol':m['energy'],'components':raw['components'],'native':raw,'native_authority':m['provenance']['native_authority'],'input':inp}
        sources[gid]=inp;metal=gid.rsplit('_',1)[1]
        ions[metal]=bare_ion_proof(verify(inp['key']).read_text(),verify(inp['parameters']).read_text(),metal,resolver,cfg['physics_sources'])
    for metal in ('La','Dy'):
        a,b=(sources[f'aquo{n}_central_{metal}'] for n in (8,9))
        require(all(a[k]['sha256']==b[k]['sha256'] for k in ('key','parameters')),'Same-metal selected key/prm differs across coordination')
    rows={}
    for cell in cells:
        cid=cell['cell_id'];gid=cell['parent_geometry_id'];inp=sources[gid];text,types,bonds=extract(cell,data['rows'][gid]['geometry'],verify(inp['xyz']).read_text())
        directory=out/'inputs'/cid;directory.mkdir(parents=True);(directory/'system.xyz').write_text(text)
        for k,n in (('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(inp[k]),directory/n)
        require(re.findall(r'(?im)^\s*parameters\s+(\S+)',(directory/'system.key').read_text())==['amoeba09.prm'],'Local parameter filename differs')
        rows[cid]={'id':cid,'parent_geometry_id':gid,'water_index':cell['water_index'],'parent_atom_indices':cell['parent_atom_indices'],'geometry':cell['geometry'],
            'input':{k:pin(directory/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))},'types':types,'bonds_zero_based':bonds,'selected_source_input':inp}
    manifest={'status':'prepared','implementation':pin(__file__),'config_fingerprint':fingerprint(cfg),'panel':cfg['panel'],'data':cfg['data'],
        'candidate_result':cfg['candidate_result'],'selected_review':cfg['selected_review'],'rows':rows,'saved_cluster_models':clusters,'analytic_bare_ions':ions,
        'maximum_native_programs':34,'workers':4,'new_native_calls':0,'geometry_rotation':False,'fitting':False}
    save(out/'MANIFEST.json',manifest);return manifest

def bound_inputs(cfg,manifest):
    checked_config(cfg);require(manifest['implementation']==cfg['implementation'] and manifest['config_fingerprint']==fingerprint(cfg) and manifest['panel']==cfg['panel'],'Prepared manifest differs')
    require(list(manifest['rows'])==list(IDS) and set(manifest['saved_cluster_models'])==set(GIDS),'Prepared denominator differs')
    panel=read(cfg['panel']);byid={c['cell_id']:c for c in panel['cells']};data=read(cfg['data'])
    for cid,row in manifest['rows'].items():
        inp=row['input'];source=row['selected_source_input'];gid=row['parent_geometry_id'];cell=byid[cid]
        require(row['geometry']==cell['geometry'] and row['parent_atom_indices']==cell['parent_atom_indices'] and source==data['baseline_model'][gid]['provenance']['effective_input'],'Prepared source map differs')
        text,types,bonds=extract(cell,data['rows'][gid]['geometry'],verify(source['xyz']).read_text())
        require(verify(inp['xyz']).read_text()==text and types==row['types'] and bonds==row['bonds_zero_based'],'Prepared literal typed geometry differs')
        for k in ('key','parameters'):require(verify(inp[k]).read_bytes()==verify(source[k]).read_bytes(),'Selected key/prm bytes changed')
    return manifest

def model_cycle(rows,clusters,ions):
    require(set(rows)==set(IDS) and set(clusters)==set(GIDS) and set(ions)=={'La','Dy'},'Cycle denominator differs')
    require(all(v['status']=='analytically_zero' and v['energy_kcal_mol']==0. and v['missing_data_substitution'] is False for v in ions.values()),'Missing or unsupported ion convention')
    if not all(v['status']=='complete' for v in rows.values()):return None
    signs={f'aquo{n}_central_{m}':(1 if (n,m) in ((9,'Dy'),(8,'La')) else -1) for n in (8,9) for m in ('La','Dy')}
    require(all(math.isfinite(v['energy_kcal_mol']) for v in rows.values()) and all(math.isfinite(v['energy_kcal_mol']) for v in clusters.values()),'Nonfinite cycle term')
    C=math.fsum(signs[g]*v['energy_kcal_mol'] for g,v in clusters.items());P=math.fsum(signs[v['parent_geometry_id']]*v['energy_kcal_mol'] for v in rows.values())-(ions['Dy']['energy_kcal_mol']-ions['La']['energy_kcal_mol'])
    return {'units':'kcal/mol','C':C,'P':P,'B':C-P,'C_minus_P_minus_B':C-P-(C-P),'new_target':False}

class Counter:
    def __init__(self):self.lock=threading.RLock();self.stop=threading.Event();self.count=0
    def reserve(self):
        with self.lock:
            require(not self.stop.is_set() and self.count<34,'Global stop or 34-program ceiling');self.count+=1;return self.count

def release_gate(cfg,prep,manifest,runtime=True):
    require(cfg['root_release'] is not None,'Root release required');r=read(cfg['root_release'])
    require(r['status']=='approved_matched_pair_model_energy' and r['implementation']==cfg['implementation'] and r['config_fingerprint']==fingerprint(cfg) and r['preparation_result']==cfg['preparation_result'] and r['panel']==cfg['panel'],'Root release binding differs')
    require(r['maximum_native_programs']==34 and r['workers']==4 and r['fit_authorized'] is False,'Root release scope differs')
    require(prep['status']=='prepared' and prep['new_native_calls']==0 and prep['manifest']==pin(verify(prep['manifest'])),'Preparation not admitted')
    if runtime:require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4,'One four-CPU Slurm task required')

def run(cfg,out):
    start=time.monotonic();cpu=time.process_time();child0=resource.getrusage(resource.RUSAGE_CHILDREN);checked_config(cfg,True)
    prep=read(cfg['preparation_result']);manifest=bound_inputs(cfg,read(prep['manifest']));release_gate(cfg,prep,manifest);native=module(cfg['native_helper'],'literal_pair_native');count=Counter()
    result={'status':'running','implementation':pin(__file__),'config':cfg,'manifest':prep['manifest'],'panel':cfg['panel'],'candidate_result':cfg['candidate_result'],'selected_review':cfg['selected_review'],
        'job_id':os.environ['SLURM_JOB_ID'],'native_program_invocations':0,'global_stop':False,'error':None,'rows':{cid:{'id':cid,'status':'not_attempted'} for cid in IDS},
        'saved_cluster_models':manifest['saved_cluster_models'],'analytic_bare_ions':manifest['analytic_bare_ions'],'model_cycle':None,'fit_targets_changed':False,'new_forces':False}
    def persist():
        with count.lock:result.update(native_program_invocations=count.count,global_stop=count.stop.is_set());save(out/'RESULT.json',result)
    def worker(cid):
        if count.stop.is_set():return
        row=manifest['rows'][cid];directory=out/'cells'/cid;record={'id':cid,'status':'started','parent_geometry_id':row['parent_geometry_id'],'geometry':row['geometry'],'input':row['input'],'receipts':[]}
        try:
            directory.mkdir(parents=True)
            for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['input'][k]),directory/n)
            record['runtime_input']={k:pin(directory/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            with count.lock:result['rows'][cid]={'id':cid,'status':'started'};save(directory/'RESULT.json',record);count.reserve();save(out/'CALL_COUNT.json',{'native_program_invocations':count.count,'maximum':34})
            native.invoke([cfg['analyze']['path'],'system.xyz','E'],directory,directory/'analyze.out',record['receipts'])
            parsed=parse_analyze(directory/'analyze.out')
            with count.lock:record.update(status='complete',energy_kcal_mol=parsed['total_energy_kcal_mol'],components=parsed['components'],native=parsed)
        except Exception as exc:
            count.stop.set()
            with count.lock:record.update(status='failed_or_unavailable',error=repr(exc));result['error']=repr(exc)
        finally:
            with count.lock:
                result['rows'][cid]=record
                if directory.exists():save(directory/'RESULT.json',record)
                persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(worker,IDS))
        require(not count.stop.is_set() and count.count==34 and all(v['status']=='complete' for v in result['rows'].values()),'Incomplete 34-row model denominator')
        result['model_cycle']=model_cycle(result['rows'],manifest['saved_cluster_models'],manifest['analytic_bare_ions']);result['status']='matched_pair_model_complete'
    except Exception as exc:count.stop.set();result.update(status='failed_or_unavailable',error=repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu,'aggregate_child_user_CPU_seconds':ch.ru_utime-child0.ru_utime,'aggregate_child_system_CPU_seconds':ch.ru_stime-child0.ru_stime,'allocated_CPUs':4,'per_receipt_child_CPU_deltas_not_summed':True};persist()
    return result

def admit_result(cfg,manifest,result):
    """Saved-only exact 34-row raw admission; root collector owns QM comparison."""
    bound_inputs(cfg,manifest)
    release_gate(cfg,read(cfg['preparation_result']),manifest,runtime=False)
    require(result['implementation']==cfg['implementation'] and result['config']==cfg and result['panel']==cfg['panel'] and read(result['manifest'])==manifest,'Runtime provenance differs')
    require(result['candidate_result']==cfg['candidate_result'] and result['selected_review']==cfg['selected_review'],'Runtime selected authority differs')
    require(result['status']=='matched_pair_model_complete' and result['error'] is None and result['global_stop'] is False and result['native_program_invocations']==34 and set(result['rows'])==set(IDS),'Runtime denominator/failure differs')
    require(result['saved_cluster_models']==manifest['saved_cluster_models'] and result['analytic_bare_ions']==manifest['analytic_bare_ions'],'Saved anchors differ')
    for cid,row in result['rows'].items():
        frozen=manifest['rows'][cid];require(row['status']=='complete' and row['parent_geometry_id']==frozen['parent_geometry_id'] and row['geometry']==frozen['geometry'] and row['input']==frozen['input'],'Native row identity differs')
        for k in ('xyz','key','parameters'):require(verify(row['runtime_input'][k]).read_bytes()==verify(frozen['input'][k]).read_bytes(),'Actual native input differs')
        parsed=parse_analyze(verify(row['native']['analyze_output']));require(parsed==row['native'] and row['energy_kcal_mol']==parsed['total_energy_kcal_mol'] and row['components']==parsed['components'],'Native energy/component replay differs')
        require(len(row['receipts'])==1,'Native receipt count differs');r=row['receipts'][0];directory=Path(row['runtime_input']['xyz']['path']).parent
        require(r['argv']==[cfg['analyze']['path'],'system.xyz','E'] and r['cwd']==str(directory) and r['returncode']==0 and r['output']==row['native']['analyze_output'] and all(r['threads'][k]=='1' for k in THREADS),'Native program receipt differs')
    cycle=model_cycle(result['rows'],manifest['saved_cluster_models'],manifest['analytic_bare_ions']);require(cycle==result['model_cycle'],'Model cycle replay differs')
    return {'rows':result['rows'],'saved_cluster_models':manifest['saved_cluster_models'],'analytic_bare_ion':manifest['analytic_bare_ions'],'model_cycle':cycle,'new_native_calls':0}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('config','prepare','run'));p.add_argument('--panel');p.add_argument('--config');p.add_argument('--output',required=True);a=p.parse_args()
    if a.stage=='config':save(a.output,default_config(a.panel));raise SystemExit(0)
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);cfg=json.loads(Path(a.config).read_text());start=time.monotonic();cpu=time.process_time()
    try:
        if a.stage=='prepare':
            require(cfg['root_release'] is None and cfg['preparation_result'] is None,'Preparation must be unreleased');prepare(cfg,out)
            r={'status':'prepared','implementation':pin(__file__),'config':pin(a.config),'manifest':pin(out/'MANIFEST.json'),'new_native_calls':0,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu};save(out/'RESULT.json',r)
        else:r=run(cfg,out)
    except Exception as exc:r={'status':'precall_or_preparation_failed','implementation':pin(__file__),'error':repr(exc)};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'output':str(out)}));raise SystemExit(0 if r['status'] in ('prepared','matched_pair_model_complete') else 1)
