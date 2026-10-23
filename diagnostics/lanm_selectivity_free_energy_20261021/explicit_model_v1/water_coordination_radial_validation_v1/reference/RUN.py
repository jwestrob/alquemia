#!/usr/bin/env python3
"""Four conditional W9 radial references; unchanged original ECP native runner."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'
GIDS=tuple('aquo9_'+s+'_'+m for s in ('minus','plus') for m in ('La','Dy'))
CELLS=tuple('ecp_pbe0__'+g for g in GIDS)
TRAINING=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central') for m in ('La','Dy'))
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
FLOAT=r'[-+]?\d+(?:\.\d+)?(?:[eEdD][-+]?\d+)?'

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,value):
    with Path(path).open('x') as h:json.dump(value,h,indent=2,allow_nan=False);h.write('\n')
    return pin(path)
def module(p):
    s=importlib.util.spec_from_file_location('unchanged_reference',verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def default_config(geometry_manifest):
    ref=D/'reference_chemistry_v1'
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'root_declaration':pin(HERE.parent/'ROOT_DECLARATION.md'),
        'geometry_manifest':pin(geometry_manifest),'original_native_panel':pin(W/'reference_chemistry_v1/prepared_v2/PANEL.json'),
        'native_adapter':pin(ref/'REFERENCE.py'),'launcher':pin(ref/'RUN.sbatch'),
        'central_references':{m:pin(W/f'reference_chemistry_v1/ecp_remaining_v1/tasks/ecp_pbe0__aquo9_central_{m}/execution/RESULT.json') for m in ('La','Dy')},
        'candidate_implementation':pin(D/'explicit_model_v1/paired_water_coordination_joint_v1/REFINE.py'),
        'reference_cell_ids':list(CELLS),'maximum_native_calls':4,'MPI_ranks':4,'threads_per_rank':1,
        'candidate_result':None,'selected_review':None,'root_release':None,'fresh_outcomes_opened':False}

def stage_metadata(text):
    t=text.split('****END OF INPUT****')[-1];lines=t.splitlines()
    grid=t.split('COSX GRID GENERATION',1)[1].split('SCF SETTINGS',1)[0] if 'COSX GRID GENERATION' in t else ''
    signature=[l.strip() for l in grid.splitlines() if re.match(r'^(GRIDX\s+[123]|General Integration Accuracy|Radial Grid Type|Angular Grid \(max|Angular grid pruning|Weight generation|Basis function cutoff|Integration weight cutoff|Partially contracted|Rotationally invariant|Angular grids for H|UseSFitting)',l)]
    return {'COSX_grid_settings_signature':signature,'final_exchange_recomputed':'Recomputing exchange energy using gridx3' in t,
        'final_exchange_change_hartree':[float(v.replace('D','E').replace('d','e')) for v in re.findall(r'Exchange energy change after final integration\s*:\s*('+FLOAT+')',t)],
        'native_selected_lines':[{'line':i+1,'text':l} for i,l in enumerate(lines) if re.search(r'Last (Energy|MAX-Density|RMS-Density|DIIS|Orbital)|SCF CONVERGED|check signals convergence|TRAH|NR .*MAcro|final integration|exchange energy|Dispersion correction|FINAL SINGLE POINT ENERGY|WARNING',l,re.I)],
        'analytic_gradient_stage':'ORCA SCF GRADIENT CALCULATION' in t,'force_derivative_qualified':False}

def load(config):
    cfg=json.loads(Path(config).read_text())
    require(cfg==default_config(cfg['geometry_manifest']['path']),'Exact source/input configuration required')
    require(cfg['native_adapter']['sha256']=='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129','Original native adapter differs')
    require(cfg['original_native_panel']['sha256']=='d5929526cf70ac957c847a5936ed229e9edf1322764b0cc287a59292e01204f9','Original reference panel differs')
    ref=module(cfg['native_adapter']);original=ref.load_panel(verify(cfg['original_native_panel']));geometry=read(cfg['geometry_manifest'])
    require(set(geometry['rows'])==set(GIDS),'Exact four fresh geometry slots required')
    require(geometry['status']=='geometry_metadata_prepared' and geometry['cell_order']==list(GIDS) and geometry['count']==4 and all(c['passed'] for c in geometry['checks']) and geometry['new_native_calls']==0 and geometry['fit_trajectory_or_outcomes_read'] is False and geometry['validation_energies_or_forces_read'] is False,'Frozen zero-call geometry admission differs')
    require(geometry['root_declaration']==cfg['root_declaration'],'Geometry declaration differs')
    verify(geometry['implementation']);verify(geometry['configuration'])
    original_typed={r['id']:r for r in read(geometry['source_manifest'])['rows']}
    centers={};cells=[]
    for metal in ('La','Dy'):
        c=next(c for c in original['cells'] if c['cell_id']=='ecp_pbe0__aquo9_central_'+metal)
        r=read(cfg['central_references'][metal]);replay=ref.collect_one(original,c,Path(cfg['central_references'][metal]['path']).parent,r['returncode'])
        require(r['status']=='complete' and all(r['checks'].values()) and all(r[k]==v for k,v in replay.items()),'Central raw reference replay failed')
        require(verify(r['input']).read_text()==ref.deck(original,c,4,3000),'Actual central method/resource deck differs')
        centers[metal]=(c,r)
    for gid in GIDS:
        row=geometry['rows'][gid];g=row['geometry'];metal=gid.rsplit('_',1)[1];center=centers[metal][0];old=center['geometry']
        require(g['geometry_id']==gid and g['metal']==metal and g['atoms']==28 and g['water_count']==9 and g['charge']==3 and g['physical_multiplicity']==(1 if metal=='La' else 6),'Fresh geometry state/inventory differs')
        symbols,xyz=ref.xyz_read(verify(g['xyz']));require(symbols==g['symbols']==old['symbols'] and xyz==g['positions_A'],'Plain coordinate bytes differ')
        typed=original_typed['aquo9_central_'+metal]
        require(xyz[0]==old['positions_A'][0] and row['input']['bonds_zero_based']==typed['bonds_zero_based'] and row['input']['types']==typed['types'],'Metal/topology/types changed')
        verify(row['input']['xyz']);verify(row['source']['central_ordinary_xyz']);verify(row['source']['central_typed_xyz'])
        require(row['source']['central_ordinary_xyz']['sha256']==old['xyz']['sha256'] and row['source']['central_typed_xyz']==typed['xyz'],'Exact central geometry lineage differs')
        shift=-.1 if '_minus_' in gid else .1
        for i in range(1,28,3):
            delta=[old['positions_A'][i][k]-old['positions_A'][0][k] for k in range(3)];distance=math.sqrt(sum(v*v for v in delta));require(distance>0,'Zero central metal-O distance')
            unit=[v/distance for v in delta]
            for j in (i,i+1,i+2):require(max(abs(xyz[j][k]-old['positions_A'][j][k]-shift*unit[k]) for k in range(3))<=2e-12,'Rigid radial translation differs')
        c=copy.deepcopy(center);c.update(cell_id='ecp_pbe0__'+gid,geometry_id=gid,geometry=g,pilot=False)
        require(c['recipe']=='ecp_pbe0' and c['charge']==3 and c['multiplicity']==1 and c['physical_multiplicity']==(1 if metal=='La' else 6) and c['ecp_core']==(46 if metal=='La' else 55) and c['explicit_electrons']==98,'Reference electronic conventions changed')
        require(ref.deck(original,c,4,3000)==ref.deck(original,center,4,3000),'Fresh versus central native deck changed')
        cells.append(c)
    for s in ('minus','plus'):require(geometry['rows'][f'aquo9_{s}_La']['geometry']['positions_A']==geometry['rows'][f'aquo9_{s}_Dy']['geometry']['positions_A'],'Paired geometry mismatch')
    return cfg,ref,original,geometry,cells,centers

def prepare(config,out):
    start,cpu=time.monotonic(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    require(all(os.environ.get(k)=='1' for k in THREADS),'Preparation must use one local thread')
    cfg,ref,original,geometry,cells,centers=load(config)
    panel=copy.deepcopy(original);panel.update(geometries=[c['geometry'] for c in cells],cells=cells,free_fragments=[],pilot_cell_ids=[],source_inventory=cfg['geometry_manifest'],source_panel=cfg['original_native_panel'])
    panelpin=save(out/'NATIVE_PANEL.json',panel);templates={};inventory=[]
    for c in cells:
        folder=out/'inputs'/c['cell_id'];folder.mkdir(parents=True)
        (folder/'endpoint.inp.template').write_text(ref.deck(panel,c,4,'@FLOOR_0.75_GRANTED_MIB_DIV_4@'))
        shutil.copyfile(verify(c['geometry']['xyz']),folder/'geometry.xyz')
        templates[c['cell_id']]={'input':pin(folder/'endpoint.inp.template'),'geometry':pin(folder/'geometry.xyz')}
        inventory.append({k:c[k] for k in ('cell_id','geometry_id','recipe','metal','charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')}|{'atoms':28})
    anchors={m:{'raw_result':cfg['central_references'][m],'energy':r['energy_hartree']*ref.HA,'forces':r['forces_kcal_mol_A'],'energy_stages':stage_metadata(verify(r['output']).read_text()),'native_checks':r['checks'],'input':r['input'],'output':r['output'],'engrad':r['engrad']} for m,(c,r) in centers.items()}
    manifest={'status':'prepared_unbound_no_calls','implementation':pin(__file__),'configuration':pin(config),'geometry_manifest':cfg['geometry_manifest'],'native_adapter':cfg['native_adapter'],'native_panel':panelpin,'declared':4,'inventory':inventory,'input_templates':templates,'anchor_reference':anchors,'orca':panel['orca'],'assets':panel['assets'],'candidate_result':None,'selected_review':None,'root_release':None,'native_calls':0,'fresh_outcomes_opened':False}
    mp=save(out/'MANIFEST.json',manifest)
    release={'status':None,'implementation':pin(__file__),'configuration':pin(config),'prepared_manifest':mp,'root_declaration':cfg['root_declaration'],'cell_ids':list(CELLS),'maximum_reference_calls':4,'ranks':4,'automatic_retry':False,'candidate_result':None,'selected_review':None,'parameters':None,'task_output_directories':{c:'@NEW_TASK_DIRECTORY@/execution' for c in CELLS},'independent_source_review':None}
    save(out/'ROOT_RELEASE.template.json',release)
    commands={c:['sbatch','--qos=standard_expand','--ntasks=4','--cpus-per-task=1',cfg['launcher']['path'],str(Path(__file__).resolve()),'run','--config',str(Path(config).resolve()),'--manifest',mp['path'],'--release','@ROOT_RELEASE.json@','--cell',c,'--out','@NEW_TASK_DIRECTORY@/execution'] for c in CELLS}
    save(out/'LAUNCH_TEMPLATES.json',{'status':'unbound_templates_only','commands':commands,'calls':0})
    rp=save(out/'RESULT.json',{'status':'prepared_unbound_no_calls','manifest':mp,'configuration':pin(config),'cells':4,'central_references_reused':2,'native_calls':0,'candidate_bound':False,'fresh_outcomes_opened':False,'local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu})
    print(json.dumps({'manifest':mp,'preparation':rp}));return rp

def gate(config,manifest,release):
    cfg,ref,original,geometry,cells,centers=load(config);mp=pin(manifest);prepared=read(mp);rp=pin(release);r=read(rp)
    require(prepared['status']=='prepared_unbound_no_calls' and prepared['configuration']==pin(config) and prepared['implementation']==pin(__file__) and prepared['geometry_manifest']==cfg['geometry_manifest'] and prepared['declared']==4,'Prepared source/inventory differs')
    require(r['status']=='approved_fresh_W9_radial_reference' and r['implementation']==pin(__file__) and r['configuration']==pin(config) and r['prepared_manifest']==mp and r['root_declaration']==cfg['root_declaration'],'Exact root release required')
    require(r['cell_ids']==list(CELLS) and r['maximum_reference_calls']==4 and r['ranks']==4 and r['automatic_retry'] is False and set(r['task_output_directories'])==set(CELLS) and len(set(r['task_output_directories'].values()))==4,'Released finite denominator differs')
    verify(r['independent_source_review']);candidate=read(r['candidate_result']);review=read(r['selected_review']);verify(review['implementation'])
    require(candidate['job_id']=='17676' and candidate['status']=='bounded_joint_coordination_water_complete' and candidate['error'] is None and candidate['candidate_frozen'] and candidate['development_checks_passed'] and candidate['implementation']==cfg['candidate_implementation'],'Terminal passing17676 required')
    selected=candidate['selected'];require(selected['score']['development_checks_passed'] and set(selected['model'])==set(TRAINING) and selected['parameters']==r['parameters'] and len(r['parameters'])==6,'Frozen passing selected-ten candidate differs')
    require(review['status']=='selected_ten_raw_review_passed' and review['result']==r['candidate_result'] and review['selected_parameters']==r['parameters'] and review['selected_cell_ids']==list(TRAINING) and review['complete']==10 and review['new_native_calls']==0 and review['error'] is None,'Independent selected-ten review differs')
    require(all(review[k] is True for k in ('all_selected_raw_replays_exact','exact_selected_score_replayed','all_development_controls_passed')),'Independent selected-ten review not passed')
    p=ref.load_panel(verify(prepared['native_panel']));require(p['cells']==cells and p['source_inventory']==cfg['geometry_manifest'] and p['source_panel']==cfg['original_native_panel'],'Prepared native cells differ')
    return cfg,ref,p,geometry,cells,prepared,{'release':rp,'manifest':mp,'candidate_result':r['candidate_result'],'selected_review':r['selected_review'],'parameters':r['parameters'],'task_output_directories':r['task_output_directories']}

def run(config,manifest,release,cell,out):
    cfg,ref,p,geometry,cells,prepared,evidence=gate(config,manifest,release);require(cell in CELLS,'Cell outside four fresh references')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','0'))==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Exact four MPI ranks required')
    out=Path(out).resolve();require(str(out)==evidence['task_output_directories'][cell] and not out.exists(),'Exact released fresh task path required; no retry')
    out.parent.mkdir(parents=True,exist_ok=True);save(Path(str(out)+'.release.json'),evidence|{'cell_id':cell,'job_id':os.environ['SLURM_JOB_ID'],'native_calls_authorized':1})
    return ref.run(prepared['native_panel']['path'],cell,out)

def collect(config,manifest,release,result_inventory,out):
    cfg,ref,p,geometry,cells,prepared,evidence=gate(config,manifest,release);paths=json.loads(Path(result_inventory).read_text());require(set(paths)==set(CELLS),'Exact four-result denominator required')
    rows={};available=0;cost=0.
    for c in cells:
        gid=c['geometry_id'];meta=geometry['rows'][gid];row={'status':'unavailable','geometry':c['geometry'],'input':meta['input'],'bodies':meta['bodies'],'reference':None};rows[gid]=row
        expected=Path(evidence['task_output_directories'][c['cell_id']])/'RESULT.json';require(Path(paths[c['cell_id']]).resolve()==expected,'Result inventory differs from released task path')
        if not expected.exists():row['reason']='missing_result';continue
        authority=pin(expected)
        try:
            saved=read(authority);cost+=saved.get('allocated_cpu_seconds',0.)
            require(saved['status']=='complete' and saved['returncode']==0 and all(saved['checks'].values()) and saved['cell_id']==c['cell_id'] and saved['panel']==prepared['native_panel'] and saved['adapter']==cfg['native_adapter'],'Native reference admission/provenance failed')
            replay=ref.collect_one(p,c,expected.parent,saved['returncode']);require(all(saved[k]==v for k,v in replay.items()),'Raw reference replay differs')
            require(len(saved['forces_kcal_mol_A'])==28 and all(math.isfinite(v) for f in saved['forces_kcal_mol_A'] for v in f),'Finite full28 force inventory required')
            native_start=json.loads((expected.parent/'START.json').read_text());require(native_start['ranks']==4 and native_start['input']==saved['input'] and verify(saved['input']).read_text()==ref.deck(p,c,4,native_start['maxcore_MiB']),'Executed native recipe/resource differs')
            stages=stage_metadata(verify(saved['output']).read_text());center=prepared['anchor_reference'][c['metal']]['energy_stages']
            changed=stages['COSX_grid_settings_signature']!=center['COSX_grid_settings_signature'] or stages['final_exchange_recomputed']!=center['final_exchange_recomputed'] or stages['analytic_gradient_stage']!=center['analytic_gradient_stage']
            row.update(status='complete',reference={'energy':saved['energy_hartree']*ref.HA,'energy_hartree':saved['energy_hartree'],'forces':saved['forces_kcal_mol_A'],'raw_result':authority,'native_checks':saved['checks'],'state':{k:c[k] for k in ('charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')},'energy_stages':stages,'changed_central_energy_stage_signature':changed,'output':saved['output'],'input':saved['input'],'engrad':saved['engrad'],'physical_state_qualified':False,'force_derivative_qualified':False,'method_accuracy_qualified':False});available+=1
        except Exception as exc:row.update(reason=repr(exc),raw_result=authority)
    return save(out,{'status':'complete' if available==4 else 'unavailable','declared':4,'available':available,'rows':rows,'anchor_reference':prepared['anchor_reference'],'binding':evidence,'result_inventory':pin(result_inventory),'allocated_native_CPU_seconds':cost,'native_calls_this_collection':0,'model_scoring_performed':False,'physical_accuracy_qualified':False})

if __name__=='__main__':
    a=argparse.ArgumentParser();sub=a.add_subparsers(dest='action',required=True)
    q=sub.add_parser('config');q.add_argument('--geometry-manifest',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--out',required=True)
    for action in ('run','collect'):
        q=sub.add_parser(action)
        for name in ('config','manifest','release','out'):q.add_argument('--'+name,required=True)
        q.add_argument('--cell' if action=='run' else '--result-inventory',required=True)
    args=a.parse_args()
    if args.action=='config':save(args.out,default_config(args.geometry_manifest))
    elif args.action=='prepare':prepare(args.config,args.out)
    elif args.action=='run':raise SystemExit(run(args.config,args.manifest,args.release,args.cell,args.out))
    else:print(json.dumps(collect(args.config,args.manifest,args.release,args.result_inventory,args.out)))
