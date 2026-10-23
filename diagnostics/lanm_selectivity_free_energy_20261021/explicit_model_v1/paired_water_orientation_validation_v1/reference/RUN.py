#!/usr/bin/env python3
"""Thin conditional four-water-cell wrapper around the unchanged native reference adapter."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
RESERVED=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/reserved_validation_v1'
FAMILIES=('aquo8',)
GIDS=tuple('validation_'+f+'_orientation_'+s+'_'+m for f in FAMILIES for s in ('minus','plus') for m in ('La','Dy'))
CELLS=tuple('ecp_pbe0__'+g for g in GIDS)
WATER_PARAMETERS=(3.906136716492464,.28782327521600637,3.5865909437359376,.2837564269122712,.20765190997435956,.3173248056679709)
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')

def require(ok,message):
    if not ok:raise ValueError(message)

def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):require(pin(p['path'])==p,'Changed evidence: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,value):
    with Path(path).open('x') as h:json.dump(value,h,indent=2,allow_nan=False);h.write('\n')
    return pin(path)

def module(p):
    s=importlib.util.spec_from_file_location('unchanged_reserved_reference',verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def default_config():
    panel=pin(RESERVED/'NATIVE_PANEL.json');p=read(panel)
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'root_declaration':pin(HERE.parent/'ROOT_DECLARATION.md'),
        'outline':pin(HERE.parents[1]/'paired_donor_validation_outline_v1/OUTLINE.md'),
        'reserved_panel':pin(RESERVED/'RESERVED_PANEL.json'),'native_panel':panel,'native_adapter':p['adapter'],
        'water_result':pin(ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_six_parameter_v1/run_v1/execution/RESULT.json'),
        'water_review':pin(ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_six_parameter_v1/result_review_v1/REVIEW.json'),
        'water_numerical_result':pin(ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_six_transfer_v1/run_v1/execution/RESULT.json'),
        'water_numerical_review':pin(ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_six_transfer_v1/actual_transfer_review_v1/REVIEW.json'),
        'slurm_wrapper':pin(ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/RUN.sbatch'),
        'water_cell_ids':list(CELLS),'ranks':4,'cpus_per_rank':1,'maximum_reference_calls':4,
        'candidate_result':None,'root_release':None,'reserved_outcomes_opened':False}

def load(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['water_cell_ids']==list(CELLS),'Source/inventory differs')
    require(cfg['ranks']==4 and cfg['cpus_per_rank']==1 and cfg['maximum_reference_calls']==4,'Resource/scope differs')
    require(cfg['candidate_result'] is None and cfg['root_release'] is None and cfg['reserved_outcomes_opened'] is False,'Source-only configuration was changed')
    for k in ('specification','root_declaration','outline','slurm_wrapper'):verify(cfg[k])
    require(cfg['native_panel']['sha256']=='75d3fc69703c8a5acca61f547e3fcf5d614cfe6ac653d2f3ecd3f34a93bcaf5a' and cfg['reserved_panel']['sha256']=='eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b','Reserved inventory changed')
    require(cfg['native_adapter']['sha256']=='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129','Native recipe changed')
    ref=module(cfg['native_adapter']);p=ref.load_panel(verify(cfg['native_panel']));reserved=read(cfg['reserved_panel'])
    require(len(p['cells'])==24 and len(p['geometries'])==24 and not p['free_fragments'] and p['adapter']==cfg['native_adapter'] and p['source_inventory']==cfg['reserved_panel'],'Original24 scope differs')
    chosen=[c for c in p['cells'] if c['cell_id'] in CELLS]
    require(len(chosen)==4 and [c['cell_id'] for c in chosen]==list(CELLS),'Water subset denominator/order differs')
    eligible={c['id']:c for c in reserved['checks']};geometries={g['geometry_id']:g for g in reserved['geometries']}
    for c in chosen:
        g=c['geometry'];require(g==geometries[c['geometry_id']] and eligible[c['geometry_id']]['passed'],'Reserved geometry eligibility changed')
        require(c['recipe']=='ecp_pbe0' and c['multiplicity']==1 and c['physical_multiplicity']==(1 if c['metal']=='La' else 6) and c['ecp_core']==(46 if c['metal']=='La' else 55),'Electronic-state convention differs')
        require(c['explicit_electrons']==sum(ref.Z[x] for x in g['symbols'])-c['charge']-c['ecp_core'] and c['explicit_electrons']%2==0,'Explicit electron count differs')
        symbols,xyz=ref.xyz_read(verify(g['xyz']));require(symbols==g['symbols'] and xyz==g['positions_A'],'Pinned XYZ differs')
        peer=next(q for q in chosen if q['geometry_id'][:-2]==c['geometry_id'][:-2] and q['metal']!=c['metal'])
        require(peer['geometry']['positions_A']==xyz and peer['charge']==c['charge'],'Matched La/Dy coordinates/charge differ')
    return ref,p,chosen

def prepare(config,out):
    cfg=json.loads(Path(config).read_text());out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    require(all(os.environ.get(k)=='1' for k in THREADS),'Preparation uses one local thread')
    ref,p,cells=load(cfg);admission=water_admission(cfg);templates={};inventory=[]
    for c in cells:
        d=out/'inputs'/c['cell_id'];d.mkdir(parents=True)
        deck=ref.deck(p,c,4,'@FLOOR_0.75_GRANTED_MIB_DIV_4@')
        require(deck.startswith('! PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad\n') and '%pal nprocs 4 end' in deck and ' Guess PModel\n MaxIter 300\n' in deck and 'N_core '+str(c['ecp_core']) in deck,'Actual frozen deck differs')
        (d/'endpoint.inp.template').write_text(deck);templates[c['cell_id']]=pin(d/'endpoint.inp.template')
        inventory.append({k:c[k] for k in ('cell_id','geometry_id','recipe','metal','charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')}|{'atoms':c['geometry']['atoms'],'xyz':c['geometry']['xyz']})
    manifest={'status':'prepared_unreleased','configuration':pin(config),'implementation':pin(__file__),'native_adapter':cfg['native_adapter'],
        'native_panel':cfg['native_panel'],'reserved_panel':cfg['reserved_panel'],'original_denominator':24,'water_denominator':4,'full_Stage_A_denominator':16,'other_original_slots_unexecuted':20,
        'inventory':inventory,'input_templates':templates,'frozen_water_admission':admission,'orca':p['orca'],'assets':p['assets'],'resource_rule':'four MPI ranks; MaxCore=floor(0.75*actual granted MiB/4); default scheduler memory/time',
        'candidate_result':None,'consumed_water_numerical_evidence':cfg['water_numerical_result'],'root_release':None,'native_calls':0,'reserved_outcomes_opened':False}
    mp=save(out/'MANIFEST.json',manifest)
    commands={c['cell_id']:['sbatch','--qos=standard_expand','--ntasks=4','--cpus-per-task=1',cfg['slurm_wrapper']['path'],str(Path(__file__).resolve()),'run','--config',str(Path(config).resolve()),'--manifest',mp['path'],'--release','@ROOT_RELEASE.json@','--cell',c['cell_id'],'--out','@NEW_TASK_DIRECTORY@/execution'] for c in cells}
    save(out/'LAUNCH_TEMPLATES.json',{'status':'unreleased_templates_only','commands':commands,'native_calls':0})
    result={'status':'prepared_unreleased','manifest':mp,'configuration':pin(config),'source':pin(__file__),'cells':4,'native_calls':0,'reserved_outcomes_opened':False,'candidate_bound':False,
        'local_threads':1,'elapsed_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu}
    save(out/'RESULT.json',result);print(json.dumps(result))

def water_admission(cfg):
    require(cfg['water_result']['sha256']=='6027e319a261eb85171fb06e30aa3e519411f6371a7905eeca0e8a79812c48fc' and cfg['water_review']['sha256']=='a3147adcd7448012437e3ac275872fa58fe65e98a114619a95ec01e17a0df187','Frozen water authorities differ')
    require(cfg['water_numerical_result']['sha256']=='e3eeb870aac22f408c8b38a5f050404cfee5a8f4c0468633c07122fbc5e07254' and cfg['water_numerical_review']['sha256']=='3f7b3a9b375cf3cdad0580a75852ce9488a664232aa1da8b1254b547ac2b4ce4','Consumed numerical authorities differ')
    water=read(cfg['water_result']);review=read(cfg['water_review']);actual=read(cfg['water_numerical_result']);nr=read(cfg['water_numerical_review'])
    require(water['status']=='bounded_six_parameter_water_complete' and water['job_id']=='17605' and water['error'] is None and water['candidate_frozen'] and water['development_checks_passed'],'Passing actual water candidate required')
    require(water['selected']['parameters']==list(WATER_PARAMETERS) and water['selected']['score']['development_checks_passed'],'Frozen water values/controls differ')
    require(review['status']=='saved_result_review_passed' and review['result']==cfg['water_result'] and review['all_selected_raw_replays_exact'] and review['exact_evaluator_replay'] and review['selected_parameters']==list(WATER_PARAMETERS),'Independent water replay unavailable')
    require(actual['status']=='consumed_six_parameter_transfer_complete' and actual['job_id']=='17612' and actual['error'] is None and actual['stage_A_passed'] and not actual['global_stop'],'Consumed transfer native evidence unavailable')
    require(nr['status']=='passed_actual_raw_numerical_and_comparison_review' and nr['actual_result']==cfg['water_numerical_result'],'Consumed numerical review unavailable')
    bound=read(actual['bound']);require(bound['candidate_result']==cfg['water_result'] and bound['parameters']==list(WATER_PARAMETERS),'Consumed input/candidate join differs')
    numerical={}
    for metal in ('La','Dy'):
        gid='aquo8_central_'+metal;row=actual['rows'][gid];selected=water['selected']['model'][gid];c=row['cartesian_check']
        require(row['status']=='complete' and all(row['model'][k]==selected[k] for k in ('energy','forces','native')),'Consumed water model changed')
        require(c['fd_passed'] and c['all_atoms_active_and_components_present'] and c['fd_total_separation_A']==1e-5 and c['fd_each_side_displacement_A']==5e-6,'Original Cartesian recipe/admission differs')
        require(abs(c['energy']-selected['energy'])<=1e-4 and c['max_component_error_kcal_mol_A']<=1e-3 and max(abs(a-b) for aa,bb in zip(c['forces'],selected['forces']) for a,b in zip(aa,bb))<=1e-3,'Consumed numerical tolerances failed')
        verify(c['output']);verify(selected['provenance']['native_authority'])
        files=row['input']['files'];require(files==bound['inputs'][gid]['files'],'Consumed water inputs differ')
        for key in ('xyz','key','parameters'):
            verify(files[key]);require(files[key]['sha256']==selected['provenance']['effective_input'][key]['sha256'],'Selected water physical input differs')
        numerical[gid]={'cartesian_output':c['output'],'native_authority':selected['provenance']['native_authority'],'input':files,'passed':True}
    return {'water_result':cfg['water_result'],'water_review':cfg['water_review'],'water_numerical_result':cfg['water_numerical_result'],'water_numerical_review':cfg['water_numerical_review'],
        'parameters':list(WATER_PARAMETERS),'consumed_water_checks':numerical,'historical_donor_failures_preserved':True,'reserved_water_accuracy_not_yet_known':True}

def gate(config,manifest,release):
    cfg=json.loads(Path(config).read_text());ref,p,cells=load(cfg);mp=pin(manifest);prepared=read(mp);rp=pin(release);r=read(rp)
    require(prepared['configuration']==pin(config) and prepared['implementation']==pin(__file__) and prepared['status']=='prepared_unreleased' and prepared['water_denominator']==4,'Preparation binding differs')
    require(r['status']=='approved_independent_water_orientation_reference' and r['implementation']==pin(__file__) and r['configuration']==pin(config) and r['prepared_manifest']==mp and r['root_declaration']==cfg['root_declaration'],'Root reference release differs')
    require(r['cell_ids']==list(CELLS) and r['maximum_reference_calls']==4 and r['ranks']==4 and r['automatic_retry'] is False,'Released finite scope differs')
    for key in ('water_result','water_review','water_numerical_result','water_numerical_review'):require(r[key]==cfg[key],'Root water authority differs: '+key)
    evidence=water_admission(cfg)
    return cfg,ref,p,cells,{'release':rp,'manifest':mp,'water_admission':evidence}

def run(config,manifest,release,cell,out):
    cfg,ref,p,cells,evidence=gate(config,manifest,release);require(cell in CELLS,'Cell outside the four reserved aquo orientations')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','0'))==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Exact four-rank allocation required')
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);require(not out.exists(),'No automatic replay/retry of existing task directory')
    save(Path(str(out)+'.release.json'),evidence|{'cell_id':cell,'job_id':os.environ['SLURM_JOB_ID'],'native_calls_authorized':1})
    return ref.run(cfg['native_panel']['path'],cell,out)

def collect(config,manifest,release,result_inventory,out):
    cfg,ref,p,cells,evidence=gate(config,manifest,release);paths=json.loads(Path(result_inventory).read_text());require(set(paths)==set(CELLS),'Full4 water result-path denominator required')
    rows={};cost=0.;available=0
    for c in cells:
        gid=c['geometry_id'];row={'status':'unavailable','geometry':c['geometry'],'reference':None};path=Path(paths[c['cell_id']]);rows[gid]=row
        if not path.exists():row['reason']='missing_result';continue
        authority=pin(path);saved=read(authority);cost+=saved.get('allocated_cpu_seconds',0.)
        try:
            require(saved['cell_id']==c['cell_id'] and saved['panel']==cfg['native_panel'] and saved['adapter']==cfg['native_adapter'],'Result provenance differs')
            require(saved['status']=='complete' and saved['returncode']==0 and all(saved['checks'].values()),'Reference native admission failed')
            again=ref.collect_one(p,c,path.parent,saved['returncode']);require(all(saved[k]==v for k,v in again.items()),'Raw reference replay differs')
            row.update(status='complete',result=authority,reference={'energy':saved['energy_hartree']*ref.HA,'energy_hartree':saved['energy_hartree'],'forces':saved['forces_kcal_mol_A'],
                'native_checks':saved['checks'],'printed_scf_criteria':saved['printed_scf_criteria'],'state':{k:c[k] for k in ('charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')},
                'output':saved['output'],'engrad':saved['engrad'],'input':saved['input'],'geometry':saved['geometry'],'physical_state_qualified':False,'method_accuracy_qualified':False});available+=1
        except Exception as exc:row.update(reason=repr(exc),result=authority)
    return save(out,{'status':'complete' if available==4 else 'unavailable','declared':4,'available':available,'rows':rows,'binding':evidence,'result_inventory':pin(result_inventory),
        'allocated_native_cpu_seconds':cost,'native_calls_this_collection':0,'model_scoring_performed':False,'physical_accuracy_qualified':False})

if __name__=='__main__':
    a=argparse.ArgumentParser();sub=a.add_subparsers(dest='action',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--out',required=True)
    for action in ('run','collect'):
        q=sub.add_parser(action)
        for name in ('config','manifest','release','out'):q.add_argument('--'+name,required=True)
        q.add_argument('--cell' if action=='run' else '--result-inventory',required=True)
    args=a.parse_args()
    if args.action=='prepare':prepare(args.config,args.out)
    elif args.action=='run':raise SystemExit(run(args.config,args.manifest,args.release,args.cell,args.out))
    else:print(json.dumps(collect(args.config,args.manifest,args.release,args.result_inventory,args.out)))
