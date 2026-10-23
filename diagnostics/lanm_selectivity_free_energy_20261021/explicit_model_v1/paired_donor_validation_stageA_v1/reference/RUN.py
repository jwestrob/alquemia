#!/usr/bin/env python3
"""Thin conditional16-cell wrapper around the unchanged native reference adapter."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
RESERVED=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/reserved_validation_v1'
FAMILIES=('aquo8','hydrated_acetate_mono','hydrated_acetate_bi','hydrated_acetamide')
GIDS=tuple('validation_'+f+'_orientation_'+s+'_'+m for f in FAMILIES for s in ('minus','plus') for m in ('La','Dy'))
CELLS=tuple('ecp_pbe0__'+g for g in GIDS)
CHECKS=tuple(prefix+f+'_central_'+m for prefix in ('hydrated_','matched_dry_') for f in ('acetate_mono','acetate_bi','acetamide') for m in ('La','Dy'))
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
        'slurm_wrapper':pin(ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/RUN.sbatch'),
        'stage_A_cell_ids':list(CELLS),'ranks':4,'cpus_per_rank':1,'maximum_reference_calls':16,
        'candidate_result':None,'root_release':None,'reserved_outcomes_opened':False}

def load(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['stage_A_cell_ids']==list(CELLS),'Source/inventory differs')
    require(cfg['ranks']==4 and cfg['cpus_per_rank']==1 and cfg['maximum_reference_calls']==16,'Resource/scope differs')
    require(cfg['candidate_result'] is None and cfg['root_release'] is None and cfg['reserved_outcomes_opened'] is False,'Source-only configuration was changed')
    for k in ('specification','root_declaration','outline','slurm_wrapper'):verify(cfg[k])
    require(cfg['native_panel']['sha256']=='75d3fc69703c8a5acca61f547e3fcf5d614cfe6ac653d2f3ecd3f34a93bcaf5a' and cfg['reserved_panel']['sha256']=='eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b','Reserved inventory changed')
    require(cfg['native_adapter']['sha256']=='eec0e60d982bda4e1a8685c6120f25c6ab74cc57e114cd3d8f370b4baccba129','Native recipe changed')
    ref=module(cfg['native_adapter']);p=ref.load_panel(verify(cfg['native_panel']));reserved=read(cfg['reserved_panel'])
    require(len(p['cells'])==24 and len(p['geometries'])==24 and not p['free_fragments'] and p['adapter']==cfg['native_adapter'] and p['source_inventory']==cfg['reserved_panel'],'Original24 scope differs')
    chosen=[c for c in p['cells'] if c['cell_id'] in CELLS]
    require(len(chosen)==16 and [c['cell_id'] for c in chosen]==list(CELLS),'Stage A denominator/order differs')
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
    ref,p,cells=load(cfg);templates={};inventory=[]
    for c in cells:
        d=out/'inputs'/c['cell_id'];d.mkdir(parents=True)
        deck=ref.deck(p,c,4,'@FLOOR_0.75_GRANTED_MIB_DIV_4@')
        require(deck.startswith('! PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad\n') and '%pal nprocs 4 end' in deck and ' Guess PModel\n MaxIter 300\n' in deck and 'N_core '+str(c['ecp_core']) in deck,'Actual frozen deck differs')
        (d/'endpoint.inp.template').write_text(deck);templates[c['cell_id']]=pin(d/'endpoint.inp.template')
        inventory.append({k:c[k] for k in ('cell_id','geometry_id','recipe','metal','charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')}|{'atoms':c['geometry']['atoms'],'xyz':c['geometry']['xyz']})
    manifest={'status':'prepared_unreleased','configuration':pin(config),'implementation':pin(__file__),'native_adapter':cfg['native_adapter'],
        'native_panel':cfg['native_panel'],'reserved_panel':cfg['reserved_panel'],'original_denominator':24,'stage_A_denominator':16,'stage_B_unexecuted':8,
        'inventory':inventory,'input_templates':templates,'orca':p['orca'],'assets':p['assets'],'resource_rule':'four MPI ranks; MaxCore=floor(0.75*actual granted MiB/4); default scheduler memory/time',
        'candidate_result':None,'consumed_numerical_review':None,'root_release':None,'native_calls':0,'reserved_outcomes_opened':False}
    mp=save(out/'MANIFEST.json',manifest)
    commands={c['cell_id']:['sbatch','--qos=standard_expand','--ntasks=4','--cpus-per-task=1',cfg['slurm_wrapper']['path'],str(Path(__file__).resolve()),'run','--config',str(Path(config).resolve()),'--manifest',mp['path'],'--release','@ROOT_RELEASE.json@','--cell',c['cell_id'],'--out','@NEW_TASK_DIRECTORY@/execution'] for c in cells}
    save(out/'LAUNCH_TEMPLATES.json',{'status':'unreleased_templates_only','commands':commands,'native_calls':0})
    result={'status':'prepared_unreleased','manifest':mp,'configuration':pin(config),'source':pin(__file__),'cells':16,'native_calls':0,'reserved_outcomes_opened':False,'candidate_bound':False,
        'local_threads':1,'elapsed_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu}
    save(out/'RESULT.json',result);print(json.dumps(result))

def gate(config,manifest,release):
    cfg=json.loads(Path(config).read_text());ref,p,cells=load(cfg);mp=pin(manifest);prepared=read(mp);rp=pin(release);r=read(rp)
    require(prepared['configuration']==pin(config) and prepared['implementation']==pin(__file__) and prepared['status']=='prepared_unreleased' and prepared['stage_A_denominator']==16,'Preparation binding differs')
    require(r['status']=='approved_paired_donor_validation_stageA_reference' and r['implementation']==pin(__file__) and r['configuration']==pin(config) and r['prepared_manifest']==mp and r['root_declaration']==cfg['root_declaration'],'Root reference release differs')
    require(r['cell_ids']==list(CELLS) and r['maximum_reference_calls']==16 and r['ranks']==4 and r['automatic_retry'] is False,'Released finite scope differs')
    candidate=read(r['candidate_result']);review=read(r['candidate_review']);checks=read(r['consumed_numerical_review'])
    require(candidate['status']=='bounded_paired_donor_complete' and candidate['job_id']=='17621' and candidate['error'] is None and candidate['candidate_frozen'] and candidate['development_checks_passed'] and candidate['selected']['score']['development_checks_passed'] and candidate['frozen_water_invariants']['passed'],'Actual passing donor candidate required')
    require(review['status']=='passed' and review['result']==r['candidate_result'] and review['development_checks_passed'],'Selected raw/input/objective review required')
    require(checks['status']=='consumed_numerical_checks_complete' and checks['numerical_checks_passed'] and checks['candidate_result']==r['candidate_result'] and checks['error'] is None and not checks['global_stop'],'Consumed numerical qualification required')
    require(set(checks['checks'])==set(CHECKS) and all(q['status']=='complete' and q['passed'] for q in checks['checks'].values()),'Every consumed numerical cell required')
    require(checks['native_program_invocations']<=12 and checks['internal_FD_energies_requested']<=1344 and checks['internal_FD_energies_completed']<=1344,'Consumed numerical scope changed')
    return cfg,ref,p,cells,{'release':rp,'manifest':mp,'candidate_result':r['candidate_result'],'candidate_review':r['candidate_review'],'consumed_numerical_review':r['consumed_numerical_review']}

def run(config,manifest,release,cell,out):
    cfg,ref,p,cells,evidence=gate(config,manifest,release);require(cell in CELLS,'Cell outside Stage A')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','0'))==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Exact four-rank allocation required')
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);require(not out.exists(),'No automatic replay/retry of existing task directory')
    save(Path(str(out)+'.release.json'),evidence|{'cell_id':cell,'job_id':os.environ['SLURM_JOB_ID'],'native_calls_authorized':1})
    return ref.run(cfg['native_panel']['path'],cell,out)

def collect(config,manifest,release,result_inventory,out):
    cfg,ref,p,cells,evidence=gate(config,manifest,release);paths=json.loads(Path(result_inventory).read_text());require(set(paths)==set(CELLS),'Full16 result-path denominator required')
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
    return save(out,{'status':'complete' if available==16 else 'unavailable','declared':16,'available':available,'rows':rows,'binding':evidence,'result_inventory':pin(result_inventory),
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
