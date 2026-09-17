#!/usr/bin/env python3
"""Frozen-input, one-model-per-protein XoxF preparation; never executes ORCA."""
from __future__ import annotations
import argparse
import collections
import copy
import hashlib
import importlib.util
import json
import math
import multiprocessing
import os
from pathlib import Path
import shutil
import sys
import traceback

for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','OPENMM_CPU_THREADS'):
    os.environ[_key]='1'

A=Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
HERE=Path(__file__).resolve().parent
CAL=A/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914'
PROTOCOL='pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'
WRAPPER=A/'diagnostics/plm_adh9_fixed_core_20260916/prepare_candidates.py'
WRAPPER_SHA='8d988aa7bd5e10924f879b464b58f4ad3403829e0b21b7db6e55bfefb304ec5a'
REUSE_IDS={'PQQSEQ_242fa05e3ffc20087d42','PQQSEQ_faa97386262eec4316fc'}
ROLE_MAP={'conserved_metal_Glu':'anchor_glutamate','conserved_metal_Asn':'anchor_asparagine',
          'catalytic_Asp':'catalytic_aspartate','Ln_associated_position':'extra_acidic_ligand_homolog'}
GATES={'minimum_confidence':.9,'minimum_CN':7,'cutoff_A':3.1,'maximum_direct_N':2,'maximum_Asp_partner_A':3.5}
SELECTOR={'protein_chain':'A','metal':{'chain':'B','resnum':1,'resname':'LA','atom':'LA','icode':''},
          'pqq':{'chain':'C','resnum':1,'resname':'PQQ','icode':''}}
read=lambda p:json.loads(Path(p).read_text())

def record(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(rec):
    if record(rec['path'])['sha256']!=rec['sha256']:raise ValueError('Pinned input changed: '+rec['path'])
    return Path(rec['path'])

def write(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m

def wrapper():
    if record(WRAPPER)['sha256']!=WRAPPER_SHA:raise ValueError('Frozen wrapper changed')
    return load('xoxf_batch_frozen_preparer',WRAPPER)

def select(rows):
    """Rank all structurally valid samples, without filtering admission first."""
    if not rows:raise ValueError('No structurally valid AF3 samples')
    return sorted(rows,key=lambda r:(-r['CN'],-r['protein_La_iptm'],-r['protein_CN'],r['sample']))[0]

def admission(row):
    return row['CN']>=7 and row['protein_La_iptm']>=.9 and row['direct_N']<=2

def by_id(rows):
    result={r['target_id']:r for r in rows}
    if len(result)!=len(rows):raise ValueError('Duplicate target identifier')
    return result

def raw_model(w,target,entry):
    """Get donor metrics without applying CN/confidence eligibility prematurely."""
    import gemmi
    cif=verify(entry['source_cif']);summary=verify(entry['summary_confidences'])
    structure=gemmi.read_structure(str(cif));w.heavy_map(structure);model=structure[0]
    chains=[c.name for c in model]
    if len(chains)!=3 or set(chains)!={'A','B','C'}:raise ValueError('Expected exactly protein A, La B and PQQ C')
    residues=list(model['A']);observed=''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in residues)
    if observed!=target['sequence'] or [r.seqid.num for r in residues]!=list(range(1,len(observed)+1)):
        raise ValueError('Fold sequence or residue numbering differs from frozen full protein')
    if any(a.element.name in {'H','D'} for c in model for r in c for a in r):raise ValueError('Raw AF3 model contains hydrogen')
    if any(not all(math.isfinite(v) for v in (a.pos.x,a.pos.y,a.pos.z)) for c in model for r in c for a in r):
        raise ValueError('Nonfinite raw source coordinates')
    sites=w.base.select_metal_sites(model,site_chain='B',site_resnum=1,site_icode='',site_resname='LA',site_atom='LA')
    if len(sites)!=1 or sites[0].element.upper()!='LA':raise ValueError('Expected exactly selected La')
    site=sites[0];pqq=w.base._prepare_nearby_pqq(model,site.atom.pos,w.fixed.DEFAULT_PQQ_MICROSTATE)
    w.base._reject_unsupported_nearby_species(model,site,pqq)
    pkey=w.fixed.residue_key(SELECTOR['pqq'],'PQQ')
    if set(pqq)!={pkey}:raise ValueError('Missing, multiple or incomplete selected PQQ')
    _,prepared=pqq[pkey]
    if prepared.microstate.formal_charge!=-3 or prepared.formula!='C14H3N2O8' or prepared.schema_id!='pdb_ccd_pqq_v1':
        raise ValueError('PQQ chemistry/schema differs from frozen method')
    contacts,_,sulfur,untyped=w.base._classify_contacts(model,site,pqq)
    s=read(summary);confidence=float(s['chain_pair_iptm'][chains.index('A')][chains.index('B')])
    if not math.isfinite(confidence):raise ValueError('Nonfinite protein–La confidence')
    return {'sample':entry['sample'],'seed':entry['seed'],'CN':len(contacts),
       'protein_CN':sum(c.residue.chain=='A' for c in contacts),'direct_N':sum(c.element.upper()=='N' for c in contacts),
       'protein_La_iptm':confidence,'protein_PQQ_iptm':float(s['chain_pair_iptm'][chains.index('A')][chains.index('C')]),
       'source_cif':record(cif),'summary_confidence':record(summary),
       'typed_contacts':[w.fixed.contact_metadata(c,frozenset()) for c in contacts]}

def map_roles(w,target,chosen):
    """Use frozen four-position homology; unique geometric Arg/Lys partner."""
    import gemmi
    model=gemmi.read_structure(chosen['source_cif']['path'])[0];chain=model['A'];by_num={r.seqid.num:r for r in chain}
    observations=target['role_mappings']
    if isinstance(observations,dict):observations=[dict(value,role=key) for key,value in observations.items()]
    indexed={r['role']:r for r in observations}
    if len(observations)!=4 or set(indexed)!=set(ROLE_MAP):raise ValueError('Four unambiguous homologous core observations required')
    roles={}
    for old,new in ROLE_MAP.items():
        row=indexed[old]
        if row['status']!='observed' or row.get('sequence_identity_verified') is not True:
            raise ValueError('Unobserved or unverified mapped role: '+old)
        pos=int(row['full_protein_position']);aa=row['amino_acid'];r=by_num[pos]
        if target['sequence'][pos-1]!=aa or gemmi.find_tabulated_residue(r.name).one_letter_code!=aa:
            raise ValueError('Mapped residue differs from frozen full sequence: '+old)
        if r.name not in w.fixed.ROLE_ALLOWED_RESNAMES[new]:raise ValueError('Unsupported mapped '+new+': '+r.name)
        roles[new]={'chain':'A','resname':r.name,'resnum':pos,'icode':''}
    catpos=roles['catalytic_aspartate']['resnum'];extra=roles['extra_acidic_ligand_homolog']['resnum']
    if extra!=catpos+2 or target['sequence'][catpos-2:catpos]!='WD':raise ValueError('Catalytic WD / D+2 homology inconsistent')
    cat=by_num[catpos];candidates=[]
    for r in chain:
        if r.name not in {'ARG','LYS'}:continue
        distances=sorted((a.pos.dist(b.pos),a.name,b.name) for a in cat if a.name in {'OD1','OD2'}
                         for b in r if b.name in w.fixed.CATION_HBOND_ATOMS[r.name])
        if distances:candidates.append({'resname':r.name,'resnum':r.seqid.num,'distance_A':distances[0][0],
                                      'asp_atom':distances[0][1],'partner_atom':distances[0][2]})
    candidates.sort(key=lambda x:x['distance_A']);near=[p for p in candidates if p['distance_A']<=3.5]
    if len(near)!=1:raise ValueError(f'Expected one catalytic Asp Arg/Lys partner within3.5A; found {len(near)}')
    p=near[0];roles['catalytic_asp_cationic_partner']={'chain':'A','resname':p['resname'],'resnum':p['resnum'],'icode':''}
    w.role_state(model,{**SELECTOR,'roles':roles})
    return roles,candidates

def choose_target(w,target,fold,root,sources):
    uid=target['target_id'];directory=root/'selection'/uid;directory.mkdir(parents=True,exist_ok=False)
    review={'target_id':uid,'models':[],'selected_sample':None,'selection_rule':'CN descending; protein-La confidence descending; protein CN descending; sample index ascending; then selected-model admission'}
    try:
        if fold is None:raise ValueError('No fold outcome for frozen target')
        if fold['sequence_sha256']!=target['sequence_sha256']:raise ValueError('Fold/full-sequence SHA mismatch')
        native_input=fold.get('native_input')
        if native_input:verify(native_input)
        entries=fold['models']
        if len(entries)!=3 or {e['sample'] for e in entries}!={0,1,2}:raise ValueError('Expected exactly samples0,1,2')
        valid=[]
        for entry in entries:
            try:
                if entry['seed']!=101 or entry['status']!='complete':raise ValueError('AF3 model missing or seed differs')
                value=raw_model(w,target,entry);value['status']='valid';value['admission_pass']=admission(value);valid.append(value)
            except Exception as exc:value={'sample':entry.get('sample'),'status':'invalid_or_missing','reason':f'{type(exc).__name__}: {exc}','input':entry}
            review['models'].append(value)
        chosen=select(valid);review['selected_sample']=chosen['sample'];review['selected_metrics']=chosen
        if not admission(chosen):raise ValueError('Selected model fails frozen CN/confidence/direct-N gate; no substitution')
        roles,partners=map_roles(w,target,chosen);review.update(roles=roles,cationic_partners_by_distance=partners)
        case=f'{uid}_AF3_sample{chosen["sample"]}';destination=root/'adapter_inputs'/f'{case}.cif';summary=destination.with_name(case+'_summary.json')
        shutil.copyfile(chosen['source_cif']['path'],destination);shutil.copyfile(chosen['summary_confidence']['path'],summary)
        candidate={'case_id':case,'target_id':uid,'rank':chosen['sample'],'gene_id':target.get('gene_id',';'.join(target.get('source_gene_ids',[])) or uid),
            'source_cif':record(destination),'summary_confidence':record(summary),'roles':roles,**SELECTOR,
           'source_provenance':{'raw_AF3_model':chosen['source_cif'],'raw_AF3_summary':chosen['summary_confidence'],**sources}}
        if native_input:candidate['source_provenance']['native_AF3_input']=native_input
        review['status']='selected';write(directory/'selection.json',review)
        candidate['source_provenance']['selection_review']=record(directory/'selection.json')
        return candidate,None
    except Exception as exc:
        reason=f'{type(exc).__name__}: {exc}';review.update(status='unsupported',reason=reason)
        write(directory/'selection.json',review)
        return None,{'target_id':uid,'status':'unsupported','reason':reason,'selection_review':record(directory/'selection.json')}

def prepare_worker(payload):
    """One fresh process per protein; original RNG/minimizer state cannot leak."""
    target,pins_path,output,approval=payload
    from geometry_checks import check
    import hydrogen_preparation as hp
    hp.OUT=Path(output);hp.APPROVAL=Path(approval)
    w=wrapper();cid=target['case_id'];case={k:target[k] for k in ('case_id','target_id','rank','source_cif')}
    try:
        pins,original,cp=w.verify_pins(pins_path)
        raw,audit=hp.capture_original(w,target,pins,original,cp);carve=read(raw['path']);geometry=check(carve['qm_fragments'])
        original_geometry=copy.deepcopy(geometry);write(audit/'original_core_geometry.json',original_geometry);recovery=None
        if geometry['status']!='PASS':carve,geometry,recovery=hp.recover(w,target,raw,audit)
        else:
            for rec in carve['outputs'].values():rec['path']=str(Path(raw['path']).parent/rec['path'])
        observation=read(audit/'minimization_observation.json');force=observation['calls'][0]['final']['hydrogen_force_rms_kJ_mol_nm']
        details={'status':'PASS','target_id':target['target_id'],'selected_sample':target['rank'],'source_cif':target['source_cif'],
          'geometry':geometry,'original_geometry':original_geometry,'original_minimizer_observation':record(audit/'minimization_observation.json'),
          'original_minimizer_force_converged':force<=1.0,'original_H_force_rms_kJ_mol_nm':force,
          'hydrogen_recovery':recovery,'preparation_mode':'recovered_same_objective' if recovery else 'original_frozen_protonation',
          'all_heavy_atoms_preserved':True,'pqq_microstate':'pqq_ox_3minus_v1','calibration_implementation_changed':bool(recovery),
          'limitation':'Geometry PASS does not mean original50-step H minimization reached force tolerance. Numerical recovery, when used, is explicit; no calibration refit.'}
        validation=Path(output)/'prepared'/f'{cid}_geometry_validation.json';write(validation,details)
        carve.update(geometry_validation=record(validation),hydrogen_geometry_audit=details,metadata_packager=record(__file__),original_preparation_carve=raw)
        manifest=Path(output)/'prepared'/f'{cid}_carve_manifest.json';write(manifest,carve)
        case.update(status='ready_for_orca',reason=None,carve_manifest=record(manifest))
        outcome={'target_id':target['target_id'],'status':'ready','reason':None,'case_id':cid,'geometry_validation':record(validation)}
        result={'case':case,'outcome':outcome,'geometry':details}
    except Exception as exc:
        reason=f'{type(exc).__name__}: {exc}'
        result={'case':None,'outcome':{'target_id':target['target_id'],'status':'preparation_failed','reason':reason,'case_id':cid},
                'geometry':{'status':'UNSUPPORTED','target_id':target['target_id'],'reason':reason},'traceback':traceback.format_exc()}
    write(Path(output)/'worker_results'/f'{target["target_id"]}.json',result)
    return result

def worker_count(request,target_count):
    cpus=int(os.environ.get('SLURM_CPUS_ON_NODE','1').split('(')[0])
    # Conservative8GiB/worker, retaining20% of physical/allocation RAM; OpenMM workers each use one CPU.
    physical=os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES')
    mem_mib=os.environ.get('SLURM_MEM_PER_NODE');available=min(physical,int(mem_mib)*1024**2) if mem_mib else physical
    memory_limit=max(1,int(available*.8)//(8*1024**3))
    cap=max(1,min(cpus,memory_limit,max(1,target_count)))
    if request!='auto':cap=min(cap,max(1,int(request)))
    return cap,{'allocated_cpus':cpus,'available_memory_bytes':available,'memory_reserve_fraction':.2,'per_worker_budget_GiB':8,'worker_count':cap,'threads_per_worker':1,'fresh_process_per_target':True}

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('inventory','inventory-sha256','fold-input-manifest','fold-input-sha256','fold-results','approval','approval-sha256','output'):
        parser.add_argument('--'+name,required=True)
    parser.add_argument('--workers',default='auto');args=parser.parse_args(argv)
    root=Path(args.output).resolve()
    sources={'inventory':{'path':str(Path(args.inventory).resolve()),'sha256':args.inventory_sha256},
             'fold_input_manifest':{'path':str(Path(args.fold_input_manifest).resolve()),'sha256':args.fold_input_sha256},
             'approval':{'path':str(Path(args.approval).resolve()),'sha256':args.approval_sha256}}
    for rec in sources.values():verify(rec)
    inventory=read(args.inventory);fold_inputs=read(args.fold_input_manifest);approval=read(args.approval)
    if inventory.get('status')!='PASS':raise ValueError('Inventory source validation is not frozen PASS')
    if not os.environ.get('SLURM_JOB_ID'):raise ValueError('Production batch preparation requires an allocated Slurm CPU job')
    targets=by_id(inventory['targets']);fold_targets=by_id(fold_inputs['targets'])
    if set(targets)!=set(fold_targets):raise ValueError('Inventory/fold input target scope differs')
    if len(targets)!=176 or not REUSE_IDS.issubset(targets):raise ValueError('Frozen approved176 target/two reuse scope differs')
    for uid,target in targets.items():
        if hashlib.sha256(target['sequence'].encode()).hexdigest()!=target['sequence_sha256']:raise ValueError('Full sequence SHA mismatch: '+uid)
        if fold_targets[uid]['sequence_sha256']!=target['sequence_sha256']:raise ValueError('Fold input protein differs: '+uid)
        if (fold_targets[uid]['action']=='reuse_complete')!=(uid in REUSE_IDS):raise ValueError('Completed reuse action mismatch')
    sources['fold_results']=record(args.fold_results);fold_results=read(args.fold_results)
    if fold_results.get('schema')!='plm.xoxf_all.fold_results.v1' and fold_results.get('schema_version')!='plm.xoxf_all.fold_results.v1':raise ValueError('Unrecognized fold-results schema')
    # Runtime outputs must explicitly pin both immutable source manifests.
    def contains_record(value,expected):
        if isinstance(value,dict):return (value.get('sha256')==expected['sha256'] and value.get('path')==expected['path']) or any(contains_record(v,expected) for v in value.values())
        return isinstance(value,list) and any(contains_record(v,expected) for v in value)
    for key in ('inventory','fold_input_manifest'):
        if not contains_record(fold_results,sources[key]):raise ValueError('Fold results do not pin frozen '+key)
    completed=by_id(fold_results['targets'])
    if set(completed)-set(targets):raise ValueError('Unexpected folded target')
    for uid,fold in completed.items():
        if fold.get('native_input')!=fold_targets[uid].get('native_input'):
            raise ValueError('Runtime native input differs from frozen fold manifest: '+uid)
    reuse_source=approval['completed_reuse_results'];verify(reuse_source);reuse_data=read(reuse_source['path'])
    reused_rows=by_id(reuse_data['rows'])
    if set(reused_rows)!=REUSE_IDS or any(r['status']!='complete' for r in reused_rows.values()):
        raise ValueError('Completed reuse source does not contain the two approved complete pairs')
    root.mkdir(parents=True,exist_ok=False)
    for name in ('selection','adapter_inputs','original_preparation','prepared','worker_results'):(root/name).mkdir()
    w=wrapper();cp=CAL/'implementation_pins.json';w.fixed.verify_pins(cp)
    cases=[];outcomes=[];reused=[];candidates=[]
    for uid,target in sorted(targets.items()):
        if uid in REUSE_IDS:
            reused.append({'target_id':uid,'source':reuse_source,'status':'reused'})
            outcomes.append({'target_id':uid,'status':'reused','reason':'Existing completed, independently checked AF3 fixed-core result','source':reuse_source})
            continue
        candidate,failure=choose_target(w,target,completed.get(uid),root,sources)
        if candidate:candidates.append(candidate)
        else:outcomes.append(failure)
    candidate_path=root/'candidate_manifest.json';write(candidate_path,{'schema_version':'plm.adh9.candidate_manifest.v1','protocol_id':PROTOCOL,
         'approval':sources['approval'],'sources':sources,'targets':candidates,'target_count':len(targets),'reused_results':reused})
    pins_path=root/'candidate_implementation_pins.json'
    pins={'schema_version':'plm.adh9.fixed_core_pins.v1','protocol_id':PROTOCOL,'wrapper':record(WRAPPER),
       'candidate_manifest':record(candidate_path),'calibration_implementation_pins':record(cp),'calibration_result':record(CAL/'result.json'),
       'candidate_gates':GATES,'authorities':{**sources,**{p.name:record(p) for p in HERE.glob('*.py')}}}
    write(pins_path,pins);w.verify_pins(pins_path)
    nworkers,resources=worker_count(args.workers,len(candidates));write(root/'resources.json',resources)
    payloads=[(target,str(pins_path),str(root),str(Path(args.approval).resolve())) for target in candidates]
    results=[]
    if payloads:
        context=multiprocessing.get_context('spawn')
        with context.Pool(processes=nworkers,maxtasksperchild=1) as pool:
            for result in pool.imap_unordered(prepare_worker,payloads,chunksize=1):
                results.append(result);print(result['outcome']['target_id'],result['outcome']['status'],flush=True)
    for result in results:
        outcomes.append(result['outcome'])
        if result['case']:cases.append(result['case'])
    if set(by_id(outcomes))!=set(targets):raise ValueError('Target accounting incomplete')
    outcomes.sort(key=lambda x:x['target_id']);cases.sort(key=lambda x:x['target_id'])
    counts=dict(collections.Counter(x['status'] for x in outcomes));unsupported=[x for x in outcomes if x['status'] in {'unsupported','preparation_failed'}]
    validation=root/'geometry_validation.json';write(validation,{'status':'PASS','scope':'All admitted cases pass geometry; failed and unsupported targets remain unscored',
        'ready_count':len(cases),'unsupported_count':len(unsupported),'cases':{r['outcome']['case_id']:r['geometry'] for r in results if r['case']},
        'approval':sources['approval']})
    outcome_path=root/'target_outcomes.json';write(outcome_path,{'schema_version':'plm.xoxf_all.preparation_outcomes.v1','target_count':len(targets),
       'counts':counts,'per_target':outcomes,'sources':sources,'orca_executed':False})
    implementation=root/'implementation_pins.json';write(implementation,{'protocol_id':PROTOCOL,'approval':sources['approval'],
       'candidate_pins':record(pins_path),'candidate_manifest':record(candidate_path),'geometry_validation':record(validation),
       'target_outcomes':record(outcome_path),'scripts':{p.name:record(p) for p in HERE.glob('*.py')},
       'original_calibration_implementation_pins':record(cp),
       'H_observer_helpers':record(A/'diagnostics/plm_adh9_af3_20260916/hydrogen_repair/repair_hydrogens.py')})
    write(root/'prepared_pairs.json',{'schema_version':'plm.adh9.prepared_pairs.v1','protocol_id':PROTOCOL,
       'candidate_manifest':record(candidate_path),'implementation_pins':record(implementation),'geometry_validation':record(validation),
       'target_outcomes':record(outcome_path),'cases':cases,'prepared_pair_count':len(cases),'target_count':len(targets),
       'reused_results':reused,'reused_count':len(reused),'unsupported_targets':unsupported,'unsupported_count':len(unsupported),'orca_executed':False})
    print(json.dumps({'prepared_pairs':str(root/'prepared_pairs.json'),'counts':counts}),flush=True)

if __name__=='__main__':main()
