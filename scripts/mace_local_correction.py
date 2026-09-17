"""Matched local/global MACE decomposition, with explicit archived DFT availability."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
import numpy as np
from affordable_common import InvalidArtifact, cache_key, energy, HA_TO_KCAL, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, check_atoms, write_xyz
from mace_global_prepare import CASES, LABELS
from mace_global_benchmark import snapshot, numerical_parent_gate, CONTRASTS
from mace_gb import MODEL as GB_MODEL

MACE_SCHEMA='alquemia.mace_local_correction.v1'
GB_SCHEMA='alquemia.mace_local_gb.v1'
POLICY='archived_and_global_H_source_mapped_core_v1'
STATES=('archived','global_H')
TOL={'total_charge_e':1e-5,'physical_mapping_A':.002}


def source_mapping(prep,core):
    """Map real source atoms, explicitly excluding each documented artificial cap."""
    old={a['id']:a for a in prep['preparation_details']['original_protein_atoms']}
    full={a['id']:a for a in prep['physical_atoms']}
    if len(old)!=len(prep['preparation_details']['original_protein_atoms']):raise InvalidArtifact('duplicate original source identity')
    manifest=read_json(verify(prep['source_preparation']))
    mapping=[{'core_index':0,'physical_id':'metal','kind':'metal'}]
    if 'atom_graph' in manifest:
        for m in manifest['atom_graph']['source_to_qm']:
            i=m['qm_index']
            if m['kind']=='source':
                s=m['source'];identity=f"{s['chain']}/{s['resnum']}/{s['insertion_code']}/{s['atom']}"
                mapping.append({'core_index':i,'physical_id':identity,'kind':'physical_source'})
            elif m['kind']=='sigma_link_H':
                if core[i][0]!='H':raise InvalidArtifact('cap is not H')
                mapping.append({'core_index':i,'physical_id':None,'kind':'synthetic_cap','cap_record':m})
            else:raise InvalidArtifact('unsupported graph atom kind')
    else:
        offset=1
        for f in manifest['qm_fragments']:
            count=f['atom_count']
            if f['kind']=='fixed_core_pqq':
                for i,a in enumerate(f['atom_records']):
                    mapping.append({'core_index':offset+i,'physical_id':'pqq/'+a['name'],'kind':'PQQ'})
            elif f['kind'] in ('fixed_core_protein_sidechain','fixed_core_cationic_sidechain'):
                match=re.fullmatch(r'([^:]+):([A-Z]+)(-?\d+)([A-Za-z]?)',f['id'])
                if not match:raise InvalidArtifact('unsupported PQQ fragment residue identity')
                chain,resname,resid,icode=match.groups();prefix=f'{chain}/{resid}/{icode}/'
                for i,a in enumerate(core[offset:offset+count]):
                    if i==count-1:
                        if a[0]!='H':raise InvalidArtifact('documented fragment cap is absent')
                        mapping.append({'core_index':offset+i,'physical_id':None,'kind':'synthetic_cap','fragment':f['id']})
                    else:
                        found=[s for identity,s in old.items() if identity.startswith(prefix) and s['resname']==resname and s['element']==a[0]
                               and np.linalg.norm(np.array(a[1:])-s['xyz_A'])<=TOL['physical_mapping_A']]
                        if len(found)!=1:raise InvalidArtifact('ambiguous/missing exact-fragment physical source match')
                        mapping.append({'core_index':offset+i,'physical_id':found[0]['id'],'kind':'physical_source'})
            else:raise InvalidArtifact('unsupported frozen fragment kind')
            offset+=count
        if offset!=len(core):raise InvalidArtifact('unaccounted frozen fragment atoms')
    mapping=sorted(mapping,key=lambda a:a['core_index'])
    if [a['core_index'] for a in mapping]!=list(range(len(core))):raise InvalidArtifact('core atom map is not complete and unique')
    physical=[a['physical_id'] for a in mapping if a['physical_id'] is not None]
    if len(set(physical))!=len(physical):raise InvalidArtifact('physical source atom mapped twice')
    changed=[];rows=list(core)
    water_before={m['id']:m['before_A'] for m in prep['water_H_moves']}
    for m in mapping:
        i=m['core_index'];identity=m['physical_id'];atom=core[i]
        if identity is None:continue
        if identity not in full:raise InvalidArtifact('core source missing from declared global system')
        original=old.get(identity,full[identity]);original_position=water_before.get(identity,original['xyz_A'])
        if identity!='metal' and original['element']!=atom[0]:raise InvalidArtifact('mapped source element mismatch')
        error=float(np.linalg.norm(np.array(atom[1:])-original_position))
        if error>TOL['physical_mapping_A']:raise InvalidArtifact('core/original source-coordinate mismatch')
        m['original_coordinate_error_A']=error
        if atom[0]=='H' and m['kind']!='PQQ':
            new=tuple(full[identity]['xyz_A']);rows[i]=(atom[0],*new)
            changed.append({'core_index':i,'physical_id':identity,'before_A':list(atom[1:]),'after_A':list(new)})
        elif np.linalg.norm(np.array(atom[1:])-full[identity]['xyz_A'])>TOL['physical_mapping_A']:
            raise InvalidArtifact('fixed heavy/PQQ coordinates incompatible with full system')
    return rows,mapping,changed


def archived_dft(endpoint,manifests):
    candidates=[]
    direct=Path(verify(endpoint['input'])).with_suffix('.out')
    if direct.is_file():candidates.append(direct)
    for path in manifests:
        m=read_json(path)
        for t in m['tasks']:
            if t['source_xyz']['sha256']==endpoint['xyz']['sha256'] and t['source_input']['sha256']==endpoint['input']['sha256']:
                candidates.append(Path(t['output_path']))
    candidates=sorted(set(candidates))
    if not candidates:return {'status':'unavailable','energy_hartree':None,'reason':'no exact archived endpoint found'}
    if len(candidates)!=1:raise InvalidArtifact('ambiguous archived DFT endpoint')
    path=candidates[0];receipt=read_json(str(path)+'.execution.json')
    if receipt['returncode']!=0 or not receipt['normal_termination'] or not receipt['scf_converged']:
        raise InvalidArtifact('archived DFT execution failed')
    artifacts=receipt['artifacts']
    for key in ('output','runtime_input','runtime_sidecar','template_input','xyz'):verify(artifacts[key])
    if artifacts['output']!=record(path) or artifacts['xyz']['sha256']!=endpoint['xyz']['sha256'] or artifacts['template_input']['sha256']!=endpoint['input']['sha256']:
        raise InvalidArtifact('DFT receipt is incompatible with exact core input')
    return {'status':'computed','energy_hartree':energy(path),'output':record(path),
            'receipt':record(str(path)+'.execution.json'),'input':endpoint['input'],'xyz':endpoint['xyz'],
            'orca_version':receipt['orca_version'],'units':'hartree'}


def specs():
    return [(case,state,metal) for case in CASES for state in STATES for metal in ('La','Ca')]


def prepare_mace(global_collection,dft_manifests,agreement,output):
    global_c=read_json(global_collection);gm=read_json(verify(global_c['manifest']))
    if gm['checkpoint_label']!='medium' or global_c['status']!='complete' or gm['schema_version']!='alquemia.mace_global_benchmark.v1':
        raise InvalidArtifact('completed primary medium global MACE panel required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,(*gm['implementation'],'mace_local_correction.py'))
    cases={};tasks=[]
    for case in CASES:
        prep=read_json(verify(gm['cases'][case]));endpoints=prep['source_audit_row']['endpoints']
        old=xyz(verify(endpoints['La']['xyz']));new,mapping,changes=source_mapping(prep,old)
        d=out/case;d.mkdir()
        baseline={metal:archived_dft(endpoints[metal],dft_manifests) for metal in ('La','Ca')}
        baseline_r=None
        if all(r['status']=='computed' for r in baseline.values()):baseline_r=(baseline['Ca']['energy_hartree']-baseline['La']['energy_hartree'])*HA_TO_KCAL
        data={'case_id':case,'policy_id':POLICY,'global_preparation':gm['cases'][case],'source_endpoints':endpoints,
              'mapping':mapping,'H_changes':changes,'archived_DFT':baseline,'archived_DFT_R_kcal_mol':baseline_r,
              'global_H_DFT':None,'global_H_DFT_status':'not_computed','evidence':LABELS[case]}
        write_new(d/'mapping.json',data);cases[case]=record(d/'mapping.json')
        for state in STATES:
            for metal in ('La','Ca'):
                endpoint=endpoints[metal]
                if state=='archived':path=verify(endpoint['xyz']);rows=xyz(path)
                else:
                    rows=[(metal if i==0 else a[0],*a[1:]) for i,a in enumerate(new)]
                    path=d/(metal+'_'+state+'.xyz');write_xyz(path,rows)
                tasks.append({'task_id':f'{case}_{state}_{metal}','case_id':case,'local_state':state,'kind':'core',
                              'variant':'primary','metal':metal,'charge':endpoint['charge'],'spin_multiplicity':1,
                              'xyz':record(path),'state':check_atoms(rows,endpoint['charge']),'mapping':cases[case],
                              'rotation_matrix':np.eye(3).tolist()})
    model=copy.deepcopy(gm['model']);model['preparation_policy']=POLICY
    m={'schema_version':MACE_SCHEMA,'protocol_id':'mace_polar_1m_local_core_H_comparison_v1',
       'agreement':record(agreement),'software':gm['software'],'model':model,'implementation':pins,
       'global_mace_collection':record(global_collection),'numerical_reference':gm['numerical_reference'],
       'DFT_source_manifests':[record(p) for p in dft_manifests],'cases':cases,'tasks':tasks,'tolerances':TOL,
       'reference':None,'S_kcal_mol':None,'calibrated_class':None,'hybrid':None,
       'run_inventory':{'new_MACE_calls':20,'new_DFT_endpoints':0},'evidence_use':'consumed_method_development'}
    for t in tasks:t['cache_key']=cache_key({'task':t,'model':model,'software':m['software'],'implementation':pins})
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest);is_gb=m['schema_version']==GB_SCHEMA
    if m['schema_version'] not in (MACE_SCHEMA,GB_SCHEMA) or m['tolerances']!=TOL:raise InvalidArtifact('unsupported local protocol/tolerances')
    sm=read_json(verify(m['software']))
    for ref in [m['agreement'],*m['implementation'].values(),sm['python'],sm['requirements'],
                *read_json(verify(sm['backend_source_inventory']))['files'],*m['DFT_source_manifests']]:verify(ref)
    parent=read_json(verify(m['global_mace_collection']));pm=read_json(verify(parent['manifest']))
    expected=copy.deepcopy(pm['model']);expected['preparation_policy']=POLICY
    if numerical_parent_gate(m)['status']!='pass':raise InvalidArtifact('matching MACE numerical gate required')
    if is_gb:
        ref=read_json(verify(m['solver_validation']));rm=read_json(verify(ref['manifest']))
        if not ref['numerical_checks_pass'] or m['model']!=GB_MODEL or m['software']!=rm['software']:
            raise InvalidArtifact('matching GB physical model/software gate required')
        source=read_json(verify(m['source_mace_collection']));lm=read_json(verify(source['manifest']))
        if source['status']!='complete' or lm['cases']!=m['cases']:raise InvalidArtifact('incomplete/mismatched local MACE source')
        validate(verify(source['manifest']))
    elif m['model']!=expected or m['software']!=pm['software']:
        raise InvalidArtifact('local MACE changes the electronic model')
    if len(m['tasks'])!=20 or set(m['cases'])!=set(CASES):raise InvalidArtifact('local inventory changed')
    percase={}
    for case,ref in m['cases'].items():
        c=read_json(verify(ref));prep=read_json(verify(c['global_preparation']))
        if c['global_preparation']!=pm['cases'][case] or c['evidence']!=LABELS[case]:raise InvalidArtifact('case geometry/evidence differs')
        endpoints=prep['source_audit_row']['endpoints'];core=xyz(verify(endpoints['La']['xyz']))
        new,mapping,moves=source_mapping(prep,core)
        if (len(mapping)!=len(c['mapping']) or moves!=c['H_changes'] or c['source_endpoints']!=endpoints):
            raise InvalidArtifact('local source mapping/H preparation changed')
        for stored,replayed in zip(c['mapping'],mapping):
            numeric='original_coordinate_error_A'
            if ({k:v for k,v in stored.items() if k!=numeric}!={k:v for k,v in replayed.items() if k!=numeric}
                    or abs(stored.get(numeric,0.)-replayed.get(numeric,0.))>1e-12):
                raise InvalidArtifact('local atom identity/coordinate mapping changed')
        for metal,r in c['archived_DFT'].items():
            if r['status']=='computed':
                verify(r['receipt']);value=energy(verify(r['output']))
                if value!=r['energy_hartree'] or r['xyz']!=endpoints[metal]['xyz']:raise InvalidArtifact('archived DFT provenance changed')
        percase[case]=(core,new,endpoints)
    for task,(case,state,metal) in zip(m['tasks'],specs()):
        if task['task_id']!=f'{case}_{state}_{metal}' or (task['case_id'],task['local_state'],task['metal'])!=(case,state,metal):
            raise InvalidArtifact('task inventory identity mismatch')
        old,new,endpoints=percase[case]
        expected_rows=[(metal if i==0 else a[0],*a[1:]) for i,a in enumerate(old if state=='archived' else new)]
        rows=xyz(verify(task['xyz']))
        if ([a[0] for a in rows]!=[a[0] for a in expected_rows] or
                not np.allclose([a[1:] for a in rows],[a[1:] for a in expected_rows],atol=1e-12,rtol=0)
                or task['charge']!=endpoints[metal]['charge'] or task['spin_multiplicity']!=1):
            raise InvalidArtifact('local paired geometry/charge mismatch')
        if state=='archived' and task['xyz']!=endpoints[metal]['xyz']:raise InvalidArtifact('archived XYZ is not byte-exact')
        check_atoms(rows,task['charge'])
        if is_gb:
            r=source['rows'][task['task_id']]
            from mace_hybrid import accepted_attempt
            source_manifest=verify(source['manifest'])
            source_task=next(t for t in lm['tasks'] if t['task_id']==task['task_id'])
            candidates=(source_manifest.parent/'execution'/task['task_id']).glob('attempt_*')
            if r not in [accepted_attempt(a,source_task,source_manifest) for a in candidates]:
                raise InvalidArtifact('local GB source has no matching successful MACE receipt')
            if task['source_density']!=r['density_coefficients'] or task['source_vacuum_energy_eV']!=r['energy_eV']:
                raise InvalidArtifact('local GB density/energy differs from MACE source')
            verify(task['source_density'])
            if task['solver']!='native' or task['platform']!='CUDA' or task['solvent_dielectric']!=78.5:raise InvalidArtifact('local GB physical settings differ')
        base={k:v for k,v in task.items() if k!='cache_key'}
        if task['cache_key']!=cache_key({'task':base,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('local scientific cache key mismatch')
    return {'status':'pass','manifest':record(manifest),'tasks':20,'new_DFT_endpoints':0}


def prepare_gb(collection,global_gb,solver_validation,agreement,output):
    c=read_json(collection);m=read_json(verify(c['manifest']));validate(verify(c['manifest']))
    if c['status']!='complete':raise InvalidArtifact('local MACE endpoints incomplete')
    g=read_json(global_gb);gm=read_json(verify(g['manifest']))
    if g['status']!='complete' or gm['source_mace_collection']!=m['global_mace_collection']:
        raise InvalidArtifact('global solvent output lacks compatible physical/model source')
    ref=read_json(solver_validation);rm=read_json(verify(ref['manifest']))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);pins=snapshot(out,(*m['implementation'],'mace_local_correction.py'))
    result=copy.deepcopy(m);tasks=[]
    for t in m['tasks']:
        t=copy.deepcopy(t);t.pop('cache_key');r=c['rows'][t['task_id']]
        t.update(solver='native',platform='CUDA',solvent_dielectric=78.5,source_density=r['density_coefficients'],source_vacuum_energy_eV=r['energy_eV'])
        tasks.append(t)
    result.update(schema_version=GB_SCHEMA,protocol_id='mace_polar_1m_local_frozen_obc2_H_comparison_v1',
                  agreement=record(agreement),model=GB_MODEL,software=rm['software'],implementation=pins,tasks=tasks,
                  source_mace_collection=record(collection),global_gb_collection=record(global_gb),solver_validation=record(solver_validation),
                  run_inventory={'new_GB_calls':20,'new_MACE_calls':0,'new_DFT_endpoints':0})
    for t in tasks:t['cache_key']=cache_key({'task':t,'model':result['model'],'software':result['software'],'implementation':pins})
    write_new(out/'manifest.json',result);return validate(out/'manifest.json')


def collect_local(manifest):
    from mace_hybrid import accepted_attempt
    mp=Path(manifest).resolve();m=read_json(mp);is_gb=m['schema_version']==GB_SCHEMA
    rows={};attempts=[]
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:
                if is_gb:verify(r['serialized_system'])
                valid.append(r)
            attempts.append({'task_id':t['task_id'],'attempt':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']]=valid[-1] if valid else {'status':'unavailable','energy_eV':None}
    source=read_json(verify(m['source_mace_collection'])) if is_gb else None
    global_c=read_json(verify(m['global_gb_collection'])) if is_gb else read_json(verify(m['global_mace_collection']))
    scores={}
    for case in CASES:
        prep=read_json(verify(m['cases'][case]));r={'states':{},'archived_DFT_R_kcal_mol':prep['archived_DFT_R_kcal_mol'],
            'global_H_DFT':None,'hybrid':None,'evidence':LABELS[case],'geometry_effect_kcal_mol':None,'environment_contribution_kcal_mol':None}
        for state in STATES:
            la,ca=[rows[f'{case}_{state}_{metal}'] for metal in ('La','Ca')]
            score={'status':'unavailable','R_kcal_mol':None,'S_kcal_mol':None,'calibrated_class':None}
            if la['status']==ca['status']=='computed':
                value=(ca['energy_eV']-la['energy_eV'])*EV_TO_KCAL
                score.update(status='computed',R_kcal_mol=value)
                if is_gb:
                    gb=ca['GB_reaction_kcal_mol']-la['GB_reaction_kcal_mol']
                    raw=source['scores'][case]['states'][state]['R_kcal_mol']
                    score.update(R_kcal_mol=raw+gb,MACE_vacuum_R_kcal_mol=raw,GB_Ca_minus_La_kcal_mol=gb)
            r['states'][state]=score
        a,b=[r['states'][s]['R_kcal_mol'] for s in STATES]
        if a is not None and b is not None:r['geometry_effect_kcal_mol']=b-a
        full=global_c['scores'][case]['R_kcal_mol']
        if b is not None and full is not None:r['environment_contribution_kcal_mol']=full-b
        if is_gb and r['environment_contribution_kcal_mol'] is not None:
            raw_core=source['scores'][case]['states']['global_H']['R_kcal_mol']
            full_score=global_c['scores'][case]
            r['environment_components']={'vacuum_kcal_mol':full_score['MACE_vacuum_R_kcal_mol']-raw_core,
                'GB_kcal_mol':full_score['GB_Ca_minus_La_kcal_mol']-r['states']['global_H']['GB_Ca_minus_La_kcal_mol']}
        if is_gb and a is not None and prep['archived_DFT_R_kcal_mol'] is not None:
            r['archived_DFT_minus_local_MACE_GB_kcal_mol']=prep['archived_DFT_R_kcal_mol']-a
        scores[case]=r
    complete=all(r['status']=='computed' for r in rows.values())
    contrasts=[]
    for higher,lower in CONTRASTS:
        r={'higher_expected':higher,'lower_expected':lower,'local_states':{}}
        for state in STATES:
            a,b=[scores[case]['states'][state]['R_kcal_mol'] for case in (higher,lower)]
            r['local_states'][state]=None if a is None or b is None else a-b
        a,b=[scores[case]['environment_contribution_kcal_mol'] for case in (higher,lower)]
        r['environment_contribution_difference_kcal_mol']=None if a is None or b is None else a-b
        contrasts.append(r)
    return {'status':'complete' if complete else 'incomplete','protocol_id':m['protocol_id'],'manifest':record(mp),
            'collection_implementation':record(__file__),'rows':rows,'attempts':attempts,'scores':scores,'contrasts':contrasts,
            'reference':None,'S_kcal_mol':None,'calibrated_class':None,'hybrid':None,'new_DFT_endpoints':0}


def report(collection,output):
    c=read_json(collection);m=read_json(verify(c['manifest']))
    if m['schema_version']!=GB_SCHEMA:raise InvalidArtifact('report requires the local solvent collection')
    lines=['# Matched local/global MACE decomposition','',f"Status: {c['status']}. Baseline unchanged.",'',
           'All values are kcal/mol, R = E(Ca) - E(La). Global minus core uses matched corrected H coordinates.',
           'Archived DFT uses its original core geometry; no matched corrected-H DFT or valid hybrid is available.','',
           '| Case | Archived DFT R | Archived MACE+GB R | Corrected-H MACE+GB R | Full−core R contribution | H effect |',
           '|---|---:|---:|---:|---:|---:|']
    for case,r in c['scores'].items():
        values=[r['archived_DFT_R_kcal_mol'],r['states']['archived']['R_kcal_mol'],r['states']['global_H']['R_kcal_mol'],
                r['environment_contribution_kcal_mol'],r['geometry_effect_kcal_mol']]
        lines.append('| '+case+' | '+' | '.join('unavailable' if v is None else f'{v:.6f}' for v in values)+' |')
    lines.extend(['','## Predeclared differences between sites','','```json',json.dumps(c['contrasts'],indent=2),'```','',
                  'Positive local differences support the predeclared qualitative ordering. A negative environment',
                  'difference erodes that ordering; it includes changed cavity, direct interactions and learned',
                  'charge response together. Components do not establish a unique cause.','',
                  'The archived DFT versus MACE+GB discrepancy includes both electronic and solvent approximations.',
                  'Do not mix archived DFT with corrected-H low-level subtraction and call it a computed hybrid.',
                  'All cases are consumed development data. Alpha structures are one observation; PQQ class',
                  'association is separate from direct affinity evidence. No threshold, aquo gauge or absolute class.'])
    with Path(output).open('x') as f:f.write('\n'.join(lines)+'\n')
    return {'status':'report_written','report':record(output)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare-mace')
    for key in ('global-collection','agreement','output'):q.add_argument('--'+key,required=True)
    q.add_argument('--dft-manifest',action='append',required=True)
    q=sub.add_parser('prepare-gb')
    for key in ('collection','global-gb','solver-validation','agreement','output'):q.add_argument('--'+key,required=True)
    q=sub.add_parser('report');q.add_argument('--collection',required=True);q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare-mace':r=prepare_mace(a.global_collection,a.dft_manifest,a.agreement,a.output)
    elif a.command=='prepare-gb':r=prepare_gb(a.collection,a.global_gb,a.solver_validation,a.agreement,a.output)
    else:r=report(a.collection,a.output)
    print(json.dumps(r,indent=2))
