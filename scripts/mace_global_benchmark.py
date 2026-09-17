"""Finite whole-protein MACE/GB development panel using the existing task runner."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import check_atoms, write_xyz, rotation, EV_TO_KCAL
from mace_global_prepare import CASES, LABELS, POLICY
from mace_gb import MODEL as GB_MODEL, TOL as GB_TOL

MACE_SCHEMA='alquemia.mace_global_benchmark.v1'
GB_SCHEMA='alquemia.mace_global_gb.v1'
TOL={'energy_kcal_mol':.01,'force_max_eV_A':.001,'total_charge_e':1e-5}
CONTRASTS=(('PQQ_4MAE','PQQ_1H4I'),('ALPHA_1F6S','GGR_1GLG'),('ALPHA_6IP9','GGR_1GLG'))


def snapshot(out,names):
    d=out/'implementation';d.mkdir()
    pins={}
    for name in set(names)|{'mace_global_benchmark.py','mace_global_prepare.py','mace_gb.py','mace_hydrogen.py','mace_response_trace.py','affordable_common.py','mace_hybrid.py'}:
        shutil.copyfile(Path(__file__).with_name(name),d/name);pins[name]=record(d/name)
    return pins


def physical_preparations(preparation):
    p=read_json(preparation)
    if p['policy_id']!=POLICY or tuple(r['case_id'] for r in p['rows'])!=CASES or p['status']!='complete':
        raise InvalidArtifact('five declared complete global preparations required')
    preservation=Path(preparation).with_name('implementation_preservation.json')
    archive=read_json(preservation)
    for name,ref in p['implementation'].items():
        if archive[name]['recorded_source']!=ref or archive[name]['preserved_copy']['sha256']!=ref['sha256']:
            raise InvalidArtifact('preparation implementation preservation mismatch')
        verify(archive[name]['preserved_copy'])
    for row in p['rows']:
        c=read_json(verify(row['preparation']))
        for key in ('source','source_preparation','forcefield','water_forcefield'):verify(c[key])
        if c['case_id']!=row['case_id'] or c['evidence']!=LABELS[row['case_id']]:
            raise InvalidArtifact('case identity/evidence mismatch')
        physical=c['physical_atoms']
        if len({a['id'] for a in physical})!=len(physical):raise InvalidArtifact('duplicate physical mapping')
        for metal in ('La','Ca'):
            endpoint=c['endpoints'][metal];rows=xyz(verify(endpoint['xyz']))
            expected=[(metal if a['element']=='M' else a['element'],*a['xyz_A']) for a in physical]
            if rows!=expected:raise InvalidArtifact('paired full coordinates do not match physical preparation')
            wanted=c['protein_charge_e']+c['cofactor_charge_e']+(3 if metal=='La' else 2)
            if endpoint['charge']!=wanted:raise InvalidArtifact('full formal charge closure failed')
            check_atoms(rows,wanted)
    return p,record(preservation)


def expected_model(parent):
    model=copy.deepcopy(parent['model'])
    for key in ('assembly','microstate','explicit_waters','response_trace'):model.pop(key,None)
    model['preparation_policy']=POLICY
    model['assembly_microstate_water_scope']='per_case_pinned_preparation'
    return model


def expected_tasks(checkpoint):
    return [(case,metal,variant) for case in CASES for variant in
            (('primary','rotate') if checkpoint=='medium' and case in ('GGR_1GLG','PQQ_4MAE') else ('primary',))
            for metal in ('La','Ca')]


def numerical_parent_gate(m):
    c=read_json(verify(m['numerical_reference']))
    good=(c['status']=='complete' and c['analytic_comparison']['core_gate'] and
          bool(c['checks']) and all(check['pass'] for check in c['checks']))
    return {'status':'pass' if good else 'fail',
            'scope':'analytic_core_and_rigid_checks; direct_scoring_does_not_use_subtractive_partition',
            'parent_partition_pass':c['partition_check']}


def prepare_mace(preparation,medium,large,agreement,output):
    p,preserved=physical_preparations(preparation)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    results={}
    for label,collection in (('medium',medium),('large',large)):
        c=read_json(collection);parent=read_json(verify(c['manifest']))
        d=out/label;d.mkdir();pins=snapshot(d,parent['implementation'])
        cases={r['case_id']:r['preparation'] for r in p['rows']};tasks=[]
        for case,metal,variant in expected_tasks(label):
            state=read_json(verify(cases[case]));endpoint=state['endpoints'][metal]
            rows=xyz(verify(endpoint['xyz']));matrix=rotation() if variant=='rotate' else np.eye(3)
            if variant=='rotate':
                center=np.array(next(a[1:] for a in rows if a[0]==metal))
                rows=[(a[0],*((np.array(a[1:])-center)@matrix.T+center)) for a in rows]
            task_id=f'{case}_{metal}_{variant}';path=d/(task_id+'.xyz');write_xyz(path,rows)
            tasks.append({'task_id':task_id,'case_id':case,'kind':'full','variant':variant,'metal':metal,
                          'charge':endpoint['charge'],'spin_multiplicity':1,'xyz':record(path),
                          'state':check_atoms(rows,endpoint['charge']),'rotation_matrix':matrix.tolist(),
                          'preparation':cases[case]})
        m={'schema_version':MACE_SCHEMA,'protocol_id':f'mace_polar_1{label[0]}_global_analytic_H_OXT_v1',
           'checkpoint_label':label,'agreement':record(agreement),'software':parent['software'],
           'preparation':record(preparation),'preserved_preparation_implementation':preserved,
           'cases':cases,'numerical_reference':record(collection),'implementation':pins,
           'model':expected_model(parent),'tasks':tasks,'tolerances':TOL,
           'reference':None,'S_kcal_mol':None,'calibrated_class':None,'evidence_use':'consumed_method_development',
           'run_inventory':{'new_MACE_energy_force_calls':len(tasks),'new_DFT_endpoints':0},'compute_budget':None}
        for t in tasks:t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':pins})
        write_new(d/'manifest.json',m);results[label]=validate(d/'manifest.json')
    return results


def validate(manifest):
    m=read_json(manifest);schema=m['schema_version'];is_gb=schema==GB_SCHEMA
    if schema not in (MACE_SCHEMA,GB_SCHEMA):raise InvalidArtifact('unsupported global schema')
    p,preserved=physical_preparations(verify(m['preparation']))
    if preserved!=m['preserved_preparation_implementation']:raise InvalidArtifact('preparation snapshot mismatch')
    software=read_json(verify(m['software']))
    for ref in [m['agreement'],*m['implementation'].values(),software['python'],software['requirements'],
                *read_json(verify(software['backend_source_inventory']))['files']]:verify(ref)
    if m['cases']!={r['case_id']:r['preparation'] for r in p['rows']} or m['tolerances']!=TOL:
        raise InvalidArtifact('case inventory/tolerances changed')
    label=m['checkpoint_label']
    if label not in ('medium','large'):raise InvalidArtifact('unsupported checkpoint label')
    expected=expected_tasks(label)
    if len(m['tasks'])!=len(expected):raise InvalidArtifact('task inventory changed')
    if is_gb:
        reference=read_json(verify(m['solver_validation']))
        if reference['status']!='complete' or not reference['numerical_checks_pass']:
            raise InvalidArtifact('passed GB solver checks required')
        validated=read_json(verify(reference['manifest']))
        if m['software']!=validated['software'] or m['model']!=GB_MODEL:
            raise InvalidArtifact('GB software or physical model changed')
        source=read_json(verify(m['source_mace_collection']));sm=read_json(verify(source['manifest']))
        if source['status']!='complete' or sm['cases']!=m['cases'] or sm['checkpoint_label']!=label:
            raise InvalidArtifact('complete matching MACE outputs required')
        validate(verify(source['manifest']))
    else:
        c=read_json(verify(m['numerical_reference']));parent=read_json(verify(c['manifest']))
        if numerical_parent_gate(m)['status']!='pass' or m['model']!=expected_model(parent) or m['software']!=parent['software']:
            raise InvalidArtifact('MACE implementation lacks matching numerical parent gate')
        verify(software['checkpoint'])
    for t,(case,metal,variant) in zip(m['tasks'],expected):
        if (t['case_id'],t['metal'],t['variant'])!=(case,metal,variant) or t['task_id']!=f'{case}_{metal}_{variant}':
            raise InvalidArtifact('task scientific identity changed')
        if t['preparation']!=m['cases'][case]:raise InvalidArtifact('task physical mapping changed')
        prepared=read_json(verify(t['preparation']));endpoint=prepared['endpoints'][metal]
        original=xyz(verify(endpoint['xyz']));rows=xyz(verify(t['xyz']))
        matrix=rotation() if variant=='rotate' else np.eye(3)
        coords=np.array([a[1:] for a in original]);center=np.array(next(a[1:] for a in original if a[0]==metal))
        wanted=(coords-center)@matrix.T+center if variant=='rotate' else coords
        if ([a[0] for a in rows]!=[a[0] for a in original] or
                not np.allclose([a[1:] for a in rows],wanted,atol=1e-12,rtol=0)
                or t['rotation_matrix']!=matrix.tolist() or t['charge']!=endpoint['charge'] or t['spin_multiplicity']!=1):
            raise InvalidArtifact('geometry or electronic state changed')
        check_atoms(rows,t['charge'])
        if is_gb:
            r=source['rows'][t['task_id']]
            st=next(x for x in sm['tasks'] if x['task_id']==t['task_id'])
            from mace_hybrid import accepted_attempt
            sourcepath=verify(source['manifest']).parent/'execution'/t['task_id']
            if r not in [accepted_attempt(a,st,verify(source['manifest'])) for a in sourcepath.glob('attempt_*')]:
                raise InvalidArtifact('GB source MACE receipt mismatch')
            if t['source_density']!=r['density_coefficients'] or t['source_vacuum_energy_eV']!=r['energy_eV']:
                raise InvalidArtifact('GB source electronic output changed')
            if t['platform']!='CUDA' or t['solver']!='native' or t['solvent_dielectric']!=78.5:
                raise InvalidArtifact('GB solver parameters changed')
        base={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':base,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('scientific cache key mismatch')
    return {'status':'pass','manifest':record(manifest),'tasks':len(m['tasks']),'new_DFT_endpoints':0}


def prepare_gb(collection,solver_validation,agreement,output):
    c=read_json(collection);m=read_json(verify(c['manifest']));validate(verify(c['manifest']))
    if c['status']!='complete':raise InvalidArtifact('MACE panel incomplete; do not fill missing charges')
    ref=read_json(solver_validation);rm=read_json(verify(ref['manifest']))
    if not ref['numerical_checks_pass']:raise InvalidArtifact('GB solver gate failed')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pins=snapshot(out,m['implementation'])
    tasks=[]
    for t in m['tasks']:
        task=copy.deepcopy(t);task.pop('cache_key');r=c['rows'][t['task_id']]
        task.update(platform='CUDA',solver='native',solvent_dielectric=78.5,
                    source_density=r['density_coefficients'],source_vacuum_energy_eV=r['energy_eV'])
        tasks.append(task)
    result=copy.deepcopy(m)
    result.update(schema_version=GB_SCHEMA,protocol_id=f'mace_polar_1{m["checkpoint_label"][0]}_global_frozen_obc2_H_OXT_v1',
                  agreement=record(agreement),software=rm['software'],implementation=pins,model=GB_MODEL,tasks=tasks,
                  source_mace_collection=record(collection),solver_validation=record(solver_validation),
                  run_inventory={'new_GB_energy_force_calls':len(tasks),'new_MACE_calls':0,'new_DFT_endpoints':0},
                  corrected_gradient=None,corrected_gradient_status='frozen_charge_GB_force_omits_charge_chain_rule')
    for task in tasks:task['cache_key']=cache_key({'task':task,'model':result['model'],'software':result['software'],'implementation':pins})
    write_new(out/'manifest.json',result);return validate(out/'manifest.json')


def collect_global(manifest):
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
    scores={};checks=[]
    parent=read_json(verify(m['source_mace_collection'])) if is_gb else None
    for case in CASES:
        score={'evidence':LABELS[case],'status':'unavailable','R_kcal_mol':None,'S_kcal_mol':None,'calibrated_class':None}
        la,ca=(rows[f'{case}_{metal}_primary'] for metal in ('La','Ca'))
        if la['status']==ca['status']=='computed':
            r=(ca['energy_eV']-la['energy_eV'])*EV_TO_KCAL
            score.update(status='computed',R_kcal_mol=r)
            if is_gb:
                # Convert each native kJ/mol endpoint once, then contrast.
                correction=ca['GB_reaction_kcal_mol']-la['GB_reaction_kcal_mol']
                raw=parent['scores'][case]['R_kcal_mol']
                score.update(R_kcal_mol=raw+correction,MACE_vacuum_R_kcal_mol=raw,GB_Ca_minus_La_kcal_mol=correction)
        scores[case]=score
        rotated=[]
        for metal in ('La','Ca'):
            name=f'{case}_{metal}_rotate'
            if name not in rows:continue
            primary=rows[f'{case}_{metal}_primary'];other=rows[name];rotated.append(other)
            if primary['status']!=other['status'] or primary['status']!='computed':
                checks.append({'name':name,'status':'unavailable','pass':False});continue
            delta=abs(primary['energy_eV']-other['energy_eV'])*EV_TO_KCAL
            f=np.load(verify(primary['forces']));g=np.load(verify(other['forces']))@rotation()
            df=float(np.max(np.abs(f-g)))
            checks.append({'name':name,'status':'computed','energy_difference_kcal_mol':delta,
                           'force_difference_max_eV_A':df,'pass':delta<=.01 and df<=.001})
        if rotated and all(r['status']=='computed' for r in [la,ca,*rotated]):
            delta=abs(rotated[1]['energy_eV']-rotated[0]['energy_eV']-(ca['energy_eV']-la['energy_eV']))*EV_TO_KCAL
            checks.append({'name':case+'_rotation_contrast','status':'computed','difference_kcal_mol':delta,'pass':delta<=.01})
    contrasts=[]
    for higher,lower in CONTRASTS:
        a,b=(scores[x]['R_kcal_mol'] for x in (higher,lower))
        delta=None if a is None or b is None else a-b
        contrasts.append({'higher_expected':higher,'lower_expected':lower,'difference_kcal_mol':delta,
                          'expected_order':None if delta is None else delta>0})
    complete=all(r['status']=='computed' for r in rows.values())
    numerical=complete and all(c['pass'] for c in checks)
    if parent is not None:numerical=numerical and parent['numerical_checks_pass']
    return {'status':'complete' if complete else 'incomplete','protocol_id':m['protocol_id'],'manifest':record(mp),
            'collection_implementation':record(__file__),'rows':rows,'attempts':attempts,'scores':scores,
            'checks':checks,'numerical_checks_pass':numerical,'contrasts':contrasts,
            'development_ordering_screen_pass':numerical and all(c['expected_order'] is True for c in contrasts),
            'reference':None,'S_kcal_mol':None,'calibrated_class':None,'corrected_gradient':None,
            'evidence_use':'consumed_method_development','independent_biological_groups':4,
            'structural_cases':5,'goal_completion_claimed':False}


def compare_report(medium_gb,large_gb,output):
    collections={label:read_json(path) for label,path in (('medium',medium_gb),('large',large_gb))}
    manifests={label:read_json(verify(c['manifest'])) for label,c in collections.items()}
    if any(m['schema_version']!=GB_SCHEMA or m['checkpoint_label']!=label for label,m in manifests.items()):
        raise InvalidArtifact('comparison requires the declared medium/large solvent collections')
    if manifests['medium']['cases']!=manifests['large']['cases']:
        raise InvalidArtifact('comparison physical preparations differ')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for case in CASES:
        r={'case_id':case,'evidence':LABELS[case]}
        for label,c in collections.items():
            score=c['scores'][case]
            r[label]={k:score.get(k) for k in ('status','MACE_vacuum_R_kcal_mol','GB_Ca_minus_La_kcal_mol','R_kcal_mol')}
        a,b=(r[label]['R_kcal_mol'] for label in ('medium','large'))
        r['large_minus_medium_R_kcal_mol']=None if a is None or b is None else b-a
        rows.append(r)
    result={'status':'complete' if all(c['status']=='complete' for c in collections.values()) else 'incomplete',
            'collections':{'medium':record(medium_gb),'large':record(large_gb)},'rows':rows,
            'numerical_checks':{label:c['numerical_checks_pass'] for label,c in collections.items()},
            'contrasts':{label:c['contrasts'] for label,c in collections.items()},
            'primary_candidate_ordering_screen_pass':collections['medium']['development_ordering_screen_pass'],
            'reference':None,'S_kcal_mol':None,'calibrated_class':None,'goal_completion_claimed':False}
    write_new(out/'comparison.json',result)
    lines=['# Whole-protein MACE / frozen-solvent development panel','',
           'Production baseline unchanged. Five consumed structures, four biological groups.',
           'Alpha structures are one qualified affinity observation; PQQ labels describe functional association.',
           'R = E(Ca) - E(La). Solvent R adds G(Ca) - G(La). All values below are kcal/mol.',
           'No aquo offset, inherited decision band, fitted threshold or absolute class.','',
           '| Case | Medium vacuum R | Medium solvent R | Large solvent R | Large−medium solvent R |',
           '|---|---:|---:|---:|---:|']
    for r in rows:
        values=[r['medium']['MACE_vacuum_R_kcal_mol'],r['medium']['R_kcal_mol'],r['large']['R_kcal_mol'],r['large_minus_medium_R_kcal_mol']]
        lines.append('| '+r['case_id']+' | '+' | '.join('unavailable' if v is None else f'{v:.6f}' for v in values)+' |')
    lines.extend(['','## Predeclared relative-order checks','','```json',json.dumps(result['contrasts'],indent=2),'```','',
                  f"Medium primary ordering screen passes: {result['primary_candidate_ordering_screen_pass']}.",
                  f"Numerical checks: {result['numerical_checks']}.",'',
                  'The ordering screen is a development result, not broad validation or goal completion.',
                  'Frozen monopoles omit dipolar/self-consistent solvent response; metal cavity parameters remain unvalidated.',
                  'No combined descriptor gradient, relaxation or entropy correction is available.',
                  'Exact inputs, exclusions, terminal/H operations and chemical states are pinned per case.',
                  'Per-task receipts retain evaluation/initialization time and memory; whole Slurm costs are reported separately.'])
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'status':result['status'],'comparison':record(out/'comparison.json'),'report':record(out/'REPORT.md')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare-mace')
    for key in ('preparation','medium-reference','large-reference','agreement','output'):q.add_argument('--'+key,required=True)
    q=sub.add_parser('prepare-gb')
    for key in ('collection','solver-validation','agreement','output'):q.add_argument('--'+key,required=True)
    q=sub.add_parser('compare')
    for key in ('medium-gb','large-gb','output'):q.add_argument('--'+key,required=True)
    a=p.parse_args()
    if a.command=='prepare-mace':r=prepare_mace(a.preparation,a.medium_reference,a.large_reference,a.agreement,a.output)
    elif a.command=='prepare-gb':r=prepare_gb(a.collection,a.solver_validation,a.agreement,a.output)
    else:r=compare_report(a.medium_gb,a.large_gb,a.output)
    print(json.dumps(r,indent=2))
