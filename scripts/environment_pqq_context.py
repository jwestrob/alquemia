"""Read-only compatibility inventory for the unchanged compact-context PQQ rule."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_hybrid import EV_TO_KCAL,check_atoms,write_xyz

PROTOCOL='native_OMOL_complete_polar_context_disulfide_v2'


def transfer_item(old,root):
    artifacts=old['baseline']['artifacts']
    if old['case_id']=='1KB0':
        release=read_json(root/'diagnostics/baseline_benchmark_20260915/RESULT.json')
        parent=next(r['source_manifest'] for r in release['scores'] if r['case']=='pqq_1kb0')
    else:
        directory=verify(artifacts['Ca']['xyz']).parent;parents=list(directory.glob('*carve_manifest.json'))
        if len(parents)!=1:raise InvalidArtifact('ambiguous transfer source preparation')
        parent=record(parents[0])
    p=read_json(verify(parent));eps={}
    for metal in ('Ca','La'):
        a=artifacts[metal]
        eps[metal]={'xyz':p['outputs'][metal+'_xyz'],'input':p['outputs'][metal+'_input'],
                    'output':a['output'],'receipt':a.get('execution',a.get('receipt')),
                    'charge':p['charge_ledger'][metal+'_total'],'energy_hartree':a['energy_hartree']}
    return {'case':old['case_id'],'parent':parent,'endpoints':eps,'group':old['sequence_accession_group'],
            'evidence_stratum':old['evidence_stratum'],'expected_class':old['expected_class'],
            'source_native_OMOL_R_kcal_mol':old['R_kcal_mol']}


def inventory(source,static_config,output,chemistry_version='v1'):
    from second_shell_context import parent_state,expansion
    if chemistry_version=='v2':
        from environment_context_chemistry import complete_expansion as expansion
    elif chemistry_version!='v1':raise InvalidArtifact('unsupported context chemistry version')
    report=read_json(source);cfg=read_json(static_config);rows=[]
    calibration=[r for r in report['canonical_rows'] if r['evaluation_role']=='calibration']
    if len(calibration)!=25:raise InvalidArtifact('canonical 25-case inventory differs')
    for old in calibration:
        artifacts=old['baseline']['artifacts'];directory=verify(artifacts['Ca']['xyz']).parent
        parents=list(directory.glob('*carve_manifest.json'))
        if len(parents)!=1:raise InvalidArtifact('ambiguous canonical source preparation')
        p=read_json(parents[0]);eps={}
        for metal in ('Ca','La'):
            a=artifacts[metal]
            eps[metal]={'xyz':a['xyz'],'input':a['input'],'output':a['output'],'receipt':a['execution'],
                        'charge':p['charge_ledger'][metal+'_total'],'energy_hartree':a['energy_hartree']}
        item={'case':old['case_id'],'parent':record(parents[0]),'endpoints':eps,'group':old['sequence_accession_group'],
              'evidence_stratum':old['evidence_stratum'],'expected_class':old['expected_class'],
              'source_native_OMOL_R_kcal_mol':old['R_kcal_mol']}
        try:
            state=parent_state(item,cfg['topology']);_,audit=expansion(state)
            rows.append({'case':item['case'],'status':'supported','source':item,'atoms_before':audit['original_atom_count'],
                         'atoms_after':audit['new_atom_count'],'added_formal_charge':audit['added_formal_charge'],
                         'added_fragments':audit['added_fragments'],'disulfide_closures':audit.get('disulfide_closures',[])})
        except (ValueError,KeyError,OSError) as exc:rows.append({'case':item['case'],'status':'unsupported','source':item,'reason':str(exc)})
    result={'source':record(source),'static_config':record(static_config),'policy':cfg['policy'],'cases':rows,
            'supported':sum(r['status']=='supported' for r in rows),'denominator':25,'new_energy_calls':0,
            'classification':None,'basis_for_new_threshold':None,'chemistry_version':chemistry_version}
    write_new(output,result);return result


def prepare(inventory,agreement,output):
    from second_shell_context import parent_state
    from environment_context_chemistry import complete_expansion,POLICY_ID
    import mace_omol as omol
    inv=read_json(inventory);cfg=read_json(verify(inv['static_config']));source=read_json(verify(inv['source']))
    root=verify(inv['source']).parents[2]
    if inv['supported']!=25 or inv['chemistry_version']!='v2':raise InvalidArtifact('all25 versioned preparations required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        x=impl/p.name;shutil.copyfile(p,x);pins[p.name]=record(x)
    static=root/'workspaces/second_shell_20260919/prepared_v2'
    old=read_json(static/'mace_collection_1202083.json');oldm=read_json(static/'mace_manifest.json')
    model=omol.model(verify(cfg['software']))
    if model!=oldm['model'] or model!=source['model']:raise InvalidArtifact('native model compatibility failed')
    items=[r['source']|{'evaluation_role':'calibration'} for r in inv['cases']]
    items += [transfer_item(r,root)|{'evaluation_role':'retrospective_structural_transfer'}
              for r in source['canonical_rows'] if r['evaluation_role']!='calibration']
    tasks=[];reused=[];cases=[]
    for item in items:
        state=parent_state(item,cfg['topology']);rows,audit=complete_expansion(state);case=item['case']
        d=out/'cases'/case;d.mkdir(parents=True);ap=d/'preparation.json';write_new(ap,audit)
        for metal in ('Ca','La'):
            charge=item['endpoints'][metal]['charge']+audit['added_formal_charge'];xp=d/(metal+'.xyz');write_xyz(xp,rows[metal])
            task={'task_id':case+'__expanded__'+metal,'case_id':case,'metal':metal,'metal_index':0,'kind':'core',
                  'variant':'expanded','energy_component':omol.COMPONENT,'energy_only':False,'charge':charge,
                  'spin_multiplicity':1,'xyz':record(xp),'state':check_atoms(rows[metal],charge),'preparation':record(ap)}
            task['cache_key']=cache_key({'task':task,'model':model,'software':cfg['software'],'implementation':pins})
            if case in ('1H4I','4MAE'):
                prior=next(t for t in oldm['tasks'] if t['case_id']==case and t['variant']=='expanded' and t['metal']==metal)
                if prior['xyz']['sha256']!=task['xyz']['sha256'] or prior['charge']!=charge:raise InvalidArtifact('transfer reuse differs')
                accepted=next(r['accepted'] for r in old['rows'] if r['case']==case and r['variant']=='expanded' and r['metal']==metal)
                if not omol.accepted_state(accepted,prior):raise InvalidArtifact('prior native endpoint invalid')
                reused.append({'task':task,'source_collection':record(static/'mace_collection_1202083.json'),'accepted':accepted})
            else:tasks.append(task)
        cases.append({'source':item,'preparation':record(ap),'status':'prepared'})
    if len(tasks)!=52 or len(cases)!=28 or len(reused)!=4:raise InvalidArtifact('finite52task/28case scope differs')
    m={'protocol_id':PROTOCOL,'stage':'environment_pqq_context','inventory':record(inventory),'agreement':record(agreement),
       'source':inv['source'],'source_static_collection':record(static/'mace_collection_1202083.json'),
       'chemistry_policy':POLICY_ID,'policy':cfg['policy'],'model':model,'software':cfg['software'],'implementation':pins,
       'tasks':tasks,'reused':reused,'cases':cases,'calibration_rule':{'minimum_gap_model_kcal_mol':.02,'bands':'calibration_class_extrema'},
       'compute_budget':None,'reference':None,'baseline_changed':False,'new_DFT_endpoints':0}
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def validate(manifest):
    import mace_omol as omol
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or len(m['tasks'])!=52 or len(m['reused'])!=4:raise InvalidArtifact('frozen panel differs')
    for p in m['implementation'].values():verify(p)
    for key in ('inventory','agreement','source','source_static_collection'):verify(m[key])
    if m['model']!=omol.model(verify(m['software'])):raise InvalidArtifact('model changed')
    for t in m['tasks']+[r['task'] for r in m['reused']]:
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('cache changed')
        if check_atoms(xyz(verify(t['xyz'])),t['charge'])!=t['state']:raise InvalidArtifact('state changed')
        verify(t['preparation'])
    return {'status':'pass','new_native_MACE_endpoints':52,'reused_context_endpoints':4,'cases':28,'DFT_endpoints':0}


def execute(manifest):
    import mace_omol as omol
    if not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('allocation required')
    validate(manifest);m=read_json(manifest);out=Path(manifest).parent/'mace_execution';out.mkdir(exist_ok=False);start=time.monotonic();pins=[]
    for t in m['tasks']:
        d=out/t['task_id'];d.mkdir();r=omol.worker(manifest,t['task_id'],d,'native')
        if r['status']!='computed' or not omol.accepted_state(r,t):raise InvalidArtifact('native MACE endpoint failed')
        pins.append(record(d/'result.json'))
    col={'manifest':record(manifest),'results':pins,'wall_seconds':time.monotonic()-start,'slurm_job_id':os.environ['SLURM_JOB_ID'],
         'allocated_cpus':int(os.environ['SLURM_CPUS_PER_TASK']),'status':'complete'}
    path=Path(manifest).parent/('collection_'+os.environ['SLURM_JOB_ID']+'.json');write_new(path,col)
    return {'status':'complete','collection':record(path)}


def collect(collection,output):
    col=read_json(collection);m=read_json(verify(col['manifest']));energies={}
    for pin in col['results']:
        r=read_json(verify(pin));energies[r['task_id']]=r['energy_eV']
    for r in m['reused']:energies[r['task']['task_id']]=r['accepted']['energy_eV']
    rows=[]
    for case in m['cases']:
        s=case['source'];cid=s['case'];prep=read_json(verify(case['preparation']))
        r=(energies[cid+'__expanded__Ca']-energies[cid+'__expanded__La'])*EV_TO_KCAL
        rows.append({'case':cid,'role':s['evaluation_role'],'group':s['group'],'expected_class':s['expected_class'],
                     'core_R_model_kcal_mol':s['source_native_OMOL_R_kcal_mol'],'context_R_model_kcal_mol':r,
                     'added_formal_charge':prep['added_formal_charge'],'atoms':prep['new_atom_count'],'decision':None})
    calibration=[r for r in rows if r['role']=='calibration'];ca=max(r['context_R_model_kcal_mol'] for r in calibration if r['expected_class']=='Ca')
    la=min(r['context_R_model_kcal_mol'] for r in calibration if r['expected_class']=='La');gap=la-ca
    passed=len(calibration)==25 and gap>m['calibration_rule']['minimum_gap_model_kcal_mol']
    if passed:
        for r in rows:
            v=r['context_R_model_kcal_mol'];r['decision']='Ca-supported' if v<=ca else 'La-supported' if v>=la else 'inconclusive'
            r['expected_region_pass']=r['decision']==r['expected_class']+'-supported'
    result={'manifest':col['manifest'],'collection':record(collection),'rows':rows,'calibration_count':25,'calibration_pass':passed,
            'gap_model_kcal_mol':gap,'Ca_band_max':ca if passed else None,'La_band_min':la if passed else None,
            'correct_calibration':sum(r.get('expected_region_pass',False) for r in calibration),
            'correct_transfer':sum(r.get('expected_region_pass',False) for r in rows if r['role']!='calibration'),
            'transfer_denominator':3,'baseline_changed':False,'all_cases_consumed':True,'broad_affinity_validated':False}
    write_new(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('inventory')
    for s in ('source','static-config','output'):q.add_argument('--'+s,required=True)
    q.add_argument('--chemistry-version',choices=('v1','v2'),default='v1')
    q=sub.add_parser('prepare')
    for s in ('inventory','agreement','output'):q.add_argument('--'+s,required=True)
    for s in ('validate','execute'):
        q=sub.add_parser(s);q.add_argument('--manifest',required=True)
    q=sub.add_parser('collect');q.add_argument('--collection',required=True);q.add_argument('--output',required=True)
    argv=sys.argv[1:]
    if argv and argv[0].startswith('--'):argv=['inventory']+argv
    a=vars(p.parse_args(argv));op=a.pop('command');r=globals()[op](**a)
    print(json.dumps(r if op!='inventory' else {k:r[k] for k in ('supported','denominator','new_energy_calls')},indent=2))
