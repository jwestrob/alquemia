"""Unchanged static second-shell physics on two archived GGR structures."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import time
from affordable_common import HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from hydration_square import endpoint
from mace_hybrid import EV_TO_KCAL,check_atoms,write_xyz

CASES=('2FW0','2FVY')
PROTOCOL='native_r2scan3c_CPCM_static_second_shell_GGR_replication_v1'


def configuration(root,agreement,output):
    root=Path(root).resolve();src=root/'workspaces/second_shell_20260919/prepared_v2/manifest.json'
    old=read_json(src);cfg=read_json(verify(old['configuration']));cases=[]
    for case in CASES:
        c=case.lower();p=root/f'workspaces/ggr_mechanism_20260915/stage_b_prepared_v1/{c}/amide_v3/repair_manifest.json'
        parent=read_json(p);eps={}
        for metal in ('Ca','La'):
            ep=parent['outputs'][metal];op=root/f'workspaces/ggr_mechanism_20260915/stage_b_tasks_v1/ggr_{c}_formamide/{metal}/ggr_{c}_GGR_amide_v3_{metal}.out'
            e=endpoint(record(op),record(str(op)+'.execution.json'),ep['xyz'],ep['input'])
            eps[metal]={'xyz':ep['xyz'],'input':ep['input'],'output':record(op),'receipt':record(str(op)+'.execution.json'),
                        'charge':ep['charge'],'energy_hartree':e['energy_hartree']}
        cases.append({'case':case,'parent':record(p),'endpoints':eps,'group':'GGR',
                      'evidence_stratum':'condition_qualified_affinity_direction'})
    new={k:cfg[k] for k in ('policy','topology','orca','execution_policy','software')}
    new.update(protocol_id=PROTOCOL,cases=cases,agreement=record(agreement),source_static_manifest=record(src),
               source_static_collection=record(src.parent/'collection_1202082.json'),
               source_static_MACE=record(src.parent/'mace_collection_1202083.json'))
    write_new(output,new);return new


def prepare(config,output):
    from second_shell_context import parent_state,expansion,POLICY
    import mace_omol as omol
    cfg=read_json(config)
    if tuple(c['case'] for c in cfg['cases'])!=CASES or cfg['policy']!=POLICY:raise InvalidArtifact('frozen scope differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        x=impl/p.name;shutil.copyfile(p,x);pins[p.name]=record(x)
    tasks=[];states=[];mtasks=[];model=omol.model(verify(cfg['software']))
    for item in cfg['cases']:
        state=parent_state(item,cfg['topology']);rows,audit=expansion(state);case=item['case']
        ap=out/(case+'_preparation.json');write_new(ap,audit);states.append({'case':case,'source':item,'preparation':record(ap)})
        for metal in ('Ca','La'):
            ep=item['endpoints'][metal];charge=ep['charge']+audit['added_formal_charge'];tid=case+'__expanded__'+metal
            d=out/'tasks'/tid;d.mkdir(parents=True);xp=d/'core.xyz';write_xyz(xp,rows[metal]);ip=d/'endpoint.inp'
            ip.write_text('! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n'+f'* xyzfile {charge} 1 core.xyz\n')
            check_atoms(rows[metal],charge)
            tasks.append({'task_id':tid,'case':case,'metal':metal,'charge':charge,'multiplicity':1,'input':record(ip),
                          'xyz':record(xp),'output_path':str(d/'endpoint.out'),'preparation':record(ap),'atom_count':len(rows[metal])})
            for variant,xpin,q in [('core',ep['xyz'],ep['charge']),('expanded',record(xp),charge)]:
                mtasks.append({'task_id':case+'__'+variant+'__'+metal,'case_id':case,'metal':metal,'metal_index':0,
                    'kind':'core','variant':variant,'energy_component':omol.COMPONENT,'energy_only':False,
                    'charge':q,'spin_multiplicity':1,'xyz':xpin,'state':check_atoms(xyz(verify(xpin)),q)})
    m={'protocol_id':PROTOCOL,'agreement':cfg['agreement'],'configuration':record(config),'policy':POLICY,
       'states':states,'tasks':tasks,'implementation':pins,'orca':cfg['orca'],'execution_policy':cfg['execution_policy'],
       'execution_resources':{'mpi_ranks':16,'concurrent_tasks':4},'reference':None,'calibrated_decision':None,
       'baseline_changed':False,'compute_budget':None,'new_DFT_endpoints':4}
    write_new(out/'manifest.json',m)
    mm={'protocol_id':'native_OMOL_static_second_shell_GGR_replication_v1','stage':'environment_replica_context',
        'agreement':cfg['agreement'],'source_DFT_manifest':record(out/'manifest.json'),'implementation':pins,
        'model':model,'software':cfg['software'],'tasks':mtasks}
    for t in mtasks:t['cache_key']=cache_key({'task':t,'model':model,'software':cfg['software'],'implementation':pins})
    write_new(out/'mace_manifest.json',mm);return validate(out/'manifest.json')


def validate(manifest):
    from second_shell_context import parent_state,expansion,POLICY
    m=read_json(manifest);cfg=read_json(verify(m['configuration']));verify(m['agreement'])
    if m['protocol_id']!=PROTOCOL or m['policy']!=POLICY or len(m['tasks'])!=4:raise InvalidArtifact('frozen scope changed')
    for p in m['implementation'].values():verify(p)
    inventory=[]
    for s in m['states']:
        rows,audit=expansion(parent_state(s['source'],cfg['topology']))
        if json.loads(json.dumps(audit))!=read_json(verify(s['preparation'])):raise InvalidArtifact('regeneration differs')
        for metal in ('Ca','La'):
            t=next(t for t in m['tasks'] if t['case']==s['case'] and t['metal']==metal)
            if xyz(verify(t['xyz']))!=rows[metal]:raise InvalidArtifact('source coordinates differ')
            check_atoms(rows[metal],t['charge'])
        inventory.append({'case':s['case'],'before_atoms':audit['original_atom_count'],'expanded_atoms':audit['new_atom_count'],
                          'added_formal_charge':audit['added_formal_charge']})
    from affordable_workflow import dry_run
    return dry_run(manifest)|{'inventory':inventory}


def execute_mace(manifest):
    import mace_omol as omol
    if not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('allocation required')
    m=read_json(manifest);out=Path(manifest).parent/'mace_execution';out.mkdir(exist_ok=False);start=time.monotonic();pins=[]
    if len(m['tasks'])!=8 or m['model']!=omol.model(verify(m['software'])):raise InvalidArtifact('native MACE scope/model differs')
    for p in m['implementation'].values():verify(p)
    for t in m['tasks']:
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('cache differs')
        d=out/t['task_id'];d.mkdir();r=omol.worker(manifest,t['task_id'],d,'native')
        if r['status']!='computed' or not omol.accepted_state(r,t):raise InvalidArtifact('native MACE endpoint failed')
        pins.append(record(d/'result.json'))
    col={'manifest':record(manifest),'results':pins,'wall_seconds':time.monotonic()-start,'slurm_job_id':os.environ['SLURM_JOB_ID'],
         'allocated_cpus':int(os.environ['SLURM_CPUS_PER_TASK']),'status':'complete'}
    write_new(Path(manifest).parent/('mace_collection_'+os.environ['SLURM_JOB_ID']+'.json'),col);return col


def collect(manifest,output):
    m=read_json(manifest);cfg=read_json(verify(m['configuration']));old=read_json(verify(cfg['source_static_collection']));rows=[]
    for s in m['states']:
        eps={}
        for metal in ('Ca','La'):
            t=next(t for t in m['tasks'] if t['case']==s['case'] and t['metal']==metal)
            eps[metal]=endpoint(record(t['output_path']),record(t['output_path']+'.execution.json'),t['xyz'],t['input'])
        item=s['source'];prep=read_json(verify(s['preparation']));b=item['endpoints']['Ca']['energy_hartree']-item['endpoints']['La']['energy_hartree']
        a=eps['Ca']['energy_hartree']-eps['La']['energy_hartree']
        rows.append({'case':s['case'],'group':'GGR','status':'complete','endpoints':eps,'core_R_hartree':b,'expanded_R_hartree':a,
                     'delta_R_kcal_mol':(a-b)*HA_TO_KCAL,'added_formal_charge':prep['added_formal_charge'],
                     'atoms_before':prep['original_atom_count'],'atoms_after':prep['new_atom_count']})
    lookup={r['case']:r for r in old['rows']+rows};matrix=[]
    for a in ('1F6S','6IP9'):
        for b in ('1GLG','2FW0','2FVY'):
            x,y=lookup[a],lookup[b]
            matrix.append({'La_like':a,'Ca_like':b,'before_kcal_mol':(x['core_R_hartree']-y['core_R_hartree'])*HA_TO_KCAL,
                           'after_kcal_mol':(x['expanded_R_hartree']-y['expanded_R_hartree'])*HA_TO_KCAL})
    result={'manifest':record(manifest),'source_static':cfg['source_static_collection'],'new_rows':rows,'matrix':matrix,
            'before_correct':sum(r['before_kcal_mol']>0 for r in matrix),'after_correct':sum(r['after_kcal_mol']>0 for r in matrix),
            'structural_directions':6,'independent_new_biological_groups':0,'reference':None,'calibrated_decision':None,'baseline_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('configuration')
    for s in ('root','agreement','output'):q.add_argument('--'+s,required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    for s in ('validate','execute-mace'):
        q=sub.add_parser(s);q.add_argument('--manifest',required=True)
    q=sub.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('command').replace('-','_');print(json.dumps(globals()[op](**a),indent=2))


if __name__=='__main__':main()
