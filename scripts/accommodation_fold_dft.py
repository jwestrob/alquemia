"""Thin frozen-fold DFT adapter around the existing manifested ORCA runner."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import datetime as dt
import json
from pathlib import Path
import re
import shutil
import statistics

from affordable_common import HA_TO_KCAL,InvalidArtifact,energy,read_json,record,verify,write_new,xyz
from result_protocol import verify_orca_execution_record
from run_orca_task_manifest import load_manifest_tasks,run_manifest,_verify_prepared_task

PROTOCOL='pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'
FILES=('accommodation_fold_dft.py','affordable_common.py','run_orca_task_manifest.py',
       'render_orca_runtime_input.py','result_protocol.py')


def recipe(path,coordinate_path):
    lines=[s.strip() for s in Path(path).read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    if len(lines)!=3 or lines[0].split()!=['!','r2SCAN-3c','NoAutostart','CPCM(Water)','DefGrid3'] or lines[1]!='%maxcore 8000':
        raise InvalidArtifact('not the preserved native baseline recipe: '+str(path))
    m=re.fullmatch(r'\* xyzfile (-?\d+) (\d+) (\S+)',lines[2])
    if not m or (Path(path).parent/m[3]).resolve()!=Path(coordinate_path).resolve():
        raise InvalidArtifact('input XYZ identity differs: '+str(path))
    return {'charge':int(m[1]),'multiplicity':int(m[2]),'method':'native_r2SCAN3c_CPCM_Water_DefGrid3_NormalSCF',
            'maxcore_MB_per_rank':8000}


def elapsed(receipt):
    return (dt.datetime.fromisoformat(receipt['finished_at_utc'])-dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds()


def checked_endpoint(input_pin,xyz_pin,output_pin,receipt_pin,manifest_pin,task_id):
    ip,xp,op,rp,mp=(verify(p) for p in (input_pin,xyz_pin,output_pin,receipt_pin,manifest_pin))
    problems=verify_orca_execution_record(output_path=op,manifest_path=mp,task_id=task_id,input_path=ip,xyz_path=xp)
    if problems:raise InvalidArtifact('; '.join(problems))
    rr=read_json(rp)
    if rr['orca_version']!='6.1.1':raise InvalidArtifact('ORCA version differs')
    raw=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+\d.EeDd]+)',op.read_text())[-1]
    return {'status':'complete','energy_hartree':energy(op),'energy_hartree_text':raw,
            'input':input_pin,'xyz':xyz_pin,'output':output_pin,'receipt':receipt_pin,
            'manifest':manifest_pin,'task_id':task_id,'wall_seconds':elapsed(rr),
            'nprocs':rr['parallelism']['nprocs'],'endpoint_rank_seconds':elapsed(rr)*rr['parallelism']['nprocs'],
            'allocation':rr['allocation'],'orca_executable':rr['orca_executable']}


def prepare(preparation,reference,agreement,output):
    prep=read_json(preparation);ref=read_json(reference);out=Path(output).resolve()
    if prep['denominator']!=250 or prep['supported']!=233 or ref['protocol_id']!=PROTOCOL:raise InvalidArtifact('scope/protocol differs')
    src=read_json(verify(prep['source_manifest']));sources={r['case_id']:r for r in src['cases']}
    archived={r['panel_id']:r for r in ref['scores']}
    oldprep=read_json(verify(ref['preparation']));parents={r['panel_id']:r for r in oldprep['targets']}
    if len(sources)!=250 or len(archived)!=25:raise InvalidArtifact('source denominator differs')
    # Complete scientific and archival audit before creating a task tree.
    reused={};cost=[];orc=None;audits=[]
    for c in prep['cases']:
        s=sources[c['case_id']]
        if c['status']!='prepared':continue
        old=archived[s['root_case_id']];cp=read_json(verify(c['core']['parent']))
        if cp['protocol_id']!=PROTOCOL or cp['fixed_core']['water_policy']!='dry_exclude_all_source_and_synthetic_waters':
            raise InvalidArtifact('core chemistry/protocol differs')
        pair={}
        for z in ('Ca','La'):
            ep=c['core']['endpoints'][z];a=old['artifacts'][z]
            method=recipe(verify(ep['input']),verify(ep['xyz']));oldmethod=recipe(verify(a['input']),verify(a['xyz']))
            if method!=oldmethod or ep['charge']!=method['charge'] or ep['multiplicity']!=method['multiplicity']:
                raise InvalidArtifact('native recipe/state differs')
            newxyz=xyz(verify(ep['xyz']));oldxyz=xyz(verify(a['xyz']))
            if [x[0] for x in newxyz]!=[x[0] for x in oldxyz]:raise InvalidArtifact('fixed core composition/order differs')
            pair[z]=newxyz
            if s['canonical_coordinate_match']:
                delta=max(abs(a[i]-b[i]) for a,b in zip(newxyz,oldxyz) for i in (1,2,3))
                if delta>1e-12:raise InvalidArtifact('canonical coordinates do not match archive')
                epold=checked_endpoint(a['input'],a['xyz'],a['output'],a['execution'],parents[s['root_case_id']]['manifest'],z)
                if epold['energy_hartree']!=old['energies_hartree'][z]:raise InvalidArtifact('archived energy differs from release')
                if orc is None:orc=record(verify(epold['orca_executable']))
                elif epold['orca_executable']!=orc:raise InvalidArtifact('archived executable mismatch')
                epold['reused']=True;epold['current_core_xyz']=ep['xyz'];epold['canonical_max_displacement_A']=delta
                reused.setdefault(c['case_id'],{})[z]=epold
                cost.append({'case_id':c['case_id'],'metal':z,**{k:epold[k] for k in ('wall_seconds','nprocs','endpoint_rank_seconds','allocation','receipt')}})
        if pair['Ca'][1:]!=pair['La'][1:] or pair['Ca'][0][1:]!=pair['La'][0][1:]:raise InvalidArtifact('paired core coordinates differ')
        audits.append({'case_id':c['case_id'],'source_core':c['core']['parent'],'atom_count':len(pair['Ca']),
                       'canonical':s['canonical_coordinate_match'],'paired_coordinates_exact':True,'native_recipe_matches_released_parent':True})
    if len(reused)!=25 or len(cost)!=50 or len(audits)!=233:raise InvalidArtifact('supported/reused count differs')
    out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir()
    for name in FILES:shutil.copy2(Path(__file__).parent/name,impl/name)
    policy={'workers':4,'nprocs':16,'task_runner':record(impl/'run_orca_task_manifest.py'),
            'runtime_renderer':record(impl/'render_orca_runtime_input.py')}
    primary=[c for c in prep['cases'] if c['status']=='prepared' and not sources[c['case_id']]['canonical_coordinate_match']]
    primary.sort(key=lambda c:c['case_id'])
    if len(primary)!=208:raise InvalidArtifact('finite new-source count differs')
    manifests=[];newtasks=[]
    for shard in range(4):
        sd=out/f'shard_{shard}';sd.mkdir();tasks=[]
        for c in primary[shard::4]:
            for z in ('Ca','La'):
                ep=c['core']['endpoints'][z];td=sd/'tasks'/c['case_id']/z;td.mkdir(parents=True)
                ip=verify(ep['input']);xp=verify(ep['xyz']);shutil.copy2(ip,td/ip.name);shutil.copy2(xp,td/xp.name)
                task={'task_id':c['case_id']+'__'+z,'case_id':c['case_id'],'metal':z,
                      'input':record(td/ip.name),'xyz':record(td/xp.name),'output_path':str(td/'endpoint.out'),
                      'source_core':c['core']['parent'],'source_endpoint':ep}
                if task['input']['sha256']!=ep['input']['sha256'] or task['xyz']['sha256']!=ep['xyz']['sha256']:raise InvalidArtifact('copy changed')
                tasks.append(task)
        if len(tasks)!=104:raise InvalidArtifact('shard size differs')
        m={'protocol_id':PROTOCOL,'agreement':record(agreement),'preparation':record(preparation),
           'execution_policy':policy,'orca':orc,'shard':shard,'tasks':tasks}
        mp=sd/'manifest.json';write_new(mp,m);manifests.append(record(mp))
        for t in tasks:newtasks.append(dict(t,manifest=record(mp)))
    result={'schema_version':'accommodation_fold_DFT_v1','protocol_id':PROTOCOL,'agreement':record(agreement),
            'preparation':record(preparation),'sources':prep['source_manifest'],'reference':record(reference),
            'reference_protocol':ref['protocol_id'],'bands':{'Ca_max':ref['calibration']['U_max_Ca_kcal_mol'],
            'La_min':ref['calibration']['L_min_La_kcal_mol']},'aquo_reporting_gauge':ref['aquo_reporting_gauge'],
            'denominator':250,'primary_denominator':225,'new_source_count':208,'new_endpoint_count':416,
            'reused_canonical_pairs':reused,'new_tasks':newtasks,'shards':manifests,'orca':orc,
            'execution_policy':policy,'core_audits':audits,'historical_endpoint_costs':cost,
            'historical_cost_summary':{'endpoint_count':50,'sum_wall_seconds':sum(r['wall_seconds'] for r in cost),
                 'median_wall_seconds':statistics.median(r['wall_seconds'] for r in cost),
                 'sum_endpoint_rank_seconds':sum(r['endpoint_rank_seconds'] for r in cost),
                 'matched_scientific_protocol':True,'matched_future_hardware':False},
            'implementation':{n:record(impl/n) for n in FILES},'production_changed':False,'threshold_refitted':False}
    write_new(out/'manifest.json',result);return dry_run(out/'manifest.json')


def dry_run(manifest):
    m=read_json(manifest);verify(m['orca']);seen=set()
    for shard in m['shards']:
        _,tasks=load_manifest_tasks(verify(shard))
        for task in tasks:
            if task['task_id'] in seen:raise InvalidArtifact('overlapping shards')
            seen.add(task['task_id']);_verify_prepared_task(task);recipe(task['input'],task['xyz'])
    if len(seen)!=416:raise InvalidArtifact('new endpoint count differs')
    return {'status':'ready','new_endpoints':416,'shards':4,'endpoints_per_shard':104,'reused_pairs':25,
            'unavailable_sources':17,'manifest':record(manifest),'historical_costs':m['historical_cost_summary']}


def decision(value,bands):
    if value is None:return 'unavailable'
    if value<=bands['Ca_max']:return 'Ca-supported'
    if value>=bands['La_min']:return 'La-supported'
    return 'inconclusive'


def counts(rows):
    calls=Counter(r['decision'] for r in rows)
    return {'denominator':len(rows),'correct':sum(r['decision']==r['expected_class']+'-supported' for r in rows),
            'wrong':sum(r['decision'] in ('Ca-supported','La-supported') and r['decision']!=r['expected_class']+'-supported' for r in rows),
            'inconclusive':calls['inconclusive'],'unavailable':calls['unavailable'],'calls':dict(calls)}


def collect(manifest,output,final_if_complete=False):
    m=read_json(manifest);prep=read_json(verify(m['preparation']));source=read_json(verify(m['sources']))
    supported={r['case_id']:r for r in prep['cases']};endpoints=dict(m['reused_canonical_pairs']);progress=Counter()
    for t in m['new_tasks']:
        op=Path(t['output_path']);rp=Path(str(op)+'.execution.json')
        if not rp.exists():ep={'status':'pending','output_path':str(op),'receipt_path':str(rp)}
        else:
            try:
                ep=checked_endpoint(t['input'],t['xyz'],record(op),record(rp),t['manifest'],t['task_id']);ep['reused']=False
            except Exception as exc:ep={'status':'failed','reason':str(exc),'receipt':record(rp),'output':record(op) if op.exists() else None}
        endpoints.setdefault(t['case_id'],{})[t['metal']]=ep;progress[ep['status']]+=1
    rows=[];groups=defaultdict(list)
    for c in source['cases']:
        row={k:c[k] for k in ('case_id','root_case_id','source_conditioning_metal','canonical_coordinate_match',
                              'primary_evaluation_pool','biological_group','expected_class','evidence_stratum')}
        row.update(preparation_status=supported[c['case_id']]['status'],preparation_reason=supported[c['case_id']].get('reason'),
                   endpoints=endpoints.get(c['case_id'],{}),R_hartree=None,R_kcal_mol=None,S_aquo_gauge_kcal_mol=None)
        if len(row['endpoints'])==2 and all(e['status']=='complete' for e in row['endpoints'].values()):
            r=row['endpoints']['Ca']['energy_hartree']-row['endpoints']['La']['energy_hartree']
            row.update(R_hartree=r,R_kcal_mol=r*HA_TO_KCAL,S_aquo_gauge_kcal_mol=(r-m['aquo_reporting_gauge']['delta_E_aquo_hartree'])*HA_TO_KCAL)
        row['decision']=decision(row['R_kcal_mol'],m['bands']);rows.append(row);groups[c['root_case_id']].append(row)
    protein=[]
    for root,rr in groups.items():
        ps={}
        for name,arm in [('La4',[r for r in rr if r['source_conditioning_metal']=='La' and not r['canonical_coordinate_match']]),
                         ('Ca5',[r for r in rr if r['source_conditioning_metal']=='Ca'])]:
            if len(arm)!=(4 if name=='La4' else 5):raise InvalidArtifact('arm size differs')
            vals=[r['R_kcal_mol'] for r in arm];value=statistics.median(vals) if all(v is not None for v in vals) else None
            ps[name]={'R_kcal_mol':value,'decision':decision(value,m['bands']),
                      'declared_members':[r['case_id'] for r in arm],'available_members':sum(v is not None for v in vals),
                      'expected_class':rr[0]['expected_class']}
        la,ca=ps['La4']['R_kcal_mol'],ps['Ca5']['R_kcal_mol'];balanced=None if la is None or ca is None else (la+ca)/2
        ps['balanced']={'R_kcal_mol':balanced,'decision':decision(balanced,m['bands']),'expected_class':rr[0]['expected_class']}
        protein.append({'root_case_id':root,'biological_group':rr[0]['biological_group'],'expected_class':rr[0]['expected_class'],
                        'pools':ps,'Ca5_median_minus_La4_median_kcal_mol':None if la is None or ca is None else ca-la})
    subsets={'all250':rows,'canonical25':[r for r in rows if r['canonical_coordinate_match']],
             'primary225':[r for r in rows if r['primary_evaluation_pool']],
             'La100':[r for r in rows if r['source_conditioning_metal']=='La' and not r['canonical_coordinate_match']],
             'Ca125':[r for r in rows if r['source_conditioning_metal']=='Ca']}
    computed=[e for rr in endpoints.values() for e in rr.values() if e['status']=='complete' and not e.get('reused')]
    result={'schema_version':'accommodation_fold_DFT_collection_v1','manifest':record(manifest),
            'preparation':m['preparation'],'sources':m['sources'],'reference':m['reference'],'bands':m['bands'],
            'new_endpoint_statuses':dict(progress),'rows':rows,'protein_summaries':protein,
            'counts':{'single_sources':{k:counts(v) for k,v in subsets.items()},
                      'ensemble_descriptors':{k:counts([p['pools'][k] for p in protein]) for k in ('La4','Ca5','balanced')}},
            'measured_new_completed_endpoint_cost':{'count':len(computed),'wall_seconds_sum':sum(e['wall_seconds'] for e in computed),
                'endpoint_rank_seconds_sum':sum(e['endpoint_rank_seconds'] for e in computed)},
            'all_evidence_consumed':True,'threshold_refitted':False,'production_changed':False,
            'ensemble_interpretation':'fixed structural robustness descriptors, not thermal ensembles or newly calibrated classifier'}
    write_new(output,result);report(result,Path(output).with_suffix('.md'))
    if final_if_complete and progress['pending']==0:
        final=Path(manifest).parent/'final_collection.json'
        try:write_new(final,result)
        except FileExistsError:pass
        else:report(result,final.with_suffix('.md'))
    return {'new_endpoint_statuses':dict(progress),'counts':result['counts'],'output':str(output)}


def report(result,output):
    lines=['# Preserved DFT baseline: fold-conditioning challenge','',
           'All250 sources remain in the denominator. Canonical25 pairs are verified archive reuse; primary225 are structural repeats of consumed reference proteins. No fitting.', '',
           '| Descriptor | Correct | Wrong | Inconclusive | Unavailable | Denominator |',
           '|---|---:|---:|---:|---:|---:|']
    for section in result['counts'].values():
        for name,c in section.items():lines.append(f"| {name} | {c['correct']} | {c['wrong']} | {c['inconclusive']} | {c['unavailable']} | {c['denominator']} |")
    lines+=['',f"New endpoints: {result['new_endpoint_statuses']}.",
            f"Measured completed endpoint costs: {result['measured_new_completed_endpoint_cost']}.",
            '', 'La4/Ca5 require every declared member; balanced is the equal mean of those arm medians. Missing members make the descriptor unavailable. Frozen original single-source bands are used without recalibration. S is aquo-gauge display only.', '']
    with Path(output).open('x') as f:f.write('\n'.join(lines))


def compare(dft,mace,output):
    d=read_json(dft);m=read_json(mace)
    if d['preparation']!=m['preparation'] or d['sources']!=m['sources']:raise InvalidArtifact('DFT/MACE source preparations differ')
    lookup={r['case_id']:r for r in d['rows']};paired=[]
    for r in m['rows']:
        old=lookup[r['case_id']]
        if (old['expected_class'],old['biological_group'])!=(r['expected_class'],r['biological_group']):raise InvalidArtifact('label/group differs')
        paired.append({'case_id':r['case_id'],'representation':r['representation'],'DFT_R_kcal_mol':old['R_kcal_mol'],
                       'DFT_decision':old['decision'],'MACE_native_R_model_kcal_mol':r['native_R_model_kcal_mol'],
                       'MACE_native_decision':r['native_decision'],'MACE_composite_R_model_kcal_mol':r['composite_R_model_kcal_mol'],
                       'MACE_composite_decision':r['composite_decision']})
    result={'DFT_collection':record(dft),'MACE_comparison':record(mace),'preparation':d['preparation'],
            'DFT_counts':d['counts'],'MACE_counts':m['counts'],'paired_rows':paired,'threshold_refitted':False,'production_changed':False}
    write_new(output,result);return {'output':str(output),'DFT_counts':d['counts'],'MACE_counts':m['counts']}


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare')
    for key in ('preparation','reference','agreement','output'):q.add_argument('--'+key,required=True)
    q=sub.add_parser('dry-run');q.add_argument('--manifest',required=True)
    q=sub.add_parser('execute');q.add_argument('--manifest',required=True)
    q=sub.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True);q.add_argument('--final-if-complete',action='store_true')
    q=sub.add_parser('compare')
    for key in ('dft','mace','output'):q.add_argument('--'+key,required=True)
    a=vars(p.parse_args());op=a.pop('op')
    if op=='execute':
        m=read_json(a['manifest']);result=run_manifest(Path(a['manifest']),orca_path=verify(m['orca']),workers=4,nprocs=16,
                           expected_runner_sha256=m['execution_policy']['task_runner']['sha256'])
    else:result=globals()[op.replace('-','_')](**a)
    print(json.dumps(result))


if __name__=='__main__':main()
