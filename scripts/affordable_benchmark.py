"""Frozen baseline/amide benchmark using the existing ORCA task runner."""
from pathlib import Path
from collections import Counter
import argparse
import json
import re
import shutil
import statistics
from datetime import datetime
from affordable_common import (read_json, record, verify, write_new, energy, paired,
                               contrast, classify_raw, cache_key, InvalidArtifact)


def audit_reference(root):
    p=root/'reference_inputs/aquo_cn8_native_r2scan3c_v2/aquo_reference.json'
    ref=read_json(p)
    if ref['status']!='complete': raise InvalidArtifact('reference incomplete')
    artifacts={}
    for name,c in ref['calculations'].items():
        files={}
        for key,hkey in [('input','input_sha256'),('xyz','xyz_sha256'),('expected_output','output_sha256')]:
            q=verify({'path':str(p.parent/c[key]),'sha256':c[hkey]});files[key]=record(q)
        if energy(Path(files['expected_output']['path']))!=c['energy_hartree']:
            raise InvalidArtifact('reference energy mismatch')
        receipt=Path(files['expected_output']['path']+'.execution.json')
        r=read_json(receipt)
        if not r['normal_termination'] or not r['scf_converged']: raise InvalidArtifact('invalid reference receipt')
        files['receipt']=record(receipt);artifacts[name]=files
    e=ref['energies_hartree']
    gap=((e['Ca_on_Ca_geometry']-e['La_on_Ca_geometry'])+(e['Ca_on_La_geometry']-e['La_on_La_geometry']))/2
    if abs(gap-ref['delta_E_aquo_hartree'])>1e-10: raise InvalidArtifact('reference algebra mismatch')
    return {'manifest':record(p),'artifacts':artifacts,'delta_E_aquo_hartree':gap,
            'role':'common_ranking_gauge_no_inherited_v3_threshold'}


def prepare(root,output):
    root=root.resolve();output=output.resolve()
    if not output.is_relative_to(root/'workspaces'): raise InvalidArtifact('workspace required')
    if output.exists(): raise InvalidArtifact('refusing existing preparation')
    ref=audit_reference(root)
    inv=read_json(root/'diagnostics/affordable_challenger_20260915/repair_inventory_verified.json')
    sources=[]
    for entry in inv:
        case=Path(entry['repaired']['path']).parent.name
        for lane in ('baseline','repaired'): sources.append((case,lane,verify(entry[lane])))
    sources.append(('pqq_1kb0','fixed_core',root/'diagnostics/pqq_q46444_1kb0_external_validation_20260915/prepared/01_1KB0/pmdh_fc_holdout_1kb0_carve_manifest.json'))
    tasks=[];cases=[]
    for case,lane,mp in sources:
        m=read_json(mp);ids=[];charges={};xps={}
        for metal in ('La','Ca'):
            o=m['outputs'];nested=o.get(metal)
            ir=nested.get('input',nested.get('orca_input')) if nested else o[metal+'_input']
            xr=nested['xyz'] if nested else o[metal+'_xyz']
            ip=verify(ir);xp=verify(xr);text=ip.read_text()
            simple=[line.strip().lower() for line in text.splitlines() if line.strip().startswith('!')]
            if len(simple)!=1 or set(simple[0].split()[1:])!={'r2scan-3c','noautostart','cpcm(water)','defgrid3'}: raise InvalidArtifact('unexpected Hamiltonian')
            if re.search(r'%basis|%pointcharges|%pal|NumGrad',text,re.I): raise InvalidArtifact('unexpected overrides')
            match=re.search(r'^\s*\*\s+xyzfile\s+(-?\d+)\s+1\s+(\S+)\s*$',text,re.M|re.I)
            if not match or Path(match[2]).name!=xp.name: raise InvalidArtifact('coordinate binding mismatch')
            charges[metal]=int(match[1]);d=output/case/lane/metal;d.mkdir(parents=True)
            dest=d/ip.name;xd=d/xp.name;shutil.copyfile(ip,dest);shutil.copyfile(xp,xd);xps[metal]=xd
            tid=f'{case}_{lane}_{metal}';ids.append(tid)
            tasks.append({'task_id':tid,'case':case,'lane':lane,'metal':metal,'charge':charges[metal],
                          'input':record(dest),'xyz':record(xd),'output_path':str(d/(ip.stem+'.out')),
                          'source_input':record(ip),'source_xyz':record(xp),'source_manifest':record(mp),'protocol_id':m['protocol_id']})
        descriptors={'core_charge_La':charges['La'],'core_charge_Ca':charges['Ca'],'coordination':m['coordination']}
        donors=m['coordination'].get('typed_donors',m['coordination'].get('protonated_typed_direct_donors_ordered_by_distance',[]))
        descriptors['donor_composition']=dict(Counter(x['donor_type'] for x in donors))
        descriptors['residue_denticity']=dict(Counter(f"{x['chain']}:{x['resname']}:{x['resnum']}" for x in donors))
        descriptors['mean_donor_distance_A']=sum(x['distance_A'] for x in donors)/len(donors) if donors else None
        group='ggr' if case.startswith('ggr') else 'aequorin' if case.startswith('aequorin') else 'parvalbumin' if case.startswith('carp') else 'pqq_1kb0'
        stratum={'ggr':'direct_same_assay_direction','aequorin':'protein_level_apparent_affinity','parvalbumin':'supporting_cross_study','pqq_1kb0':'canonical_PQQ_class_transfer'}[group]
        row={'case':case,'lane':lane,'protocol_id':m['protocol_id'],'source_manifest':record(mp),'tasks':ids,
             'paired_invariants':paired(xps['La'],xps['Ca'],charges['La'],charges['Ca']),
             'descriptors':descriptors,'biological_group':group,'evidence_stratum':stratum,
             'evaluation_role':'frozen_original_evaluation' if lane!='repaired' else 'chemistry_repair_development'}
        row['cache_key']=cache_key(row);cases.append(row)
    release=root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json'
    result={'schema_version':'alquemia.baseline_benchmark.v1','status':'prepared','tasks':tasks,'cases':cases,
            'reference':ref,'release':record(release),
            'agreement':record(root/'diagnostics/baseline_benchmark_20260915/AGREEMENT.md'),
            'execution_policy':{n:record(root/'scripts'/f) for n,f in [('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py')]},
            'implementation':record(__file__),
            'orca':record('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca')}
    write_new(output/'manifest.json',result)
    return result


def collect(manifest):
    from run_orca_task_manifest import load_manifest_tasks, _completed_attempt_is_valid
    from affordable_common import digest
    m,loaded=load_manifest_tasks(manifest);verify(m['reference']['manifest']);verify(m['release'])
    for artifacts in m['reference']['artifacts'].values():
        for artifact in artifacts.values(): verify(artifact)
    lookup={t['task_id']:t for t in loaded}; endpoints={};rows=[]
    for t in m['tasks']:
        op=Path(t['output_path']);rp=Path(str(op)+'.execution.json');r={'status':'not_run','energy_hartree':None}
        if op.exists():
            try:
                if not rp.exists(): raise InvalidArtifact('running_or_missing_receipt')
                if not _completed_attempt_is_valid(rp,op,manifest_sha256=digest(manifest),task=lookup[t['task_id']],runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer']): raise InvalidArtifact('receipt validation failed')
                receipt=read_json(rp)
                if receipt['orca_executable']!=m['orca'] or receipt['task_id']!=t['task_id']: raise InvalidArtifact('executable or task mismatch')
                for artifact in receipt['artifacts'].values(): verify(artifact)
                for key in ('source_input','source_xyz','source_manifest'): verify(t[key])
                seconds=(datetime.fromisoformat(receipt['finished_at_utc'])-datetime.fromisoformat(receipt['started_at_utc'])).total_seconds()
                r={'status':'complete','energy_hartree':energy(op),'output':record(op),'receipt':record(rp),
                   'wall_seconds':seconds,'mpi_ranks':receipt['parallelism']['nprocs'],
                   'assigned_rank_seconds':seconds*receipt['parallelism']['nprocs']}
            except (InvalidArtifact,ValueError) as exc: r.update(status='unavailable',reason=str(exc))
        endpoints[t['task_id']]=r
    for c in m['cases']:
        row=dict(c);ep={metal:endpoints[tid] for metal,tid in zip(('La','Ca'),c['tasks'])};row['endpoints']=ep
        row.update(status='unavailable',score=None,decision='unavailable',environment_status='excluded_failed_validation')
        if all(e['status']=='complete' for e in ep.values()):
            score=contrast(ep['Ca']['energy_hartree'],ep['La']['energy_hartree'],m['reference']['delta_E_aquo_hartree'])
            row.update(status='complete',score=score,decision=classify_raw(score['R_kcal_mol'],read_json(verify(m['release'])),c['protocol_id']))
            if c['case'].startswith('ggr') and c['lane']=='baseline': row['frozen_direction_test']='pass' if score['S_kcal_mol']<0 else 'fail'
            if c['case']=='pqq_1kb0': row['frozen_class_transfer_test']='pass' if row['decision']=='Ca-supported' else 'fail'
        rows.append(row)
    differences=[]
    for case in sorted({r['case'] for r in rows if r['lane']=='repaired'}):
        pair={r['lane']:r for r in rows if r['case']==case}
        delta=None;cost_ratio=None
        if all(r['status']=='complete' for r in pair.values()):
            delta=pair['repaired']['score']['R_kcal_mol']-pair['baseline']['score']['R_kcal_mol']
            assigned={lane:sum(e['assigned_rank_seconds'] for e in row['endpoints'].values()) for lane,row in pair.items()}
            cost_ratio=assigned['repaired']/assigned['baseline']
        differences.append({'case':case,'repair_delta_R_kcal_mol':delta,'repair_to_original_assigned_rank_seconds_ratio':cost_ratio})
    vectors=[]
    for group in ('aequorin','parvalbumin'):
        for lane in ('baseline','repaired'):
            members=[row for row in rows if row['biological_group']==group and row['lane']==lane]
            values=[None if row['score'] is None else row['score']['S_kcal_mol'] for row in members]
            summary={'biological_group':group,'lane':lane,'ordered_sites':[row['case'] for row in members],
                     'S_vector_kcal_mol':values,'independent_site_labels':False}
            if group=='aequorin' and all(value is not None for value in values):
                summary.update(median=statistics.median(values),minimum=min(values),maximum=max(values),
                               negative_count=sum(value<0 for value in values),positive_count=sum(value>0 for value in values))
            vectors.append(summary)
    events=manifest.parent/'budget_events.jsonl'
    cost_events=[json.loads(line) for line in events.read_text().splitlines()] if events.exists() else []
    costs=[event for event in cost_events if 'allocated_core_seconds' in event]
    return {'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
            'collector':record(__file__),'execution_cost_records':costs,'ordered_vectors':vectors,
            'manifest':record(manifest),'rows':rows,'repair_differences':differences,'completed_endpoints':sum(e['status']=='complete' for e in endpoints.values()),'total_endpoints':len(endpoints)}


def report(collection):
    r=read_json(collection);verify(r['manifest'])
    lines=['# Baseline benchmarking results','',
           f"Status: {r['status']}; {r['completed_endpoints']}/{r['total_endpoints']} endpoints complete.",'',
           'S uses the verified symmetric CN8 ranking gauge. Larger values are more La-like on each protocol scale. Generic and repaired protocols have no calibrated classification bands. Only exact fixed-core PQQ inherits its released bands.','',
           '| Case | Representation | S (kcal/mol) | Decision / frozen test |',
           '|---|---|---:|---|']
    for row in r['rows']:
        score='unavailable' if row['score'] is None else f"{row['score']['S_kcal_mol']:.6f}"
        decision=row.get('frozen_direction_test',row['decision'])
        if 'frozen_class_transfer_test' in row: decision+='; '+row['frozen_class_transfer_test']
        lines.append(f"| {row['case']} | {row['lane']} | {score} | {decision} |")
    lines+=['','## Chemistry repair effect','','Differences below are repaired minus original R; the aquo offset cancels.','']
    for d in r['repair_differences']:
        value=d['repair_delta_R_kcal_mol']
        lines.append(f"- {d['case']}: {'unavailable' if value is None else f'{value:.6f} kcal/mol'}")
    lines+=['','## Interpretation limits','',
            'GGR supplies one direct same-assay Ca-favoring control. Its frozen sign criterion applies to the original v2 record; repaired v3 is a development comparison. Aequorin EF1/EF3/EF4 is an ordered vector against protein-level evidence, not three labeled independent sites. Parvalbumin CD/EF is supporting cross-study evidence. 1KB0 tests PQQ class transfer, not measured relative affinity. Structural sites are grouped by protein. No flexible classifier or threshold was fitted. All unavailable endpoints remain in the denominator.',
            '', 'Typed donor composition, denticity, core charges and mean donor distance are in the JSON ledger. These are descriptive comparators without a fitted decision rule. Canonical calibration already separates by composition; this panel cannot alone prove incremental DFT value.',
            '', 'Baseline defaults and historical results are unchanged. The APBS and global challengers are excluded. No environmental or response correction is supplied.',
            '', '## Measured execution cost','']
    for c in r['execution_cost_records']:
        lines.append(f"- Job {c['slurm_job_id']}: {c['elapsed_seconds']:.3f} workflow seconds, {c['allocated_cpus']} allocated CPUs, {c['allocated_core_seconds']:.3f} allocated core-seconds. Scheduler totals and peak memory are recorded separately after termination.")
    if not r['execution_cost_records']:lines.append('Execution cost is not final yet.')
    lines+=['',f"Unrounded energies and pinned receipts: `{Path(collection).resolve()}`."]
    return '\n'.join(lines)+'\n'


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('prepare');q.add_argument('--root',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    q=s.add_parser('report');q.add_argument('--collection',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.root,a.output)
    elif a.op=='collect':r=collect(a.manifest);write_new(a.output,r)
    else:
        with a.output.open('x') as f:f.write(report(a.collection))
        return
    print(r['status'])
if __name__=='__main__':main()
