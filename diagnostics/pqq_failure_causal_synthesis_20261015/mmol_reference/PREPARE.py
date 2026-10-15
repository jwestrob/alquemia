#!/usr/bin/env python3
"""Metadata/hash/copy only: bind exact MMOL193 target, native sources and fixtures."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[3]
EXAA=ROOT/'workspaces/exaa_donor_response_20261014/reference_v1/manifest.json'
EXAA_RESULT=ROOT/'workspaces/exaa_donor_response_20261014/full_comparison_v1/execution/RESULT.json'
PROTOCOL='MMOL1770_Ca0_saved_common3_native_r2SCAN3c_vacuum_v1'
ORDER=('origin','adaptive_Ca','adaptive_La')

def pin(path):
    path=Path(path).resolve();h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(path),'sha256':h.hexdigest()}

def checked(r):
    if pin(r['path'])['sha256']!=r['sha256']:raise ValueError('Changed input: '+r['path'])
    return Path(r['path'])

def write(p,data):
    with p.open('x') as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('inventory',type=Path);p.add_argument('inventory_sha');p.add_argument('workspace',type=Path);a=p.parse_args()
    invpin={'path':str(a.inventory.resolve()),'sha256':a.inventory_sha};inv=json.loads(checked(invpin).read_text())
    old=json.loads(EXAA.read_text());fixture=json.loads(EXAA_RESULT.read_text())
    if fixture['manifest']!=pin(EXAA) or fixture['available_cells']!=8:raise ValueError('Completed exact ExaA fixture required')
    row,=[r for r in inv['rows'] if r['role']=='target']
    if row['source_case_id']!='mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-0' or row['sample']!=0:raise ValueError('Wrong target')
    if row['state']['atom_count']!=193 or row['state']['cap_count']!=13:raise ValueError('Wrong declared source inventory')
    w=a.workspace.resolve();w.mkdir(parents=True,exist_ok=False);impl=w/'implementation';impl.mkdir()
    copies={name:checked(old['implementation'][name]['snapshot']) for name in ('run_orca_task_manifest.py','render_orca_runtime_input.py','_base_render_orca_runtime_input.py','result_protocol.py','OLD_REFERENCE.py','COLLECT_V2.py')}
    copies['EXAA_REFERENCE.py']=checked(old['implementation']['REFERENCE.py']['snapshot'])
    copies.update({name:Path(__file__).with_name(name) for name in ('PREPARE.py','REFERENCE.py','COUNTED_COLLECT.py','COMPARE.py','RUN.sbatch','COMPARE.sbatch')})
    implementation={}
    for name,src in copies.items():
        dst=impl/name;shutil.copyfile(src,dst);implementation[name]={'source':pin(src),'snapshot':pin(dst)}
    source={k:row[k] for k in ('case_id','source_case_id','role','sample','selection_id','prepared_admission','source_preparation','state_signature','maps','native_proposals','proposal_manifest','endpoints')}
    source.update(atom_count=193,cap_count=13,water_inventory=[])
    for key in ('prepared_admission','source_preparation','state_signature','proposal_manifest'):checked(source[key])
    tasks=[]
    for metal in ('Ca','La'):
        for candidate in ORDER:
            cell=row['matrix'][metal][candidate];xyz=checked(cell['xyz']);checked(cell['MACE']);checked(cell['MACE_request'])
            tid=f'sample0_{metal}_{candidate}_vacuum';td=w/'tasks'/tid;td.mkdir(parents=True)
            shutil.copyfile(xyz,td/'core.xyz');q=cell['charge'];mult=cell['multiplicity']
            if q!=(-3 if metal=='Ca' else -2) or mult!=1 or cell['atom_count']!=193:raise ValueError('Wrong declared endpoint state')
            (td/'endpoint.inp').write_text(f'! r2SCAN-3c NoAutostart DefGrid3 TightSCF\n* xyzfile {q} {mult} core.xyz\n')
            tasks.append({'task_id':tid,'sample':0,'candidate':candidate,'metal':metal,'medium':'vacuum','atom_count':193,
                'charge':q,'multiplicity':mult,'source_xyz':cell['xyz'],'input':pin(td/'endpoint.inp'),'xyz':pin(td/'core.xyz'),
                'output_path':str(td/'endpoint.out'),'saved_MACE':cell['MACE'],'saved_MACE_eV':cell['MACE_eV'],
                'saved_MACE_request':cell['MACE_request'],'saved_frozen_ALPB':cell['native_result']})
    for r in old['runtime_files']:checked(r)
    for r in fixture['collections']:checked(r)
    m={'schema':'mmol_saved_common3_dft_reference.v1','protocol_id':PROTOCOL,'created_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
       'source_collection':invpin,'sources':[source],'tasks':tasks,'implementation':implementation,'orca':old['orca'],
       'runtime_files':old['runtime_files'],'expected_version':'6.1.1','expected_basis':'def2-mTZVPP',
       'expected_auxiliary_basis':'def2-mTZVPP/J','La_ECP_core_electrons':46,
       'execution_policy':{'task_runner':pin(impl/'run_orca_task_manifest.py'),'runtime_renderer':pin(impl/'render_orca_runtime_input.py'),
           'workers':1,'mpi_ranks':12,'memory':'actual shared allocation with existing25percentreserve; no override'},
       'parser_fixtures':old['parser_fixtures'],'fixture_manifest':pin(EXAA),'collection_fixtures':fixture['collections'],
       'fixture_join':pin(EXAA_RESULT),'automatic_retry':False,'candidate_order':list(ORDER),
       'scout_tasks':['sample0_Ca_origin_vacuum','sample0_La_origin_vacuum'],
       'continuation_policy':'four proposals require explicit reviewed two-origin receipt for this exact manifest',
       'scientific_scope':'target-only same-coordinate vacuum common3 reference; peer source1 inventory-only; no fit/classifier/optimization',
       'preparation_lineage':pin(EXAA),'historical_operational_selection':row['stored_score']['geometry_choices'],
       'historical_score_result':row['score_result'],
       'plan':pin(ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/MMOL_REFERENCE_PLAN.md')}
    write(w/'manifest.json',m);print(json.dumps(pin(w/'manifest.json')))

if __name__=='__main__':main()
