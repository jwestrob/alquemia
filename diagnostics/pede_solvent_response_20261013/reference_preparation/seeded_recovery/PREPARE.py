#!/usr/bin/env python3
"""Metadata/byte-only preparation of one conditional same-state MORead recovery."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[4]
PRIMARY=ROOT/'workspaces/pede_solvent_response_20261013/reference_sample1_Ca_v1'
DIAG=Path(__file__).parent

def pin(p):
    p=Path(p).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}

def main():
    p=argparse.ArgumentParser();p.add_argument('workspace',type=Path);a=p.parse_args()
    m=json.loads((PRIMARY/'manifest.json').read_text());pm=pin(PRIMARY/'manifest.json')
    source_collection=PRIMARY/'continuation_v2/sample1_Ca_water_COLLECTION.json'
    c=json.loads(source_collection.read_text());row=next(r for r in c['cells'] if r['task_id']=='sample1_Ca_water')
    assert row['status']=='complete_diagnostic' and all(row['checks'].values())
    gbw,=[r for r in row['retained_artifacts'] if Path(r['path']).name=='endpoint.runtime.gbw']
    assert pin(gbw['path'])==gbw
    target=next(t for t in m['tasks'] if t['task_id']=='sample1_Ca_vacuum')
    source=next(t for t in m['tasks'] if t['task_id']=='sample1_Ca_water')
    assert target['xyz']['sha256']==source['xyz']['sha256'] and target['charge']==source['charge']==-3 and target['multiplicity']==source['multiplicity']==1
    w=a.workspace.resolve();w.mkdir(parents=True,exist_ok=False);impl=w/'implementation';impl.mkdir()
    implementation={}
    for name,r in m['implementation'].items():
        if name in ('PREPARE.py','RUN.sbatch'):continue
        src=Path(r['snapshot']['path']);assert pin(src)==r['snapshot'];shutil.copyfile(src,impl/name)
        implementation[name]={'source':pin(src),'snapshot':pin(impl/name)}
    for name,src in [('SEEDED.py',DIAG/'SEEDED.py'),('PREPARE.py',Path(__file__)),('CHECK.sbatch',DIAG/'CHECK.sbatch'),('EXECUTE.sbatch',DIAG/'EXECUTE.sbatch'),
                     ('COLLECT_V2.py',DIAG.parent/'recovery/COLLECT_V2.py')]:
        shutil.copyfile(src,impl/name);implementation[name]={'source':pin(src),'snapshot':pin(impl/name)}
    td=w/'tasks/sample1_Ca_vacuum_seeded';td.mkdir(parents=True)
    shutil.copyfile(target['xyz']['path'],td/'core.xyz');shutil.copyfile(gbw['path'],td/'initial.gbw')
    assert pin(td/'core.xyz')['sha256']==target['xyz']['sha256'] and pin(td/'initial.gbw')['sha256']==gbw['sha256']
    text='! r2SCAN-3c NoAutostart DefGrid3 TightSCF\n%scf\n Guess MORead\n MOInp "initial.gbw"\n MaxIter 125\nend\n* xyzfile -3 1 core.xyz\n'
    (td/'endpoint.inp').write_text(text)
    t=dict(target,task_id='sample1_Ca_vacuum_seeded',input=pin(td/'endpoint.inp'),xyz=pin(td/'core.xyz'),output_path=str(td/'endpoint.out'),initial_gbw=pin(td/'initial.gbw'))
    m.update(schema='pede_conditional_same_state_MORead.v1',protocol_id='PedE_sample1_Ca_vacuum_same_state_water_MORead_recovery_v1',tasks=[t],implementation=implementation,
             primary_manifest=pm,primary_task_id=target['task_id'],primary_task=target,source_task=source,
             seed_source={'collection':pin(source_collection),'execution':row['execution'],'output':row['output'],'gbw':gbw},
             actual_parser_recovery=pin(ROOT/'workspaces/pede_solvent_response_20261013/reference_collection_recovery_v2/execution/EXECUTION.json'),
             moread_fixture=pin(ROOT/'workspaces/metal_environment_response_20260926/lady_repaired_la_seeded_v1/La_A/endpoint.out'),
             root_plan=pin(ROOT/'diagnostics/pede_solvent_response_20261013/SEEDED_RECOVERY_PLAN.md'),
             initialization_only=True,automatic_fallback=False,primary_denominator_unchanged=True)
    m['execution_policy']['task_runner']=pin(impl/'run_orca_task_manifest.py');m['execution_policy']['runtime_renderer']=pin(impl/'render_orca_runtime_input.py')
    (w/'manifest.json').write_text(json.dumps(m,indent=2)+'\n');print(json.dumps(pin(w/'manifest.json')))

if __name__=='__main__':main()
