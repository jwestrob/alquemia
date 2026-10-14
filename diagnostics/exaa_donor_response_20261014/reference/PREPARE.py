#!/usr/bin/env python3
"""Bind saved coordinates and executable bytes; numerical preflight runs in Slurm."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[3]
OLD = ROOT/'workspaces/pede_solvent_response_20261013/reference_v1/manifest.json'
PROTOCOL = 'ExaA_saved_common_La_path_native_r2SCAN3c_vacuum_v1'

def pin(path):
    path=Path(path).resolve(); h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return {'path':str(path),'sha256':h.hexdigest()}

def checked(r):
    if pin(r['path'])['sha256']!=r['sha256']: raise ValueError('Changed input: '+r['path'])
    return Path(r['path'])

def write(p,data):
    with p.open('x') as f: json.dump(data,f,indent=2); f.write('\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('inventory',type=Path);p.add_argument('inventory_sha');p.add_argument('workspace',type=Path)
    a=p.parse_args(); invpin={'path':str(a.inventory.resolve()),'sha256':a.inventory_sha}
    inv=json.loads(checked(invpin).read_text());old=json.loads(OLD.read_text())
    w=a.workspace.resolve();w.mkdir(parents=True,exist_ok=False);impl=w/'implementation';impl.mkdir()
    copies={name:checked(old['implementation'][name]['snapshot']) for name in (
        'run_orca_task_manifest.py','render_orca_runtime_input.py','_base_render_orca_runtime_input.py','result_protocol.py')}
    copies.update({'OLD_REFERENCE.py':checked(old['implementation']['REFERENCE.py']['snapshot']),
        'COLLECT_V2.py':ROOT/'diagnostics/pede_solvent_response_20261013/reference_preparation/recovery/COLLECT_V2.py',
        **{name:Path(__file__).with_name(name) for name in ('PREPARE.py','REFERENCE.py','RUN.sbatch')}})
    implementation={}
    for name,src in copies.items():
        dst=impl/name;shutil.copyfile(src,dst);implementation[name]={'source':pin(src),'snapshot':pin(dst)}
    sources=[];tasks=[]
    for sample in (0,3):
        row,=[r for r in inv['rows'] if r['case_id'].endswith(f'sample-{sample}')]
        source={k:row[k] for k in ('case_id','role','selection_id','prepared_admission','source_preparation','state_signature','maps','native_proposals')}
        source.update(sample=sample,source_case_id=row['case_id'],atom_count=192,cap_count=14,water_inventory=[],
                      endpoints={z:row['matrix'][z]['origin']['xyz'] for z in ('Ca','La')})
        for key in ('prepared_admission','source_preparation','state_signature'): checked(source[key])
        sources.append(source)
        for metal in ('Ca','La'):
            for candidate in ('origin','adaptive_La'):
                cell=row['matrix'][metal][candidate];xyz=checked(cell['xyz'])
                tid=f'sample{sample}_{metal}_{candidate}_vacuum';td=w/'tasks'/tid;td.mkdir(parents=True)
                shutil.copyfile(xyz,td/'core.xyz');charge=cell['charge'];mult=cell['multiplicity']
                (td/'endpoint.inp').write_text(f'! r2SCAN-3c NoAutostart DefGrid3 TightSCF\n* xyzfile {charge} {mult} core.xyz\n')
                tasks.append({'task_id':tid,'sample':sample,'candidate':candidate,'metal':metal,'medium':'vacuum',
                    'charge':charge,'multiplicity':mult,'source_xyz':cell['xyz'],'input':pin(td/'endpoint.inp'),
                    'xyz':pin(td/'core.xyz'),'output_path':str(td/'endpoint.out'),
                    'saved_MACE':cell['MACE'],'saved_MACE_eV':cell['MACE_eV'],
                    'saved_frozen_ALPB':cell['native_frozen_ALPB_result']})
    for r in old['runtime_files']: checked(r)
    fixture=next(t['output_path'] for t in old['tasks'] if t['task_id']=='sample1_La_water')
    m={'schema':'exaa_common_path_dft_reference.v1','protocol_id':PROTOCOL,
       'created_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'source_collection':invpin,
       'sources':sources,'tasks':tasks,'implementation':implementation,'orca':old['orca'],
       'runtime_files':old['runtime_files'],'expected_version':'6.1.1','expected_basis':'def2-mTZVPP',
       'expected_auxiliary_basis':'def2-mTZVPP/J','La_ECP_core_electrons':46,
       'execution_policy':{'task_runner':pin(impl/'run_orca_task_manifest.py'),
            'runtime_renderer':pin(impl/'render_orca_runtime_input.py'),'workers':1,'mpi_ranks':12,
            'memory':'actual shared allocation with25percentreserve; no override'},
       'parser_fixtures':old['parser_fixtures']+[pin(fixture)],'automatic_retry':False,
       'scout_tasks':['sample3_La_origin_vacuum','sample3_La_adaptive_La_vacuum'],
       'scientific_scope':'Same-coordinate vacuum finite-motion diagnostic; no complete pool, classifier or force qualification',
       'preparation_lineage':pin(OLD),'plan':pin(ROOT/'diagnostics/exaa_donor_response_20261014/PLAN.md')}
    write(w/'manifest.json',m);print(json.dumps(pin(w/'manifest.json')))

if __name__=='__main__':main()
