#!/usr/bin/env python3
"""Read-only checks on the real acquired/prepared extension; never executes ORCA."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import run_orca_task_manifest as runner


def verify(record):
    p=Path(record['path'])
    assert p.is_file(),str(p)
    assert hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256'],str(p)
    return p


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    start=time.monotonic();j=json.loads(args.manifest.read_text());checks=[]
    for r in j['new_structure_inventory']:verify(r['source'])
    checks.append('19 acquired source hashes match')
    row=next(r for r in j['new_structure_inventory'] if r['pdb_id']=='8GY2')
    prep=json.loads(verify(row['preparation']).read_text());assert prep['status']=='prepared_unscored'
    target=prep['target'];manifest_path=verify(target['manifest']);m,tasks=runner.load_manifest_tasks(manifest_path)
    assert len(tasks)==2
    for task in tasks:
        runner._verify_prepared_task(task)
        assert not task['output'].exists(),str(task['output'])
    checks.append('existing runner loads and verifies two real prepared tasks; no endpoint outputs exist')
    for name in ('normalized','normalization_manifest','protonated','protonation_manifest','heavy_coordinate_check'):
        verify(target[name])
    h=json.loads(verify(target['heavy_coordinate_check']).read_text());assert h['passes']
    assert m['fixed_core']['final_fragment_heavy_coordinate_check']['maximum_displacement_A']==0
    p=json.loads(verify(target['protonation_manifest']).read_text())
    assert p['repaired_missing_atom_count']==p['repaired_missing_terminal_atom_count']==0
    checks.append('source-heavy coordinates unchanged; no missing heavy atoms reconstructed')
    ca=verify(m['outputs']['Ca_xyz']).read_text().splitlines()[2:]
    la=verify(m['outputs']['La_xyz']).read_text().splitlines()[2:]
    assert len(ca)==len(la)==71
    differences=[i for i,(a,b) in enumerate(zip(ca,la)) if a!=b]
    assert len(differences)==1
    i=differences[0];assert ca[i].split()[0]=='Ca' and la[i].split()[0]=='La'
    assert ca[i].split()[1:]==la[i].split()[1:]
    assert m['charge_ledger']['expected_total_charges']=={'Ca':-2,'La':-1}
    assert m['charge_ledger']['multiplicity']==1
    checks.append('71-atom pairs differ only in metal element; exact paired coordinates and charge closure')
    assert not row['eligible_direct_affinity'] and row['eligible_PQQ_association_transfer']
    assert row['nearest_existing_sequence']['identity_relative_longer']<.5
    assert j['ready_matched_DFT_MACE_count']==j['new_direct_affinity_groups']==0
    assert all(not r['eligible_direct_affinity'] for r in j['new_structure_inventory'])
    checks.append('association/affinity evidence kept separate; unsupported matched features remain unavailable')
    assert j['new_energy_evaluations']==j['new_classifier_fits']==j['new_cluster_jobs']==0
    checks.append('no new scoring/fitting/cluster execution represented')
    result={'passed':True,'checks':checks,'wall_seconds':time.monotonic()-start,'integration_scope':'real-source integrity and preparation only; no scientific energy executable run'}
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
