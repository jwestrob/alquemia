#!/usr/bin/env python3
"""Metadata-only freeze of eight byte-identical saved-origin DFT inputs."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[3]
ORCA = Path('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca')
SOURCE = ROOT/'workspaces/ca_coverage_recovery_20261012/two_ca_scoring_c5q88_v1/NATIVE_COLLECTION.json'
PROTOCOL = 'PedE_fixed_origin_native_r2SCAN3c_DefGrid3_TightSCF_CPCM_response_v1'

def pin(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return {'path': str(path), 'sha256': h.hexdigest()}

def checked(record):
    assert pin(record['path'])['sha256'] == record['sha256'], record['path']
    return Path(record['path'])

def main():
    p = argparse.ArgumentParser(); p.add_argument('workspace', type=Path); a = p.parse_args()
    w = a.workspace.resolve(); w.mkdir(parents=True, exist_ok=False)
    impl = w/'implementation'; impl.mkdir()
    copies = {'run_orca_task_manifest.py': ROOT/'scripts/run_orca_task_manifest.py',
              'render_orca_runtime_input.py': ROOT/'scripts/metal_environment_runtime.py',
              '_base_render_orca_runtime_input.py': ROOT/'scripts/render_orca_runtime_input.py',
              'result_protocol.py': ROOT/'scripts/result_protocol.py',
              'PREPARE.py': Path(__file__),
              'REFERENCE.py': Path(__file__).with_name('REFERENCE.py'),
              'RUN.sbatch': Path(__file__).with_name('RUN.sbatch')}
    implementation = {}
    for name, src in copies.items():
        dst = impl/name; shutil.copyfile(src, dst)
        implementation[name] = {'source': pin(src), 'snapshot': pin(dst)}
    runtime_paths = [ORCA]
    for stem in ('orca_startup', 'orca_guess', 'orca_leanscf', 'orca_prop', 'orca_propint'):
        runtime_paths += [ORCA.parent/stem, ORCA.parent/(stem+'_mpi')]
    runtime_paths += [ORCA.parent/'lib'/x for x in ('liborca_tools_6_1_1.so.6', 'libstdc++.so.6.0.32', 'libgcc_s.so.1')]
    runtime = [pin(x) for x in runtime_paths]
    rows = json.loads(SOURCE.read_text())['cases']; sources = []; tasks = []
    for sample in (1, 3):
        sid = f'q88jh5-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-{sample}'
        row, = [r for r in rows if r['source_case_id'] == sid]
        records = {k: row[k] for k in ('prepared_admission', 'source_preparation', 'state_signature')}
        for rec in records.values(): checked(rec)
        admission = json.loads(Path(records['prepared_admission']['path']).read_text())
        state = json.loads(Path(records['state_signature']['path']).read_text())
        assert admission['status'] == 'prepared' and admission['atom_count'] == 192 and admission['cap_count'] == 14
        assert state['water_inventory'] == []
        source = {'sample': sample, 'case_id': row['case_id'], 'source_case_id': sid,
                  'selection_id': row['selection_id'], **records,
                  'atom_count': 192, 'cap_count': 14, 'water_inventory': [], 'endpoints': {}}
        for metal, charge in (('Ca', -3), ('La', -2)):
            origin = row['matrix'][metal]['origin']
            assert origin['status'] == 'complete' and origin['charge'] == charge and origin['multiplicity'] == 1
            xyz = checked(origin['xyz']); source['endpoints'][metal] = origin['xyz']
            for medium in ('water', 'vacuum'):
                tid = f'sample{sample}_{metal}_{medium}'; td = w/'tasks'/tid; td.mkdir(parents=True)
                shutil.copyfile(xyz, td/'core.xyz')
                header = '! r2SCAN-3c NoAutostart DefGrid3 TightSCF' + (' CPCM(Water)' if medium == 'water' else '')
                (td/'endpoint.inp').write_text(header+f'\n* xyzfile {charge} 1 core.xyz\n')
                tasks.append({'task_id': tid, 'sample': sample, 'metal': metal, 'medium': medium,
                              'charge': charge, 'multiplicity': 1, 'source_xyz': origin['xyz'],
                              'input': pin(td/'endpoint.inp'), 'xyz': pin(td/'core.xyz'),
                              'output_path': str(td/'endpoint.out')})
        sources.append(source)
    assert sources[0]['state_signature']['sha256'] == sources[1]['state_signature']['sha256']
    # Bounded hash-only cache lookup: these are known native-DFT manifest families.
    cache_paths = [ROOT/x for x in (
        'workspaces/environment_pqq_dft_20260919/prepared_v2/manifest.json',
        'workspaces/accommodation_torsion_20260920/prepared_v3/dft/manifest.json',
        'workspaces/coordination_preparation_20260919/dft_v1/manifest.json',
        'workspaces/accommodation_fold_DFT_20260920/run_v1/manifest.json')]
    hashes = {t['xyz']['sha256'] for t in tasks}; cache = []
    for path in cache_paths:
        if not path.exists(): cache.append({'path': str(path), 'status': 'absent'}); continue
        raw = path.read_text()
        cache.append({'manifest': pin(path), 'exact_xyz_hash_mentions': sorted(x for x in hashes if x in raw),
                      'scope': 'this manifest only, no recursive archive search'})
    manifest = {'schema': 'pede_fixed_origin_dft_reference.v1', 'protocol_id': PROTOCOL,
                'created_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                'source_collection': pin(SOURCE), 'sources': sources, 'tasks': tasks,
                'implementation': implementation, 'orca': pin(ORCA), 'runtime_files': runtime,
                'expected_version': '6.1.1', 'expected_basis': 'def2-mTZVPP',
                'expected_auxiliary_basis': 'def2-mTZVPP/J', 'La_ECP_core_electrons': 46,
                'execution_policy': {'task_runner': pin(impl/'run_orca_task_manifest.py'),
                                     'runtime_renderer': pin(impl/'render_orca_runtime_input.py'),
                                     'workers': 1, 'mpi_ranks': 12, 'memory': 'actual shared allocation, 25% reserve; no submission override'},
                'scout_order': ['sample1_La_water', 'sample1_La_vacuum'],
                'parser_fixtures': [pin(ROOT/'workspaces/environment_pqq_dft_20260919/prepared_v2/tasks/q88jh5-pqq-la_model__expanded__La/endpoint.out'),
                                    pin(ROOT/'workspaces/metal_environment_response_20260926/scf_check_v1/La_A_strict/endpoint.out')],
                'cache_audit': cache, 'automatic_retry': False,
                'interpretation': 'new fixed-origin diagnostic; no force/density-response or classifier qualification from normal termination alone'}
    (w/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(pin(w/'manifest.json')))

if __name__ == '__main__': main()
