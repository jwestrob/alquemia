#!/usr/bin/env python3
"""Metadata/hash audit only: no molecular imports or numerical analysis."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
verified = {}


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(value):
    if isinstance(value, dict):
        if 'path' in value and 'sha256' in value:
            path = value['path']
            if path not in verified:
                verified[path] = pin(path)
            assert verified[path]['sha256'] == value['sha256'], path
        else:
            for child in value.values():
                verify(child)
    elif isinstance(value, list):
        for child in value:
            verify(child)


def read(path):
    return json.loads(Path(path).read_text())


def main():
    sources = {
        'inventory': pin(ROOT / 'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json'),
        'MACE_collection': pin(ROOT / 'workspaces/exaa_metal_position_20261014/scoring_v1/execution/FINAL_COLLECTION.json'),
        'component_collection': pin(ROOT / 'workspaces/ca_coverage_recovery_20261012/two_ca_solvent_q9_v1/FINAL_COLLECTION.json'),
        'SCC_manifest': pin(ROOT / 'workspaces/ca_coverage_recovery_20261012/two_ca_solvent_q9_v1/manifest.json'),
        'runtime': pin(ROOT / 'workspaces/orca_parameter_audit_20261012/context_gfn2_v1/RUNTIME.json'),
        'native_manifest': pin(ROOT / 'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/manifest.json'),
        'native_build': pin(ROOT / 'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/canonical_v1/BUILD.json'),
    }
    inventory = read(sources['inventory']['path'])
    scores = {r['source_case_id']: r for r in read(sources['MACE_collection']['path'])['rows']}
    components = {r['source_case_id']: r for r in read(sources['component_collection']['path'])['rows']}
    native_manifest = read(sources['native_manifest']['path'])
    runtime = read(sources['runtime']['path'])
    build = read(sources['native_build']['path'])
    verify(sources)
    # Verify runtime/build pins as data. No executable is loaded or invoked.
    verify(runtime)
    verify(native_manifest['implementation'])
    verify(native_manifest['provenance'])
    verify(native_manifest['build'])
    verify(build)
    rows = []
    for source in inventory['rows']:
        source_id = source['case_id']
        for metal in ('Ca', 'La'):
            for candidate in ('origin', 'adaptive_La'):
                cell = source['matrix'][metal][candidate]
                original = scores[source_id]['matrix'][metal][candidate]['original']
                assert original['status'] == 'reused'
                for name in ('xyz', 'MACE', 'MACE_request'):
                    assert original[name] == cell[name], (source_id, metal, candidate, name)
                    verify(cell[name])
                assert original['energy_eV'] == cell['MACE_eV']
                native = read(cell['native_frozen_ALPB_result']['path'])
                inputs = native['input']
                assert native['status'] == 'available'
                assert native['manifest'] == sources['native_manifest']
                assert inputs['xyz'] == cell['xyz']
                assert inputs['state_signature'] == source['state_signature']
                assert inputs['maps'] == source['maps']
                assert inputs['source_preparation'] == source['source_preparation']
                assert inputs['atom_count'] == cell['atom_count'] == 192
                assert inputs['charge'] == cell['charge'] == original['charge']
                assert inputs['multiplicity'] == cell['multiplicity'] == original['multiplicity'] == 1
                assert inputs['valence_electrons'] == cell['GFN2_valence_electrons'] == 588
                comp = components[source_id]['components'][metal + '__' + candidate]
                assert comp['accepted'] is True
                assert len(native['checks']) == len(native['vacuum_precision']) == 4
                assert all(c['passed'] and c['informative'] for c in native['checks'])
                assert all(c['informative'] for c in native['vacuum_precision'])
                endpoints = {}
                for medium in ('vacuum', 'water'):
                    endpoint = inputs['endpoints'][medium]
                    point = read(endpoint['point']['path'])
                    assert endpoint['status'] == 'available'
                    assert endpoint['point'] == comp['inputs'][medium]
                    assert endpoint['point_manifest'] == sources['SCC_manifest']
                    assert endpoint['executed_xyz']['sha256'] == cell['xyz']['sha256']
                    assert native[medium]['xyz']['sha256'] == cell['xyz']['sha256']
                    assert native[medium]['charges']['sha256'] == endpoint['charges']['sha256']
                    assert all(point['backends'][b]['status'] == 'available' for b in ('xtb', 'tblite'))
                    assert point['backends']['xtb']['solver_evidence']['convergence_declared'] is True
                    assert point['backends']['tblite']['solver_evidence']['library_returned_without_error'] is True
                    assert len(point['backends']['xtb']['charges_e']) == 192
                    assert len(point['backends']['tblite']['charges_e']) == 192
                    verify(endpoint)
                    verify(native[medium])
                    execution = read(endpoint['execution']['path'])
                    assert execution['returncode'] == 0
                    endpoints[medium] = {
                        **endpoint,
                        'original_absolute_gate_accepted': point['accepted'],
                        'original_absolute_gate_energy_hartree': point['energy_hartree'],
                        'raw_SCC_returncode': execution['returncode'],
                        'raw_SCC_command': execution['command'],
                        'backend_status': {b: point['backends'][b]['status'] for b in ('xtb', 'tblite')},
                        'charge_sum_residual_e': {b: point['backends'][b]['charge_sum_residual_e'] for b in ('xtb', 'tblite')},
                    }
                for name in ('cavity', 'born_matrix_with_HB'):
                    assert native['water'][name]['sha256'] == native['vacuum'][name]['sha256']
                verify(comp['verified_output_evidence'])
                verify(cell['native_frozen_ALPB_result'])
                verify(source['state_signature'])
                verify(source['maps'])
                verify(source['source_preparation'])
                rows.append({
                    'source_case_id': source_id, 'role': source['role'],
                    'metal': metal, 'candidate': candidate, 'status': 'available', 'missing': [],
                    'xyz': cell['xyz'], 'MACE': cell['MACE'], 'MACE_request': cell['MACE_request'],
                    'atom_count': cell['atom_count'], 'charge': cell['charge'], 'multiplicity': cell['multiplicity'],
                    'valence_electrons': cell['GFN2_valence_electrons'],
                    'state_signature': source['state_signature'], 'maps': source['maps'],
                    'source_preparation': source['source_preparation'],
                    'native_result': cell['native_frozen_ALPB_result'],
                    'G_qvac_hartree': native['G_qvac_hartree'],
                    'native_components_hartree': native['vacuum']['components_Eh'],
                    'native_diagonal_checks': native['checks'], 'vacuum_precision': native['vacuum_precision'],
                    'native_replay': {m: native[m] for m in ('water', 'vacuum')},
                    'SCC_endpoints': endpoints,
                    'component_v2_acceptance': {k: v for k, v in comp.items() if k != 'verified_output_evidence'},
                    'same_cavity_and_Born_matrix_bytes': True,
                    'new_transplanted_xyz': scores[source_id]['matrix'][metal][candidate]['transplanted']['xyz'],
                    'original_charge_reuse_valid_for_transplanted_geometry': False,
                })
    assert len(rows) == 8
    pointers = {
        'joint_search': pin(ROOT / 'scripts/adaptive_metal_proposals.py'),
        'frozen_joint_search': pin(ROOT / 'workspaces/motion_envelope_20260923/pilot34_searches_v1/implementation/adaptive_metal_proposals.py'),
        'kinematics': pin(ROOT / 'scripts/mace_site_kinematics.py'),
        'historical_result': pin(ROOT / 'diagnostics/adaptive_metal_20260922/FINAL_RESULT.md'),
        'historical_manifest': pin(ROOT / 'workspaces/adaptive_metal_20260922/common_pool_v1/manifest.json'),
    }
    result = {
        'status': 'eight_original_cells_exactly_available', 'audit_kind': 'metadata_and_byte_hashes_only',
        'declared_original_cells': 8, 'available_original_cells': 8, 'missing_original_cells': [],
        'sources': sources, 'rows': rows,
        'runtime': {k: runtime[k] for k in ('xtb', 'xtb_parameter_file', 'tblite_python', 'tblite_launcher_path', 'driver')},
        'SCC_settings': read(sources['SCC_manifest']['path'])['settings'],
        'native_model': native_manifest['native_model'], 'native_tolerances': native_manifest['tolerances'],
        'native_build': build['binary'], 'qualified_engine_build': build['qualified_engine_build'],
        'static_future_reuse_pointers': pointers,
        'preserved_gate_distinction': 'Original per-point absolute parity remains failed; accepted component-v2 pairs and native replay checks are separately named evidence.',
        'transplanted_scope': 'Original charges/results are controls only; no SCC or solvent evaluation of transplanted coordinates is performed by this audit.',
        'scientific_calls': 0, 'new_geometry': False, 'verified_artifacts': list(verified.values()),
    }
    (OUT / 'EVIDENCE.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'evidence': pin(OUT / 'EVIDENCE.json'), 'verified_artifacts': len(verified)}))


if __name__ == '__main__':
    main()
