#!/usr/bin/env python3
"""Check saved metadata only: no molecular imports, geometry or scientific math."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED = {
    'Hans8DQ2': '85770b1d42be2ebe92051d6596206ad25a09011ca5e9c5a40b5a3e0ed5f2cec0',
    'Hans_8FNR': '8aa2b5a94c56a3f4052ce172d3f990272f2b8c6ba8c98996a1418c69da90fbe2',
    'Mex_8FNS': '4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930',
}


def pin(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return {'path': str(path), 'sha256': digest.hexdigest()}


def read(path):
    return json.loads(Path(path).read_text())


def pins(value):
    # Walk the declared inventory only; never recursively expand linked files.
    if isinstance(value, dict):
        if 'path' in value and 'sha256' in value:
            yield {'path': value['path'], 'sha256': value['sha256']}
        for child in value.values():
            yield from pins(child)
    elif isinstance(value, list):
        for child in value:
            yield from pins(child)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checks = []
    result = {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'failed',
              'scoring_ready': False, 'molecular_calls': 0, 'scientific_allocations': 0,
              'method': 'standard-library file hashes and saved metadata equality only',
              'checks': checks, 'unique_declared_file_pins': 0}

    def check(name, passed):
        checks.append({'name': name, 'passed': bool(passed)})
        if not passed:
            raise ValueError(name)

    try:
        result.update(inventory=pin(args.inventory), contract=pin(args.contract),
                      checker=pin(__file__))
        inv, contract = read(args.inventory), read(args.contract)
        check('inventory_protocol', inv['protocol'] == 'natural_transfer_source_contract_v1')
        check('contract_inventory_binding', contract['inventory'] == result['inventory'])
        declared = {}
        for item in list(pins(inv)) + list(pins(contract)):
            check('consistent_pin:' + item['path'],
                  declared.get(item['path'], item['sha256']) == item['sha256'])
            declared[item['path']] = item['sha256']
        for path, sha in declared.items():
            check('file_hash:' + path, pin(path)['sha256'] == sha)
        result['unique_declared_file_pins'] = len(declared)
        check('all_three_ordered_sources', [s['source_id'] for s in inv['sources']] == list(EXPECTED)
              and contract['ordered_sources'] == list(EXPECTED))
        check('source_only_scope', inv['scope']['source_only'] is True
              and inv['scope']['scoring_ready'] is False
              and inv['scope']['scientific_numerical_analysis'] is False
              and inv['scope']['molecular_calls'] == inv['scope']['scientific_allocations'] == 0)
        for source in inv['sources']:
            sid = source['source_id']
            inputs = read(source['inputs']['path'])
            rec = source['recorded_inventory']
            local = source['local_model']
            whole = source['whole_electronic_state']
            check(sid + ':input_identity', source['inputs']['sha256'] == EXPECTED[sid]
                  and inputs['source_id'] == sid)
            check(sid + ':configuration_namespace', source['configuration']['selected'] == 'A'
                  and source['configuration']['labels_are_biological_chains'] is False
                  and source['prepared_chain'] == 'A'
                  and source['prepared_assembly'] == inputs['assembly'])
            check(sid + ':local_endpoints', source['local_electronic_states'] == inputs['configurations']['A']['endpoints']
                  and local['configuration_A'] == inputs['configurations']['A']
                  and local['configuration_B_historical'] == inputs['configurations']['B'])
            check(sid + ':charges', rec['charge_accounting'] == inputs['charge_accounting']
                  and whole['charge_e'] == inputs['charge_accounting']['total_e'])
            check(sid + ':no_aggregate_state_default', all(whole[k] is None for k in
                  ('multiplicity', 'spin_2S', 'physical_state_qualified'))
                  and all(v is None for k, v in source['unresolved'].items() if k != 'reason'))
            check(sid + ':water_and_spectators', rec['water_count'] == inputs['waters']['total_source']
                  and rec['QM_water_residues'] == inputs['waters']['QM_selected_residues']
                  and rec['spectators'] == inputs['spectator_occupancy'])
            for key in ('physical_sources', 'export', 'core_mapping', 'boundary_mapping'):
                check(sid + ':' + key, source['artifacts'][key] == inputs[key])
            core = read(inputs['core_mapping']['A']['path'])
            check(sid + ':target_mapping', rec['target']['mapping_row'] == core[0]
                  and core[0]['id'] == 'target_EF3' and core[0]['qm_index'] == 0
                  and core[0]['source_site']['chain'] == 'A' and core[0]['source_site']['resid'] == 203
                  and rec['target']['substitution_oxidation_state'] == 3)
            assembly = source['deposited_assembly']
            excerpt = assembly['literal_evidence']
            lines = Path(assembly['source']['path']).read_text().splitlines(keepends=True)
            check(sid + ':literal_assembly_evidence',
                  ''.join(lines[excerpt['first_line'] - 1:excerpt['last_line']]).rstrip('\n')
                  == excerpt['text'].rstrip('\n') and assembly['prepared_partner'] is None)
        expected_pairs = [('Hans8DQ2', 'Mex_8FNS'), ('Hans_8FNR', 'Mex_8FNS')]
        check('two_separate_balanced_contrasts',
              [(c['left_source'], c['right_source']) for c in contract['cross_source_contrasts']] == expected_pairs
              and [(c['Hans_source'], c['Mex_source']) for c in inv['exchange_definitions']['balanced_contrasts']] == expected_pairs)
        check('no_scoring_release', contract['scoring_ready'] is False and all(contract[k] is None for k in
              ('whole_system_energy_recipe', 'whole_system_electronic_states', 'partner_cycle_preparation'))
              and inv['historical_Mex_aggregate_hypotheses']['selected_hypothesis'] is None)
        result['status'] = 'metadata_fidelity_passed'
    except Exception as exc:
        result['error'] = f'{type(exc).__name__}: {exc}'
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ('status', 'scoring_ready', 'unique_declared_file_pins')}))
    return 0 if result['status'] == 'metadata_fidelity_passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
