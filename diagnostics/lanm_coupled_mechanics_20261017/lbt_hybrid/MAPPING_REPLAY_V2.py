"""Strict map replay with a bounded allowance only for derived cap Jacobians."""
import copy
import math
import struct
import sys

POLICY = {'revision': 'derived_cap_Jacobian_binary64_replay_v2',
          'max_ULP_distance': 8, 'absolute_maximum': 8*sys.float_info.epsilon,
          'float_format': 'IEEE754_binary64', 'non_Jacobian_leaves': 'exact_type_and_value',
          'zero_sign_and_value': 'exact', 'scientific_tolerances_changed': False}
JACOBIANS = ('jacobian_retained', 'jacobian_omitted')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def bits(value):
    raw = struct.unpack('>Q', struct.pack('>d', value))[0]
    return 0xffffffffffffffff-raw if raw >> 63 else raw+0x8000000000000000


def strict_equal(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(strict_equal(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(strict_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, float):
        return math.isfinite(a) and math.isfinite(b) and a.hex() == b.hex()
    return a == b


def jacobian_paths(mapping):
    require(type(mapping) is dict and set(mapping) == {'core', 'caps', 'field', 'physical_atoms'}, 'Map schema differs')
    require(all(type(mapping[k]) is list for k in mapping), 'Map list type differs')
    require(len(mapping['caps']) == 2, 'Two physical caps required')
    caps = mapping['caps']; core = mapping['core']
    indices = [i for i, r in enumerate(core) if r.get('kind') == 'cap']
    require(indices == [c['qm_index'] for c in caps], 'Cap/core index correspondence differs')
    paths = set()
    for index, cap in enumerate(caps):
        require(type(cap['qm_index']) is int and strict_equal(cap, core[cap['qm_index']]), 'Duplicated cap/core row differs')
        for key in JACOBIANS:
            matrix = cap[key]
            require(type(matrix) is list and len(matrix) == 3
                    and all(type(r) is list and len(r) == 3 for r in matrix), 'Cap Jacobian shape differs')
            require(all(type(v) is float and math.isfinite(v) for r in matrix for v in r), 'Cap Jacobian type or finiteness differs')
            for i in range(3):
                for j in range(3):
                    paths.add(('caps', index, key, i, j))
                    paths.add(('core', cap['qm_index'], key, i, j))
    return paths


def compare(archived, regenerated):
    """No coordinate tolerance; no mutation; return every allowed float difference."""
    require(sys.float_info.radix == 2 and sys.float_info.mant_dig == 53, 'Binary64 Python floats required')
    allowed = jacobian_paths(archived)
    require(jacobian_paths(regenerated) == allowed, 'Allowed cap locations differ')
    differences = []

    def visit(a, b, path):
        label = '/'.join(map(str, path))
        require(type(a) is type(b), 'Map leaf/container type differs: '+label)
        if isinstance(a, dict):
            require(a.keys() == b.keys(), 'Map keys differ: '+label)
            for k in sorted(a):
                visit(a[k], b[k], path+(k,))
        elif isinstance(a, list):
            require(len(a) == len(b), 'Map list length differs: '+label)
            for i, (x, y) in enumerate(zip(a, b)):
                visit(x, y, path+(i,))
        elif isinstance(a, float):
            require(math.isfinite(a) and math.isfinite(b), 'Nonfinite map float: '+label)
            if a.hex() == b.hex():
                return
            require(path in allowed, 'Non-Jacobian map float differs: '+label)
            require(a != 0. and b != 0., 'Zero Jacobian coefficient/sign differs: '+label)
            delta = abs(b-a); ulps = abs(bits(a)-bits(b))
            require(delta <= POLICY['absolute_maximum'] and ulps <= POLICY['max_ULP_distance'],
                    'Derived Jacobian differs beyond bounded roundoff: '+label)
            differences.append({'path': list(path), 'archived': a, 'regenerated': b,
                                'archived_hex': a.hex(), 'regenerated_hex': b.hex(),
                                'absolute_difference': delta, 'ULP_distance': ulps})
        else:
            require(a == b, 'Non-Jacobian map value differs: '+label)

    visit(archived, regenerated, ())
    return {'pass': True, 'policy': POLICY, 'exact_all_leaves': not differences,
            'derived_Jacobian_differences': differences, 'differing_leaves': len(differences),
            'max_absolute_difference': max((r['absolute_difference'] for r in differences), default=0.),
            'max_ULP_distance': max((r['ULP_distance'] for r in differences), default=0),
            'coordinates_identities_charges_types_exact': True}


def negative_controls(archived):
    """Small exact-data mutations; never energy or model calls."""
    require(compare(archived, archived)['exact_all_leaves'], 'Exact archived replay fails')
    ci = archived['caps'][0]['qm_index']
    controls = []

    def reject(name, change):
        altered = copy.deepcopy(archived)
        change(altered)
        try:
            compare(archived, altered)
        except ValueError as exc:
            controls.append({'name': name, 'pass': True, 'expected_rejection': str(exc)})
        else:
            raise ValueError('Negative control admitted: '+name)

    def coordinate(container, index):
        return lambda m: m[container][index]['xyz_A'].__setitem__(0, math.nextafter(m[container][index]['xyz_A'][0], math.inf))

    reject('one_ULP_physical_coordinate_change', coordinate('physical_atoms', 0))
    reject('one_ULP_field_coordinate_change', coordinate('field', 0))
    reject('one_ULP_real_QM_coordinate_change', coordinate('core', 1))
    def cap_coordinate(m):
        for row in (m['caps'][0], m['core'][ci]):
            row['xyz_A'][0] = math.nextafter(row['xyz_A'][0], math.inf)
    reject('one_ULP_cap_coordinate_change', cap_coordinate)
    reject('one_ULP_external_charge_change', lambda m: m['field'][0].__setitem__('charge_e', math.nextafter(m['field'][0]['charge_e'], math.inf)))
    reject('physical_atom_identity_change', lambda m: m['physical_atoms'][0].__setitem__('id', m['physical_atoms'][0]['id']+'_changed'))
    reject('core_atom_identity_change', lambda m: m['core'][1].__setitem__('source_id', m['core'][1]['source_id']+'_changed'))
    def mutate_both(m, f):
        for row in (m['caps'][0], m['core'][ci]):
            row['jacobian_retained'][0][0] = f(row['jacobian_retained'][0][0])
    reject('Jacobian_float_type_change', lambda m: mutate_both(m, str))
    def beyond(v):
        for _ in range(POLICY['max_ULP_distance']+1):
            v = math.nextafter(v, math.inf)
        return v
    reject('Jacobian_nine_ULP_change', lambda m: mutate_both(m, beyond))
    reject('Jacobian_nonroundoff_change', lambda m: mutate_both(m, lambda v: v+1e-8))
    reject('Jacobian_nonfinite_change', lambda m: mutate_both(m, lambda v: math.inf))
    reject('duplicated_cap_mapping_inconsistency', lambda m: m['caps'][0]['jacobian_retained'][0].__setitem__(0, math.nextafter(m['caps'][0]['jacobian_retained'][0][0], math.inf)))
    return {'pass': True, 'controls': controls, 'new_molecular_calls': 0}


def validate_diagnostic(record, read):
    result = read(record)
    require(result['status'] == 'complete' and result['slurm_job_id'] == '12984'
            and result['compared_origins'] == result['declared_origins'] == 4
            and result['new_molecular_calls'] == result['model_loads'] == 0, 'Actual mapping diagnosis unavailable')
    expected = [('7CCO', 'baseline'), ('7CCO', 'expanded'), ('7CCN', 'baseline'), ('7CCN', 'expanded')]
    require([(r['source_id'], r['region']) for r in result['rows']] == expected, 'Actual mapping diagnosis denominator differs')
    checks = []
    for row in result['rows']:
        require(row['status'] == 'compared' and row['repeated_dictionary_equality']
                and row['repeated_typed_leaf_difference_count'] == 0, 'Actual mapping repeat unavailable')
        archived = read(row['archived_mapping'])
        regenerated = [read(p) for p in row['regenerated_mappings']]
        require(len(regenerated) == 2 and strict_equal(*regenerated), 'Actual repeated map differs')
        check = compare(archived, regenerated[0])
        expected_count = 32 if (row['source_id'], row['region']) == ('7CCN', 'baseline') else 0
        require(check['differing_leaves'] == row['typed_leaf_difference_count'] == expected_count
                and check['max_ULP_distance'] == row['max_ulp_distance']
                and check['max_absolute_difference'] == row['max_numeric_absolute_difference'], 'Actual diagnosed residual differs')
        checks.append({'source_id': row['source_id'], 'region': row['region'], 'replay': check,
                       'archived_mapping': row['archived_mapping'], 'regenerated_mapping': row['regenerated_mappings'][0]})
    return {'pass': True, 'diagnostic': record, 'checks': checks,
            'negative_controls': negative_controls(read(result['rows'][2]['archived_mapping']))}
