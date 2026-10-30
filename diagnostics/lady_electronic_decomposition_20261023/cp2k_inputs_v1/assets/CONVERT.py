#!/usr/bin/env python3
"""Literal saved-library conversion only. No native imports, jobs, or calculations.

Run from any directory with the exact Python recorded in PROVENANCE.json.
Existing outputs must be absent or byte-identical. The four geometries are read
and pinned, never rewritten. All decimal coefficient tokens remain strings.
"""
from pathlib import Path
from decimal import Decimal
import collections
import hashlib
import json
import re
import sys

HERE = Path(__file__).resolve().parent
DIAG = HERE.parents[1]
ROOT = DIAG.parents[1]
INPUTS = DIAG / 'inputs'
CP2K = ROOT / 'workspaces/lady_electronic_decomposition_20261023/cp2k_build_v1/cp2k-2026.2'
ANGULAR = 'SPDFG'


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verified(record):
    path = Path(record['path'])
    assert pin(path) == record, f'Pin mismatch: {path}'
    return path


def save(name, text):
    path = HERE / name
    data = text.encode()
    if path.exists():
        assert path.read_bytes() == data, f'Refuse to overwrite different output: {path}'
    else:
        with path.open('xb') as stream:
            stream.write(data)
    return pin(path)


def orca_metal(path, metal):
    lines = [x.strip() for x in path.read_text().splitlines() if x.strip()]
    assert lines.pop(0) == f'NewECP {metal}'
    core = int(lines.pop(0).split()[1])
    assert lines.pop(0) == 'lmax g'
    ecp = {}
    for channel in 'spdfg':
        title, n = lines.pop(0).split()
        assert title == channel
        rows = []
        for i in range(int(n)):
            index, exponent, coefficient, power = lines.pop(0).split()
            assert int(index) == i + 1 and power == '2'
            assert Decimal(exponent) > 0
            rows.append([power, exponent, coefficient])
        ecp[channel] = rows
    assert lines.pop(0) == 'end'
    assert lines.pop(0) == f'NewGTO {metal}'
    shells = []
    while lines[0] != 'end':
        channel, n = lines.pop(0).split()
        rows = []
        for i in range(int(n)):
            index, exponent, coefficient = lines.pop(0).split()
            assert int(index) == i + 1 and Decimal(exponent) > 0
            rows.append([exponent, coefficient])
        shells.append({'l': ANGULAR.index(channel), 'rows': rows})
    assert lines == ['end']
    assert core == {'La': 46, 'Dy': 55}[metal]
    # This converter is intentionally restricted to the archived zero-local ECPs.
    assert ecp['g'] == [['2', '1.0000000', '0.0000000']]
    return {'core': core, 'local_original_channel': 'g', 'channels': ecp, 'shells': shells}


def parse_cp2k_basis(text):
    lines = [x.strip() for x in text.splitlines() if x.strip() and not x.startswith('#')]
    out = {}
    while lines:
        element, name = lines.pop(0).split()
        count = int(lines.pop(0))
        shells = []
        for _ in range(count):
            n, lo, hi, primitives, contractions = map(int, lines.pop(0).split())
            assert n == 1 and lo == hi and contractions == 1
            shells.append({'l': lo, 'rows': [lines.pop(0).split() for _ in range(primitives)]})
        out[element] = {'name': name, 'shells': shells}
    return out


def parse_bse_orca(text):
    names = {'HYDROGEN': 'H', 'CARBON': 'C', 'OXYGEN': 'O'}
    lines = [x.strip() for x in text.splitlines() if x.strip() and not x.startswith('!')]
    assert lines.pop(0) == '$DATA' and lines.pop() == '$END'
    out = {}
    while lines:
        element = names[lines.pop(0)]
        shells = []
        while lines and lines[0] not in names:
            channel, count = lines.pop(0).split()
            rows = []
            for i in range(int(count)):
                index, exponent, coefficient = lines.pop(0).split()
                assert int(index) == i + 1
                rows.append([exponent, coefficient])
            shells.append({'l': ANGULAR.index(channel), 'rows': rows})
        out[element] = shells
    return out


def pattern(shells):
    groups = collections.defaultdict(list)
    for s in shells:
        groups[s['l']].append(len(s['rows']))
    return '/'.join(''.join(map(str, groups[l])) for l in sorted(groups))


def main():
    parent = json.loads((INPUTS / 'PROVENANCE.json').read_text())
    downloads = json.loads((HERE / 'BSE_DOWNLOADS.json').read_text())
    for row in downloads['records']:
        assert pin(row['file'])['sha256'] == row['sha256']
    v1 = json.loads((HERE / 'def2-TZVP_HCO_BSE_v1.json').read_text())
    v0 = json.loads((HERE / 'def2-TZVP_HCO_BSE_v0.json').read_text())
    bse_cp2k = parse_cp2k_basis((HERE / 'def2-TZVP_HCO_BSE_v1.cp2k').read_text())
    bse_orca = parse_bse_orca((HERE / 'def2-TZVP_HCO_BSE_v1.orca').read_text())
    basis, metals = {}, {}
    for z, element in [('1', 'H'), ('6', 'C'), ('8', 'O')]:
        assert v1['elements'][z]['electron_shells'] == v0['elements'][z]['electron_shells']
        shells = []
        for shell in v1['elements'][z]['electron_shells']:
            assert len(shell['angular_momentum']) == len(shell['coefficients']) == 1
            shells.append({'l': shell['angular_momentum'][0],
                           'rows': [list(x) for x in zip(shell['exponents'], shell['coefficients'][0], strict=True)]})
        assert shells == bse_cp2k[element]['shells'] == bse_orca[element]
        basis[element] = shells
    for metal in ['La', 'Dy']:
        source = verified(parent['files']['metal_assets'][metal]['basis']['copy'])
        metals[metal] = orca_metal(source, metal)
        basis[metal] = metals[metal]['shells']
    expected_patterns = {'H': '311/1', 'C': '62111/411/11/1', 'O': '62111/411/11/1',
                         'La': '411111/4211/311/11', 'Dy': '411111/4211/311/11'}
    assert {k: pattern(v) for k, v in basis.items()} == expected_patterns

    text = ['# Literal archived La/Dy orbitals; documented BSE def2-TZVP H/C/O.',
            '# No exponent/contraction rounding, refitting, decontraction or rescaling.']
    for element, shells in basis.items():
        text += [f'{element} LADY-ARCHIVED', str(len(shells))]
        for shell in shells:
            l = shell['l']
            text.append(f'1 {l} {l} {len(shell["rows"])} 1')
            text += [' '.join(row) for row in shell['rows']]
        text.append('')
    basis_text = '\n'.join(text) + '\n'
    parsed = parse_cp2k_basis(basis_text)
    assert {k: v['shells'] for k, v in parsed.items()} == basis

    # Retain official all-electron charge-distribution records literally.
    all_text = (CP2K / 'data/ALL_POTENTIALS').read_text()
    all_records = {}
    potential = ['# H/C/O ALL records copied literally from this CP2K release.',
                 '# Metal rows: radial n, Gaussian exponent, coefficient; local g is UL.']
    for element in ['H', 'C', 'O']:
        match = re.search(rf'^{element} ALLELECTRON ALL\n[^\n]+\n[^\n]+\n', all_text, re.M)
        assert match
        all_records[element] = match.group()
        potential.append(match.group().rstrip('\n'))
    for metal, obj in metals.items():
        potential += [f'ECP-LADY-{metal}{obj["core"]}', f'{metal} nelec {obj["core"]}', f'{metal} ul']
        potential += [' '.join(row) for row in obj['channels']['g']]
        for channel in 'spdf':
            potential.append(f'{metal} {channel.upper()}')
            potential += [' '.join(row) for row in obj['channels'][channel]]
        potential.append('END')
    potential_text = '\n'.join(potential) + '\n'
    # Independent token reconstruction of the exact emitted metal sections.
    for metal, obj in metals.items():
        block = potential_text.split(f'ECP-LADY-{metal}{obj["core"]}\n')[1].split('END\n')[0].splitlines()
        assert block.pop(0) == f'{metal} nelec {obj["core"]}'
        channels = {}
        for channel, label in [('g', 'ul'), ('s', 'S'), ('p', 'P'), ('d', 'D'), ('f', 'F')]:
            assert block.pop(0) == f'{metal} {label}'
            channels[channel] = [block.pop(0).split() for _ in obj['channels'][channel]]
        assert not block and channels == obj['channels']

    geometry_checks = {}
    counts = {element: sum(2*s['l']+1 for s in shells) for element, shells in basis.items()}
    for key, system in parent['systems'].items():
        xyz = verified(system['xyz'])
        out = verified(system['original_output'])
        inp = verified(system['archived_native_input_copy'])
        atoms = [row.split()[0] for row in xyz.read_text().splitlines()[2:]]
        assert len(atoms) == system['atoms']
        assert '! PBE0 D4 def2-TZVP def2/J' in inp.read_text()
        native = out.read_text()
        for element in set(atoms):
            assert re.search(rf'Group\s+\d+ Type {element}\s+:.*pattern \{{{re.escape(expected_patterns[element])}\}}', native)
        ao_count = sum(counts[x] for x in atoms)
        assert re.search(rf'Number of basis functions\s+\.\.\.\s+{ao_count}\b', native)
        assert system['effective_multiplicity'] == 1
        geometry_checks[key] = {'geometry': pin(xyz), 'native_input': pin(inp), 'native_output': pin(out),
                               'atom_counts': dict(collections.Counter(atoms)), 'spherical_AOs': ao_count,
                               'archived_patterns_and_AO_count_match': True,
                               'charge': system['charge'], 'effective_multiplicity': 1,
                               'explicit_electrons': system['explicit_electrons']}

    files = {'basis': save('BASIS_LADY', basis_text), 'potentials': save('POTENTIAL_LADY', potential_text)}
    source_files = ['src/atom_types.F', 'src/qs_kind_types.F', 'src/aobasis/basis_set_types.F',
                    'src/libgrpp_integrals.F', 'src/grpp/grpp_potential.c', 'data/ECP_POTENTIALS',
                    'data/ALL_POTENTIALS']
    result = {'status': 'prepared_literal_cp2k_basis_and_ecp_no_native_validation',
              'implementation': pin(__file__), 'python': pin(sys.executable),
              'parent_inputs': pin(INPUTS / 'PROVENANCE.json'), 'files': files,
              'download_authority': pin(HERE / 'BSE_DOWNLOADS.json'),
              'cp2k_reader_authorities': {name: pin(CP2K / name) for name in source_files},
              'metal_source_assets': {m: parent['files']['metal_assets'][m]['basis']['copy'] for m in metals},
              'metal_records': metals, 'basis_records': basis, 'all_electron_records': all_records,
              'geometry_checks': geometry_checks, 'spherical_AOs_by_element': counts,
              'HCO_equivalence': {'status': 'documented_named_library_equivalence',
                'literal_ORCA_runtime_coefficients_available': False,
                'BSE_versions_0_and_1_electron_shells_identical': True,
                'BSE_ORCA_CP2K_exports_tokens_identical': True,
                'archived_ORCA_name_patterns_and_counts_match': True,
                'basis_reference': 'Weigend and Ahlrichs 2005, DOI 10.1039/B508541A',
                'ORCA_official_library_documentation': 'https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html',
                'limit': 'Saved outputs do not print H/C/O primitive coefficients; no byte-level ORCA library extraction claimed.'},
              'ECP_convention': {'row_order': ['radial_n', 'exponent', 'coefficient'],
                'radial_factor': 'r^(n-2) exp(-exponent*r^2)', 'n_is_preserved': 2,
                'local_channel': 'original zero g moved literally to UL; S/P/D/F remain ordered nonlocal projectors',
                'ECP_SEMI_LOCAL_required': True,
                'cross_format_example': 'ORCA manual SD(10,MDF) V example matches CP2K data/ECP_POTENTIALS V rows including n=2'},
              'runtime_requirements': {'basis_name_all_elements': 'LADY-ARCHIVED',
                'potentials': {'H': 'ALL', 'C': 'ALL', 'O': 'ALL', 'La': 'ECP-LADY-La46', 'Dy': 'ECP-LADY-Dy55'},
                'metal_ECP_SEMI_LOCAL': True, 'orbital_normalization': 2,
                'representation': 'spherical Gaussian orbitals; unchanged contractions'},
              'not_converted': 'ORCA def2/J and metal AuxJ are RIJCOSX numerical auxiliaries, not orbital basis; CP2K numerical settings are separately declared.',
              'molecular_calls': 0, 'jobs': 0,
              'limitations': ['Source and formatting checks only; no CP2K input execution, molecular energy or EDA validation.',
                              'Effective Dy singlet is the archived frozen-f model, not physical spin/SOC qualification.']}
    save('PROVENANCE.json', json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'files': files,
                      'provenance': pin(HERE / 'PROVENANCE.json'), 'AO_counts': counts}))


if __name__ == '__main__':
    main()
