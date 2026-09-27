"""Prepare an additive QM/MM parameter ledger; never evaluates molecular energies.

This is an explicitly unqualified research Hamiltonian, not a production scorer.
Protein force objects remain pinned in their native OpenMM XML representation.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from affordable_common import read_json, record, verify, write_new

PROTOCOL = 'nikasha_1h4i_additive_mechanics_ledger_v1'


def lj_from_amber(radius, epsilon):
    """Amber Rmin/2 in angstrom, epsilon kcal -> sigma nm, epsilon kJ."""
    return {'sigma_nm': 0.2 * radius / 2**(1/6), 'epsilon_kJ_mol': 4.184 * epsilon}


def metal_parameters(path):
    section = False
    result = {}
    for line in Path(path).read_text().splitlines():
        if line.strip() == 'NONBON':
            section = True
            continue
        fields = line.split()
        if section and fields and fields[0] in ('Ca2+', 'La3+'):
            metal = {'Ca2+': 'Ca', 'La3+': 'La'}[fields[0]]
            radius, epsilon = map(float, fields[1:3])
            result[metal] = dict(lj_from_amber(radius, epsilon), atom_type=fields[0],
                                Rmin_half_A=radius, epsilon_kcal_mol=epsilon,
                                source_line=line, status='unqualified_QMMM_candidate')
    if set(result) != {'Ca', 'La'}:
        raise ValueError('exact Ca2+/La3+ NONBON entries required')
    # Prevent accidental use of a C4-fitted family or an undocumented input.
    if Path(path).name != 'frcmod.ions234lm_126_tip3p':
        raise ValueError('only explicitly declared pure 12-6 TIP3P family supported')
    return result


def graph_distances(adjacency, start, depth=3):
    distances = {start: 0}
    frontier = {start}
    for step in range(1, depth + 1):
        frontier = {j for i in frontier for j in adjacency[i] if j not in distances}
        distances.update({j: step for j in frontier})
    return distances


def build(source_manifest, protein_result, pqq_file, metal_file):
    c = read_json(source_manifest)
    pr = read_json(protein_result)
    pins = pr['artifact_pins']
    atoms = read_json(verify(pins['atoms.json']))
    exceptions = read_json(verify(pins['exceptions.json']))
    bonded = read_json(verify(pins['bonded_supports.json']))
    verify(pins['system.xml'])
    ff = verify(pr['forcefield'])
    coulomb14 = float(ET.parse(ff).find('.//NonbondedForce').attrib['coulomb14scale'])
    core = read_json(verify(c['core_mapping']))
    boundary = read_json(verify(c['boundary_mapping']))
    env = read_json(verify(c['environments']['A']['atoms']))
    pqq = read_json(pqq_file)
    qm = {a['source_index'] for a in core if a['kind'] == 'protein_source'}
    if len(core) != 54 or len(qm) != 23 or len(atoms) != 9113:
        raise ValueError('this finite ledger supports exact current 1H4I 54-QM preparation only')
    if pr['source'] != c['source']:
        raise ValueError('source pin differs')
    physical = []
    qtilde = {a['source_index']: a['charge_e'] for a in env}
    removed = set().union(*(set(l['removed_source_indices']) for l in boundary['ledgers']))
    if set(qtilde) != set(range(len(atoms))) - removed:
        raise ValueError('embedding complement differs')
    qm_by_source = {a['source_index']: a for a in core if a['kind'] == 'protein_source'}
    for i, a in enumerate(atoms):
        if a['source_index'] != i:
            raise ValueError('source indexing changed')
        row = dict(a, physical_index=i, region='QM' if i in qm else 'MM',
                   original_ff_charge_e=a['charge_e'], charge_e=None if i in qm else qtilde.get(i, 0.0),
                   qm_index=qm_by_source[i]['qm_index'] if i in qm else None)
        if i in qm:
            # Use exact archived electronic geometry, not reserialized protein coords.
            row['xyz_A'] = qm_by_source[i]['xyz_A']
        physical.append(row)
    pqq_by_qm = {a['qm_index']: a for a in pqq['atoms']}
    for a in core:
        if a['kind'] in ('cap', 'protein_source'):
            continue
        row = dict(a, physical_index=len(physical), region='QM', charge_e=None)
        if a['kind'] == 'metal':
            row['endpoint_lj'] = metal_parameters(metal_file)
        else:
            if a['qm_index'] not in pqq_by_qm:
                raise ValueError('unsupported real QM atom')
            pp = pqq_by_qm[a['qm_index']]
            if pp['partial_charge_e'] is not None or pp['xyz_A'] != a['xyz_A']:
                raise ValueError('PQQ charge/coordinate artifact differs')
            row.update(lj_from_amber(pp['Rmin_half_A'], pp['epsilon_kcal_mol']))
            row['gaff2_type'] = pp['gaff2_type']
        physical.append(row)
    adjacency = defaultdict(set)
    for t in bonded:
        if t['force'] == 'HarmonicBondForce':
            a,b = t['atoms']; adjacency[a].add(b); adjacency[b].add(a)
    distances = {}
    exrows = []
    for e in exceptions:
        a,b = e['atoms']
        if a not in distances:
            distances[a] = graph_distances(adjacency, a)
        d = distances[a].get(b)
        if d not in (1,2,3):
            raise ValueError('exception without supported 1-2/1-3/1-4 topology')
        scale = coulomb14 if d == 3 else 0.0
        expected = atoms[a]['charge_e'] * atoms[b]['charge_e'] * scale
        if abs(expected - e['chargeprod_e2']) > 1e-12:
            raise ValueError('native exception charge scaling differs')
        region = 'QM' if a in qm and b in qm else ('cross' if a in qm or b in qm else 'MM')
        chargeprod = None if region != 'MM' else physical[a]['charge_e'] * physical[b]['charge_e'] * scale
        exrows.append(dict(e, region=region, bond_distance=d,
                           original_chargeprod_e2=e['chargeprod_e2'],
                           chargeprod_e2=chargeprod, coulomb_scale=scale,
                           action='omit' if region == 'QM' else 'replace_generic_pair'))
    caps = [a for a in core if a['kind'] == 'cap']
    id_map = {a['id']: a['physical_index'] for a in physical}
    for a in caps:
        if a['retained_source_id'] not in id_map or a['omitted_source_id'] not in id_map:
            raise ValueError('cap source mapping absent')
    selected = [t for t in bonded if any(i not in qm for i in t['atoms'])]
    report = {
        'protocol_id': PROTOCOL, 'status': 'parameter_ledger_only_no_energies',
        'inputs': {k: record(v) for k,v in dict(source_manifest=source_manifest,
                    protein_result=protein_result,pqq_cross_lj=pqq_file,metal_12_6=metal_file,
                    implementation=__file__).items()},
        'native_system': pins['system.xml'], 'forcefield':pr['forcefield'],
        'boundary_mapping':c['boundary_mapping'], 'source_states':c['source_states'],
        'environments':c['environments'],
        'counts': {'physical_atoms':len(physical), 'real_QM':sum(a['region']=='QM' for a in physical),
                   'MM':sum(a['region']=='MM' for a in physical),'caps':len(caps),
                   'zero_charge_boundary_MM':len(removed-qm),
                   'bonded_retained':dict(Counter(t['force'] for t in selected)),
                   'bonded_omitted_QM':dict(Counter(t['force'] for t in bonded if all(i in qm for i in t['atoms']))),
                   'exception_regions':dict(Counter(e['region'] for e in exrows))},
        'MM_charge_e':sum(a['charge_e'] for a in physical if a['region']=='MM'),
        'pair_policy': {'mixing':'Lorentz-Berthelot: sigma arithmetic; epsilon geometric',
          'boundary':'finite nonperiodic NoCutoff; no reaction field; no long-range dispersion correction',
          'MM_MM':'LJ plus redistributed-charge Coulomb; native exceptions replace generic pairs',
          'real_QM_MM':'LJ only; protein exceptions replace generic pairs; PQQ/metal have no covalent MM exclusions',
          'QM_QM':'no FF nonbonded', 'caps':'no independent FF particle, charge, LJ or bonded term',
          'C4':'absent', 'QM_MM_classical_Coulomb':'absent'},
        'components': {'0':'embedded r2SCAN-3c with mapped caps and qtilde; DoEQ false',
                       '1':'retained protein harmonic bonds', '2':'retained protein angles',
                       '3':'retained protein torsions', '4':'retained protein CMAP',
                       '5':'MM-MM LJ', '6':'MM-MM Coulomb', '7':'real QM-MM LJ'},
        'qualification': {'native_ORCA_equivalence':False,'full_hybrid_force_checks':False,
                          'metal_LJ_QMMM_validated':False,'molecular_energy_force_calls':0},
        'limitations':['Capped electronic boundary approximates missing covalent chemistry; no exact cancellation claim',
                       'Pure 12-6 hydration-fit metal LJ is an explicit unqualified QM/MM candidate',
                       'Fixed protein charges respond by movement, not electronic polarization',
                       'Finite dry system: no solvent-consistent whole-protein relaxation model; no optimization authorized',
                       'Original embedded electronic rigid-transform gate remains failed']}
    return report,physical,exrows,selected,caps


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('source-manifest','protein-result','pqq-cross-lj','metal-parameters','output'):
        p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args()
    report,particles,exceptions,bonded,caps=build(a.source_manifest,a.protein_result,a.pqq_cross_lj,a.metal_parameters)
    a.output.mkdir(parents=True,exist_ok=False)
    for name,data in [('particles.json',particles),('exceptions.json',exceptions),('bonded_terms.json',bonded),('caps.json',caps)]:
        write_new(a.output/name,data)
    report['artifacts']={f.name:record(f) for f in a.output.iterdir()}
    write_new(a.output/'LEDGER.json',report)
    print(json.dumps({'status':report['status'],'counts':report['counts'],'MM_charge_e':report['MM_charge_e']},indent=2))

if __name__=='__main__':main()
