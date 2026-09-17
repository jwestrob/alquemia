"""Prepare real AMOEBA protein frameworks without evaluating energies or forces.

The selected quantum metal is retained in the physical ledger, explicitly
unparameterized. A framework pass is never a complete hybrid-model pass.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz

PROTOCOL = 'frozen_physical_AMOEBA2018_framework_capability_v1'
CASES = ('GGR_1GLG', 'ALPHA_1F6S', 'ALPHA_6IP9')
SETTINGS = dict(nonbondedMethod='NoCutoff', constraints=None, rigidWater=False,
                polarization='mutual', mutualInducedTargetEpsilon=1e-5,
                mutualInducedMaxIterations=60, removeCMMotion=False,
                energy_evaluations=0, force_evaluations=0, coordinate_changes=0)


def atom_id(atom):
    return f'{atom.residue.chain.id}/{atom.residue.id}/{atom.residue.insertionCode.strip()}/{atom.name}'


def topology(prep):
    """Replay the archived source graph and copy the normalized coordinates."""
    from openmm import app, unit
    if prep['case_id'] not in CASES or prep['preparation_details']['terminal_additions']:
        raise InvalidArtifact('unsupported case or added terminal topology')
    physical = prep['physical_atoms']
    by_id = {a['id']: a for a in physical}
    if len(by_id) != len(physical):
        raise InvalidArtifact('duplicate physical identity')
    if set(a['kind'] for a in physical) - {'protein_source', 'retained_site_water', 'selected_metal'}:
        raise InvalidArtifact('unsupported cofactor/physical atom kind')
    metals = [a for a in physical if a['kind'] == 'selected_metal']
    if len(metals) != 1 or metals[0]['id'] != 'metal':
        raise InvalidArtifact('expected one explicitly unparameterized QM metal')
    for metal in ('Ca', 'La'):
        endpoint = xyz(verify(prep['endpoints'][metal]['xyz']))
        expected = [(metal if a['kind'] == 'selected_metal' else a['element'], *a['xyz_A']) for a in physical]
        if [a[0] for a in endpoint] != [a[0] for a in expected] or not np.allclose(
                [a[1:] for a in endpoint], [a[1:] for a in expected], atol=1e-9, rtol=0):
            raise InvalidArtifact('paired endpoint coordinates or species differ')
    pdb = app.PDBFile(str(verify(prep['source'])))
    model = app.Modeller(pdb.topology, pdb.positions)
    original = prep['preparation_details']['original_protein_atoms']
    wanted = {a['id'] for a in original}
    model.delete([a for a in model.topology.atoms() if atom_id(a) not in wanted])
    atoms = list(model.topology.atoms())
    ids = [atom_id(a) for a in atoms]
    if ids != [a['id'] for a in original] or set(ids) != {
            a['id'] for a in physical if a['kind'] == 'protein_source'}:
        raise InvalidArtifact('source protein inventory differs')
    source_xyz = np.asarray(model.positions.value_in_unit(unit.angstrom))
    for a, old in zip(atoms, original):
        if a.element.symbol != old['element'] or not np.allclose(
                source_xyz[a.index], old['xyz_A'], atol=1e-12, rtol=0):
            raise InvalidArtifact('source coordinates or elements differ')
        if a.element.symbol != 'H' and old['xyz_A'] != by_id[atom_id(a)]['xyz_A']:
            raise InvalidArtifact('source heavy coordinates changed')
    bonds = sorted(tuple(sorted((atom_id(a), atom_id(b)))) for a, b in model.topology.bonds())
    expected_bonds = sorted(tuple(sorted((b['atom_a_id'], b['atom_b_id'])))
                            for b in prep['preparation_details']['bonds'])
    if bonds != expected_bonds:
        raise InvalidArtifact('source connectivity/disulfides differ')
    # These waters already exist in the frozen physical preparation. No atoms,
    # positions or hydrogen geometries are generated here.
    waters = defaultdict(list)
    for a in physical:
        if a['kind'] == 'retained_site_water':
            waters[a['id'].rsplit('/', 1)[0]].append(a)
    if len(waters) != len(prep['explicit_waters']):
        raise InvalidArtifact('water inventory differs')
    chains = {c.id: c for c in model.topology.chains()}
    water_bonds = []
    for key, rows in waters.items():
        chain, resid, icode = key.split('/')
        if chain not in chains:
            chains[chain] = model.topology.addChain(chain)
        residue = model.topology.addResidue('HOH', chains[chain], resid, icode)
        named = {r['id'].rsplit('/', 1)[1]: r for r in rows}
        if set(named) != {'O', 'H1', 'H2'} or {n: r['element'] for n, r in named.items()} != {
                'O': 'O', 'H1': 'H', 'H2': 'H'}:
            raise InvalidArtifact('unsupported retained water topology')
        wa = {n: model.topology.addAtom(n, app.element.get_by_symbol(r['element']), residue)
              for n, r in named.items()}
        for n in ('H1', 'H2'):
            model.topology.addBond(wa['O'], wa[n])
            water_bonds.append(sorted((atom_id(wa['O']), atom_id(wa[n]))))
    ids = [atom_id(a) for a in model.topology.atoms()]
    if set(ids) | {'metal'} != set(by_id) or len(ids) + 1 != len(physical):
        raise InvalidArtifact('unaccounted physical atoms')
    positions = np.asarray([by_id[i]['xyz_A'] for i in ids])
    return model.topology, positions, ids, bonds, water_bonds, metals


def value(obj):
    from openmm import unit
    if unit.is_quantity(obj):
        return value(obj.value_in_unit_system(unit.md_unit_system))
    if isinstance(obj, (list, tuple)):
        return [value(x) for x in obj]
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return list(obj)


def prepare_case(pin, ff_paths, output):
    import openmm as mm
    from openmm import app, unit
    start, cpu = time.monotonic(), time.process_time()
    prep = read_json(verify(pin))
    result = dict(case_id=prep['case_id'], source_preparation=pin,
                  status='preparation_failed', energy=None, calibrated_class=None,
                  full_model_status='unvalidated_QM_source_boundary_and_damping',
                  new_energy_calls=0, new_force_calls=0, settings=SETTINGS)
    output.mkdir(exist_ok=False)
    try:
        top, pos, ids, bonds, water_bonds, metals = topology(prep)
        ff = app.ForceField(*map(str, ff_paths))
        unmatched = ff.getUnmatchedResidues(top)
        result['unmatched_residues'] = [dict(chain=r.chain.id, resid=r.id,
            insertion_code=r.insertionCode, name=r.name) for r in unmatched]
        if unmatched:
            raise InvalidArtifact('AMOEBA does not cover the unchanged framework topology')
        templates = ff.getMatchingTemplates(top)
        system = ff.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None,
            rigidWater=False, polarization='mutual', mutualInducedTargetEpsilon=1e-5,
            mutualInducedMaxIterations=60, removeCMMotion=False)
        if system.getNumParticles() != len(ids) or system.getNumConstraints():
            raise InvalidArtifact('unexpected particles or constraints')
        mp = next(f for f in system.getForces() if isinstance(f, mm.AmoebaMultipoleForce))
        gk = next(f for f in system.getForces() if isinstance(f, mm.AmoebaGeneralizedKirkwoodForce))
        parameters = []
        for i, pid in enumerate(ids):
            m = value(mp.getMultipoleParameters(i))
            g = value(gk.getParticleParameters(i))
            if abs(m[0] - g[0]) > 1e-12:
                raise InvalidArtifact('multipole/GK source charge mismatch')
            parameters.append(dict(id=pid, index=i, multipole_md_units=m, gk_md_units=g,
                covalent_maps={str(k): list(mp.getCovalentMap(i, k)) for k in range(8)}))
        charge = math.fsum(p['multipole_md_units'][0] for p in parameters)
        if abs(charge - prep['protein_charge_e']) > 1e-5:
            raise InvalidArtifact('framework formal charge changed')
        xml = output/'framework.xml'
        with xml.open('x') as f:
            f.write(mm.XmlSerializer.serialize(system))
        mapping = dict(physical_atoms=prep['physical_atoms'], system_atom_ids=ids,
            positions_A=pos.tolist(), protein_bonds=bonds, retained_water_bonds=water_bonds,
            unparameterized_QM_sources=metals, parameters=parameters,
            residue_templates=[dict(chain=r.chain.id, resid=r.id, insertion_code=r.insertionCode,
                source_name=r.name, template=t.name) for r, t in zip(top.residues(), templates)],
            assembly=prep['assembly'], microstate=prep['microstate'],
            explicit_waters=prep['explicit_waters'], evidence=prep['evidence'],
            evidence_use=prep['evidence_use'])
        write_new(output/'mapping.json', mapping)
        result.update(status='framework_parameterized_no_energy', framework=record(xml),
            mapping=record(output/'mapping.json'), physical_atoms=len(prep['physical_atoms']),
            framework_atoms=len(ids), unparameterized_atoms=len(metals),
            protein_charge_e=charge, explicit_waters=prep['explicit_waters'],
            force_classes=dict(Counter(type(f).__name__ for f in system.getForces())),
            gk_settings=dict(solute_dielectric=gk.getSoluteDielectric(),
                solvent_dielectric=gk.getSolventDielectric(), include_cavity=gk.getIncludeCavityTerm(),
                probe_radius_nm=value(gk.getProbeRadius()), surface_area_factor=value(gk.getSurfaceAreaFactor())),
            coordinate_max_change_A=0.0, added_atoms=0, removed_physical_atoms=0)
    except Exception as error:
        result['failure_reason'] = f'{type(error).__name__}: {error}'
        (output/'failure.txt').write_text(traceback.format_exc())
    result['receipt'] = dict(wall_seconds=time.monotonic()-start,
        process_cpu_seconds=time.process_time()-cpu, peak_process_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        gpu_seconds=0, scheduler_allocation=None)
    write_new(output/'result.json', result)
    return result


def prepare(inputs, plan, output):
    import openmm as mm
    from openmm import app
    from openmm.app.internal import amoebaforces
    if len(inputs) != 3 or [read_json(p)['case_id'] for p in inputs] != list(CASES):
        raise InvalidArtifact('exactly the three declared inputs in declared order required')
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    impl = root/'implementation'; impl.mkdir()
    for source in (Path(__file__), Path(__file__).with_name('affordable_common.py')):
        shutil.copyfile(source, impl/source.name)
    data = Path(app.__file__).parent/'data'
    ff_paths = [data/'amoeba2018.xml', data/'amoeba2018_gk.xml']
    config = dict(protocol=PROTOCOL, settings=SETTINGS, inputs=[record(p) for p in inputs],
        plan=record(plan), forcefields=[record(p) for p in ff_paths], openmm_version=mm.__version__,
        forcefield_builder=record(amoebaforces.__file__),
        forcefield_parser=record(Path(app.__file__).parent/'forcefield.py'),
        implementation=[record(p) for p in sorted(impl.glob('*.py'))], executable=str(Path(sys.executable).resolve()))
    write_new(root/'manifest.json', dict(**config, cache_key=cache_key(config)))
    results = [prepare_case(p, ff_paths, root/case) for p, case in zip(config['inputs'], CASES)]
    write_new(root/'result.json', dict(protocol=PROTOCOL, manifest=record(root/'manifest.json'),
        cases=results, framework_passes=sum(r['status']=='framework_parameterized_no_energy' for r in results),
        denominator=3, full_hybrid_qualified=False, numerical_score=None, new_energy_calls=0,
        new_force_calls=0, baseline_changed=False))
    return read_json(root/'result.json')


def validate(path):
    m = read_json(path)
    key = m.pop('cache_key')
    if cache_key(m) != key:
        raise InvalidArtifact('manifest key differs')
    for rec in m['inputs']+m['forcefields']+m['implementation']+[m['plan'],m['forcefield_builder'],m['forcefield_parser']]:
        verify(rec)
    return dict(status='pins_verified', protocol=m['protocol'], tasks=len(m['inputs']), energy_calls=0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--inputs', nargs=3, required=True)
    prep.add_argument('--plan', required=True)
    prep.add_argument('--output', required=True)
    dry = sub.add_parser('dry-run'); dry.add_argument('--manifest', required=True)
    report = sub.add_parser('report'); report.add_argument('--result', required=True)
    a = p.parse_args()
    if a.command == 'prepare':
        r = prepare(a.inputs, a.plan, a.output)
        print(json.dumps({k:r[k] for k in ('framework_passes','denominator','full_hybrid_qualified','new_energy_calls')}))
    elif a.command == 'dry-run':
        print(json.dumps(validate(a.manifest)))
    else:
        print(json.dumps(read_json(a.result), indent=2))


if __name__ == '__main__':
    main()
