"""Real GGR artifact/geometry tests; quantum gradient tests require real outputs.

Coordinate finite differences below test geometric Jacobians only. They do not
evaluate DFT gradients, create invented energies, or validate a response model.
"""
import copy
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import gemmi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, energy, paired, read_json, record, verify, write_new, xyz
from affordable_response import bounded_response, paired_sensitivity, read_engrad
from ggr_sensitivity import (AMPLITUDES, collect, consistency, endpoint_output_validation,
    input_state, materialize_motion, moved_source, physical_velocity, prepare, prepare_half,
    source_key, validate_motion, validate_representation)

BASE = ROOT / 'workspaces/ggr_mechanism_20260915'
PREP = BASE / 'stage_a_prepared_v1'
REPRESENTATIONS = {'extended': PREP/'ggr_1glg_nma/preparation_manifest.json',
                   'connected': PREP/'ggr_1glg_connected/preparation_manifest.json'}
AGREEMENT = ROOT/'diagnostics/ggr_mechanism_plan_20260915/AGREEMENT.md'
PLAN = ROOT/'diagnostics/ggr_mechanism_plan_20260915/PLAN.md'
OLD_REPAIR = ROOT/'workspaces/affordable_challenger_20260915/verified_repairs/ggr_1glg_GGR/repair_manifest.json'
OLD_RUN = ROOT/'workspaces/baseline_benchmark_20260915/run_v1/ggr_1glg_GGR/repaired'


class PhysicalGGR(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not all(p.exists() for p in REPRESENTATIONS.values()):
            raise unittest.SkipTest('real approved Stage A preparations unavailable')
        cls.models = {name: validate_representation(name, path) for name, path in REPRESENTATIONS.items()}

    def test_native_centers_reconstruct_without_geometry_change(self):
        for name, (r, g, coords, metal, mk, specs) in self.models.items():
            atoms, jacobian = materialize_motion(r, coords, metal, {})
            np.testing.assert_allclose([a[1:] for a in atoms],
                                      [a[1:] for a in xyz(verify(r['outputs']['La']['xyz']))], atol=1e-9, rtol=0)
            self.assertEqual(len(atoms), 58 if name == 'extended' else 111)
            self.assertTrue(np.all(jacobian == 0))

    def test_geometric_jacobians_in_physical_units(self):
        for r, g, coords, metal, mk, specs in self.models.values():
            for coordinate, spec in specs.items():
                velocity = physical_velocity(coords, spec)
                _, jacobian = materialize_motion(r, coords, metal, velocity)
                step = 1e-6  # Angstrom or radians; purely geometric, not electronic.
                displacements = []
                for sign in (-1, 1):
                    c, m = moved_source(coords, metal, spec, sign*step)
                    atoms, _ = materialize_motion(r, c, m, {})
                    displacements.append(np.asarray([a[1:] for a in atoms]))
                np.testing.assert_allclose(jacobian, (displacements[1]-displacements[0])/(2*step), atol=1e-8, rtol=0)
                caps = [row['qm_index'] for row in r['atom_graph']['source_to_qm'] if row['kind'] == 'sigma_link_H']
                self.assertTrue(np.all(jacobian[caps] == 0))

    def test_crankshaft_preserves_real_peptide_and_anchor_bond_lengths(self):
        for r, g, coords, metal, mk, specs in self.models.values():
            spec = specs['peptide']
            moving = spec['moving_source_keys']
            a, b = spec['fixed_anchor_keys']
            for sign in (-1, 1):
                displaced, newmetal = moved_source(coords, metal, spec, sign*AMPLITUDES['peptide'])
                np.testing.assert_array_equal(newmetal, metal)
                for key in coords:
                    if key not in moving:
                        np.testing.assert_array_equal(displaced[key], coords[key])
                for i, left in enumerate(moving):
                    for right in moving[i+1:]:
                        self.assertAlmostEqual(np.linalg.norm(coords[left]-coords[right]),
                                               np.linalg.norm(displaced[left]-displaced[right]), places=12)
                for anchor, child in [(a, moving[0]), (b, moving[2])]:
                    self.assertAlmostEqual(np.linalg.norm(coords[anchor]-coords[child]),
                                           np.linalg.norm(displaced[anchor]-displaced[child]), places=12)

    def test_nominal_and_half_membership_trust_region(self):
        for r, g, coords, metal, mk, specs in self.models.values():
            for coordinate, spec in specs.items():
                for sign in (-1, 1):
                    for scale in (1, .5):
                        atoms, validation = validate_motion(r, g, coords, metal, mk, spec, sign*scale*AMPLITUDES[coordinate])
                        self.assertLessEqual(validation['maximum_source_heavy_atom_displacement_A'], .05)
                        self.assertEqual(validation['coordination_membership'], 'unchanged')
                        self.assertEqual(validation['explicit_water_inventory'], [])

    def test_large_motion_rejected_without_alternative_coordinate(self):
        r,g,c,m,mk,s = self.models['extended']
        with self.assertRaisesRegex(InvalidArtifact, 'trust region'):
            validate_motion(r,g,c,m,mk,s['metal'],.1)

    def test_manifest_recipe_missing_outputs_and_no_false_half_step(self):
        from run_orca_task_manifest import load_manifest_tasks
        from affordable_workflow import dry_run
        with tempfile.TemporaryDirectory(prefix='stage_c_software_test_', dir=BASE) as directory:
            path = Path(directory)/'nominal'
            manifest = prepare(REPRESENTATIONS, path, AGREEMENT, PLAN)
            mp = path/'manifest.json'
            self.assertEqual(len(manifest['tasks']), 20)
            self.assertEqual(len(load_manifest_tasks(mp)[1]), 20)
            self.assertEqual(dry_run(mp)['tasks'], 20)
            self.assertEqual(sum(t['task_type'] == 'analytic_gradient' for t in manifest['tasks']), 4)
            self.assertEqual(len({t['cache_key'] for t in manifest['tasks']}), 20)
            for task in manifest['tasks']:
                state = input_state(verify(task['input']), gradient=task['task_type']=='analytic_gradient', tight=True)
                self.assertEqual(state['charge'], task['charge'])
            for name in self.models:
                subset = {t['task_id']: t for t in manifest['tasks'] if t['representation'] == name}
                for suffix in ('center', 'metal_minus', 'metal_plus', 'peptide_minus', 'peptide_plus'):
                    la, ca = (subset[f'{name}_{suffix}_{metal}'] for metal in ('La','Ca'))
                    paired(verify(la['xyz']), verify(ca['xyz']), 0, -1)
            comparison = collect(mp)
            self.assertTrue(all(row['status']=='unavailable' for row in comparison['comparisons']))
            self.assertEqual(comparison['half_step_blocks'], [])
            cp = Path(directory)/'missing_results.json'
            write_new(cp, comparison)
            self.assertEqual(prepare_half(mp, cp, Path(directory)/'half')['status'], 'no_failed_valid_blocks')
            self.assertIsNone(comparison['relaxation_correction_kcal_mol'])
            self.assertIsNone(comparison['entropy_correction_kcal_mol'])


class ArchivedEndpointChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repair = read_json(OLD_REPAIR)

    def task(self, metal):
        return {**self.repair['outputs'][metal], 'metal':metal}

    def output(self, metal):
        return OLD_RUN/metal/f'ggr_1glg_GGR_amide_v3_{metal}.out'

    def test_real_native_ecp_and_charge_state(self):
        for metal in ('La','Ca'):
            result = endpoint_output_validation(self.task(metal), self.output(metal), require_tight=False)
            self.assertEqual(result['ecp_core_electrons'], 46 if metal=='La' else 0)
            self.assertEqual(result['explicit_electrons']%2,0)

    def test_corrupted_real_output_ecp_and_charge_fail(self):
        original = self.output('La').read_text()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'explicitly_corrupted_actual_GGR_output.out'
            for corrupted in [original.replace('replacing 46 core electrons','replacing 45 core electrons'),
                              original.replace('Total Charge           Charge          ....    0',
                                               'Total Charge           Charge          ....    1')]:
                self.assertNotEqual(corrupted, original)
                path.write_text(corrupted)
                with self.assertRaises(InvalidArtifact):
                    endpoint_output_validation(self.task('La'), path, require_tight=False)

    def test_real_single_point_cannot_masquerade_as_tight_gradient(self):
        with self.assertRaises(InvalidArtifact):
            endpoint_output_validation(self.task('La'), self.output('La'), require_gradient=True)

    def test_actual_energy_algebra_and_scaled_half_tolerance_only(self):
        # Algebra regression using two archived energies; these are NOT a real
        # displacement experiment or an invented scientific gradient fixture.
        eca, ela = (energy(self.output(m)) for m in ('Ca','La'))
        difference = (eca-ela)*HA_TO_KCAL
        result = consistency(eca, ela, ela, difference/.04, .02)
        self.assertAlmostEqual(result['odd_energy_kcal_mol'], difference/2)
        self.assertAlmostEqual(result['predicted_odd_energy_kcal_mol'], difference/2)
        self.assertAlmostEqual(result['residual_kcal_mol'],0.,places=9)
        nominal = consistency(eca,eca,eca,0.,.02)
        half = consistency(eca,eca,eca,0.,.01,half=True)
        self.assertEqual(nominal['tolerance_kcal_mol']/.02,half['tolerance_kcal_mol']/.01)

    def test_curvature_cannot_enable_response(self):
        self.assertIsNone(bounded_response(None,None)['relaxation_energy_kcal_mol'])


class RealGradientIntegration(unittest.TestCase):
    def test_actual_missing_gradient_component_rejected(self):
        manifest = BASE/'stage_c_tasks_v1/manifest.json'
        if not manifest.exists():
            self.skipTest('approved Stage C manifest not prepared')
        tasks = [t for t in read_json(manifest)['tasks'] if t['task_type']=='analytic_gradient'
                 and Path(str(t['output_path'])+'.execution.json').exists()]
        if not tasks:
            self.skipTest('no completed approved analytic-gradient endpoint')
        task = tasks[0]
        original = Path(task['output_path']).read_text()
        with tempfile.TemporaryDirectory() as directory:
            corrupted = Path(directory)/'explicitly_corrupted_missing_CPCM_gradient.out'
            corrupted.write_text(original.replace('CPCM gradient','REMOVED_COMPONENT'))
            with self.assertRaisesRegex(InvalidArtifact,'analytic-gradient component'):
                endpoint_output_validation(task,corrupted,require_gradient=True)

    def test_actual_paired_gradient_and_corrupted_charge_rejected(self):
        manifest = BASE/'stage_c_tasks_v1/manifest.json'
        if not manifest.exists():
            self.skipTest('approved Stage C manifest not prepared')
        result = collect(manifest)
        keys = ['extended_center_La','extended_center_Ca']
        if any(k not in result['gradients'] for k in keys):
            self.skipTest('both real extended-center gradients not yet complete')
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory)/f'{metal}.json' for metal in ('La','Ca')]
            for key,path in zip(keys,paths):
                write_new(path,result['gradients'][key])
            pair = paired_sensitivity(*paths)
            self.assertEqual(pair['quantity'],'grad_E_Ca_minus_E_La')
            self.assertIsNone(pair['relaxation_correction_kcal_mol'])
            corrupt = copy.deepcopy(result['gradients'][keys[1]])
            original = verify(corrupt['artifacts']['input']).read_text()
            badinput = Path(directory)/'explicitly_corrupted_Ca_charge.inp'
            badinput.write_text(original.replace('* xyzfile -1 1','* xyzfile 0 1'))
            self.assertNotEqual(badinput.read_text(),original)
            corrupt['artifacts']['input'] = record(badinput)
            badrecord = Path(directory)/'explicitly_corrupted_Ca_gradient.json'
            write_new(badrecord,corrupt)
            with self.assertRaisesRegex(InvalidArtifact,'charge/multiplicity invariant'):
                paired_sensitivity(paths[0],badrecord)

    def test_real_engrad_layout_and_coordinate_validation(self):
        paths = [p for p in sorted(BASE.glob('stage_c*/**/endpoint.engrad'))
                 if (p.parent/'endpoint.out.execution.json').exists()]
        if not paths:
            self.skipTest('no executed approved Stage C analytic-gradient artifact yet')
        for path in paths:
            parsed = read_engrad(path)
            atoms = xyz(path.parent/'core.xyz')
            self.assertEqual(parsed['atom_count'],len(atoms))
            self.assertTrue(np.isfinite(parsed['gradient_Ha_per_bohr']).all())
            np.testing.assert_array_equal(parsed['atomic_numbers'],[gemmi.Element(a[0]).atomic_number for a in atoms])
            np.testing.assert_allclose(parsed['coordinates_bohr']*BOHR_TO_A,
                                      [a[1:] for a in atoms],atol=1e-6,rtol=0)
            self.assertAlmostEqual(parsed['energy_Ha'],energy(path.parent/'endpoint.out'),places=8)
        manifests = {p.parent.parent/'manifest.json' for p in paths}
        for manifest in manifests:
            result = collect(manifest)
            self.assertTrue(result['gradients'])
            for gradient in result['gradients'].values():
                self.assertEqual(gradient['state_validation']['orca_version'],'6.1.1')
                self.assertIsNone(gradient['relaxation_correction_kcal_mol'])

    def test_actual_engrad_format_and_corrupted_copy(self):
        paths = [p for p in sorted(BASE.glob('stage_c*/**/endpoint.engrad'))
                 if (p.parent/'endpoint.out.execution.json').exists()]
        if not paths:
            self.skipTest('no executed approved Stage C analytic-gradient artifact yet')
        path = paths[0]
        original = read_engrad(path)
        tokens = [token for line in path.read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith('#') for token in line.split()]
        with tempfile.TemporaryDirectory() as directory:
            scalar = Path(directory)/'actual_gradient_scalar_layout.engrad'
            scalar.write_text('\n'.join(tokens)+'\n')
            np.testing.assert_array_equal(read_engrad(scalar)['gradient_Ha_per_bohr'],original['gradient_Ha_per_bohr'])
            corrupted = Path(directory)/'explicitly_corrupted_actual_gradient.engrad'
            corrupted.write_text('\n'.join(tokens[:-1])+'\n')
            with self.assertRaises(InvalidArtifact):
                read_engrad(corrupted)


if __name__ == '__main__':
    unittest.main()
