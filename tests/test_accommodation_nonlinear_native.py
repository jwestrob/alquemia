"""Real archived gradients/receipts/algebra only; no manufactured native result."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import accommodation_nonlinear_native as native
from affordable_common import BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, read_json, record, verify, xyz
from affordable_response import extract, read_engrad
from hydration_square import endpoint
from mace_site_kinematics import Kinematics


class NativeAdapter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.design = ROOT/'workspaces/accommodation_torsion_20260920/prepared_v3/design.json'
        cls.controls = ROOT/'workspaces/second_shell_20260919/prepared_v2/manifest.json'
        cls.torsion = ROOT/'workspaces/accommodation_torsion_20260920/prepared_v3/dft/manifest.json'
        cls.result = read_json(ROOT/'workspaces/accommodation_torsion_20260920/partial_DFT_result_v3.json')
        cls.water = ROOT/'workspaces/hydration_basin_20260919/validation_v1/dft/manifest.json'

    def test_actual_expanded_origins_state_precision_and_scf(self):
        rows = native.origins(self.design, self.controls, self.torsion)
        self.assertEqual(len(rows), 8)
        for r in rows:
            self.assertLessEqual(r['executed_vs_source_max_abs_A'], 1e-12)
            if r['status'] == 'complete':
                self.assertEqual(r['SCF_tolerances']['TolE'], 1e-6)
                self.assertEqual(r['SCF_tolerances']['TolG'], 5e-5)
            else:
                self.assertIsNone(r['result'])
                self.assertIn('reason', r)
        control = [r for r in rows if r['case_id'] == '1H4I']
        self.assertTrue(all(r['status'] == 'complete' for r in control))
        self.assertEqual([r['result']['energy_hartree'] for r in control], [-4916.025198080411, -4269.95639325921])
        self.assertTrue(all(r['executed_vs_source_max_abs_A'] == 0 for r in control))

    def test_new_recipe_adds_engrad_without_tightscf(self):
        for c in read_json(self.design)['cases']:
            for z in ('Ca', 'La'):
                e = c['origins'][z]; text = native.recipe(e['charge'], e['multiplicity'])
                self.assertEqual(text.splitlines()[0], native.BASE_HEADER+' EnGrad')
                self.assertNotIn('TightSCF', text)
                with tempfile.TemporaryDirectory() as d:
                    path = Path(d)/'endpoint.inp'; path.write_text(text)
                    native.check_recipe(path, e['charge'], e['multiplicity'], True)
                    # Explicit corrupted copy of this real-state input.
                    path.write_text(text.replace(' EnGrad', ' TightSCF EnGrad'))
                    with self.assertRaises(InvalidArtifact): native.check_recipe(path, e['charge'], e['multiplicity'], True)

    def test_actual_native_drivers_and_gradient_units(self):
        manifest = read_json(self.water)
        for z in ('Ca', 'La'):
            t = next(t for t in manifest['tasks'] if t['task_id'] == '1F6S__11__'+z+'__soft__p')
            out = Path(t['output_path']); receipt = read_json(Path(str(out)+'.execution.json'))
            drivers = native.analytic_drivers(out.read_text(), z)
            self.assertTrue(all(drivers.values()))
            self.assertEqual('ECP' in drivers, z == 'La')
            # This real hydration fixture deliberately used TightSCF. Reuse its
            # parser evidence, never its numerical recipe for the new pilot.
            self.assertEqual(native.scf_tolerances(out.read_text())['TolE'], 1e-8)
            gp = verify(receipt['artifacts']['engrad'])
            raw = read_engrad(gp)
            parsed = extract(gp, out, verify(t['input']), verify(t['xyz']))
            np.testing.assert_array_equal(parsed['gradient_kcal_mol_per_A'], raw['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A)
            actual = endpoint(record(out), record(Path(str(out)+'.execution.json')), t['xyz'], t['input'])
            self.assertAlmostEqual(parsed['energy_hartree'], actual['energy_hartree'], places=8)
            corrupted = out.read_text().replace('CPCM gradient', 'REMOVED gradient')
            with self.assertRaises(InvalidArtifact): native.analytic_drivers(corrupted, z)
            with self.assertRaises(InvalidArtifact): native.analytic_drivers(out.read_text()+'\nnumerical differentiation\n', z)

    def test_projection_replays_real_cartesian_gradient_at_current_q(self):
        d = read_json(ROOT/'workspaces/accommodation_response_20260920/prepared_v2/design.json')
        col = read_json(ROOT/'workspaces/accommodation_response_20260920/prepared_v2/gate/collection_1203465.json')
        r = next(r for r in col['rows'] if r['metal'] == 'Ca' and r['medium'] == 'alpb')
        mapping = d['maps']['1H4I__Ca']; kin = Kinematics(read_json(verify(mapping))['context'])
        i = next(i for i, m in enumerate(kin.modes) if m['id'] == 'A/177/chi3')
        q = np.zeros(len(kin.modes)); g = np.asarray(r['gradient_kcal_mol_A'])
        got = native.project_gradient(mapping, q, g, [i])
        self.assertAlmostEqual(got['active_gradient_kcal_mol_rad'][0], r['projected_gradient_kcal_per_unit'][i], places=10)
        q[i] = .1
        moved = native.project_gradient(mapping, q, g, [i])['active_gradient_kcal_mol_rad'][0]
        # Pure coordinate-chain-rule check with the same real archived gradient;
        # this is not a claimed molecular gradient at a displaced state.
        dq = np.zeros_like(q); dq[i] = 1e-6
        direct = float(np.sum(g*(kin.evaluate(q+dq)[1]-kin.evaluate(q-dq)[1]))/2e-6)
        self.assertAlmostEqual(moved, direct, places=6)
        self.assertGreater(abs(moved-got['active_gradient_kcal_mol_rad'][0]), .001)
        self.assertTrue(all(not r['active'] for r in got['all_mode_derivatives'] if r['id'].startswith('metal_')))

    def test_actual_native_work_and_ca_minus_la_sign(self):
        points = {r['task_id']: r for r in self.result['points']}
        for profile in self.result['profiles']:
            if profile['case_id'] != 'PQQSEQ_83440678cbbd658047c9' or profile['role'] != 'extra_acidic_ligand_homolog': continue
            for r in profile['rows']:
                if abs(r['angle_radian']) != .2: continue
                works = {}
                for z, p in r['endpoints'].items():
                    e = points[p['task_id']]['DFT_hartree']; e0 = points[profile['case_id']+'__origin__'+z]['DFT_hartree']
                    x = native.work_comparison(e, e0, p['work_kcal_mol'])
                    self.assertEqual(x['native_work_kcal_mol'], p['DFT_work_kcal_mol'])
                    self.assertAlmostEqual(x['model_minus_native_work_kcal_mol']['composite'], p['work_kcal_mol']['composite']-p['DFT_work_kcal_mol'])
                    missing = native.work_comparison(e, None, p['work_kcal_mol'])
                    self.assertIsNone(missing['native_work_kcal_mol'])
                    self.assertIsNone(missing['model_minus_native_work_kcal_mol'])
                    works[z] = x['native_work_kcal_mol']
                self.assertEqual(works['Ca']-works['La'], r['DFT_Ca_minus_La_work_kcal_mol'])
                self.assertEqual(np.sign(works['Ca']-works['La']), np.sign(r['angle_radian']))


if __name__ == '__main__': unittest.main()
