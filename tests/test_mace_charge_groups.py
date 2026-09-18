"""Real source-group and saved-tensor tests; no molecular energy inference."""
from pathlib import Path
import sys
import tempfile
import types
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, record
from mace_charge_groups import prepared, validate, kernel_parent_gate
from mace_group_constraints import GroupRestoration, configure
from mace_hybrid import accepted_attempt, EV_TO_KCAL

PREP = ROOT / 'workspaces/mace_charge_groups_20260918/groups_v1/preparation.json'
TRACE = ROOT / 'workspaces/mace_response_trace_20260916/medium_v1/collection_job_1200676.json'


@unittest.skipUnless(PREP.exists() and TRACE.exists(), 'pinned real preparations/traces unavailable')
class GroupTests(unittest.TestCase):
    @unittest.skipUnless(PREP.parent.parent.joinpath('model_v3/execution/GGR_1GLG_La_one_group/attempt_001/result.json').exists(),
                         'actual native-limit pilot evaluations unavailable')
    def test_actual_native_limit_energy_density_and_force(self):
        mp = PREP.parent.parent / 'model_v3/manifest.json'; m = read_json(mp)
        old = read_json(ROOT / 'workspaces/mace_global_benchmark_20260916/mace_v1/medium/collection_job_1200701.json')
        for metal in ('Ca', 'La'):
            tid = 'GGR_1GLG_' + metal + '_one_group'
            task = next(t for t in m['tasks'] if t['task_id'] == tid)
            actual = accepted_attempt(mp.parent / 'execution' / tid / 'attempt_001', task, mp)
            self.assertIsNotNone(actual)
            native = old['rows'][tid.replace('one_group', 'primary')]
            self.assertLessEqual(abs(actual['energy_eV'] - native['energy_eV']) * EV_TO_KCAL, .01)
            np.testing.assert_allclose(np.load(verify(actual['density_coefficients'])),
                                       np.load(verify(native['density_coefficients'])), atol=1e-9, rtol=0)
            np.testing.assert_allclose(np.load(verify(actual['forces'])),
                                       np.load(verify(native['forces'])), atol=.001, rtol=0)

    @unittest.skipUnless(PREP.parent.parent.joinpath('report_v1/result.json').exists(),
                         'complete actual model/solvent report unavailable')
    def test_actual_report_energy_algebra_and_replay(self):
        from mace_charge_groups import report
        root = PREP.parent.parent
        saved = read_json(root / 'report_v1/result.json')
        self.assertEqual(saved['status'], 'complete')
        self.assertEqual(len(saved['contrasts']), 7)
        self.assertEqual(len(saved['grouping_checks']), 3)
        for case, variants in saved['scores'].items():
            for variant, score in variants.items():
                ca, la = [case + '_' + metal + '_' + variant for metal in ('Ca', 'La')]
                molecular = saved['MACE']['rows']; solvent = saved['solvent']['rows']
                vacuum = (molecular[ca]['energy_eV'] - molecular[la]['energy_eV']) * EV_TO_KCAL
                correction = solvent[ca]['GB_reaction_kcal_mol'] - solvent[la]['GB_reaction_kcal_mol']
                for endpoint in (ca, la):
                    self.assertEqual(solvent[endpoint]['GB_reaction_kcal_mol'],
                                     solvent[endpoint]['GB_reaction_kJ_mol'] / 4.184)
                    self.assertFalse(solvent[endpoint]['direct_Coulomb_included'])
                self.assertAlmostEqual(score['R_candidate_model_kcal'], vacuum + correction, places=9)
                self.assertIsNone(score['calibrated_class'])
        for contrast in saved['contrasts']:
            hi, lo = [saved['scores'][contrast[k]]['primary']['R_candidate_model_kcal']
                      for k in ('higher_expected', 'lower_expected')]
            self.assertEqual(contrast['margin_model_kcal'], hi - lo)
            self.assertEqual(contrast['pass'], hi - lo > .02)
        for check in saved['grouping_checks']:
            scores = saved['scores'][check['case']]
            delta = scores['connected']['R_candidate_model_kcal'] - scores['primary']['R_candidate_model_kcal']
            self.assertEqual(check['connected_minus_primary_model_kcal'], delta)
            self.assertEqual(check['pass'], abs(delta) <= 2.)
        with tempfile.TemporaryDirectory() as d:
            report(root / 'model_v3/manifest.json', verify(saved['native_reference']),
                   Path(d) / 'replay', root / 'solvent_v1/manifest.json')
            replay = read_json(Path(d) / 'replay/result.json')
            for key in ('scores', 'contrasts', 'checks', 'grouping_checks', 'raw_pass_count',
                        'qualified_pass_count', 'numerical_checks_pass', 'representation_checks_pass'):
                self.assertEqual(replay[key], saved[key])
            # With actual molecular results but no supplied solvent, every
            # candidate score must remain unavailable, never vacuum or zero.
            report(root / 'model_v3/manifest.json', verify(saved['native_reference']),
                   Path(d) / 'missing_solvent')
            missing = read_json(Path(d) / 'missing_solvent/result.json')
            self.assertEqual(missing['status'], 'incomplete')
            self.assertEqual(missing['qualified_pass_count'], 0)
            for variants in missing['scores'].values():
                for score in variants.values():
                    self.assertIsNone(score['R_candidate_model_kcal'])
                    self.assertIsNone(score['GB_Ca_minus_La_kcal'])

    @unittest.skipUnless(PREP.parent.parent.joinpath('model_v2/manifest.json').exists(), 'real pilot manifest unavailable')
    def test_actual_manifest_execution_gate_and_corrupted_charge(self):
        import copy
        import json
        path = PREP.parent.parent / 'model_v2/manifest.json'
        self.assertEqual(validate(path)['tasks'], 24)
        self.assertEqual(kernel_parent_gate(read_json(path))['status'], 'pass')
        broken = copy.deepcopy(read_json(path))
        broken['tasks'][2]['group_charges_e'][0] += 1
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'corrupted_real_manifest.json'; p.write_text(json.dumps(broken))
            with self.assertRaises(InvalidArtifact):
                validate(p)

    def test_actual_complete_group_inventory(self):
        p, cases, cfg = prepared(PREP)
        self.assertEqual(len(cases), 7)
        for name, (row, physical) in cases.items():
            g = row['primary']; n = len(physical['physical_atoms'])
            self.assertEqual(sorted(sum(g['group_members'], [])), list(range(n)))
            delta = np.array(g['endpoint_group_charges_e']['La']) - g['endpoint_group_charges_e']['Ca']
            self.assertEqual(np.flatnonzero(delta).tolist(), [g['site_group_index']])
            self.assertEqual(delta.sum(), 1)
            for metal in ('Ca', 'La'):
                self.assertEqual(sum(g['endpoint_group_charges_e'][metal]), physical['endpoints'][metal]['charge'])
            if name.startswith('GGR'):
                self.assertTrue(set(g['selected_residues']) <= set(row['connected']['selected_residues']))

    def test_actual_trace_one_group_identity_and_corruption(self):
        import torch
        archive = read_json(TRACE)
        manifest = read_json(verify(archive['manifest']))
        charges = {t['task_id']: t['charge'] for t in manifest['tasks']}
        for key, endpoint in archive['rows'].items():
            trace = read_json(verify(endpoint['charge_trace']))
            arrays = np.load(verify(trace['arrays']))
            n = len(arrays['initial_weights_output'])
            observer = GroupRestoration([0] * n, [charges[key]])
            local = sum(arrays[name].reshape(n, 2, -1) for name in arrays.files
                        if name.startswith('local_source_') and name.endswith('_output'))
            for stage in range(3):
                if stage == 0:
                    density, weights = local, arrays['initial_weights_output']
                else:
                    update = arrays[f'update_{stage-1}_output']
                    density = arrays[f'update_{stage-1}_local_charges'].reshape(n, 2, -1) + update[:, :-2].reshape(n, 2, -1)
                    weights = update[:, -2:]
                got = observer(torch.tensor(density), torch.tensor(weights)).numpy()
                np.testing.assert_allclose(got, arrays[f'stage_{stage}_restored_density_internal'], atol=1e-9, rtol=0)
                np.testing.assert_array_equal(got[:, :, 1:], density[:, :, 1:])
            self.assertEqual(len(observer.receipt()['stages']), 3)
            # Explicitly corrupted real weights test the unsupported state;
            # they are not computed charges or a successful scientific fixture.
            broken = torch.tensor(weights).clone(); broken.zero_()
            with self.assertRaises(InvalidArtifact):
                GroupRestoration([0] * n, [charges[key]])(torch.tensor(density), broken)

    def test_pinned_native_method_can_be_adapted_without_forward(self):
        import torch
        software = read_json(ROOT / 'workspaces/mace_hybrid_20260916/software_v1/software_manifest.json')
        checkpoint = verify(software['checkpoint'])
        model = torch.load(checkpoint, map_location='cpu', weights_only=False)
        parameters = {name: p._version for name, p in model.named_parameters()}
        import inspect
        source = record(inspect.getsourcefile(type(model).forward))
        row = read_json(verify(read_json(PREP)['cases']['GGR_1GLG']))['primary']
        with tempfile.TemporaryDirectory() as d:
            calc = types.SimpleNamespace(models=[model], use_compile=False)
            observer, receipt = configure(calc, row['atom_group_indices'],
                                           row['endpoint_group_charges_e']['Ca'], source['sha256'], Path(d))
            self.assertEqual(receipt['native_restoration_blocks_replaced'], 2)
            self.assertEqual(observer.trace, [])
            self.assertEqual(parameters, {name: p._version for name, p in model.named_parameters()})
            import mace.modules.extensions as extensions
            self.assertIs(model.forward.__func__.__globals__, extensions.__dict__)
            from mace_realspace_compat import configure as realspace
            from mace_analytic import configure as analytic
            realspace(calc)
            analytic(calc, 256, 4096, 256)
            self.assertEqual(model.forward.__func__.__name__, 'inference')
            self.assertEqual(observer.trace, [])


if __name__ == '__main__':
    unittest.main()
