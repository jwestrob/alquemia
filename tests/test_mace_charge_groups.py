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

PREP = ROOT / 'workspaces/mace_charge_groups_20260918/groups_v1/preparation.json'
TRACE = ROOT / 'workspaces/mace_response_trace_20260916/medium_v1/collection_job_1200676.json'


@unittest.skipUnless(PREP.exists() and TRACE.exists(), 'pinned real preparations/traces unavailable')
class GroupTests(unittest.TestCase):
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
