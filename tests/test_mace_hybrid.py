"""Integrity/algebra tests on pinned real 1H4I artifacts, not model-output mocks."""
import copy
import json
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import mace_hybrid as mh
from affordable_common import InvalidArtifact, cache_key, energy, read_json, verify, xyz


@unittest.skipUnless((ROOT/'workspaces/mace_hybrid_20260916/pilot_v1/manifest.json').exists(),
                     'requires pinned real prepared 1H4I pilot artifacts')
class RealPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mp = ROOT/'workspaces/mace_hybrid_20260916/pilot_v1/manifest.json'
        cls.m = read_json(cls.mp)
        cls.tasks = {t['task_id']: t for t in cls.m['tasks']}

    def test_complete_manifest_pins(self):
        self.assertEqual(mh.dry_run(self.mp)['tasks'], 12)

    def test_paired_physical_coordinates_and_parity(self):
        for variant in ('primary', 'repeat', 'translate', 'rotate'):
            tasks = [self.tasks[f'1h4i_full_{metal}_{variant}'] for metal in ('La', 'Ca')]
            rows = [xyz(verify(t['xyz'])) for t in tasks]
            self.assertEqual(len(rows[0]), 9141)
            self.assertEqual(tasks[0]['charge']-tasks[1]['charge'], 1)
            self.assertEqual([a[1:] for a in rows[0]], [a[1:] for a in rows[1]])
            self.assertEqual(sum(a[0] != b[0] for a,b in zip(*rows)), 1)
            for t, atoms in zip(tasks, rows):
                self.assertEqual(mh.check_atoms(atoms, t['charge'])['all_electron_count'] % 2, 0)

    def test_real_core_parity_corruption_is_rejected(self):
        t = self.tasks['1h4i_qm33_Ca']
        with self.assertRaises(InvalidArtifact):
            mh.check_atoms(xyz(verify(t['xyz'])), t['charge']+1)

    def test_caps_have_no_full_atom_mapping(self):
        mappings = read_json(verify(self.m['atom_mappings']))
        for rows in mappings.values():
            self.assertTrue(any(a['kind'] == 'cap' for a in rows))
            for a in rows:
                self.assertEqual(a['physical_index'] is None, a['kind'] == 'cap')

    def test_rigid_transform_and_repeat_inventory(self):
        a = np.array([r[1:] for r in xyz(verify(self.tasks['1h4i_full_La_primary']['xyz']))])
        b = np.array([r[1:] for r in xyz(verify(self.tasks['1h4i_full_La_repeat']['xyz']))])
        self.assertTrue(np.array_equal(a, b))
        c = np.array([r[1:] for r in xyz(verify(self.tasks['1h4i_full_La_translate']['xyz']))])
        np.testing.assert_allclose(c-a, np.broadcast_to([10., -7., 3.], a.shape), atol=2e-14, rtol=0)
        matrix = mh.rotation()
        np.testing.assert_allclose(matrix @ matrix.T, np.eye(3), atol=1e-15)
        self.assertAlmostEqual(np.linalg.det(matrix), 1.)
        self.assertNotEqual(self.tasks['1h4i_full_La_primary']['cache_key'], self.tasks['1h4i_full_La_repeat']['cache_key'])

    def test_archive_energy_and_sign_replay(self):
        c = read_json(verify(self.m['archived_endpoints']))
        rows = {r['task_id']: r for r in c['rows']}
        for row in rows.values():
            self.assertEqual(energy(verify(row['output'])), row['energy_hartree'])
        contrasts = {}
        for size in (33, 36):
            contrasts[size] = (rows[f'1h4i_qm{size}_Ca']['energy_hartree'] - rows[f'1h4i_qm{size}_La']['energy_hartree'])*mh.HA_TO_KCAL
        self.assertAlmostEqual(contrasts[36]-contrasts[33], 61.738603357, places=6)
        self.assertAlmostEqual(mh.EV_TO_KCAL, 1.602176634e-19 * 6.02214076e23 / 4184, places=13)
        self.assertIsNone(self.m['reference'])

    def test_missing_attempt_is_never_cached(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            self.assertIsNone(mh.accepted_attempt(p, self.tasks['1h4i_qm33_Ca'], self.mp))
            # Explicit malformed receipt; no invented successful scientific output.
            (p/'receipt.json').write_text(json.dumps({'result': None}))
            self.assertIsNone(mh.accepted_attempt(p, self.tasks['1h4i_qm33_Ca'], self.mp))

    def test_model_settings_change_cache_identity(self):
        base = self.m['model']; original = cache_key(base)
        for field, value in [('dtype', 'float32'), ('pbc_handling', 'pbc'),
                             ('assembly', 'different_assembly'), ('microstate', 'different_microstate'),
                             ('explicit_waters', ['corrupted_inventory'])]:
            changed = copy.deepcopy(base); changed[field] = value
            self.assertNotEqual(cache_key(changed), original)

    @unittest.skipUnless((ROOT/'workspaces/mace_hybrid_20260916/pilot_v2/manifest.json').exists(),
                         'requires real versioned interface-repair preparation')
    def test_interface_repair_preserves_scientific_inputs(self):
        mp = ROOT/'workspaces/mace_hybrid_20260916/pilot_v2/manifest.json'
        newer = read_json(mp)
        self.assertEqual(mh.dry_run(mp)['tasks'], 12)
        for field in ('protocol_id', 'software', 'agreement', 'source_states', 'tolerances', 'run_inventory'):
            self.assertEqual(newer[field], self.m[field])
        model = dict(newer['model']); model.pop('execution_adapter')
        self.assertEqual(model, self.m['model'])
        for before, after in zip(self.m['tasks'], newer['tasks']):
            if before['variant'] == 'rotate':
                # Independently prepared rotations across NumPy environments can
                # differ at binary roundoff; the physical-source tolerance is 1e-12 A.
                old = xyz(verify(before['xyz'])); new = xyz(verify(after['xyz']))
                self.assertEqual([r[0] for r in old], [r[0] for r in new])
                np.testing.assert_allclose([r[1:] for r in old], [r[1:] for r in new], atol=1e-12, rtol=0)
            else:
                self.assertEqual(before['xyz']['sha256'], after['xyz']['sha256'])
            np.testing.assert_allclose(before['rotation_matrix'], after['rotation_matrix'], atol=1e-15, rtol=0)
            exclude = ('xyz', 'cache_key', 'rotation_matrix')
            self.assertEqual({k:v for k,v in before.items() if k not in exclude},
                             {k:v for k,v in after.items() if k not in exclude})
            self.assertNotEqual(before['cache_key'], after['cache_key'])

    @unittest.skipUnless((ROOT/'workspaces/mace_hybrid_20260916/pilot_v3/manifest.json').exists(),
                         'requires exact-byte technical retry manifest')
    def test_revised_manifest_reuses_every_input_byte(self):
        mp = ROOT/'workspaces/mace_hybrid_20260916/pilot_v3/manifest.json'
        newer = read_json(mp)
        self.assertEqual(mh.dry_run(mp)['tasks'], 12)
        for old, new in zip(self.m['tasks'], newer['tasks']):
            self.assertEqual({k:v for k,v in old.items() if k != 'cache_key'},
                             {k:v for k,v in new.items() if k != 'cache_key'})
            self.assertNotEqual(old['cache_key'], new['cache_key'])

    @unittest.skipUnless(importlib.util.find_spec('mace'), 'requires isolated MACE installation')
    def test_real_checkpoint_adapter_preserves_all_parameters_and_buffers(self):
        import torch
        from mace.calculators import mace_polar
        from mace_realspace_compat import configure, require_isolated
        torch.set_num_threads(1)
        c = mace_polar(model=str(verify(self.m['model']['checkpoint'])), device='cpu', default_dtype='float64')
        before = {k:v.clone() for k,v in c.models[0].state_dict().items()}
        metadata = configure(c)
        after = c.models[0].state_dict()
        self.assertEqual(before.keys(), after.keys())
        self.assertTrue(all(torch.equal(before[k], after[k]) for k in before))
        self.assertEqual(len(metadata['restored_dimension_metadata']), 4)
        self.assertEqual(metadata['features'], 'realspace')
        self.assertEqual(metadata['energy'], 'realspace')
        require_isolated(torch.tensor([[False, False, False]]))
        with self.assertRaises(ValueError):
            require_isolated(torch.tensor([[True, False, False]]))
        with self.assertRaises(ValueError):
            require_isolated(torch.tensor([[False, False, False]]), True)

    @unittest.skipUnless(importlib.util.find_spec('mace'), 'requires isolated MACE installation')
    def test_adapter_matches_original_kernels_on_real_computed_density(self):
        mp = ROOT/'workspaces/mace_hybrid_20260916/pilot_v3/manifest.json'
        if not mp.exists():
            self.skipTest('real interface-repair campaign unavailable')
        c = mh.collect(mp)
        if not c['core_partition']:
            self.skipTest('four real MACE core results have not completed')
        import torch
        from mace.calculators import mace_polar
        from mace_realspace_compat import configure
        torch.set_num_threads(1)
        calc = mace_polar(model=str(verify(self.m['model']['checkpoint'])), device='cpu', default_dtype='float64')
        configure(calc)
        features = calc.models[0].electric_potential_descriptor
        energy_module = calc.models[0].coulomb_energy
        for name in mh.TASK_IDS:
            coords = [a[1:] for a in xyz(verify(self.tasks[name]['xyz']))]
            positions = torch.tensor(coords, dtype=torch.float64, requires_grad=True)
            density = torch.tensor(np.load(verify(c['rows'][name]['density_coefficients'])), dtype=torch.float64)
            batch = torch.zeros(len(coords), dtype=torch.long)
            pbc = torch.zeros((1,3), dtype=torch.bool)
            kwargs = dict(k_vectors=torch.empty((0,3), dtype=torch.float64),
                          k_norm2=torch.empty(0, dtype=torch.float64), k_vector_batch=torch.empty(0, dtype=torch.long),
                          k0_mask=torch.empty(0, dtype=torch.float64), node_positions=positions,
                          batch=batch, volume=torch.ones(1, dtype=torch.float64), pbc=pbc)
            cache = features.precompute_geometry(**kwargs, force_pbc_evaluator=False)
            wrapped_features = features.forward_dynamic(cache=cache, source_feats=density[:,None,:], pbc=pbc)
            original_features = features.realspace_features(source_feats=density, node_positions=positions, batch=batch)[0]
            torch.testing.assert_close(wrapped_features, original_features, rtol=0, atol=0)
            wrapped_energy = energy_module(source_feats=density, **kwargs, force_pbc_evaluator=False)
            original_energy = energy_module.realspace_energy(source_feats=density, positions=positions, batch=batch)
            torch.testing.assert_close(wrapped_energy, original_energy, rtol=0, atol=0)
            grad_wrapped = torch.autograd.grad(wrapped_energy.sum(), positions, retain_graph=True)[0]
            grad_original = torch.autograd.grad(original_energy.sum(), positions)[0]
            torch.testing.assert_close(grad_wrapped, grad_original, rtol=0, atol=0)
        cp = c['core_partition']
        archive = {r['task_id']:r for r in read_json(verify(self.m['archived_endpoints']))['rows']}
        delta = lambda values, size, key: values[f'1h4i_qm{size}_Ca'][key]-values[f'1h4i_qm{size}_La'][key]
        expected = ((delta(archive,36,'energy_hartree')-delta(archive,33,'energy_hartree'))*mh.HA_TO_KCAL
                    -(delta(c['rows'],36,'energy_eV')-delta(c['rows'],33,'energy_eV'))*mh.EV_TO_KCAL)
        self.assertAlmostEqual(cp['hybrid_partition_shift_kcal_mol'], expected, places=8)


if __name__ == '__main__':
    unittest.main()
