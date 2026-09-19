"""Real pinned water arrangements; these checks run no scientific executable."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import hydration_occupancy as h
from hydration_network import write_xyz


class OccupancyFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp = ROOT/'workspaces/hydration_occupancy_20260918/states_v1/manifest.json'
        if not mp.exists():
            raise unittest.SkipTest('real enumerated occupancy preparations absent')
        cls.manifest = h.read_json(mp)
        cls.preps = {p['case']: p for p in
                     (h.read_json(h.verify(pin)) for pin in cls.manifest['preparations'])}

    def test_real_subsets_preserve_source_graph_and_neutral_water_accounting(self):
        for case, p in self.preps.items():
            n = 2 if case == '1F6S' else 3
            rows = [t for t in self.manifest['tasks'] if t['case'] == case and t['seed'] == 'source']
            self.assertEqual({t['pattern'] for t in rows}, h.pattern_set(n)-{'1'*n})
            for t in rows:
                atoms = h.xyz(h.verify(t['xyz']))
                self.assertEqual(len(atoms), p['atom_count']-3*(n-t['pattern'].count('1')))
                self.assertEqual(t['charge'], -1 if t['metal'] == 'La' else -2)
                partner = next(s for s in rows if s['pattern'] == t['pattern'] and s['metal'] != t['metal'])
                other = h.xyz(h.verify(partner['xyz']))
                self.assertEqual(atoms[0][1:], other[0][1:])
                self.assertEqual(atoms[1:], other[1:])

    def test_full_cached_endpoints_match_identical_physical_context(self):
        path = ROOT/'workspaces/hydration_network_20260918/proposal_dft_v1/manifest.json'
        if not path.exists():
            self.skipTest('completed DFT proposal fixture missing')
        for task in h.read_json(path)['tasks']:
            n = 2 if task['case'] == '1F6S' else 3
            indices = h.check_state(dict(task, pattern='1'*n), self.preps[task['case']])
            self.assertEqual(len(indices), 70 if n == 2 else 76)

    def test_corrupted_real_coordinate_or_method_cannot_reuse_cache(self):
        t = copy.deepcopy(self.manifest['tasks'][0]); p = self.preps[t['case']]
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory); ip = directory/'endpoint.inp'
            ip.write_text(h.METHOD+f"\n* xyzfile {t['charge']} 1 core.xyz\n")
            t['input'] = h.record(ip)
            h.check_state(t, p)
            damaged = h.xyz(h.verify(t['xyz']))
            damaged[1] = (damaged[1][0], damaged[1][1]+.1, *damaged[1][2:])
            xp = directory/'corrupted_real.xyz'; write_xyz(xp, damaged, 'explicit corrupted real fixture')
            with self.assertRaises(h.InvalidArtifact):
                h.check_state(dict(t, xyz=h.record(xp)), p)
            ip.write_text(ip.read_text().replace('TightSCF', 'LooseSCF'))
            with self.assertRaises(h.InvalidArtifact):
                h.check_state(dict(t, input=h.record(ip)), p)

    def test_omitted_real_occupancy_task_fails_validation(self):
        from hydration_proposal_opt import validate
        path = ROOT/'workspaces/hydration_occupancy_20260918/mace_v1/manifest.json'
        if not path.exists():
            self.skipTest('real proposal fixture unavailable')
        m = h.read_json(path); m['tasks'] = m['tasks'][1:]
        with tempfile.TemporaryDirectory() as directory:
            damaged = Path(directory)/'missing_real_task.json'; h.write_new(damaged, m)
            with self.assertRaises(h.InvalidArtifact):
                validate(damaged)

    def test_partial_real_collection_keeps_missing_states_unavailable(self):
        path = ROOT/'workspaces/hydration_occupancy_20260918/dft_v1/preexecution_collection.json'
        reference = ROOT/'workspaces/hydration_network_20260918/water_reference_v2/reference_1201831.json'
        if not path.exists() or not reference.exists():
            self.skipTest('actual partially populated experiment fixture unavailable')
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)/'analysis.json'
            result = h.analyze(path, reference, out)
            self.assertEqual(result['status'], 'incomplete')
            data = h.read_json(out)
            self.assertIsNone(data['occupancy_probabilities'])
            for case in data['cases']:
                self.assertFalse(case['complete_occupancy_table'])
                self.assertEqual(case['water_addition_edges'], [])
                self.assertEqual(len(case['fixed_count_electronic_contrasts']), 1)
                self.assertFalse(any('exchange_ledger' in s for s in case['states']))

    def test_real_saved_gradient_projects_physical_water_translation(self):
        import numpy as np
        path = ROOT/'workspaces/hydration_network_20260918/proposal_dft_v1/collection_numeric_parser_v2.json'
        if not path.exists():
            self.skipTest('native gradient fixture unavailable')
        c = h.read_json(path); m = h.read_json(h.verify(c['manifest']))
        t = m['tasks'][0]; row = next(r for r in c['rows'] if r['task_id'] == t['task_id'])
        g = h.read_engrad(h.verify(row['gradient']))['gradient_Ha_per_bohr']
        coords = np.array([a[1:] for a in h.xyz(h.verify(t['xyz']))])
        water = next(w for w in t['groups'] if w['role'] == 'variable')
        displacement = np.array([.001, -.002, .003])
        changed = coords.copy(); changed[water['indices']] += displacement
        self.assertAlmostEqual(float(g[water['indices']].sum(axis=0)@displacement),
                               float(np.sum(g*(changed-coords))), places=12)

    def test_tampered_parsed_energy_is_rejected_against_real_output(self):
        path = ROOT/'workspaces/hydration_occupancy_20260918/dft_v1/preexecution_collection.json'
        reference = ROOT/'workspaces/hydration_network_20260918/water_reference_v2/reference_1201831.json'
        if not path.exists():
            self.skipTest('real partially populated experiment fixture unavailable')
        c = h.read_json(path)
        next(r for r in c['rows'] if r['status']=='complete')['result']['energy_hartree'] += 1
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'corrupted_real_collection.json'; h.write_new(p, c)
            with self.assertRaises(h.InvalidArtifact):
                h.analyze(p, reference, Path(directory)/'result.json')

    def test_experimental_site_alignment_is_a_proper_reversible_rigid_transform(self):
        import numpy as np
        from hydration_site_proposals import protein_atoms, fit
        a, b = (protein_atoms(self.preps[c]) for c in ('1F6S', '6IP9'))
        keys = sorted(a.keys() & b.keys())
        x, y = np.array([a[k] for k in keys]), np.array([b[k] for k in keys])
        r, t, rmsd = fit(x, y); inverse, shift, _ = fit(y, x)
        self.assertAlmostEqual(float(np.linalg.det(r)), 1, places=12)
        np.testing.assert_allclose((x@r+t)@inverse+shift, x, atol=1e-12)
        translated = np.array([9., -5., 7.])
        rr, tt, rmsd2 = fit(x+translated, y)
        np.testing.assert_allclose((x+translated)@rr+tt, x@r+t, atol=1e-12)
        self.assertAlmostEqual(rmsd, rmsd2, places=12)

    def test_completed_real_addition_cycles_and_bulk_reference_cancellation(self):
        path = ROOT/'diagnostics/hydration_occupancy_20260918/RESULT.json'
        if not path.exists():
            self.skipTest('actual completed occupancy integration output unavailable')
        result = h.read_json(path)
        self.assertEqual(result['status'], 'complete')
        self.assertIsNone(result['occupancy_probabilities'])
        self.assertEqual(sum(len(c['states']) for c in result['cases']), 12)
        for case in result['cases']:
            states = {s['pattern']: s for s in case['states']}
            edges = {(e['from'], e['to']): e for e in case['water_addition_edges']}
            n = len(case['variable_water_order'])
            path_a = ['1'*i+'0'*(n-i) for i in range(n+1)]
            path_b = ['0'*(n-i)+'1'*i for i in range(n+1)]
            totals = [sum(edges[a,b]['delta_R_water_addition_kcal_mol'] for a,b in zip(p,p[1:]))
                      for p in (path_a, path_b)]
            self.assertAlmostEqual(totals[0], totals[1], places=8)
            direct = (states['1'*n]['R_hartree']-states['0'*n]['R_hartree'])*h.HA_TO_KCAL
            self.assertAlmostEqual(totals[0], direct, places=8)
            for edge in edges.values():
                ledger = edge['exchange_ledger']
                required_difference = (ledger['Ca']['missing_bound_contribution_for_neutral_exchange_kcal_mol']
                                       -ledger['La']['missing_bound_contribution_for_neutral_exchange_kcal_mol'])
                self.assertAlmostEqual(required_difference, -edge['delta_R_water_addition_kcal_mol'], places=7)


if __name__ == '__main__':
    unittest.main()
