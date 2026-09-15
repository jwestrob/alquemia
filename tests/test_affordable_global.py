"""Saved scientific artifacts only; these are preparation/algebra, not xTB runs."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, xyz
from affordable_global import load_primary, direct_decomposition
from affordable_global_collect import mulliken
from affordable_global_execute import layout
from run_orca_task_manifest import load_manifest_tasks, _verify_prepared_task

RESULT = ROOT/'workspaces/affordable_challenger_20260915/solver_completion/solver_result.json'
PREPARED = ROOT/'workspaces/global_representation_20260915/prepared_v1/global_manifest.json'


@unittest.skipUnless(RESULT.exists() and PREPARED.exists(), 'pinned local scientific artifacts unavailable')
class GlobalPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.states, cls.components, _ = load_primary(RESULT)
        cls.prepared = read_json(PREPARED)

    def test_saved_direct_contrast_reproduces_both_partitions(self):
        for n in (33, 36):
            left, right = f'1h4i_qm{n}_La', f'1h4i_qm{n}_Ca'
            d = direct_decomposition(self.states[left], self.states[right])
            expected = self.components[right]['direct_core_environment_kcal_mol']-self.components[left]['direct_core_environment_kcal_mol']
            self.assertAlmostEqual(d['direct_score_kcal_mol'], expected, places=9)
            self.assertAlmostEqual(sum(d['per_environment_residue_kcal_mol'].values()), expected, places=9)
            self.assertAlmostEqual(sum(a['direct_score_kcal_mol'] for a in d['per_core_atom']), expected, places=9)

    def test_corrupted_pair_rejected(self):
        # Explicitly corrupted copy of a real state, never scientific evidence.
        ca = deepcopy(self.states['1h4i_qm33_Ca'])
        ca['environment_atoms'][0]['charge_e'] += 1
        with self.assertRaises(InvalidArtifact):
            direct_decomposition(self.states['1h4i_qm33_La'], ca)

    def test_complete_physical_pair_and_source_coordinates(self):
        m = self.prepared
        endpoints = [xyz(verify(t['xyz'])) for t in m['tasks']]
        self.assertEqual(len(endpoints[0]), 9141)
        self.assertEqual([a[1:] for a in endpoints[0]], [a[1:] for a in endpoints[1]])
        differences = [(a[0], b[0]) for a,b in zip(*endpoints) if a[0] != b[0]]
        self.assertEqual(differences, [('La','Ca')])
        self.assertEqual([t['charge'] for t in m['tasks']], [-8,-9])
        self.assertTrue(all(t['all_electron_count_for_parity_only']%2 == 0 for t in m['tasks']))
        self.assertLess(m['source_coordinate_max_difference_A'], 1e-9)
        self.assertEqual(m['synthetic_caps'], [])
        self.assertEqual(m['cut_bonds'], [])
        self.assertIsNone(m['score'])
        self.assertIsNone(m['settings']['reference'])

    def test_asp303_covalent_connection_and_pqq_graph(self):
        m = self.prepared
        idx = {a['id']: i for i,a in enumerate(m['physical_atoms'])}
        bonds = {tuple(b) for b in m['bonds_zero_based']}
        self.assertIn(tuple(sorted((idx['A/303/ /CA'],idx['A/303/ /CB']))), bonds)
        self.assertIn(tuple(sorted((idx['A/303/ /C'],idx['A/304/ /N']))), bonds)
        self.assertEqual(len(bonds), len(m['bonds_zero_based']))
        self.assertEqual(len(m['pqq_source_name_to_global_index']), 27)

    def test_existing_runner_accepts_hashed_manifest(self):
        _, tasks = load_manifest_tasks(PREPARED)
        self.assertEqual(len(tasks), 2)
        for task in tasks:
            _verify_prepared_task(task)
        self.assertEqual(verify(self.prepared['implementation']), ROOT/'scripts/affordable_global.py')

    def test_population_parser_on_six_archived_quantum_outputs(self):
        pilot = read_json(ROOT/'workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json')
        for task in pilot['tasks']:
            population = mulliken(task['output_path'], verify(task['xyz']), task['charge'])
            self.assertAlmostEqual(population['printed_sum_e'], task['charge'], places=6)
            self.assertEqual(len(population['charge_e']), len(xyz(verify(task['xyz']))))

    def test_layout_fits_observed_cluster_hardware(self):
        # Real sinfo hardware capacities; no chemical calculations are mocked.
        for cpus, memory_MiB, expected in [(64,773914,4),(112,1546754,8),(344,8256990,16)]:
            plan = layout(cpus, memory_MiB*1024**2)
            self.assertEqual(plan['nprocs_per_endpoint'], expected)
            self.assertLessEqual(2*expected*64*1024**3, memory_MiB*1024**2*.75)
            self.assertLessEqual(2*expected, cpus)


if __name__ == '__main__':
    unittest.main()
