"""Real prepared fold sources and actual archived energies; no solver mocks."""
import copy
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import accommodation_proposals as engine
import accommodation_fold_proposals as folds
from affordable_common import InvalidArtifact, read_json, record, verify, xyz


class FoldProposalPreflight(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = ROOT / 'workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/manifest.json'
        if not cls.path.exists(): raise unittest.SkipTest('real frozen fold manifest not prepared')
        cls.m = read_json(cls.path)

    def test_exact_primary_population_and_one_actual_missing_origin(self):
        result = engine.validate(self.path)
        self.assertEqual(result['starts'], 416)
        self.assertEqual(result['q0_available'], 415)
        self.assertEqual(len(self.m['unavailable_cases']), 17)
        for t in self.m['tasks']:
            if t['q0_status'] != 'available': continue
            for low in t['q0']['low'].values():
                self.assertEqual(xyz(verify(t['xyz'])), xyz(verify(low['task']['xyz'])))

    def test_corrupted_real_scope_rejected(self):
        corrupted = copy.deepcopy(self.m); corrupted['maximum_new_GFN2_singlepoints'] += 2
        with self.assertRaisesRegex(InvalidArtifact, 'finite fold scope'): folds.validate_scope(corrupted)
        corrupted = copy.deepcopy(self.m); corrupted['cases'].pop()
        with self.assertRaisesRegex(InvalidArtifact, 'population differs'): folds.validate_scope(corrupted)

    def test_real_archived_composite_algebra_replays_every_available_pair(self):
        old = read_json(verify(self.m['fold_comparison']))
        index = {r['case_id']: r for r in old['rows'] if r['representation'] == 'context'}
        n = 0
        for c in self.m['cases']:
            pair = [next(t for t in self.m['tasks'] if t['case_id'] == c['case_id'] and t['metal'] == z) for z in ('Ca', 'La')]
            if any(t['q0_status'] != 'available' for t in pair): continue
            score = engine.score(*[t['q0']['components'] for t in pair])
            self.assertEqual(score['composite_R_model_kcal_mol'], index[c['case_id']]['composite_R_model_kcal_mol'])
            n += 1
        self.assertEqual(n, 207)

    def test_actual_prelaunch_missing_proposals_preserve_baselines_and_denominators(self):
        if any(self.path.parent.glob('proposals/*/result.json')):
            self.skipTest('real run started; original prelaunch receipt retained')
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'selected.json'; result = engine.collect(self.path, out)
            self.assertEqual((result['case_denominator'], result['endpoint_denominator'], result['prepared_endpoint_denominator']), (225, 450, 416))
            self.assertEqual(result['available_cases'], 0)
            dft = ROOT / 'workspaces/accommodation_fold_DFT_20260920/partial_comparison_v2/dft_collection.json'
            counts = folds.compare(self.path, out, dft, Path(tmp) / 'comparison.json')
            self.assertEqual(counts['single']['all225']['proposal_new']['unavailable'], 225)
            self.assertGreater(counts['single']['all225']['context_composite']['correct'], 0)
            self.assertEqual(counts['all100_triples']['proposal_new']['unavailable'], 100)
            self.assertEqual(counts['pools']['La4']['proposal_new']['unavailable'], 25)
            self.assertEqual(counts['all100_triples']['context_composite']['correct'], 94)

    def test_original_completed_thirty_case_selection_regression(self):
        mp = ROOT / 'workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json'
        old = read_json(mp.parent / 'final_1204171.json')
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'replay.json'; engine.collect(mp, out); new = read_json(out)
            self.assertEqual(new['selection_counts'], old['selection_counts'])
            # Collection adds source-group metadata; every original field must replay exactly.
            for before, after in zip(old['cases'], new['cases']):
                self.assertEqual({k: after[k] for k in before}, before)
                self.assertLessEqual(set(after) - set(before), {'root_case_id', 'biological_group',
                                                              'source_conditioning_metal', 'canonical_coordinate_match'})
            self.assertEqual(new['counts'], old['counts'])

    def test_actual_proposal_input_staging_in_four_disjoint_shards(self):
        original = read_json(ROOT / 'workspaces/accommodation_nonlinear_20260920/proposals_v1/proposal_GFN2/manifest.json')
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp); tasks = []
            for i, source in enumerate(original['all_tasks'][:12]):
                task = dict(source); folder = dest / ('shard_' + str(i % 4)) / 'tasks' / task['task_id']
                folder.mkdir(parents=True)
                for key, name in [('input', 'endpoint.inp'), ('xyz', 'core.xyz')]:
                    path = folder / name; shutil.copyfile(verify(task[key]), path); task[key] = record(path)
                task['output_path'] = str(folder / 'endpoint.out'); tasks.append(task)
            master = engine.write_low_shards(dest, {**original, 'all_tasks': tasks, 'tasks': tasks, 'reused': {}}, 4)
            self.assertEqual(master['tasks'], [])
            self.assertEqual(len(master['task_manifests']), 12)
            seen = set()
            for pin in master['shards']:
                path = verify(pin); sm = read_json(path)
                self.assertEqual(read_json(path.parent / 'PREFLIGHT.json')['status'], 'dry_run_pass')
                for task in sm['tasks']:
                    self.assertNotIn(task['task_id'], seen); seen.add(task['task_id'])
                    self.assertTrue(Path(task['output_path']).is_relative_to(path.parent))
                    self.assertEqual(master['task_manifests'][task['task_id']], pin)

    def test_receipt_indirection_uses_actual_archived_shard(self):
        old = read_json(verify(self.m['fold_comparison'])); tasks = {}
        for pin in old['collections']:
            mp = read_json(verify(pin))['manifest']; sm = read_json(verify(mp))
            task = next(t for t in sm['all_tasks'] if t['task_id'] not in sm['reused'] and engine.completed(verify(mp), t['task_id']) is not None)
            tasks[task['task_id']] = mp
        master = {'task_manifests': tasks, 'reused': {}}
        for tid, pin in tasks.items():
            receipt = engine.low_completed(master, None, tid)
            self.assertEqual(receipt['manifest'], pin)
            self.assertEqual(receipt['task_id'], tid)


if __name__ == '__main__': unittest.main()
