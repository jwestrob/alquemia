"""Real pinned fragment inputs and native output accounting; no invented energies."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, write_new
from interaction_decomposition import (prepare, prepare_order_retry, parse_native_terms,
                                        parse_eda, inspect_fragment_input, collect, report)

WORK = ROOT / 'workspaces/density_embedding_20260916'


class InteractionArtifacts(unittest.TestCase):
    def require(self, path):
        if not path.exists(): self.skipTest('real artifact unavailable: ' + str(path))
        return path

    def test_retry_preserves_every_token_except_block_order(self):
        source = self.require(WORK / 'eda_v1/manifest.json')
        old = read_json(source)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'retry'
            result = prepare_order_retry(source, out)
            self.assertEqual(result['status'], 'dry_run_pass')
            new = read_json(out / 'manifest.json')
            for a,b in zip(old['tasks'],new['tasks']):
                before, after = verify(a['input']).read_text(), verify(b['input']).read_text()
                self.assertEqual(sorted(before.splitlines()), sorted(after.splitlines()))
                self.assertLess(after.index('* xyzfile'), after.index('%Frag'))
                self.assertEqual(a['xyz']['sha256'], b['xyz']['sha256'])
                self.assertEqual(a['fragment1_indices'], b['fragment1_indices'])
                self.assertEqual(a['fragment2_indices'], b['fragment2_indices'])
            with self.assertRaises(InvalidArtifact): prepare_order_retry(source, out)

    def test_reproducible_preparation_matches_actual_submitted_science(self):
        frozen = read_json(self.require(WORK / 'eda_order_retry_v1/manifest.json'))
        sm = self.require(ROOT / 'workspaces/global_electrostatic_20260916/states_v1/states_manifest.json')
        states = read_json(sm)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'prepared'
            result = prepare(sm, verify(states['endpoint_manifest']),
                             verify(frozen['agreement']), verify(frozen['parent_authorization']), out)
            self.assertEqual(result['status'], 'dry_run_pass')
            for a,b in zip(frozen['tasks'], read_json(out / 'manifest.json')['tasks']):
                self.assertEqual(a['input']['sha256'], b['input']['sha256'])
                self.assertEqual(a['xyz']['sha256'], b['xyz']['sha256'])

    def test_archived_composite_energies_close_and_failed_output_rejected(self):
        frozen = read_json(self.require(WORK / 'eda_v1/manifest.json'))
        for task in frozen['tasks']:
            for name in ('source_adduct_task', 'source_fragment1_task'):
                terms = parse_native_terms(task[name]['output_path'])
                self.assertLess(abs(terms['native_closure_kcal_mol']), .01)
            with self.assertRaises(InvalidArtifact): parse_eda(task['output_path'])

    def test_native_generated_ghosts_are_distinct_from_physical_atoms(self):
        m = read_json(self.require(WORK / 'eda_order_retry_v1/manifest.json'))
        for task in m['tasks']:
            for number, count in ((1, 47), (2, 7)):
                p = self.require(Path(task['output_path']).parent / f'endpoint.runtime_frag{number}.inp')
                info = inspect_fragment_input(p, task, number)
                self.assertEqual(info['physical_atom_count'], count)
                self.assertEqual(info['ghost_atom_count'], 54-count)
                self.assertLess(info['max_shape_serialization_error_A'], 2e-6)
                # Explicit corruption of actual native ghost/nucleus mapping.
                with tempfile.TemporaryDirectory() as tmp:
                    bad = Path(tmp) / 'corrupted.inp'
                    bad.write_text(p.read_text().replace(':(', '(', 1))
                    with self.assertRaises(InvalidArtifact): inspect_fragment_input(bad, task, number)

    def test_actual_failed_execution_has_no_decomposition_or_fallback(self):
        result = collect(self.require(WORK / 'eda_v1/manifest.json'))
        self.assertEqual(result['status'], 'failed_native_execution')
        self.assertIsNone(result['comparison'])
        self.assertIsNone(result['affinity_score'])
        self.assertFalse(result['attribution_to_original_partition_supported'])
        self.assertEqual(len(result['rows']), 2)
        for row in result['rows']:
            self.assertEqual(row['execution_returncode'], 55)
            self.assertEqual(row['status'], 'unavailable')
            self.assertEqual(sum(x['converged_scf_markers'] for x in row['scf_inventory']), 0)

    def test_failed_real_attempt_report_does_not_invent_components(self):
        result = collect(self.require(WORK / 'eda_v1/manifest.json'))
        with tempfile.TemporaryDirectory() as tmp:
            saved, rendered = Path(tmp) / 'result.json', Path(tmp) / 'report.md'
            write_new(saved, result)
            report(saved, rendered)
            text = rendered.read_text()
            self.assertIn('failed_native_execution', text)
            self.assertIn('No complete two-metal decomposition is available', text)
            self.assertNotIn('| orbital |', text)
            with self.assertRaises(FileExistsError): report(saved, rendered)


if __name__ == '__main__': unittest.main()
