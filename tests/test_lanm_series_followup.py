"""Real prepared-source and archived-output checks; no new molecular evaluations."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import lanm_series_followup as lanm
from affordable_common import InvalidArtifact, read_json, verify, xyz

MANIFEST = ROOT/'workspaces/lanm_series_followup_20260923/prepared_v2/manifest.json'
SOURCES = ROOT/'workspaces/lanm_series_followup_20260923/sources_v2/sources.json'
UNRUN_MANIFEST = ROOT/'workspaces/lanm_series_followup_20260923/prepared_v1/manifest.json'


class RealLanMPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = read_json(MANIFEST)
        cls.sources = read_json(SOURCES)
        cls.params = read_json(verify(cls.m['parameter_export']))

    def test_full_real_inventory_and_contained_native_inputs(self):
        r = lanm.validate(MANIFEST)
        self.assertEqual((r['supported_sites'], r['prepared_MACE'], r['initial_native_tasks']), (6,12,24))
        self.assertEqual(r['new_endpoint_calls'], 0)

    def test_distinct_physical_and_effective_states(self):
        for e in self.m['endpoints']:
            st = lanm.state(xyz(verify(e['xyz'])), e['charge'], e['metal'], self.params)
            self.assertEqual(st['physical_multiplicity'], 6 if e['metal']=='Dy' else 1)
            self.assertEqual(st['effective_GFN_multiplicity'], 1)
            self.assertEqual(st['GFN_valence_electron_count'], 150 if e['case_id'].startswith('Hans') else 136)
        # Deliberately corrupt one real fixture's charge: never a scientific state.
        e = next(e for e in self.m['endpoints'] if e['metal']=='Dy')
        with self.assertRaisesRegex(InvalidArtifact, 'parity'):
            lanm.state(xyz(verify(e['xyz'])), e['charge']+1, e['metal'], self.params)

    def test_real_water_identity_and_exact_source_heavy_atoms(self):
        expected = {'Mex_EF1':[326,345], 'Mex_EF2':[322,328], 'Mex_EF3':[320,331]}
        for case in self.sources['cases']:
            repair = read_json(verify(case['repair_manifest']))
            self.assertEqual([r['resnum'] for r in repair['explicit_water_inventory']], expected.get(case['case_id'], []))
            self.assertEqual(repair['heavy_coordinate_max_displacement_A'], 0.)
        original = ROOT/'workspaces/benchmark_set_20260915/evidence/lanthanide/8FNS.cif'
        actual = SOURCES.parent/'mex_protonated.pdb'
        self.assertEqual(lanm.heavy_atoms(original), lanm.heavy_atoms(actual))

    def test_altloc_selection_uses_existing_policy_without_moving_atoms(self):
        original = ROOT/'workspaces/benchmark_set_20260915/evidence/lanthanide/8FNS.cif'
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'selected.pdb'
            r = lanm.select_source_conformers(original,p)
            self.assertEqual(lanm.heavy_atoms(original),lanm.heavy_atoms(p))
            self.assertIn(52,[x['resnum'] for x in r['selections']])
        failed = read_json(ROOT/'workspaces/lanm_series_followup_20260923/sources_v1/sources.json')
        self.assertEqual(failed['supported'],3)
        self.assertGreater(failed['source_heavy_check']['max_displacement_A'],1.5)

    def test_private_Nd_selector_does_not_change_production(self):
        import carve_generic
        old = carve_generic.METAL_ELEMENTS
        adapter = lanm.private_carver()
        self.assertIn('ND', adapter.METAL_ELEMENTS)
        self.assertEqual(carve_generic.METAL_ELEMENTS,old)
        self.assertNotIn('ND',old)

    def test_unrun_endpoints_remain_null_all_denominators_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            result = lanm.collect(UNRUN_MANIFEST,Path(d)/'collection.json')
            self.assertEqual(result['complete_endpoints'],0)
            self.assertEqual(len(result['endpoints']),12)
            self.assertEqual([x['site'] for x in result['ordered_vector']],['EF1','EF2','EF3'])
            self.assertTrue(all(x['D_composite_kcal_mol'] is None for x in result['ordered_vector']))
            continuation = lanm.prepare_continuation(UNRUN_MANIFEST,ROOT/'diagnostics/lanm_series_followup_20260923/PLAN.md',Path(d)/'continuation')
            self.assertEqual(continuation['prepared'],0)
            self.assertEqual(len(continuation['unavailable']),24)

    def test_archived_native_restart_recipe_state_and_balanced_algebra(self):
        from compact_solvation import diagnostics
        archive = read_json(ROOT/'workspaces/native_xtb_restart_20260922/recovery_v1/collection_1210185.json')
        real = archive['rows'][0]
        m = read_json(verify(real['actual']['manifest']))
        task = next(t for t in m['tasks'] if t['task_id']==real['actual']['task_id'])
        audit = diagnostics(real['actual'],task)
        self.assertEqual(audit['metal_reference_occupation'],[1.,1.,1.])
        self.assertIn('INITIAL GUESS: XTBRESTART',verify(real['actual']['output']).read_text())
        self.assertEqual(lanm.recipe(-2,task['medium'],True),verify(task['input']).read_text())
        energies = {(r['geometry'],r['medium']): r['actual']['energy_hartree']
                    for r in archive['rows'] if r['seed_kind']=='self'}
        values = [energies[g,s] for g in ('old_adaptive','new_template') for s in ('vacuum','alpb')]
        # Pure exchange algebra using four real PQQ energies, not invented LanM/Dy data.
        expected = (values[1]-values[0])-(values[3]-values[2])
        self.assertEqual(lanm.balanced_difference(*values),expected)
        self.assertEqual(lanm.balanced_difference(*values[2:],*values[:2]),-expected)
        self.assertIsNone(lanm.balanced_difference(None,*values[1:]))


if __name__ == '__main__':
    unittest.main()
