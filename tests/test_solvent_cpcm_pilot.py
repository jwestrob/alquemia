"""Actual archived/native fixtures; unavailable CPCM is never fabricated."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import solvent_cpcm_pilot as pilot
import solvent_cpcm_matched as matched
from compact_solvation import completed
from affordable_common import InvalidArtifact, HA_TO_KCAL, read_json, write_new


class CPCMPilot(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = ROOT/'workspaces/solvent_cpcm_20260922/pilot_v1/manifest.json'
        if not cls.manifest.exists(): raise unittest.SkipTest('real prepared pilot unavailable')
        cls.m = read_json(cls.manifest)

    def test_actual_finite_manifest(self):
        self.assertEqual(pilot.validate(self.manifest)['tasks'], 8)
        self.assertEqual(self.m['new_DFT_calls'], 0)
        self.assertEqual(self.m['new_MACE_calls'], 0)

    def test_real_origin_pair_and_algebra(self):
        for c in self.m['cases']:
            e = c['endpoints']
            self.assertEqual(e['La']['charge'], e['Ca']['charge']+1)
            totals = {}
            for z in ('Ca', 'La'):
                s = e[z]
                totals[z] = s['MACE_energy_eV']*pilot.EV_TO_KCAL + (
                    s['GFN2']['alpb']['endpoint']['energy_hartree']-
                    s['GFN2']['vacuum']['endpoint']['energy_hartree'])*HA_TO_KCAL
            self.assertAlmostEqual(totals['Ca']-totals['La'], c['old_result']['R0']['composite_R_model_kcal_mol'], places=8)

    def test_archived_alpb_cannot_pass_as_cpcm(self):
        c = self.m['cases'][0]['endpoints']['Ca']
        with self.assertRaisesRegex(InvalidArtifact, 'actual CPCM cavity not confirmed'):
            pilot.cpcm_audit(c['GFN2']['alpb']['endpoint'], c['GFN2']['alpb']['task'], c)

    def test_explicitly_corrupted_real_method_rejected(self):
        damaged = copy.deepcopy(self.m); damaged['settings']['epsilon'] = 4.0
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'CORRUPTED_real_manifest.json'; write_new(p, damaged)
            with self.assertRaisesRegex(InvalidArtifact, 'frozen solvent protocol differs'):
                pilot.validate(p)

    def test_actual_unrun_collection_keeps_missing_not_zero(self):
        if any(Path(t['output_path']).exists() for t in self.m['tasks']):
            self.skipTest('pilot has started; retained prelaunch incomplete receipt tests this state')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'collection.json'; pilot.collect(self.manifest, p); c=read_json(p)
            self.assertEqual(c['complete_endpoints'], 0)
            self.assertEqual(c['available_pairs'], 0)
            self.assertTrue(all(r['CPCM_R'] is None and r['new_decision'] is None for r in c['rows']))
            self.assertTrue(all(r['released_ALPB_R'] is not None for r in c['rows']))

    def test_actual_scientific_cpcm_result_when_present(self):
        paths=sorted(self.manifest.parent.glob('collection_*.json'))
        if not paths: self.skipTest('complete primary-phase collection is not yet available')
        c=read_json(paths[-1])
        for r in c['rows']:
            self.assertIsNone(r['new_decision'])
            if r['status']=='available':
                self.assertEqual(r['CPCM_R'], r['components']['Ca']['composite_energy_kcal_mol']-r['components']['La']['composite_energy_kcal_mol'])
            else: self.assertIsNone(r['CPCM_R'])

    def test_actual_ordinary_cpcm_boundary(self):
        t=next(t for t in self.m['tasks'] if t['case_id']=='1H4I' and t['metal']=='Ca')
        e=completed(self.manifest,t['task_id'])
        if not e: self.skipTest('actual first CPCM endpoint unavailable')
        s=next(c for c in self.m['cases'] if c['case_id']=='1H4I')['endpoints']['Ca']
        audit=matched.audit_ordinary(e,t,s)
        self.assertEqual(audit['effective_solver'],'ordinary_ORCA_SCF')
        self.assertLess(audit['CPCM']['printed_dielectric_energy_hartree'],0)
        self.assertEqual(audit['CPCM']['radii_A']['Ca'],2.772)
        with self.assertRaisesRegex(InvalidArtifact,'effective mixer changed'):
            pilot.cpcm_audit(e,t,s)

    def test_actual_matched_vacuum_manifest_and_restart(self):
        mp=self.manifest.parent.parent/'matched_vacuum_v1/manifest.json'
        if not mp.exists():self.skipTest('actual matched vacuum preparation unavailable')
        self.assertEqual(matched.validate(mp)['fresh_CPCM'],0)
        m=read_json(mp)
        t=next(t for t in m['tasks'] if t['case_id']=='4MAE' and t['metal']=='Ca')
        e=completed(mp,t['task_id'])
        if not e:self.skipTest('actual ordinary vacuum integration not complete')
        s=next(c for c in self.m['cases'] if c['case_id']=='4MAE')['endpoints']['Ca']
        audit=matched.audit_ordinary(e,t,s,restart=True)
        self.assertEqual(audit['restart']['restart_status'],'orbital_restart_confirmed')
        self.assertEqual(audit['TolE_hartree'],1e-6)


if __name__=='__main__': unittest.main()
