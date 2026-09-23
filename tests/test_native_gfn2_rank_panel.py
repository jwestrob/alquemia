"""Pinned 28-source preparation and real-result algebra, without invented energies."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import native_gfn2_rank_panel as p

MANIFEST=ROOT/'workspaces/native_gfn2_rank_panel_20260923/run_v2/manifest.json'


class NativeRankPanel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=p.read_json(MANIFEST)

    def test_all_28_sources_and_336_exact_inputs(self):
        m=self.manifest
        self.assertEqual(p.validate(MANIFEST)['tasks'],336)
        self.assertEqual(len(m['cases']),28)
        self.assertEqual(sum(c['role']=='calibration' for c in m['cases']),25)
        self.assertEqual(m['inventory']['exact_pilot_reuse'],[])
        self.assertEqual(m['inventory']['SCF_limits'],{'explicit_500':284,'default_125':52})
        for t in m['tasks']:
            for field in ('input','xyz'):
                self.assertEqual(p.verify(t[field]).read_bytes(),p.verify(t['original_task'][field]).read_bytes())

    def test_archived_real_energy_and_pool_replay(self):
        cases={c['case_id']:c for c in self.manifest['cases']}
        for c in self.manifest['cases']:
            self.assertEqual(p.choose_rows(c['archived_matrix'],p.CANDIDATES),c['archived_pool'])
        for t in self.manifest['tasks']:
            receipt=p.read_json(p.verify(t['archived']['receipt']))
            self.assertTrue(receipt['normal_termination'])
            self.assertTrue(receipt['scf_converged'])
            self.assertEqual(receipt['parallelism']['nprocs'],8)
            self.assertEqual(p.read_json(p.verify(t['native']['result']))['energy_eV'],t['native']['energy_eV'])
            field='GFN2_vacuum_hartree' if t['medium']=='vacuum' else 'GFN2_ALPB_hartree'
            self.assertEqual(cases[t['case_id']]['archived_matrix'][t['metal']][t['candidate']]['components'][field],t['archived']['energy_hartree'])

    def test_corrupted_actual_manifest_is_rejected(self):
        for change in ('rank','membership','state'):
            m=copy.deepcopy(self.manifest)
            if change=='rank':m['execution_resources']['mpi_ranks']=4
            elif change=='membership':m['tasks']=m['tasks'][:-1];m['all_tasks']=m['tasks']
            else:m['tasks'][0]['charge']+=1;m['all_tasks']=m['tasks']
            with tempfile.TemporaryDirectory() as d:
                path=Path(d)/'explicitly_corrupted_actual_manifest.json';p.write_new(path,m)
                with self.assertRaises(p.InvalidArtifact):p.validate(path)

    def test_actual_incomplete_collection_retains_denominator(self):
        # This file was collected from the real empty execution directories before launch.
        path=MANIFEST.parent/'INCOMPLETE_COLLECTION.json'
        if not path.exists():self.skipTest('actual prelaunch collection has not been written')
        c=p.read_json(path);self.assertEqual(c['denominator'],336);self.assertEqual(c['complete'],0)
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/'comparison.json';p.compare(path,output);r=p.read_json(output)
            self.assertEqual(r['case_denominator'],28);self.assertEqual(r['complete'],0)
            self.assertFalse(r['all_numerical_gates_pass'])
            for case in r['cases']:
                self.assertEqual(case['archived8_pool']['status'],'available')
                self.assertEqual(case['rank1_pool']['status'],'unavailable')
                self.assertIsNone(case['variants']['operational']['rank1_R'])

    def test_actual_new_scientific_results_when_available(self):
        path=MANIFEST.parent/'COMPARISON.json'
        if not path.exists():self.skipTest('336-cell molecular execution not completed')
        r=p.read_json(path);self.assertEqual(r['case_denominator'],28);self.assertEqual(r['cell_denominator'],336)
        self.assertFalse(r['native_gradient_qualified']);self.assertFalse(r['reference_refitted'])
        for c in r['cases']:
            for v in c['variants'].values():
                if v['rank1_R'] is not None:
                    self.assertAlmostEqual(v['rank1_R']-v['archived8_R'],v['delta_R_kcal_mol'],places=12)
                    self.assertEqual(v['passes_pool_R_gate'],abs(v['delta_R_kcal_mol'])<=.2)


if __name__=='__main__':unittest.main()
