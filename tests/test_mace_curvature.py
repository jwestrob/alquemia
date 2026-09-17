"""Real archived physical displacements; no fabricated scientific output."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL,InvalidArtifact,read_json,verify,xyz
from mace_curvature import dft_source,differences,validate,collect_curvature
DFT=ROOT/'workspaces/ggr_mechanism_20260915/report_v1/collection_c.json'
MANIFEST=ROOT/'workspaces/mace_curvature_20260916/mace_v1/medium/manifest.json'


@unittest.skipUnless(DFT.exists() and MANIFEST.exists(),'real executed DFT and prepared MACE inputs required')
class CurvatureTests(unittest.TestCase):
    def test_reuses_exact_paired_DFT_inputs(self):
        self.assertEqual(validate(MANIFEST)['tasks'],20)
        d,m=dft_source(DFT,True); new=read_json(MANIFEST)
        for old,t in zip(m['tasks'],new['tasks']):self.assertEqual(old['xyz'],t['xyz'])
        for r in ('extended','connected'):
            for state in ('center','metal_minus','metal_plus','peptide_minus','peptide_plus'):
                a,b=[xyz(verify(next(t for t in new['tasks'] if t['task_id']==f'{r}_{state}_{metal}')['xyz'])) for metal in ('La','Ca')]
                self.assertEqual(a[0][1:],b[0][1:]);self.assertEqual(a[1:],b[1:])

    def test_actual_even_odd_units_and_negative_curvature_preserved(self):
        c,m=dft_source(DFT)
        for row in c['comparisons']:
            rep,coord,h=row['representation'],row['coordinate'],row['amplitude']
            for metal in ('La','Ca'):
                ids=[f'{rep}_{coord}_minus_{metal}',f'{rep}_center_{metal}',f'{rep}_{coord}_plus_{metal}']
                vals=np.array([c['energies'][i]['energy_hartree'] for i in ids]);got=differences((vals-vals[1])*HA_TO_KCAL,h)
                expected=row['checks'][metal]
                json.dumps(got,allow_nan=False)
                self.assertAlmostEqual(got['even_kcal_mol'],expected['even_curvature_sensitive_energy_kcal_mol'],places=10)
                self.assertAlmostEqual(got['odd_kcal_mol'],expected['odd_energy_kcal_mol'],places=10)
                if coord=='peptide' and metal=='La':self.assertLess(got['secant_curvature'],0)

    def test_partial_collection_keeps_corrections_unavailable(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'manifest.json';p.write_bytes(MANIFEST.read_bytes())
            c=collect_curvature(p)
            self.assertEqual(c['status'],'incomplete');self.assertIsNone(c['relaxation_correction_kcal_mol'])
            self.assertTrue(all(v['energy_eV'] is None for v in c['rows'].values()))

    def test_changed_real_charge_or_model_rejected(self):
        for modification in ('charge','model'):
            m=copy.deepcopy(read_json(MANIFEST))
            if modification=='charge':m['tasks'][0]['charge']+=1
            else:m['model']['solvent']='water'
            with tempfile.TemporaryDirectory() as directory:
                p=Path(directory)/'manifest.json';p.write_text(json.dumps(m))
                with self.assertRaises(InvalidArtifact):validate(p)

    def test_completed_real_screen_does_not_enable_relaxation(self):
        path=ROOT/'workspaces/mace_curvature_20260916/report_v1/result.json'
        if not path.exists():self.skipTest('actual MACE/GB curvature calls have not completed')
        r=read_json(path)
        for label in ('medium','large'):
            self.assertTrue(all(r['gates'][label].values()))
            self.assertEqual(len(r['rows'][label]),4)
            for row in r['rows'][label]:
                if row['coordinate']=='peptide':
                    self.assertLess(row['terms']['La']['MACE_plus_GB']['secant_curvature'],0)
        row=next(x for x in r['rows']['medium'] if x['representation']=='connected' and x['coordinate']=='metal')
        self.assertAlmostEqual(row['checks']['R']['curvature_even_error_kcal_mol'],.00320889408045,places=10)
        self.assertEqual(r['response_status'],'response_model_not_validated')
        self.assertIsNone(r['relaxation_correction_kcal_mol']);self.assertIsNone(r['entropy_correction_kcal_mol'])


if __name__=='__main__':unittest.main()
