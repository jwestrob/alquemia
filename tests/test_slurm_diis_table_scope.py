"""Actual 1220300 logs: paired orbital rows must not become SCF iterations."""
import json
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from slurm_native_health_watch import diis_table,inspect_text
W=ROOT/'workspaces/metal_environment_response_20260926/dy_frozen_f_direction_v1'

class RealDIISSectionScope(unittest.TestCase):
    def test_actual_orbital_output_does_not_extend_diis_table(self):
        for name in ('plus','minus','half_plus','half_minus'):
            with self.subTest(endpoint=name):
                text=(W/name/'endpoint.out').read_text()
                prefix=text[:text.index('Initializing SOSCF')]
                early,active=diis_table(prefix)
                complete,active_complete=diis_table(text)
                self.assertTrue(active)
                self.assertFalse(active_complete)
                self.assertEqual(complete,early)
                self.assertEqual([r['iteration'] for r in complete],list(range(1,8)))
                self.assertTrue(all(r['energy_hartree']<-1000 for r in complete))
                self.assertEqual(inspect_text(text)['alerts'],[])

    def test_actual_postscf_orbital_section_without_scf_header_is_not_diis(self):
        text=(W/'plus/endpoint.out').read_text()
        rows,active=diis_table(text[text.index('ORBITAL ENERGIES'):])
        self.assertEqual(rows,[])
        self.assertFalse(active)

    def test_saved_bad_receipt_contrasts_with_actual_scoped_rows(self):
        old=json.loads((W/'health_receipt.json').read_text())
        row=old['snapshot']['endpoints']['plus']
        self.assertGreater(row['last_diis_iteration'],7)
        corrected=inspect_text((W/'plus/endpoint.out').read_text())
        self.assertEqual(corrected['last_diis_iteration'],7)
        self.assertNotEqual(row['last_diis_rows'],corrected['last_diis_rows'])

if __name__=='__main__':unittest.main()
