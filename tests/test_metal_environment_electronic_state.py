"""Actual archived outputs and states; no fabricated quantum result."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json
from metal_environment_electronic_state import describe,parse,electronic_evidence
from mace_omol_vacuum import parse_endpoint

class ExplicitStateTests(unittest.TestCase):
    def test_actual_ca_la_outputs_unchanged(self):
        w=ROOT/'workspaces/metal_environment_response_20260926/reference_scout_v1'
        m=read_json(w/'manifest.json')
        for t in m['tasks']:
            with self.subTest(task=t['task_id']):
                old=parse_endpoint(t,t['output_path'],t['engrad_path'],permanent_field=t['environment']!='isolated')
                new=parse(t,t['output_path'],t['engrad_path'],embedded=t['environment']!='isolated')
                self.assertEqual(old['energy_hartree'],new['energy_hartree'])
                self.assertEqual(old['gradient_kcal_mol_per_A'],new['gradient_kcal_mol_per_A'])
                self.assertFalse(new['electronic_state_qualified'])
                self.assertEqual(new['executed_hftype'], 'RHF')
                self.assertEqual(new['spin_evidence']['stability_status'], 'not_run_or_not_printed')
                self.assertIsNone(new['spin_evidence']['s2_last'])
                self.assertIsNone(new['spin_evidence']['population_analyses']['mulliken'])

    def test_actual_lanm_state_counts(self):
        a=read_json(ROOT/'diagnostics/metal_environment_response_20260926/lanm_state_feasibility/STATE_AUDIT.json')
        for row in a['states']:
            with self.subTest(source=row['source'],metal=row['metal']):
                d=describe(row['xyz'],row['metal'],row['charge'],row['multiplicity'])
                for k in ('all_electron_count','explicit_electrons','ecp_core_electrons','alpha_electrons','beta_electrons'):
                    self.assertEqual(d[k],row[k])
                if row['metal']=='Dy':
                    self.assertFalse(d['f_in_core'])
                    with self.assertRaises(InvalidArtifact):describe(row['xyz'],'Dy',row['charge'],1)
                    with self.assertRaises(InvalidArtifact):describe(row['xyz'],'Dy',row['charge']+1,6)

    def test_real_archived_open_shell_diagnostics_format(self):
        # Real ORCA6.1.1 Er output tests syntax only, not Er or Dy capability,
        # convergence, chemistry or scientific validation of the new protocol.
        fixture=ROOT/'legacy/qmmm/aquo_panel/Er'
        text=(fixture/'opt_aquo_Er_cn8.out').read_text()
        lines=(fixture/'aquo_Er_cn8.xyz').read_text().splitlines()
        symbols=[line.split()[0] for line in lines[2:2+int(lines[0])]]
        e=electronic_evidence(text,symbols,4)
        self.assertEqual(e['s2_last'],3.752366)
        self.assertEqual(e['spin_pure_s2'],3.75)
        self.assertEqual(e['population_analyses']['mulliken']['atoms'][0]['spin_population'],3.017912)
        self.assertEqual(e['population_analyses']['loewdin']['atoms'][0]['spin_population'],3.016722)
        self.assertAlmostEqual(e['population_analyses']['mulliken']['spin_sum'],3.0,places=5)
        self.assertEqual(e['stability_status'],'not_run_or_not_printed')
        self.assertFalse(e['spin_localization_qualified'])
        # Explicitly corrupted copy of this real output must not silently map
        # a partial population list to atom coordinates.
        corrupt=text.replace('   0 Er:    1.909484    3.017912','   1 Er:    1.909484    3.017912',1)
        self.assertNotEqual(corrupt,text)
        with self.assertRaises(InvalidArtifact):electronic_evidence(corrupt,symbols,4)

    @unittest.skip('No executed native r2SCAN-3c Dy EnGrad in this campaign yet; scout must establish output support')
    def test_native_dy_execution(self):
        pass

if __name__=='__main__':unittest.main()
