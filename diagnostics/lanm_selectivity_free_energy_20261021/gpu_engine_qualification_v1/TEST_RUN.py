#!/usr/bin/env python3
"""Saved-output qualification of the parser and its refusal paths; no GPU."""
import copy
import importlib.util
import io
import json
from pathlib import Path
import re
import time
import unittest
import xml.etree.ElementTree as ET

DIRECTORY = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('gpu_qualification', DIRECTORY/'RUN.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
config = json.loads((DIRECTORY/'CONFIG_TEMPLATE.json').read_text())
controls = runner.load(config['cpu_controls_result'])
ROWS = {r['id']: r for r in controls['rows'] if r['status'] == 'complete'}
TEXTS = {key: runner.verify(row['native']['gradient_output']).read_text()
         for key, row in ROWS.items()}


def suite_xml():
    """Synthetic reporter fixture follows pinned Catch2 XML writer, not a native result."""
    root = ET.Element('Catch2TestRun', {'filters': '[polpair]', 'xml-format-version': '2'})
    for name in sorted(runner.SUITE_NAMES):
        case = ET.SubElement(root, 'TestCase', {'name': name, 'tags': '[ff][amoeba][polpair]'})
        ET.SubElement(case, 'OverallResult', {'success': 'true'})
    ET.SubElement(root, 'OverallResults', {'successes': '100', 'failures': '0', 'expectedFailures': '0'})
    ET.SubElement(root, 'OverallResultsCases', {'successes': '4', 'failures': '0', 'expectedFailures': '0'})
    return root


def encoded(root):
    return ET.tostring(root, encoding='unicode')


class GradientTests(unittest.TestCase):
    def test_all_four_actual_native_outputs_exact_replay(self):
        self.assertEqual(set(ROWS), set(runner.CONTROL_IDS))
        for key in runner.CONTROL_IDS:
            row = ROWS[key]
            with self.subTest(control=key):
                parsed = runner.parse_gradient(TEXTS[key], row['input']['atoms'])
                self.assertEqual(parsed['total_energy_kcal_mol'], row['native']['total_energy_kcal_mol'])
                self.assertEqual(parsed['force_kcal_mol_A'], row['native']['force_kcal_mol_A'])
                comparison = runner.compare(parsed, row['native'], row['input']['atoms'])
                self.assertTrue(comparison['passed'])
                self.assertEqual(comparison['force_RMS_denominator'], 3*row['input']['atoms'])
                self.assertEqual(comparison['force_component_RMS_error_kcal_mol_A'], 0)

    def setUp(self):
        self.row = ROWS['acetate_bi_central_La']
        self.text = TEXTS['acetate_bi_central_La']
        self.atoms = self.row['input']['atoms']
        self.parsed = runner.parse_gradient(self.text, self.atoms)

    def test_tinker9_does_not_need_CPU_specific_breakdown_header(self):
        text = self.text.replace('Cartesian Gradient Breakdown over Individual Atoms :', '')
        self.assertEqual(runner.parse_gradient(text, self.atoms), self.parsed)

    def test_zero_exit_native_exception_with_otherwise_complete_output_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Native failure'):
            runner.parse_gradient(self.text+'\n Terminating with uncaught exception : test\n', self.atoms)

    def test_native_nonconvergence_rejected(self):
        with self.assertRaises(ValueError):
            runner.parse_gradient(self.text+'\n Induced Dipoles not converged\n', self.atoms)

    def test_nonfinite_and_overflow_rejected(self):
        for bad in ('NaN', 'Infinity', '1e309'):
            with self.subTest(value=bad), self.assertRaises(ValueError):
                text = re.sub(r'(Total Potential Energy\s*:\s*)'+runner.NUMBER,
                              lambda m: m[1]+bad, self.text)
                runner.parse_gradient(text, self.atoms)

    def test_missing_index_duplicate_index_and_missing_summary_rejected(self):
        bad_texts = [re.sub(r'^\s*Anlyt\s+8\s+.*$', '', self.text, flags=re.M),
                     re.sub(r'^(\s*Anlyt\s+)8(\s+)', r'\g<1>7\2', self.text, flags=re.M),
                     re.sub(r'^\s*Numer\s+Total Gradient Norm Value.*$', '', self.text, flags=re.M)]
        for text in bad_texts:
            with self.subTest(kind=text[-100:]), self.assertRaises(ValueError):
                runner.parse_gradient(text, self.atoms)

    def test_duplicate_energy_wrong_units_and_nonnumeric_suffix_rejected(self):
        for text in (self.text+'\n Total Potential Energy : 0.0 Kcal/mole\n',
                     self.text.replace('Kcal/mole', 'Hartree'),
                     self.text.replace('-493.63525889', '-493.63525889xyz')):
            with self.subTest(tail=text[-50:]), self.assertRaises(ValueError):
                runner.parse_gradient(text, self.atoms)

    def test_finite_difference_gate_from_actual_output_mutation(self):
        text = re.sub(r'^(\s*Numer\s+1\s+)('+runner.NUMBER+')',
                      lambda m: m[1]+str(float(m[2])+0.002), self.text, flags=re.M)
        parsed = runner.parse_gradient(text, self.atoms)
        self.assertFalse(parsed['finite_difference']['passed'])
        self.assertFalse(runner.compare(parsed, self.row['native'], self.atoms)['passed'])

    def test_energy_gate(self):
        parsed = copy.deepcopy(self.parsed)
        parsed['total_energy_kcal_mol'] += 0.0002
        comparison = runner.compare(parsed, self.row['native'], self.atoms)
        self.assertFalse(comparison['energy_passed'])
        self.assertFalse(comparison['passed'])

    def test_force_gate_and_denominator(self):
        parsed = copy.deepcopy(self.parsed)
        parsed['force_kcal_mol_A'][0][0] += 0.002
        comparison = runner.compare(parsed, self.row['native'], self.atoms)
        self.assertFalse(comparison['force_passed'])
        self.assertFalse(comparison['passed'])
        self.assertAlmostEqual(comparison['force_component_RMS_error_kcal_mol_A'],
                               0.002/(3*self.atoms)**0.5, places=12)

    def test_missing_native_atom_and_nonfinite_reference_rejected(self):
        parsed = copy.deepcopy(self.parsed)
        parsed['force_kcal_mol_A'].pop()
        with self.assertRaises(ValueError):
            runner.compare(parsed, self.row['native'], self.atoms)
        reference = copy.deepcopy(self.row['native'])
        reference['force_kcal_mol_A'][0][0] = float('nan')
        with self.assertRaises(ValueError):
            runner.compare(self.parsed, reference, self.atoms)


class SuiteTests(unittest.TestCase):
    def test_exact_four_case_report_success(self):
        self.assertEqual(runner.parse_suite(encoded(suite_xml()), '')['test_case_count'], 4)

    def test_native_failure_in_console_and_captured_output(self):
        with self.assertRaises(ValueError):
            runner.parse_suite(encoded(suite_xml()), 'Terminating with uncaught exception : bad')
        doc = suite_xml()
        ET.SubElement(doc.find('TestCase'), 'StdOut').text = 'Polarization not converged'
        with self.assertRaises(ValueError):
            runner.parse_suite(encoded(doc), '')

    def test_missing_duplicate_or_additional_case(self):
        for mode in ('missing', 'duplicate', 'extra'):
            doc = suite_xml()
            if mode == 'missing':
                doc.remove(doc.find('TestCase'))
            elif mode == 'duplicate':
                doc.findall('TestCase')[0].set('name', doc.findall('TestCase')[1].get('name'))
            else:
                ET.SubElement(doc, 'TestCase', {'name': 'Other'})
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                runner.parse_suite(encoded(doc), '')

    def test_failed_case_assertion_expected_failure_or_empty_assertions(self):
        for field, value in (('failures', '1'), ('expectedFailures', '1'), ('successes', '0')):
            doc = suite_xml()
            doc.find('OverallResults').set(field, value)
            with self.subTest(field=field), self.assertRaises(ValueError):
                runner.parse_suite(encoded(doc), '')
        doc = suite_xml()
        doc.find('TestCase/OverallResult').set('success', 'false')
        with self.assertRaises(ValueError):
            runner.parse_suite(encoded(doc), '')

    def test_wrong_filter_missing_totals_and_malformed_XML(self):
        doc = suite_xml()
        doc.set('filters', '[ff]')
        with self.assertRaises(ValueError):
            runner.parse_suite(encoded(doc), '')
        doc = suite_xml()
        doc.remove(doc.find('OverallResultsCases'))
        with self.assertRaises(ValueError):
            runner.parse_suite(encoded(doc), '')
        with self.assertRaises(ET.ParseError):
            runner.parse_suite('<Catch2TestRun>', '')


if __name__ == '__main__':
    started = time.monotonic()
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.loadTestsFromModule(__import__('__main__'))
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    print(stream.getvalue())
    print(json.dumps({'status': 'passed' if result.wasSuccessful() else 'failed',
                      'tests_run': result.testsRun, 'failures': len(result.failures),
                      'errors': len(result.errors), 'elapsed_seconds': time.monotonic()-started,
                      'adapter': runner.pin(DIRECTORY/'RUN.py'),
                      'actual_CPU_result': config['cpu_controls_result'],
                      'native_program_invocations': 0, 'GPU_calls': 0,
                      'suite_XML_fixture': 'synthetic format fixture; not native execution'}))
    raise SystemExit(0 if result.wasSuccessful() else 1)
