#!/usr/bin/env python3
"""Finite source/fixture controls. No engine imports or native calculations."""
import importlib.util
from pathlib import Path
import re
import unittest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'


def get(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


p=get(HERE/'PATH.py','render');r=get(HERE/'RUN.py','runner')
h=get(HERE.parent/'gpu_engine_qualification_v1/RUN.py','qualified')
fixture=W/'explicit_model_v1/prepared/water_central_La'
PRM=(fixture/'amoeba09.prm').read_text();KEY=(fixture/'system.key').read_text();XYZ=(fixture/'system.xyz').read_text()


class Controls(unittest.TestCase):
    def test_actual_water_type_class_and_native_values(self):
        resolved=p.endpoint_resolution(p.parse_parameters(PRM,KEY),[36,37])
        self.assertEqual(resolved['36']['La']['nonmetal_class'],34)
        self.assertEqual(resolved['37']['La']['nonmetal_class'],35)
        self.assertAlmostEqual(resolved['36']['La']['Rmin_pair_A'],3.6985300908646654,places=13)
        self.assertAlmostEqual(resolved['36']['Dy']['epsilon_kcal_mol'],.37332330078584813,places=13)

    def test_explicit_override_is_full_pair_Rmin(self):
        resolved=p.resolve_pair(p.parse_parameters(PRM,KEY),500,403)
        self.assertEqual(resolved['Rmin_pair_A'],3.6318)
        self.assertEqual(resolved['epsilon_kcal_mol'],.228)
        self.assertEqual(resolved['thole'],.1518)

    def test_endpoint_and_midpoint_pair_parameters(self):
        _,a=p.render_key(PRM,KEY,[403],0);_,b=p.render_key(PRM,KEY,[403],1);text,m=p.render_key(PRM,KEY,[403],.5)
        self.assertAlmostEqual(m['actual_pair_parameters']['403']['Rmin_pair_A'],3.4895188293023635,places=13)
        self.assertAlmostEqual(m['actual_pair_parameters']['403']['thole'],.1454,places=13)
        self.assertEqual(a['alpha_A3'],1.6597);self.assertEqual(b['alpha_A3'],.728)
        self.assertIn('polpair 403 502',text)

    def test_rendered_water_distinct_pair_indexes(self):
        key,_=p.render_key(PRM,KEY,[36,37],.5)
        self.assertIn('vdwpr 34 502',key);self.assertIn('polpair 36 502',key)
        self.assertIn('vdwpr 35 502',key);self.assertIn('polpair 37 502',key)
        self.assertTrue(key.startswith(KEY.rstrip()+'\n\n'))

    def test_coordinates_nonmetal_lines_topology_preserved(self):
        new,old=p.render_xyz(XYZ)
        self.assertEqual(new.splitlines()[2:],XYZ.splitlines()[2:])
        self.assertEqual(new.splitlines()[1].split()[2:5],old[0][2:5])
        self.assertEqual(new.splitlines()[1].split()[5:],[str(p.TYPE)])

    def test_asymmetric_metal_bond_rejected(self):
        lines=XYZ.splitlines();lines[2]+=' 1'
        with self.assertRaisesRegex(ValueError,'bond'):p.render_xyz('\n'.join(lines)+'\n')

    def test_existing_class_alias_rejected(self):
        with self.assertRaisesRegex(ValueError,'class already'):p.parse_parameters(PRM,KEY+'\natom 900 502 C "alias" 6 12.0 0\n')

    def test_unsupported_convention_rejected(self):
        with self.assertRaisesRegex(ValueError,'conventions'):p.parse_parameters(PRM,KEY+'\nradiustype SIGMA\n')

    def test_direct_damping_rejected(self):
        with self.assertRaisesRegex(ValueError,'direct Thole'):p.parse_parameters(PRM,KEY+'\npolarize 500 1.6597 0.25 0.6\n')

    def test_no_fractional_charge_or_mass(self):
        for lam in p.LAMBDAS:
            text,meta=p.render_key(PRM,KEY,[36,37],lam)
            self.assertEqual(meta['charge'],3);self.assertEqual(meta['fixed_mass_u'],138.91)
            self.assertIn('multipole 502 0 0 3.000000000000000',text)
        self.assertEqual(2*len(p.LAMBDAS)+2,20)

    def test_analytic_saved_native_output(self):
        # Existing CPU Y/N output provides actual positive completion evidence.
        import json
        result=json.loads((W/'explicit_model_v1/matched_dry_v1/run/RESULT.json').read_text())
        row=result['rows'][0];native=row['native'];text=Path(native['gradient_output']['path']).read_text()
        parsed=r.parse_analytic(text,row['input']['atoms'],h)
        self.assertEqual(parsed['total_energy_kcal_mol'],native['total_energy_kcal_mol'])
        self.assertEqual(parsed['force_kcal_mol_A'],native['force_kcal_mol_A'])

    def test_native_zero_exit_failure_text_rejected(self):
        with self.assertRaisesRegex(ValueError,'failure'):r.parse_analytic('Uncaught exception: device error\n',4,h)

    def test_truncated_gradient_rejected(self):
        text='Total Potential Energy : 1.0 Kcal/mole\n Anlyt 1 1.0 2.0 3.0 3.74\n'
        with self.assertRaisesRegex(ValueError,'Incomplete'):r.parse_analytic(text,4,h)

    def test_parity_and_force_denominator(self):
        x={'total_energy_kcal_mol':1.,'force_kcal_mol_A':[[1.,2.,3.]]}
        y={'total_energy_kcal_mol':1.,'force_kcal_mol_A':[[1.01,2.,3.]]}
        result=r.compare(x,y);self.assertFalse(result['passed']);self.assertEqual(result['force_RMS_denominator'],3)


if __name__=='__main__':unittest.main()
