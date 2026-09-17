"""Real archived inputs and malformed copies; scientific integrations run in Slurm."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, xyz
from mace_gb import validate, collect_gb, core_gate, build_system, EV_TO_KJ
from mace_hybrid import EV_TO_KCAL, dry_run
MANIFEST=ROOT/'workspaces/mace_gb_20260916/pilot_v1/manifest.json'
COMPLETED=ROOT/'workspaces/mace_gb_20260916/pilot_v2/manifest.json'


@unittest.skipUnless(MANIFEST.exists(),'requires real prepared GB pilot')
class FrozenGBTests(unittest.TestCase):
    @unittest.skipUnless((COMPLETED.parent/'collection_job_1200700.json').exists(),'requires actually executed GB pilot')
    def test_archived_solver_sign_units_and_missing_calibration(self):
        c=collect_gb(COMPLETED)
        self.assertTrue(c['numerical_checks_pass']);self.assertEqual(len(c['rows']),19)
        self.assertAlmostEqual(c['direct']['medium']['GB_Ca_minus_La_kcal_mol'],-84.0550594094293,places=8)
        self.assertAlmostEqual(c['checkpoint_disagreement_kcal_mol']['descriptor_R_kcal_mol'],106.38302595552523,places=8)
        self.assertIsNone(c['S_kcal_mol']);self.assertIsNone(c['calibrated_class']);self.assertIsNone(c['corrected_gradient'])
        for r in c['rows'].values():
            self.assertEqual(r['extracted_force_groups'],[0]);self.assertFalse(r['direct_Coulomb_included'])
            self.assertAlmostEqual(r['energy_eV']*EV_TO_KCAL,r['GB_reaction_kcal_mol'],places=9)

    def test_real_manifest_and_units(self):
        result=dry_run(MANIFEST)
        self.assertEqual(result['tasks'],19)
        self.assertAlmostEqual(EV_TO_KCAL,EV_TO_KJ/4.184,places=13)

    def test_scientific_changes_and_corrupt_density_are_rejected(self):
        m=read_json(MANIFEST)
        mutations=[lambda d:d['model'].update(solvent_dielectric=40.),
                   lambda d:d['tolerances'].update(energy_kcal_mol=1.),
                   lambda d:d['tasks'][0].update(source_vacuum_energy_eV=0.),
                   lambda d:d['tasks'][0]['source_density'].update(sha256='0'*64),
                   lambda d:d['tasks'][0].update(cache_key='0'*64)]
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json'
            for mutate in mutations:
                altered=copy.deepcopy(m);mutate(altered);p.write_text(json.dumps(altered))
                with self.assertRaises(InvalidArtifact):validate(p)

    def test_missing_partial_outputs_remain_unavailable(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json';p.write_bytes(MANIFEST.read_bytes())
            c=collect_gb(p)
            self.assertEqual(c['status'],'incomplete');self.assertFalse(c['numerical_checks_pass'])
            self.assertIsNone(c['checkpoint_disagreement_kcal_mol'])
            for d in c['direct'].values():self.assertIsNone(d['descriptor_R_kcal_mol'])
            self.assertEqual(core_gate(p)['status'],'fail')

    def test_only_reaction_energy_is_in_selected_force_group(self):
        try:
            import openmm as mm
            from openmm import unit
        except ImportError:
            self.skipTest('OpenMM unavailable; no fabricated solver output')
        task=read_json(MANIFEST)['tasks'][0]
        atoms=xyz(verify(task['xyz']));q=np.load(verify(task['source_density']))[:,0]
        for solver in ('native','custom'):
            system=build_system(atoms,q,solver,78.5)
            forces=list(system.getForces());self.assertEqual(len(forces),2)
            gb,nb=forces
            self.assertEqual(gb.getForceGroup(),0);self.assertEqual(nb.getForceGroup(),1)
            self.assertEqual(nb.getReactionFieldDielectric(),1.)
            self.assertEqual(nb.getNonbondedMethod(),mm.NonbondedForce.NoCutoff)
            for i,charge in enumerate(q):
                self.assertEqual(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge),charge)
                params=gb.getParticleParameters(i)
                value=params[0].value_in_unit(unit.elementary_charge) if solver=='native' else params[0]
                self.assertEqual(value,charge)
            if solver=='native':self.assertEqual(gb.getSurfaceAreaEnergy().value_in_unit(unit.kilojoule_per_mole/unit.nanometer**2),0.)


if __name__=='__main__':unittest.main()
