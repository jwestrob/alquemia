"""Geometry/algebra tests use the archived real 1H4I protein parent; no FF calls."""
import copy
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import scaffold_accommodation as scaffold
from affordable_common import read_json,verify,InvalidArtifact


class ActualParentGeometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        result=read_json(ROOT/'diagnostics/scaffold_environment_20260922/RESULT.json')
        inv=read_json(verify(result['inventory_artifact']))
        row=next(r for r in inv['rows'] if r['case_id']=='1H4I')
        cls.row=row
        cls.atoms=read_json(verify(row['artifacts']['parent_atoms.json']))
        for a in cls.atoms:
            a['name']=a['source_id'].rsplit('/',1)[1]
            a['residue_id']=a['source_id'].rsplit('/',1)[0]
        cls.positions=np.asarray([a['xyz_A'] for a in cls.atoms])
        terms=read_json(verify(row['artifacts']['term_supports.json']))
        cls.bonds=next(f['supports'] for f in terms if f['type']=='HarmonicBondForce')
        cls.mapping=read_json(verify(row['artifacts']['context_mapping.json']))
        core=read_json(verify(row['endpoints']['Ca']['mapping']))['core']
        cls.core_heavy=[p for p,a in zip(core['positions_A'],core['source_atom_metadata'])
                        if a and a['element']!='H']
        cls.mobile,cls.residues=scaffold.mobile_atoms(cls.atoms,cls.bonds,cls.positions,
                                     cls.core_heavy,cls.mapping['local_parent_indices'])
        # Freeze actual first-shell anchor glutamate O atoms, as declared role targets.
        cls.donors=[i for i,a in enumerate(cls.atoms) if a['source_id'] in ('A/177//OE1','A/177//OE2')]
        if len(cls.donors)!=2:raise AssertionError('actual donor fixture mapping changed')
        cls.free=sorted(set(cls.mobile)-set(cls.donors))
        cls.constraints=scaffold.Constraints(cls.positions,cls.positions,cls.bonds,cls.free)

    def test_complete_residues_and_real_peptide_neighbours(self):
        selected=set(self.mobile)
        for residue in self.residues:
            self.assertTrue(all(i in selected for i,a in enumerate(self.atoms) if a['residue_id']==residue))
        self.assertTrue(set(self.mapping['local_parent_indices'])<=selected)
        self.assertLess(len(selected),len(self.atoms))

    def test_actual_bond_jacobian(self):
        residual,jac=self.constraints.residual_jacobian(self.positions)
        self.assertLess(np.max(abs(residual)),1e-12)
        column=3*len(self.free)//2; index,axis=self.free[column//3],column%3
        plus=self.positions.copy();minus=self.positions.copy();plus[index,axis]+=1e-6;minus[index,axis]-=1e-6
        fd=(self.constraints.residual_jacobian(plus)[0]-self.constraints.residual_jacobian(minus)[0])/2e-6
        self.assertLess(np.max(abs(fd-jac[:,column].toarray().ravel())),1e-7)

    def test_projection_and_retraction_real_coordinate_displacement(self):
        # Algebraic vector from actual source coordinates, not a fabricated scientific force.
        displacement=self.positions[self.free]-self.positions[self.free].mean(axis=0)
        vector,info=self.constraints.project(self.positions,displacement.ravel())
        _,jac=self.constraints.residual_jacobian(self.positions)
        self.assertLess(np.max(abs(jac@vector)),1e-5)
        vector=vector.reshape((-1,3));vector*=.005/np.max(np.linalg.norm(vector,axis=1))
        trial=self.positions.copy();trial[self.free]+=vector
        repaired,info=self.constraints.retract(trial)
        self.assertTrue(np.array_equal(repaired[self.constraints.fixed],self.positions[self.constraints.fixed]))
        self.assertLessEqual(info['max_bond_error_A'],1e-8)

    def test_real_source_stereochemistry_and_no_hidden_forces(self):
        geometry=scaffold.Geometry(self.atoms,self.bonds,self.positions,self.positions,self.free,[])
        check=geometry.check(self.positions)
        self.assertTrue(check['pass'])
        self.assertGreater(len(geometry.guards['chirality']),100)
        self.assertGreater(len(geometry.guards['peptides']),100)
        self.assertFalse(scaffold.SETTINGS['forcefield_energy_added_to_score'])

    def test_corrupted_real_fixed_target_rejected(self):
        initial=self.positions.copy()
        fixed=self.constraints.fixed
        pair=next((a,b) for a,b in self.bonds if a in set(fixed) and b in set(fixed))
        initial[pair[0],0]+=.1
        with self.assertRaises(InvalidArtifact):scaffold.Constraints(self.positions,initial,self.bonds,self.free)

    def test_corrupted_real_bond_and_domain_are_not_admitted(self):
        geometry=scaffold.Geometry(self.atoms,self.bonds,self.positions,self.positions,self.free,[])
        invalid=self.positions.copy();i=next(i for i in self.free if self.atoms[i]['element']!='H')
        invalid[i,0]+=1.
        self.assertFalse(geometry.check(invalid)['pass'])


class DeclaredEightParents(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=read_json(ROOT/'workspaces/collective_scaffold_20260922/parents_v1/manifest.json')

    def test_all24_actual_target_constraints_and_contexts(self):
        self.assertEqual(len(self.manifest['cases']),8)
        count=0
        for case in self.manifest['cases']:
            data=scaffold.load_parent(case)
            self.assertEqual(case['status'],'prepared')
            self.assertEqual(set(case['targets']),set(scaffold.TARGETS))
            for target in scaffold.TARGETS:
                _,initial,constraints,geometry=scaffold.target_problem(case,target,data)
                self.assertTrue(geometry.check(initial)['pass'])
                self.assertTrue(set(data['donors'])<=set(constraints.fixed))
                self.assertEqual(len(data['positions']),case['parent_atom_count'])
                count+=1
        self.assertEqual(count,24)

    def test_parent_proton_and_energy_model_are_unchanged(self):
        for case in self.manifest['cases']:
            self.assertEqual(case['energy_evaluations'],0)
            self.assertEqual(case['force_evaluations'],0)
            nonprotein=read_json(verify(case['nonprotein_inventory']))
            self.assertFalse(nonprotein['FF_parameters_assigned'])
            self.assertEqual(case['origins']['La']['charge']-case['origins']['Ca']['charge'],1)
        self.assertEqual(scaffold.SETTINGS['precision'],'double')
        self.assertEqual(scaffold.SETTINGS['maximum_accepted_iterations'],200)

    def test_corrupted_real_target_coordinate_rejected(self):
        case=self.manifest['cases'][0];data=scaffold.load_parent(case)
        wrong=copy.deepcopy(case);wrong['targets']['origin']['donor_positions_A'][0][0]+=.1
        with self.assertRaises(InvalidArtifact):scaffold.target_problem(wrong,'origin',data)

    def test_achiral_CH2_sites_do_not_invent_stereocenters(self):
        for cid,aid in [('4MAE','A/387//CE'),('q9z4j7-pqq-la_model','A/440//CG')]:
            case=next(c for c in self.manifest['cases'] if c['case_id']==cid)
            data=scaffold.load_parent(case)
            guards=scaffold.covalent_guards(data['atoms'],data['bonds'],data['positions'])
            centers={data['atoms'][r['center']]['source_id'] for r in guards['chirality']}
            self.assertNotIn(aid,centers)
            for r in guards['chirality']:
                a=data['atoms'][r['center']]
                self.assertTrue((a['name']=='CA' and a['resname']!='GLY') or
                                (a['name']=='CB' and a['resname'] in ('ILE','THR')))

    def test_real_stereocenter_handedness_and_substituents_protected(self):
        case=self.manifest['cases'][0];data,_,_,geometry=scaffold.target_problem(case,'origin')
        center=next(r for r in geometry.guards['chirality'] if data['atoms'][r['center']]['resname']=='THR' and
                    data['atoms'][r['center']]['name']=='CB')
        a,b=center['neighbours'][:2];corrupt=data['positions'].copy();corrupt[[a,b]]=corrupt[[b,a]]
        self.assertIn(center['center'],geometry.check(corrupt)['chirality_failures'])
        wrong=copy.deepcopy(data['atoms']);wrong[a]['name']='CORRUPTED_REAL_NEIGHBOR'
        with self.assertRaises(InvalidArtifact):scaffold.covalent_guards(wrong,data['bonds'],data['positions'])


if __name__=='__main__':unittest.main()
