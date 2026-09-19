"""Real archived forces: mapping/parser/algebra tests, no molecular evaluations."""
import copy
import json
from pathlib import Path
import re
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import response_probe_archive as probe


class ArchivedResponseProbeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=ROOT/'workspaces/response_probe_20260919/features_v2.json'
        if not path.exists():raise unittest.SkipTest('real prepared force-probe archive unavailable')
        cls.data=json.loads(path.read_text()); cls.rows={r['case_id']:r for r in cls.data['rows']}
        cls.pqq=cls.data['rows'][0]; cls.direct=cls.rows['GGR_extended']

    def test_paired_actual_coordinates_and_corrupted_charge(self):
        for row in self.rows.values():
            _,coordinates=probe.paired_coordinates(row['states'])
            self.assertTrue(np.isfinite(coordinates).all())
        corrupted=copy.deepcopy(self.pqq['states']);corrupted['La']['charge']+=1
        with self.assertRaisesRegex(probe.InvalidArtifact,'charge'):
            probe.paired_coordinates(corrupted)

    def test_force_sign_units_and_physical_donor_extraction(self):
        for row in (self.pqq,self.direct):
            _,coordinates=probe.paired_coordinates(row['states'])
            actual={m:-np.load(probe.verify(row['force_artifacts'][m]))*probe.EV_TO_KCAL for m in ('Ca','La')}
            mapped={m:probe.map_gradient(g,row['physical_mapping']) for m,g in actual.items()}
            difference=[mapped['Ca'][d['key']]-mapped['La'][d['key']] for d in row['donors']]
            np.testing.assert_allclose(difference,row['donor_delta_gradient_kcal_mol_A'],atol=1e-12,rtol=0)
            # Endpoint reversal negates signed radial load while preserving tangential fraction.
            result=probe.observables(-np.array(difference),[d['xyz_A'] for d in row['donors']],coordinates[0])
            self.assertAlmostEqual(result[probe.PRIMARY[0]],-row['features'][probe.PRIMARY[0]],places=12)
            self.assertAlmostEqual(result[probe.PRIMARY[1]],row['features'][probe.PRIMARY[1]],places=12)

    def test_arbitrary_rigid_transform_of_real_donor_vectors(self):
        rotation=np.array([[.36,-.48,.8],[.8,.6,0],[-.48,.64,.6]])
        np.testing.assert_allclose(rotation@rotation.T,np.eye(3),atol=1e-15)
        for row in self.rows.values():
            _,coordinates=probe.paired_coordinates(row['states'])
            positions=np.array([d['xyz_A'] for d in row['donors']]);g=np.array(row['donor_delta_gradient_kcal_mol_A'])
            result=probe.observables(g@rotation,positions@rotation+[50,-31,12],coordinates[0]@rotation+[50,-31,12])
            for key in probe.PRIMARY:self.assertAlmostEqual(result[key],row['features'][key],places=10)

    def test_actual_cap_chain_rule_against_geometry_derivatives(self):
        for row in (self.pqq,self.direct):
            preparation=probe.read_json(probe.verify(row['preparation']))
            _,coordinates=probe.paired_coordinates(row['states'])
            mapper=probe.canonical_mapping if row['target']=='PQQ_functional_class' else probe.generic_mapping
            mapping,physical,_,_=mapper(preparation,coordinates)
            for cap in (m for m in mapping if m['kind']=='cap'):
                a=np.array(physical[cap['retained']]['xyz_A']);b=np.array(physical[cap['omitted']]['xyz_A'])
                # Use the pinned chemistry's exact length, not rounded cap coordinates.
                length=1.09 if row['target']=='PQQ_functional_class' else next(
                    c['length_A'] for c in preparation['atom_graph']['source_to_qm']
                    if c['qm_index']==cap['qm_index'])
                def cap_position(x,y):return x+length*(y-x)/np.linalg.norm(y-x)
                h=1e-5
                for label in ('retained','omitted'):
                    numerical=np.zeros((3,3))
                    for axis in range(3):
                        offset=np.eye(3)[axis]*h
                        numerical[:,axis]=(cap_position(a+offset,b)-cap_position(a-offset,b))/(2*h) if label=='retained' else (
                            cap_position(a,b+offset)-cap_position(a,b-offset))/(2*h)
                    np.testing.assert_allclose(numerical,cap['J_'+label],atol=2e-9,rtol=0)

    def test_missing_or_corrupted_real_mapping_fails(self):
        preparation=probe.read_json(probe.verify(self.pqq['preparation']))
        _,coordinates=probe.paired_coordinates(self.pqq['states'])
        corrupted=copy.deepcopy(preparation);corrupted['qm_fragments'][0]['atom_records'][0]['xyz_A'][0]+=.02
        with self.assertRaises(probe.InvalidArtifact):probe.canonical_mapping(corrupted,coordinates)
        corrupted=copy.deepcopy(preparation);corrupted['qm_fragments'][-1]['atom_records'][-1]['origin']='unsupported_corrupted_origin'
        with self.assertRaises(probe.InvalidArtifact):probe.canonical_mapping(corrupted,coordinates)
        invalid=np.array(self.pqq['donor_delta_gradient_kcal_mol_A']);invalid[0,0]=float('nan')
        with self.assertRaises(probe.InvalidArtifact):
            probe.observables(invalid,[d['xyz_A'] for d in self.pqq['donors']],coordinates[0])

    def test_real_DFT_pair_Hamiltonian_charge_and_analytic_method(self):
        methods=[]
        for row in self.rows.values():
            if row['target']!='direct_site_affinity_direction':continue
            for metal,artifacts in row['DFT_artifacts'].items():
                text=probe.verify(artifacts['input']).read_text()
                match=re.search(r'^\* xyzfile (-?\d+) (\d+) core.xyz$',text,re.M)
                self.assertIsNotNone(match)
                self.assertEqual(int(match[1]),row['states'][metal]['charge'])
                self.assertEqual(int(match[2]),1)
                methods.append(re.sub(r'^\* xyzfile.*$','',text,flags=re.M).strip())
        self.assertEqual(set(methods),{'! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF EnGrad'})

    def test_complete_PQQ_and_declared_donors_no_nitrogen_loss(self):
        for row in self.rows.values():
            self.assertTrue(all(d['distance_A']<=3.2 for d in row['donors']))
            self.assertEqual(len(row['donors']),len({d['key'] for d in row['donors']}))
            if row['target']=='PQQ_functional_class':
                self.assertEqual(row['cofactor']['atom_count'],27)
                self.assertTrue(any(d['donor_type']=='cofactor_N' for d in row['donors']))
        self.assertTrue(any(3.1<d['distance_A']<=3.2 for d in self.pqq['donors']))

    def test_saved_folds_and_DFT_offset_replay(self):
        path=ROOT/'workspaces/response_probe_20260919/evaluation_v2.json'
        if not path.exists():self.skipTest('real completed grouped fit unavailable')
        result=probe.read_json(path)
        self.assertLess(result['DFT_S_to_R_max_logit_replay_error'],1e-7)
        for arm in result['results'].values():
            self.assertEqual(len(arm['predictions']),25)
            for fold in arm['folds']:
                held={r['case_id'] for r in self.rows.values() if r.get('validation_group')==fold['held_group']}
                self.assertFalse(held & set(fold['model']['training_cases']))
            self.assertEqual(len(arm['consumed_transfers']),3)

    def test_direct_metadata_matches_authoritative_ledger(self):
        ledger=probe.read_json(probe.verify(self.data['sources']['old_features']))
        byid={r['case_id']:r for r in ledger['rows']}
        for name,row in self.rows.items():
            if row['target']!='direct_site_affinity_direction':continue
            source=byid['GGR_1GLG' if name.startswith('GGR') else name]
            for field in ('label','biological_group','evidence_stratum'):
                self.assertEqual(row[field],source[field])
        # The actually produced wrong v1 reporting metadata must be rejected.
        w=ROOT/'workspaces/response_probe_20260919'
        with self.assertRaisesRegex(probe.InvalidArtifact,'reporting metadata'):
            probe.transfer(w/'features_v1.json',w/'evaluation_v2.json',
                           ROOT/'diagnostics/response_probe_20260919/TRANSFER_PLAN.md',w/'must_not_exist.json')

    def test_metadata_repair_preserves_all_numerical_features_and_fits(self):
        w=ROOT/'workspaces/response_probe_20260919'
        before=probe.read_json(w/'features_v1.json')
        for a,b in zip(before['rows'],self.data['rows']):
            for key in ('features','donors','physical_mapping','donor_delta_gradient_kcal_mol_A','DFT_response'):
                self.assertEqual(a.get(key),b.get(key))
        self.assertEqual(probe.read_json(w/'evaluation_v1.json')['models'],probe.read_json(w/'evaluation_v2.json')['models'])

    def test_transfer_difference_matches_saved_model_algebra(self):
        w=ROOT/'workspaces/response_probe_20260919'; result=probe.read_json(w/'transfer_v1.json')
        fitted=probe.read_json(w/'evaluation_v2.json')
        self.assertEqual(result['new_fits'],0)
        for comparison in result['rows']:
            for arm,value in comparison['arms'].items():
                model=fitted['models'][arm]; logits=[]
                for case in (comparison['alpha'],comparison['GGR']):
                    row=self.rows[case]; values={**row['features'],'DFT_R_kcal_mol':row['DFT_R_kcal_mol']}
                    standard=(np.array([values[k] for k in model['features']])-model['mean'])/model['scale']
                    standard[~np.array(model['active'])]=0
                    logits.append(float(np.dot(standard,model['coefficients'])+model['intercept']))
                self.assertAlmostEqual(logits[0]-logits[1],value['alpha_minus_GGR_logit'],places=11)
                self.assertIsNone(value['absolute_class']);self.assertIsNone(value['probability'])


if __name__=='__main__':unittest.main()
