"""Real archived fixtures; software/algebra checks, not new scientific evidence."""
import copy
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import site_classifier as classifier


class ArchivedClassifierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT/'workspaces/site_classifier_20260918/features_v1/features.json'
        if not path.exists():
            raise unittest.SkipTest('real extracted benchmark fixture unavailable')
        cls.data = json.loads(path.read_text())
        cls.rows = cls.data['rows']

    def test_geometry_rigid_transform_on_real_proteins(self):
        rotation = np.array([[.36,-.48,.8],[.8,.6,0],[-.48,.64,.6]])
        np.testing.assert_allclose(rotation@rotation.T,np.eye(3),atol=1e-15)
        for case in ('a0a3f2yly8-pqq-la_model','GGR_1GLG'):
            row = next(r for r in self.rows if r['case_id']==case)
            original = classifier.load_pin(row['preparation']); moved = copy.deepcopy(original)
            for atom in moved['physical_atoms']:
                atom['xyz_A'] = (np.array(atom['xyz_A'])@rotation + [40,-20,11]).tolist()
            actual,_,_,_ = classifier.geometry(moved)
            for name in classifier.STRUCTURE:
                self.assertAlmostEqual(actual[name],row['features'][name],places=10)

    def test_corrupted_real_coordinates_fail(self):
        original = classifier.load_pin(self.rows[0]['preparation'])
        corrupted = copy.deepcopy(original)
        corrupted['physical_atoms'][0]['xyz_A'][0] = float('nan')
        with self.assertRaises(classifier.InvalidArtifact):classifier.geometry(corrupted)

    def test_logistic_gradient_on_real_feature_values(self):
        rows = [r for r in self.rows if r['target']=='PQQ_functional_class' and r['role']=='calibration']
        x = np.array([[r['features'][k] for k in classifier.ARMS['DFT_structure_MACE']] for r in rows])
        weights = classifier.balanced_weights(rows); mean = weights@x
        scale = np.sqrt(weights@((x-mean)**2)); x = (x-mean)/np.where(scale>1e-12,scale,1.)
        y = np.array([r['label']=='La' for r in rows],dtype=float)
        theta = np.zeros(x.shape[1]+1); _,analytic = classifier.objective(theta,x,y,weights)
        numerical = []
        for i in range(len(theta)):
            plus = theta.copy();minus = theta.copy();plus[i]+=1e-5;minus[i]-=1e-5
            numerical.append((classifier.objective(plus,x,y,weights)[0]-classifier.objective(minus,x,y,weights)[0])/2e-5)
        np.testing.assert_allclose(analytic,numerical,atol=2e-10,rtol=1e-8)

    def test_structure_replicates_do_not_multiply_group_weight(self):
        rows = [r for r in self.rows if r['target']=='direct_site_affinity_direction']
        duplicated = rows + [copy.deepcopy(r) for r in rows if r['biological_group']=='GGR_MglB']
        def totals(data):
            weights = classifier.balanced_weights(data)
            return {g:sum(w for w,r in zip(weights,data) if r['biological_group']==g)
                    for g in {r['biological_group'] for r in data}}
        for group,value in totals(rows).items():
            self.assertAlmostEqual(value,totals(duplicated)[group],places=14)
        with self.assertRaises(classifier.InvalidArtifact):
            classifier.balanced_weights([r for r in rows if r['label']=='Ca'])

    def test_duplicates_and_family_members_stay_grouped(self):
        groups,pairs = classifier.group_rows(self.rows)
        for a in self.rows:
            for b in self.rows:
                if a['biological_group']==b['biological_group'] and a['target']==b['target']:
                    self.assertEqual(groups[a['case_id']],groups[b['case_id']])
        for pair in pairs:
            if pair['identity_relative_longer']>=.5:
                self.assertEqual(groups[pair['a']],groups[pair['b']])

    def test_saved_real_fit_replay_and_scope_guard(self):
        work = ROOT/'workspaces/site_classifier_20260918'
        result = work/'evaluation_v1/result.json'
        if not result.exists():self.skipTest('actual completed classifier fit unavailable')
        data = json.loads((work/'features_v2/features.json').read_text())
        fitted = json.loads(result.read_text()); rows = {r['case_id']:r for r in data['rows']}
        for target,report in fitted['results'].items():
            for fold in report['folds']:
                self.assertFalse(set(fold['held_cases']) & set(fold['training_cases']))
            for arm,report in report['arms'].items():
                model = fitted['models'][target][arm]
                for saved in report['training_replay']:
                    row = rows[saved['case_id']]
                    replay = classifier.predict(model,row)
                    self.assertAlmostEqual(replay['logit'],saved['logit'],places=14)
                    self.assertEqual(replay['class'],saved['class'])
                    unlabelled = {k:v for k,v in row.items() if k!='label'}
                    self.assertIsNone(classifier.predict(model,unlabelled)['correct'])
                wrong = copy.deepcopy(row);wrong['DFT_protocol']='explicitly_corrupted_real_protocol'
                with self.assertRaises(classifier.InvalidArtifact):classifier.predict(model,wrong)


if __name__ == '__main__':unittest.main()
