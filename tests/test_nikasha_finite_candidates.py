"""Finite-candidate admission against completed molecular artifacts; no fake energies."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import nikasha_finite_candidates as finite
from affordable_common import InvalidArtifact,read_json,record,write_new


class FiniteCandidates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ip=ROOT/'diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json'
        if not ip.exists():raise unittest.SkipTest('actual pinned common8 unavailable')
        cls.inputs=read_json(ip)
        cls.tmp=tempfile.TemporaryDirectory(prefix='finite_candidate_test_',dir=ROOT/'workspaces')
        cls.root=Path(cls.tmp.name)
        cases=[]
        for row in cls.inputs['cases']:
            _,_,_,tasks=finite.sources(cls.inputs,row['case_id'])
            kin=finite.Kinematics(finite.pool.pinned(tasks['Ca']['mapping'])['context'])
            cases.append({'case_id':row['case_id'],'status':'prepared',
                          'candidates':[{'id':'replayed_origin','full_q':[0.]*len(kin.modes)}]})
        cls.spec={'branch':'solvent_guided','inputs':record(ip),'agreement':cls.inputs['agreement'],
                  'coordinate_limits':{'maximum_angle_radian':.8,'maximum_heavy_displacement_A':.8},
                  'maximum_candidates_per_case':4,'cases':cases}
        cls.sp=cls.root/'spec.json';write_new(cls.sp,cls.spec)
        finite.prepare(cls.sp,cls.root/'identity')
        cls.mp=cls.root/'identity/manifest.json'

    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    def test_real_origins_deduplicate_and_replay_unchanged_contrasts(self):
        result=finite.validate(self.mp)
        self.assertEqual(result['new_MACE_cells'],0);self.assertEqual(result['new_GFN2_calls'],0)
        out=self.root/'identity_collection.json'
        finite.pool.collect(self.mp,out)
        data=read_json(out);self.assertEqual(data['available'],8)
        for c in data['cases']:
            self.assertEqual(c['pool'],c['prior_pool'])
            self.assertEqual(c['aliases']['replayed_origin']['representative'],'origin')

    def test_corrupted_real_coordinates_and_unselected_motion_rejected(self):
        cid=self.spec['cases'][0]['case_id'];_,_,_,ts=finite.sources(self.inputs,cid)
        original=self.spec['cases'][0]['candidates'][0]
        bad=copy.deepcopy(original);bad['full_q'][0]=.1
        self.assertNotIn(0,ts['Ca']['active_indices'])
        with self.assertRaises(InvalidArtifact):finite.physical_candidate(bad,ts,self.spec['coordinate_limits'])
        bad=copy.deepcopy(original);bad['full_q'][ts['Ca']['active_indices'][0]]=.9
        with self.assertRaises(InvalidArtifact):finite.physical_candidate(bad,ts,self.spec['coordinate_limits'])
        bad=copy.deepcopy(original)
        atoms=finite.xyz(finite.verify(ts['Ca']['xyz']));atoms[1]=(atoms[1][0],atoms[1][1]+.01,*atoms[1][2:])
        path=self.root/'explicitly_corrupted_real_coordinate.xyz';finite.write_xyz(path,atoms)
        bad['coordinate']=record(path)
        with self.assertRaises(InvalidArtifact):finite.physical_candidate(bad,ts,self.spec['coordinate_limits'])

    def test_unavailable_proposal_does_not_report_baseline_as_branch_success(self):
        spec=copy.deepcopy(self.spec);spec['cases'][0].update(status='unavailable',reason='explicit_missing_proposal_fixture',candidates=[])
        sp=self.root/'missing.json';write_new(sp,spec);finite.prepare(sp,self.root/'missing')
        out=self.root/'missing_collection.json';finite.pool.collect(self.root/'missing/manifest.json',out)
        data=read_json(out);self.assertEqual(data['available'],7)
        self.assertEqual(data['cases'][0]['pool']['status'],'unavailable')
        self.assertEqual(data['cases'][0]['prior_pool']['status'],'available')

    def test_source_recipe_and_archived_energy_corruption_rejected(self):
        for key in ('GFN2_maxiter','prior_energy'):
            m=read_json(self.mp)
            if key=='GFN2_maxiter':m[key]=125
            else:m['cases'][0]['matrix']['Ca']['origin']['components']['MACE_eV']+=.01
            path=self.root/('corrupt_'+key+'.json');write_new(path,m)
            with self.assertRaises(InvalidArtifact):finite.validate(path)

    def test_repeated_and_missing_case_rejected(self):
        for rows in (self.spec['cases'][:-1],self.spec['cases']+[self.spec['cases'][0]]):
            spec={**self.spec,'cases':rows};p=self.root/('denominator_'+str(len(rows))+'.json');write_new(p,spec)
            with self.assertRaises(InvalidArtifact):finite.load_spec(p)


if __name__=='__main__':unittest.main()
