"""Actual prepared La triples and pinned archived origins; no fabricated science."""
from pathlib import Path
from collections import Counter
import copy
import itertools
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new,xyz,paired
import union_triple_preparation as u
import consistent_context as context

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'workspaces/union_triple_preparation_20260923/run_v3'


class TriplePreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (RUN/'REUSE.json').exists():raise unittest.SkipTest('actual completed preparation/reuse inventory absent')
        cls.m=u.validate_selection(RUN/'SELECTION.json');cls.p=read_json(RUN/'PREPARATION.json');cls.r=read_json(RUN/'REUSE.json')
        cls.index={(r['selection_id'],r['case_id']):r for r in cls.p['cases']}

    def test_all100_exact_subsets_and_six_old_missing(self):
        source=read_json(verify(self.m['source_manifest']));groups={}
        for r in source['cases']:
            if r['source_conditioning_metal']=='La' and not r['canonical_coordinate_match']:
                groups.setdefault(r['root_case_id'],[]).append(r['case_id'])
        expected={(g,tuple(sorted(t))) for g,ids in groups.items() for t in itertools.combinations(ids,3)}
        actual={(t['protein_id'],tuple(sorted(t['members']))) for t in self.p['triples']}
        self.assertEqual(actual,expected);self.assertEqual(len(actual),100)
        old=read_json(verify(self.m['inventory']))
        self.assertEqual({tuple(sorted(t['members'])) for t in self.p['triples'] if t['status']!='prepared'},
            {tuple(sorted(t['sources'])) for t in old['triples'] if t['missing_sources']})
        self.assertEqual(Counter(t['status'] for t in self.p['triples']),{'prepared':94,'prior_preparation_unavailable':6})

    def test_real_deduplicated_population_and_pair_invariants(self):
        self.assertEqual(len(self.m['groups']),34);self.assertEqual(len(self.index),126)
        self.assertEqual(sum(not r['is_separate_stress_probe'] for r in self.p['cases']),125)
        for r in self.p['cases']:
            self.assertEqual(r['status'],'prepared');ep=r['representations']['context']['endpoints']
            paired(verify(ep['La']['xyz']),verify(ep['Ca']['xyz']),ep['La']['charge'],ep['Ca']['charge'])
            self.assertEqual(xyz(verify(ep['La']['xyz']))[1:],xyz(verify(ep['Ca']['xyz']))[1:])
            anchor=self.index[r['selection_id'],r['state_anchor_source']]
            self.assertEqual(r['state_key'],anchor['state_key']);self.assertEqual(r['protein_key'],anchor['protein_key'])
            audit=read_json(verify(r['representations']['context']['preparation']))
            self.assertTrue(audit['core_source_coordinates_unchanged']);self.assertTrue(audit['water_inventory_unchanged'])
            self.assertFalse(audit['new_protons_or_waters'])

    def test_exact_vs_changed_graph_and_xyz_matches_actual_tenfold(self):
        old={r['case_id']:r for r in read_json(verify(self.m['tenfold_preparation']))['cases']}
        counts=Counter()
        for r in self.p['cases']:
            prev=old[r['case_id']];ep=r['representations']['context']['endpoints']
            same=r['state_key']==prev['state_key'] and all(context.reusable_state(ep[z],prev['representations']['context']['endpoints'][z]) for z in ('Ca','La'))
            self.assertEqual(same,r['tenfold_comparison']['exact_reuse_eligible'])
            if not r['is_separate_stress_probe']:counts[same]+=1
        self.assertEqual(counts,{True:70,False:55})
        self.assertEqual(sum(t.get('all_members_exact_tenfold',False) for t in self.p['triples']),61)

    def test_stress_probe_uses_only_La_selection_and_not_denominator(self):
        stress=next(r for r in self.p['cases'] if r['is_separate_stress_probe'])
        self.assertEqual(stress['case_id'],u.STRESS);self.assertEqual(stress['source']['source_conditioning_metal'],'Ca')
        self.assertEqual((stress['atom_count'],stress['cap_count'],stress['new_La_charge']),(169,13,-2))
        g=next(g for g in self.m['groups'] if g['selection_id']==stress['selection_id'])
        self.assertTrue(all('__conditioned_La__' in c for c in g['members']))
        self.assertTrue(all(u.STRESS not in t['members'] for t in self.p['triples']))
        self.assertEqual(stress['tenfold_comparison']['atom_count'],193)

    def test_actual_origin_receipts_forces_and_explicit_missing_values(self):
        self.assertEqual(len(self.r['rows']),126)
        for r in self.r['rows']:
            prepared=self.index[r['selection_id'],r['case_id']]
            for z,ep in r['native_endpoints'].items():
                actual=read_json(verify(ep['receipt']));self.assertEqual(actual['energy_eV'],ep['energy_eV'])
                self.assertEqual(actual['forces'],ep['forces']);force=np.load(verify(ep['forces']),allow_pickle=False)
                self.assertEqual(force.shape,(prepared['atom_count'],3));self.assertTrue(np.isfinite(force).all())
            self.assertEqual(set(r['missing_native_endpoints']),{'Ca','La'}-set(r['native_endpoints']))
            for z,media in r['solvent_endpoints'].items():
                for name,ep in media.items():
                    for k in ('manifest','receipt','output'):verify(ep[k])
                    self.assertIsInstance(ep['energy_hartree'],float)
            self.assertIsNone(r['classification']);self.assertEqual(r['candidate_pool_reuse'],'not_audited_requires_exact_physical_maps_and_settings')
        self.assertEqual(self.r['summary']['source_selection_pairs'],125)
        self.assertEqual(self.r['new_molecular_calls'],0);self.assertIsNone(self.r['new_reference'])

    def test_corrupted_real_subset_deletion_is_rejected(self):
        d=read_json(verify(self.m['inventory']));d['triples'].pop()
        with tempfile.TemporaryDirectory(dir=ROOT/'workspaces') as td:
            path=Path(td)/'corrupted_missing_triple.json';write_new(path,d)
            with self.assertRaisesRegex(InvalidArtifact,'all100'):u.selection(path,Path(td)/'must_not_prepare')

    def test_corrupted_real_frozen_membership_is_rejected(self):
        d=copy.deepcopy(self.m);d['groups'][0]['members'].pop()
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'corrupted_membership.json';write_new(path,d)
            with self.assertRaisesRegex(InvalidArtifact,'frozen triple selection'):u.validate_selection(path)


if __name__=='__main__':unittest.main()
