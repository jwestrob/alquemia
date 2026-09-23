"""Actual saved source/proposal geometry, never synthetic molecular output."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import omitted_context_contacts as audit
from affordable_common import read_json,verify,xyz

RESULT=ROOT/'workspaces/omitted_context_contacts_20260923/RESULT_v1.json'

class ActualOmittedContacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=read_json(RESULT);cls.inputs=read_json(verify(cls.data['inputs']))
        cls.parent=read_json(verify(cls.inputs['audit']))

    def test_exact_fixed_population_and_honest_graph_coverage(self):
        expected={r['pair_id'] for r in self.parent['rows'] if not r['pool_reuse']}
        self.assertEqual(len(expected),55);self.assertEqual(expected,{r['pair_id'] for r in self.data['rows']})
        self.assertEqual(len({r['case_id'] for r in self.data['rows']}),41)
        for r in self.data['rows']:
            for v in r['representations'].values():
                self.assertEqual(v['status'],'available')
                self.assertEqual(v['graph']['raw_heavy_max_difference_A'],0)
                self.assertEqual(v['graph']['prepared_only_heavy_atoms'],[])
                self.assertEqual(v['PQQ_heavy_atoms_fully_mapped_and_fixed'],24)

    def test_actual_coordinate_distances_not_only_cached_counts(self):
        # Existing A8 sample1 and A0A3 sample2 are members of the declared55.
        targets=[next(r for r in self.data['rows'] if r['protein_id'].startswith(name) and 'sample-'+sample in r['case_id'])
                 for name,sample in [('a8r3s4','1'),('a0a3','2')]]
        for r in targets:
            for v in r['representations'].values():
                graph=audit.SourceGraph(verify(v['graph']['raw_source']),verify(self.parent['config']['topology']))
                source_by_label={audit.label(meta):graph.atoms[k] for k,meta in graph.meta.items()}
                mapping=read_json(verify(v['source_preparation']))['mapping']['source_to_qm']
                query_indices={audit.label(x['source']):x['qm_index'] for x in mapping if 'source' in x}
                for name,value in v['contacts'].items():
                    nearest=value['moving_donor_subset']['nearest'];atoms=xyz(verify(v['coordinates'][name]['coordinate']))
                    q=np.array(atoms[query_indices[nearest['query']]][1:]);o=source_by_label[nearest['outside']].pos
                    actual=float(np.linalg.norm(q-[o.x,o.y,o.z]))
                    self.assertAlmostEqual(actual,nearest['distance_A'],places=12)
                    self.assertNotIn(nearest['outside'],query_indices)

    def test_three_bond_exclusion_on_actual_glutamate_graph(self):
        source=next(r['source'] for r in self.parent['rows'] if not r['pool_reuse'])
        state=audit.parent_state(source['original_core'],self.parent['config']['topology'],require_endpoint_receipts=False)
        graph=state['graph'];res=next(r for r in graph.residues.values() if r.canonical_resname=='GLU')
        adjacency={k:set() for k in graph.atoms}
        for a,b in graph.edges:adjacency[a].add(b);adjacency[b].add(a)
        seen=audit.within_three(graph.key(res.key,'CD'),adjacency,graph)
        self.assertIn(graph.key(res.key,'CA'),seen)
        self.assertNotIn(graph.key(res.key,'N'),seen)

    def test_threshold_lists_and_outside_stationarity_use_real_results(self):
        fixed=read_json(RESULT.parent/'OUTSIDE_FIXED_CHECK.json')
        self.assertEqual(fixed['checked_context_geometries'],330)
        self.assertEqual(fixed['maximum_unscored_physical_heavy_displacement_A'],0)
        for r in self.data['rows']:
            for v in r['representations'].values():
                for geom in v['contacts'].values():
                    for subset in geom.values():
                        self.assertEqual(subset['below_3_5_A'],len(subset['contacts_below_3_5_A']))
                        self.assertEqual(subset['new_below_2_0_A'],len(subset['new_severe_contacts']))
                        self.assertTrue(all(c['distance_A']<3.5 for c in subset['contacts_below_3_5_A']))
                        self.assertTrue(all(c['distance_A']<2<=c['origin_distance_A'] for c in subset['new_severe_contacts']))
        self.assertEqual(self.data['new_molecular_calls'],0)
        self.assertFalse(self.data['scores_or_classifications_changed'])

if __name__=='__main__':unittest.main()
