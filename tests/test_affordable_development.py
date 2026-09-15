"""Real pinned artifact regressions; no synthetic scientific energies or proteins."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import (record, verify, read_json, energy, contrast, classify_raw,
                               corrected, paired, cache_key, InvalidArtifact)
from affordable_peptide import SourceGraph

BASE = ROOT / 'diagnostics/pqq_pmdh_fixed_core_calibration_20260914'
GGR = ROOT / 'diagnostics/nonpqq_direct_site_benchmark_20260915/prepared/ggr_1glg/GGR/ggr_1glg_GGR_carve_manifest.json'
TOPOLOGY = Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/site-packages/openmm/app/data/amber19/protein.ff19SB.xml')


class ArchivedRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.release = read_json(BASE / 'result.json')
        cls.holdout = read_json(BASE / 'reserved_crystal_holdout/result/holdout_result.json')

    def test_all_released_energies_and_bands(self):
        for payload in (self.release, self.holdout):
            for row in payload['scores']:
                endpoints = {m: energy(verify(a['output'])) for m, a in row['artifacts'].items()}
                for m in endpoints:
                    self.assertAlmostEqual(endpoints[m], row['artifacts'][m]['energy_hartree'], places=10)
                score = contrast(endpoints['Ca'], endpoints['La'], self.release['aquo_reporting_gauge']['delta_E_aquo_hartree'])
                self.assertAlmostEqual(score['R_kcal_mol'], row['R_kcal_mol'], places=7)
                self.assertAlmostEqual(score['S_kcal_mol'], row['S_aquo_gauge_kcal_mol'], places=7)
                self.assertIn(classify_raw(score['R_kcal_mol'], self.release, payload['protocol_id']), ('Ca-supported','Ln-supported'))

    def test_reference_and_threshold_do_not_transfer(self):
        row = self.release['scores'][0]
        score = contrast(row['artifacts']['Ca']['energy_hartree'],row['artifacts']['La']['energy_hartree'])
        self.assertIsNone(score['S_kcal_mol'])
        self.assertEqual(classify_raw(score['R_kcal_mol'], self.release, 'new_model'), 'uncalibrated_protocol')
        self.assertIsNone(corrected(score, None, None)['S_env_kcal_mol'])

    def test_endpoint_contribution_sign_with_actual_energies(self):
        a, b = self.holdout['scores']
        # An algebra identity using actual archived differences, not an environmental result.
        delta_ca = (b['artifacts']['Ca']['energy_hartree']-a['artifacts']['Ca']['energy_hartree'])*627.509474
        delta_la = (b['artifacts']['La']['energy_hartree']-a['artifacts']['La']['energy_hartree'])*627.509474
        result = corrected({'S_kcal_mol':a['S_aquo_gauge_kcal_mol']},delta_ca,delta_la)
        self.assertAlmostEqual(result['S_env_kcal_mol'],b['S_aquo_gauge_kcal_mol'],places=7)

    def test_corrupted_real_output_is_not_converged(self):
        op = verify(self.holdout['scores'][0]['artifacts']['Ca']['output'])
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'explicitly_corrupted_real_output.out'
            p.write_text(op.read_text().replace('ORCA TERMINATED NORMALLY','REMOVED BY MALFORMED-INPUT TEST'))
            with self.assertRaises(InvalidArtifact): energy(p)
            p.write_text(op.read_text().replace('SCF CONVERGED','SCF NOT CONVERGED'))
            with self.assertRaises(InvalidArtifact): energy(p)

    def test_complete_configuration_cache(self):
        config = {'source':record(GGR), 'method':'r2SCAN-3c', 'cavity': 'common_chain',
                  'charge_model':'MBIS', 'assembly':'A', 'microstate':'frozen', 'protocol':'development'}
        original=cache_key(config)
        for key in config:
            changed = dict(config); changed[key] = {'explicit_test_mutation':config[key]}
            self.assertNotEqual(original,cache_key(changed))


class RealPeptideGraph(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=read_json(GGR)
        cls.graph=SourceGraph(verify(cls.manifest['source_structure']),TOPOLOGY)

    def test_ggr_retains_actual_nitrogen_and_source_coordinates(self):
        g=self.graph
        r=g.locate({'chain':'A','resnum':140,'resname':'GLN'})
        chosen=g.peptide(r.key)
        n=g.amide_next[r.key]
        self.assertIn(n,chosen)
        self.assertEqual(g.meta[n]['resnum'],141)
        atoms, manifest=g.materialize(chosen)
        self.assertEqual(sorted(x[0] for x in atoms),['C','H','H','H','N','O'])
        self.assertEqual(len(manifest['cut_bonds_and_caps']),2)
        for row in manifest['source_to_qm']:
            if row['kind']=='source':
                key=(row['source']['chain_index'],row['source']['residue_index'],row['source']['atom'])
                self.assertEqual(tuple(row['xyz_A']), tuple(g.atoms[key].pos))

    def test_real_proline_preserves_ring(self):
        g=self.graph
        carbon=next(k for k,n in g.amide_next.items() if g.residues[n[:2]].resname=='PRO')
        chosen=g.peptide(carbon)
        n=g.amide_next[carbon]
        for name in ('N','CA','CB','CG','CD'):
            self.assertIn((*n[:2],name),chosen)
        atoms,m=g.materialize(chosen)
        self.assertEqual(sum(x[0]=='N' for x in atoms),1)
        self.assertTrue(any(b['a']['atom']=='CD' and b['b']['atom']=='N' or b['b']['atom']=='CD' and b['a']['atom']=='N' for b in m['retained_bonds']))

    def test_overlapping_peptides_and_sidechain_merge(self):
        g=self.graph
        r=g.locate({'chain':'A','resnum':140,'resname':'GLN'})
        next_r=g.amide_next[r.key][:2]
        selected=g.peptide(r.key)|g.peptide(next_r)
        selected.update(g.key(r.key,n) for n in ('CB','CG','CD','OE1','NE2'))
        atoms,m=g.materialize(selected)
        source=[str(x['source']) for x in m['source_to_qm'] if x['kind']=='source']
        self.assertEqual(len(source),len(set(source)))
        self.assertEqual(len(atoms),len({tuple(a[1:]) for a in atoms}))
        self.assertIn(g.key(r.key,'CA'),g.close_overlaps(selected))
        self.assertIn(g.key(next_r,'CA'),g.close_overlaps(selected))

    def test_missing_neighbor_fails_in_corrupted_real_structure(self):
        g=self.graph
        r=g.locate({'chain':'A','resnum':140,'resname':'GLN'})
        with tempfile.TemporaryDirectory() as d:
            st=g.structure.clone()
            n=g.amide_next[r.key]
            st[0][n[0]][n[1]].remove_atom('N','*')
            p=Path(d)/'explicitly_corrupted_missing_N.pdb'; st.write_pdb(str(p))
            modified=SourceGraph(p,TOPOLOGY)
            with self.assertRaises(InvalidArtifact): modified.peptide(r.key)

    def test_numbering_and_insertion_codes_are_not_bonds(self):
        g=self.graph
        r=g.locate({'chain':'A','resnum':140,'resname':'GLN'})
        with tempfile.TemporaryDirectory() as d:
            st=g.structure.clone(); n=g.amide_next[r.key]
            st[0][n[0]][n[1]].seqid.num=7141
            st[0][n[0]][n[1]].seqid.icode='B'
            p=Path(d)/'real_structure_identifier_test.pdb'; st.write_pdb(str(p))
            modified=SourceGraph(p,TOPOLOGY)
            self.assertEqual(modified.amide_next[r.key],n)
            self.assertEqual(modified.meta[n]['insertion_code'],'B')

    def test_incompatible_altlocs_and_broken_real_peptide_fail(self):
        g=self.graph;r=g.locate({'chain':'A','resnum':140,'resname':'GLN'});n=g.amide_next[r.key]
        with tempfile.TemporaryDirectory() as d:
            st=g.structure.clone()
            st[0][r.key[0]][r.key[1]]['C'][0].altloc='A'
            st[0][n[0]][n[1]]['N'][0].altloc='B'
            p=Path(d)/'corrupted_incompatible_conformers.pdb';st.write_pdb(str(p))
            with self.assertRaises(InvalidArtifact): SourceGraph(p,TOPOLOGY).peptide(r.key)
            st=g.structure.clone();st[0][n[0]][n[1]]['N'][0].pos.x+=10
            p=Path(d)/'corrupted_broken_peptide.pdb';st.write_pdb(str(p))
            with self.assertRaises(InvalidArtifact): SourceGraph(p,TOPOLOGY).peptide(r.key)

    def test_all_six_repaired_pair_invariants(self):
        inventory=read_json(ROOT/'diagnostics/affordable_challenger_20260915/repair_inventory_verified.json')
        self.assertEqual(len(inventory),6)
        for r in inventory:
            m=read_json(verify(r['repaired']))
            out=m['outputs']
            paired(verify(out['La']['xyz']),verify(out['Ca']['xyz']),out['La']['charge'],out['Ca']['charge'])
            self.assertEqual(m['heavy_coordinate_max_displacement_A'],0)
            self.assertEqual(out['La']['all_electron_count_for_parity']%2,0)
            self.assertEqual(out['Ca']['all_electron_count_for_parity']%2,0)

    def test_physical_cap_jacobians_on_real_bond(self):
        from affordable_response import cap_jacobians, bounded_response
        g=self.graph;r=g.locate({'chain':'A','resnum':140,'resname':'GLN'})
        n=g.amide_next[r.key]
        x=np.array(tuple(g.atoms[n].pos));y=np.array(tuple(g.atoms[g.key(n[:2],'CA')].pos))
        ja,jb=cap_jacobians(x,y,1.01)
        def cap(a,b): return a+1.01*(b-a)/np.linalg.norm(b-a)
        for axis in range(3):
            step=np.eye(3)[axis]*1e-6
            np.testing.assert_allclose(ja[:,axis],(cap(x+step,y)-cap(x-step,y))/2e-6,atol=1e-8)
            np.testing.assert_allclose(jb[:,axis],(cap(x,y+step)-cap(x,y-step))/2e-6,atol=1e-8)
        np.testing.assert_allclose(ja+jb,np.eye(3),atol=1e-14)
        self.assertIsNone(bounded_response(None,None)['relaxation_energy_kcal_mol'])


class PilotIntegration(unittest.TestCase):
    def test_actual_solver_accounting_is_charged_to_recovery(self):
        from affordable_solver import prior_allocated_cost
        p=ROOT/'workspaces/affordable_challenger_20260915/solver_terminal_receipt.json'
        if not p.exists(): self.skipTest('terminal solver accounting not yet available')
        cost,rows=prior_allocated_cost([p])
        self.assertEqual(cost,688)
        with self.assertRaises(InvalidArtifact): prior_allocated_cost([p,p])

    def test_missing_mbis_charge_table_in_corrupted_real_output(self):
        from affordable_environment import mbis_charges
        m=read_json(ROOT/'workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json')
        t=m['tasks'][0];source=Path(t['output_path'])
        if not source.exists(): self.skipTest('real MBIS endpoint unavailable')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_charge_header.out'
            import re
            p.write_text(re.sub(r'ATOM\s+CHARGE\s+POPULATION\s+SPIN','REMOVED_CHARGE_HEADER',source.read_text()))
            with self.assertRaises(InvalidArtifact): mbis_charges(p,verify(t['xyz']),t['charge'])

    def test_single_point_is_not_an_analytic_gradient(self):
        from affordable_response import extract
        row=read_json(BASE/'reserved_crystal_holdout/result/holdout_result.json')['scores'][0]['artifacts']['La']
        with self.assertRaises(InvalidArtifact):
            extract(Path('/no_gradient_fixture'),verify(row['output']),verify(row['input']),verify(row['xyz']))

    def test_real_environment_pairs_and_common_partition_cavity(self):
        from affordable_state import validate_skeleton_pair
        base=ROOT/'workspaces/affordable_challenger_20260915/environment'
        hashes=[]
        for radius in ('33','36'):
            la,ca=[read_json(base/f'1h4i_qm{radius}_{metal}'/'skeleton.json') for metal in ('La','Ca')]
            hashes.append(validate_skeleton_pair(la,ca)['physical_boundary_hash'])
            ca['environment_atoms'][0]['charge_e']+=.01  # corrupted real prepared state
            with self.assertRaises(InvalidArtifact): validate_skeleton_pair(la,ca)
        self.assertEqual(hashes[0],hashes[1])

    def test_pilot_inputs_are_finite_and_frozen(self):
        from affordable_workflow import dry_run
        r=dry_run(ROOT/'workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json')
        self.assertEqual(r['tasks'],6)
        self.assertFalse(r['budget_enforced'])

    def test_historical_budget_does_not_stop_agreed_tasks(self):
        from affordable_workflow import dry_run
        source=ROOT/'workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json'
        # Keep the real coordinate/input records; only corrupt the budget.
        m=read_json(source);m['budget']['max_endpoint_evaluations_including_retries']=5
        p=source.parent/'explicitly_corrupted_budget_test.json'
        try:
            p.write_text(__import__('json').dumps(m))
            self.assertEqual(dry_run(p)['tasks'],6)
            self.assertFalse(dry_run(p)['budget_enforced'])
        finally: p.unlink(missing_ok=True)

    def test_actual_mbis_endpoint_integration(self):
        from affordable_environment import mbis_charges
        m=read_json(ROOT/'workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json')
        if not all(Path(t['output_path']).exists() for t in m['tasks']):
            self.skipTest('scientific integration unrun: manifested ORCA pilot not complete')
        for t in m['tasks']:
            q=mbis_charges(t['output_path'],t['xyz']['path'],t['charge'])
            self.assertLess(abs(q['sum_e']-t['charge']),1e-4)

    def test_actual_apbs_identity_integration(self):
        p=ROOT/'workspaces/affordable_challenger_20260915/solver_completion/identity/1h4i_qm33_La/result.json'
        if not p.exists(): self.skipTest('scientific integration unrun: APBS identity result unavailable')
        r=read_json(p)
        self.assertLessEqual(abs(r['components']['delta_U_kcal_mol']),.01)

    def test_split_inputs_preserve_real_charging_blocks(self):
        from affordable_solver import split_charging_input
        import re
        p=ROOT/'workspaces/affordable_challenger_20260915/solver_recovery/1h4i_qm33_La/primary/calculation/transfer.in'
        original=p.read_text();parts=split_charging_input(original)
        self.assertEqual(len(parts),6)
        for name,text in parts.items():
            block=re.search(r'^elec name '+name+r'\n(.*?)^end\s*$',original,re.M|re.S)[1]
            self.assertIn('elec name '+name+'\n'+block+'end\n',text)
            self.assertEqual(text.split('elec name ')[0],original.split('elec name ')[0])
            self.assertEqual(text.count('print elecEnergy'),1)

    def test_real_completion_states_match_frozen_schedule(self):
        from affordable_solver import completion_preflight
        p=ROOT/'workspaces/affordable_challenger_20260915/solver_completion/completion_manifest.json'
        if not p.exists(): self.skipTest('completion preparation not yet available')
        m=completion_preflight(p)
        self.assertEqual(len(m['entries']),18);self.assertEqual(len(m['tasks']),102)
        self.assertEqual(sum('cached_result' in e for e in m['entries']),1)
        schedule=read_json(verify(m['previous_schedule']))
        self.assertEqual({e['label']:e['state_cache_key'] for e in m['entries']},{e['label']:e['state_cache_key'] for e in schedule['entries']})
        old=read_json(ROOT/'workspaces/affordable_challenger_20260915/solver_recovery/1h4i_qm33_La/primary/calculation/apbs_manifest.json')
        repeated=read_json(p.parent/'repeat/1h4i_qm33_La/calculation/apbs_manifest.json')
        self.assertAlmostEqual(old['direct_coulomb_kcal_mol'],repeated['direct_coulomb_kcal_mol'],places=10)
        for name in old['pqr']: self.assertEqual(old['pqr'][name]['sha256'],repeated['pqr'][name]['sha256'])


if __name__=='__main__': unittest.main()
