"""Real common8 parent/geometry replay; no molecular energies or forces."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
import numpy as np
import openmm as mm
from openmm import app, unit
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, xyz
import scaffold_parent_prepare as parent
import scaffold_environment_inventory as old

class RealParentPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=read_json(ROOT/'workspaces/collective_scaffold_20260922/parents_v1/manifest.json')
        cls.rows=cls.manifest['cases']

    def test_all_declared_sources_targets_and_no_molecular_evaluation(self):
        original=read_json(verify(self.manifest['inputs']))
        self.assertEqual([r['case_id'] for r in self.rows],[r['case_id'] for r in original['cases']])
        self.assertEqual(self.manifest['prepared_cases'],8)
        self.assertEqual(self.manifest['prepared_targets'],24)
        for k in ('energy_evaluations','force_evaluations','molecular_searches','new_terminal_repairs'):
            self.assertEqual(self.manifest[k],0)
        for row in self.rows:
            self.assertEqual(set(row['targets']),set(parent.TARGETS))
            self.assertEqual(row['origins']['La']['charge']-row['origins']['Ca']['charge'],1)
            self.assertEqual(read_json(verify(row['context_maps']['Ca'])),read_json(verify(row['context_maps']['La'])))

    def test_every_source_protein_atom_and_hydrogen_is_preserved(self):
        for row in self.rows:
            source=app.PDBFile(str(verify(row['source'])))
            pos=np.asarray(source.positions.value_in_unit(unit.angstrom))
            actual={old.atom_id(a):(a.element.symbol,pos[a.index]) for a in source.topology.atoms() if a.residue.name in old.STANDARD}
            atoms=read_json(verify(row['parent']['atoms']));p=np.asarray(read_json(verify(row['parent']['positions_A'])))
            for a in atoms:
                if a['source_id'] not in actual:
                    self.assertIn(row['case_id'],('1H4I','4MAE'));self.assertEqual(a['name'],'OXT');continue
                sym,q=actual.pop(a['source_id']);self.assertEqual(sym,a['element'])
                np.testing.assert_array_equal(p[a['parent_index']],q)
            self.assertFalse(actual)
            self.assertEqual(read_json(verify(row['parent']['source_atom_ids'])),[a['source_id'] for a in atoms])
            state=read_json(verify(row['parent']['state']))
            self.assertFalse(state['source_H_moved']);self.assertFalse(state['global_archived_H_normalization_reused'])

    def test_parameterized_system_and_archived_crystal_reuse(self):
        for row in self.rows:
            system=mm.XmlSerializer.deserialize(verify(row['parent']['system']).read_text())
            self.assertEqual(system.getNumParticles(),row['parent_atom_count']);self.assertEqual(system.getNumConstraints(),0)
            self.assertEqual(bool(row['parent']['archived_reuse']),row['case_id'] in ('1H4I','4MAE'))
            nb=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
            self.assertEqual(nb.getNonbondedMethod(),mm.NonbondedForce.NoCutoff)
            q=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
            self.assertAlmostEqual(q,round(q),places=8)

    def test_all_targets_reconstruct_exact_archived_context_and_preserve_parent_bonds(self):
        for row in self.rows:
            top=read_json(verify(row['parent']['topology']));mapping=read_json(verify(row['context_parent_mapping']))
            g=read_json(verify(row['context_maps']['Ca']))['context'];q0=np.asarray(read_json(verify(row['parent']['positions_A'])))
            mobile=read_json(verify(row['mobile_set']));b=np.asarray([x['indices'] for x in top['bonds']]);length=np.linalg.norm(q0[b[:,0]]-q0[b[:,1]],axis=1)
            for target in row['targets'].values():
                p=np.asarray(read_json(verify(target['parent_positions_A'])))
                np.testing.assert_allclose(np.linalg.norm(p[b[:,0]]-p[b[:,1]],axis=1),length,atol=1e-12,rtol=0)
                np.testing.assert_array_equal(p[mobile['frozen_parent_indices']],q0[mobile['frozen_parent_indices']])
                context=parent.reconstruct_context(p,mapping,g)
                actual=np.asarray([a[1:] for a in xyz(verify(target['actual_archived_coordinate']))])
                np.testing.assert_allclose(context,actual,atol=1e-12,rtol=0)
                ca,la=[xyz(verify(target['context_coordinates'][z])) for z in ('Ca','La')]
                self.assertEqual(ca[0][0],'Ca');self.assertEqual(la[0][0],'La');self.assertEqual(ca[1:],la[1:])
                np.testing.assert_allclose([a[1:] for a in ca],actual,atol=1e-12,rtol=0)

    def test_mobile_shell_is_q0_complete_residues_and_one_actual_neighbor_layer(self):
        for row in self.rows:
            top=read_json(verify(row['parent']['topology']));mobile=read_json(verify(row['mobile_set']))
            p=np.asarray(read_json(verify(row['parent']['positions_A'])));core=np.asarray(mobile['canonical_core_heavy_positions_A'])
            atoms=top['atoms'];expected=set()
            for r in top['residues']:
                ids=[i for i in r['parent_indices'] if atoms[i]['element']!='H']
                if np.linalg.norm(p[ids,None]-core[None],axis=2).min()<=8:expected.add(r['residue_index'])
            expected|={atoms[i]['residue_index'] for i in mobile['core_parent_indices']};seed=expected.copy()
            for b in top['bonds']:
                a,c=[atoms[i] for i in b['indices']]
                if a['residue_index']!=c['residue_index'] and {a['name'],c['name']}=={'C','N'}:
                    if a['residue_index'] in seed:expected.add(c['residue_index'])
                    if c['residue_index'] in seed:expected.add(a['residue_index'])
            self.assertEqual(sorted(expected),mobile['mobile_residue_indices'])
            self.assertEqual([a['parent_index'] for a in atoms if a['residue_index'] in expected],mobile['mobile_parent_indices'])

    def test_role_targets_fixed_inventory_and_no_nonprotein_parameters(self):
        for row in self.rows:
            donors=read_json(verify(row['donors']));names={}
            for d in donors['atoms']:names.setdefault(d['role'],[]).append(d['source_id'].split('/')[-1])
            self.assertEqual(names['anchor_glutamate'],['OE1','OE2']);self.assertEqual(names['anchor_asparagine'],['OD1'])
            if 'extra_acidic_ligand_homolog' in names:self.assertEqual(names['extra_acidic_ligand_homolog'],['OD1','OD2'])
            self.assertNotIn('catalytic_aspartate',names)
            inventory=read_json(verify(row['nonprotein_inventory']));self.assertFalse(inventory['FF_parameters_assigned']);self.assertEqual(inventory['source_water_count'],0)
            self.assertEqual(len([a for a in inventory['source_atoms'] if a['resname']=='PQQ']),24)
            for t in row['targets'].values():
                p=np.asarray(read_json(verify(t['parent_positions_A'])))
                np.testing.assert_array_equal(t['donor_positions_A'],p[[d['parent_index'] for d in donors['atoms']]])

    def test_corrupted_real_map_source_and_caps_fail_explicitly(self):
        row=self.rows[0];g=read_json(verify(row['context_maps']['Ca']))['context'];top=read_json(verify(row['parent']['topology']));p=np.asarray(read_json(verify(row['parent']['positions_A'])))
        bad=deepcopy(g);l=next(l for l in bad['core_links'] if l[1]=='cap');bad['core_links'].remove(l)
        with self.assertRaisesRegex(InvalidArtifact,'cut bond/cap mismatch'):parent.context_map(bad,top,p)
        bad=deepcopy(g);i=next(i for i,m in enumerate(bad['source_atom_metadata']) if m and m['resname'] in old.STANDARD);bad['positions_A'][i][0]+=.01
        with self.assertRaisesRegex(InvalidArtifact,'source coordinate mismatch'):parent.context_map(bad,top,p)

if __name__=='__main__':unittest.main()
