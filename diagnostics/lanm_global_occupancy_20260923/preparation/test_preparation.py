"""Actual full-chain fixture checks; no manufactured scientific outputs."""
import unittest
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import read_json, verify, xyz, record
from lanm_global_occupancy_prepare import SOURCES, OCCUPANCIES
from lanm_series_followup import heavy_atoms

MANIFEST = ROOT/'workspaces/lanm_global_occupancy_20260923/prepared_v1/manifest.json'

class PreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = read_json(MANIFEST)

    def test_exact_scope_and_actual_pins(self):
        m = self.m
        self.assertEqual([s['source_id'] for s in m['sources']], list(SOURCES))
        self.assertEqual([s['state_id'] for s in m['states']],
                         [sid+'__'+occ for sid in SOURCES for occ in OCCUPANCIES])
        self.assertTrue(all(s['status'] == 'prepared' for s in m['states']))
        for p in [m['agreement'],m['preparation_agreement'],m['parameter_export'],
                  *m['implementation'].values(),*m['forcefields'].values()]: verify(p)
        self.assertEqual([m[k] for k in ('molecular_calls','new_protonations','geometry_optimizations')],[0,0,0])

    def test_full_source_heavy_coordinates_and_peptide_topology(self):
        for s in self.m['sources']:
            raw = heavy_atoms(verify(s['crystal']))
            pp = heavy_atoms(verify(s['protonated_source']))
            self.assertEqual({k:v for k,v in raw.items() if k[0]=='A'},
                             {k:v for k,v in pp.items() if k[0]=='A'})
            data = read_json(verify(s['atoms'])); atoms = data['atoms']
            for a in atoms:
                if a['element'] != 'H':
                    np.testing.assert_allclose(a['xyz_A'],a['source']['xyz_A'],rtol=0,atol=1e-12)
            n = s['protein_atoms']; graph = {i:set() for i in range(n)}
            for b in data['bonds']:
                i,j=b['indices']
                if b['kind']=='protein_covalent':graph[i].add(j);graph[j].add(i)
            seen={0}; stack=[0]
            while stack:
                for j in graph[stack.pop()]-seen:seen.add(j);stack.append(j)
            self.assertEqual(len(seen),n)
            self.assertEqual(s['peptide_bonds'],s['protein_residues']-1)
            self.assertEqual(s['missing_heavy_added'],0)
            self.assertEqual(s['caps_added'],0)

    def test_hydrogens_and_fixed_water_inventory(self):
        for s in self.m['sources']:
            data=read_json(verify(s['atoms'])); atoms=data['atoms']; degree=[0]*len(atoms)
            for b in data['bonds']:
                i,j=b['indices'];degree[i]+=1;degree[j]+=1
                if 'H' in (atoms[i]['element'],atoms[j]['element']):
                    dist=np.linalg.norm(np.array(atoms[i]['xyz_A'])-atoms[j]['xyz_A'])
                    self.assertAlmostEqual(dist,b['forcefield_equilibrium_A'],places=12)
            self.assertTrue(all(degree[i]==1 for i,a in enumerate(atoms) if a['element']=='H'))
            for state in [r for r in self.m['states'] if r['source_id']==s['source_id']]:
                mapping=read_json(verify(state['mapping']))
                self.assertEqual(mapping['atoms'][:len(atoms)],atoms)
                self.assertEqual(mapping['bonds'],data['bonds'])
                self.assertEqual(state['source_water_inventory'],s['water_inventory'])
                self.assertGreaterEqual(s['minimum_all_atom_distance_A'],.5)

    def test_exact_La_Dy_pair_charge_spin_and_parameters(self):
        params=read_json(verify(self.m['parameter_export']))
        for state in self.m['states']:
            mapping=read_json(verify(state['mapping']));atoms=mapping['atoms']
            ca=xyz(verify(state['endpoints']['La']['xyz']));dy=xyz(verify(state['endpoints']['Dy']['xyz']))
            self.assertEqual(len(ca),state['n_atoms'])
            self.assertEqual([a[1:] for a in ca],[a[1:] for a in dy])
            self.assertEqual(len(state['metal_indices']),state['n'])
            self.assertEqual(state['charge'],state['protein_charge']+3*state['n'])
            for metal,rows in [('La',ca),('Dy',dy)]:
                ep=state['endpoints'][metal]
                self.assertEqual(ep['physical_multiplicity'],1 if metal=='La' else 1+5*state['n'])
                self.assertEqual(ep['all_electron_count']%2,(ep['physical_multiplicity']-1)%2)
                self.assertEqual(ep['native_valence_electron_count']%2,0)
                self.assertEqual(ep['native_effective_multiplicity'],1)
                self.assertEqual(params['element'][metal]['refocc'],[1.,1.,1.])
                self.assertEqual(params['element'][metal]['shells'],['5d','6s','6p'])
                self.assertEqual([i for i,a in enumerate(rows) if a[0]==metal],state['metal_indices'])
            self.assertEqual([a[0] for i,a in enumerate(ca) if i not in state['metal_indices']],
                             [a[0] for i,a in enumerate(dy) if i not in state['metal_indices']])

    def test_real_EF4_source_identity_and_removed_ions(self):
        expected={'Hans_8DQ2':'Na','Hans_8FNR':'Dy','Mex_8FNS':'Nd'}
        for s in self.m['sources']:
            self.assertEqual(s['sites']['EF4']['source']['element'],expected[s['source_id']])
        for state in self.m['states']:
            self.assertEqual(len(state['removed_source_ions']),4-state['n'])
            self.assertEqual(state['EF4_Na_replacement_hypothesis'],
                             state['source_id']=='Hans_8DQ2' and state['n']==4)
            self.assertFalse(read_json(verify(state['mapping']))['coordination_bonds_included'])

if __name__=='__main__':unittest.main()
