"""Real alpha-lactalbumin network geometry tests; no quantum backend calls."""
from pathlib import Path
import sys
import unittest
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import hydration_network as h


class NetworkFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config=ROOT/'diagnostics/hydration_network_20260918/CONFIG.json'
        if not cls.config.exists():raise unittest.SkipTest('pinned local fixture config missing')
        cfg=h.read_json(cls.config);source=h.read_json(h.verify(cfg['source_config']))
        if not all(Path(i['preparation']['path']).exists() for i in source['parents']):
            raise unittest.SkipTest('real alpha source structures unavailable')
        cls.cfg,cls.source,cls.parents,cls.union=h.discover(cls.config)
        cls.reference=h.xyz(h.verify(cls.source['water_reference_geometry']))

    def test_omitted_polar_neighbors_are_retained_with_complete_amides(self):
        for data in self.parents:
            atoms,mapping,ledger,groups=h.materialize(data,self.union,self.reference)
            present={h.source_id(a['source']) for a in mapping['source_to_qm'] if a['kind']=='source'}
            for contact in data['contacts']:
                self.assertIn(h.source_id(contact['neighbor']),present)
            cuts=mapping['cut_bonds_and_caps']
            self.assertFalse(any({c['retained']['atom'],c['omitted']['atom']}=={'C','N'} for c in cuts))
            self.assertEqual(sum(a['formal_charge'] for a in ledger),-4)
            self.assertEqual(len(present),len([a for a in mapping['source_to_qm'] if a['kind']=='source']))
            for a in mapping['source_to_qm']:
                if a['kind']=='source' and a['source']['element'] not in ('H','D'):
                    np.testing.assert_array_equal(a['source_xyz_A'],a['xyz_A'])

    def test_replica_composition_matches_and_outer_waters_are_distinct(self):
        counts=[]
        for data in self.parents:
            atoms,mapping,ledger,groups=h.materialize(data,self.union,self.reference)
            nonwater=[a for a in mapping['source_to_qm'] if a['kind']!='source' or a['source']['canonical_resname']!='HOH']
            counts.append(h.square.Counter(atoms[a['qm_index']][0] for a in nonwater))
            expected=(2,1) if data['item']['case']=='1F6S' else (3,2)
            self.assertEqual(tuple(sum(w['role']==role for w in groups) for role in ('variable','frozen_outer')),expected)
        self.assertEqual(counts[0],counts[1])

    def test_radial_seed_is_rigid_and_rotates_with_the_real_structure(self):
        data=self.parents[0];atoms,mapping,ledger,groups=h.materialize(data,self.union,self.reference)
        seed=h.radial_seed(atoms,groups,mapping)
        from scipy.spatial.transform import Rotation
        rotation=Rotation.from_rotvec([.3,-.6,.2]).as_matrix();translation=np.array([12.,-8.,3.])
        transformed=[(a[0],*(rotation@np.array(a[1:])+translation)) for a in atoms]
        transformed_seed=h.radial_seed(transformed,groups,mapping)
        np.testing.assert_allclose([a[1:] for a in transformed_seed],
                                   [rotation@np.array(a[1:])+translation for a in seed],atol=1e-12)
        mobile={i for w in groups if w['role']=='variable' for i in w['hydrogen_indices']}
        for i,a in enumerate(atoms):
            if i not in mobile:self.assertEqual(a,seed[i])
        for w in groups:
            before=np.array([atoms[i][1:] for i in w['indices']]);after=np.array([seed[i][1:] for i in w['indices']])
            np.testing.assert_allclose(np.linalg.norm(before[:,None]-before[None,:],axis=2),
                                       np.linalg.norm(after[:,None]-after[None,:],axis=2),atol=1e-13)

    def test_native_constraints_free_only_selected_water_hydrogens(self):
        data=self.parents[1];atoms,mapping,ledger,groups=h.materialize(data,self.union,self.reference)
        text=h.constrained_input(atoms,groups,-1,80)
        mobile={i for w in groups if w['role']=='variable' for i in w['hydrogen_indices']}
        self.assertEqual(len(mobile),6)
        for i in range(len(atoms)):
            self.assertEqual(f'{{ C {i} C }}' in text,i not in mobile)
        self.assertNotIn('NumGrad',text);self.assertNotIn('Freq',text)

    def test_both_real_structures_are_prepared_without_case_filter(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'prepared'
            h.prepare(str(self.config),str(output))
            m=h.read_json(output/'manifest.json')
            self.assertEqual(len(m['tasks']),8)
            self.assertEqual({t['case'] for t in m['tasks']},{'1F6S','6IP9'})
            self.assertIsNone(m['case_filter'])

    def test_tabulated_water_reference_units_and_equilibrium(self):
        import hydration_reference as r
        data=h.read_json(ROOT/'diagnostics/hydration_network_20260918/WATER_THERMOCHEMISTRY.json')
        terms=r.thermochemical_terms(data)
        self.assertAlmostEqual(terms['vapor_pressure_bar'],.0316674874006291,places=12)
        self.assertGreater(terms['harmonic_gas_ZPE_kcal_mol'],13)
        self.assertLess(terms['harmonic_gas_ZPE_kcal_mol'],14)
        self.assertAlmostEqual(terms['gas_thermal_enthalpy_increment_kcal_mol']*4.184,9.904)
        # At coexistence the extra liquid pressure correction must vanish.
        coexist=dict(data,liquid_pressure_bar=terms['vapor_pressure_bar'])
        self.assertEqual(r.thermochemical_terms(coexist)['liquid_pressure_term_kcal_mol'],0)
        with self.assertRaises(h.InvalidArtifact):r.thermochemical_terms(dict(data,temperature_K=310))

    def test_real_water_exchange_ledger_cancels_reference_in_metal_difference(self):
        import hydration_state_analysis as a
        rp=ROOT/'workspaces/hydration_network_20260918/water_reference_v2/reference_1201831.json'
        if not rp.exists():self.skipTest('computed gas-water reference unavailable')
        reference=h.read_json(rp)
        old=h.read_json(ROOT/'workspaces/hydration_square_20260918/repaired_v1/collection_1201801.json')
        rows={r['state_id']:r for r in old['rows']}
        for square in old['squares']:
            full,deleted=rows[square['full_state']],rows[square['deleted_state']]
            ledger={metal:a.exchange_ledger(full['endpoints'][metal]['energy_hartree'],
                    deleted['endpoints'][metal]['energy_hartree'],1,reference) for metal in ('Ca','La')}
            delta=ledger['Ca']['missing_bound_contribution_for_neutral_exchange_kcal_mol']-ledger['La']['missing_bound_contribution_for_neutral_exchange_kcal_mol']
            self.assertAlmostEqual(delta,-square['delta_S_water_addition_kcal_mol'],places=8)
            self.assertIsNone(ledger['Ca']['formation_free_energy_kcal_mol'])
            self.assertIsNone(ledger['La']['bound_state_correction_kcal_mol'])

    def test_real_completed_gas_optimization_trace_has_no_invented_final_gradient(self):
        base=ROOT/'workspaces/hydration_network_20260918/water_reference_v2/gas'
        if not (base/'endpoint.out').exists():self.skipTest('computed gas optimization unavailable')
        text=(base/'endpoint.out').read_text();initial=h.xyz(base/'water.xyz')
        trace=h.optimization_trace(text,initial)
        self.assertEqual(len(trace['analytic_gradient_evaluations']),4)
        self.assertIsNone(trace['final_endpoint_gradient'])
        final=h.xyz(base/'endpoint.runtime.xyz')
        np.testing.assert_allclose(trace['last_printed_coordinates_A'],[a[1:] for a in final],atol=6e-7)
        self.assertNotEqual(trace['analytic_gradient_evaluations'][-1]['energy_hartree'],-76.418935329701)
        # Explicitly corrupted copy of the real printed output, not scientific data.
        damaged=text.replace('CARTESIAN GRADIENT','CARTESIAN GRADIENT_CORRUPTED')
        with self.assertRaises(h.InvalidArtifact):h.optimization_trace(damaged,initial)

    def test_omitted_seed_cannot_be_treated_as_completed_search(self):
        import hydration_state_analysis as analysis
        import json
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'prepared';h.prepare(str(self.config),str(output),case='1F6S')
            mp=output/'manifest.json';manifest=h.read_json(mp)
            # Malformed fixture: real prepared manifest with one requested start lost.
            manifest['tasks']=[t for t in manifest['tasks'] if t['seed']=='source']
            mp.write_text(json.dumps(manifest))
            cp=output/'missing_seed_collection.json'
            cp.write_text(json.dumps({'manifest':h.record(mp),'rows':[dict(task_id=t['task_id'],case=t['case'],
                metal=t['metal'],pattern=t['pattern'],seed=t['seed'],status='unavailable') for t in manifest['tasks']]}))
            result=output/'analysis.json';analysis.analyze([str(cp)],str(result))
            data=h.read_json(result);self.assertEqual(data['status'],'incomplete')
            self.assertIn('Ca:incomplete_predefined_seeds',data['cases'][0]['states'][0]['failures'])
            self.assertEqual(data['cases'][0]['fixed_count_electronic_contrasts'],[])

    def test_water_rotation_projection_matches_real_gradient_chain_rule(self):
        import hydration_mace as m
        from scipy.spatial.transform import Rotation
        mp=ROOT/'workspaces/hydration_network_20260918/mace_proposal_v2/manifest.json'
        if not mp.exists():self.skipTest('actual DFT checkpoint fixtures unavailable')
        task=h.read_json(mp)['tasks'][0];d=h.read_json(h.verify(task['dft']))
        coords=np.array(d['coordinates_A']);gradient=np.array(d['gradient_Ha_per_bohr'])
        groups=[w for w in task['water_groups'] if w['role']=='variable']
        projected=m.rotational_gradient(coords,gradient,groups)
        direction=np.array([.2,-.4,.3]);eps=1e-5
        w=groups[0];ids=w['hydrogen_indices'];o=coords[w['oxygen_index']]
        plus=(coords[ids]-o)@Rotation.from_rotvec(eps*direction).as_matrix().T+o
        minus=(coords[ids]-o)@Rotation.from_rotvec(-eps*direction).as_matrix().T+o
        numerical=float(np.sum(gradient[ids]*(plus-minus)/(2*eps)))
        self.assertAlmostEqual(float(projected[0]@direction),numerical,places=10)
        # Translation of actual coordinates cannot change rotations about water O.
        np.testing.assert_allclose(m.rotational_gradient(coords+[8.,-3.,10.],gradient,groups),projected,atol=1e-14)

    def test_exact_water_rotation_jacobian_on_real_gradient(self):
        from hydration_proposal_opt import rotate_waters
        from hydration_mace import rotational_gradient
        mp=ROOT/'workspaces/hydration_network_20260918/mace_proposal_v2/manifest.json'
        if not mp.exists():self.skipTest('actual DFT checkpoint fixtures unavailable')
        task=h.read_json(mp)['tasks'][0];d=h.read_json(h.verify(task['dft']))
        initial=np.array(d['coordinates_A']);gradient=np.array(d['gradient_Ha_per_bohr'])
        groups=[w for w in task['water_groups'] if w['role']=='variable']
        for scale in (0.,1e-7,1.,3.):
            parameters=np.tile([.2,-.4,.3],len(groups))*scale
            coords,jacobians=rotate_waters(initial,groups,parameters)
            torque=rotational_gradient(coords,gradient,groups)
            analytic=np.einsum('wij,wi->wj',jacobians,torque).reshape(-1)
            numerical=[]
            for k in range(len(parameters)):
                offset=np.zeros_like(parameters);offset[k]=1e-5
                plus,_=rotate_waters(initial,groups,parameters+offset);minus,_=rotate_waters(initial,groups,parameters-offset)
                numerical.append(np.sum(gradient*(plus-minus))/(2e-5))
            np.testing.assert_allclose(analytic,numerical,atol=1e-10,rtol=1e-7)
            mobile={i for w in groups for i in w['hydrogen_indices']}
            fixed=[i for i in range(len(initial)) if i not in mobile]
            np.testing.assert_array_equal(coords[fixed],initial[fixed])

    def test_real_proposal_transfer_preserves_original_core_and_water_count(self):
        mp=ROOT/'workspaces/hydration_network_20260918/core_transfer_v1/manifest.json'
        if not mp.exists():self.skipTest('real proposed water transfer unavailable')
        m=h.read_json(mp)
        for t in m['tasks']:
            before=h.xyz(h.verify(t['original_xyz']));after=h.xyz(h.verify(t['xyz']))
            self.assertEqual(len(before),40 if t['case']=='1F6S' else 43)
            self.assertEqual([a[0] for a in before],[a[0] for a in after])
            for i in range(len(before)):
                if i not in t['water_H_indices']:self.assertEqual(before[i],after[i])
            parent=h.read_json(h.verify(t['parent']))
            for w in h.square.water_groups(parent,before):
                oh=[np.linalg.norm(np.array(after[i][1:])-after[w['oxygen_index']][1:]) for i in w['hydrogen_indices']]
                oxygen=next(np.array(a[1:]) for a in self.reference if a[0]=='O')
                reference_oh=[np.linalg.norm(np.array(a[1:])-oxygen) for a in self.reference if a[0]=='H']
                np.testing.assert_allclose(oh,reference_oh,atol=2e-9,rtol=0)
            original=h.verify(parent['outputs'][t['metal']]['input']).read_text().splitlines()[0]
            self.assertEqual(h.verify(t['input']).read_text().splitlines()[0],original)

    def test_real_analytic_output_component_parser_ignores_gradient_progress_ellipsis(self):
        mp=ROOT/'workspaces/hydration_network_20260918/proposal_dft_v1/manifest.json'
        if not mp.exists():self.skipTest('real DFT adjudication unavailable')
        for t in h.read_json(mp)['tasks']:
            p=Path(t['output_path']);receipt=Path(str(p)+'.execution.json')
            if not receipt.exists():self.skipTest('actual DFT endpoint not complete')
            e=h.square.endpoint(h.record(p),h.record(receipt),t['xyz'],t['input'])
            self.assertIsNotNone(e['components_hartree']['gCP'])
            self.assertAlmostEqual(sum(e['components_hartree'][k] for k in ('SCF','gCP','dispersion')),e['energy_hartree'],places=10)


if __name__=='__main__':unittest.main()
