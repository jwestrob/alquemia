"""Original-autograd comparisons on real computed densities and source geometry."""
import copy
import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import read_json, verify, xyz

COLLECTION = ROOT/'workspaces/mace_hybrid_20260916/pilot_v3/collection_job_1200308.json'


@unittest.skipUnless(COLLECTION.exists() and importlib.util.find_spec('mace'),
                     'requires real completed MACE cores and isolated installation')
class RealKernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import torch
        from mace.calculators import mace_polar
        from mace_realspace_compat import configure
        cls.torch = torch
        torch.set_num_threads(1)
        cls.collection = read_json(COLLECTION)
        cls.manifest = read_json(verify(cls.collection['manifest']))
        cls.calc = mace_polar(model=str(verify(cls.manifest['model']['checkpoint'])), device='cpu', default_dtype='float64')
        configure(cls.calc)
        cls.reference_features = copy.deepcopy(cls.calc.models[0].electric_potential_descriptor.realspace_features)
        cls.reference_energy = copy.deepcopy(cls.calc.models[0].coulomb_energy.realspace_energy)
        cls.fixtures = []
        for t in cls.manifest['tasks']:
            if t['kind'] == 'core':
                r = cls.collection['rows'][t['task_id']]
                cls.fixtures.append((t['task_id'], np.array([a[1:] for a in xyz(verify(t['xyz']))]),
                                     np.load(verify(r['density_coefficients']))))

    def evaluate(self, module, coords, density, feature, batch=None):
        torch = self.torch
        r = torch.tensor(coords, dtype=torch.float64, requires_grad=True)
        q = torch.tensor(density, dtype=torch.float64, requires_grad=True)
        b = torch.zeros(len(coords), dtype=torch.long) if batch is None else batch
        if feature:
            values = module(source_feats=q, node_positions=r, batch=b)[0]
        else:
            values = module(source_feats=q, positions=r, batch=b)
        gradients = torch.autograd.grad(values.square().sum(), (q,r))
        return [values.detach(), *[g.detach() for g in gradients]]

    def test_real_features_energies_charge_and_position_derivatives(self):
        from mace_blocked import configure
        for tile in (17, 64):
            # Restore original modules before configuring another execution tile.
            self.calc.models[0].electric_potential_descriptor.realspace_features = copy.deepcopy(self.reference_features)
            self.calc.models[0].coulomb_energy.realspace_energy = copy.deepcopy(self.reference_energy)
            if hasattr(self.calc.models[0], '_alquemia_blocked_kernel'):
                del self.calc.models[0]._alquemia_blocked_kernel
            configure(self.calc, tile)
            for name, coords, density in self.fixtures:
                for feature in (True, False):
                    with self.subTest(tile=tile, task=name, feature=feature):
                        old = self.reference_features if feature else self.reference_energy
                        new = (self.calc.models[0].electric_potential_descriptor.realspace_features if feature
                               else self.calc.models[0].coulomb_energy.realspace_energy)
                        expected = self.evaluate(old, coords, density, feature)
                        actual = self.evaluate(new, coords, density, feature)
                        for a,b in zip(actual, expected):
                            self.torch.testing.assert_close(a,b,atol=1e-8,rtol=1e-10)

    def test_batch_excludes_interactions_between_real_core_systems(self):
        # Actual original fixtures, grouped as separate systems, with no moved atoms.
        torch = self.torch
        coords = np.concatenate([f[1] for f in self.fixtures])
        density = np.concatenate([f[2] for f in self.fixtures])
        batch = torch.cat([torch.full((len(f[1]),),i,dtype=torch.long) for i,f in enumerate(self.fixtures)])
        from mace_blocked import configure
        self.calc.models[0].electric_potential_descriptor.realspace_features = copy.deepcopy(self.reference_features)
        self.calc.models[0].coulomb_energy.realspace_energy = copy.deepcopy(self.reference_energy)
        if hasattr(self.calc.models[0], '_alquemia_blocked_kernel'):
            del self.calc.models[0]._alquemia_blocked_kernel
        configure(self.calc, 64)
        new = self.calc.models[0].coulomb_energy.realspace_energy
        for a,b in zip(self.evaluate(new,coords,density,False,batch),
                       self.evaluate(self.reference_energy,coords,density,False,batch)):
            torch.testing.assert_close(a,b,atol=1e-8,rtol=1e-10)

    def test_saved_pair_storage_is_linear_on_real_input(self):
        from mace_blocked import potential
        torch = self.torch
        _, coords, density = self.fixtures[-1]
        r = torch.tensor(coords, dtype=torch.float64, requires_grad=True)
        q = torch.tensor(density[:,0], dtype=torch.float64, requires_grad=True)
        b = torch.zeros(len(coords), dtype=torch.long)
        widths = self.reference_features.total_width_factors
        saved = []
        def pack(t):
            saved.append(tuple(t.shape)); return t
        with torch.autograd.graph.saved_tensors_hooks(pack, lambda t:t):
            output = potential(q,r,b,widths,1,17)
        self.assertEqual(saved, [(len(q),), (len(q),3), (len(q),), tuple(widths.shape)])
        torch.autograd.grad(output.sum(),(q,r))

    def test_local_chunking_preserves_checkpoint_state(self):
        from mace_local_memory import configure
        state = {k:v.clone() for k,v in self.calc.models[0].state_dict().items()}
        result = configure(self.calc,128,17)
        self.assertEqual(len(result['modules']),2)
        self.assertEqual(len(result['product_modules']),2)
        after = self.calc.models[0].state_dict()
        self.assertEqual(state.keys(),after.keys())
        self.assertTrue(all(self.torch.equal(state[k],after[k]) for k in state))

    def test_edge_accumulation_derivatives_without_saved_messages(self):
        from mace_local_memory import AccumulateEdges
        torch = self.torch
        _, coords, density = self.fixtures[0]
        # Real core neighbor graph supplies duplicate receiver indices.
        from mace.data.neighborhood import get_neighborhood
        edges = get_neighborhood(coords,6.)[0]
        sender,receiver = [torch.tensor(x,dtype=torch.long) for x in edges]
        q = torch.tensor(density,dtype=torch.float64,requires_grad=True)
        base = torch.zeros_like(q,requires_grad=True)
        expected = base.index_add(0,receiver,q[sender])
        saved = []
        def pack(t): saved.append(tuple(t.shape)); return t
        values = q[sender]
        with torch.autograd.graph.saved_tensors_hooks(pack,lambda t:t):
            actual = AccumulateEdges.apply(base,receiver,values)
        self.assertEqual(saved,[tuple(receiver.shape)])
        torch.testing.assert_close(actual,expected,rtol=0,atol=0)
        for a,b in zip(torch.autograd.grad(actual.square().sum(),(base,q),retain_graph=True),
                       torch.autograd.grad(expected.square().sum(),(base,q))):
            torch.testing.assert_close(a,b,rtol=0,atol=0)


if __name__ == '__main__':
    unittest.main()
