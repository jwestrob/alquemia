"""Analytic fields tested against independent autodiff on real archived cores."""
import copy
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz
from mace_hybrid import rotation
from mace_rotation import rotate_vectors
COLLECTION=ROOT/'workspaces/mace_hybrid_20260916/pilot_v3/collection_job_1200308.json'

@unittest.skipUnless(COLLECTION.exists() and importlib.util.find_spec('mace'),'requires actual completed cores and isolated MACE installation')
class AnalyticMultipoleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import torch
        from mace.calculators import mace_polar
        from mace_realspace_compat import configure
        cls.torch=torch;torch.set_num_threads(1);torch.set_default_dtype(torch.float64)
        cls.c=read_json(COLLECTION);cls.m=read_json(verify(cls.c['manifest']))
        cls.calc=mace_polar(model=str(verify(cls.m['model']['checkpoint'])),device='cpu',default_dtype='float64')
        configure(cls.calc)
        cls.old_state={k:v.clone() for k,v in cls.calc.models[0].state_dict().items()}
        cls.fixtures=[]
        for t in cls.m['tasks']:
            if t['kind']=='core':
                cls.fixtures.append((t['task_id'],np.array([a[1:] for a in xyz(verify(t['xyz']))]),
                    np.load(verify(cls.c['rows'][t['task_id']]['density_coefficients']))))
        from mace_analytic import configure as analytic
        analytic(cls.calc,tile=17,edge_tile=128,node_tile=17)

    def reference(self,q,r,b,widths):
        """Different route: autodiff scalar kernel in Cartesian space, no f'/f'' formulas."""
        torch=self.torch
        from mace_analytic import EPSILON,COULOMB_FACTOR
        pairs=((b[:,None]==b[None,:])&~torch.eye(len(b),dtype=torch.bool)).nonzero()
        i,j=pairs.T
        i=i.repeat_interleave(len(widths));j=j.repeat_interleave(len(widths))
        w=widths.repeat(len(pairs));d=r[i]-r[j]
        def radial(vector,width):
            norm=torch.linalg.vector_norm(vector)
            return torch.erf(norm/(2*width))/(norm+EPSILON)
        values=torch.vmap(radial)(d,w)
        gradfun=torch.func.jacrev(radial,argnums=0)
        gradient=torch.vmap(gradfun)(d,w)
        hessian=torch.vmap(torch.func.jacrev(gradfun,argnums=0))(d,w)
        scalar=q[j,0]*values-(q[j,1:]*gradient).sum(-1)
        vector=q[j,0,None]*gradient-torch.einsum('ijk,ik->ij',hessian,q[j,1:])
        result=q.new_zeros((len(q)*len(widths),4))
        dest=i*len(widths)+torch.arange(len(widths)).repeat(len(pairs))
        result.index_add_(0,dest,torch.cat((scalar[:,None],vector),dim=1))
        return result.reshape(len(q),len(widths),4)*COULOMB_FACTOR

    def test_four_real_cores_fields_and_first_derivatives(self):
        torch=self.torch
        from mace_analytic import fields,cartesian
        widths=self.calc.models[0].electric_potential_descriptor.realspace_features.total_width_factors
        for name,coords,density in self.fixtures:
            for tile in (17,64):
                with self.subTest(core=name,tile=tile):
                    r=torch.tensor(coords,requires_grad=True);q=cartesian(torch.tensor(density)).detach().requires_grad_(True)
                    b=torch.zeros(len(q),dtype=torch.long)
                    expected=self.reference(q,r,b,widths)
                    derivatives=torch.autograd.grad(expected.square().sum(),(q,r))
                    actual=fields(q,r,b,widths,tile)
                    actual_derivatives=torch.autograd.grad(actual.square().sum(),(q,r))
                    for a,e in zip((actual,*actual_derivatives),(expected,*derivatives)):
                        torch.testing.assert_close(a,e,atol=1e-8,rtol=1e-10)

    def test_real_mixed_batch_exclusions(self):
        torch=self.torch
        from mace_analytic import fields,cartesian
        r=torch.tensor(np.concatenate([f[1] for f in self.fixtures]))
        q=cartesian(torch.tensor(np.concatenate([f[2] for f in self.fixtures])))
        b=torch.cat([torch.full((len(f[1]),),i,dtype=torch.long) for i,f in enumerate(self.fixtures)])
        w=self.calc.models[0].electric_potential_descriptor.realspace_features.total_width_factors
        torch.testing.assert_close(fields(q,r,b,w,31),self.reference(q,r,b,w),atol=1e-8,rtol=1e-10)

    def test_energy_features_and_forces_rotate_on_real_densities(self):
        torch=self.torch;R=rotation()
        e=self.calc.models[0].coulomb_energy.realspace_energy
        f=self.calc.models[0].electric_potential_descriptor.realspace_features
        for name,coords,density in self.fixtures:
            b=torch.zeros(len(coords),dtype=torch.long)
            results=[]
            for rr,qq in ((coords,density),(coords@R.T,rotate_vectors(density,R))):
                r=torch.tensor(rr,requires_grad=True);q=torch.tensor(qq)
                energy=e(source_feats=q,positions=r,batch=b)
                force=-torch.autograd.grad(energy.sum(),r)[0]
                features=f(source_feats=q,node_positions=r,batch=b)[0]
                results.append((energy.detach().numpy(),force.detach().numpy(),features.detach().numpy()))
            old,new=results
            for a,z in ((new[0],old[0]),(new[1]@R,old[1]),(new[2],rotate_vectors(old[2],R,f.num_radial))):
                np.testing.assert_allclose(a,z,atol=1e-8,rtol=1e-10)

    def test_saved_storage_linear_and_parameters_unchanged(self):
        torch=self.torch
        from mace_analytic import fields,cartesian
        _,coords,density=self.fixtures[-1]
        q=cartesian(torch.tensor(density)).detach().requires_grad_(True)
        r=torch.tensor(coords,requires_grad=True);b=torch.zeros(len(q),dtype=torch.long)
        w=self.calc.models[0].electric_potential_descriptor.realspace_features.total_width_factors
        saved=[]
        def pack(t):saved.append(t.shape);return t
        with torch.autograd.graph.saved_tensors_hooks(pack,lambda t:t):fields(q,r,b,w,17)
        self.assertEqual(saved,[q.shape,r.shape,b.shape,w.shape])
        current=self.calc.models[0].state_dict()
        self.assertEqual(current.keys(),self.old_state.keys())
        for name,value in current.items():self.assertTrue(torch.equal(value,self.old_state[name]),name)

if __name__=='__main__':unittest.main()
