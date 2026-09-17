"""Blocked analytic monopole/dipole fields of the pinned regularized kernel.

Opt-in replacement for fixed-axis finite displacements. Scalar/vector fields
are derivatives of erf(r/(2w))/(r+epsilon); widths, normalization and self terms
are retained. Backward recomputes one tile and differentiates its analytic
expression, including charge and dipole response. No numerical derivatives.
"""
import math
from types import MethodType
import torch
from torch.autograd.function import once_differentiable
from graph_longrange.utils import FIELD_CONSTANT

KERNEL_ID='graph044_analytic_multipoles_blocked_v1'
EPSILON=1e-6
COULOMB_FACTOR=FIELD_CONSTANT/(4*math.pi)


def block_fields(ri,rj,source,widths,allowed):
    """Output [receivers,widths,(potential,Cartesian gradient)], no self pairs."""
    d=ri[:,None,:]-rj[None,:,:]
    r=torch.sqrt(torch.where(allowed,(d*d).sum(-1),1.))
    rw=r[...,None]; den=rw+EPSILON
    t=rw/(2*widths); erf=torch.erf(t)
    first_erf=torch.exp(-t*t)/(math.sqrt(math.pi)*widths)
    f=erf/den
    first=first_erf/den-erf/den.square()
    second=-(rw/(2*widths.square()))*first_erf/den-2*first_erf/den.square()+2*erf/den.pow(3)
    a=first/rw; b=(second-a)/rw.square()
    f=torch.where(allowed[...,None],f,0.)
    a=torch.where(allowed[...,None],a,0.)
    b=torch.where(allowed[...,None],b,0.)
    q,p=source[:,0],source[:,1:]
    pd=(d*p[None,:,:]).sum(-1)[...,None]
    scalar=(f*q[None,:,None]-a*pd).sum(1)
    vector=((a*q[None,:,None]-b*pd)[...,None]*d[:,:,None,:]
            -a[...,None]*p[None,:,None,:]).sum(1)
    return torch.cat((scalar[...,None],vector),dim=-1)*COULOMB_FACTOR


class MultipoleFields(torch.autograd.Function):
    @staticmethod
    def forward(ctx, source, positions, batch, widths, tile):
        if source.dtype!=torch.float64 or positions.dtype!=torch.float64:
            raise ValueError('analytic multipoles require float64')
        if source.shape!=(len(positions),4) or positions.shape!=(len(source),3) or batch.shape!=(len(source),):
            raise ValueError('invalid multipole input shapes')
        if widths.requires_grad or tile<1 or not bool((widths>0).all()):
            raise ValueError('unsupported widths/tile')
        identity=torch.cat((batch[:,None].to(positions.dtype),positions),dim=1)
        if len(torch.unique(identity,dim=0))!=len(positions):
            raise ValueError('coincident physical atom centers unsupported')
        ctx.save_for_backward(source,positions,batch,widths);ctx.tile=tile
        result=source.new_zeros((len(source),len(widths),4))
        for i in range(0,len(source),tile):
            for j in range(0,len(source),tile):
                ii=torch.arange(i,min(i+tile,len(source)),device=batch.device)
                jj=torch.arange(j,min(j+tile,len(source)),device=batch.device)
                allowed=(batch[ii,None]==batch[None,jj])&(ii[:,None]!=jj[None,:])
                result[i:i+tile]+=block_fields(positions[i:i+tile],positions[j:j+tile],source[j:j+tile],widths,allowed)
        return result

    @staticmethod
    @once_differentiable
    def backward(ctx, gradient):
        source,positions,batch,widths=ctx.saved_tensors;tile=ctx.tile
        ds=torch.zeros_like(source);dr=torch.zeros_like(positions)
        for i in range(0,len(source),tile):
            for j in range(0,len(source),tile):
                ii=torch.arange(i,min(i+tile,len(source)),device=batch.device)
                jj=torch.arange(j,min(j+tile,len(source)),device=batch.device)
                allowed=(batch[ii,None]==batch[None,jj])&(ii[:,None]!=jj[None,:])
                with torch.enable_grad():
                    ri=positions[i:i+tile].detach().requires_grad_(True)
                    rj=positions[j:j+tile].detach().requires_grad_(True)
                    q=source[j:j+tile].detach().requires_grad_(True)
                    fields=block_fields(ri,rj,q,widths,allowed)
                    gi,gj,gq=torch.autograd.grad(fields,(ri,rj,q),gradient[i:i+tile].detach())
                dr[i:i+tile]+=gi;dr[j:j+tile]+=gj;ds[j:j+tile]+=gq
        return ds,dr,None,None,None


def fields(source,positions,batch,widths,tile=256):
    return MultipoleFields.apply(source,positions,batch,widths,tile)


def cartesian(source):
    if source.shape[1]==1:
        return torch.cat((source,source.new_zeros((len(source),3))),dim=1)
    return source[:,[0,3,1,2]]


def configure(calc,tile=256,edge_tile=4096,node_tile=256):
    if not isinstance(tile,int) or tile<1:
        raise ValueError('positive tile required')
    model=calc.models[0]
    f=model.electric_potential_descriptor.realspace_features
    e=model.coulomb_energy.realspace_energy
    if (f.density_max_l>1 or f.projection_max_l>1 or e.density_max_l>1 or
        hasattr(model,'_alquemia_blocked_kernel')):
        raise ValueError('unsupported/already configured multipole architecture')
    def features_l0(self,source_feats,positions,batch):
        v=fields(cartesian(source_feats),positions,batch,self.total_width_factors,tile)
        return self.l0_factors*v[:,:,0]
    def features_l1(self,source_feats,positions,batch):
        v=fields(cartesian(source_feats),positions,batch,self.total_width_factors,tile)
        scalar=self.l0_factors*v[:,:,0]
        # Saved finite-difference factor contains 1/offset; analytic gradient does not.
        vector=(self.l1_factors*self.offset)[None,:,None]*v[:,:,[2,3,1]]
        return torch.cat((scalar,vector.reshape(len(source_feats),-1)),dim=1)
    def energy(self,source_feats,positions,batch):
        q=cartesian(source_feats)
        v=fields(q,positions,batch,q.new_tensor([self.density_smearing_width]),tile)[:,0,:]
        node=0.5*(q*v).sum(-1)
        if self.include_self_interaction:
            node=node+0.5*(source_feats*self.self_interaction(source_feats)).sum(-1)
        result=q.new_zeros((int(batch.max())+1,));result.index_add_(0,batch,node)
        return result
    f.call_density_0_feats_0=MethodType(features_l0,f)
    f.call_density_1_feats_1=MethodType(features_l1,f)
    e.forward=MethodType(energy,e)
    original=model.forward
    def inference(self,data,*args,**kwargs):
        if args:
            raise ValueError('analytic inference requires named forward controls')
        for key in ('training','compute_hessian','compute_edge_forces','compute_atomic_stresses'):
            if kwargs.get(key,False):raise ValueError(f'analytic inference does not support {key}')
        return original(data,**kwargs)
    model.forward=MethodType(inference,model);model._alquemia_blocked_kernel=KERNEL_ID
    from mace_local_memory import configure as local
    memory=local(calc,edge_tile,node_tile)
    return {'kernel_id':KERNEL_ID,'tile_atoms':tile,'local_memory':memory,
            'epsilon_A':EPSILON,'precision':'float64','analytic_derivative_order':1,
            'finite_displacements_used':False,'physical_model_changed':True,
            'weights_changed':False,'saved_pair_activations':False,'pair_work_order':'quadratic'}
