"""Exact edge batching for the pinned OMOL nonlinear interaction; energy only."""
from __future__ import annotations
import types
from affordable_common import InvalidArtifact

ADAPTER='omol_nonlinear_exact_edge_batches_v1'
BLOCK='RealAgnosticResidualNonLinearInteractionBlock'


def forward(self,node_attrs,node_feats,edge_attrs,edge_feats,edge_index,cutoff=None,
            lammps_class=None,lammps_natoms=(0,0),first_layer=False):
    import torch
    if torch.is_grad_enabled() or lammps_class is not None or tuple(lammps_natoms)!=(0,0):
        raise InvalidArtifact('exact edge adapter supports no-grad, non-LAMMPS inference only')
    n=node_feats.shape[0];count=edge_index.shape[1];chunk=self._omol_edge_chunk_size
    sc=self.skip_tp(node_feats)
    feats=self.linear_up(node_feats)
    residual=self.linear_res(feats)
    source=self.source_embedding(node_attrs);target=self.target_embedding(node_attrs)
    message=feats.new_zeros((n,self.conv_tp.irreps_out.dim));density=feats.new_zeros((n,1))
    processed=0
    for start in range(0,count,chunk):
        stop=min(start+chunk,count);sender=edge_index[0,start:stop];receiver=edge_index[1,start:stop]
        ef=torch.cat((edge_feats[start:stop],source[sender],target[receiver]),dim=-1)
        weights=self.conv_tp_weights(ef)
        ed=torch.tanh(self.density_fn(ef)**2)
        if cutoff is not None:
            weights=weights*cutoff[start:stop];ed=ed*cutoff[start:stop]
        density.index_add_(0,receiver,ed)
        message.index_add_(0,receiver,self.conv_tp(feats[sender],edge_attrs[start:stop],weights))
        processed+=stop-start
    if processed!=count:raise InvalidArtifact('edge inventory was not fully processed')
    self._omol_edge_calls.append({'atoms':n,'edges':count,'processed_edges':processed,'chunk_size':chunk})
    message=self.linear_1(message)/(density*self.beta+self.alpha)
    message=message+residual
    message=self.equivariant_nonlin(message)
    message=self.linear_2(message)
    return self.reshape(message),sc


def install(model,chunk_size):
    if chunk_size not in (1024,2048):raise InvalidArtifact('undeclared edge chunk size')
    if type(model).__name__!='ScaleShiftMACE' or len(model.interactions)!=3:
        raise InvalidArtifact('unsupported native OMOL architecture')
    for block in model.interactions:
        if (type(block).__name__!=BLOCK or hasattr(block,'conv_fusion')
                or getattr(block,'cueq_config',None) or getattr(block,'oeq_config',None)):
            raise InvalidArtifact('unsupported interaction/backend for exact edge batching')
        if hasattr(block,'_omol_edge_chunk_size'):raise InvalidArtifact('adapter already installed')
    for block in model.interactions:
        block._omol_edge_chunk_size=chunk_size;block._omol_edge_calls=[]
        block.forward=types.MethodType(forward,block)


def receipt(model,chunk_size):
    calls=[b._omol_edge_calls for b in model.interactions]
    if any(len(c)!=1 for c in calls):raise InvalidArtifact('each interaction must execute exactly once')
    layers=[c[0] for c in calls]
    if len({(r['atoms'],r['edges'],r['processed_edges'],r['chunk_size']) for r in layers})!=1:
        raise InvalidArtifact('edge adapter changed graph across layers')
    return {'id':ADAPTER,'chunk_size':chunk_size,'layers':layers,'gradient_support':False}
