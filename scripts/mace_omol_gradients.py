"""Exact, separately qualified OMOL autograd batching with activation recomputation."""
from __future__ import annotations
from types import MethodType
from affordable_common import InvalidArtifact

ADAPTER = 'omol_exact_checkpointed_edge_product_autograd_v3'
CHUNK = 1024


def install_mask(model):
    """The categorical charge feature is a coordinate-independent constant zero."""
    import torch
    from mace_omol_ablation import install
    install(model)
    emb = model.joint_embedding.embedders['total_charge']
    original = emb.forward

    def forward(this, indices):
        # Only the original constant-zero categorical feature is outside autograd.
        # All geometry-dependent operations retain their full derivative graph.
        with torch.no_grad():
            return original(indices)
    emb.forward = MethodType(forward, emb)


def interaction(self, node_attrs, node_feats, edge_attrs, edge_feats, edge_index,
                cutoff=None, lammps_class=None, lammps_natoms=(0,0), first_layer=False):
    import torch
    from torch.utils.checkpoint import checkpoint
    if lammps_class is not None or tuple(lammps_natoms)!=(0,0):
        raise InvalidArtifact('checkpoint adapter requires a non-LAMMPS graph')
    n = node_feats.shape[0]; count = edge_index.shape[1]
    sc = self.skip_tp(node_feats); feats = self.linear_up(node_feats)
    residual = self.linear_res(feats)
    source = self.source_embedding(node_attrs); target = self.target_embedding(node_attrs)
    message = feats.new_zeros((n,self.conv_tp.irreps_out.dim)); density = feats.new_zeros((n,1))

    def kernel(full_feats, full_source, full_target, attrs, radial, indices, envelope, start, stop):
        self._alquemia_kernel_calls += 1
        sender, receiver = indices[0,start:stop], indices[1,start:stop]
        ef = torch.cat((radial[start:stop],full_source[sender],full_target[receiver]),dim=-1)
        weights = self.conv_tp_weights(ef); ed = torch.tanh(self.density_fn(ef)**2)
        if envelope is not None:
            weights = weights*envelope[start:stop]; ed = ed*envelope[start:stop]
        return self.conv_tp(full_feats[sender],attrs[start:stop],weights), ed

    for start in range(0,count,CHUNK):
        stop = min(start+CHUNK,count); receiver = edge_index[1,start:stop]
        args = (feats,source,target,edge_attrs,edge_feats,edge_index,cutoff,start,stop)
        messages, ed = (checkpoint(kernel,*args,use_reentrant=False,preserve_rng_state=False)
                        if torch.is_grad_enabled() else kernel(*args))
        # IndexAddBackward retains each source message tensor even though its
        # derivative only gathers the destination gradient. ScatterAddBackward
        # retains the expanded index view instead. Both sum the same messages
        # onto the same receivers, in the same edge order; no graph is detached.
        density.scatter_add_(0,receiver[:,None].expand_as(ed),ed)
        message.scatter_add_(0,receiver[:,None].expand_as(messages),messages)
    self._alquemia_forward_calls.append({'atoms':n,'edges':count,'chunks':(count+CHUNK-1)//CHUNK})
    message = self.linear_1(message)/(density*self.beta+self.alpha)
    message = self.equivariant_nonlin(message+residual)
    return self.reshape(self.linear_2(message)), sc


def product(self, node_feats, sc, node_attrs):
    import torch
    from torch.utils.checkpoint import checkpoint
    count = node_feats.shape[0]
    if count != node_attrs.shape[0] or (sc is not None and count != sc.shape[0]):
        raise InvalidArtifact('product atom inventory differs')

    def kernel(full_feats, full_sc, full_attrs, start, stop):
        self._alquemia_kernel_calls += 1
        return self._alquemia_original_forward(node_feats=full_feats[start:stop],
            sc=None if full_sc is None else full_sc[start:stop],node_attrs=full_attrs[start:stop])
    outputs = []
    for start in range(0,count,CHUNK):
        stop = min(start+CHUNK,count); args = (node_feats,sc,node_attrs,start,stop)
        outputs.append(checkpoint(kernel,*args,use_reentrant=False,preserve_rng_state=False)
                       if torch.is_grad_enabled() else kernel(*args))
    self._alquemia_forward_calls.append({'atoms':count,'chunks':(count+CHUNK-1)//CHUNK})
    return torch.cat(outputs,dim=0)


def install(model):
    if type(model).__name__!='ScaleShiftMACE' or len(model.interactions)!=3 or len(model.products)!=3:
        raise InvalidArtifact('unsupported checkpointed OMOL architecture')
    for block in model.interactions:
        if (type(block).__name__!='RealAgnosticResidualNonLinearInteractionBlock'
                or hasattr(block,'conv_fusion') or getattr(block,'cueq_config',None)
                or getattr(block,'oeq_config',None) or hasattr(block,'_alquemia_forward_calls')):
            raise InvalidArtifact('unsupported/already adapted interaction')
    for block in model.products:
        if (type(block).__name__!='EquivariantProductBasisBlock' or getattr(block,'cueq_config',None)
                or getattr(block,'oeq_config',None) or hasattr(block,'_alquemia_forward_calls')):
            raise InvalidArtifact('unsupported/already adapted product')
    for block in model.interactions:
        block._alquemia_forward_calls=[]; block._alquemia_kernel_calls=0
        block.forward=MethodType(interaction,block)
    for block in model.products:
        block._alquemia_forward_calls=[]; block._alquemia_kernel_calls=0
        block._alquemia_original_forward=block.forward
        block.forward=MethodType(product,block)


def receipt(model, gradients):
    groups={}
    for name, blocks in [('interactions',model.interactions),('products',model.products)]:
        rows=[]
        for b in blocks:
            if len(b._alquemia_forward_calls)!=1:
                raise InvalidArtifact('one full model forward required')
            row=b._alquemia_forward_calls[0]
            if b._alquemia_kernel_calls < row['chunks']:
                raise InvalidArtifact('incomplete checkpoint execution')
            rows.append({**row,'kernel_invocations_including_recomputations':b._alquemia_kernel_calls})
        groups[name]=rows
    return {'id':ADAPTER,'chunk_size':CHUNK,'use_reentrant':False,'preserve_rng_state':False,'early_stop':False,
            'aggregation':'scatter_add_expanded_receiver_index',
            'analytic_gradients_requested':gradients,'coordinate_graph_preserved':True,**groups}
