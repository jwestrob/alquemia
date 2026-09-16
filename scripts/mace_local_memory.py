"""Checkpointed edge execution of the released nonlinear MACE interaction.

Same operations/weights as mace.modules.blocks.RealAgnosticResidualNonLinear-
InteractionBlock (MACE0.3.16, MIT license); only edge batching and activation
retention change. All neighbors contribute before density normalization and
node nonlinearities. No neighborhood truncation or independent cluster scoring.
"""
from types import MethodType
import torch
from torch.utils.checkpoint import checkpoint, set_checkpoint_early_stop

LOCAL_MEMORY_ID = 'mace0316_edge_product_checkpoint_v5'


class AccumulateEdges(torch.autograd.Function):
    """index_add whose backward needs indices, never the original messages.

    torch.index_add retains its source tensor even though its derivative is just
    gathering the outgoing adjoint. That defeats checkpointing large messages.
    """
    @staticmethod
    def forward(ctx, accumulated, indices, values):
        ctx.save_for_backward(indices)
        result = accumulated.clone()
        result.index_add_(0,indices,values)
        return result

    @staticmethod
    def backward(ctx, gradient):
        indices, = ctx.saved_tensors
        return gradient, None, gradient.index_select(0,indices)


def configure(calc, edge_tile, node_tile=None):
    if not isinstance(edge_tile, int) or edge_tile < 1:
        raise ValueError('edge tile must be a positive integer')
    patched = []
    for name, module in calc.models[0].named_modules():
        if type(module).__name__ != 'RealAgnosticResidualNonLinearInteractionBlock':
            continue
        if hasattr(module, 'conv_fusion') or hasattr(module, '_alquemia_edge_tile'):
            raise ValueError('fused/already-patched local interaction unsupported')

        def forward(self, node_attrs, node_feats, edge_attrs, edge_feats, edge_index,
                    cutoff=None, lammps_class=None, lammps_natoms=(0,0), first_layer=False):
            if lammps_class is not None:
                raise ValueError('blocked local inference does not support LAMMPS')
            num_nodes = len(node_feats)
            sc = self.skip_tp(node_feats)
            up = self.linear_up(node_feats)
            residual = self.linear_res(up)
            source = self.source_embedding(node_attrs)
            target = self.target_embedding(node_attrs)

            def edge_messages(up, source, target, attrs, radial, indices, envelope):
                sender, receiver = indices
                augmented = torch.cat((radial,source[sender],target[receiver]),dim=-1)
                weights = self.conv_tp_weights(augmented)
                density = torch.tanh(self.density_fn(augmented)**2)
                if envelope is not None:
                    weights = weights * envelope
                    density = density * envelope
                return self.conv_tp(up[sender],attrs,weights), density

            message = None
            density = up.new_zeros((num_nodes,1))
            for start in range(0, edge_index.shape[1], self._alquemia_edge_tile):
                stop = start+self._alquemia_edge_tile
                indices = edge_index[:,start:stop]
                # TorchScript wraps the checkpoint early-stop exception as an
                # opaque RuntimeError. Recompute the complete block instead.
                with set_checkpoint_early_stop(False):
                    values, local_density = checkpoint(
                        edge_messages,up,source,target,edge_attrs[start:stop],edge_feats[start:stop],
                        indices,None if cutoff is None else cutoff[start:stop],
                        use_reentrant=False, preserve_rng_state=False)
                if message is None:
                    message = values.new_zeros((num_nodes,values.shape[-1]))
                message = AccumulateEdges.apply(message,indices[1],values)
                density = AccumulateEdges.apply(density,indices[1],local_density)
            if message is None:
                raise ValueError('no edges in selected real protein system')
            message = self.linear_1(message) / (density*self.beta+self.alpha)
            message = self.equivariant_nonlin(message+residual)
            return self.reshape(self.linear_2(message)), sc

        module._alquemia_edge_tile = edge_tile
        module.forward = MethodType(forward,module)
        patched.append(name)
    if not patched or len(patched) != len(calc.models[0].interactions):
        raise ValueError('unsupported checkpoint interaction architecture')
    products, sparse_products = [], []
    if node_tile is not None:
        if not isinstance(node_tile,int) or node_tile < 1:
            raise ValueError('positive node tile required')
        def wrap_product(original):
            def forward(self, node_feats, sc, node_attrs):
                pieces = []
                for start in range(0,len(node_feats),node_tile):
                    stop = start+node_tile
                    with set_checkpoint_early_stop(False):
                        pieces.append(checkpoint(original,node_feats[start:stop],
                            None if sc is None else sc[start:stop],node_attrs[start:stop],
                            use_reentrant=False,preserve_rng_state=False))
                return torch.cat(pieces,dim=0)
            return forward
        for name,module in calc.models[0].named_modules():
            if type(module).__name__ == 'EquivariantProductBasisBlock':
                module.forward = MethodType(wrap_product(module.forward),module)
                products.append(name)
        if len(products) != len(calc.models[0].products):
            raise ValueError('unsupported checkpoint product architecture')
        def wrap_sparse(original):
            def forward(self,x1,x2):
                pieces = []
                for start in range(0,len(x1),node_tile):
                    stop = start+node_tile
                    with set_checkpoint_early_stop(False):
                        pieces.append(checkpoint(original,x1[start:stop],x2[start:stop],
                            use_reentrant=False,preserve_rng_state=False))
                return torch.cat(pieces,dim=0)
            return forward
        for name,module in calc.models[0].named_modules():
            if type(module).__name__ == 'SparseUvuTensorProduct':
                module.forward = MethodType(wrap_sparse(module.forward),module)
                sparse_products.append(name)
    return {'implementation_id': LOCAL_MEMORY_ID, 'edge_tile': edge_tile,
            'modules': patched, 'node_tile': node_tile, 'product_modules': products,
            'sparse_product_modules': sparse_products,
            'checkpoint_backward': True, 'neighbor_cutoff_changed': False}
