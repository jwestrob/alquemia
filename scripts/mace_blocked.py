"""Bounded-memory equivalent of graph-longrange 0.4.4 realspace pair sums.

Retains the released displaced-multipole representation and all exclusions.
First derivatives w.r.t. positions AND charges are evaluated analytically in
blocks, including their propagation through the learned charge-response layers.
No global N*N arrays or pair activations are retained. Inference only.
"""
import math
from types import MethodType

import torch
from torch.autograd.function import once_differentiable
from graph_longrange.utils import FIELD_CONSTANT

KERNEL_ID = 'graph044_blocked_analytic_backward_v1'
COULOMB_FACTOR = FIELD_CONSTANT / (4 * math.pi)
EPSILON = 1e-6  # Exact released kernel denominator regularization, in Angstrom.


def pair_terms(positions, batch, multiplicity, widths, i, j, tile, derivative=False):
    """Rows are receiver sites; columns are source sites (including dipole sites)."""
    ri, rj = positions[i:i+tile], positions[j:j+tile]
    difference = ri[:, None, :] - rj[None, :, :]
    distance = torch.linalg.vector_norm(difference, dim=-1)
    indices_i = torch.arange(i, i+len(ri), device=positions.device) // multiplicity
    indices_j = torch.arange(j, j+len(rj), device=positions.device) // multiplicity
    allowed = ((batch[i:i+tile, None] == batch[None, j:j+tile]) &
               (indices_i[:, None] != indices_j[None, :]))
    x = distance[..., None] * 0.5 / widths
    erf = torch.erf(x)
    denominator = distance[..., None] + EPSILON
    kernel = torch.where(allowed[..., None], erf / denominator, 0.)
    if not derivative:
        return kernel
    radial = ((torch.exp(-x*x) / (math.sqrt(math.pi)*widths))*denominator - erf) / denominator.square()
    radial = torch.where(allowed[..., None], radial, 0.)
    direction = difference / distance.clamp_min(torch.finfo(distance.dtype).tiny)[..., None]
    return kernel, radial, direction


class PairPotential(torch.autograd.Function):
    @staticmethod
    def forward(ctx, charges, positions, batch, widths, multiplicity, tile):
        if charges.dtype != torch.float64 or positions.dtype != torch.float64:
            raise ValueError('validated blocked kernel requires float64')
        if widths.requires_grad or multiplicity not in (1, 4) or tile < 1:
            raise ValueError('unsupported width gradients, multiplicity or tile size')
        if positions.shape != (len(charges), 3) or batch.shape != charges.shape:
            raise ValueError('incompatible pair-kernel tensor shapes')
        ctx.save_for_backward(charges, positions, batch, widths)
        ctx.multiplicity, ctx.tile = multiplicity, tile
        result = charges.new_zeros((len(charges), len(widths)))
        for i in range(0, len(charges), tile):
            for j in range(0, len(charges), tile):
                kernel = pair_terms(positions, batch, multiplicity, widths, i, j, tile)
                result[i:i+tile] += torch.einsum('ijw,j->iw', kernel, charges[j:j+tile])
        return result * COULOMB_FACTOR

    @staticmethod
    @once_differentiable
    def backward(ctx, gradient):
        charges, positions, batch, widths = ctx.saved_tensors
        multiplicity, tile = ctx.multiplicity, ctx.tile
        dq, dr = torch.zeros_like(charges), torch.zeros_like(positions)
        for i in range(0, len(charges), tile):
            gi = gradient[i:i+tile] * COULOMB_FACTOR
            for j in range(0, len(charges), tile):
                kernel, radial, direction = pair_terms(positions, batch, multiplicity, widths, i, j, tile, True)
                dq[j:j+tile] += torch.einsum('iw,ijw->j', gi, kernel)
                coefficient = torch.einsum('iw,ijw->ij', gi, radial) * charges[None, j:j+tile]
                pair_gradient = coefficient[..., None] * direction
                dr[i:i+tile] += pair_gradient.sum(dim=1)
                dr[j:j+tile] -= pair_gradient.sum(dim=0)
        return dq, dr, None, None, None, None


def potential(charges, positions, batch, widths, multiplicity, tile):
    return PairPotential.apply(charges, positions, batch, widths, multiplicity, tile)


def displaced(module, source_feats, positions, batch):
    """Same offsets, charge signs and spherical component order as upstream."""
    extended_positions = positions.repeat_interleave(4, dim=0)
    extended_positions[1::4] += module.x
    extended_positions[2::4] += module.y
    extended_positions[3::4] += module.z
    charges = torch.zeros_like(extended_positions[:, 0])
    charges[1::4] = source_feats[:, 3] / module.offset
    charges[2::4] = source_feats[:, 1] / module.offset
    charges[3::4] = source_feats[:, 2] / module.offset
    charges[0::4] = source_feats[:, 0] - (charges[1::4]+charges[2::4]+charges[3::4])
    return charges, extended_positions, batch.repeat_interleave(4)


def configure(calc, tile=512, edge_tile=None, node_tile=None):
    """Replace pair execution only; leave feature assembly and self terms intact."""
    if not isinstance(tile, int) or tile < 1:
        raise ValueError('positive integer tile required')
    model = calc.models[0]
    f = model.electric_potential_descriptor.realspace_features
    e = model.coulomb_energy.realspace_energy
    if f.density_max_l > 1 or f.projection_max_l > 1 or e.density_max_l > 1:
        raise ValueError('only released monopole/dipole realspace representation supported')
    if hasattr(model, '_alquemia_blocked_kernel'):
        raise ValueError('blocked kernel already configured')

    def features_l0(self, source_feats, positions, batch):
        values = potential(source_feats[:, 0], positions, batch, self.total_width_factors, 1, tile)
        return self.l0_factors * values

    def features_l1(self, source_feats, positions, batch):
        q, r, b = displaced(self, source_feats, positions, batch)
        scalar = potential(q, r, b, self.total_width_factors, 4, tile)
        result = source_feats.new_zeros((len(batch), 4*self.num_radial))
        result[:, :self.num_radial] = self.l0_factors * scalar[0::4]
        result[:, self.num_radial::3] = self.l1_factors * (scalar[2::4]-scalar[0::4])
        result[:, self.num_radial+1::3] = self.l1_factors * (scalar[3::4]-scalar[0::4])
        result[:, self.num_radial+2::3] = self.l1_factors * (scalar[1::4]-scalar[0::4])
        return result

    def energy(self, source_feats, positions, batch):
        if self.density_max_l == 0:
            q, r, b, multiplicity = source_feats[:, 0], positions, batch, 1
        else:
            q, r, b = displaced(self, source_feats, positions, batch)
            multiplicity = 4
        widths = q.new_tensor([self.density_smearing_width])
        values = potential(q, r, b, widths, multiplicity, tile)[:, 0]
        energies = q.new_zeros((int(batch.max())+1,))
        energies.index_add_(0, b, 0.5*q*values)
        if self.include_self_interaction:
            self_fields = self.self_interaction(source_feats)
            energies.index_add_(0, batch, 0.5*torch.einsum('nb,nb->n', source_feats, self_fields))
        return energies

    f.call_density_0_feats_0 = MethodType(features_l0, f)
    f.call_density_1_feats_1 = MethodType(features_l1, f)
    e.forward = MethodType(energy, e)
    original_forward = model.forward

    def inference(self, data, *args, **kwargs):
        for key in ('training', 'compute_hessian', 'compute_edge_forces', 'compute_atomic_stresses'):
            if kwargs.get(key, False):
                raise ValueError(f'blocked first-derivative inference does not support {key}')
        return original_forward(data, *args, **kwargs)

    model.forward = MethodType(inference, model)
    model._alquemia_blocked_kernel = KERNEL_ID
    local = None
    if edge_tile is not None:
        from mace_local_memory import configure as configure_local
        local = configure_local(calc, edge_tile, node_tile)
    return {'kernel_id': KERNEL_ID, 'tile_sites': tile, 'local_memory': local, 'pair_work_order': 'quadratic',
            'saved_pair_activations': False, 'precision': 'float64', 'analytic_derivative_order': 1,
            'unchanged_offsets_A': {'features': f.offset, 'energy': e.offset}, 'epsilon_A': EPSILON}
