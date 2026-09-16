"""Isolated PolarMACE 0.3.16 / graph-longrange 0.4.4 interface repair.

No electrostatic kernel or learned parameter is changed. The released calculator
ignores pbc_handling, its forward passes obsolete backend keywords and a singleton
density channel, and it constructs an unused reciprocal grid for isolated atoms.
Only nonperiodic, uncompiled inference is supported by this adapter.
"""
from contextlib import contextmanager
from types import MethodType

ADAPTER_ID = 'polar0316_graph044_isolated_interface_v2'


def require_isolated(pbc, force_pbc_evaluator=False):
    if force_pbc_evaluator or bool(pbc.any()):
        raise ValueError('isolated realspace adapter rejects periodic evaluation')


def configure(calc):
    import torch
    import mace.modules.extensions as extensions

    if calc.use_compile or len(calc.models) != 1:
        raise ValueError('adapter requires one uncompiled PolarMACE model')
    model = calc.models[0]
    if hasattr(model, '_alquemia_realspace_adapter'):
        raise ValueError('adapter already installed')
    features = model.electric_potential_descriptor
    coulomb = model.coulomb_energy
    restored_dimensions = []
    for name, module in model.named_modules():
        if type(module).__name__ == 'GTOSelfInteractionBlock':
            # New backend caches an integer already specified by the checkpoint's
            # irreducible representation. Retain all saved overlap coefficients.
            expected = module.features_irreps.dim
            if hasattr(module, 'features_dim') and module.features_dim != expected:
                raise ValueError('self-interaction dimension disagrees with checkpoint')
            if not hasattr(module, 'features_dim'):
                module.features_dim = expected
                restored_dimensions.append(name)
    # Checkpoint predates these dispatch attributes. Use the maintained backend's
    # public setters; realspace does not use the new periodic output permutation.
    features.set_pbc_handling('realspace')
    coulomb.set_pbc_handling('realspace')
    original_precompute = features.precompute_geometry
    original_dynamic = features.forward_dynamic
    original_coulomb = coulomb.forward
    original_forward = model.forward
    calls = {'omitted_unused_kgrids': 0}

    def precompute(self, *, force_pbc_evaluator=False, **kwargs):
        require_isolated(kwargs['pbc'], force_pbc_evaluator)
        return original_precompute(**kwargs)

    def dynamic(self, *, cache, source_feats, pbc):
        require_isolated(pbc)
        if source_feats.ndim != 3 or source_feats.shape[1] != 1:
            raise ValueError('expected singleton density channel [N,1,multipoles]')
        return original_dynamic(cache=cache, source_feats=source_feats[:, 0, :])

    def energy(self, *, force_pbc_evaluator=False, **kwargs):
        require_isolated(kwargs['pbc'], force_pbc_evaluator)
        return original_coulomb(**kwargs)

    def unused_grid(cutoff, cell_vectors, r_cell_vectors):
        # Both backend realspace dispatchers ignore all four reciprocal tensors.
        # Empty tensors express absence of a grid, never a zero energy correction.
        calls['omitted_unused_kgrids'] += 1
        return (cell_vectors.new_empty((0, 3)), cell_vectors.new_empty((0,)),
                torch.empty((0,), dtype=torch.long, device=cell_vectors.device),
                cell_vectors.new_empty((0,)))

    @contextmanager
    def without_unused_grid():
        original_grid = extensions.compute_k_vectors_flat
        extensions.compute_k_vectors_flat = unused_grid
        try:
            yield
        finally:
            extensions.compute_k_vectors_flat = original_grid

    def forward(self, data, *args, **kwargs):
        require_isolated(data['pbc'], kwargs.get('use_pbc_evaluator', False))
        if features.pbc_handling != 'realspace' or coulomb.pbc_handling != 'realspace':
            raise ValueError('electrostatic dispatcher changed after adapter setup')
        with without_unused_grid():
            return original_forward(data, *args, **kwargs)

    features.precompute_geometry = MethodType(precompute, features)
    features.forward_dynamic = MethodType(dynamic, features)
    coulomb.forward = MethodType(energy, coulomb)
    model.forward = MethodType(forward, model)
    model._alquemia_realspace_adapter = ADAPTER_ID
    return {'adapter_id': ADAPTER_ID, 'features': features.pbc_handling,
            'energy': coulomb.pbc_handling, 'counters': calls,
            'restored_dimension_metadata': restored_dimensions,
            'kernel_changes': False, 'learned_parameter_changes': False}
