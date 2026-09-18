"""Task-local group constraints for the pinned native POLAR restoration blocks.

This changes the scientific model. It is not an optimization of native POLAR.
The installed implementation and all learned weights remain untouched.
"""
from __future__ import annotations
import inspect
import textwrap
import types
from affordable_common import InvalidArtifact, record, write_new

ADAPTER = 'polar_formal_group_Fukui_restoration_v1'
WEIGHT_RATIO_MIN = 1e-10
CHARGE_TOLERANCE = 1e-5


class GroupRestoration:
    def __init__(self, atom_groups, group_charges):
        if (not atom_groups or not group_charges or
                any(type(g) is not int for g in atom_groups) or
                set(atom_groups) != set(range(len(group_charges))) or
                any(type(q) is not int for q in group_charges)):
            raise InvalidArtifact('invalid fixed formal-charge group inventory')
        self.atom_groups = list(atom_groups)
        self.group_charges = list(group_charges)
        self.trace = []

    def __call__(self, density, weights):
        import torch
        n, channels, _ = density.shape
        if (n != len(self.atom_groups) or channels != 2 or weights.shape != (n, 2)
                or density.dtype != torch.float64 or weights.dtype != torch.float64
                or not torch.isfinite(density).all() or not torch.isfinite(weights).all()):
            raise InvalidArtifact('unsupported group-restoration tensor/state')
        ids = torch.tensor(self.atom_groups, dtype=torch.long, device=density.device)
        target = torch.tensor(self.group_charges, dtype=density.dtype, device=density.device)[:, None] / 2
        shape = (len(self.group_charges), 2)
        sums = torch.zeros(shape, dtype=torch.float64, device=density.device).index_add(0, ids, weights)
        absolute = torch.zeros_like(sums).index_add(0, ids, weights.abs())
        if (absolute == 0).any() or (sums.abs() <= WEIGHT_RATIO_MIN * absolute).any():
            raise InvalidArtifact('near-cancelling or zero group Fukui weights; no clipping permitted')
        raw = torch.zeros_like(sums).index_add(0, ids, density[:, :, 0])
        correction = weights / sums[ids] * (target - raw)[ids]
        result = density.clone()
        result[:, :, 0] = density[:, :, 0] + correction
        restored = torch.zeros_like(sums).index_add(0, ids, result[:, :, 0])
        error = float((restored - target).abs().max().detach())
        if error > CHARGE_TOLERANCE:
            raise InvalidArtifact('group spin-channel closure failed')
        self.trace.append({'stage': len(self.trace), 'group_raw_charge_per_spin_e': raw.detach().cpu().tolist(),
                           'group_restored_charge_per_spin_e': restored.detach().cpu().tolist(),
                           'maximum_channel_charge_error_e': error,
                           'minimum_weight_ratio': float((sums.abs() / absolute).min().detach()),
                           'charge_correction_L1_e': float(correction.abs().sum().detach())})
        return result

    def receipt(self):
        if len(self.trace) != 3:
            raise InvalidArtifact('expected exactly initial plus two native restoration stages')
        return {'adapter_id': ADAPTER, 'atom_group_indices': self.atom_groups,
                'group_charges_e': self.group_charges, 'stages': self.trace,
                'weight_ratio_minimum': WEIGHT_RATIO_MIN, 'charge_tolerance_e': CHARGE_TOLERANCE,
                'weights_changed': False, 'energy_readouts_changed': False,
                'native_global_model_reproduced_only_in_one_group_limit': True}


def configure(calc, atom_groups, group_charges, native_source_sha256, output):
    """Replace only the two explicit normalization blocks in a copied method.

    Both full source and copied method are pinned and saved. This avoids editing
    site-packages or replacing the physical graph's batch IDs with charge groups.
    """
    model = calc.models[0]
    if len(calc.models) != 1 or hasattr(model, '_alquemia_group_restore'):
        raise InvalidArtifact('one unmodified checkpoint instance required')
    if hasattr(model, '_alquemia_realspace_adapter') or hasattr(model, '_alquemia_blocked_kernel'):
        raise InvalidArtifact('install group constraints before the existing realspace/memory wrappers')
    original = type(model).forward
    path = inspect.getsourcefile(original)
    native = record(path)
    if native['sha256'] != native_source_sha256:
        raise InvalidArtifact('native POLAR source differs from inspected version')
    source = textwrap.dedent(inspect.getsource(original))
    start = source.index('    fukui_norm = scatter_sum(')
    end = source.index('    # print("spin_charge_density", spin_charge_density)', start)
    source = source[:start] + (
        '    spin_charge_density = self._alquemia_group_restore(spin_charge_density, fukui_sources)\n'
    ) + source[end:]
    start = source.index('        fukui_norm2 = scatter_sum(')
    end = source.index('\n    total_energy = e0 + inter_e', start)
    source = source[:start] + (
        '        spin_charge_density = self._alquemia_group_restore(spin_charge_density, current_fukui_sources)\n'
    ) + source[end:]
    if source.count('self._alquemia_group_restore(') != 2 or 'Q_p_S' in source or 'Q_m_S' in source:
        raise InvalidArtifact('native restoration block extraction changed unexpectedly')
    # Keep the native module globals live: the existing realspace wrapper
    # temporarily redirects compute_k_vectors_flat there. A copied namespace
    # would bypass that memory fix and allocate an unused reciprocal grid.
    source = source.replace('def forward(', 'def _alquemia_group_forward(', 1)
    path = output / 'group_forward.py'
    path.write_text(source)
    namespace = original.__globals__
    exec(compile(source, str(path), 'exec'), namespace)
    observer = GroupRestoration(atom_groups, group_charges)
    model._alquemia_group_restore = observer
    model.forward = types.MethodType(namespace['_alquemia_group_forward'], model)
    receipt = {'adapter_id': ADAPTER, 'native_source': native, 'copied_forward': record(path),
               'adapter_source': record(__file__), 'learned_parameters_modified': False,
               'native_restoration_blocks_replaced': 2, 'physical_graph_batch_unchanged': True}
    write_new(output / 'group_adapter.json', receipt)
    return observer, receipt
