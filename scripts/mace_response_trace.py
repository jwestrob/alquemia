"""Read-only charge/field tracing of pinned real MACE primary evaluations."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new

SCHEMA = 'alquemia.mace_response_trace.v1'
TRACE_ID = 'polar0316_charge_update_observer_v1'


class ChargeTrace:
    def __init__(self, model):
        self.arrays = {}
        self.handles = []
        self.steps = model.num_recursion_steps
        self.layers = len(model.lr_source_maps)
        self.model_type = type(model).__name__
        for i, module in enumerate(model.lr_source_maps):
            self.attach(module, f'local_source_{i}')
        self.attach(model.fukui_source_map, 'initial_weights')
        for i, module in enumerate(model.field_dependent_charges_maps):
            self.attach(module, f'update_{i}', ('local_charges', 'potential_features'))
        self.attach(model.local_electron_energy, 'electron_readout',
                    ('charges_0', 'charges_induced'))

    def attach(self, module, name, inputs=()):
        def hook(_module, _args, kwargs, output):
            key = name + '_output'
            if key in self.arrays:
                raise InvalidArtifact(f'trace module called repeatedly: {name}')
            self.arrays[key] = output.detach().cpu().numpy().copy()
            for field in inputs:
                self.arrays[name + '_' + field] = kwargs[field].detach().cpu().numpy().copy()
            # Returning None preserves the module's exact output tensor/graph.
        self.handles.append(module.register_forward_hook(hook, with_kwargs=True))

    def finish(self, output, total_charge, final_density):
        for handle in self.handles:
            handle.remove()
        n = len(final_density)
        local = sum(self.arrays[f'local_source_{i}_output'].reshape(n, 2, -1)
                    for i in range(self.layers))
        stages = []
        previous = local
        for stage in range(self.steps + 1):
            if stage == 0:
                provisional = local
                weights = self.arrays['initial_weights_output']
            else:
                i = stage - 1
                incoming = self.arrays[f'update_{i}_local_charges'].reshape(n, 2, -1)
                raw = self.arrays[f'update_{i}_output']
                provisional = incoming + raw[:, :-2].reshape(n, 2, -1)
                weights = raw[:, -2:]
            if weights.shape != (n, 2):
                raise InvalidArtifact('unexpected charge-restoration weight shape')
            sums = weights.sum(axis=0)
            absolute = np.abs(weights).sum(axis=0)
            normalized = weights / np.where(sums == 0, 1., sums)
            raw_charge = provisional[:, :, 0].sum(axis=0)
            correction = normalized * (total_charge / 2 - raw_charge)
            restored = provisional.copy()
            restored[:, :, 0] += correction
            expected = (self.arrays[f'update_{stage}_local_charges'] if stage < self.steps
                        else self.arrays['electron_readout_charges_induced']).reshape(n, 2, -1)
            residual = float(np.max(np.abs(restored - expected)))
            if residual > 1e-8:
                raise InvalidArtifact(f'charge trace reconstruction mismatch {residual}')
            self.arrays[f'stage_{stage}_restored_density_internal'] = expected.copy()
            self.arrays[f'stage_{stage}_restoration_correction_e'] = correction
            self.arrays[f'stage_{stage}_normalized_weights'] = normalized
            stages.append({'stage': stage, 'raw_charge_per_spin_e': raw_charge.tolist(),
                           'restored_charge_per_spin_e': expected[:, :, 0].sum(axis=0).tolist(),
                           'raw_weight_sum': sums.tolist(), 'raw_weight_L1': absolute.tolist(),
                           'weight_sum_cancellation_ratio': (np.abs(sums) / np.where(absolute == 0, 1., absolute)).tolist(),
                           'normalized_weight_L1': np.abs(normalized).sum(axis=0).tolist(),
                           'normalized_weight_max_abs': np.abs(normalized).max(axis=0).tolist(),
                           'restoration_charge_L1_e': np.abs(correction).sum(axis=0).tolist(),
                           'restored_charge_min_e': float(expected[:, :, 0].sum(axis=1).min()),
                           'restored_charge_max_e': float(expected[:, :, 0].sum(axis=1).max()),
                           'total_charge_update_L1_e': float(np.abs(expected[:, :, 0] - previous[:, :, 0]).sum()),
                           'reconstruction_max_abs': residual})
            previous = expected
        scalar_error = float(np.max(np.abs(previous[:, :, 0].sum(axis=1) - final_density[:, 0])))
        if scalar_error > 1e-8:
            raise InvalidArtifact('trace final density does not match saved result')
        path = Path(output) / 'charge_trace_arrays.npz'
        np.savez_compressed(path, **self.arrays)
        summary = {'trace_id': TRACE_ID, 'model_type': self.model_type,
                   'layers': self.layers, 'recursion_steps': self.steps,
                   'atoms': n, 'stages': stages, 'final_scalar_error_e': scalar_error,
                   'arrays': record(path), 'spin_policy': 'closed-shell multiplicity 1',
                   'observer_changes_energy_or_gradient': False}
        write_new(Path(output) / 'charge_trace.json', summary)
        return record(Path(output) / 'charge_trace.json')


def prepare(source_collection, agreement, output):
    from mace_hybrid import dry_run, TASK_IDS
    c = read_json(source_collection)
    sm = read_json(verify(c['manifest']))
    dry_run(verify(c['manifest']))
    if c['status'] != 'complete' or sm['schema_version'] != 'alquemia.mace_analytic.v1':
        raise InvalidArtifact('complete analytic primary reference required')
    m = copy.deepcopy(sm)
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    impl = out / 'implementation'
    impl.mkdir()
    pins = {}
    for name in (*sm['implementation'], 'mace_response_trace.py'):
        shutil.copyfile(Path(__file__).with_name(name), impl / name)
        pins[name] = record(impl / name)
    names = TASK_IDS + [f'1h4i_full_{metal}_primary' for metal in ('La', 'Ca')]
    tasks = [copy.deepcopy(next(t for t in sm['tasks'] if t['task_id'] == name)) for name in names]
    m.update(schema_version=SCHEMA, protocol_id=sm['protocol_id'] + '_charge_trace_v1',
             implementation=pins, tasks=tasks, agreement=record(agreement),
             trace_reference=record(source_collection),
             trace_tolerances={'energy_eV': 1e-6, 'force_eV_A': 1e-6, 'density': 1e-8},
             run_inventory={'distinct_mace_energy_force_calls': 6, 'new_DFT_endpoints': 0})
    m['model']['response_trace'] = TRACE_ID
    for t in tasks:
        t.pop('cache_key')
        t['cache_key'] = cache_key({'task': t, 'model': m['model'], 'software': m['software'], 'implementation': pins})
    write_new(out / 'manifest.json', m)
    return dry_run(out / 'manifest.json')


def validate(m):
    from mace_hybrid import dry_run, TASK_IDS
    source = read_json(verify(m['trace_reference']))
    sm = read_json(verify(source['manifest']))
    dry_run(verify(source['manifest']))
    expected_model = copy.deepcopy(sm['model'])
    expected_model['response_trace'] = TRACE_ID
    if (m['model'] != expected_model or m['software'] != sm['software'] or
            m['tolerances'] != sm['tolerances'] or
            m['protocol_id'] != sm['protocol_id'] + '_charge_trace_v1'):
        raise InvalidArtifact('trace changes physical model or software')
    if m['trace_tolerances'] != {'energy_eV': 1e-6, 'force_eV_A': 1e-6, 'density': 1e-8}:
        raise InvalidArtifact('trace equivalence tolerances changed')
    names = TASK_IDS + [f'1h4i_full_{metal}_primary' for metal in ('La', 'Ca')]
    if [t['task_id'] for t in m['tasks']] != names:
        raise InvalidArtifact('trace task inventory changed')
    for t in m['tasks']:
        original = next(x for x in sm['tasks'] if x['task_id'] == t['task_id'])
        if {k:v for k,v in t.items() if k != 'cache_key'} != {k:v for k,v in original.items() if k != 'cache_key'}:
            raise InvalidArtifact('trace changes physical input')


def add_analysis(c, m):
    source = read_json(verify(m['trace_reference']))
    tol = m['trace_tolerances']
    checks = []
    for name, new in c['rows'].items():
        old = source['rows'][name]
        if new['status'] != 'computed':
            checks.append({'task_id': name, 'pass': False, 'status': 'unavailable'})
            continue
        de = new['energy_eV'] - old['energy_eV']
        df = float(np.max(np.abs(np.load(verify(new['forces'])) - np.load(verify(old['forces'])))))
        dd = float(np.max(np.abs(np.load(verify(new['density_coefficients'])) - np.load(verify(old['density_coefficients'])))))
        trace = read_json(verify(new['charge_trace']))
        verify(trace['arrays'])
        checks.append({'task_id': name, 'energy_delta_eV': de,
                       'force_delta_eV_A': df, 'density_delta': dd,
                       'pass': abs(de) <= tol['energy_eV'] and df <= tol['force_eV_A'] and dd <= tol['density'],
                       'charge_trace': new['charge_trace']})
    c['trace_equivalence'] = checks


def core_gate(manifest):
    from mace_hybrid import collect, TASK_IDS
    c = collect(manifest)
    checks = [r for r in c['trace_equivalence'] if r['task_id'] in TASK_IDS]
    return {'status': 'pass' if len(checks) == 4 and all(r['pass'] for r in checks) else 'failed', 'checks': checks}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('source-collection', 'agreement', 'output'):
        p.add_argument('--' + key, required=True)
    a = p.parse_args()
    print(json.dumps(prepare(a.source_collection, a.agreement, a.output), indent=2))


if __name__ == '__main__':
    main()
