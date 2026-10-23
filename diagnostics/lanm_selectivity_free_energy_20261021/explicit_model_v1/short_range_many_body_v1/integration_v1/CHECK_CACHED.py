"""Exercise the new Hamiltonian on immutable real baseline exports, without native calls.

Zero-correction parity uses inherited export bounds. Fixed nonzero probes are
implementation checks, not selected/fitted candidates. Frozen A coordinate
gradients are explicitly partial, never physical total forces.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import torch


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def load_kernel(path):
    spec = importlib.util.spec_from_file_location('srmb_kernel_under_test', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def scalar(x):
    return float(x.detach())


def run_one(path, mod):
    raw = json.loads(Path(path).read_text())
    op, geom = raw['native_operator'], raw['geometry']
    if not raw['compatible'] or op is None:
        raise ValueError('fixture lacks an admitted compatible baseline')
    tensor = lambda x: torch.tensor(x, dtype=torch.float64)
    A_raw = tensor(op['A0_Aminus3'])
    snapshots = raw['native_qualification']['snapshots']
    bound = min(s['numerical_bounds']['arithmetic'] for s in snapshots)
    asym = scalar(torch.max(torch.abs(A_raw-A_raw.T)))
    if asym > bound:
        raise ValueError('raw A asymmetry exceeds inherited numerical export bound')
    # Explicit, recorded projection of admitted floating-point asymmetry only.
    A, symmetry_receipt = mod.symmetrize_native_operator(A_raw,
        max_absolute_asymmetry=bound,
        admission_authority={'fixture': pin(path), 'field': 'native_qualification'})
    xyz = tensor(geom['xyz_A']).requires_grad_()
    geometry = mod.Geometry(geom['metal_atom'], tuple(geom['atomic_numbers']),
                            tuple(mod.Ligand.from_dict(x) for x in geom['ligands']))
    base = mod.Baseline(A, tensor(op['Ed_e_per_A2']), tensor(op['Ep_e_per_A2']),
        tuple(op['site_atom_indices']), op['K_kcal_A_per_mol_e2'],
        tensor(op['nonpolar_energy_kcal']), tensor(op['permanent_mu_eA']),
        authority=raw['baseline_id'])
    zero = mod.Parameters()
    out = mod.evaluate(xyz, geometry, zero, base)
    checks = {}

    def check(name, error, tolerance):
        error, tolerance = float(error), float(tolerance)
        checks[name] = {'error': error, 'tolerance': tolerance, 'passed': error <= tolerance}
        if not checks[name]['passed']:
            raise ValueError(f'{name}: {error} > {tolerance}')

    ub = max(s['numerical_bounds']['native_vs_dense_u'] for s in snapshots)
    vb = max(s['numerical_bounds']['native_vs_dense_v'] for s in snapshots)
    check('native_u', scalar(torch.max(torch.abs(out['u']-tensor(op['u_native_eA'])))), ub)
    check('native_v', scalar(torch.max(torch.abs(out['v']-tensor(op['v_native_eA'])))), vb)
    pe = raw['native_qualification']['polarization_energy']
    energy_bound = pe['print_half_unit'] + .5*base.K*scalar(base.Ep0.abs().sum())*ub
    check('native_total_energy', abs(scalar(out['energy'])-raw['baseline_observables']['energy_kcal_mol']), energy_bound)
    check('native_polarization', abs(scalar(out['polarization_energy'])-pe['native_contraction']), energy_bound)
    B = torch.eye(3, dtype=torch.float64).repeat(len(base.site_atom_indices), 1)
    native_mu = base.permanent_mu_eA + .5*B.T@(tensor(op['u_native_eA'])+tensor(op['v_native_eA']))
    check('native_total_moment', scalar(torch.max(torch.abs(out['mu_eA']-native_mu))), .5*len(base.site_atom_indices)*(ub+vb))
    check('zero_V', abs(scalar(out['delta']['V'])), 0)
    check('zero_h', scalar(out['delta']['h'].abs().max()), 0)
    check('zero_H', scalar(out['delta']['H'].abs().max()), 0)
    if out['coordinate_derivative_scope'] != 'partial_frozen_baseline_not_total_force':
        raise ValueError('frozen cache misrepresented as total-force provider')

    # Probe magnitudes fixed before execution; these do not seek quantum agreement.
    kinds = sorted({l.kind for l in geometry.ligands})
    probes = {
        'mechanical': mod.Parameters(pairs={k: mod.Pair(va=.01) for k in kinds}),
        'drive': mod.Parameters(pairs={k: mod.Pair(hf=.0001) for k in kinds}),
        'stiffness': mod.Parameters(pairs={k: mod.Pair(kff=.0001) for k in kinds}),
    }
    if len(geometry.ligands) >= 2:
        triple_kinds = {tuple(sorted((a.kind, b.kind)))
            for i, a in enumerate(geometry.ligands) for b in geometry.ligands[i+1:]}
        probes.update({
            'triple_mechanical': mod.Parameters(triples={k: mod.Triple(vf=.01) for k in triple_kinds}),
            'triple_drive': mod.Parameters(triples={k: mod.Triple(h_neighbor=.0001) for k in triple_kinds}),
            'triple_stiffness': mod.Parameters(triples={k: mod.Triple(k_cross=.0001) for k in triple_kinds}),
        })
    probe_records = {}
    for label, parameters in probes.items():
        field = torch.zeros(3, dtype=torch.float64, requires_grad=True)
        p = mod.evaluate(xyz, geometry, parameters, base, field)
        gradient = torch.autograd.grad(p['energy'], field, create_graph=True)[0]
        hessian = torch.stack([torch.autograd.grad(gradient[i], field, retain_graph=True)[0] for i in range(3)])
        check(label+'_conjugate_moment', scalar(torch.max(torch.abs(-gradient/base.K-p['mu_eA']))), 1e-11)
        check(label+'_conjugate_susceptibility', scalar(torch.max(torch.abs(-hessian/base.K-p['chi_A3']))), 1e-10)
        if label.endswith('mechanical'):
            check(label+'_preserves_moment', scalar(torch.max(torch.abs(p['mu_eA']-out['mu_eA']))), 0)
            check(label+'_preserves_susceptibility', scalar(torch.max(torch.abs(p['chi_A3']-out['chi_A3']))), 0)
        if label.endswith('drive'):
            check(label+'_preserves_susceptibility', scalar(torch.max(torch.abs(p['chi_A3']-out['chi_A3']))), 0)
        probe_records[label] = {
            'energy_shift_kcal_mol': scalar(p['energy']-out['energy']),
            'moment_shift_eA': (p['mu_eA']-out['mu_eA']).detach().tolist(),
            'susceptibility_shift_A3': (p['chi_A3']-out['chi_A3']).detach().tolist(),
            'minimum_eigenvalue': scalar(p['eigenvalues'][0]),
        }
    return {'fixture': pin(path), 'id': raw['id'], 'baseline_id': raw['baseline_id'],
        'symmetry_projection': symmetry_receipt,
        'checks': checks, 'probes': probe_records, 'coordinate_derivative_scope': out['coordinate_derivative_scope']}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kernel', required=True)
    ap.add_argument('--fixture', action='append', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    torch.set_num_threads(1)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=False)
    start, cpu = time.monotonic(), time.process_time()
    kernel_pin = pin(args.kernel)
    result = {'status': 'failed', 'implementation': pin(__file__), 'kernel': kernel_pin,
        'rows': [], 'error': None, 'new_native_calls': 0, 'new_quantum_calls': 0,
        'coefficients_fitted': False, 'physical_accuracy_qualified': False,
        'full_real_coordinate_forces_qualified': False}
    try:
        module = load_kernel(args.kernel)
        for fixture in args.fixture:
            result['rows'].append(run_one(fixture, module))
        if pin(args.kernel) != kernel_pin:
            raise ValueError('kernel changed during execution')
        result['status'] = 'passed_cached_baseline_and_new_term_integration'
    except Exception as exc:
        result['error'] = repr(exc)
    result['cost'] = {'wall_seconds': time.monotonic()-start,
                      'CPU_seconds': time.process_time()-cpu, 'threads': 1}
    (out/'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'rows': len(result['rows']), 'error': result['error']}))
    return 0 if result['error'] is None else 1


if __name__ == '__main__':
    raise SystemExit(main())
