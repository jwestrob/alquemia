"""Source-only native/local force composition; no executable dispatch or solve."""
import torch


def combine(kernel, xyz, geometry, parameters, solved, native, K):
    """Caller must bind native/solved to the same immutable state and source.

    native dictionary: other_energy, baseline_gradient[N,3], supplied_u/v[3S].
    The native gradient already includes permanent/nonpolar derivatives and
    fixed-u/v derivatives of A0, Ed0 and Ep0. It is NOT a frozen-A gradient.
    """
    if not xyz.requires_grad:
        raise ValueError('xyz must require gradients for the explicit local derivative')
    u, v = solved['u'].detach(), solved['v'].detach()
    uread=torch.as_tensor(native['supplied_u'],dtype=xyz.dtype)
    vread=torch.as_tensor(native['supplied_v'],dtype=xyz.dtype)
    if not torch.equal(u,uread) or not torch.equal(v,vread):
        raise ValueError('native readback differs from exact supplied corrected dual solution')
    d=kernel.correction(xyz,geometry,native['site_atom_indices'],parameters)
    # Envelope derivative at fixed independent u,v. Differentiate this scalar,
    # not the reduced frozen-baseline polarization energy.
    explicit=d['V']+0.5*K*(u@d['H']@v-(u+v)@d['h'])
    local_gradient=torch.autograd.grad(explicit,xyz)[0]
    baseline_gradient=torch.as_tensor(native['baseline_gradient'],dtype=xyz.dtype)
    if baseline_gradient.shape!=xyz.shape or not bool(torch.isfinite(baseline_gradient).all()):
        raise ValueError('invalid native baseline contraction gradient')
    gradient=baseline_gradient+local_gradient
    energy=float(native['other_energy'])+float(d['V'].detach())+float(solved['polarization_energy'].detach())
    return {'energy_kcal':energy,'gradient_kcal_per_A':gradient,'force_kcal_per_A':-gradient,
            'local_explicit_gradient_kcal_per_A':local_gradient,
            'native_baseline_gradient_kcal_per_A':baseline_gradient,
            'scope':'complete_vacuum_chain_rule_pending_native_qualification'}
