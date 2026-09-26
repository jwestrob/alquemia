"""Finite, fixed-charge electrostatic inputs and their coordinate chain rule.

All quantities are atomic units: coordinates bohr, charges e, phi Eh/e,
field Eh/(e bohr), energy Eh, forces Eh/bohr. The potential is zero at
infinity; an optional common constant gauge shifts phi but never the field.
No periodic images, dielectric, charge response, MM internal energy, van der
Waals terms, link-atom mapping, or intrinsic electronic forces are implied.
"""
import hashlib
import json

import numpy as np

BOHR_ANGSTROM = 0.529177210903


def _inputs(qm_bohr, mm_bohr, mm_charges):
    x, y, q = (np.asarray(a, dtype=np.float64) for a in
               (qm_bohr, mm_bohr, mm_charges))
    if x.ndim != 2 or x.shape[1] != 3 or y.ndim != 2 or y.shape[1] != 3:
        raise ValueError("Coordinates require shape (N,3), in bohr")
    if q.shape != (len(y),) or not all(np.isfinite(a).all() for a in (x, y, q)):
        raise ValueError("Finite coordinates and one finite charge per MM atom required")
    return x, y, q


def _pairs(x, y, q, chunk_size):
    if not isinstance(chunk_size, int) or chunk_size < 1:
        raise ValueError("chunk_size must be a positive integer")
    for start in range(0, len(y), chunk_size):
        d = x[:, None, :] - y[None, start:start + chunk_size, :]
        r2 = np.sum(d*d, axis=2)
        if np.any(r2 == 0):
            raise ValueError("Overlapping QM/MM centers: boundary/exclusions must be resolved")
        invr = 1 / np.sqrt(r2)
        yield start, d, invr, q[start:start + chunk_size]


def potential_field(qm_bohr, mm_bohr, mm_charges, *, gauge=0.0, chunk_size=512):
    """Return phi_i=sum_j q_j/r_ij + gauge and F_i=sum_j q_j*d_ij/r_ij^3."""
    x, y, q = _inputs(qm_bohr, mm_bohr, mm_charges)
    if not np.isfinite(gauge):
        raise ValueError("Finite gauge required")
    phi = np.full(len(x), gauge, dtype=float)
    field = np.zeros_like(x)
    for _, d, ir, charges in _pairs(x, y, q, chunk_size):
        phi += np.sum(ir * charges, axis=1)
        field += np.sum(d * (ir**3 * charges)[:, :, None], axis=1)
    return phi, field


def chain_rule_forces(qm_bohr, mm_bohr, mm_charges, denergy_dphi,
                      denergy_dfield, *, chunk_size=512):
    """Coordinate forces due ONLY to E's dependence on phi and field.

    Add the returned QM forces to intrinsic model forces evaluated at fixed
    phi/field. Add MM forces to the separately accounted MM/cross terms.
    Inputs are derivatives of the complete candidate electronic energy,
    not learned charges/dipoles unless their derivative equivalence is proven.
    Fixed MM charges have no coordinate response beyond their moving centers.
    Physical cap/source chain-rule projection remains the caller's responsibility.
    """
    x, y, q = _inputs(qm_bohr, mm_bohr, mm_charges)
    a, b = np.asarray(denergy_dphi, float), np.asarray(denergy_dfield, float)
    if a.shape != (len(x),) or b.shape != x.shape or not all(
            np.isfinite(v).all() for v in (a, b)):
        raise ValueError("Finite dE/dphi (N,) and dE/dfield (N,3) required")
    fx, fy = np.zeros_like(x), np.zeros_like(y)
    for start, d, ir, charges in _pairs(x, y, q, chunk_size):
        bd = np.sum(b[:, None, :] * d, axis=2)
        gradient = charges[None, :, None] * (
            -a[:, None, None] * d * ir[:, :, None]**3
            + b[:, None, :] * ir[:, :, None]**3
            - 3*d * (bd * ir**5)[:, :, None])
        fx -= gradient.sum(axis=1)
        fy[start:start + len(charges)] += gradient.sum(axis=0)
    return fx, fy


def charged_gauge_shift(total_charge_e, potential_shift_au):
    """Expected Eh shift for a charge-inclusive electronic energy: Q*delta_phi.

    Requires electron AND nuclear coupling in the declared model. This identity
    does not establish that an executable implements this convention.
    """
    values = np.asarray([total_charge_e, potential_shift_au], float)
    if not np.isfinite(values).all():
        raise ValueError("Finite charge and gauge shift required")
    return float(values.prod())


def evaluation_fingerprint(*, state, geometry, environment, model, protocol):
    """Hash complete caller-supplied JSON identities, rejecting missing scopes.

    State includes element ordering, charge, physical multiplicity, backend spin
    and ECP convention. Geometry includes ordered physical/cap maps and coordinates
    or their content hashes. Environment includes coordinates, charges, parameters,
    exclusions, phi/field overrides, gauge and solvent settings. Model includes
    code/checkpoint hashes and units; protocol includes all numerical settings.
    The function cannot detect an omitted setting inside these records: adapters
    must construct these fields, not pass a protein ID or an XYZ hash alone.
    """
    payload = dict(state=state, geometry=geometry, environment=environment,
                   model=model, protocol=protocol)
    if any(not isinstance(v, dict) or not v for v in payload.values()):
        raise ValueError("All five fingerprint scopes must be nonempty dictionaries")
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                     allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()
