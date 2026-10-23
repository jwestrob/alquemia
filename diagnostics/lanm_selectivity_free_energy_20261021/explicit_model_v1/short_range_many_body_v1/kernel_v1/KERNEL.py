"""Finite local energy/response prototype. No fitted parameters or native calls.

Coordinates A; energy kcal/mol; fields e/A^2; dipoles e A; A/H in A^-3.
Only small dense CPU float64 operators are supported. See SPEC.md.
"""
from dataclasses import dataclass, field
from typing import Callable
import math
import torch
import torch.nn.functional as functional


@dataclass(frozen=True)
class Ligand:
    id: str
    kind: str
    atom_indices: tuple
    donor_atoms: tuple
    frame_atoms_by_donor: dict

    @classmethod
    def from_dict(cls, row):
        return cls(row['id'], row['kind'], tuple(row['atom_indices']),
                   tuple(row['donor_atoms']),
                   {int(k): tuple(v) for k, v in row['frame_atoms_by_donor'].items()})


@dataclass(frozen=True)
class Geometry:
    metal_atom: int | None
    atomic_numbers: tuple
    ligands: tuple

    def validate(self, xyz, sites):
        if xyz.dtype != torch.float64 or xyz.device.type != 'cpu' or xyz.shape != (len(self.atomic_numbers), 3):
            raise ValueError('CPU float64 xyz_A[N,3] required')
        if not bool(torch.isfinite(xyz).all()):
            raise ValueError('nonfinite coordinates')
        if self.metal_atom is not None and not 0 <= self.metal_atom < len(self.atomic_numbers):
            raise ValueError('metal atom index')
        if len(set(sites)) != len(sites) or any(i < 0 or i >= len(xyz) for i in sites):
            raise ValueError('invalid response-site atom order')
        used, ids = (set() if self.metal_atom is None else {self.metal_atom}), set()
        for ligand in self.ligands:
            if ligand.kind not in ('water', 'acetate', 'amide') or ligand.id in ids:
                raise ValueError('ligand kind/id')
            ids.add(ligand.id)
            atoms = set(ligand.atom_indices)
            if len(atoms) != len(ligand.atom_indices) or atoms & used or any(i < 0 or i >= len(xyz) for i in atoms):
                raise ValueError('intact ligand atoms overlap/invalid')
            used |= atoms
            expected = 2 if ligand.kind == 'acetate' else 1
            if len(ligand.donor_atoms) != expected or len(set(ligand.donor_atoms)) != expected:
                raise ValueError('one water/amide O or two same-fragment acetate O required')
            if set(ligand.frame_atoms_by_donor) != set(ligand.donor_atoms):
                raise ValueError('donor/frame map mismatch')
            for donor in ligand.donor_atoms:
                frame = ligand.frame_atoms_by_donor[donor]
                if donor not in atoms or donor not in sites or self.atomic_numbers[donor] != 8:
                    raise ValueError('donor must be an oxygen response site in its intact ligand')
                if not frame or len(set(frame)) != len(frame) or donor in frame or not set(frame) <= atoms:
                    raise ValueError('invalid covalent frame atoms')
                if ligand.kind == 'water' and (len(frame) != 2 or any(self.atomic_numbers[j] != 1 for j in frame)):
                    raise ValueError('water frame is its two equivalent H atoms')


@dataclass(frozen=True)
class Pair:
    v0: float = 0.0       # kcal/mol: compact energy
    vr: float = 0.0       # kcal/mol: independent radial slope basis
    va: float = 0.0       # kcal/mol: donor-frame angular energy
    hr: float = 0.0       # e/A^2: metal-to-donor drive
    hf: float = 0.0       # e/A^2: donor-frame drive
    krr: float = 0.0      # A^-3: radial donor stiffness
    kff: float = 0.0      # A^-3: frame donor stiffness


@dataclass(frozen=True)
class Triple:
    v0: float = 0.0       # kcal/mol: intact-ligand nonadditivity
    va: float = 0.0       # kcal/mol: inter-ligand angle
    vf: float = 0.0       # kcal/mol: symmetric directed ligand-frame orientation
    h_neighbor: float = 0.0   # e/A^2: cross-directed drive
    k_cross: float = 0.0      # A^-3: symmetric cross-site response


@dataclass(frozen=True)
class Parameters:
    pairs: dict = field(default_factory=dict)     # ligand kind -> Pair
    triples: dict = field(default_factory=dict)   # sorted two-kind tuple -> Triple
    cutoff_A: float = 6.0
    radial_origin_A: float = 2.5
    contact_A: float = 0.1
    frame_softness: float = 0.1  # dimensionless; prevents singular cancelling bonds

    def validate(self):
        kinds = {'water', 'acetate', 'amide'}
        if set(self.pairs) - kinds or any(len(k) != 2 or tuple(sorted(k)) != k or set(k) - kinds for k in self.triples):
            raise ValueError('parameter kind keys; triple keys must be sorted')
        if not (self.cutoff_A > self.radial_origin_A > 0 and self.contact_A > 0 and self.frame_softness > 0):
            raise ValueError('positive fixed geometric scales required')
        for values in list(self.pairs.values()) + list(self.triples.values()):
            if any(not math.isfinite(float(v)) for v in vars(values).values()):
                raise ValueError('nonfinite coefficient')


@dataclass
class Baseline:
    A0: torch.Tensor
    Ed0: torch.Tensor
    Ep0: torch.Tensor
    site_atom_indices: tuple
    K: float                    # literal electric/dielec from native authority
    U_other0: torch.Tensor | None = None  # all baseline energy except polarization
    permanent_mu_eA: torch.Tensor | None = None
    coordinate_derivatives_complete: bool = False
    authority: str = 'manufactured'


def compact_weight(r2, cutoff_A):
    """(1-r^2/Rc^2)^3_+: C2, compact, finite at contact."""
    return torch.clamp(1.0-r2/(cutoff_A*cutoff_A), min=0.0)**3


def _embedded(vector, site, dim):
    return functional.pad(vector, (3*site, dim-3*site-3))


def correction(xyz, geometry, sites, parameters):
    """Pre-relaxation pair/triple DeltaV,h,H. No fragments split into O pairs."""
    geometry.validate(xyz, sites)
    parameters.validate()
    dim = 3*len(sites)
    site_map = {atom: i for i, atom in enumerate(sites)}
    zero = xyz.sum()*0.0
    V = zero
    h = torch.zeros(dim, dtype=xyz.dtype) + zero
    H = torch.zeros((dim, dim), dtype=xyz.dtype) + zero
    frames = []
    counts = {'intact_pairs': len(geometry.ligands), 'intact_triples': 0, 'donor_pair_contributions': 0, 'cross_ligand_donor_contributions': 0}
    if geometry.metal_atom is None:
        counts['intact_pairs'] = 0
        return {'V': V, 'h': h, 'H': H, 'counts': counts}
    delta2 = parameters.contact_A**2
    for ligand in geometry.ligands:
        p = parameters.pairs.get(ligand.kind, Pair())
        local = []
        for donor in ligand.donor_atoms:
            R = xyz[donor]-xyz[geometry.metal_atom]
            r2 = torch.dot(R, R)
            n = R/torch.sqrt(r2+delta2)
            bonds = xyz[list(ligand.frame_atoms_by_donor[donor])]-xyz[donor]
            directions = bonds/torch.sqrt((bonds*bonds).sum(dim=1, keepdim=True)+delta2)
            bisector = directions.sum(dim=0)
            b = bisector/torch.sqrt(torch.dot(bisector, bisector)+parameters.frame_softness**2)
            w = compact_weight(r2, parameters.cutoff_A)
            a = torch.dot(n, b)
            rho = (r2-parameters.radial_origin_A**2)/parameters.radial_origin_A**2
            en, eb = _embedded(n, site_map[donor], dim), _embedded(b, site_map[donor], dim)
            V = V + w*(p.v0+p.vr*rho+p.va*a)
            h = h + w*(p.hr*en+p.hf*eb)
            H = H + w*(p.krr*torch.outer(en, en)+p.kff*torch.outer(eb, eb))
            local.append((w, n, en, site_map[donor], a))
            counts['donor_pair_contributions'] += 1
        frames.append((ligand.kind, local))
    for i, (kind_i, donors_i) in enumerate(frames):
        for kind_j, donors_j in frames[i+1:]:
            p = parameters.triples.get(tuple(sorted((kind_i, kind_j))), Triple())
            counts['intact_triples'] += 1
            for wi, ni, ei, si, ai in donors_i:
                for wj, nj, ej, sj, aj in donors_j:
                    w = wi*wj
                    c = torch.dot(ni, nj)
                    # Embed neighbor direction at the receiving donor site.
                    cross_i, cross_j = _embedded(nj, si, dim), _embedded(ni, sj, dim)
                    V = V + w*(p.v0+p.va*c+p.vf*0.5*(ai+aj))
                    h = h + p.h_neighbor*w*(cross_i+cross_j)
                    H = H + p.k_cross*w*(torch.outer(ei, ej)+torch.outer(ej, ei))
                    counts['cross_ligand_donor_contributions'] += 1
    return {'V': V, 'h': h, 'H': H, 'counts': counts}


def _validate_baseline(base, dim):
    for name, shape in (('A0', (dim, dim)), ('Ed0', (dim,)), ('Ep0', (dim,))):
        a = getattr(base, name)
        if a.dtype != torch.float64 or a.device.type != 'cpu' or a.shape != shape or not bool(torch.isfinite(a).all()):
            raise ValueError('invalid baseline '+name)
    if dim > 600:
        raise ValueError('small dense prototype limited to 200 response sites')
    if not math.isfinite(base.K) or base.K <= 0:
        raise ValueError('literal positive native K required')
    if not torch.equal(base.A0, base.A0.T):
        raise ValueError('A0 must be exactly symmetric; no implicit symmetrization')
    if base.permanent_mu_eA is not None and (base.permanent_mu_eA.shape != (3,) or not bool(torch.isfinite(base.permanent_mu_eA).all())):
        raise ValueError('permanent dipole shape/values')


def symmetrize_native_operator(raw_A, *, max_absolute_asymmetry, admission_authority):
    """Explicit upstream representation repair under a supplied inherited bound.

    Caller preserves raw array/pins; the returned receipt must accompany use.
    This is never called automatically by evaluate and establishes no new gate.
    """
    if not admission_authority or not math.isfinite(max_absolute_asymmetry) or max_absolute_asymmetry < 0:
        raise ValueError('exact admission authority and nonnegative inherited bound required')
    if raw_A.dtype != torch.float64 or raw_A.device.type != 'cpu' or raw_A.ndim != 2 or raw_A.shape[0] != raw_A.shape[1] or not bool(torch.isfinite(raw_A).all()):
        raise ValueError('invalid raw operator')
    asym = float(torch.max(torch.abs(raw_A-raw_A.T)))
    if asym > max_absolute_asymmetry:
        raise ValueError('raw operator exceeds inherited symmetry admission')
    symmetric = 0.5*(raw_A+raw_A.T)
    return symmetric, {'policy': 'explicit_arithmetic_symmetric_part',
                       'admission_authority': admission_authority,
                       'inherited_max_absolute_asymmetry': max_absolute_asymmetry,
                       'raw_max_absolute_asymmetry': asym,
                       'max_absolute_change': float(torch.max(torch.abs(symmetric-raw_A))),
                       'frobenius_change': float(torch.linalg.vector_norm(symmetric-raw_A))}


def evaluate(xyz, geometry, parameters, baseline, field_native=None, *, baseline_mode='frozen'):
    """Evaluate one state. field_native is uniform e/A^2 in BOTH native channels.

    baseline_mode='provider' is an explicit provider's completeness assertion,
    not an inference from requires_grad. Frozen-cache gradients remain partial.
    """
    if baseline_mode not in ('frozen', 'provider'):
        raise ValueError('baseline mode')
    if baseline_mode == 'provider' and (not baseline.coordinate_derivatives_complete or baseline.U_other0 is None or baseline.permanent_mu_eA is None):
        raise ValueError('total derivative provider must supply complete A/fields/other energy/permanent moment')
    dim = 3*len(baseline.site_atom_indices)
    _validate_baseline(baseline, dim)
    delta = correction(xyz, geometry, baseline.site_atom_indices, parameters)
    G = baseline.A0+delta['H']
    if not torch.equal(G, G.T):
        raise ValueError('response operator symmetry lost')
    eigenvalues = torch.linalg.eigvalsh(G)
    if not bool(torch.isfinite(eigenvalues).all()) or not bool(eigenvalues[0] > 0):
        raise ValueError('full response operator is not numerically SPD')
    chol, info = torch.linalg.cholesky_ex(G)
    if int(info) != 0:
        raise ValueError('full response Cholesky failed')
    B = torch.eye(3, dtype=xyz.dtype).repeat(dim//3, 1)
    F = torch.zeros(3, dtype=xyz.dtype) if field_native is None else field_native
    if F.dtype != torch.float64 or F.device.type != 'cpu' or F.shape != (3,) or not bool(torch.isfinite(F).all()):
        raise ValueError('field_native must be finite CPU float64[3]')
    ed = baseline.Ed0+delta['h']+B@F
    ep = baseline.Ep0+delta['h']+B@F
    solved = torch.cholesky_solve(torch.column_stack((ed, ep, B)), chol)
    u, v = solved[:, 0], solved[:, 1]
    polar = -0.5*baseline.K*torch.dot(u, ep)
    other = xyz.sum()*0 if baseline.U_other0 is None else baseline.U_other0
    permanent = torch.zeros(3, dtype=xyz.dtype) if baseline.permanent_mu_eA is None else baseline.permanent_mu_eA
    energy = other-baseline.K*torch.dot(permanent, F)+delta['V']+polar
    mu_ind = 0.5*B.T@(u+v)
    chi = B.T@solved[:, 2:]
    residual = torch.max(torch.abs(G@solved-torch.column_stack((ed, ep, B))))
    scale = max(1.0, float(torch.max(torch.abs(G)).detach())*float(torch.max(torch.abs(solved)).detach())*dim,
                float(torch.max(torch.abs(torch.column_stack((ed, ep, B))))))
    if float(residual.detach()) > 100*dim*torch.finfo(torch.float64).eps*scale:
        raise ValueError('dense solve residual qualification failed')
    return {'energy': energy, 'polarization_energy': polar, 'u': u, 'v': v,
            'mu_induced_eA': mu_ind, 'mu_eA': mu_ind+permanent, 'chi_A3': chi,
            'G': G, 'delta': delta, 'eigenvalues': eigenvalues, 'solve_residual': residual,
            'energy_scope': 'total_declared_baseline' if baseline.U_other0 is not None else 'polarization_plus_local_correction',
            'moment_scope': 'total_with_declared_permanent_origin' if baseline.permanent_mu_eA is not None else 'induced_only',
            'coordinate_derivative_scope': 'total_by_provider_contract' if baseline_mode == 'provider' else 'partial_frozen_baseline_not_total_force',
            'baseline_authority': baseline.authority}


def evaluate_provider(xyz, geometry, parameters, provider: Callable, field_native=None, cell=None):
    """Provider(xyz, cell) must differentiate all baseline coordinate/cell paths."""
    return evaluate(xyz, geometry, parameters, provider(xyz, cell), field_native, baseline_mode='provider')


def field_au_to_native(field_au, units):
    """Exact historical energy-consistent conversion; retain literal constants."""
    return field_au*units['hartree_kcal_mol']/(units['electric']*units['bohr_A'])


def response_to_source_au(mu_eA, chi_A3, units):
    return {'mu_e_bohr': mu_eA/units['bohr_A'],
            'chi_source_au': chi_A3*units['hartree_kcal_mol']/(units['electric']*units['bohr_A']**2),
            'chi_geometric_bohr3': chi_A3/units['bohr_A']**3}
