"""Sparse lower-triangle real operator; no fixture or native execution API."""
from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree

PARENT = Path(__file__).resolve().parent.parent / 'real_bound_feasibility_v1/BOUND.py'
PARENT_SHA = '2f6c5474dcc82308a65b8a4ae96d2d401630fcd783eb64cd376d4893f1386929'


def parent():
    if hashlib.sha256(PARENT.read_bytes()).hexdigest() != PARENT_SHA:
        raise ValueError('real-tensor source pin mismatch')
    spec = importlib.util.spec_from_file_location('pinned_real_tensor', PARENT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def assemble(xyz, types, alpha, thole, sides, *, cutoff, ewald,
             pair_overrides=None, stiffness=None, shift=0.,
             max_assembly_bytes, batch_pairs=4096, mutual_scales_all_one=True):
    """Return lower CSC triangle and an unshifted site-block upper bound.

    stiffness maps zero-based active site to normalized isotropic beta=exp(sigma)-1.
    All sites here are active positive-alpha unknowns. No hidden zero-alpha
    filtering, covalent scaling, PME term, or SPD admission is supplied.
    """
    xyz = np.asarray(xyz, dtype=np.float64)
    raw_types = np.asarray(types)
    alpha, thole, sides = (np.asarray(x, dtype=np.float64)
                           for x in (alpha, thole, sides))
    n = len(xyz)
    if (xyz.shape != (n, 3) or n == 0 or raw_types.shape != (n,)
            or not np.issubdtype(raw_types.dtype, np.integer)
            or alpha.shape != (n,) or thole.shape != (n,) or sides.shape != (3,)):
        raise ValueError('invalid active-site arrays')
    if not all(np.isfinite(x).all() for x in (xyz, alpha, thole, sides)):
        raise ValueError('nonfinite input')
    if (np.any(alpha <= 0) or np.any(thole < 0) or np.any(sides <= 0)
            or not np.isfinite([cutoff, ewald, shift]).all()
            or not 0 < cutoff < min(sides)/2 or ewald < 0 or shift < 0
            or mutual_scales_all_one is not True):
        raise ValueError('unsupported real-operator domain')
    if (not isinstance(batch_pairs, int) or not 1 <= batch_pairs <= 65536
            or not isinstance(max_assembly_bytes, int) or max_assembly_bytes <= 0):
        raise ValueError('invalid assembly budget')
    overrides = {}
    for key, gamma in (pair_overrides or {}).items():
        if (len(key) != 2 or any(not isinstance(t, (int, np.integer)) for t in key)
                or not np.isfinite(gamma) or gamma < 0):
            raise ValueError('invalid pair override')
        canonical = tuple(sorted(int(t) for t in key))
        if canonical in overrides:
            raise ValueError('duplicate symmetric pair override')
        overrides[canonical] = float(gamma)
    beta = np.zeros(n)
    for site, value in (stiffness or {}).items():
        if (not isinstance(site, (int, np.integer)) or not 0 <= site < n
                or not np.isfinite(value) or value <= -1):
            raise ValueError('invalid normalized isotropic stiffness')
        beta[site] = value
    # The pair-query allocation precedes this exact-count budget; no factor fill
    # is included. Scope is documented so this is not a process-memory promise.
    pairs = cKDTree(xyz % sides, boxsize=sides).query_pairs(cutoff, output_type='ndarray')
    npairs = len(pairs)
    entries = 9*npairs + 3*n
    if 3*n >= np.iinfo(np.int32).max or entries >= np.iinfo(np.int32).max:
        raise ValueError('outside declared int32 sparse index domain')
    # COO arrays, lower CSC and conversion/sorting workspace, pair indices and
    # a conservative bounded tensor batch; actual backend fill is separate.
    estimated_bytes = 3*entries*(8+4+4) + pairs.nbytes + 256*min(batch_pairs, npairs) + 128*n
    if estimated_bytes > max_assembly_bytes:
        raise MemoryError('declared sparse assembly storage ceiling exceeded')
    rows, cols = np.empty(entries, np.int32), np.empty(entries, np.int32)
    data = np.empty(entries, np.float64)
    tself = 4*ewald**3/(3*np.sqrt(np.pi))
    diagonal = 1-alpha*tself+beta
    d = np.arange(3*n, dtype=np.int32)
    rows[:3*n], cols[:3*n], data[:3*n] = d, d, np.repeat(diagonal-shift, 3)
    row_norm = np.zeros(n)
    original = parent()
    ai = np.repeat(np.arange(3), 3)
    bj = np.tile(np.arange(3), 3)
    offset = 3*n
    for start in range(0, npairs, batch_pairs):
        chunk = pairs[start:start+batch_pairs]
        i, j = chunk[:, 0], chunk[:, 1]
        if np.any(i >= j):
            raise ValueError('pair enumeration is not unique upper site pairs')
        disp = xyz[j]-xyz[i]
        disp -= sides*np.rint(disp/sides)
        gamma = np.minimum(thole[i], thole[j])
        gamma = np.where(gamma == 0, np.maximum(thole[i], thole[j]), gamma)
        for (ta, tb), value in overrides.items():
            mask = ((raw_types[i] == ta) & (raw_types[j] == tb)) | ((raw_types[i] == tb) & (raw_types[j] == ta))
            gamma[mask] = value
        tensors = original.tensor(disp, alpha[i], alpha[j], gamma, ewald)
        count = 9*len(chunk)
        # Upper tensor T[a,b] maps to lower block entry (j,b),(i,a).
        rows[offset:offset+count] = (3*j[:, None]+bj).ravel()
        cols[offset:offset+count] = (3*i[:, None]+ai).ravel()
        data[offset:offset+count] = -tensors.reshape(-1)
        offset += count
        norms = np.max(np.abs(np.linalg.eigvalsh(tensors)), axis=1)
        np.add.at(row_norm, i, norms)
        np.add.at(row_norm, j, norms)
    if offset != entries or not np.isfinite(data).all():
        raise ValueError('incomplete or nonfinite assembly')
    lower = coo_matrix((data, (rows, cols)), shape=(3*n, 3*n)).tocsc()
    lower.sort_indices()
    if lower.nnz != entries or not lower.has_canonical_format:
        raise ValueError('duplicate or noncanonical sparse entries')
    info = {'sites': n, 'dimension': 3*n, 'pairs': npairs, 'lower_entries': entries,
            'triangle': 'lower', 'normalization': 'sqrt(alpha) D sqrt(alpha)',
            'shift': float(shift), 'unshifted_site_upper': float(np.max(diagonal+row_norm)),
            'unshifted_site_lower': float(np.min(diagonal-row_norm)),
            'estimated_assembly_bytes': int(estimated_bytes),
            'native_pair_source_sha256': PARENT_SHA,
            'mutual_scales_all_one': True, 'PME_included': False,
            'stability_admitted': False, 'native_operator_parity_admitted': False}
    return lower, info
