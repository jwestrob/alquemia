"""Detector v0 (reachable donor capacity) and baselines B2/B3 on one chain. Constants frozen in PLAN.md."""
import zlib

import numpy as np
from scipy.spatial import cKDTree

BAND = (2.2, 2.9)
REACH = {'SER': [(1.35, 1.55, 0.5)], 'THR': [(1.35, 1.55, 0.5)],
         'ASP': [(2.25, 2.55, 1.0)], 'ASN': [(2.25, 2.55, 0.5)],
         'GLU': [(2.70, 3.90, 1.0)], 'GLN': [(2.70, 3.90, 0.5)],
         'HIS': [(2.40, 2.60, 0.25), (3.50, 3.75, 0.25)], 'TYR': [(5.55, 5.85, 0.5)]}
BB_W = 0.5
MAX_ANGLE = 120.0
GRID_R = 9.0
NMS_R = 4.0
CLASH = {'bbNCAC': 3.0, 'bbO': 2.2, 'CB': 2.8}


def virtual_cb(n, ca, c):
    b = ca - n
    cc = c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def prepare(residues):
    """residues: list of (resname, num, icode, atoms) -> arrays."""
    names, anchor, ca, bbO, bbC, bbNCAC, heavy = [], [], [], [], [], [], []
    for (rn, num, ic, at) in residues:
        names.append(rn)
        ca.append(at['CA'])
        if 'CB' in at:
            cb = at['CB']
        elif rn != 'GLY' and 'N' in at and 'C' in at:
            cb = virtual_cb(at['N'], at['CA'], at['C'])
        else:
            cb = at['CA'] if rn == 'GLY' or 'CB' not in at else at['CB']
        anchor.append(cb)
        bbO.append(at.get('O', np.full(3, np.nan)))
        bbC.append(at.get('C', np.full(3, np.nan)))
        for k in ('N', 'CA', 'C'):
            if k in at:
                bbNCAC.append(at[k])
        heavy.extend(at.values())
    return {'names': np.array(names), 'anchor': np.array(anchor), 'ca': np.array(ca), 'bbO': np.array(bbO),
            'bbC': np.array(bbC), 'bbNCAC': np.array(bbNCAC), 'heavy': np.array(heavy),
            'is_gly': np.array([n == 'GLY' for n in names])}


def lattice_near(centres, r):
    lo = np.floor(centres.min(0) - r)
    hi = np.ceil(centres.max(0) + r)
    axes = [np.arange(lo[i], hi[i] + 1) for i in range(3)]
    g = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
    t = cKDTree(centres)
    d, _ = t.query(g, distance_upper_bound=r)
    return g[np.isfinite(d)]


def candidate_points(P):
    don = np.array([n in REACH for n in P['names']])
    okO = np.isfinite(P['bbO'][:, 0])
    g1 = lattice_near(P['anchor'][don & ~P['is_gly']], GRID_R) if (don & ~P['is_gly']).any() else np.zeros((0, 3))
    g2 = lattice_near(P['bbO'][okO], BAND[1])
    g = np.unique(np.concatenate([g1, g2]), axis=0)
    keep = np.ones(len(g), bool)
    for arr, r in ((P['bbNCAC'], CLASH['bbNCAC']), (P['bbO'][okO], CLASH['bbO']),
                   (P['anchor'][~P['is_gly']], CLASH['CB'])):
        if len(arr):
            d, _ = cKDTree(arr).query(g, distance_upper_bound=r)
            keep &= ~np.isfinite(d)
    return g[keep]


def _pairs(tree, pts, r):
    lst = tree.query_ball_point(pts, r)
    lens = np.fromiter((len(x) for x in lst), int, len(lst))
    pi = np.repeat(np.arange(len(pts)), lens)
    ri = np.fromiter((j for x in lst for j in x), int, lens.sum())
    return pi, ri


def v0_scores(P, pts):
    """Return S, n_anion, tiebreak for every point."""
    n_pts, nres = len(pts), len(P['names'])
    don_idx = np.array([i for i, n in enumerate(P['names']) if n in REACH and not P['is_gly'][i]], int)
    pi_all, ri_all, w_all = [], [], []
    S_an = np.zeros(n_pts)
    if len(don_idx):
        tree = cKDTree(P['anchor'][don_idx])
        pi, rj = _pairs(tree, pts, max(r[1] for v in REACH.values() for r in v) + BAND[1])
        ri = don_idx[rj]
        cb = P['anchor'][ri]
        v = pts[pi] - cb
        D = np.linalg.norm(v, axis=1)
        ax = cb - P['ca'][ri]
        cosang = np.einsum('ij,ij->i', ax, v) / (np.linalg.norm(ax, axis=1) * np.maximum(D, 1e-9))
        ang_ok = cosang >= np.cos(np.radians(MAX_ANGLE))
        w = np.zeros(len(pi))
        names = P['names'][ri]
        for rn, shells in REACH.items():
            m = names == rn
            if not m.any():
                continue
            for (rmin, rmax, wt) in shells:
                lo = max(0.0, rmin - BAND[1], BAND[0] - rmax)
                ok = m & ang_ok & (D >= lo) & (D <= rmax + BAND[1])
                w = np.where(ok, np.maximum(w, wt), w)
        sel = w > 0
        pi_all.append(pi[sel]); ri_all.append(ri[sel]); w_all.append(w[sel])
        an = sel & ((names == 'ASP') | (names == 'GLU'))
        S_an = np.bincount(pi[an], minlength=n_pts).astype(float)
    okO = np.where(np.isfinite(P['bbO'][:, 0]))[0]
    tO = cKDTree(P['bbO'][okO])
    pi, rj = _pairs(tO, pts, BAND[1])
    ri = okO[rj]
    O, Cc = P['bbO'][ri], P['bbC'][ri]
    d = np.linalg.norm(pts[pi] - O, axis=1)
    u1, u2 = Cc - O, pts[pi] - O
    cosang = np.einsum('ij,ij->i', u1, u2) / (np.linalg.norm(u1, axis=1) * np.maximum(np.linalg.norm(u2, axis=1), 1e-9))
    sel = (d >= BAND[0]) & (cosang <= np.cos(np.radians(120.0)))
    pi_all.append(pi[sel]); ri_all.append(ri[sel]); w_all.append(np.full(sel.sum(), BB_W))
    pi = np.concatenate(pi_all); ri = np.concatenate(ri_all); w = np.concatenate(w_all)
    S = np.zeros(n_pts); tb = np.full(n_pts, np.inf)
    if len(pi):
        key = pi.astype(np.int64) * nres + ri
        order = np.lexsort((-w, key))
        key, w2, pi2, ri2 = key[order], w[order], pi[order], ri[order]
        first = np.ones(len(key), bool); first[1:] = key[1:] != key[:-1]
        pi2, ri2, w2 = pi2[first], ri2[first], w2[first]
        S = np.bincount(pi2, weights=w2, minlength=n_pts)
        dist = np.linalg.norm(pts[pi2] - P['anchor'][ri2], axis=1)
        cnt = np.bincount(pi2, minlength=n_pts)
        tb = np.where(cnt > 0, np.bincount(pi2, weights=dist, minlength=n_pts) / np.maximum(cnt, 1), np.inf)
    return S, S_an, tb


def b2_scores(P, pts):
    ac = np.array([i for i, n in enumerate(P['names']) if n in ('ASP', 'GLU')], int)
    n = np.zeros(len(pts)); tb = np.full(len(pts), np.inf)
    if len(ac):
        pi, rj = _pairs(cKDTree(P['anchor'][ac]), pts, 8.0)
        n = np.bincount(pi, minlength=len(pts)).astype(float)
        dist = np.linalg.norm(pts[pi] - P['anchor'][ac[rj]], axis=1)
        tb = np.where(n > 0, np.bincount(pi, weights=dist, minlength=len(pts)) / np.maximum(n, 1), np.inf)
    return n, tb


def nms(pts, order, r=NMS_R):
    tree = cKDTree(pts)
    sup = np.zeros(len(pts), bool)
    kept = []
    for i in order:
        if sup[i]:
            continue
        kept.append(i)
        sup[tree.query_ball_point(pts[i], r)] = True
    return np.array(kept, int)


def rank(score, tb, pts):
    m = np.where(score > 0)[0]
    order = m[np.lexsort((tb[m], -score[m]))]
    return nms(pts, order)


def b3_points(P, unit_id):
    h = P['heavy']
    g = lattice_near(h, 5.0)
    t = cKDTree(h)
    d, _ = t.query(g, distance_upper_bound=5.0)
    g = g[np.isfinite(d) & (d >= 3.0)]
    rng = np.random.default_rng([42, zlib.crc32(unit_id.encode())])
    order = rng.permutation(len(g))
    return g, nms(g, order)
