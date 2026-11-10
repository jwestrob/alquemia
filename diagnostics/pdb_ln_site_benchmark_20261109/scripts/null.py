"""Composition-matched null (PLAN.md section 7).
Usage: python -I null.py WS SCRIPTS AUDIT_CIF NPROC  -> WS/out/null.jsonl"""
import json
import os
import random
import sys
from multiprocessing import Pool

import numpy as np
from scipy.spatial import cKDTree

WS, SCR, ACIF, NPROC = sys.argv[1:5]
sys.path.insert(0, SCR)
import common as C  # noqa: E402
import detect as Dt  # noqa: E402

C.CIF_DIRS[:] = [os.path.join(WS, 'cif'), ACIF]
NPERM = 200


def stats(P, sites):
    pts = Dt.candidate_points(P)
    S, _, _ = Dt.v0_scores(P, pts)
    mx = float(S.max()) if len(S) else 0.0
    if len(sites) and len(pts):
        d, _ = cKDTree(np.array(sites)).query(pts)
        m = d <= 4.0
        loc = float(S[m].max()) if m.any() else 0.0
    else:
        loc = 0.0
    return mx, loc


def work(u):
    st = C.read(u['entry'])
    res = C.chain_residues(st[0], u['chain'])
    P = Dt.prepare(res)
    sites = [s['xyz'] for s in u['sites'] if s['n_don'] >= 3]
    obs = stats(P, sites)
    names = P['names'].copy()
    pos = np.where((names != 'GLY') & (names != 'PRO'))[0]
    expo = cKDTree(P['ca']).query_ball_point(P['anchor'][pos], 10.0, return_length=True)
    q = np.quantile(expo, [1 / 3, 2 / 3])
    b = np.digitize(expo, q)
    rng = np.random.default_rng([42, len(res)])
    null = []
    for _ in range(NPERM):
        nn = names.copy()
        for k in range(3):
            idx = pos[b == k]
            nn[idx] = names[idx][rng.permutation(len(idx))]
        P2 = dict(P)
        P2['names'] = nn
        null.append(stats(P2, sites))
    null = np.array(null)
    return {'uid': u['uid'], 'c95': u['c95'], 'c30': u['c30'], 'obs_max': obs[0], 'obs_site': obs[1],
            'pct_max': float((null[:, 0] < obs[0]).mean() + 0.5 * (null[:, 0] == obs[0]).mean()),
            'pct_site': float((null[:, 1] < obs[1]).mean() + 0.5 * (null[:, 1] == obs[1]).mean()),
            'null_max_p95': float(np.percentile(null[:, 0], 95)), 'null_site_p95': float(np.percentile(null[:, 1], 95)),
            'null_max_mean': float(null[:, 0].mean()), 'null_site_mean': float(null[:, 1].mean())}


if __name__ == '__main__':
    U = json.load(open(os.path.join(WS, 'out', 'units.json')))
    part = sorted({u['c95'] for u in U if u['tier'] == 'partner' and any(s['n_don'] >= 3 for s in u['sites'])})
    random.Random(42).shuffle(part)
    chosen = []
    for c in part:
        su = sorted([u for u in U if u['tier'] == 'stripped' and u['c95'] == c and any(s['n_don'] >= 3 for s in u['sites'])],
                    key=lambda u: u['uid'])
        if su:
            chosen.append(su[0])
        if len(chosen) == 50:
            break
    with Pool(int(NPROC)) as pool, open(os.path.join(WS, 'out', 'null.jsonl'), 'w') as fh:
        for r in pool.imap_unordered(work, chosen):
            fh.write(json.dumps(r) + '\n'); fh.flush()
            print(r['uid'], r['pct_max'], r['pct_site'], flush=True)
