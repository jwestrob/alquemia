"""Metrics, frozen reading, capacity curve, diagnostics (PLAN.md sections 4-6).
Usage: python -I stats.py WS OUTDIR"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

WS, OUT = sys.argv[1:3]
units = {u['uid']: u for u in json.load(open(os.path.join(WS, 'out', 'units.json')))}
scores = {}
errs = []
for line in open(os.path.join(WS, 'out', 'scores.jsonl')):
    r = json.loads(line)
    if 'error' in r:
        errs.append(r['uid'])
    else:
        scores[r['uid']] = r
split = dict(l.split() for l in open(os.path.join(WS, 'out', 'splits.tsv')))
M = ('v0', 'b1', 'b2', 'b3')
KS = (1, 3, 5)
recs = []
for uid, sc in scores.items():
    u = units[uid]
    sm = {s['sid']: s for s in u['sites']}
    for s in sc['sites']:
        m = sm[s['sid']]
        recs.append({'tier': u['tier'], 'kind': u['kind'], 'uid': uid, 'entry': u['entry'], 'c95': u['c95'],
                     'c30': u['c30'], 'split': split[u['c30']], 'n_don': m['n_don'], 'n_chain': m['n_chain'],
                     'cls': m['cls'], 'cluster5': m['cluster5'], 'rmsd': m['rmsd'], 'metal': m['metal_at_site'],
                     **{f'r_{k}': s['rank'][k] for k in M}, 'S_true': s['S_true'], 'an_true': s['an_true'],
                     'clash': s['clash_true'], 'near': s['n_near_res'], 'reach': s['n_near_reach'],
                     'near_res': s['near_res'], 'S_top1': s['S_top1'], 'sid': s['sid']})
L = []


def p(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    L.append(s)


def protein_vals(rs, f):
    """site -> uid -> entry -> protein(c95); returns {c95: (c30, value)}."""
    a = defaultdict(list)
    for r in rs:
        a[(r['c95'], r['entry'], r['uid'])].append(f(r))
    b = defaultdict(list)
    for (c, e, u), v in a.items():
        b[(c, e)].append(np.mean(v))
    c_ = defaultdict(list)
    for (c, e), v in b.items():
        c_[c].append(np.mean(v))
    c30 = {r['c95']: r['c30'] for r in rs}
    return {c: (c30[c], float(np.mean(v))) for c, v in c_.items()}


RNG_RESAMPLES = None


def boot(pv_list):
    """pv_list: list of {c95:(c30,val)} with identical keys; returns point and CI per entry + diff 0-1."""
    keys = sorted(pv_list[0])
    cl = sorted({pv_list[0][k][0] for k in keys})
    ci = {c: i for i, c in enumerate(cl)}
    idx = np.array([ci[pv_list[0][k][0]] for k in keys])
    V = np.array([[pv[k][1] for k in keys] for pv in pv_list])
    sums = np.zeros((len(pv_list), len(cl)))
    for j in range(len(pv_list)):
        np.add.at(sums[j], idx, V[j])
    cnt = np.bincount(idx, minlength=len(cl)).astype(float)
    rng = np.random.default_rng(42)
    draws = rng.integers(0, len(cl), size=(2000, len(cl)))
    W = np.stack([np.bincount(d, minlength=len(cl)) for d in draws])  # 2000 x ncl
    est = (W @ sums.T) / (W @ cnt)[:, None]
    point = V.mean(1)
    out = [(point[j], np.percentile(est[:, j], 2.5), np.percentile(est[:, j], 97.5)) for j in range(len(pv_list))]
    diff = est[:, 0] - est[:, 1] if len(pv_list) > 1 else None
    return out, diff, len(keys), len(cl)


def recall_table(rs, label):
    res = {}
    if not rs:
        p(f'{label}: no sites'); return res
    for k in KS:
        pvs = [protein_vals(rs, lambda r, m=m: float(r[f'r_{m}'] is not None and r[f'r_{m}'] < k)) for m in M]
        out, diff, npro, ncl = boot([pvs[0], pvs[2], pvs[1], pvs[3]])
        d_pt = out[0][0] - out[1][0]
        dl, dh = np.percentile(diff, 2.5), np.percentile(diff, 97.5)
        res[k] = {'v0': out[0], 'b2': out[1], 'b1': out[2], 'b3': out[3], 'diff': (d_pt, dl, dh)}
        p(f'{label} k={k} sites={len(rs)} proteins={npro} c30={ncl} | ' +
          ' '.join(f'{n}={o[0]:.3f}[{o[1]:.3f},{o[2]:.3f}]' for n, o in zip(('v0', 'b2', 'b1', 'b3'), out)) +
          f' | v0-b2={d_pt:+.3f}[{dl:+.3f},{dh:+.3f}]')
    return res


p('scored units', len(scores), 'errors', len(errs), errs[:5])
p('B1 errors', sum(1 for s in scores.values() if s['b1_err']), 'B1 units with <5 centres',
  sum(1 for s in scores.values() if s['n_cent']['b1'] < 5), 'of', len(scores))
p('v0 centres per unit median', np.median([s['n_cent']['v0'] for s in scores.values()]),
  'b1', np.median([s['n_cent']['b1'] for s in scores.values()]),
  'unit sec total', round(sum(s['sec'] for s in scores.values())))

# ---------------- benchmark size
for t in ('stripped', 'partner'):
    rs = [r for r in recs if r['tier'] == t]
    p(f'\n== {t}: sites {len(rs)} units {len({r["uid"] for r in rs})} entries {len({r["entry"] for r in rs})} '
      f'c95 {len({r["c95"] for r in rs})} c30 {len({r["c30"] for r in rs})}')
    p(' donors', sorted(Counter(min(r['n_don'], 8) for r in rs).items()))
    p(' class', dict(Counter(r['cls'] for r in rs)), 'class|>=3', dict(Counter(r['cls'] for r in rs if r['n_don'] >= 3)),
      'cluster5', sum(r['cluster5'] for r in rs))
    if t == 'partner':
        for k in ('ca', 'metal_free', 'other_metal'):
            q = [r for r in rs if r['kind'] == k]
            p(f' kind {k}: sites {len(q)} entries {len({r["entry"] for r in q})} c95 {len({r["c95"] for r in q})} '
              f'c30 {len({r["c30"] for r in q})} >=3 {sum(r["n_don"] >= 3 for r in q)}')
        p(' metal at transferred site', dict(Counter(r['metal'] for r in rs)),
          'rmsd median', np.median([r['rmsd'] for r in rs]), '>2A', sum(r['rmsd'] > 2 for r in rs))

# ---------------- recall
allres = {}
for t, kinds in (('partner', (None, 'ca', 'metal_free', 'other_metal')), ('stripped', (None,))):
    for kd in kinds:
        for mind in (3, 2):
            for cls in ('all', 'single-chain'):
                for sp in ('all', 'test'):
                    rs = [r for r in recs if r['tier'] == t and (kd is None or r['kind'] == kd) and r['n_don'] >= mind
                          and (cls == 'all' or r['cls'] == cls) and (sp == 'all' or r['split'] == 'test')]
                    lab = f"{t}{'/' + kd if kd else '/pooled' if t == 'partner' else ''} >={mind} {cls} split={sp}"
                    if sp == 'test' and (kd is not None or mind == 2):
                        continue
                    allres[lab] = recall_table(rs, lab)

# ---------------- frozen reading
p('\n== FROZEN READING (partner pooled, >=3 donors, k=5)')
for cls in ('all', 'single-chain'):
    r = allres[f'partner/pooled >=3 {cls} split=all'][5]
    v, d = r['v0'], r['diff']
    useful = v[0] >= 0.5 and d[0] >= 0.10 and d[1] > 0
    p(f' [{cls}] v0 recall {v[0]:.3f} [{v[1]:.3f},{v[2]:.3f}]; v0-B2 {d[0]:+.3f} [{d[1]:+.3f},{d[2]:+.3f}] -> '
      + ('useful' if useful else 'not better than counting acidic residues')
      + ('' if cls == 'all' else '  (co-reported, lead addition)'))

# ---------------- capacity curve
p('\n== capacity curve (stripped tier, all v0 centres; x = Asp/Glu of ONE chain that can reach the centre)')
bins = list(range(8)) + ['>=8']
for mode in ('any Ln', 'single-chain site'):
    per = defaultdict(lambda: np.zeros((9, 2)))
    for uid, sc in scores.items():
        if units[uid]['tier'] != 'stripped':
            continue
        c30 = units[uid]['c30']
        for S, an, anyln, scp in sc['cap']:
            b = min(an, 8)
            if mode == 'any Ln':
                pos = anyln
            else:
                if anyln and not scp:
                    continue
                pos = scp
            per[c30][b] += (1, pos)
    cl = sorted(per)
    A = np.stack([per[c] for c in cl])  # ncl x 9 x 2
    tot = A.sum(0)
    rng = np.random.default_rng(42)
    draws = rng.integers(0, len(cl), size=(2000, len(cl)))
    W = np.stack([np.bincount(d, minlength=len(cl)) for d in draws])
    est = np.einsum('dc,cbk->dbk', W, A)
    frac = est[:, :, 1] / np.maximum(est[:, :, 0], 1)
    p(f' [{mode}] c30={len(cl)}')
    for i, b in enumerate(bins):
        p(f'   n_anion={b}: centres {int(tot[i, 0])} with Ln {int(tot[i, 1])} fraction {tot[i, 1] / max(tot[i, 0], 1):.4f} '
          f'[{np.percentile(frac[:, i], 2.5):.4f},{np.percentile(frac[:, i], 97.5):.4f}]')

# ---------------- diagnostics
p('\n== diagnostics at true site positions')
for t in ('stripped', 'partner'):
    rs = [r for r in recs if r['tier'] == t and r['n_don'] >= 3]
    p(f' {t} >=3: clash at true position {sum(r["clash"] for r in rs)}/{len(rs)}; '
      f'same-chain O/N residues reachable {sum(r["reach"] for r in rs)}/{sum(r["near"] for r in rs)}; '
      f'S_true median {np.median([r["S_true"] for r in rs]):.2f}, S_top1 median {np.median([r["S_top1"] or 0 for r in rs]):.2f}; '
      f'S_true>=S_top1 {sum((r["S_true"] >= (r["S_top1"] or 0)) for r in rs)}')
    for cls in ('single-chain', 'symmetry-bridging', 'inter-chain'):
        q = [r for r in rs if r['cls'] == cls]
        if q:
            p(f'   {cls}: n={len(q)} v0 hit@5 raw {np.mean([r["r_v0"] is not None and r["r_v0"] < 5 for r in q]):.3f} '
              f'b2 {np.mean([r["r_b2"] is not None and r["r_b2"] < 5 for r in q]):.3f} '
              f'median same-chain donors {np.median([r["n_chain"] for r in q])}')
miss = [r for r in recs if r['tier'] == 'partner' and r['n_don'] >= 3 and not (r['r_v0'] is not None and r['r_v0'] < 5)]
p(f'\n== partner >=3 misses at k=5: {len(miss)} sites; 10 sampled (seed 42)')
sel = np.random.default_rng(42).choice(len(miss), size=min(10, len(miss)), replace=False)
for i in sel:
    r = miss[i]
    p(f" {r['uid']} site={r['sid']} kind={r['kind']} cls={r['cls']} n_don={r['n_don']} n_chain={r['n_chain']} "
      f"rmsd={r['rmsd']} metal={r['metal']} clash={r['clash']} S_true={r['S_true']} an={r['an_true']} S_top1={r['S_top1']} "
      f"rank_v0={r['r_v0']} rank_b2={r['r_b2']} near={r['near_res']}")
open(os.path.join(OUT, 'stats.log'), 'w').write('\n'.join(L) + '\n')
with open(os.path.join(OUT, 'site_ranks.tsv'), 'w') as fh:
    cols = ['tier', 'kind', 'uid', 'sid', 'c95', 'c30', 'split', 'n_don', 'n_chain', 'cls', 'cluster5', 'rmsd', 'metal',
            'r_v0', 'r_b1', 'r_b2', 'r_b3', 'S_true', 'an_true', 'clash', 'near', 'reach', 'S_top1']
    fh.write('\t'.join(cols) + '\n')
    for r in recs:
        fh.write('\t'.join('' if r[c] is None else str(r[c]) for c in cols) + '\n')
