"""Step 3: frozen paired statistics (Part 2) and benchmark power (Part 3).
Usage: python -I stats.py WS OUTDIR"""
import collections
import csv
import json
import math
import os
import sys

import numpy as np

WS, OUT = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(os.path.join(WS, 'out', 'pairs.tsv')), delimiter='\t'))
for r in rows:
    r['protein'] = r['uniprot'] or r['cluster']
    r['native_ca'] = int(r['ca_site_res_CaBinding']) >= 2
    r['excluded'] = bool(r['exclusion'])
    r['reverse'] = int(r['ln_site_res_LnBinding']) >= 1
    for k in ('ln_cn', 'ca_cn', 'ln_chel', 'ca_chel', 'ln_rule_La', 'ca_rule_La', 'ca_protein_donors'):
        r[k] = float(r[k])
    for k in ('ln_mean_d', 'ca_mean_d'):
        r[k] = float(r[k]) if r[k] not in ('', 'None') else float('nan')

STRATA = {
    'PRIMARY native-Ca, >=4 protein donors (Ca deposit)': lambda r: r['native_ca'] and not r['excluded'] and not r['reverse'] and r['ca_protein_donors'] >= 4,
    'native-Ca, any donor count': lambda r: r['native_ca'] and not r['excluded'] and not r['reverse'],
    'native-Ca primary, light Ln (La-Eu)': lambda r: r['native_ca'] and not r['excluded'] and not r['reverse'] and r['ca_protein_donors'] >= 4 and r['ln_light'] == '1',
    'native-Ca primary, heavy Ln (Gd-Lu)': lambda r: r['native_ca'] and not r['excluded'] and not r['reverse'] and r['ca_protein_donors'] >= 4 and r['ln_light'] == '0',
    'UNVERIFIED, >=4 protein donors': lambda r: not r['native_ca'] and not r['excluded'] and not r['reverse'] and r['ca_protein_donors'] >= 4,
    'UNVERIFIED, any donor count': lambda r: not r['native_ca'] and not r['excluded'] and not r['reverse'],
    'excluded by frozen keyword rule (sensitivity)': lambda r: r['excluded'],
}


def per_protein(rr):
    """protein -> dict(p, q, dCN, dMO, dChel, n_ln_units, n_ca_units)."""
    out = {}
    byp = collections.defaultdict(list)
    for r in rr:
        byp[r['protein']].append(r)
    for prot, pr in byp.items():
        ln_units = collections.defaultdict(dict)     # (ln_pdb, ln_key) -> {ln_site: call}
        ca_units = collections.defaultdict(dict)
        deltas = collections.defaultdict(lambda: collections.defaultdict(list))  # ln instance -> partner feats
        for r in pr:
            ln_units[(r['ln_pdb'], r['ln_key'])][r['ln_site']] = r['ln_rule_La']
            ca_units[(r['ca_pdb'], r['ca_key'])][r['ca_site']] = r['ca_rule_La']
            deltas[(r['ln_pdb'], r['ln_key'], r['ln_site'])]['x'].append(r)
        p = np.mean([np.mean(list(v.values())) for v in ln_units.values()])
        q = np.mean([np.mean(list(v.values())) for v in ca_units.values()])
        unit_d = collections.defaultdict(list)
        for (lp, lk, ls), d in deltas.items():
            xs = d['x']
            unit_d[(lp, lk)].append((xs[0]['ln_cn'] - np.mean([x['ca_cn'] for x in xs]),
                                     xs[0]['ln_mean_d'] - np.nanmean([x['ca_mean_d'] for x in xs]),
                                     xs[0]['ln_chel'] - np.mean([x['ca_chel'] for x in xs])))
        dd = np.array([np.mean(v, axis=0) for v in unit_d.values()])
        out[prot] = {'p': float(p), 'q': float(q), 'dCN': float(dd[:, 0].mean()), 'dMO': float(np.nanmean(dd[:, 1])),
                     'dChel': float(dd[:, 2].mean()), 'n_ln_units': len(ln_units), 'n_ca_units': len(ca_units),
                     'n_pairs': len(pr), 'metals': sorted({r['ln_metal'] for r in pr}),
                     'ln_entries': sorted({r['ln_pdb'] for r in pr}), 'ca_entries': sorted({r['ca_pdb'] for r in pr}),
                     'ln_cn_mean': float(np.mean([r['ln_cn'] for r in pr])), 'ca_cn_mean': float(np.mean([r['ca_cn'] for r in pr]))}
    return out


def boot(pp, key_fn, B=2000, seed=42):
    keys = sorted(pp)
    rng = np.random.default_rng(seed)
    vals = np.array([key_fn(pp[k]) for k in keys], float)
    est = float(np.nanmean(vals))
    bs = [np.nanmean(vals[rng.integers(0, len(keys), len(keys))]) for _ in range(B)]
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return est, float(lo), float(hi)


summary = {'n_pairs_total': len(rows), 'strata': {}}
for name, sel in STRATA.items():
    rr = [r for r in rows if sel(r)]
    if not rr:
        summary['strata'][name] = {'n_proteins': 0}
        continue
    pp = per_protein(rr)
    p = boot(pp, lambda d: d['p'])
    q = boot(pp, lambda d: d['q'])
    pq = boot(pp, lambda d: d['p'] - d['q'])
    s = {'n_proteins': len(pp), 'n_pairs': len(rr), 'n_ln_entries': len({r['ln_pdb'] for r in rr}),
         'n_ca_entries': len({r['ca_pdb'] for r in rr}),
         'p': p, 'q': q, 'p_minus_q': pq,
         'dCN': boot(pp, lambda d: d['dCN']), 'dMeanMO_A': boot(pp, lambda d: d['dMO']), 'dChel': boot(pp, lambda d: d['dChel']),
         'pooled_unweighted': {'p': float(np.mean([r['ln_rule_La'] for r in {(r['ln_pdb'], r['ln_site']): r for r in rr}.values()])),
                               'q': float(np.mean([r['ca_rule_La'] for r in {(r['ca_pdb'], r['ca_site']): r for r in rr}.values()]))},
         'proteins': pp}
    if name.startswith('PRIMARY'):
        n, d, (e, lo, hi) = len(pp), pq[0], pq
        if n >= 10 and d >= 0.40 and lo > 0:
            verdict = 'IMPRINT SUPPORTED'
        elif n >= 10 and d <= 0.10:
            verdict = 'IMPRINT REFUTED'
        else:
            verdict = 'INCONCLUSIVE'
        s['frozen_reading'] = verdict
    summary['strata'][name] = s
summary['reverse_pairs'] = [{k: r[k] for k in ('ln_pdb', 'ln_site', 'ca_pdb', 'ca_site', 'ln_rule_La', 'ca_rule_La')}
                            for r in rows if r['reverse']]


# ---------------- Part 3: power ----------------
def binom_sf(k, n, p):  # P(X >= k)
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def mcnemar_min_k(alpha=0.05):  # one-directional discordance, exact two-sided
    k = 1
    while 2 * 0.5 ** k > alpha:
        k += 1
    return k


def mdd_paired_best_case(n, alpha=0.05, power=0.8):
    """Smallest per-unit excess error rate d detectable when every discordant unit favours one method."""
    k = mcnemar_min_k(alpha)
    if n < k:
        return None
    lo, hi = 0.0, 1.0
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if binom_sf(k, n, m) < power else (lo, m)
    return hi


def fisher_two_sided(a, n1, b, n2):
    tot = a + b
    N = n1 + n2

    def h(x):
        return math.comb(n1, x) * math.comb(n2, tot - x) / math.comb(N, tot)
    p0 = h(a)
    return sum(h(x) for x in range(max(0, tot - n2), min(n1, tot) + 1) if h(x) <= p0 * (1 + 1e-9))


def mdd_unpaired(n, p1, alpha=0.05, power=0.8):
    """Smallest drop from accuracy p1 detectable by Fisher's exact test, n units per arm (exact power)."""
    rej = {}
    for p2 in np.arange(p1, -0.0001, -0.005):
        pw = 0.0
        for a in range(n + 1):
            pa = math.comb(n, a) * p1 ** a * (1 - p1) ** (n - a)
            if pa < 1e-12:
                continue
            for b in range(n + 1):
                pb = math.comb(n, b) * p2 ** b * (1 - p2) ** (n - b)
                if pb < 1e-12:
                    continue
                key = (a, b)
                if key not in rej:
                    rej[key] = fisher_two_sided(a, n, b, n) <= alpha
                if rej[key]:
                    pw += pa * pb
        if pw >= power:
            return round(float(p1 - p2), 3)
    return None


power = {
    'assumptions': 'alpha 0.05 two-sided; 80% power; exact tests; paired = both discriminators scored on the same units '
                   '(exact McNemar, best case: every discordant unit favours one method); unpaired = Fisher exact with '
                   'baseline accuracy fixed at the observed rule value. Balanced-accuracy change = half the arm change.',
    'La_family_unit_n3': {
        'min_two_sided_p_sign_test': 2 * 0.5 ** 3, 'min_two_sided_p_fisher_3v3': fisher_two_sided(3, 3, 0, 3),
        'paired_mdd': mdd_paired_best_case(3), 'note': 'no outcome on 3 units reaches p<=0.05: power 0 for any effect'},
    'La_site_unit_n22_ignoring_clustering': {'paired_mdd': mdd_paired_best_case(22), 'unpaired_mdd_from_21of22': mdd_unpaired(22, 21 / 22)},
    'Ca_site_unit_n125': {'paired_mdd': mdd_paired_best_case(125), 'unpaired_mdd_from_122of125': mdd_unpaired(125, 122 / 125)},
    'Ca_deposit_unit_n34': {'paired_mdd': mdd_paired_best_case(34), 'unpaired_mdd_from_0.976': mdd_unpaired(34, 0.976)},
    'mcnemar_min_discordant_units': mcnemar_min_k(),
}
summary['power'] = power
json.dump(summary, open(os.path.join(OUT, 'stats_summary.json'), 'w'), indent=1, default=float)

# compact per-protein table for the report
with open(os.path.join(OUT, 'per_protein.tsv'), 'w') as fh:
    fh.write('stratum\tprotein\tmetals\tn_ln_entries\tn_ca_entries\tp_ln_called_La\tq_ca_called_La\tln_CN_mean\tca_CN_mean\tdCN\tdMeanMO_A\tdChel\n')
    for name, s in summary['strata'].items():
        for prot, d in sorted(s.get('proteins', {}).items()):
            fh.write(f"{name}\t{prot}\t{','.join(d['metals'])}\t{len(d['ln_entries'])}\t{len(d['ca_entries'])}\t{d['p']:.3f}\t{d['q']:.3f}\t"
                     f"{d['ln_cn_mean']:.2f}\t{d['ca_cn_mean']:.2f}\t{d['dCN']:.2f}\t{d['dMO']:.3f}\t{d['dChel']:.2f}\n")
for name, s in summary['strata'].items():
    if s.get('n_proteins'):
        print(f"{name}: n_prot={s['n_proteins']} pairs={s['n_pairs']} p={s['p'][0]:.3f} q={s['q'][0]:.3f} "
              f"p-q={s['p_minus_q'][0]:.3f} [{s['p_minus_q'][1]:.3f},{s['p_minus_q'][2]:.3f}] "
              f"dCN={s['dCN'][0]:.2f} [{s['dCN'][1]:.2f},{s['dCN'][2]:.2f}] dMO={s['dMeanMO_A'][0]:.3f} "
              f"[{s['dMeanMO_A'][1]:.3f},{s['dMeanMO_A'][2]:.3f}] dChel={s['dChel'][0]:.2f} [{s['dChel'][1]:.2f},{s['dChel'][2]:.2f}] "
              f"pooled p={s['pooled_unweighted']['p']:.3f} q={s['pooled_unweighted']['q']:.3f} {s.get('frozen_reading', '')}")
    else:
        print(name, 'n_prot=0')
print('reverse pairs', len(summary['reverse_pairs']))
print(json.dumps(power, indent=1, default=float))
