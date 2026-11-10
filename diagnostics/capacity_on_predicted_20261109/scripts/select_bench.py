"""Part 2 sample: benchmark proteins for the predicted tier (PLAN.md section 2).
Usage: python -I select_bench.py BENCH_WS BENCH_SCRIPTS AUDIT_CIF OUT_DIR
Writes OUT_DIR/bench_sample.tsv, OUT_DIR/bench_sample.faa, OUT_DIR/bench_eligible.tsv"""
import json
import os
import random
import sys
from collections import defaultdict

import gemmi

BWS, BSCR, ACIF, OUT = sys.argv[1:5]
sys.path.insert(0, BSCR)
import common as C  # noqa: E402

C.CIF_DIRS[:] = [os.path.join(BWS, 'cif'), ACIF]
MAXLEN = 1000


def entity_seq(st, chain):
    ent = None
    for en in st.entities:
        if en.entity_type == gemmi.EntityType.Polymer and any(
                sc in en.subchains for sc in {r.subchain for r in st[0][chain]}):
            ent = en
            break
    if ent is None:
        return None
    s = []
    for name in ent.full_sequence:
        t = gemmi.find_tabulated_residue(gemmi.Entity.first_mon(name))
        s.append(t.one_letter_code.upper() if t and t.one_letter_code.strip() else 'X')
    return ''.join(s)


units = json.load(open(os.path.join(BWS, 'out', 'units.json')))
elig = []
for u in units:
    if u['tier'] != 'stripped':
        continue
    sites = [s for s in u['sites'] if s['cls'] == 'single-chain' and s['n_don'] >= 3]
    if not sites:
        continue
    st = C.read(u['entry'])
    seq = entity_seq(st, u['chain'])
    if seq is None or len(seq) > MAXLEN:
        continue
    elig.append({'uid': u['uid'], 'entry': u['entry'], 'chain': u['chain'], 'c95': u['c95'], 'c30': u['c30'],
                 'n_sites': len(sites), 'len': len(seq), 'n_X': seq.count('X'), 'seq': seq})
by = defaultdict(list)
for e in elig:
    by[e['c30']].append(e)
cl = sorted(by)
rng = random.Random(42)
rng.shuffle(cl)
pick = []
for c in cl[:50]:
    cand = sorted(by[c], key=lambda e: e['uid'])
    pick.append(rng.choice(cand))
with open(os.path.join(OUT, 'bench_eligible.tsv'), 'w') as fh:
    fh.write('uid\tc95\tc30\tn_sites\tlen\n')
    for e in elig:
        fh.write(f"{e['uid']}\t{e['c95']}\t{e['c30']}\t{e['n_sites']}\t{e['len']}\n")
with open(os.path.join(OUT, 'bench_sample.tsv'), 'w') as fh, open(os.path.join(OUT, 'bench_sample.faa'), 'w') as fa:
    fh.write('name\tuid\tentry\tchain\tc95\tc30\tn_sites\tlen\tn_X\n')
    for e in pick:
        name = f"B_{e['entry']}_{e['chain']}"
        fh.write(f"{name}\t{e['uid']}\t{e['entry']}\t{e['chain']}\t{e['c95']}\t{e['c30']}\t{e['n_sites']}\t{e['len']}\t{e['n_X']}\n")
        fa.write(f">{name}\n{e['seq']}\n")
print('eligible units', len(elig), 'c95', len({e['c95'] for e in elig}), 'c30', len(cl), 'picked', len(pick),
      'sites', sum(e['n_sites'] for e in pick), 'aa', sum(e['len'] for e in pick), 'with X', sum(e['n_X'] > 0 for e in pick))
