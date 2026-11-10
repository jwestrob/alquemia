"""Step 2: frozen first-shell features for every Ca/Ln ion, same-site Ln/Ca pairing, UniProt curation.
Usage (gemmi env): python -I pair_sites.py WS COD_SCRIPTS_DIR
Writes WS/out/sites.tsv, WS/out/pairs.tsv, WS/out/pair_log.json."""
import ast
import itertools
import json
import os
import re
import sys

import numpy as np
import gemmi

WS, COD = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rcsb as r  # noqa: E402
r.set_cache(os.path.join(WS, 'api'))

sys.path.insert(0, COD)
import protein_shells as ps  # noqa: E402  (frozen, imported unchanged)

LN = ('La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu')
LIGHT = {'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu'}
ps.TARGET = ('Ca',) + LN                       # adaptation 1 (PLAN)
ps.cshm.all_shapes = lambda *a, **k: {}         # adaptation 2: shape/arrangement unused here
ps.go.arrangement = lambda *a, **k: {}
ps.BV_R0 = {**{m: ps.BV_R0['La'] for m in LN}, 'Ca': ps.BV_R0['Ca']}

# features() and its constants, extracted verbatim from la_ca_diagnostics.py
src = open(os.path.join(COD, 'la_ca_diagnostics.py')).read()
tree = ast.parse(src)
ns = {'np': np, 'math': __import__('math'), 'PRIMARY': 3.0}
for node in tree.body:
    if isinstance(node, ast.Assign) and any(getattr(t, 'id', '') in ('WATER', 'NEUTRAL') for t in node.targets):
        exec(compile(ast.Module([node], []), 'la_ca', 'exec'), ns)
    if isinstance(node, ast.FunctionDef) and node.name == 'features':
        exec(compile(ast.Module([node], []), 'la_ca', 'exec'), ns)
features = ns['features']


def rule(f):
    return f['cn'] >= 9 or f['n_chel'] >= 2


cand = json.load(open(os.path.join(WS, 'out', 'candidates.json')))
BB = ('N', 'CA', 'C', 'O')
LBL = re.compile(r'^([^:]+):([A-Z0-9]{1,3}?)(-?\d+):')


def load(pdb):
    path = os.path.join(WS, 'cif', pdb + '.cif')
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_hydrogens()
    sites = {s['site']: s for s in ps.process(path, None, None)}
    model = st[0]
    chains = {}
    for ch in model:
        poly = ch.get_polymer()
        if len(poly) < 20:
            continue
        seq = ''.join(gemmi.find_tabulated_residue(x.name).one_letter_code.upper()
                      if gemmi.find_tabulated_residue(x.name) else 'X' for x in poly)
        chains[ch.name] = {'res': [x for x in poly], 'seq': seq,
                           'idx': {x.seqid.num: i for i, x in enumerate(poly) if x.seqid.icode == ' '}}
    ions = []
    for ch in model:
        for res in ch:
            for at in res:
                alt = at.altloc if at.altloc != '\0' else ''
                if at.element.name == 'Ca' and alt in ('', 'A'):
                    ions.append((f"{ch.name}:{res.name}{res.seqid.num}:{at.name}" + (f":{alt}" if alt else ''),
                                 np.array(at.pos.tolist())))
    return {'st': st, 'sites': sites, 'chains': chains, 'ca_ions': ions}


def coord_res(site):
    out = []
    for d in site['donors']:
        if d['d'] > 3.0 or not d['polymer'] or '@sym' in d['label']:
            continue
        m = LBL.match(d['label'])
        if m and m.group(2) in ps.AA:
            out.append((m.group(1), int(m.group(3)), m.group(2)))
    return sorted(set(out))


def align_map(a, b):
    al = gemmi.align_string_sequences(list(a), list(b), [])
    ident = al.match_count / max(1, min(len(a), len(b)))
    mp, i, j = {}, 0, 0
    for op in re.findall(r'(\d+)([MID])', al.cigar_str()):
        n, t = int(op[0]), op[1]
        for _ in range(n):
            if t == 'M':
                mp[i] = j; i += 1; j += 1
            elif t == 'I':
                i += 1
            else:
                j += 1
    return ident, mp


def kabsch(P, Q):  # find R,t minimising |R Q + t - P|
    pc, qc = P.mean(0), Q.mean(0)
    H = (Q - qc).T @ (P - pc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = pc - R @ qc
    rmsd = float(np.sqrt(np.mean(np.sum((Q @ R.T + t - P) ** 2, axis=1))))
    return R, t, rmsd


def bb_atoms(chain, i):
    out = {}
    for k in range(i - 2, i + 3):
        if 0 <= k < len(chain['res']):
            res = chain['res'][k]
            for a in BB:
                at = res.find_atom(a, '*')
                if at:
                    out[(k, a)] = np.array(at.pos.tolist())
    return out


# ---------- UniProt / entity metadata ----------
def entity_meta(pdb):
    q = ('{entry(entry_id:"%s"){polymer_entities{rcsb_polymer_entity_container_identifiers{entity_id auth_asym_ids uniprot_ids} '
         'rcsb_polymer_entity{pdbx_description} rcsb_entity_source_organism{scientific_name} '
         'rcsb_polymer_entity_align{reference_database_name reference_database_accession '
         'aligned_regions{entity_beg_seq_id ref_beg_seq_id length}}}}}') % pdb
    d = r.graphql(q)['data']['entry']
    out = {}
    for pe in d['polymer_entities'] or []:
        ids = pe['rcsb_polymer_entity_container_identifiers']
        aln = [a for a in pe.get('rcsb_polymer_entity_align') or [] if a['reference_database_name'] == 'UniProt']
        meta = {'desc': (pe.get('rcsb_polymer_entity') or {}).get('pdbx_description') or '',
                'org': ';'.join(o.get('scientific_name') or '' for o in pe.get('rcsb_entity_source_organism') or []),
                'aln': aln}
        for c in ids.get('auth_asym_ids') or []:
            out[c] = meta
    return out


def to_uniprot(meta, label_seq):
    for a in meta['aln']:
        for reg in a['aligned_regions']:
            if reg['entity_beg_seq_id'] <= label_seq < reg['entity_beg_seq_id'] + reg['length']:
                return a['reference_database_accession'], reg['ref_beg_seq_id'] + label_seq - reg['entity_beg_seq_id']
    return None, None


AUTO = {'ECO:0000256', 'ECO:0000259', 'ECO:0007829', 'ECO:0000313', 'ECO:0008006'}
LN_LIG = re.compile(r'^(La|Ce|Pr|Nd|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu)\(3\+\)$|lanthanide', re.I)
EXCL = re.compile(r'lanthan|lanmodulin|\bLBT\b|xoxf|exaf|pedh|de novo|designed|synthetic', re.I)


def uniprot_info(acc):
    u = r.uniprot(acc)
    if '_error' in u:
        return None
    reviewed = 'Swiss-Prot' in u.get('entryType', '')
    name = (((u.get('proteinDescription') or {}).get('recommendedName') or {}).get('fullName') or {}).get('value', '')
    if not name:
        sub = (u.get('proteinDescription') or {}).get('submissionNames') or [{}]
        name = ((sub[0].get('fullName') or {}).get('value')) or ''
    cof = ';'.join(c.get('name', '') for cm in u.get('comments', []) if cm.get('commentType') == 'COFACTOR'
                   for c in cm.get('cofactors', []))
    ca_pos, ln_pos = set(), set()
    for f in u.get('features', []):
        if f.get('type') != 'Binding site':
            continue
        lig = (f.get('ligand') or {}).get('name', '')
        ev = {e.get('evidenceCode') for e in f.get('evidences', [])}
        curated = reviewed or bool(ev - AUTO)
        lo, hi = f['location']['start'].get('value'), f['location']['end'].get('value')
        if lo is None or hi is None or not curated:
            continue
        if lig == 'Ca(2+)':
            ca_pos |= set(range(lo, hi + 1))
        elif LN_LIG.search(lig):
            ln_pos |= set(range(lo, hi + 1))
    return {'reviewed': reviewed, 'name': name, 'cofactor': cof, 'ca_pos': ca_pos, 'ln_pos': ln_pos}


def curate(pdb, data, site, meta):
    """Return (accession, n_ca_annotated, n_ln_annotated, n_mapped, exclusion_text)."""
    accs, nca, nln, nmap, texts = [], 0, 0, 0, []
    for ch, num, rn in coord_res(site):
        m = meta.get(ch)
        c = data['chains'].get(ch)
        if not m or not c or num not in c['idx']:
            continue
        res = c['res'][c['idx'][num]]
        if res.label_seq is None:
            continue
        acc, pos = to_uniprot(m, res.label_seq)
        texts += [m['desc'], m['org']]
        if not acc:
            continue
        ui = uniprot_info(acc)
        if not ui:
            continue
        nmap += 1
        accs.append(acc)
        texts += [ui['name'], ui['cofactor']]
        nca += pos in ui['ca_pos']
        nln += pos in ui['ln_pos']
    excl = sorted({t for t in texts if t and EXCL.search(t)} |
                  ({'synthetic construct'} if any('synthetic construct' in (t or '').lower() for t in texts) else set()))
    acc = max(set(accs), key=accs.count) if accs else None
    return acc, nca, nln, nmap, excl


# ---------- main ----------
site_rows, pair_rows, log = [], [], {'clusters': {}}
cache = {}
ALN = {}


def get(pdb):
    if pdb not in cache:
        cache[pdb] = load(pdb)
        cache[pdb]['meta'] = entity_meta(pdb)
    return cache[pdb]


for cl, v in sorted(cand['clusters'].items()):
    lg = {'ln_sites': 0, 'ln_sites_no_protein_residue': 0, 'pairs': 0}
    for pdb in v['ln'] + v['ca']:
        try:
            d = get(pdb)
        except Exception as e:  # noqa: BLE001
            lg.setdefault('load_errors', []).append(f'{pdb}:{e}')
            continue
        for sid, s in d['sites'].items():
            f = features(s, protein=True)
            site_rows.append({'cluster': cl, 'pdb': pdb, 'site': sid, 'metal': s['metal'], 'res': s['resolution'],
                              'occ': s['metal_occ'], 'n_protein_donors': s['n_protein_donors'], 'cn': f['cn'],
                              'n_chel': f['n_chel'], 'n_carb': f['n_carb'], 'n_water': f['n_water'],
                              'mean_d': f['mean_d'], 'rule_La': int(rule(f)) if f['cn'] else ''})
    for lpdb in v['ln']:
        if lpdb not in cache:
            continue
        L = cache[lpdb]
        for lsid, ls in L['sites'].items():
            if ls['metal'] not in LN:
                continue
            lg['ln_sites'] += 1
            cr = coord_res(ls)
            if not cr:
                lg['ln_sites_no_protein_residue'] += 1
                continue
            lpos = None
            for ch in L['st'][0]:
                for res in ch:
                    for at in res:
                        alt = at.altloc if at.altloc != '\0' else ''
                        if f"{ch.name}:{res.name}{res.seqid.num}:{at.name}" + (f":{alt}" if alt else '') == lsid:
                            lpos = np.array(at.pos.tolist())
            lchains = sorted({c for c, _, _ in cr if c in L['chains']})
            if not lchains or lpos is None:
                continue
            for cpdb in v['ca']:
                if cpdb not in cache or cpdb == lpdb:
                    continue
                C = cache[cpdb]
                opts = {}
                for lc in lchains:
                    opts[lc] = []
                    for cc, cdat in C['chains'].items():
                        ak = (lpdb, lc, cpdb, cc)
                        if ak not in ALN:
                            ALN[ak] = align_map(L['chains'][lc]['seq'], cdat['seq'])
                        ident, mp = ALN[ak]
                        if ident >= 0.95:
                            opts[lc].append((cc, mp, round(ident, 4)))
                if any(not o for o in opts.values()):
                    continue
                best = None
                for combo in itertools.islice(itertools.product(*[opts[lc] for lc in lchains]), 64):
                    P, Q = [], []
                    for lc, (cc, mp, ident) in zip(lchains, combo):
                        lch, cch = L['chains'][lc], C['chains'][cc]
                        for (_, num, _) in [x for x in cr if x[0] == lc]:
                            if num not in lch['idx']:
                                continue
                            for (k, a), xyz in bb_atoms(lch, lch['idx'][num]).items():
                                if k in mp:
                                    at = cch['res'][mp[k]].find_atom(a, '*')
                                    if at:
                                        P.append(xyz); Q.append(np.array(at.pos.tolist()))
                    # deduplicate atoms shared by overlapping windows
                    if len(P) < 8:
                        continue
                    PQ = {tuple(np.round(p, 3)): q for p, q in zip(P, Q)}
                    P2 = np.array([np.array(k) for k in PQ]); Q2 = np.array(list(PQ.values()))
                    if len(P2) < 8:
                        continue
                    R, t, rmsd = kabsch(P2, Q2)
                    for cid, cxyz in C['ca_ions']:
                        dd = float(np.linalg.norm(R @ cxyz + t - lpos))
                        if best is None or dd < best[0]:
                            best = (dd, cid, rmsd, len(P2), [c[0] for c in combo], [c[2] for c in combo])
                if best is None or best[0] > 1.5 or best[1] not in C['sites']:
                    continue
                cs = C['sites'][best[1]]
                fl, fc = features(ls, protein=True), features(cs, protein=True)
                acc, nca, nln_c, nmap_c, excl_c = curate(cpdb, C, cs, C['meta'])
                acc_l, nca_l, nln_l, nmap_l, excl_l = curate(lpdb, L, ls, L['meta'])
                lg['pairs'] += 1
                pair_rows.append({
                    'cluster': cl, 'ln_pdb': lpdb, 'ln_site': lsid, 'ln_metal': ls['metal'],
                    'ln_light': int(ls['metal'] in LIGHT), 'ln_res': ls['resolution'], 'ln_occ': ls['metal_occ'],
                    'ca_pdb': cpdb, 'ca_site': best[1], 'ca_res': cs['resolution'], 'ca_occ': cs['metal_occ'],
                    'dist_A': round(best[0], 3), 'sup_rmsd_A': round(best[2], 3), 'sup_atoms': best[3],
                    'chain_ident': min(best[5]),
                    'ln_key': ';'.join(str(n) for _, n, _ in cr),
                    'ca_key': ';'.join(str(n) for _, n, _ in coord_res(cs)),
                    'ln_cn': fl['cn'], 'ln_chel': fl['n_chel'], 'ln_carb': fl['n_carb'], 'ln_water': fl['n_water'],
                    'ln_mean_d': round(fl['mean_d'], 4) if fl['mean_d'] else '', 'ln_rule_La': int(rule(fl)),
                    'ln_protein_donors': ls['n_protein_donors'],
                    'ca_cn': fc['cn'], 'ca_chel': fc['n_chel'], 'ca_carb': fc['n_carb'], 'ca_water': fc['n_water'],
                    'ca_mean_d': round(fc['mean_d'], 4) if fc['mean_d'] else '', 'ca_rule_La': int(rule(fc)),
                    'ca_protein_donors': cs['n_protein_donors'],
                    'uniprot': acc or acc_l or '', 'ca_site_res_mapped': nmap_c, 'ca_site_res_CaBinding': nca,
                    'ln_site_res_LnBinding': nln_l, 'exclusion': ' | '.join(sorted(set(excl_c) | set(excl_l))),
                })
    log['clusters'][cl] = lg
    print(cl, lg, flush=True)

os.makedirs(os.path.join(WS, 'out'), exist_ok=True)
for name, rows in (('sites.tsv', site_rows), ('pairs.tsv', pair_rows)):
    with open(os.path.join(WS, 'out', name), 'w') as fh:
        if rows:
            cols = list(rows[0].keys())
            fh.write('\t'.join(cols) + '\n')
            for x in rows:
                fh.write('\t'.join(str(x[c]) for c in cols) + '\n')
json.dump(log, open(os.path.join(WS, 'out', 'pair_log.json'), 'w'), indent=1)
print('sites', len(site_rows), 'pairs', len(pair_rows))
