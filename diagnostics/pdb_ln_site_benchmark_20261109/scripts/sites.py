"""Label lanthanide sites in every Ln entry (full crystal with symmetry images).
Writes WS/out/ln_sites.tsv (one row per Ln ion, first conformer) and WS/out/ln_ions.json (all Ln
positions for capacity labelling).
Usage: python -I sites.py WS SCRIPTS AUDIT_CIF"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import gemmi

WS, SCR, ACIF = sys.argv[1:4]
sys.path.insert(0, SCR)
import common as C  # noqa: E402

C.CIF_DIRS[:] = [os.path.join(WS, 'cif'), ACIF]
PRIMARY = 3.0  # protein_shells PRIMARY; first-shell O/N cutoff
CARBOX = {('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2')}
NEUTRAL_SC = {('ASN', 'OD1'), ('GLN', 'OE1'), ('SER', 'OG'), ('THR', 'OG1'), ('TYR', 'OH')}
HISN = {('HIS', 'ND1'), ('HIS', 'NE2')}
disc = json.load(open(os.path.join(WS, 'out', 'discovery.json')))


def dtype(res, atom):
    if atom == 'O' or atom == 'OXT':
        return 'bbO'
    if (res, atom) in CARBOX:
        return 'carbox'
    if (res, atom) in NEUTRAL_SC:
        return 'scO'
    if (res, atom) in HISN:
        return 'hisN'
    return 'otherON'


rows, ions = [], {}
for e in sorted(disc['ln_info']):
    st = C.read(e)
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    ent_of = {}
    for ch in model:
        for en in disc['ln_info'][e]['entities']:
            if ch.name in en['chains']:
                ent_of[ch.name] = en
    ca_b = [a.b_iso for ch in model for r in ch if r.name in C.AA for a in r if a.name == 'CA' and C.alt_ok(a)]
    medb = float(np.median(ca_b)) if ca_b else float('nan')
    lns = []
    for ch in model:
        for r in ch:
            for a in r:
                if a.element.name in C.LN_EL and C.alt_ok(a):
                    lns.append((ch, r, a))
    ions[e] = [{'site': f'{ch.name}:{r.name}{r.seqid.num}:{a.name}', 'xyz': a.pos.tolist()} for ch, r, a in lns]
    for ch, r, a in lns:
        site = f'{ch.name}:{r.name}{r.seqid.num}:{a.name}'
        don, other_ln5, cn = [], 0, 0
        for mk in ns.find_atoms(a.pos, '\0', radius=5.0):
            cra = mk.to_cra(model)
            b = cra.atom
            if not C.alt_ok(b):
                continue
            p = st.cell.find_nearest_pbc_position(a.pos, b.pos, mk.image_idx)
            d = p.dist(a.pos)
            if d < 0.1:
                continue
            if b.element.name in C.LN_EL:
                other_ln5 += 1
                continue
            if b.element.name not in ('O', 'N') or d > PRIMARY:
                continue
            cn += 1
            if cra.residue.name not in C.AA or cra.residue.entity_type != gemmi.EntityType.Polymer:
                continue
            # Ln position expressed in the donor chain's own (unmoved) frame
            translated = abs(b.pos.dist(a.pos) - d) > 0.01
            img = int(mk.image_idx) if mk.image_idx else (999 if translated else 0)
            if img:
                q = st.cell.find_nearest_pbc_position(b.pos, a.pos, mk.image_idx, True)
                if abs(q.dist(b.pos) - d) > 0.05:
                    q = None
            else:
                q = a.pos
            don.append({'chain': cra.chain.name, 'img': img, 'res': cra.residue.name,
                        'num': cra.residue.seqid.num, 'atom': b.name, 'type': dtype(cra.residue.name, b.name),
                        'd': d, 'q': None if q is None else q.tolist()})
        groups = defaultdict(list)
        for x in don:
            groups[(x['chain'], x['img'])].append(x)
        n_same = n_other = n_sym = 0
        if groups:
            # attributed group: most donors; ties -> image 0 first, then chain name
            gkey = sorted(groups, key=lambda k: (-len(groups[k]), k[1] != 0, k[0]))[0]
            for k, v in groups.items():
                if k == gkey:
                    n_same += len(v)
                elif k[1] != gkey[1]:
                    n_sym += len(v)
                else:
                    n_other += len(v)
            if len(groups) == 1:
                cls = 'single-chain'
            elif any(k[1] != gkey[1] for k in groups):
                cls = 'symmetry-bridging'
            else:
                cls = 'inter-chain'
            q = groups[gkey][0]['q']
        else:
            gkey, cls, q = (None, None), 'no-protein', None
        res_keys = {(x['chain'], x['img'], x['num']) for x in don}
        carb = Counter((x['chain'], x['img'], x['num']) for x in don if x['type'] == 'carbox')
        en = ent_of.get(gkey[0], {})
        rows.append({'pdb': e, 'site': site, 'metal': a.element.name, 'occ': round(a.occ, 3),
                     'b': round(a.b_iso, 2), 'b_rel': round(a.b_iso / medb, 3) if medb == medb and medb > 0 else '',
                     'res': disc['ln_info'][e]['res'], 'cn': cn, 'n_prot_don': len(don),
                     'n_prot_res': len(res_keys), 'n_same_chain': n_same, 'n_other_chain_au': n_other,
                     'n_symmetry': n_sym, 'cls': cls, 'cluster5': int(other_ln5 > 0),
                     'n_chel_carbox': sum(1 for v in carb.values() if v >= 2),
                     'types': ';'.join(f'{k}:{v}' for k, v in sorted(Counter(x['type'] for x in don).items())),
                     'chain': gkey[0] or '', 'img': gkey[1] if gkey[1] is not None else '',
                     'qx': '' if q is None else round(q[0], 3), 'qy': '' if q is None else round(q[1], 3),
                     'qz': '' if q is None else round(q[2], 3),
                     'entity': en.get('entity', ''), 'c95': en.get('c95') or '', 'c30': en.get('c30') or '',
                     'donors': ';'.join(f"{x['chain']}{'@' + str(x['img']) if x['img'] else ''}:{x['res']}{x['num']}:{x['atom']}:{x['d']:.2f}" for x in don)})
    print(e, len(lns), flush=True)
cols = list(rows[0])
with open(os.path.join(WS, 'out', 'ln_sites.tsv'), 'w') as fh:
    fh.write('\t'.join(cols) + '\n')
    for x in rows:
        fh.write('\t'.join(str(x[c]) for c in cols) + '\n')
json.dump(ions, open(os.path.join(WS, 'out', 'ln_ions.json'), 'w'))
print('ions', len(rows), Counter(x['cls'] for x in rows if x['n_prot_don'] >= 2))
