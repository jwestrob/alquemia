"""Benchmark discovery: all Ln X-ray entries (<=2.5 A), their protein entities' 95 %/30 % clusters,
and lanthanide-free partner deposits (X-ray <=3.0 A) of the same 95 % clusters.
Reuses the audit's rcsb.py client (imported by path; unchanged). Single-threaded, polite.
Usage: python -I discover_bench.py WS AUDIT_SCRIPTS_DIR"""
import json
import os
import sys

WS, AUD = sys.argv[1], sys.argv[2]
sys.path.insert(0, AUD)
import rcsb as r  # noqa: E402

r.set_cache(os.path.join(WS, 'api'))
LN = 'LA CE PR ND SM EU GD TB DY HO ER TM YB LU'.split()
LNLIKE = set(LN) | {'Y1', 'YT3', 'SC', 'U1', 'TH', 'AM', 'CM'}  # Y3+/Sc3+/actinide ions excluded from partners
METAL_IONS = set('CA MG ZN MN FE FE2 CU CU1 CU3 NI CO CO3 NA K CD HG SR BA CS RB LI TL PB AG AU AU3 PT PD IR OS RU RH '
                 'MN3 CR V GA IN SN SB BI AL BE W MO 3CO 3NI CD1 PT4 ZN2 OHG HG2 MG8 NA2 4MO 6MO'.split())
xr = r.terminal('exptl.method', 'exact_match', 'X-RAY DIFFRACTION')
res25 = r.terminal('rcsb_entry_info.resolution_combined', 'less_or_equal', 2.5)
res30 = r.terminal('rcsb_entry_info.resolution_combined', 'less_or_equal', 3.0)
cid = 'rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id'
ln_ids = sorted(r.search(r.AND(xr, res25, r.terminal(cid, 'in', LN))))


def entry_info(ids):
    out = {}
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        q = ('{entries(entry_ids:%s){rcsb_id rcsb_entry_info{resolution_combined} '
             'nonpolymer_entities{nonpolymer_comp{chem_comp{id}}} '
             'polymer_entities{rcsb_id entity_poly{rcsb_entity_polymer_type} '
             'rcsb_polymer_entity_container_identifiers{entity_id auth_asym_ids} '
             'rcsb_polymer_entity_group_membership{group_id aggregation_method similarity_cutoff}}}}') % json.dumps(chunk)
        d = r.graphql(q)
        for e in d['data']['entries'] or []:
            ents = []
            for pe in e['polymer_entities'] or []:
                if (pe.get('entity_poly') or {}).get('rcsb_entity_polymer_type') != 'Protein':
                    continue
                g = {}
                for m in pe.get('rcsb_polymer_entity_group_membership') or []:
                    if m['aggregation_method'] == 'sequence_identity':
                        g[int(m['similarity_cutoff'])] = m['group_id']
                ids_ = pe['rcsb_polymer_entity_container_identifiers']
                ents.append({'entity': ids_['entity_id'], 'chains': ids_.get('auth_asym_ids') or [],
                             'c95': g.get(95), 'c30': g.get(30)})
            comps = sorted({(n.get('nonpolymer_comp') or {}).get('chem_comp', {}).get('id')
                            for n in e.get('nonpolymer_entities') or []} - {None})
            res_ = (e['rcsb_entry_info'] or {}).get('resolution_combined') or [99]
            out[e['rcsb_id']] = {'res': min(res_), 'entities': ents, 'comps': comps}
    return out


ln_info = entry_info(ln_ids)
c95 = sorted({en['c95'] for v in ln_info.values() for en in v['entities'] if en['c95']})
gm = 'rcsb_polymer_entity_group_membership.group_id'
cand = set()
for i in range(0, len(c95), 100):
    hits = r.search(r.AND(xr, res30, r.terminal(gm, 'in', c95[i:i + 100])), 'polymer_entity')
    cand |= {x.split('_')[0] for x in hits}
cand = sorted(cand - set(ln_ids))
p_info = entry_info(cand)


def kind(comps):
    cs = set(comps)
    if cs & LNLIKE:
        return 'ln_like'
    if 'CA' in cs:
        return 'ca'
    if cs & METAL_IONS:
        return 'other_metal'
    return 'metal_free'


partners = {}
for e, v in p_info.items():
    k = kind(v['comps'])
    if k == 'ln_like':
        continue
    for en in v['entities']:
        if en['c95'] in c95:
            partners.setdefault(en['c95'], {}).setdefault(e, {'res': v['res'], 'kind': k})
out = {'n_ln_entries': len(ln_ids), 'n_ln_entries_with_info': len(ln_info), 'n_c95': len(c95),
       'n_partner_candidates': len(cand), 'ln_info': ln_info,
       'partners': {c: partners.get(c, {}) for c in c95},
       'partner_info': {e: p_info[e] for e in sorted({e for c in partners.values() for e in c})}}
json.dump(out, open(os.path.join(WS, 'out', 'discovery.json'), 'w'), indent=0)
from collections import Counter  # noqa: E402
kc = Counter(v['kind'] for c in partners.values() for v in c.values())
print('Ln entries', len(ln_ids), 'info', len(ln_info), 'c95', len(c95), 'c95 with partner', len(partners),
      'partner entries', len(out['partner_info']), 'kinds (cluster-entry pairs)', dict(kc))
