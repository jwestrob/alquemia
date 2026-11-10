"""Step 1: find Ln entries, their 95 % clusters and Ca partner entries; apply the download cap.
Writes WS/out/candidates.json. Usage: python -I discover.py WS"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rcsb as r  # noqa: E402

WS = sys.argv[1]
r.set_cache(os.path.join(WS, 'api'))
LN = 'LA CE PR ND SM EU GD TB DY HO ER TM YB LU'.split()
CAP, PER = 1000, 8
xr = r.terminal('exptl.method', 'exact_match', 'X-RAY DIFFRACTION')
res = r.terminal('rcsb_entry_info.resolution_combined', 'less_or_equal', 2.5)
cid = 'rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id'
ln_ids = sorted(r.search(r.AND(xr, res, r.terminal(cid, 'in', LN))))
ca_ids = set(r.search(r.AND(xr, res, r.terminal(cid, 'exact_match', 'CA'))))


def entry_info(ids):
    out = {}
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        q = ('{entries(entry_ids:%s){rcsb_id rcsb_entry_info{resolution_combined} '
             'polymer_entities{rcsb_id entity_poly{rcsb_entity_polymer_type} '
             'rcsb_polymer_entity_group_membership{group_id aggregation_method similarity_cutoff}}}}') % json.dumps(chunk)
        d = r.graphql(q)
        for e in d['data']['entries'] or []:
            cl = set()
            for pe in e['polymer_entities'] or []:
                if (pe.get('entity_poly') or {}).get('rcsb_entity_polymer_type') != 'Protein':
                    continue
                for g in pe.get('rcsb_polymer_entity_group_membership') or []:
                    if g['aggregation_method'] == 'sequence_identity' and int(g['similarity_cutoff']) == 95:
                        cl.add(g['group_id'])
            res_ = (e['rcsb_entry_info'] or {}).get('resolution_combined') or [99]
            out[e['rcsb_id']] = {'res': min(res_), 'clusters': sorted(cl)}
    return out


ln_info = entry_info(ln_ids)
clusters = sorted({c for v in ln_info.values() for c in v['clusters']})
# partner candidates: Ca entries with a polymer entity in any of these clusters (verified below)
gm = 'rcsb_polymer_entity_group_membership.group_id'
part = set()
for i in range(0, len(clusters), 200):
    part |= set(x.split('_')[0] for x in r.search(r.AND(r.terminal(gm, 'in', clusters[i:i + 200])), 'polymer_entity'))
part = sorted((part & ca_ids) - set())
ca_info = entry_info(part)
by_cl = {}
for c in clusters:
    lns = sorted([e for e, v in ln_info.items() if c in v['clusters']], key=lambda e: (ln_info[e]['res'], e))
    cas = sorted([e for e, v in ca_info.items() if c in v['clusters'] and e not in lns], key=lambda e: (ca_info[e]['res'], e))
    if lns and cas:
        by_cl[c] = {'ln': lns, 'ca': cas}
order = sorted(by_cl)
random.Random(42).shuffle(order)
need_all = len({e for c in by_cl for e in by_cl[c]['ln'][:PER] + by_cl[c]['ca'][:PER]})
chosen, files = [], set()
for c in order:
    add = set(by_cl[c]['ln'][:PER] + by_cl[c]['ca'][:PER])
    if len(files | add) > CAP:
        continue
    files |= add
    chosen.append(c)
json.dump({'n_ln_entries': len(ln_ids), 'n_ca_entries': len(ca_ids), 'n_ln_clusters': len(clusters),
           'n_clusters_with_ca_partner': len(by_cl), 'files_needed_all_clusters_capped_per_cluster': need_all,
           'n_clusters_chosen': len(chosen), 'n_files': len(files), 'sampled': len(chosen) < len(by_cl),
           'clusters': {c: {'ln': by_cl[c]['ln'][:PER], 'ca': by_cl[c]['ca'][:PER],
                            'n_ln_total': len(by_cl[c]['ln']), 'n_ca_total': len(by_cl[c]['ca'])} for c in chosen},
           'entry_res': {**{e: ln_info[e]['res'] for e in ln_info}, **{e: ca_info[e]['res'] for e in ca_info}}},
          open(os.path.join(WS, 'out', 'candidates.json'), 'w'), indent=1)
print('Ln entries', len(ln_ids), 'clusters', len(clusters), 'with Ca partner', len(by_cl),
      'files needed', need_all, 'chosen clusters', len(chosen), 'files', len(files))
