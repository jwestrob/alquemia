#!/usr/bin/env python3
"""Saved geometry only; execution, including numerical parsing, requires Slurm."""
import csv, hashlib, json, math, os, resource, socket, sys, time, traceback
from pathlib import Path

def pin(path):
    p = Path(path).resolve(); return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def checked(p):
    assert pin(p['path']) == p, 'Input hash changed: ' + p['path']; return Path(p['path'])
def read(p): return json.loads(checked(p).read_text())
def write(path, value): Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
def key(s): return tuple(s.get(k, '') for k in ('chain','resnum','insertion_code','resname','atom','element'))
def xyz(p):
    lines = checked(p).read_text().splitlines(); a = [s.split() for s in lines[2:] if s.strip()]
    assert int(lines[0]) == len(a) == 204
    return [s[0] for s in a], [tuple(map(float, s[1:4])) for s in a], lines[2:]

def analyze(r):
    a = read(r['geometry_admission']); cp = read(r['source_preparation']); state = read(r['state_signature'])
    maps = {z: read(r['maps'][z]) for z in ('Ca','La')}
    assert r['status'] == a['status'] == 'prepared' and a['pass'] and maps['Ca'] == maps['La']
    n = r['named_Asn_indices_zero_based']; moved = set(r['moved_indices_zero_based'])
    assert moved == {n[x] for x in ('OD1','ND2','HD21','HD22')}
    entries = {v['qm_index']: v for v in cp['mapping']['source_to_qm']}
    assert set(entries) == set(range(1,204))
    source_qm = {key(v['source']): i for i,v in entries.items() if 'source' in v}
    for name,i in n.items():
        s = entries[i]['source']; assert (s['chain'],s['resnum'],s['resname'],s['atom']) == ('A',290,'ASN',name)
    kin = maps['Ca']['context']; physical_qm = {i: source_qm[key(s)] for i,s in enumerate(kin['source_atom_metadata']) if s and key(s) in source_qm}
    bonds = {frozenset((physical_qm[i],physical_qm[j])) for i,j in kin['bonds'] if i in physical_qm and j in physical_qm}
    parents = {tuple(h): tuple(p) for h,p in state['H_covalent_parents']}
    for h in ('HD21','HD22'):
        assert parents[key(entries[n[h]]['source'])] == key(entries[n['ND2']]['source'])
        assert frozenset((n[h],n['ND2'])) in bonds
    pair = {}
    for z in ('Ca','La'):
        c = r['matrix'][z]['amide_flip']; assert c['geometry_admission'] == r['geometry_admission']
        assert c['charge'] == state['endpoint_charges'][z] and c['multiplicity'] == 1 and c['atom_count'] == 204
        before = xyz(c['original_xyz']); after = xyz(c['xyz']); assert before[0] == after[0] and before[0][0] == z
        assert all(before[2][i] == after[2][i] for i in range(204) if i not in moved)
        pair[z] = (before,after)
    for t in (0,1):
        assert pair['Ca'][t][1] == pair['La'][t][1] and pair['Ca'][t][0][1:] == pair['La'][t][0][1:]
    symbols,before,_ = pair['Ca'][0]; _,after,_ = pair['Ca'][1]
    distances = {stage: {name: math.dist(coords[0],coords[n[name]]) for name in ('OD1','ND2')} for stage,coords in (('origin',before),('amide_flip',after))}
    caps = set(a['cap_indices']); contacts = []
    for i in sorted(moved):
        for j in range(204):
            if j in moved or j in caps or 'H' not in (symbols[i],symbols[j]) or frozenset((i,j)) in bonds: continue
            contacts.append({'moving_index':i,'other_index':j,'moving':entries[i], 'other':entries.get(j,{'kind':'metal','qm_index':0}),
                             'before_A':math.dist(before[i],before[j]),'after_A':math.dist(after[i],after[j])})
    contacts.sort(key=lambda v:(v['after_A'],v['moving_index'],v['other_index']))
    shortened = [c for c in contacts if c['after_A'] < c['before_A']]
    return {'source_case_id':r['source_case_id'],'case_key':r['case_key'],'record_id':r['record_id'],'conditioning':r['conditioning'],'sample':r['sample'],
            'status':'complete','paired_Ca_La_geometry_exact':True,'metal_donor_distances_A':distances,'H_parents_verified':{'HD21':'ND2','HD22':'ND2'},
            'rigid_amide_max_error_A':a['rigid_amide_max_error_A'],'covalent_max_error_A':a['covalent_max_error_A'],
            'shortest_included_nonbonded_H_pair':contacts[0] if contacts else None,
            'shortest_shortened_included_nonbonded_H_pair':shortened[0] if shortened else None,
            'saved_nearest_contacts':a['contacts']['nearest_contacts'],'saved_nearest_omitted_contacts':a['contacts']['nearest_omitted_source_contacts'],
            'inputs':{k:r[k] for k in ('maps','state_signature','source_preparation','geometry_admission','matrix')}}

def main():
    assert os.environ.get('SLURM_JOB_ID','').isdigit(), 'Actual Slurm allocation required'
    assert int(os.environ.get('SLURM_CPUS_ON_NODE','0').split('(')[0]) >= 1
    mp = {'path':str(Path(sys.argv[1]).resolve()),'sha256':sys.argv[2]}; m = read(mp)
    assert pin(__file__) == m['implementation']['ANALYZE.py'] and pin(sys.executable) == m['python']
    checked(m['implementation']['RUN.sbatch']); checked(m['method'])
    out = Path(m['output']); out.mkdir(parents=True,exist_ok=False); start=time.monotonic(); cpu=time.process_time(); rows=[]
    receipt={'status':'failed','manifest':mp,'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'new_energy_calls':0}
    try:
        p = read(m['preparation']); assert p['status'] == 'complete' and p['prepared_cells'] == 20
        assert [r['source_case_id'] for r in p['rows']] == m['source_ids'] and len(p['rows']) == 10
        for r in p['rows']:
            try: rows.append(analyze(r))
            except Exception: rows.append({'source_case_id':r['source_case_id'],'status':'unavailable','error':traceback.format_exc()})
        result={'status':'complete' if all(r['status']=='complete' for r in rows) else 'partial','manifest':mp,'declared_sources':10,'rows':rows,'method':m['method']}
        write(out/'RESULT.json',result)
        with (out/'TABLE.csv').open('w',newline='') as f:
            w=csv.writer(f); w.writerow(['source_case_id','conditioning','sample','status','origin_OD1_A','origin_ND2_A','flip_OD1_A','flip_ND2_A'])
            for r in rows:
                d=r.get('metal_donor_distances_A',{}); w.writerow([r['source_case_id'],r.get('conditioning'),r.get('sample'),r['status'],*[d.get(s,{}).get(n) for s in ('origin','amide_flip') for n in ('OD1','ND2')]])
        receipt['status']=result['status']; receipt['result']=pin(out/'RESULT.json')
    except Exception: receipt['error']=traceback.format_exc(); raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__ == '__main__': sys.exit(main())
