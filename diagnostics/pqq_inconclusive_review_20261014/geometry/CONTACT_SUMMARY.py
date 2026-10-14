#!/usr/bin/env python3
"""Correct only the saved inter-residue contact descriptor and compact its peers."""
import argparse, csv, hashlib, importlib.util, json, os, resource, socket, sys, time, traceback
from pathlib import Path

def read(p): return json.loads(Path(p).read_text())
def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path']) == p, p
    return Path(p['path'])
def write(p, x):
    with Path(p).open('x') as f: json.dump(x, f, indent=2, allow_nan=False); f.write('\n')

def main(config):
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURM_JOB_NODELIST'):
        raise RuntimeError('actual Slurm allocation required before scientific imports')
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
        os.environ[k] = '1'
    cfg = read(config)
    for p in (cfg['driver'], cfg['plan'], cfg['manifest'], cfg['result']): verify(p)
    m = read(cfg['manifest']['path']); result = read(cfg['result']['path'])
    assert result['manifest'] == cfg['manifest'] and len(result['rows']) == 48
    for p in [m['driver'], m['helper'], *m['code'].values(), *m['inputs']]: verify(p)
    out = Path(cfg['output']); out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    receipt = {'job_id': os.environ['SLURM_JOB_ID'], 'hostname': socket.gethostname(),
               'config': pin(config), 'status': 'running', 'new_energy_calls': 0, 'new_folds': 0}
    try:
        import numpy as np, gemmi
        root = Path(m['driver']['path']).parents[3]
        sys.path.insert(0, str(root/'scripts'))
        from mace_site_kinematics import Kinematics
        import pqq_microstates as pqq
        spec = importlib.util.spec_from_file_location('saved_adapter', m['driver']['path'])
        a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
        spec = importlib.util.spec_from_file_location('saved_helper', m['helper']['path'])
        b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
        b.np, b.gemmi, b.pqq = np, gemmi, pqq
        rows, table = [], []
        for c, r in zip(m['cases'], result['rows']):
            assert c['case_key'] == r['case_key'] and r['status'] == 'complete'
            data = read(c['mapping']['path'])['context']; kin = Kinematics(data)
            meta = data['source_atom_metadata']; raw, _, _ = a.raw_map(c, b)
            mapped = {b.key(v) for v in meta if v}
            edges = {i: set() for i in range(len(meta))}
            for i, j in data['bonds']: edges[i].add(j); edges[j].add(i)
            heavy = [i for i,v in enumerate(meta) if v and v['element'] not in ('H','D','Ca','La')]
            omitted = [(k,v) for k,v in raw.items() if k not in mapped and v['element'] not in ('H','D','Ca','La')]
            omitted_xyz = np.asarray([v['xyz'] for k,v in omitted])
            donor_residues = {b.selkey(v,'')[:3] for v in c['roles'].values() if v['resname'] in ('ASP','GLU','ASN')}
            origin = None; contacts = {}; flat = {'case_key': r['case_key'], 'case_id': r['case_id'],
                'protein': r['biological_group'], 'conditioning': r['source_conditioning_metal'],
                'target': r['target'], 'origin_CN': r['local']['origin']['typed_declared_role_CN_3p1A'],
                'metal_pLDDT': r['raw_metrics']['confidence']['CIF_atom_plddt']['active_metal']['mean'],
                'metal_to_PQQ_PAE_A': r['raw_metrics']['confidence'].get('directed_PAE_A',{}).get('metal_to_PQQ',{}).get('mean')}
            for name in ('origin','adaptive_Ca','adaptive_La'):
                q = np.zeros(len(data['modes'])) if name == 'origin' else np.asarray(read(c['proposals'][name.removeprefix('adaptive_')]['path'])['proposal']['full_q'])
                physical, _, _, _ = kin.evaluate(q); sy, xyz = b.xyz_read(c['geometries'][name])
                physical[kin.metal] = xyz[0]
                for i, kind, *rest in data['core_links']:
                    if kind == 'source': physical[rest[0]] = xyz[i]
                if name == 'origin': origin = physical.copy()
                pairs = []
                for idx, i in enumerate(heavy):
                    excluded = edges[i] | {j for k in edges[i] for j in edges[k]}
                    for j in heavy[idx+1:]:
                        if j in excluded or b.key(meta[i])[:3] == b.key(meta[j])[:3]: continue
                        pairs.append((float(np.linalg.norm(physical[i]-physical[j])), i, j))
                ranked = [{'distance_A': d, 'atoms': [meta[i],meta[j]]} for d,i,j in sorted(pairs)[:12]]
                # Standard Asp/Glu/Asn sidechains have no inter-residue covalent
                # bonds. Exclude all backbone atoms and the complete same residue.
                donor_ids = [i for i in heavy if b.key(meta[i])[:3] in donor_residues and meta[i]['atom'] not in ('N','CA','C','O','OXT')]
                near = []
                for i in donor_ids:
                    distances = np.linalg.norm(omitted_xyz-physical[i], axis=1)
                    eligible = [j for j,(k,v) in enumerate(omitted) if k[:3] != b.key(meta[i])[:3]]
                    closest = sorted(eligible, key=lambda j: (distances[j], omitted[j][0]))[:12]
                    for j in closest:
                        k,v = omitted[j]
                        if k[:3] == b.key(meta[i])[:3]: continue
                        near.append((float(distances[j]), i, k, v['resname']))
                nearest = [{'distance_A': d, 'donor': meta[i], 'omitted_atom': {'chain':k[0],'resnum':k[1],'icode':k[2],'atom':k[3],'resname':rn},
                            'donor_displacement_A':float(np.linalg.norm(physical[i]-origin[i]))} for d,i,k,rn in sorted(near)[:12]]
                frame_residual = r['local']['origin']['raw_source_mapped_heavy_displacements']['max_A']
                contacts[name] = {'inter_residue_mapped_heavy': ranked,
                    'omitted_source_to_standard_donor_sidechain': nearest if frame_residual <= 1e-10 else None,
                    'source_coordinate_frame_residual_A': frame_residual,
                    'omitted_scope': 'all retained standard Asp/Glu/Asn sidechain heavy atoms; origin and proposals, same residue/backbone excluded; source membership exact',
                    'omitted_status': 'available' if frame_residual <= 1e-10 else 'unavailable_coordinate_frame'}
                if name != 'origin':
                    z = name.removeprefix('adaptive_'); e = r['optimizer_evidence'][name]
                    flat[z+'_selected'] = r['selected'][z]
                    flat[z+'_boundary'] = e['boundary_flag']
                    flat[z+'_maximum_heavy_motion_A'] = e['final_geometry']['maximum_heavy_displacement_A']
                    flat[z+'_selected_residual_kcal_mol_A'] = e['residual']['selected_max_load']
                    flat[z+'_omitted_angular_residual_kcal_mol_A'] = e['residual']['omitted_angular_max_load']
            rows.append({'case_key':r['case_key'], 'case_id':r['case_id'], 'target':r['target'], 'contacts':contacts})
            table.append(flat)
        final = {'prior_result': cfg['result'], 'config': pin(config), 'terminal': True, 'rows': rows,
                 'supersedes_only':'closest_mapped_nonbonded_heavy_pairs/nonbonded_scope; v1 included fixed PQQ covalent bonds absent from kinematic graph',
                 'interpretation':'inter-residue distances, no invented clash cutoff; proximity/strain hypothesis does not prove missing energy or incorrect class'}
        write(out/'RESULT.json',final)
        with (out/'PEERS.tsv').open('x') as f:
            w=csv.DictWriter(f,fieldnames=list(table[0]),delimiter='\t');w.writeheader();w.writerows(table)
        receipt.update(status='complete', result=pin(out/'RESULT.json'), table=pin(out/'PEERS.tsv'))
    except Exception as exc:
        receipt.update(status='failed', reason=str(exc), traceback=traceback.format_exc())
        raise
    finally:
        u=resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(wall_seconds=time.monotonic()-start, CPU_seconds=u.ru_utime+u.ru_stime, peak_RSS_KiB=u.ru_maxrss)
        write(out/'EXECUTION.json',receipt)

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);a=p.parse_args();main(a.config)
