"""One saved-only force-direction decomposition of the admitted donor fit."""
import argparse
import datetime
import hashlib
import json
import math
import mmap
import os
from pathlib import Path
import re
import time

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'

def pin(path):
    p = Path(path).resolve()
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return {'path': str(p), 'sha256': h.hexdigest()}

def root_field(path, key):
    # Read only the selected terminal record from the large immutable fit log.
    marker = ('  "' + key + '": ').encode()
    with path.open('rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as m:
        start = m.find(b'\n' + marker) + 1
        assert start > 0, key
        start += len(marker)
        next_key = re.search(rb'\n  "[^"\n]+": ', m[start:])
        end = start + next_key.start() if next_key else len(m) - 2
        return json.loads(m[start:end].strip().removesuffix(b','))

def sub(a, b): return [x-y for x, y in zip(a, b)]
def dot(a, b): return sum(x*y for x, y in zip(a, b))
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def main(out):
    start, cpu = time.monotonic(), time.process_time()
    assert all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    declaration = json.loads(Path(__file__).with_name('DECLARATION.json').read_text())
    fit_path = BASE / 'paired_donor_fit_v1/run_v1/execution/RESULT.json'
    training_path = BASE / 'matched_dry_v1/training_v1/RESULT.json'
    fit_pin, training_pin = pin(fit_path), pin(training_path)
    assert fit_pin['sha256'] == declaration['inputs']['selected_donor_result_sha256']
    assert training_pin['sha256'] == declaration['inputs']['training_result_sha256']
    selected = root_field(fit_path, 'selected')
    assert selected['status'] == 'complete'
    data = json.loads(training_path.read_text())
    model, score = selected['model'], selected['score']
    fixtures, max_rms_error, max_component_error = [], 0., 0.
    for kind, values in score['per_fixture'].items():
        for name, saved in values.items():
            row = data['rows'][name+'_La']
            positions = row['geometry']['positions_A']
            axis = row['geometry']['radial_direction']
            assert abs(dot(axis, axis)-1) < 1e-10
            def residual(metal):
                return [sub(a,b) for a,b in zip(model[name+'_'+metal]['forces'], data['rows'][name+'_'+metal]['reference']['forces'])]
            errors = [sub(a,b) for a,b in zip(residual('Dy'),residual('La'))] if kind == 'paired' else residual(kind)
            bodies = row['bodies']; categories = sorted({b['category'] for b in bodies})
            contributions = dict(radial_translation=0., transverse_translation=0., rotation=0.)
            body_records, raw = [], []
            for category in categories:
                group = [b for b in bodies if b['category'] == category]
                for body in group:
                    ids = body['indices']
                    center = [sum(positions[i][k] for i in ids)/len(ids) for k in range(3)]
                    net = [sum(errors[i][k] for i in ids) for k in range(3)]
                    torques = [cross(sub(positions[i],center),errors[i]) for i in ids]
                    torque = [sum(t[k] for t in torques) for k in range(3)]
                    radial = dot(net,axis); transverse = sub(net,[radial*x for x in axis])
                    nkeys = 1 if category == 'ion' else 2
                    weight = 1/(len(categories)*len(group)*nkeys*3)
                    contributions['radial_translation'] += weight * radial**2
                    contributions['transverse_translation'] += weight * dot(transverse,transverse)
                    raw.extend(net)
                    if nkeys == 2:
                        contributions['rotation'] += weight * dot(torque,torque)
                        raw.extend(torque)
                    body_records.append({'category':category,'indices':ids,'net_error':net,'torque_over_1A_error':torque,'radial_error':radial,'transverse_error_norm':math.sqrt(dot(transverse,transverse)),'torque_over_1A_norm':math.sqrt(dot(torque,torque))})
            rms_error = abs(sum(contributions.values())-saved['rigid_RMS_kcal_mol_A']**2)
            component_error = max(abs(a-b) for a,b in zip(raw,saved['raw_body_components_kcal_mol_A']))
            assert len(raw) == len(saved['raw_body_components_kcal_mol_A'])
            assert rms_error < 1e-8 and component_error < 1e-8
            max_rms_error = max(max_rms_error,rms_error); max_component_error = max(max_component_error,component_error)
            fixtures.append({'kind':kind,'fixture':name,'environment':'dry' if name.startswith('matched_dry_') else 'hydrated','saved_RMS':saved['rigid_RMS_kcal_mol_A'],'saved_max_component':saved['maximum_body_component_error_kcal_mol_A'],'saved_passed':saved['passed'],'RMS_squared_contributions':contributions,'bodies':body_records})
    assert len(fixtures) == 36
    aggregates = {}
    for kind in ('La','Dy','paired'):
        for environment in ('dry','hydrated'):
            rows = [r for r in fixtures if r['kind']==kind and r['environment']==environment]
            c = {key:sum(r['RMS_squared_contributions'][key] for r in rows)/len(rows) for key in contributions}
            total = sum(c.values())
            aggregates[kind+'/'+environment] = {'fixtures':len(rows),'failed':sum(not r['saved_passed'] for r in rows),'mean_RMS_squared_components':c,'fraction_of_squared_error':{k:v/total for k,v in c.items()},'pooled_RMS':math.sqrt(total)}
    result = {'status':'saved_force_direction_decomposition_complete','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'implementation':pin(__file__),'declaration':pin(Path(__file__).with_name('DECLARATION.json')),'fit':fit_pin,'training':training_pin,'fixtures':fixtures,'aggregates':aggregates,'reconstruction':{'fixtures':36,'max_RMS_squared_error':max_rms_error,'max_body_component_error':max_component_error},'interpretation':declaration['interpretation'],'new_native_calls':0,'threads':1,'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}}
    out.mkdir(parents=True,exist_ok=False)
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines = ['# Saved donor force-direction diagnostic','','All36 saved rigid-force results reconstruct from atomic forces. No new molecular calls or fitting.','', '| Subset | Failed / total | Pooled RMS | Radial share | Transverse share | Rotation share |','|---|---:|---:|---:|---:|---:|']
    for name,a in aggregates.items():
        p=a['fraction_of_squared_error'];lines.append(f"| {name} | {a['failed']}/{a['fixtures']} | {a['pooled_RMS']:.3f} | {p['radial_translation']:.1%} | {p['transverse_translation']:.1%} | {p['rotation']:.1%} |")
    lines += ['', 'RMS units are kcal/mol/Å. Shares partition category-balanced squared rigid-body error; torque uses the existing1Å conversion. The translation axis is the declared ligand radial direction, including when decomposing water-body forces. These are descriptive partitions, not new acceptance tests.', '', '| Dry fixture | Metal comparison | Ligand radial error | Ligand transverse norm | Ligand torque/1Å norm |','|---|---|---:|---:|---:|']
    for r in fixtures:
        if r['environment']=='dry':
            ligand=next(b for b in r['bodies'] if b['category']=='ligand')
            lines.append(f"| {r['fixture']} | {r['kind']} | {ligand['radial_error']:.3f} | {ligand['transverse_error_norm']:.3f} | {ligand['torque_over_1A_norm']:.3f} |")
    lines += ['',declaration['interpretation'],'','No reserved validation was opened. Original reference uncertainties and failed development requirements remain.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'status':result['status'],'reconstruction':result['reconstruction'],'aggregates':aggregates,'cost':result['cost']},indent=2))

if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();main(args.out)
