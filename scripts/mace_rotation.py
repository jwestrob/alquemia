"""Rotation attribution on pinned MACE cores; uses the existing manifest runner."""
from __future__ import annotations
import argparse
import copy
from pathlib import Path
import shutil
import json
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import TASK_IDS, EV_TO_KCAL, rotation, write_xyz

SCHEMA = 'alquemia.mace_rotation.v1'
PROTOCOL = 'mace_polar_1m_realspace_rotation_attribution_v1'


def rotate_vectors(values, matrix, scalar_channels=1):
    """Backend spherical components are [y,z,x], after the scalar channels."""
    a = np.asarray(values)
    out = a.copy()
    vectors = a[..., scalar_channels:].reshape(*a.shape[:-1], -1, 3)
    cartesian = vectors[..., [2, 0, 1]] @ np.asarray(matrix).T
    out[..., scalar_channels:] = cartesian[..., [1, 2, 0]].reshape(*a.shape[:-1], -1)
    return out


def prepare(original, blocked, agreement, output):
    from mace_hybrid import dry_run
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    refs = {name: record(path) for name, path in [('original', original), ('blocked', blocked)]}
    sources = {name: read_json(verify(ref)) for name, ref in refs.items()}
    if any(c['status'] != 'complete' and any(c['rows'][t]['status'] != 'computed' for t in TASK_IDS)
           for c in sources.values()):
        raise InvalidArtifact('completed original cores required')
    source_manifests = {name: read_json(verify(c['manifest'])) for name,c in sources.items()}
    for c in sources.values(): dry_run(verify(c['manifest']))
    inputs = out/'inputs'; inputs.mkdir()
    core_tasks = {t['task_id']: t for t in source_manifests['original']['tasks'] if t['kind']=='core'}
    rotated = {}
    for name,t in core_tasks.items():
        atoms = xyz(verify(t['xyz'])); center = np.array(next(a[1:] for a in atoms if a[0] in ('La','Ca')))
        coords = (np.array([a[1:] for a in atoms])-center) @ rotation().T + center
        path = inputs/f'{name}_rotate.xyz'
        write_xyz(path, [(a[0], *p) for a,p in zip(atoms,coords)])
        rotated[name] = record(path)
    for backend in ('original','blocked'):
        dest=out/backend; dest.mkdir(); impl=dest/'implementation'; impl.mkdir()
        m=copy.deepcopy(source_manifests[backend])
        pins={}
        for name in ('mace_hybrid.py','mace_rotation.py','affordable_common.py','mace_realspace_compat.py','mace_blocked.py','mace_local_memory.py'):
            shutil.copyfile(Path(__file__).with_name(name),impl/name); pins[name]=record(impl/name)
        m.update(schema_version=SCHEMA,protocol_id=PROTOCOL,implementation=pins,agreement=record(agreement),
                 rotation_source=refs[backend], frozen_density_source=refs['original'], backend=backend,
                 tasks=[],run_inventory={'distinct_mace_energy_force_calls':4,'new_DFT_endpoints':0})
        m['model']['rotation_density_source']=refs['original']
        for name in TASK_IDS:
            t=copy.deepcopy(core_tasks[name]); t.pop('cache_key')
            t.update(xyz=rotated[name],variant='rotate',rotation_matrix=rotation().tolist())
            t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':pins})
            m['tasks'].append(t)
        write_new(dest/'manifest.json',m); dry_run(dest/'manifest.json')
    return {'status':'prepared','campaign':str(out),'new_core_calls':8}


def validate(m):
    if m['protocol_id'] != PROTOCOL or m['backend'] not in ('original','blocked'):
        raise InvalidArtifact('unknown rotation protocol/backend')
    source=read_json(verify(m['rotation_source'])); sm=read_json(verify(source['manifest']))
    frozen=read_json(verify(m['frozen_density_source']))
    verify(frozen['manifest'])
    model=copy.deepcopy(m['model']); model.pop('rotation_density_source')
    if model != sm['model'] or m['software'] != sm['software']:
        raise InvalidArtifact('rotation changes physical model/software')
    if m['model']['rotation_density_source'] != m['frozen_density_source']:
        raise InvalidArtifact('density provenance mismatch')
    if [t['task_id'] for t in m['tasks']] != TASK_IDS:
        raise InvalidArtifact('rotation must contain exactly the four original cores')
    for t in m['tasks']:
        old=next(x for x in sm['tasks'] if x['task_id']==t['task_id'])
        if any(t[k] != old[k] for k in ('charge','spin_multiplicity','kind','metal','state')):
            raise InvalidArtifact('rotation changes scientific state')
        before=xyz(verify(old['xyz'])); after=xyz(verify(t['xyz']))
        center=np.array(next(a[1:] for a in before if a[0] in ('La','Ca')))
        expected=(np.array([a[1:] for a in before])-center) @ rotation().T + center
        if ([a[0] for a in before] != [a[0] for a in after] or
            not np.allclose(expected,np.array([a[1:] for a in after]),rtol=0,atol=1e-12) or
            not np.array_equal(np.array(t['rotation_matrix']),rotation()) or t['variant']!='rotate'):
            raise InvalidArtifact('unapproved coordinate transformation')
        verify(source['rows'][t['task_id']]['forces'])
        verify(source['rows'][t['task_id']]['density_coefficients'])
        verify(frozen['rows'][t['task_id']]['density_coefficients'])


def probe(calc, m, task, output):
    """Fixed actual density, lab-fixed versus co-rotated auxiliary stencil."""
    import torch
    frozen=read_json(verify(m['frozen_density_source'])); fm=read_json(verify(frozen['manifest']))
    old=next(t for t in fm['tasks'] if t['task_id']==task['task_id'])
    coords=np.array([a[1:] for a in xyz(verify(old['xyz']))])
    rotated=np.array([a[1:] for a in xyz(verify(task['xyz']))])
    density=np.load(verify(frozen['rows'][task['task_id']]['density_coefficients']))
    device=next(calc.models[0].parameters()).device
    tensor=lambda a:torch.as_tensor(a,dtype=torch.float64,device=device)
    r,rr,q,qr=map(tensor,(coords,rotated,density,rotate_vectors(density,rotation())))
    batch=torch.zeros(len(r),dtype=torch.long,device=device)
    energy=calc.models[0].coulomb_energy.realspace_energy
    feature=calc.models[0].electric_potential_descriptor.realspace_features
    arrays={}; answer={'density_source':frozen['rows'][task['task_id']]['density_coefficients'],
                       'energy_offset_A':energy.offset,'feature_offset_A':feature.offset}
    with torch.no_grad():
        for label,module in [('energy',energy),('features',feature)]:
            def evaluate(mod, positions, charges):
                if label=='energy': value=mod(source_feats=charges,positions=positions,batch=batch)
                else: value=mod(source_feats=charges,node_positions=positions,batch=batch)[0]
                return value.detach().cpu().numpy()
            base=evaluate(module,r,q); lab=evaluate(module,rr,qr)
            corot=copy.deepcopy(module)
            for axis in ('x','y','z'):
                vector=getattr(corot,axis)
                vector.copy_(vector @ tensor(rotation()).T)
            moving=evaluate(corot,rr,q)
            expected=base if label=='energy' else rotate_vectors(base,rotation(),feature.num_radial)
            for key,val in [('base',base),('lab_rotated',lab),('lab_expected',expected),('corotated_frame',moving)]:
                arrays[f'{label}_{key}']=val
            row={'lab_covariance_max_abs':float(np.max(np.abs(lab-expected))),
                 'corotated_frame_max_abs':float(np.max(np.abs(moving-base))),
                 'corotated_frame_pass':bool(np.allclose(moving,base,atol=1e-8,rtol=1e-10))}
            if label=='energy':
                row.update(base_eV=float(base.item()),lab_rotation_delta_kcal_mol=float((lab-base).item())*EV_TO_KCAL)
            else:
                row['scalar_lab_max_abs']=float(np.max(np.abs(lab[:,:feature.num_radial]-expected[:,:feature.num_radial])))
                row['vector_lab_max_abs']=float(np.max(np.abs(lab[:,feature.num_radial:]-expected[:,feature.num_radial:])))
            answer[label]=row
    path=Path(output)/'rotation_probe_arrays.npz'; np.savez(path,**arrays)
    answer['arrays']=record(path)
    return answer


def collect_rotation(manifest):
    from mace_hybrid import accepted_attempt
    mp=Path(manifest); m=read_json(mp); source=read_json(verify(m['rotation_source']))
    rows={}; checks=[]; attempts=[]
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:
                verify(r['rotation_probe']['arrays']); valid.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        name=t['task_id']; rows[name]=valid[-1] if valid else {'status':'unavailable'}
        if not valid: continue
        new,old=rows[name],source['rows'][name]
        de=(new['energy_eV']-old['energy_eV'])*EV_TO_KCAL
        df=float(np.max(np.abs(np.load(verify(new['forces']))@rotation()-np.load(verify(old['forces'])))))
        dd=float(np.max(np.abs(np.load(verify(new['density_coefficients']))-rotate_vectors(np.load(verify(old['density_coefficients'])),rotation()))))
        checks.append({'task_id':name,'energy_delta_kcal_mol':de,'force_max_delta_eV_A':df,
                       'density_covariance_max_abs':dd,
                       'pass':abs(de)<=m['tolerances']['energy_kcal_mol'] and df<=m['tolerances']['force_max_eV_A']})
    pair_checks=[]
    for size in (33,36):
        ca,la=f'1h4i_qm{size}_Ca',f'1h4i_qm{size}_La'
        if all(rows[n]['status']=='computed' for n in (ca,la)):
            delta=((rows[ca]['energy_eV']-rows[la]['energy_eV'])-(source['rows'][ca]['energy_eV']-source['rows'][la]['energy_eV']))*EV_TO_KCAL
            pair_checks.append({'partition':f'qm{size}','paired_R_delta_kcal_mol':delta,'pass':abs(delta)<=m['tolerances']['energy_kcal_mol']})
    return {'status':'complete' if len(checks)==4 else 'incomplete','protocol_id':PROTOCOL,
            'manifest':record(mp),'backend':m['backend'],'rows':rows,'checks':checks,'pair_checks':pair_checks,'attempts':attempts,
            'reference':None,'S_kcal_mol':None,'decision':'uncalibrated_rotation_diagnostic'}


def compare(original,blocked):
    from mace_hybrid import dry_run
    for path in (original,blocked): dry_run(path)
    a,b=collect_rotation(original),collect_rotation(blocked)
    am,bm=read_json(original),read_json(blocked)
    tol=bm['kernel_validation_tolerances']; checks=[]
    for name in TASK_IDS:
        ta=next(t for t in am['tasks'] if t['task_id']==name); tb=next(t for t in bm['tasks'] if t['task_id']==name)
        if any(ta[k]!=tb[k] for k in ('xyz','charge','spin_multiplicity','rotation_matrix')):
            raise InvalidArtifact('original/blocked comparison uses different inputs')
        x,y=a['rows'][name],b['rows'][name]
        if any(r['status']!='computed' for r in (x,y)):
            checks.append({'task_id':name,'status':'unavailable','pass':False});continue
        de=y['energy_eV']-x['energy_eV']
        df=float(np.max(np.abs(np.load(verify(y['forces']))-np.load(verify(x['forces'])))))
        dd=float(np.max(np.abs(np.load(verify(y['density_coefficients']))-np.load(verify(x['density_coefficients'])))))
        checks.append({'task_id':name,'energy_delta_eV':de,'force_max_delta_eV_A':df,'density_max_delta':dd,
                       'pass':abs(de)<=tol['energy_eV'] and df<=tol['force_eV_A'] and dd<=tol['density']})
    return {'status':'pass' if all(r['pass'] for r in checks) else 'failed','tolerances':tol,
            'checks':checks,'original':a,'blocked':b}


def report_rotation(collection, output):
    c=read_json(collection); m=read_json(verify(c['manifest']))
    if c['protocol_id'] != PROTOCOL:
        raise InvalidArtifact('not a rotation collection')
    lines=['# MACE rotation attribution', '', f"Status: {c['status']}; backend: {c['backend']}.",
           f"Protocol: `{PROTOCOL}`. No calibrated score/class. Baseline unchanged.", '',
           f"Successful core calls: {sum(r['status']=='computed' for r in c['rows'].values())}/{len(m['tasks'])}.",
           'New DFT calls: 0. Unrotated MACE cores reused.', '',
           '| Core | Rotation energy change (kcal/mol) | Maximum force change (eV/Angstrom) | Physical check |',
           '|---|---:|---:|---|']
    for r in c['checks']:
        lines.append(f"| {r['task_id']} | {r['energy_delta_kcal_mol']:.12g} | {r['force_max_delta_eV_A']:.12g} | {r['pass']} |")
    lines.extend(['', '## Paired contrasts', '', '```json', json.dumps(c['pair_checks'],indent=2), '```',
                  '', '## Frozen-density representation checks', ''])
    for name,r in c['rows'].items():
        if r['status']=='computed':
            lines.extend([f'### {name}', '', '```json', json.dumps(r['rotation_probe'],indent=2), '```', ''])
    with Path(output).open('x') as f:f.write('\n'.join(lines)+'\n')
    return {'status':'report_written','report':record(output)}


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare')
    for key in ('original','blocked','agreement','output'):q.add_argument('--'+key,required=True)
    q=sub.add_parser('compare')
    for key in ('original','blocked','output'):q.add_argument('--'+key,required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.original,a.blocked,a.agreement,a.output)
    else:r=compare(a.original,a.blocked);write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k not in ('original','blocked','checks')},indent=2))
    if r['status']=='failed':raise SystemExit(1)

if __name__=='__main__':main()
