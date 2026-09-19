"""Test native MACE as a proposal model against existing DFT water rotations."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL, BOHR_TO_A, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms

STAGE='hydration_orientation_proposal'
TOL={'resolved_energy_change_kcal_mol':.5,'minimum_descent_fraction':.875,
     'minimum_seed_order_fraction':1.,'minimum_median_torque_cosine':.5,
     'undefined_torque_norm_kcal_mol_radian':1e-6}


def prepare(manifests, software, inventory, agreement, output):
    from hydration_network import optimization_trace, geometry_checks, write_xyz
    import mace_omol as omol
    _,out,m=omol.common(inventory,software,agreement,output,STAGE)
    for name in ('hydration_mace.py','hydration_network.py','hydration_square.py','affordable_peptide.py'):
        target=out/'implementation'/name;shutil.copyfile(Path(__file__).with_name(name),target)
        m['implementation'][name]=record(target)
    snapshots=out/'dft_checkpoints';snapshots.mkdir()
    tasks=[]
    for mp in manifests:
        source=read_json(mp)
        if source['protocol_id']!='native_r2scan3c_cpcm_water_network_orientations_v1':
            raise InvalidArtifact('unsupported DFT source protocol')
        for t in source['tasks']:
            if set(t['pattern'])!={'1'}:raise InvalidArtifact('only original full-water searches allowed')
            # Read once, freeze bytes. Only completed printed gradient blocks are used.
            text=Path(t['output_path']).read_text()
            marker='Norm of the Cartesian gradient'
            stop=text.rfind(marker)
            if stop<0:raise InvalidArtifact('no completed DFT analytic gradient')
            text=text[:text.find('\n',stop)+1]
            op=snapshots/(t['task_id']+'.out');op.write_text(text)
            atoms=xyz(verify(t['xyz']));trace=optimization_trace(text,atoms)
            if len(trace['analytic_gradient_evaluations'])<3:raise InvalidArtifact('third gradient not yet available')
            for index in (0,2):
                point=trace['analytic_gradient_evaluations'][index]
                rows=[(a[0],*p) for a,p in zip(atoms,point['coordinates_A'])]
                frozen_error,water_error=geometry_checks(atoms,rows,t['groups'],t['mobile_indices'],enforce=False)
                tid=t['task_id']+f'__iteration_{index+1}'
                xp=out/(tid+'.xyz');write_xyz(xp,rows,tid)
                dp=out/(tid+'.dft.json');write_new(dp,dict(point,source_snapshot=record(op),
                    source_manifest=record(mp),source_task_id=t['task_id'],iteration=index+1,
                    full_optimization_completed=False,coordinate_print_resolution_A=1e-6,
                    frozen_coordinate_error_A=frozen_error,water_internal_error_A=water_error,
                    exact_constraint_gate_pass=frozen_error<=2e-5 and water_error<=2e-4))
                tasks.append({'task_id':tid,'case_id':t['case'],'metal':t['metal'],'metal_index':0,
                    'seed':t['seed'],'iteration':index+1,'source_task_id':t['task_id'],
                    'kind':'core','variant':'primary','energy_component':omol.COMPONENT,
                    'charge':t['charge'],'spin_multiplicity':1,'xyz':record(xp),
                    'state':check_atoms(rows,t['charge']),'energy_only':False,
                    'dft':record(dp),'water_groups':t['groups']})
    m.update(tasks=tasks,source_manifests=[record(mp) for mp in manifests],tolerances=TOL,
             purpose='water_orientation_proposal_utility_not_affinity',new_quantum_calls=0)
    return omol.seal(out,m)


def validate(manifest):
    import mace_omol as omol
    m=read_json(manifest)
    if m['stage']!=STAGE or m['model']!=omol.model(verify(m['software'])) or m['tolerances']!=TOL:
        raise InvalidArtifact('proposal model or frozen tolerances changed')
    verify(m['agreement'])
    for p in m['implementation'].values():verify(p)
    for p in m['source_manifests']:verify(p)
    expected={(case,metal,seed,iteration) for case in ('1F6S','6IP9') for metal in ('Ca','La')
              for seed in ('source','radial_away') for iteration in (1,3)}
    observed={(t['case_id'],t['metal'],t['seed'],t['iteration']) for t in m['tasks']}
    if observed!=expected or len(m['tasks'])!=16:raise InvalidArtifact('finite16-task water comparison changed')
    for t in m['tasks']:
        d=read_json(verify(t['dft']));verify(d['source_snapshot']);rows=xyz(verify(t['xyz']))
        if check_atoms(rows,t['charge'])!=t['state'] or t['spin_multiplicity']!=1:
            raise InvalidArtifact('chemical state changed')
        if not np.allclose([a[1:] for a in rows],d['coordinates_A'],atol=1e-10,rtol=0):
            raise InvalidArtifact('MACE geometry differs from evaluated DFT geometry')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('proposal cache identity changed')
    return {'status':'pass','tasks':len(m['tasks']),'manifest':record(manifest)}


def rotational_gradient(coords, gradient, groups):
    coords=np.asarray(coords);gradient=np.asarray(gradient)
    return np.array([np.cross(coords[w['hydrogen_indices']]-coords[w['oxygen_index']],
                            gradient[w['hydrogen_indices']]).sum(axis=0)
                     for w in groups if w['role']=='variable'])


def collect(manifest):
    validate(manifest);mp=Path(manifest);m=read_json(mp);rows=[];attempts=[];lookup={}
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            if (a/'receipt.json').exists():attempts.append(record(a/'receipt.json'))
            r=accepted_attempt(a,t,mp)
            if r is not None:valid.append((a,r))
        if not valid:
            rows.append({'task_id':t['task_id'],'status':'unavailable'});continue
        a,r=valid[-1]
        # The shared executor returns the validated result dictionary.
        d=read_json(verify(t['dft']));coords=np.array(d['coordinates_A'])
        forces=np.load(verify(r['forces']))
        if forces.shape!=coords.shape or not np.isfinite(forces).all():raise InvalidArtifact('bad native forces')
        qm=rotational_gradient(coords,np.array(d['gradient_Ha_per_bohr'])*HA_TO_KCAL/BOHR_TO_A,t['water_groups'])
        ml=rotational_gradient(coords,-forces*EV_TO_KCAL,t['water_groups'])
        norms=[float(np.linalg.norm(v)) for v in (qm,ml)]
        cosine=float(np.sum(qm*ml)/np.prod(norms)) if min(norms)>TOL['undefined_torque_norm_kcal_mol_radian'] else None
        row={'task_id':t['task_id'],'status':'computed','result':record(a/'result.json'),
             'dft_energy_hartree':d['energy_hartree'],'mace_energy_eV':r['energy_eV'],
             'dft_rotational_gradient_kcal_mol_radian':qm.tolist(),'mace_rotational_gradient_kcal_mol_radian':ml.tolist(),
             'torque_norms_kcal_mol_radian':norms,'torque_cosine':cosine,
             'evaluation_seconds':r['evaluation_seconds'],'wall_seconds':r['wall_seconds']}
        row['source_geometry_audit']={k:d[k] for k in ('frozen_coordinate_error_A','water_internal_error_A','exact_constraint_gate_pass')}
        rows.append(row);lookup[(t['case_id'],t['metal'],t['seed'],t['iteration'])]=row
    complete=len(lookup)==16;descents=[];seed_orders=[];checks=None
    if complete:
        def comparison(a,b,identity):
            dq=(b['dft_energy_hartree']-a['dft_energy_hartree'])*HA_TO_KCAL
            dm=(b['mace_energy_eV']-a['mace_energy_eV'])*EV_TO_KCAL
            resolved=abs(dq)>=TOL['resolved_energy_change_kcal_mol']
            return dict(identity,dft_change_kcal_mol=dq,mace_change_kcal_mol=dm,
                        error_kcal_mol=dm-dq,resolved=resolved,correct_sign=(dq*dm>0) if resolved else None)
        for case in ('1F6S','6IP9'):
            for metal in ('Ca','La'):
                for seed in ('source','radial_away'):
                    descents.append(comparison(lookup[(case,metal,seed,1)],lookup[(case,metal,seed,3)],
                                               {'case':case,'metal':metal,'seed':seed}))
                seed_orders.append(comparison(lookup[(case,metal,'source',1)],lookup[(case,metal,'radial_away',1)],
                                              {'case':case,'metal':metal,'direction':'radial_minus_source'}))
        def fraction(data):
            resolved=[x for x in data if x['resolved']]
            return sum(x['correct_sign'] for x in resolved)/len(resolved) if resolved else None
        cosines=[r['torque_cosine'] for r in rows if r['torque_cosine'] is not None]
        median=float(np.median(cosines)) if cosines else None
        df,sf=fraction(descents),fraction(seed_orders)
        checks={'descent_fraction':df,'seed_order_fraction':sf,'median_torque_cosine':median,
                'defined_torque_cases':len(cosines),'descent_pass':df is not None and df>=TOL['minimum_descent_fraction'],
                'seed_order_pass':sf is not None and sf>=TOL['minimum_seed_order_fraction'],
                'torque_pass':median is not None and median>=TOL['minimum_median_torque_cosine']}
        checks['proposal_utility_pass']=all(checks[k] for k in ('descent_pass','seed_order_pass','torque_pass'))
    return {'manifest':record(mp),'status':'complete' if complete else 'incomplete','rows':rows,'attempt_receipts':attempts,
            'orientation_changes':descents,'starting_seed_orders':seed_orders,'checks':checks,
            'baseline_changed':False,'new_quantum_calls':0,'occupancy_probabilities':None,'biological_accuracy_gain':None,
            'scope':'native vacuum MACE proposal proxy versus CPCM DFT; fixed composition, one consumed biological group'}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare');q.add_argument('--manifests',nargs='+',required=True)
    for name in ('software','inventory','agreement','output'):q.add_argument('--'+name,required=True)
    q=sub.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('op')
    if op=='prepare':result=prepare(**a)
    else:result=collect(a['manifest']);write_new(a['output'],result)
    print(json.dumps({k:result[k] for k in ('status','tasks','checks') if k in result},indent=2))
