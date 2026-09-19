"""Native DFT adjudication of exact rigid-water MACE proposals."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL, BOHR_TO_A, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from affordable_response import read_engrad
from hydration_network import geometry_checks
from hydration_mace import rotational_gradient
from hydration_square import endpoint

PROTOCOL='native_r2scan3c_cpcm_mace_proposed_water_orientations_v1'
METHOD='! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF EnGrad'


def prepare(collection,agreement,output):
    from hydration_proposal_opt import collect
    saved=read_json(collection);mp=verify(saved['manifest']);actual=collect(mp)
    if actual!=saved or not saved['all_optimization_eligible']:
        raise InvalidArtifact('all actual predefined MACE searches must converge before DFT adjudication')
    parent=read_json(mp);source=read_json(verify(parent['source_manifests'][0]));out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    source_tasks={t['task_id']:t for t in parent['tasks']};tasks=[]
    for row in saved['rows']:
        t=source_tasks[row['task_id']];initial=xyz(verify(t['xyz']));proposed=xyz(verify(row['proposed_xyz']))
        frozen,water=geometry_checks(initial,proposed,t['water_groups'],t['mobile_indices'])
        if frozen>1e-9 or water>1e-9:raise InvalidArtifact('MACE proposal lost exact geometry constraints')
        d=out/'tasks'/t['task_id'];d.mkdir(parents=True);xp=d/'core.xyz';xp.write_bytes(verify(row['proposed_xyz']).read_bytes())
        ip=d/'endpoint.inp';ip.write_text(METHOD+f"\n* xyzfile {t['charge']} 1 core.xyz\n")
        task={'task_id':t['task_id'],'case':t['case_id'],'metal':t['metal'],'charge':t['charge'],'multiplicity':1,
              'input':record(ip),'xyz':record(xp),'output_path':str(d/'endpoint.out'),'engrad_path':str(d/'endpoint.engrad'),
              'groups':t['water_groups'],'mobile_indices':t['mobile_indices'],'initial_xyz':t['xyz'],
              'selected_seed':row['selected_seed'],'mace_result':row['result'],'source_task':t}
        task['cache_key']=cache_key(task|{'protocol_id':PROTOCOL});tasks.append(task)
    for case in ('1F6S','6IP9'):
        pair={t['metal']:t for t in tasks if t['case']==case};la,ca=pair['La'],pair['Ca']
        if la['charge']-ca['charge']!=1 or la['groups']!=ca['groups']:raise InvalidArtifact('paired state inventory differs')
        a,b=xyz(verify(la['xyz'])),xyz(verify(ca['xyz']))
        fixed=[i for i in range(1,len(a)) if i not in la['mobile_indices']]
        if any(a[i]!=b[i] for i in fixed):raise InvalidArtifact('paired protein/oxygen geometry differs')
    impl=out/'implementation';impl.mkdir();pins={}
    for name in ('hydration_adjudicate.py','hydration_network.py','hydration_mace.py','hydration_square.py','affordable_peptide.py','mace_hybrid.py'):
        p=impl/name;shutil.copyfile(Path(__file__).with_name(name),p);pins[name]=record(p)
    manifest={'protocol_id':PROTOCOL,'tasks':tasks,'agreement':record(agreement),'orca':source['orca'],
              'execution_policy':source['execution_policy'],'mace_collection':record(collection),
              'initial_checkpoints':parent['proposal_collection'],'implementation':pins,
              'paired_geometry_policy':'same fixed protein/oxygen/water inventory; variable H optimized per metal',
              'new_endpoint_count':len(tasks),'baseline_changed':False}
    write_new(out/'manifest.json',manifest)
    from affordable_workflow import dry_run
    return dry_run(out/'manifest.json')


def collect(manifest,output):
    m=read_json(manifest);rows=[]
    starts=read_json(verify(m['initial_checkpoints']));start_m=read_json(verify(starts['manifest']))
    start_tasks={t['task_id']:t for t in start_m['tasks']}
    for t in m['tasks']:
        try:
            op=Path(t['output_path']);rp=Path(str(op)+'.execution.json')
            result=endpoint(record(op),record(rp),t['xyz'],t['input'])
            rec=read_json(rp);gp=verify(rec['artifacts']['engrad']);gradient=read_engrad(gp)
            atoms=xyz(verify(t['xyz']));coords=np.array([a[1:] for a in atoms])
            import gemmi
            if not np.array_equal(gradient['atomic_numbers'],[gemmi.Element(a[0]).atomic_number for a in atoms]):
                raise InvalidArtifact('native gradient atom ordering differs')
            if abs(gradient['energy_Ha']-result['energy_hartree'])>1e-7 or not np.allclose(gradient['coordinates_bohr']*BOHR_TO_A,coords,atol=2e-6,rtol=0):
                raise InvalidArtifact('native gradient energy/geometry mismatch')
            torque=rotational_gradient(coords,gradient['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A,t['groups'])
            comparisons={}
            for r in starts['rows']:
                s=start_tasks[r['task_id']]
                if s['case_id']==t['case'] and s['metal']==t['metal'] and s['iteration']==1:
                    comparisons[s['seed']]=(result['energy_hartree']-r['dft_energy_hartree'])*HA_TO_KCAL
            if set(comparisons)!={'source','radial_away'}:raise InvalidArtifact('missing both original DFT seed energies')
            row={'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'status':'complete','result':result,
                 'gradient':record(gp),'rotational_gradient_kcal_mol_radian':torque.tolist(),
                 'rotational_gradient_max_kcal_mol_radian':float(np.max(abs(torque))),
                 'DFT_change_from_initial_seed_kcal_mol':comparisons,
                 'lower_DFT_energy_than_both_initial_seeds':all(e<0 for e in comparisons.values()),
                 'selected_MACE_seed':t['selected_seed']}
        except (OSError,ValueError,KeyError) as exc:
            row={'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'status':'unavailable','failure':str(exc)}
        rows.append(row)
    cases=[]
    for case in ('1F6S','6IP9'):
        pair={r['metal']:r for r in rows if r['case']==case}
        if any(r['status']!='complete' for r in pair.values()):continue
        r=pair['Ca']['result']['energy_hartree']-pair['La']['result']['energy_hartree']
        cases.append({'case':case,'R_hartree':r,'orientation_delta_R_from_source_kcal_mol':
            pair['Ca']['DFT_change_from_initial_seed_kcal_mol']['source']-pair['La']['DFT_change_from_initial_seed_kcal_mol']['source']})
    result={'manifest':record(manifest),'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
            'rows':rows,'cases':cases,'baseline_changed':False,'calibrated_score':None,'occupancy_probabilities':None,
            'native_optimization_comparison':'pending actual converged geometry-qualified comparator',
            'all_proposals_lower_both_DFT_starts':len(rows)==4 and all(r.get('lower_DFT_energy_than_both_initial_seeds',False) for r in rows)}
    write_new(output,result);return {'status':result['status'],'cases':cases,'all_proposals_lower_both_DFT_starts':result['all_proposals_lower_both_DFT_starts']}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('prepare')
    for name in ('collection','agreement','output'):q.add_argument('--'+name,required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('op');print(json.dumps(globals()[op](**a),indent=2))
