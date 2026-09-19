"""Native checks of predetermined finite-water-integral representatives."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from hydration_basin_coordinates import WaterCoordinates
from hydration_water_motion import METHOD, checked_gradient
from hydration_square import endpoint
from mace_hybrid import write_xyz
from water_basin_sampling import SETTINGS, collect as collect_sampling, integral, representatives

PROTOCOL='native_r2scan3c_finite_water_basin_validation_v1'


def prepare(collection,agreement,output,center_id=None):
    source=read_json(collection)
    actual=collect_sampling(verify(source['manifest']));actual_rows={r['task_id']:r for r in actual['rows']}
    selected=[r for r in source['rows'] if center_id is None or r['task_id']==center_id]
    if not selected or any(r['status']!='computed' or r!=actual_rows.get(r['task_id']) for r in selected):
        raise InvalidArtifact('actual completed sampling centers required')
    sm=read_json(verify(source['manifest']));sp=read_json(verify(sm['preparation']))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';shutil.copytree(Path(verify(source['manifest'])).parent/'implementation',impl)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    pins={p.name:record(p) for p in impl.glob('*.py')}
    selection={};definitions=[]
    for row in selected:
        name=row['task_id'];r=read_json(verify(row['result']));c=sp['centers'][name]
        arrays=np.load(verify(r['samples']));chosen=representatives(arrays['q'],arrays['logp'],arrays['changes'])
        for wanted,actual in zip(chosen,r['representatives']):
            if (wanted['kind']!=actual['kind'] or wanted['sample_index']!=actual['sample_index']
                    or wanted['q']!=actual['q'] or wanted['predicted_change_kcal_mol']!=actual['predicted_change_kcal_mol']
                    or not np.isclose(wanted['normalized_importance_weight'],actual['normalized_importance_weight'],rtol=1e-12,atol=1e-15)
                    or not np.isclose(wanted['extent'],actual['extent'],rtol=0,atol=1e-12)):
                raise InvalidArtifact('selection changed after sampling')
        chosen=r['representatives'] # retain the actual compute-node selection/weights
        selection[name]={'center':c,'sampling_result':row['result'],'representatives':chosen}
        for pick in chosen:definitions.append({'task_id':name+'__'+pick['kind'],'center_id':name,**pick})
    pp=out/'preparation.json';write_new(pp,{'protocol_id':PROTOCOL,'source_collection':record(collection),
        'agreement':record(agreement),'selection':selection,'tasks':definitions,'settings':SETTINGS})
    tasks=[]
    for t in definitions:
        c=selection[t['center_id']]['center'];frame=WaterCoordinates(xyz(verify(c['common_xyz'])),c['groups'])
        td=out/'tasks'/t['task_id'];td.mkdir(parents=True);xp=td/'core.xyz';write_xyz(xp,frame.rows(t['q']))
        ip=td/'endpoint.inp';ip.write_text(METHOD+f"\n* xyzfile {c['charge']} 1 core.xyz\n")
        row={**t,'case':c['case'],'metal':c['metal'],'charge':c['charge'],'multiplicity':1,
            'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),
            'engrad_path':str(td/'endpoint.engrad'),'preparation':record(pp)}
        row['cache_key']=cache_key(row|{'protocol_id':PROTOCOL});tasks.append(row)
    old=read_json(verify(sp['validations'][-1]));dft=read_json(verify(read_json(verify(old['DFT_collection']))['manifest']))
    manifest={'protocol_id':PROTOCOL,'tasks':tasks,'preparation':record(pp),'implementation':pins,
        'agreement':record(agreement),'orca':dft['orca'],'execution_policy':dft['execution_policy'],
        'execution_resources':{'mpi_ranks':16,'concurrent_tasks':4},'baseline_changed':False}
    write_new(out/'manifest.json',manifest)
    from affordable_workflow import dry_run
    return dry_run(out/'manifest.json')


def collect(manifest,output):
    m=read_json(manifest);p=read_json(verify(m['preparation']));rows=[]
    for t in m['tasks']:
        try:
            c=p['selection'][t['center_id']]['center'];actual=xyz(verify(t['xyz']))
            frame=WaterCoordinates(xyz(verify(c['common_xyz'])),c['groups']);wanted=frame.rows(t['q'])
            if ([a[0] for a in actual]!=[a[0] for a in wanted] or
                    not np.allclose([a[1:] for a in actual],[a[1:] for a in wanted],atol=1e-12,rtol=0)):
                raise InvalidArtifact('changed physical validation coordinates')
            r=endpoint(record(t['output_path']),record(t['output_path']+'.execution.json'),t['xyz'],t['input'])
            gp=read_json(verify(r['receipt']))['artifacts']['engrad'];g=checked_gradient(gp,r,actual)
            change=(r['energy_hartree']-c['DFT_result']['energy_hartree'])*HA_TO_KCAL
            rows.append({'task_id':t['task_id'],'center_id':t['center_id'],'kind':t['kind'],'status':'complete',
                'result':r,'gradient':gp,'physical_gradient':frame.gradient(t['q'],g).tolist(),
                'actual_change_kcal_mol':change,'predicted_change_kcal_mol':t['predicted_change_kcal_mol'],
                'prediction_error_kcal_mol':t['predicted_change_kcal_mol']-change,
                'normalized_importance_weight':t['normalized_importance_weight'],'extent':t['extent']})
        except (OSError,ValueError,KeyError) as exc:
            rows.append({'task_id':t['task_id'],'status':'unavailable','failure':str(exc)})
    r={'manifest':record(manifest),'rows':rows,'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete'}
    write_new(output,r);return {'status':r['status'],'completed':sum(r['status']=='complete' for r in rows),'tasks':len(rows)}


def report(sampling,validation,output):
    sc=read_json(sampling)
    if sc['status']!='complete':raise InvalidArtifact('actual complete sampling required')
    sm=read_json(verify(sc['manifest']));sp=read_json(verify(sm['preparation']))
    native_rows=[]
    for path in validation:
        dc=read_json(path)
        if dc['status']!='complete':raise InvalidArtifact('actual complete native results required')
        dm=read_json(verify(dc['manifest']));dp=read_json(verify(dm['preparation']))
        prior=read_json(verify(dp['source_collection']))
        if prior['manifest']!=sc['manifest']:raise InvalidArtifact('different sampling/validation experiment')
        for row in dc['rows']:
            for pin in (row['result']['output'],row['result']['receipt'],row['gradient']):verify(pin)
        native_rows.extend(dc['rows'])
    if len({r['task_id'] for r in native_rows})!=len(native_rows):raise InvalidArtifact('duplicated native validation')
    rows=[];by_name={}
    for s in sc['rows']:
        name=s['task_id'];r=read_json(verify(s['result']));c=sp['centers'][name];a=np.load(verify(r['samples']))
        numerical=integral(a['q'],a['logp'],a['seeds'],a['changes'])
        checks=[d for d in native_rows if d['center_id']==name]
        if len(checks)!=4:raise InvalidArtifact('prescribed native validation incomplete')
        weights=np.array([d['normalized_importance_weight'] for d in checks]);weights/=weights.sum()
        errors=np.array([d['prediction_error_kcal_mol'] for d in checks]);mae=float(weights@abs(errors));maximum=float(np.max(abs(errors)))
        native_pass=mae<=SETTINGS['DFT_validation_weighted_absolute_error_kcal_mol'] and maximum<=SETTINGS['DFT_validation_maximum_error_kcal_mol']
        row={'center_id':name,'case':c['case'],'metal':c['metal'],'pattern':c['pattern'],'integrals':numerical,
            'representative_weighted_absolute_error_kcal_mol':mae,'maximum_native_error_kcal_mol':maximum,
            'native_representative_gate_pass':bool(native_pass),'native_checks':checks,
            'native_anchor_stationary':c['native_stationary'],'actual_MACE_evaluations':r['actual_MACE_evaluations'],
            'anchor_energy_hartree':c['DFT_result']['energy_hartree'],'DFT_basin_integral_qualified':False,
            'qualification_limit':'finite representative checks cannot establish full DFT basin or global domain completeness',
            'sampling_result':s['result'],'common_physical_frame_H_permutations':c['common_physical_frame_H_permutations'],
            'common_frame_shape_roundoff_max_A':c['common_frame_shape_roundoff_max_A']}
        rows.append(row);by_name[name]=row
    pairs=[]
    for prefix in ('1F6S__11','6IP9__110'):
        ca=by_name[prefix+'__Ca'];la=by_name[prefix+'__La'];domains=[]
        for c,l in zip(ca['integrals'],la['integrals']):
            domains.append({'translation_limit_A':c['translation_limit_A'],'rotation_limit_radian':c['rotation_limit_radian'],
                'cheap_conditional_Ca_minus_La_Fconfig_kcal_mol':c['conditional_F_config_kcal_mol']-l['conditional_F_config_kcal_mol'],
                'both_numerical_gates_pass':c['numerical_gate_pass'] and l['numerical_gate_pass']})
        pairs.append({'case':prefix,'domains':domains,'both_native_representative_gates_pass':ca['native_representative_gate_pass'] and la['native_representative_gate_pass'],
            'occupancy_score':None,'absolute_decision':None})
    result={'status':'complete','sampling':record(sampling),'validation':[record(v) for v in validation],'settings':SETTINGS,
        'centers':rows,'pairs':pairs,'baseline_changed':False,'occupancy_probabilities':None,'absolute_entropy':None,
        'DFT_reweighted_integral':None,'evidence_use':'consumed_single_biological_comparison',
        'implementation':record(__file__)}
    out=Path(output);out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    lines=['# Finite water basins: actual sampling and native checks','',
        'These are conditional cheap-potential configurational integrals on declared common physical domains. '
        'They are not occupancy, absolute entropy, calibrated discrimination scores, or unbiased DFT reweighting.','',
        '| Center | Native representative weighted error | Maximum error | Native gate | MACE calls |',
        '|---|---:|---:|---|---:|']
    for r in rows:lines.append(f"| {r['center_id']} | {r['representative_weighted_absolute_error_kcal_mol']:.5f} | {r['maximum_native_error_kcal_mol']:.5f} | {r['native_representative_gate_pass']} | {r['actual_MACE_evaluations']} |")
    lines+=['','| Center | Domain Å/rad | Conditional Fconfig | ESS | Scramble difference | Outer10% weight | Numerical gate |',
        '|---|---|---:|---:|---:|---:|---|']
    for r in rows:
        for d in r['integrals']:lines.append(f"| {r['center_id']} | {d['translation_limit_A']}/{d['rotation_limit_radian']} | {d['conditional_F_config_kcal_mol']:.5f} | {d['effective_samples']:.1f} | {d['scramble_difference_kcal_mol']:.5f} | {d['outer_10percent_weight_fraction']:.4f} | {d['numerical_gate_pass']} |")
    lines+=['','| Pair | Domain Å/rad | Cheap Ca−La Fconfig contribution | Both numerical gates | Both native gates |',
        '|---|---|---:|---|---|']
    for p in pairs:
        for d in p['domains']:lines.append(f"| {p['case']} | {d['translation_limit_A']}/{d['rotation_limit_radian']} | {d['cheap_conditional_Ca_minus_La_Fconfig_kcal_mol']:.5f} | {d['both_numerical_gates_pass']} | {p['both_native_representative_gates_pass']} |")
    lines+=['','All energies in kcal/mol. Configurational values use the explicit common measure in the agreement; '
        'its additive gauge cancels in these equal-water-count pairs. Kinetic/internal/quantum/nonpolar terms are not supplied. '
        'A failure of broad-domain target agreement does not negate the previously measured local relaxation improvement.']
    (out/'TABLES.md').write_text('\n'.join(lines)+'\n')
    return {'status':'complete','native_gate_passes':sum(r['native_representative_gate_pass'] for r in rows),'centers':len(rows)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in [('prepare',('collection','agreement','output')),('collect',('manifest','output')),('report',('sampling','validation','output'))]:
        q=sub.add_parser(op)
        for name in fields:q.add_argument('--'+name,required=True,**({'action':'append'} if op=='report' and name=='validation' else {}))
        if op=='prepare':q.add_argument('--center-id')
    a=vars(p.parse_args());op=a.pop('op');print(json.dumps(globals()[op](**a),indent=2))
