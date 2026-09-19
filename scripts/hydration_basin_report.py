"""Report native-checked water response separately from unavailable occupancy."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new, xyz
from hydration_square import endpoint
from hydration_basin_coordinates import WaterCoordinates
from hydration_basin import SETTINGS


def thermal_extent(curvature,temperature):
    from scipy.constants import R
    if not curvature['numerical_pass'] or not curvature['positive_definite']:
        return {'status':'unavailable','reason':'unstable_or_numerically_unsupported_curvature','covariance':None}
    k=np.array(curvature['symmetric_curvature'])
    if np.min(np.linalg.eigvalsh(k))<=0:raise InvalidArtifact('inconsistent positive curvature claim')
    covariance=(R*temperature/4184)*np.linalg.solve(k,np.eye(len(k)));waters=[]
    for i in range(len(k)//6):
        block=covariance[6*i:6*i+6,6*i:6*i+6]
        translation=float(np.sqrt(np.trace(block[:3,:3])));rotation=float(np.sqrt(np.trace(block[3:,3:])))
        waters.append({'translation_RMS_A':translation,'rotation_RMS_radian':rotation,
                       'RMS_within_declared_domain':translation<=SETTINGS['translation_limit_A'] and rotation<=SETTINGS['rotation_limit_radian']})
    return {'status':'conditional_classical_harmonic_model_diagnostic','temperature_K':temperature,
            'covariance':covariance.tolist(),'coordinate_units':'per water: Angstrom x3, radian x3',
            'waters':waters,'all_RMS_within_declared_domain':all(w['RMS_within_declared_domain'] for w in waters),
            'physical_fluctuation_estimate':None,'harmonic_free_energy_kcal_mol':None}


def geometry_audit(center, coordinates):
    atoms=xyz(verify(coordinates));origin=center.get('origin_xyz',center['xyz'])
    original=xyz(verify(origin));frame=WaterCoordinates(original,center['groups'])
    pos=np.array([a[1:] for a in atoms]);old=frame.initial
    if [a[0] for a in atoms]!=frame.symbols or not np.array_equal(pos[frame.fixed],old[frame.fixed]):
        raise InvalidArtifact('changed fixed atoms or atom ordering')
    rows=[]
    for w in frame.groups:
        ids=w['indices'];i=w['oxygen_index'];hydrogen=[j for j in ids if frame.symbols[j]=='H']
        if not np.allclose(np.linalg.norm(pos[ids,None]-pos[None,ids],axis=2),
                           np.linalg.norm(old[ids,None]-old[None,ids],axis=2),atol=2e-12,rtol=0):
            raise InvalidArtifact('water internal geometry changed')
        neighbors=[j for j,z in enumerate(frame.symbols) if j not in ids and j!=0 and z!='H']
        closest=min(neighbors,key=lambda j:np.linalg.norm(pos[i]-pos[j]))
        rows.append({'source':w['source'],'metal_O_initial_A':float(np.linalg.norm(old[i]-old[0])),
                     'metal_O_final_A':float(np.linalg.norm(pos[i]-pos[0])),
                     'O_displacement_A':float(np.linalg.norm(pos[i]-old[i])),
                     'OH_lengths_A':[float(np.linalg.norm(pos[j]-pos[i])) for j in hydrogen],
                     'nearest_other_heavy_atom_index':closest,'nearest_other_heavy_atom_element':frame.symbols[closest],
                     'nearest_other_heavy_atom_distance_A':float(np.linalg.norm(pos[i]-pos[closest]))})
    return {'fixed_coordinates_exact':True,'rigid_water_shapes_preserved':True,'origin_xyz':origin,
            'maximum_total_atom_displacement_A':float(np.max(np.linalg.norm(pos-old,axis=1))),'waters':rows}


def report(stage_a, validation, accounting, water_reference, output):
    out=Path(output).resolve()
    if out.exists():raise InvalidArtifact('report output already exists')
    a=read_json(stage_a)
    if a['status']!='complete':raise InvalidArtifact('complete coupled collection required')
    bulk=read_json(water_reference);bm=read_json(verify(bulk['manifest']))
    temperature=read_json(verify(bm['thermochemistry']))['temperature_K']
    initial={};steps=[];pairs=[];checks=[];pins=[]
    for generation,path in enumerate(validation):
        v=read_json(path);verify(v['implementation']);p=read_json(verify(v['preparation']))
        dc=read_json(verify(v['DFT_collection']));dm=read_json(verify(dc['manifest']))
        if v['status']!='complete' or dc['status']!='complete':raise InvalidArtifact('partial native validation')
        if generation==0 and p.get('proposals_only'):raise InvalidArtifact('initial direction checks required')
        if generation>0 and not p.get('proposals_only'):raise InvalidArtifact('expected successive proposal checks')
        defs={t['task_id']:t for t in dm['tasks']};dr={r['task_id']:r for r in dc['rows']}
        pins.append(record(path));checks.extend(v['direction_checks']);current={}
        for row in v['rows']:
            if row['kind'] not in ('proposal','minimum'):continue
            name=row['center_id'];c=p['selection'][name]['center'];t=defs[row['task_id']];d=dr[row['task_id']]
            r=endpoint(d['result']['output'],d['result']['receipt'],t['xyz'],t['input'])
            if r!=d['result']:raise InvalidArtifact('native result changed')
            if generation==0:initial[name]=c['DFT_result']['energy_hartree']
            if name not in initial:raise InvalidArtifact('new center outside initial four')
            prior=c.get('generation',0)
            if prior!=generation:raise InvalidArtifact('nonconsecutive recenter generation')
            change=(r['energy_hartree']-c['DFT_result']['energy_hartree'])*HA_TO_KCAL
            if abs(change-row['actual_change_kcal_mol'])>1e-8:raise InvalidArtifact('reported native change differs')
            step={**row,'generation':generation,'case':c['case'],'pattern':c['pattern'],'metal':c['metal'],
                  'native_energy_hartree':r['energy_hartree'],'native_result':r,
                  'total_change_from_initial_kcal_mol':(r['energy_hartree']-initial[name])*HA_TO_KCAL,
                  'xyz':t['xyz'],'water_count':c['water_count'],'charge':c['charge'],
                  'geometry':geometry_audit(c,t['xyz']),
                  'harmonic_free_energy_kcal_mol':None,'occupancy_probability':None}
            br=read_json(verify(p['selection'][name]['basin_result']))
            step['thermal_extent']=thermal_extent(br['minimum']['curvature'],temperature)
            steps.append(step);current[name]=step
        groups=sorted({(r['case'],r['pattern']) for r in current.values()})
        for case,pattern in groups:
            metal={r['metal']:r for r in current.values() if (r['case'],r['pattern'])==(case,pattern)}
            if set(metal)!={'Ca','La'}:
                pairs.append({'generation':generation,'case':case,'pattern':pattern,'status':'missing_paired_endpoint',
                              'response_score_contribution_kcal_mol':None});continue
            ca,la=metal['Ca'],metal['La']
            if ca['water_count']!=la['water_count']:raise InvalidArtifact('unbalanced paired water inventory')
            pairs.append({'generation':generation,'case':case,'pattern':pattern,'status':'computed_component',
                'R_hartree':ca['native_energy_hartree']-la['native_energy_hartree'],
                'R_kcal_mol':(ca['native_energy_hartree']-la['native_energy_hartree'])*HA_TO_KCAL,
                'response_score_contribution_kcal_mol':ca['total_change_from_initial_kcal_mol']-la['total_change_from_initial_kcal_mol'],
                'last_step_prediction_contrast_error_kcal_mol':ca['energy_error_kcal_mol']-la['energy_error_kcal_mol'],
                'both_step_energy_checks_pass':ca['energy_pass'] and la['energy_pass'],
                'both_native_stationary':ca['stationary_within_tolerance'] and la['stationary_within_tolerance'],
                'absolute_calibrated_decision':None})
    # A native-checked stationary endpoint need not be recalculated merely
    # because its partner needs another step. Make this exact reuse explicit.
    latest={s['center_id']:s for s in steps};final_pairs=[]
    for case,pattern in sorted({(s['case'],s['pattern']) for s in latest.values()}):
        metal={s['metal']:s for s in latest.values() if (s['case'],s['pattern'])==(case,pattern)}
        if set(metal)!={'Ca','La'}:raise InvalidArtifact('missing initial paired native state')
        ca,la=metal['Ca'],metal['La']
        final_pairs.append({'case':case,'pattern':pattern,'status':'latest_measured_native_geometry_component',
            'source_generations':{m:s['generation'] for m,s in metal.items()},
            'native_results':{m:s['native_result'] for m,s in metal.items()},
            'both_native_stationary':all(s['stationary_within_tolerance'] for s in metal.values()),
            'both_checked_interior_minima':all(s['kind']=='minimum' and s['energy_pass'] and s['stationary_within_tolerance'] for s in metal.values()),
            'response_score_contribution_kcal_mol':ca['total_change_from_initial_kcal_mol']-la['total_change_from_initial_kcal_mol'],
            'R_kcal_mol':(ca['native_energy_hartree']-la['native_energy_hartree'])*HA_TO_KCAL,
            'absolute_calibrated_decision':None})
    missing=['validated stationary coupled basin spectra and their physical extent',
             'bound internal water vibration and its consistent reference',
             'non-electrostatic solvent contribution',
             'competing orientation basins and all occupancy-state corrections']
    result={'status':'complete','stage_a':record(stage_a),'validations':pins,'accounting':record(accounting),
            'water_reference':record(water_reference),'thermochemistry_reference':bm['thermochemistry'],
            'steps':steps,'pairs':pairs,'latest_native_pairs':final_pairs,'direction_checks':checks,
            'thermochemistry':{'status':'unavailable','missing':missing,'occupancy_probabilities':None},
            'baseline_changed':False,'classification_accuracy_claim':None,'evidence_use':'consumed_single_biological_comparison',
            'implementation':record(__file__)}
    out.mkdir(parents=True);write_new(out/'result.json',result)
    lines=['# Coupled water response: actual native checks','',
           'All energy changes are kcal/mol. A negative change lowers the native DFT energy. '
           'These are electronic response components; occupancy remains unavailable.','',
           '| Round | Center | Native change this step | Total change | Prediction error | Energy check | Native gradient max |',
           '|---:|---|---:|---:|---:|---|---:|']
    for s in steps:
        lines.append(f"| {s['generation']} | {s['center_id']} | {s['actual_change_kcal_mol']:.5f} | "
                     f"{s['total_change_from_initial_kcal_mol']:.5f} | {s['energy_error_kcal_mol']:.5f} | "
                     f"{'pass' if s['energy_pass'] else 'fail'} | {s['DFT_gradient_max']:.4f} |")
    lines+=['','Gradient maxima use the physical coordinates: kcal/mol/Å for translations, '
            'kcal/mol/radian for rotations; the reported maximum uses the declared 1 Å/radian numerical scale.','',
            '| Round | Site | Total response contribution to Ca−La | Paired step prediction error |',
            '|---:|---|---:|---:|']
    for p in pairs:
        if p['status']!='computed_component':
            lines.append(f"| {p['generation']} | {p['case']} {p['pattern']} | unavailable | unavailable |")
        else:
            lines.append(f"| {p['generation']} | {p['case']} {p['pattern']} | "
                         f"{p['response_score_contribution_kcal_mol']:.5f} | {p['last_step_prediction_contrast_error_kcal_mol']:.5f} |")
    lines+=['','Larger Ca−La is more La-like within this protocol. The aquo offset cancels '
            'from these response contributions; historical decision bands do not transfer.','',
            '## Latest measured geometries, with explicit endpoint reuse','',
            '| Site | Ca round | La round | Total Ca−La response contribution | Both checked interior minima |',
            '|---|---:|---:|---:|---|']
    for p in final_pairs:
        lines.append(f"| {p['case']} {p['pattern']} | {p['source_generations']['Ca']} | "
                     f"{p['source_generations']['La']} | {p['response_score_contribution_kcal_mol']:.5f} | "
                     f"{'yes' if p['both_checked_interior_minima'] else 'no'} |")
    lines+=['','These use actual native energies. Earlier qualified endpoints are reused '
            'explicitly when only their partner receives another step. There is no baseline substitution.','',
            '| Center | Direction | DFT even energy | Predicted even energy | Error | Frozen tolerance | Result |',
            '|---|---|---:|---:|---:|---:|---|']
    for c in checks:
        lines.append(f"| {c['center_id']} | {c['kind']} | {c['DFT_even_kcal_mol']:.6f} | "
                     f"{c['predicted_even_kcal_mol']:.6f} | {c['even_error']:.6f} | {c['even_tolerance']:.6f} | "
                     f"{'pass' if c['pass'] else 'fail'} |")
    lines+=['','## Conditional thermal extent of the cheap basin model','',
            '| Round | Center | Largest translation RMS (Å) | Largest rotation RMS (rad) | All RMS within local domain |',
            '|---:|---|---:|---:|---|']
    for s in steps:
        e=s['thermal_extent']
        if e['status']=='unavailable':continue
        lines.append(f"| {s['generation']} | {s['center_id']} | "
                     f"{max(w['translation_RMS_A'] for w in e['waters']):.4f} | "
                     f"{max(w['rotation_RMS_radian'] for w in e['waters']):.4f} | "
                     f"{'yes' if e['all_RMS_within_declared_domain'] else 'no'} |")
    lines+=['','The domain remains 0.20 Å / 0.35 rad. These are conditional harmonic-model '
            'diagnostics, not validated physical fluctuations or entropy estimates.','',
            '## Occupancy remains unavailable','']+['- '+s+'.' for s in missing]
    lines+=['','No new biological accuracy claim, default change, threshold fit or PQQ rescore.','']
    (out/'TABLES.md').write_text('\n'.join(lines))
    return {'status':'complete','steps':len(steps),'paired_rows':len(pairs),'output':str(out)}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('stage-a','accounting','water-reference','output'):p.add_argument('--'+name,required=True)
    p.add_argument('--validation',action='append',required=True)
    print(json.dumps(report(**vars(p.parse_args())),indent=2))
