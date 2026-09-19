"""Small research classifiers from pinned, already computed La/Ca artifacts.

No energy backend is imported or launched. This module supplies inventory,
feature extraction, grouped evaluation and explicit saved-model scoring.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import resource
import shutil
import time

import gemmi
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.special import expit

from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz

PROTOCOL = 'site_label_ridge_logistic_archived_features_v1'
EV_KCAL = 23.06054783061903
STRUCTURE = ['contact_count', 'acid_oxygen_fraction', 'denticity_excess_fraction',
             'mean_contact_distance_A', 'local_probe_accessible_fraction']
MACE = ['MACE_response_le6_model_kcal', 'MACE_response_6to12_model_kcal']
ARMS = {'DFT': ['DFT_S_kcal_mol'], 'DFT_structure': ['DFT_S_kcal_mol', *STRUCTURE],
        'DFT_structure_MACE': ['DFT_S_kcal_mol', *STRUCTURE, *MACE]}
SETTINGS = {'contact_radius_A': 3.2, 'metal_radius_A': 1.8, 'probe_radius_A': 1.4,
            'sphere_points': 1024, 'homology_identity_longer': .5,
            'alignment_match': 2., 'alignment_mismatch': -1.,
            'alignment_gap_open': -10., 'alignment_gap_extend': -.5,
            'ridge_lambda': 1., 'logit_boundary': 0., 'MACE_closure_model_kcal': .01}


def load_pin(pin):
    return read_json(verify(pin))


def create_output(path):
    out = Path(path).resolve()
    out.mkdir(parents=True, exist_ok=False)
    return out


def implementation(out):
    destination = out/'implementation.py'
    shutil.copyfile(__file__, destination)
    return record(destination)


def inventory(root, agreement, output):
    root = Path(root).resolve(); out = create_output(output)
    def get(path):
        p = root / path
        sources.append(record(p))
        return read_json(p)
    sources = [record(agreement)]
    reference = get('workspaces/mace_omol_20260917/masked_calibration_v1/reference.json')
    canonical = get('workspaces/mace_omol_20260917/charge_ablation_canonical_report_v1/result.json')
    accuracy = get('workspaces/mace_pqq_utility_20260918/accuracy_v1/result.json')
    labels = {r['case_id']: r for r in reference['scores']}
    scores = {r['case_id']: r for r in accuracy['rows']}
    rows = []; unavailable = []
    for row in canonical['scores']:
        case = row['case_id']; label = labels[case]
        if row['status'] != 'computed':
            unavailable.append({'case_id': case, 'target': 'PQQ_functional_class',
                                'label': label['expected_class'], 'reason': row['preparation_status'],
                                'role': label['evaluation_role']})
            continue
        rows.append({'case_id': case, 'target': 'PQQ_functional_class',
                     'label': label['expected_class'], 'role': label['evaluation_role'],
                     'biological_group': label['sequence_accession_group'],
                     'family': 'PQQ_eight_blade_ADH', 'preparation': row['preparation'],
                     'endpoints': {m: row['endpoints'][m]['bound'] for m in ('Ca', 'La')},
                     'MACE_score': label['two_call_model_kcal'],
                     'DFT_score': scores[case]['DFT_S_kcal_mol'],
                     'DFT_protocol': accuracy['DFT_protocol'], 'MACE_protocol': reference['protocol_id'],
                     'baseline_class': scores[case]['DFT_class'],
                     'evidence_stratum': label['evidence_stratum'],
                     'disconnected': label['disconnected_atom_model_eV']})
    # Keep the direct-site stratum separate. No site labels are inferred from
    # protein-level aequorin or cross-readout domain observations.
    development = get('workspaces/mace_omol_20260917/charge_ablation_report_v1/result.json')
    ggr_reports = {name: get('workspaces/mace_omol_20260917/ggr_masked_'+name.lower()+'_report_v1/result.json')
                   for name in ('2FW0', '2FVY')}
    baseline = get('diagnostics/baseline_benchmark_20260915/RESULT.json')
    alpha_scores = get('workspaces/benchmark_set_20260915/ready_tasks_v4/collection_1199508.json')
    ggr_scores = get('diagnostics/ggr_mechanism_plan_20260915/COMPARISON_AB.json')
    get('diagnostics/benchmark_set_20260915/EVIDENCE.json')
    dft = {'GGR_1GLG': next(r['score']['S_kcal_mol'] for r in baseline['scores']
                           if r['case'] == 'ggr_1glg_GGR' and r['lane'] == 'repaired')}
    for name in ('1F6S', '6IP9'):
        dft['ALPHA_'+name] = next(r['score']['S_kcal_mol'] for r in alpha_scores['rows']
                                  if r['case'] == 'alacta_'+name.lower()+'_v2_strong_site')
    for name in ('2FW0', '2FVY'):
        dft['GGR_'+name] = next(r['score']['S_kcal_mol'] for r in ggr_scores['rows']
                               if r['case'] == 'ggr_'+name.lower()+'_formamide')
    disconnected = reference['scores'][0]['disconnected_atom_model_eV']
    for case in ('GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY', 'ALPHA_1F6S', 'ALPHA_6IP9'):
        if case in development['scores']:
            row = development['scores'][case]
            preparation = row['preparation']
            endpoints = {m: row['endpoints'][m]['bound'] for m in ('Ca', 'La')}
            score = row['R_mask_model_kcal']
        else:
            row = ggr_reports[case.split('_')[1]]
            preparation = row['preparation']
            endpoints = {m: row['rows'][m+'_bound_primary'] for m in ('Ca', 'La')}
            score = row['R_mask_model_kcal']
        group = 'GGR_MglB' if case.startswith('GGR') else 'bovine_alpha_lactalbumin'
        rows.append({'case_id': case, 'target': 'direct_site_affinity_direction',
                     'label': 'Ca' if case.startswith('GGR') else 'La', 'role': 'calibration',
                     'biological_group': group, 'family': group, 'preparation': preparation,
                     'endpoints': endpoints, 'MACE_score': score, 'DFT_score': dft[case],
                     'DFT_protocol': 'generic_peptide_amide_vertical_native_r2scan3c_v3',
                     'MACE_protocol': reference['protocol_id'],
                     'baseline_class': None, 'evidence_stratum': 'direct_direction_condition_qualified',
                     'disconnected': disconnected})
    result = {'protocol': PROTOCOL, 'agreement': record(agreement), 'settings': SETTINGS,
              'arms': ARMS, 'sources': sources, 'rows': rows, 'unavailable': unavailable,
              'exclusions': ['aequorin: site affinity labels unresolved',
                             'parvalbumin: supporting cross-study evidence only',
                             'Khoury domains: La ITC versus Ca folding proxy',
                             'LanM: no matched archived masked-MACE features in this inventory',
                             'other ledger cases: unprepared or unresolved; no new labels inferred'],
              'new_energy_evaluations': 0, 'implementation': implementation(out)}
    write_new(out/'inventory.json', result)
    return {'inventory': record(out/'inventory.json'), 'eligible_rows': len(rows),
            'unavailable': unavailable}


def sequence(atoms):
    aliases = {'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'CYX': 'CYS', 'CYM': 'CYS',
               'ASH': 'ASP', 'GLH': 'GLU', 'LYN': 'LYS'}
    residues = {}
    for a in atoms:
        if a['id'].split('/')[-1] != 'CA' or a['element'] != 'C':
            continue
        res = aliases.get(a.get('resname'), a.get('resname'))
        info = gemmi.find_tabulated_residue(res or '')
        letter = info.one_letter_code
        if not letter or not letter.isupper() or letter == 'X':
            raise InvalidArtifact('unsupported protein residue '+str(res))
        key = a['id'].rsplit('/', 1)[0]
        if key in residues:
            raise InvalidArtifact('duplicate source alpha carbon')
        residues[key] = letter
    if not residues:
        raise InvalidArtifact('protein sequence absent')
    return ''.join(residues.values())


def geometry(preparation):
    atoms = preparation['physical_atoms']
    metals = [i for i,a in enumerate(atoms) if a['kind'] == 'selected_metal']
    if len(metals) != 1:
        raise InvalidArtifact('exactly one selected metal required')
    metal = metals[0]; positions = np.array([a['xyz_A'] for a in atoms], dtype=float)
    if not np.isfinite(positions).all():
        raise InvalidArtifact('nonfinite physical coordinates')
    relative = positions - positions[metal]; distances = np.linalg.norm(relative, axis=1)
    donors = []
    for i,a in enumerate(atoms):
        if i == metal or distances[i] > SETTINGS['contact_radius_A'] + 1e-10:
            continue
        name = a['id'].split('/')[-1]; res = a.get('resname')
        eligible = (a['element'] == 'O' or
                    a['element'] == 'N' and (res in ('HIS','HID','HIE','HIP') and name in ('ND1','NE2') or res == 'PQQ') or
                    a['element'] == 'S' and res in ('CYS','CYX','CYM','MET'))
        if eligible:
            donors.append(i)
    if not donors:
        raise InvalidArtifact('no supported contacts')
    residues = {atoms[i]['id'].rsplit('/',1)[0] for i in donors}
    acid = sum(atoms[i].get('resname') in ('ASP','GLU') and
               atoms[i]['id'].split('/')[-1] in ('OD1','OD2','OE1','OE2') for i in donors)
    # A source-defined local frame makes this finite sphere sampling invariant
    # to rigid rotation/translation. It is not a bulk-solvent connectivity test.
    heavy = [i for i,a in enumerate(atoms) if i != metal and a['element'] != 'H' and distances[i] > 1e-8]
    heavy.sort(key=lambda i: (round(float(distances[i]),8), atoms[i]['id']))
    e1 = relative[heavy[0]] / distances[heavy[0]]
    e2 = None
    for i in heavy[1:]:
        v = relative[i] - np.dot(relative[i],e1)*e1
        if np.linalg.norm(v) > 1e-6:
            e2 = v/np.linalg.norm(v); break
    if e2 is None:
        raise InvalidArtifact('local frame is collinear')
    frame = np.array([e1,e2,np.cross(e1,e2)])
    n = SETTINGS['sphere_points']; k = np.arange(n)
    z = 1-2*(k+.5)/n; phi = k*math.pi*(3-math.sqrt(5))
    unit = np.column_stack((np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z))
    radius = SETTINGS['metal_radius_A']+SETTINGS['probe_radius_A']
    points = radius*(unit @ frame); accessible = np.ones(n,dtype=bool)
    for i,a in enumerate(atoms):
        if i == metal:
            continue
        excluded = float(a['radius_A'])+SETTINGS['probe_radius_A']
        if not math.isfinite(excluded) or excluded <= 0:
            raise InvalidArtifact('invalid source atomic radius')
        if distances[i] < radius+excluded:
            accessible &= np.sum((points-relative[i])**2,axis=1) >= excluded**2
    features = dict(zip(STRUCTURE, [float(len(donors)), acid/len(donors),
                    (len(donors)-len(residues))/len(donors), float(distances[donors].mean()),
                    float(accessible.mean())]))
    detail = {'selected_metal_index': metal,
              'contacts': [{'id':atoms[i]['id'],'resname':atoms[i].get('resname'),
                            'element':atoms[i]['element'],'distance_A':float(distances[i])} for i in donors],
              'exposure_definition':'local_probe_occlusion_not_bulk_solvent_connectivity'}
    return features, detail, distances, metal


def native_response(row, prep, distances, metal):
    values = {}; totals = {}; pins = []
    for m in ('Ca','La'):
        endpoint = row['endpoints'][m]
        if endpoint['status'] != 'computed' or not endpoint['parameter_versions_unchanged']:
            raise InvalidArtifact('incomplete or changed MACE calculation')
        if endpoint['charge_feature_adapter']['id'] != 'omol_zero_charge_embedding_before_joint_projection_v1':
            raise InvalidArtifact('wrong MACE charge feature protocol')
        artifact = endpoint['native_readout']['arrays']; p = verify(artifact); pins.append(artifact)
        with np.load(p, allow_pickle=False) as a:
            value = np.asarray(a['node_energy_eV']) + np.asarray(a['embedding_energy_eV'])
        if value.shape != distances.shape or not np.isfinite(value).all():
            raise InvalidArtifact('MACE readout shape or value differs')
        values[m] = value; totals[m] = endpoint['energy_eV']
        if abs(math.fsum(map(float,value))-totals[m])*EV_KCAL > SETTINGS['MACE_closure_model_kcal']:
            raise InvalidArtifact('native readout does not reproduce total')
        coordinates = xyz(verify(prep['endpoints'][m]['xyz']))
        expected = [(m if i == metal else a['element'], *a['xyz_A']) for i,a in enumerate(prep['physical_atoms'])]
        if coordinates != expected:
            raise InvalidArtifact('paired/source atom order or coordinates differ')
        manifest = load_pin(endpoint['manifest'])
        task = next(t for t in manifest['tasks'] if t['task_id'] == endpoint['task_id'])
        if task['preparation'] != row['preparation'] or task['metal_index'] != metal or task['position'] != 'bound':
            raise InvalidArtifact('readout is not mapped to source preparation')
        if task['xyz'] != prep['endpoints'][m]['xyz']:
            # A byte-identical copied input is allowed; a different atom list is not.
            if task['xyz']['sha256'] != prep['endpoints'][m]['xyz']['sha256']:
                raise InvalidArtifact('endpoint coordinates differ from feature source')
    if prep['endpoints']['La']['charge']-prep['endpoints']['Ca']['charge'] != 1:
        raise InvalidArtifact('paired physical charges differ unexpectedly')
    offset = row['disconnected']['Ca']-row['disconnected']['La']
    response = values['Ca']-values['La']; response[metal] -= offset; response *= EV_KCAL
    total = (totals['Ca']-totals['La']-offset)*EV_KCAL
    closure = math.fsum(map(float,response))-total
    if abs(closure)>.01 or abs(total-row['MACE_score'])>.01:
        raise InvalidArtifact('paired readout/reference closure fails')
    return {MACE[0]:float(response[distances <= 6.+1e-10].sum()),
            MACE[1]:float(response[(distances > 6.+1e-10)&(distances <= 12.+1e-10)].sum())}, {
                'readouts':pins,'full_score_model_kcal':total,'closure_error_model_kcal':closure,
                'response_beyond12_model_kcal':float(response[distances > 12.+1e-10].sum())}


def features(inventory, output):
    started = time.monotonic(); cpu = time.process_time()
    inv = read_json(inventory); verify(inv['agreement']); out = create_output(output)
    if inv['settings'] != SETTINGS or inv['arms'] != ARMS:
        raise InvalidArtifact('declared feature/model settings differ')
    for ref in inv['sources']:verify(ref)
    rows = []; failures = []
    for row in inv['rows']:
        try:
            prep = load_pin(row['preparation']); g, detail, distances, metal = geometry(prep)
            r, native = native_response(row,prep,distances,metal)
            seq = sequence(prep['physical_atoms'])
            rows.append({k:v for k,v in row.items() if k not in ('endpoints','disconnected')} | {
                'features': {'DFT_S_kcal_mol':row['DFT_score'], **g, **r}, 'sequence':seq,
                'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
                'geometry_audit':detail,'native_audit':native,'feature_status':'complete'})
        except (InvalidArtifact,KeyError,ValueError,StopIteration) as error:
            failures.append({'case_id':row['case_id'],'target':row['target'],
                             'reason':str(error),'status':'unsupported'})
    result = {'protocol':PROTOCOL,'inventory':record(inventory),'agreement':inv['agreement'],
              'settings':SETTINGS,'arms':ARMS,'rows':rows,'failures':failures,
              'unavailable':inv['unavailable'],'new_energy_evaluations':0,
              'wall_seconds':time.monotonic()-started,'CPU_seconds':time.process_time()-cpu,
              'max_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'software':{'numpy':np.__version__,'scipy':scipy.__version__,'gemmi':gemmi.__version__},
              'implementation':implementation(out)}
    write_new(out/'features.json',result)
    return {'features':record(out/'features.json'),'complete_rows':len(rows),'failures':failures,
            'wall_seconds':result['wall_seconds']}


def group_rows(rows):
    parent = list(range(len(rows)))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    scoring = gemmi.AlignmentScoring()
    # Gemmi accepts integer scores; multiplying all declared scores by two
    # preserves the selected affine alignment exactly.
    scoring.match = 4; scoring.mismatch = -2; scoring.gapo = -20; scoring.gape = -1
    pairs = []
    for i,a in enumerate(rows):
        for j,b in enumerate(rows[:i]):
            if a['target'] != b['target']:continue
            same = a['biological_group'] == b['biological_group']
            identity = None
            if a['target'] == 'PQQ_functional_class':
                alignment = gemmi.align_string_sequences(list(a['sequence']),list(b['sequence']),[],scoring)
                identity = alignment.match_count/max(len(a['sequence']),len(b['sequence']))
                same |= identity >= SETTINGS['homology_identity_longer']
                pairs.append({'a':a['case_id'],'b':b['case_id'],
                              'identity_relative_longer':identity,'joined':same,
                              'cigar':alignment.cigar_str()})
            if same:parent[root(i)] = root(j)
    groups = {}
    for i,row in enumerate(rows):groups.setdefault(root(i),[]).append(row['case_id'])
    labels = {case:'group:'+min(cases) for cases in groups.values() for case in cases}
    return labels,pairs


def balanced_weights(rows):
    classes = {r['label'] for r in rows}
    if classes != {'Ca','La'}:raise InvalidArtifact('training fold contains only one class')
    counts = Counter((r['label'],r['biological_group']) for r in rows)
    groups = {c:len({r['biological_group'] for r in rows if r['label']==c}) for c in classes}
    return np.array([.5/groups[r['label']]/counts[(r['label'],r['biological_group'])] for r in rows])


def objective(theta, x, y, weights):
    z = x@theta[:-1]+theta[-1]
    residual = weights*(expit(z)-y)
    loss = float(np.dot(weights,np.logaddexp(0,z)-y*z)+.5*SETTINGS['ridge_lambda']*np.dot(theta[:-1],theta[:-1]))
    grad = np.r_[x.T@residual+SETTINGS['ridge_lambda']*theta[:-1],residual.sum()]
    return loss,grad


def fit(rows, names):
    weights = balanced_weights(rows)
    compatibility = {key:rows[0][key] for key in ('target','DFT_protocol','MACE_protocol')}
    if any(any(r[key]!=value for key,value in compatibility.items()) for r in rows):
        raise InvalidArtifact('incompatible scientific target or energy protocols')
    x = np.array([[r['features'][k] for k in names] for r in rows],dtype=float)
    if not np.isfinite(x).all():raise InvalidArtifact('missing/nonfinite classifier feature')
    y = np.array([r['label']=='La' for r in rows],dtype=float)
    mean = weights@x; scale = np.sqrt(weights@((x-mean)**2)); active = scale>1e-12
    safe_scale = np.where(active,scale,1.); standard = (x-mean)/safe_scale; standard[:,~active] = 0
    result = minimize(objective,np.zeros(len(names)+1),args=(standard,y,weights),
                      jac=True,method='L-BFGS-B',options={'gtol':1e-10,'ftol':1e-14,'maxiter':2000})
    loss,gradient = objective(result.x,standard,y,weights)
    if not result.success or np.max(np.abs(gradient))>1e-7:
        raise InvalidArtifact('classifier optimizer did not converge: '+str(result.message))
    return {'features':names,'compatibility':compatibility,'mean':mean.tolist(),'scale':safe_scale.tolist(),'active':active.tolist(),
            'coefficients':result.x[:-1].tolist(),'intercept':float(result.x[-1]),
            'objective':loss,'max_abs_gradient':float(np.max(np.abs(gradient))),
            'iterations':int(result.nit),'training_cases':[r['case_id'] for r in rows],
            'training_biological_groups':sorted({r['biological_group'] for r in rows}),
            'ridge_lambda':SETTINGS['ridge_lambda'],'interpretation':'uncalibrated_classifier_logit'}


def predict(model, row):
    if any(row.get(key)!=value for key,value in model['compatibility'].items()):
        raise InvalidArtifact('incompatible prediction target or energy protocol')
    x = np.array([row['features'][k] for k in model['features']],dtype=float)
    if not np.isfinite(x).all():raise InvalidArtifact('missing/nonfinite prediction feature')
    x = (x-np.array(model['mean']))/np.array(model['scale']); x[~np.array(model['active'])] = 0
    logit = float(np.dot(x,model['coefficients'])+model['intercept'])
    return {'logit':logit,'class':'La' if logit>=0 else 'Ca',
            'correct':(('La' if logit>=0 else 'Ca')==row['label']) if row.get('label') in ('Ca','La') else None}


def evaluate(features, output):
    start = time.monotonic(); cpu = time.process_time()
    data = read_json(features); verify(data['agreement']); out = create_output(output)
    if data['settings'] != SETTINGS:raise InvalidArtifact('settings differ')
    labels,pairs = group_rows(data['rows']); rows = [r | {'validation_group':labels[r['case_id']]} for r in data['rows']]
    results = {}; models = {}
    for target in sorted({r['target'] for r in rows}):
        selected = [r for r in rows if r['target']==target]
        train = [r for r in selected if r['role']=='calibration']; transfer = [r for r in selected if r['role']!='calibration']
        folds = []; predictions = {arm:[] for arm in ARMS}
        for group in sorted({r['validation_group'] for r in train}):
            held = [r for r in train if r['validation_group']==group]
            used = [r for r in train if r['validation_group']!=group]
            fold = {'held_group':group,'held_cases':[r['case_id'] for r in held],
                    'training_cases':[r['case_id'] for r in used],'models':{}}
            for arm,names in ARMS.items():
                try:
                    model = fit(used,names); fold['models'][arm] = model
                    predictions[arm] += [{'case_id':r['case_id'],'label':r['label'],'group':group,
                                          'status':'predicted',**predict(model,r)} for r in held]
                except InvalidArtifact as error:
                    fold['models'][arm] = {'status':'untrainable','reason':str(error)}
                    predictions[arm] += [{'case_id':r['case_id'],'label':r['label'],'group':group,
                                          'status':'unavailable','class':None,'logit':None,'correct':None,
                                          'reason':str(error)} for r in held]
            folds.append(fold)
        arms = {}; models[target] = {}
        for arm,names in ARMS.items():
            model = fit(train,names); models[target][arm] = model
            predicted = predictions[arm]
            arms[arm] = {'out_of_group':predicted,'predicted':sum(r['status']=='predicted' for r in predicted),
                         'correct':sum(r['correct'] is True for r in predicted),'total':len(train),
                         'by_class':{c:{'total':sum(r['label']==c for r in predicted),
                                        'predicted':sum(r['label']==c and r['status']=='predicted' for r in predicted),
                                        'correct':sum(r['label']==c and r['correct'] is True for r in predicted)} for c in ('Ca','La')},
                         'training_replay':[{'case_id':r['case_id'],**predict(model,r)} for r in train],
                         'consumed_transfers':[{'case_id':r['case_id'],'label':r['label'],
                              'same_group_in_training':r['validation_group'] in {t['validation_group'] for t in train},
                              **predict(model,r)} for r in transfer]}
        results[target] = {'cases':len(selected),'training_cases':len(train),
                           'biological_groups':len({r['biological_group'] for r in train}),
                           'validation_groups':len({r['validation_group'] for r in train}),
                           'folds':folds,'arms':arms,'prospectively_blind':False}
    result = {'protocol':PROTOCOL,'features':record(features),'agreement':data['agreement'],
              'settings':SETTINGS,'results':results,'models':models,'group_assignments':labels,
              'sequence_pairs':pairs,'unavailable':data['unavailable'],'feature_failures':data['failures'],
              'baseline_changed':False,'new_energy_evaluations':0,
              'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,
              'implementation':implementation(out)}
    write_new(out/'result.json',result)
    for target,arms in models.items():
        for arm,model in arms.items():
            write_new(out/(target+'__'+arm+'.json'),{'protocol':PROTOCOL,'target':target,'arm':arm,
                     'settings':SETTINGS,'feature_source':record(features),'model':model,
                     'qualification':'research_only_consumed_data; see grouped results',
                     'agreement':data['agreement']})
    lines=['# Archived-feature classifier experiment','',
           'All cases are consumed development. No baseline/reference change or new energy calculation.','',
           '| Target | Model | Correct held out | Predicted / total | Groups held out |',
           '|---|---|---:|---:|---:|']
    for target,r in results.items():
        for arm,a in r['arms'].items():
            lines.append(f"| {target} | {arm} | {a['correct']} | {a['predicted']}/{a['total']} | {r['validation_groups']} |")
    lines += ['', 'Untrainable folds remain unavailable. Training replay is not an accuracy gain. '
               'PQQ testing is within one broad fold family. Reused crystals are marked when their '
               'homology group is in training. The direct-site head has only two biological groups.',
               '',f"Feature failures: {len(data['failures'])}; preparation-unavailable cases: {len(data['unavailable'])}.",
               'Detailed folds, all coefficients, row predictions, group identities and source hashes are in result.json.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'result':record(out/'result.json'),'report':record(out/'REPORT.md'),
            'summary':{t:{a:{k:v for k,v in arm.items() if k in ('correct','predicted','total','by_class')}
                          for a,arm in r['arms'].items()} for t,r in results.items()}}


def score(model, features, case, output):
    saved = read_json(model); data = read_json(features)
    if saved['protocol'] != PROTOCOL or saved['settings'] != data['settings']:
        raise InvalidArtifact('incompatible feature/model protocol')
    row = next(r for r in data['rows'] if r['case_id']==case)
    if row['target'] != saved['target']:raise InvalidArtifact('wrong scientific target')
    result = {'case_id':case,'target':saved['target'],'model':record(model),'features':record(features),
              'qualification':saved['qualification'],
              'out_of_group_evaluation':saved.get('out_of_group_evaluation'),
              **predict(saved['model'],row)}
    write_new(output,result)
    return result


def export(result, output):
    """Export already fitted coefficients with their actual evaluation status."""
    data = read_json(result); out = create_output(output); verify(data['features'])
    artifacts = []
    for target,arms in data['models'].items():
        evaluated = data['results'][target]
        for arm,model in arms.items():
            a = evaluated['arms'][arm]
            if a['predicted'] < a['total']:
                status = 'not_evaluable_missing_labelled_groups'
            elif a['correct'] < evaluated['arms']['DFT']['correct']:
                status = 'held_out_regression_do_not_promote'
            else:
                status = 'research_fidelity_retained_accuracy_gain_not_established'
            p = out/(target+'__'+arm+'.json')
            write_new(p,{'protocol':data['protocol'],'target':target,'arm':arm,
                         'settings':data['settings'],'feature_source':data['features'],
                         'model':model,'qualification':status,'evaluation_source':record(result),
                         'out_of_group_evaluation':{k:a[k] for k in ('predicted','correct','total','by_class')},
                         'validation_groups':evaluated['validation_groups'],
                         'prospectively_blind':False,'production_default':False,
                         'agreement':data['agreement']})
            artifacts.append(record(p))
    write_new(out/'manifest.json',{'source':record(result),'models':artifacts,
                                 'new_fits':0,'implementation':implementation(out)})
    return {'manifest':record(out/'manifest.json'),'models':len(artifacts),'new_fits':0}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='operation',required=True)
    p = sub.add_parser('inventory')
    for name in ('root','agreement','output'):p.add_argument('--'+name,required=True)
    p = sub.add_parser('features')
    for name in ('inventory','output'):p.add_argument('--'+name,required=True)
    p = sub.add_parser('evaluate')
    for name in ('features','output'):p.add_argument('--'+name,required=True)
    p = sub.add_parser('score')
    for name in ('model','features','case','output'):p.add_argument('--'+name,required=True)
    p = sub.add_parser('export')
    for name in ('result','output'):p.add_argument('--'+name,required=True)
    args = vars(parser.parse_args()); operation = args.pop('operation')
    print(json.dumps(globals()[operation](**args),indent=2))
