"""Pinned physical-displacement test of MACE as an inexpensive curvature source."""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact, cache_key, energy, HA_TO_KCAL, read_json, record, verify, write_new, xyz
from affordable_response import read_engrad
from mace_hybrid import EV_TO_KCAL, accepted_attempt, check_atoms
from mace_global_benchmark import snapshot, numerical_parent_gate
from mace_gb import MODEL as GB_MODEL

MACE_SCHEMA = 'alquemia.mace_curvature.v1'
GB_SCHEMA = 'alquemia.mace_curvature_gb.v1'
POLICY = 'exact_archived_ggr_physical_displacements_v1'
TOL = {'total_charge_e': 1e-5, 'vacuum_odd_absolute_kcal_mol': .02,
       'vacuum_odd_relative': .05, 'curvature_even_absolute_kcal_mol': .005,
       'curvature_even_relative': .25, 'anchored_absolute_kcal_mol': .02,
       'anchored_relative_to_even': .25}
UNAVAILABLE = {'response_status': 'response_model_not_validated',
               'relaxation_correction_kcal_mol': None, 'entropy_correction_kcal_mol': None,
               'calibrated_class': None}


def expected_ids():
    return {f'{r}_{label}_{metal}' for r in ('extended', 'connected')
            for label in ('center', 'metal_minus', 'metal_plus', 'peptide_minus', 'peptide_plus')
            for metal in ('La', 'Ca')}


def dft_source(collection, full_receipt_check=False):
    c = read_json(collection); mp = verify(c['manifest']); m = read_json(mp)
    if (m['protocol_id'] != 'ggr_local_sensitivity_tightscf_dev_v1' or
            {t['task_id'] for t in m['tasks']} != expected_ids() or len(m['tasks']) != 20 or
            len(c['comparisons']) != 4 or any(r['status'] != 'pass' for r in c['comparisons'])):
        raise InvalidArtifact('requires all exact executed GGR Stage C tasks and passed physical checks')
    for d in m['directions']:
        if (d['representation'] not in ('extended','connected') or
                d['coordinate'] not in ('metal','peptide') or
                d['amplitude'] != (.02 if d['coordinate']=='metal' else math.pi/180)):
            raise InvalidArtifact('physical direction or amplitude changed')
    if full_receipt_check:
        from ggr_sensitivity import executed
        _, executed_rows = executed(mp)
        if any(any(c['energies'][key].get(field) != value for field,value in row.items())
               for key,row in executed_rows.items()):
            raise InvalidArtifact('DFT collection differs from verified actual execution')
    for t in m['tasks']:
        e = c['energies'][t['task_id']]
        if e['status'] != 'complete' or energy(verify(e['output'])) != e['energy_hartree']:
            raise InvalidArtifact('DFT energy unavailable or changed')
        verify(t['input']); verify(t['xyz'])
        receipt = read_json(verify(e['receipt']))
        for ref in receipt['artifacts'].values(): verify(ref)
        if t['task_type'] == 'analytic_gradient':
            g = c['gradients'][t['task_id']]
            for ref in g['artifacts'].values(): verify(ref)
            raw = read_engrad(verify(g['artifacts']['engrad']))
            if abs(raw['energy_Ha']-e['energy_hartree']) > 1e-8:
                raise InvalidArtifact('DFT gradient energy mismatch')
    return c, m


def prepare_mace(dft, parent, agreement, output):
    c, dm = dft_source(dft, full_receipt_check=True)
    pc = read_json(parent); pm = read_json(verify(pc['manifest']))
    if pc['status'] != 'complete' or pm['checkpoint_label'] not in ('medium','large'):
        raise InvalidArtifact('completed pinned medium/large parent required')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    pins = snapshot(out, (*pm['implementation'], 'mace_curvature.py', 'affordable_response.py'))
    model = copy.deepcopy(pm['model']); model['preparation_policy'] = POLICY
    tasks = []
    for old in dm['tasks']:
        tasks.append({'task_id': old['task_id'], 'kind': 'core', 'variant': 'primary',
                      'xyz': old['xyz'], 'metal': old['metal'], 'charge': old['charge'],
                      'spin_multiplicity': old['multiplicity'], 'representation': old['representation'],
                      'coordinate': old['coordinate'], 'amplitude': old['amplitude'],
                      'state': check_atoms(xyz(verify(old['xyz'])), old['charge'])})
    m = {'schema_version': MACE_SCHEMA, 'protocol_id': f'mace_polar_1{pm["checkpoint_label"][0]}_ggr_curvature_v1',
         'agreement': record(agreement), 'DFT_collection': record(dft), 'DFT_manifest': c['manifest'],
         'parent_collection': record(parent), 'numerical_reference': pm['numerical_reference'],
         'checkpoint_label': pm['checkpoint_label'], 'model': model, 'software': pm['software'],
         'implementation': pins, 'tasks': tasks, 'tolerances': TOL,
         'run_inventory': {'new_MACE_calls': 20, 'new_GB_calls': 0, 'new_DFT_endpoints': 0},
         'evidence_use': 'consumed_method_development', **UNAVAILABLE}
    for t in tasks:
        t['cache_key'] = cache_key({'task':t, 'model':model, 'software':m['software'], 'implementation':pins})
    write_new(out/'manifest.json', m)
    return validate(out/'manifest.json')


def actual_collection(path):
    c = read_json(path); mp = verify(c['manifest']); m = read_json(mp)
    if c['status'] != 'complete': raise InvalidArtifact('actual complete collection required')
    for t in m['tasks']:
        candidates = [accepted_attempt(a,t,mp) for a in (mp.parent/'execution'/t['task_id']).glob('attempt_*')]
        if c['rows'][t['task_id']] not in candidates:
            raise InvalidArtifact('collection row lacks matching successful execution receipt')
    return c, m


def validate(manifest):
    m = read_json(manifest); is_gb = m['schema_version'] == GB_SCHEMA
    if m['schema_version'] not in (MACE_SCHEMA,GB_SCHEMA) or m['tolerances'] != TOL:
        raise InvalidArtifact('unsupported curvature protocol/tolerances')
    c, dm = dft_source(verify(m['DFT_collection']))
    if m['DFT_manifest'] != c['manifest']: raise InvalidArtifact('DFT reference changed')
    sm = read_json(verify(m['software']))
    for ref in [m['agreement'], *m['implementation'].values(), sm['python'], sm['requirements'],
                *read_json(verify(sm['backend_source_inventory']))['files']]: verify(ref)
    pc = read_json(verify(m['parent_collection'])); pm = read_json(verify(pc['manifest']))
    expected = copy.deepcopy(pm['model']); expected['preparation_policy'] = POLICY
    if numerical_parent_gate(m)['status'] != 'pass' or m['checkpoint_label'] != pm['checkpoint_label']:
        raise InvalidArtifact('matching MACE numerical parent required')
    if is_gb:
        reference = read_json(verify(m['solver_validation'])); rm = read_json(verify(reference['manifest']))
        if not reference['numerical_checks_pass'] or m['software'] != rm['software'] or m['model'] != GB_MODEL:
            raise InvalidArtifact('matching GB physics/software validation required')
        source, source_m = actual_collection(verify(m['source_mace_collection']))
        validate(verify(source['manifest']))
        if (source_m['DFT_collection'] != m['DFT_collection'] or source_m['checkpoint_label'] != m['checkpoint_label']):
            raise InvalidArtifact('GB source preparation/model mismatch')
    elif m['model'] != expected or m['software'] != pm['software']:
        raise InvalidArtifact('MACE physical model changed')
    if len(m['tasks']) != 20 or {t['task_id'] for t in m['tasks']} != expected_ids():
        raise InvalidArtifact('curvature task inventory changed')
    old = {t['task_id']:t for t in dm['tasks']}
    for t in m['tasks']:
        original = old[t['task_id']]
        for key in ('xyz','metal','charge','coordinate','representation','amplitude'):
            if t[key] != original[key]: raise InvalidArtifact('physical endpoint differs from frozen DFT')
        if t['spin_multiplicity'] != 1 or t['kind'] != 'core': raise InvalidArtifact('electronic state/scope changed')
        check_atoms(xyz(verify(t['xyz'])), t['charge'])
        if is_gb:
            r = source['rows'][t['task_id']]
            if (t['source_density'] != r['density_coefficients'] or t['source_vacuum_energy_eV'] != r['energy_eV'] or
                    t['solver'] != 'native' or t['platform'] != 'CUDA' or t['solvent_dielectric'] != 78.5):
                raise InvalidArtifact('GB source/solver changed')
        base = {k:v for k,v in t.items() if k != 'cache_key'}
        if t['cache_key'] != cache_key({'task':base,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('scientific cache key mismatch')
    return {'status':'pass','tasks':20,'new_DFT_endpoints':0,'manifest':record(manifest)}


def prepare_gb(collection, solver_validation, agreement, output):
    c, m = actual_collection(collection); validate(verify(c['manifest']))
    ref = read_json(solver_validation); rm = read_json(verify(ref['manifest']))
    if not ref['numerical_checks_pass']: raise InvalidArtifact('solver not validated')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    pins = snapshot(out, m['implementation']); result = copy.deepcopy(m)
    for t in result['tasks']:
        t.pop('cache_key'); r = c['rows'][t['task_id']]
        t.update(solver='native',platform='CUDA',solvent_dielectric=78.5,
                 source_density=r['density_coefficients'],source_vacuum_energy_eV=r['energy_eV'])
    result.update(schema_version=GB_SCHEMA,protocol_id=f'mace_polar_1{m["checkpoint_label"][0]}_ggr_curvature_obc2_v1',
                  agreement=record(agreement),model=GB_MODEL,software=rm['software'],implementation=pins,
                  source_mace_collection=record(collection),solver_validation=record(solver_validation),
                  run_inventory={'new_GB_calls':20,'new_MACE_calls':0,'new_DFT_endpoints':0})
    for t in result['tasks']:
        t['cache_key']=cache_key({'task':t,'model':result['model'],'software':result['software'],'implementation':pins})
    write_new(out/'manifest.json',result)
    return validate(out/'manifest.json')


def collect_curvature(manifest):
    mp=Path(manifest).resolve(); m=read_json(mp); rows={}; attempts=[]
    for t in m['tasks']:
        valid=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp)
            if r is not None:
                if m['schema_version']==GB_SCHEMA: verify(r['serialized_system'])
                valid.append(r)
            attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']]=valid[-1] if valid else {'status':'unavailable','energy_eV':None}
    return {'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete',
            'manifest':record(mp),'protocol_id':m['protocol_id'],'rows':rows,'attempts':attempts,**UNAVAILABLE}


def differences(values, h):
    """Values are [minus, center, plus] in one consistent energy unit."""
    minus,center,plus=map(float,values)
    if not np.isfinite(values).all() or not math.isfinite(h) or h<=0:
        raise InvalidArtifact('nonfinite energies or nonpositive physical amplitude')
    odd=((plus-center)-(minus-center))/2
    even=((plus-center)+(minus-center))/2
    return {'minus_change_kcal_mol':minus-center,'plus_change_kcal_mol':plus-center,
            'odd_kcal_mol':odd,'even_kcal_mol':even,'secant_gradient':odd/h,'secant_curvature':2*even/h**2}


def assess(medium, large, output):
    out=Path(output).resolve(); out.mkdir(parents=True, exist_ok=False); allrows={}; sources={}
    for label,path in (('medium',medium),('large',large)):
        gb,gm=actual_collection(path); validate(verify(gb['manifest']))
        if gm['checkpoint_label']!=label or gm['schema_version']!=GB_SCHEMA:
            raise InvalidArtifact('checkpoint or solver collection mismatch')
        mc,mm=actual_collection(verify(gm['source_mace_collection']))
        dc,dm=dft_source(verify(gm['DFT_collection']))
        sources[label]={'gb':record(path),'mace':gm['source_mace_collection'],'DFT':gm['DFT_collection']}
        rows=[]
        for direction in dm['directions']:
            rep,coordinate,h=direction['representation'],direction['coordinate'],direction['amplitude']
            terms={}; projections={}
            for metal in ('La','Ca'):
                ids=[f'{rep}_{coordinate}_minus_{metal}',f'{rep}_center_{metal}',f'{rep}_{coordinate}_plus_{metal}']
                # Subtract in the original units before conversion to retain precision.
                def changes(vals,factor):
                    v=np.array(vals,dtype=float);return differences((v-v[1])*factor,h)
                dv=[dc['energies'][i]['energy_hartree'] for i in ids]
                mv=[mc['rows'][i]['energy_eV'] for i in ids]
                bv=[gb['rows'][i]['GB_reaction_kcal_mol'] for i in ids]
                sv=[mc['rows'][i]['energy_components_eV']['interaction_energy'] for i in ids]
                solvent_changes=(np.array(mv)-mv[1])*EV_TO_KCAL+(np.array(bv)-bv[1])
                terms[metal]={'DFT_CPCM':changes(dv,HA_TO_KCAL),'MACE_vacuum':changes(mv,EV_TO_KCAL),
                              'GB_reaction':changes(bv,1.),'MACE_plus_GB':differences(solvent_changes,h),
                              'MACE_short':changes(sv,EV_TO_KCAL)}
                jac=np.array(direction['qm_jacobian'])
                dftg=np.array(dc['gradients'][ids[1]]['gradient_kcal_mol_per_A'])
                force=np.load(verify(mc['rows'][ids[1]]['forces']))
                projections[metal]={'DFT_CPCM':float(np.sum(dftg*jac)),
                                    'MACE_vacuum':-float(np.sum(force*jac))*EV_TO_KCAL}
            terms['R']={key:{name:terms['Ca'][key][name]-terms['La'][key][name] for name in terms['Ca'][key]}
                        for key in terms['La']}
            projections['R']={key:projections['Ca'][key]-projections['La'][key] for key in projections['La']}
            checks={}
            for metal in ('La','Ca','R'):
                d=terms[metal]['DFT_CPCM']; cheap=terms[metal]['MACE_plus_GB']; vac=terms[metal]['MACE_vacuum']
                numerical=vac['odd_kcal_mol']-h*projections[metal]['MACE_vacuum']
                tol=max(TOL['vacuum_odd_absolute_kcal_mol'],TOL['vacuum_odd_relative']*abs(h*projections[metal]['MACE_vacuum']))
                even_error=cheap['even_kcal_mol']-d['even_kcal_mol']
                even_tol=max(TOL['curvature_even_absolute_kcal_mol'],TOL['curvature_even_relative']*abs(d['even_kcal_mol']))
                anchor={sign:sgn*h*projections[metal]['DFT_CPCM']+cheap['even_kcal_mol']-d[sign+'_change_kcal_mol']
                        for sign,sgn in (('minus',-1),('plus',1))}
                anchor_tol=max(TOL['anchored_absolute_kcal_mol'],TOL['anchored_relative_to_even']*abs(d['even_kcal_mol']))
                checks[metal]={'vacuum_gradient_residual_kcal_mol':numerical,'vacuum_gradient_tolerance_kcal_mol':tol,
                               'vacuum_gradient_pass':abs(numerical)<=tol,'curvature_even_error_kcal_mol':even_error,
                               'curvature_even_tolerance_kcal_mol':even_tol,'curvature_pass':abs(even_error)<=even_tol,
                               'anchored_prediction_error_kcal_mol':anchor,'anchored_tolerance_kcal_mol':anchor_tol,
                               'anchored_prediction_pass':all(abs(v)<=anchor_tol for v in anchor.values())}
            rows.append({'representation':rep,'coordinate':coordinate,'amplitude':h,
                         'coordinate_units':direction['physical_coordinate']['units'],
                         'energy_units':'kcal_mol','terms':terms,'projected_gradients':projections,'checks':checks})
        allrows[label]=rows
    if sources['medium']['DFT']!=sources['large']['DFT']: raise InvalidArtifact('different DFT source panels')
    gates={label:{name:all(c[name] for r in rows for c in r['checks'].values())
                  for name in ('vacuum_gradient_pass','curvature_pass','anchored_prediction_pass')}
           for label,rows in allrows.items()}
    result={'status':'complete','sources':sources,'implementation':record(__file__),'rows':allrows,'gates':gates,
            'tolerances':TOL,'new_DFT_calls':0,'interpretation':'directional_curvature_approximation_screen_only',**UNAVAILABLE}
    write_new(out/'result.json',result)
    lines=['# MACE local curvature screen','',
           'Exact archived physical GGR displacements; zero new DFT. Medium primary, large sensitivity.',
           'Candidate curvature is MACE vacuum + frozen-monopole OBC-II, with charges reevaluated at each geometry.',
           'No combined analytic gradient, relaxation correction, stable-basin or biological-accuracy claim.','',
           '| Model | Representation | Motion | Metal | DFT even | MACE+GB even | Difference | Curvature gate |',
           '|---|---|---|---|---:|---:|---:|---|']
    for label,rows in allrows.items():
        for r in rows:
            for metal in ('La','Ca','R'):
                d=r['terms'][metal];ch=r['checks'][metal]
                lines.append(f'| {label} | {r["representation"]} | {r["coordinate"]} | {metal} | '
                             f'{d["DFT_CPCM"]["even_kcal_mol"]:.6f} | {d["MACE_plus_GB"]["even_kcal_mol"]:.6f} | '
                             f'{ch["curvature_even_error_kcal_mol"]:.6f} | {ch["curvature_pass"]} |')
    lines.extend(['','Even energies are kcal/mol. Negative curvature is retained, not clamped.',
                  'Peptide-path curvature includes the curved coordinate and is not a Cartesian eigenvalue.','',
                  '```json',json.dumps(gates,indent=2),'```','',
                  'All projections, secant curvatures, anchored residuals and vacuum/solvent/short components are in result.json.',
                  'A passing direction does not validate a coupled scaffold model. response_model_not_validated.'])
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare-mace')
    for k in ('dft','parent','agreement','output'):a.add_argument('--'+k,required=True)
    a=sub.add_parser('prepare-gb')
    for k in ('collection','solver-validation','agreement','output'):a.add_argument('--'+k,required=True)
    a=sub.add_parser('report')
    for k in ('medium','large','output'):a.add_argument('--'+k,required=True)
    d=vars(p.parse_args()); command=d.pop('command')
    result={'prepare-mace':prepare_mace,'prepare-gb':prepare_gb,'report':assess}[command](**d)
    print(json.dumps(result if command!='report' else result['gates'],indent=2))


if __name__=='__main__':main()
