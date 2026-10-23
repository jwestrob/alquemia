#!/usr/bin/env python3
"""Independent saved-only replay of the fixed 78-cell native moment result."""
import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import time

HERE = Path(__file__).resolve().parent
BRANCH = HERE.parent
DEBYE = 4.803206802
THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
NUM = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'
checks = 0

def check(value, message):
    global checks
    if not value:
        raise ValueError(message)
    checks += 1

def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    check(pin(p['path']) == p, 'Changed pin: ' + p['path'])
    return Path(p['path'])

def read(p):
    return json.loads(verify(p).read_text())

def numbers(values):
    v = [float(x.replace('D', 'E').replace('d', 'e')) for x in values]
    check(all(math.isfinite(x) for x in v), 'Nonfinite numeric data')
    return v

def one(pattern, text):
    matches = re.findall(pattern, text)
    check(len(matches) == 1, 'Raw scalar denominator: ' + pattern)
    return matches[0]

def vector_close(a, b, tol=1e-12):
    check(len(a) == len(b) == 3 and max(abs(x-y) for x,y in zip(a,b)) <= tol, 'Vector replay mismatch')

def row_replay(row, frozen, config):
    check(row['status'] == 'complete' and row['input'] == frozen['input'], 'Incomplete or changed row')
    for k in ('xyz', 'key', 'parameters'):
        check(verify(row['runtime_input'][k]).read_bytes() == verify(frozen['input'][k]).read_bytes(), 'Runtime input differs')
    native = row['moments']
    text = verify(native['analyze_output']).read_text()
    check('Tinker' in text and 'Version 26.2' in text, 'Native version absent')
    check(not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b', text), 'Native failure marker')
    energy = numbers([one(r'Total Potential Energy\s*:\s*('+NUM+')', text)])[0]
    check(energy == native['energy_kcal_mol'], 'Energy replay differs')
    check(text.count('Energy Component Breakdown :') == 1, 'Component denominator differs')
    components = {}
    for line in text.split('Energy Component Breakdown :',1)[1].splitlines():
        match = re.fullmatch(r'\s*([A-Za-z][A-Za-z /()\-]+?)\s+('+NUM+r')\s+(\d+)\s*', line)
        if match:
            check(match[1].strip() not in components, 'Duplicate component')
            components[match[1].strip()] = {'energy_kcal_mol': numbers([match[2]])[0], 'interactions': int(match[3])}
    closure = energy-math.fsum(v['energy_kcal_mol'] for v in components.values())
    check(components == native['components'] and abs(closure) <= 1e-6 and closure == native['component_closure_kcal_mol'], 'Component closure/replay failed')
    matches = re.findall(r'Induced Dipoles\s*:\s*Iterations\s+(\d+)\s+RMS Residual\s+(\S+)', text)
    check(len(matches) == len(re.findall(r'Induced Dipoles\s*:\s*Iterations',text)) == 2, 'Induction denominator differs')
    residuals = numbers([v for _,v in matches]); iterations = [int(i) for i,_ in matches]
    check(min(iterations) >= 1 and min(residuals) >= 0 and max(residuals) <= 1e-8, 'Induction gate failed')
    check(native['induction']['solutions'] == 2 and native['induction']['maximum_iterations'] == max(iterations) and native['induction']['maximum_printed_RMS_residual_Debye'] == max(residuals), 'Induction replay differs')
    q = numbers([one(r'Total Electric Charge\s*:\s*(\S+)\s+Electrons', text)])[0]
    mu = numbers(one(r'Dipole X,Y,Z-Components\s*:\s*([^\n]+)', text).split())
    check(len(mu) == 3 and abs(q-frozen['charge']) <= 1e-5 and q == native['charge_e'], 'Moment charge/shape mismatch')
    vector_close(mu,native['total_COM_Debye'],0)
    lines = verify(native['induced_output']).read_text().splitlines()
    sites = [line.split() for line in lines[1:] if line.strip()]
    check(int(lines[0].split()[0]) == len(sites) == frozen['atoms'], 'Induced site count differs')
    check(all(len(v)==5 for v in sites) and [int(v[0]) for v in sites] == list(range(1,len(sites)+1)) and [v[1] for v in sites] == frozen['native_names'], 'Induced site identity differs')
    u = [numbers(v[2:]) for v in sites]
    check(u == native['induced_Debye'], 'Induced array replay differs')
    summed = [math.fsum(v[j] for v in u) for j in range(3)]
    lab = [mu[j] + DEBYE*frozen['charge']*frozen['COM_A'][j] for j in range(3)]
    vector_close(summed,native['induced_sum_Debye'],0); vector_close(lab,native['total_lab_Debye'],0)
    replay_error = None
    if frozen['energy_replay_authority'] is not None:
        old = verify(frozen['energy_replay_authority']).read_text()
        old_energy = numbers([one(r'Total Potential Energy\s*:\s*('+NUM+')',old)])[0]
        check(old_energy == frozen['energy_replay_kcal_mol'], 'Old raw energy lineage differs')
        replay_error = energy-old_energy
        check(abs(replay_error) <= 1e-4, 'Old energy replay failed')
    if frozen['kind'] == 'ion':
        check(abs(energy) <= 1e-8 and max(map(abs,summed)) <= 0.5e-6, 'Ion zero gate failed')
    check(len(row['receipts']) == 1, 'Receipt count differs')
    receipt = row['receipts'][0]
    folder = Path(row['runtime_input']['xyz']['path']).parent
    check(receipt['argv']==[config['analyze']['path'],'system.xyz','EM'] and receipt['cwd']==str(folder) and receipt['returncode']==0 and receipt['output']==native['analyze_output'] and all(receipt['threads'][k]=='1' for k in THREADS), 'Actual invocation differs')
    named = one(r'Induced Dipoles Written To\s*:\s*(\S+)', text)
    check(named in ('system.uind','system.001u') and Path(native['induced_output']['path']) == folder/named and [n for n in ('system.uind','system.001u') if (folder/n).exists()] == [named], 'Native induced filename/denominator differs')
    symbols = {}
    for key in ('parameters','key'):
        for line in verify(frozen['input'][key]).read_text().splitlines():
            a=line.split()
            if a and a[0].lower()=='atom':
                a=shlex.split(line,comments=True); symbols[int(a[1])]=a[3][:3]
    check(frozen['native_names']==[symbols[t] for t in frozen['types']], 'Native katom name mapping differs')
    check(math.isfinite(receipt['elapsed_seconds']) and receipt['elapsed_seconds'] >= 0, 'Invalid elapsed cost')
    return {'energy_kcal_mol':energy,'induced_Debye':summed,'lab_Debye':lab,'maximum_residual_Debye':max(residuals),'maximum_iterations':max(iterations),'energy_replay_error_kcal_mol':replay_error,'native_elapsed_seconds':receipt['elapsed_seconds']}

def main(result_path):
    check(all(os.environ.get(k)=='1' for k in THREADS),'One-thread environment required')
    source_review=read(pin(BRANCH/'independent_review_v1/RECOVERY_REVIEW.json'))
    check(source_review['status']=='passed_recovery_source_and_preparation_review','Source review not passed')
    result=read(pin(result_path));config=result['config'];manifest=read(result['manifest'])
    check(result['status']=='zero_field_model_redistribution_complete' and result['error'] is None and result['global_stop'] is False and result['native_program_invocations']==74 and result['inherited_native_program_invocations']==4 and result['total_native_program_invocations']==78,'Native terminal not complete')
    check(result['implementation']==source_review['implementation']==config['implementation'] and result['manifest']==source_review['manifest'],'Frozen source/manifest mismatch')
    for v in config.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}: verify(v)
    for v in config['physics_sources'].values():verify(v)
    release=read(config['root_release']);prep=read(config['preparation_result'])
    check(config['preparation_result']==source_review['preparation_result'] and prep['manifest']==result['manifest'],'Actual preparation differs')
    fp=hashlib.sha256(json.dumps({k:v for k,v in config.items() if k not in ('root_release','preparation_result')},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    check(release['status']=='approved_zero_field_model_redistribution_recovery' and release['implementation']==result['implementation'] and release['config_fingerprint']==fp==manifest['config_fingerprint'] and release['preparation_result']==config['preparation_result'] and release['maximum_native_programs']==74 and release['maximum_total_native_programs']==78 and release['previous_result']==config['previous_result']==result['previous_result'] and release['workers']==4 and release['finite_field_authorized'] is False and release['fit_authorized'] is False,'Release does not bind actual scope')
    check(result['candidate_result']==manifest['candidate_result']==config['candidate_result'] and result['selected_review']==config['selected_review']==manifest['selected_review'],'Selected authority differs')
    check(set(result['rows'])==set(manifest['rows']) and len(result['rows'])==78,'Full 78 slots required')
    counts=dict(Counter(r['kind'] for r in manifest['rows'].values()))
    check(counts=={'cluster':4,'dimer':34,'water':34,'shell':4,'ion':2},'Inventory differs')
    prior=read(config['previous_result'])
    check(prior['native_program_invocations']==4 and len(manifest['reused_rows'])==4 and prior['status']=='failed_or_unavailable' and prior['global_stop'] is True, 'Original failure denominator differs')
    for cid,row in manifest['reused_rows'].items():
        check(result['rows'][cid]==row and row['inherited_native_authority']['result']==config['previous_result'] and row['inherited_native_authority']['producer_implementation']==prior['implementation'], 'Reused producer differs')
        original=prior['rows'][cid]
        check(original['status']=='failed_or_unavailable' and row['original_failure']==original['error'] and row['runtime_input']==original['runtime_input'] and row['receipts']==original['receipts'] and row['moments']['analyze_output']==original['actual_output'], 'Reused raw authority changed')
    check(set(manifest['pending_ids'])==set(result['rows'])-set(manifest['reused_rows']) and len(manifest['pending_ids'])==74, 'Pending74 identity differs')
    check(all('inherited_native_authority' not in result['rows'][c] for c in manifest['pending_ids']), 'New row claims inherited producer')
    rows={cid:row_replay(r,manifest['rows'][cid],config) for cid,r in result['rows'].items()}
    check(sum(r['energy_replay_error_kcal_mol'] is not None for r in rows.values())==38,'Old energy denominator differs')
    parents={}
    for n in (8,9):
        for metal in ('La','Dy'):
            g=f'aquo{n}_central_{metal}';ion='ion_'+metal;shell=g+'__shell'
            dimers=[g+f'__dimer_{i:02d}' for i in range(1,n+1)]
            waters=[g+f'__water_{i:02d}' for i in range(1,n+1)]
            coeffs={'C':{g:1,ion:-1,shell:-1},'P':{**dict.fromkeys(dimers,1),ion:-n,**dict.fromkeys(waters,-1)},'B':{g:1,**dict.fromkeys(dimers,-1),ion:n-1,shell:-1,**dict.fromkeys(waters,1)}}
            parents[g]={}
            for kind,weights in coeffs.items():
                stored=result['model_cycles']['parents'][g][kind]
                check(weights==stored['terms']==manifest['permanent_cancellation'][g][kind]['terms'],'Native cycle coefficients differ')
                u=[math.fsum(c*rows[x]['induced_Debye'][j] for x,c in weights.items()) for j in range(3)]
                lab=[math.fsum(c*rows[x]['lab_Debye'][j] for x,c in weights.items()) for j in range(3)]
                ub=math.fsum(abs(c)*manifest['rows'][x]['atoms']*0.5e-6 for x,c in weights.items())
                mb=math.fsum(abs(c)*0.0005 for c in weights.values());cb=math.fsum(abs(c)*2e-9 for c in weights.values())
                check(ub==stored['induced_print_bound_Debye_per_component'] and mb==stored['M_print_bound_Debye_per_component'] and cb==stored['COM_arithmetic_allowance_Debye_per_component'],'Printing bounds differ')
                vector_close(u,stored['induced_Debye']);vector_close(lab,stored['corrected_M_Debye'])
                check(stored['consistent'] and all(abs(a-b)<=ub+mb+cb for a,b in zip(u,lab)),'Origin/permanent consistency failed')
                parents[g][kind]={'induced_Debye':u,'corrected_M_Debye':lab,'induced_print_bound_Debye_per_component':ub,'consistency_tolerance_Debye_per_component':ub+mb+cb,'M_minus_induced_Debye':[a-b for a,b in zip(lab,u)]}
            vector_close(parents[g]['C']['induced_Debye'],[a+b for a,b in zip(parents[g]['P']['induced_Debye'],parents[g]['B']['induced_Debye'])])
    signs={'aquo8_central_La':1,'aquo8_central_Dy':-1,'aquo9_central_La':-1,'aquo9_central_Dy':1}
    contrasts={}
    for kind in ('C','P','B'):
        value=[math.fsum(c*parents[g][kind]['induced_Debye'][j] for g,c in signs.items()) for j in range(3)]
        bound=math.fsum(p[kind]['induced_print_bound_Debye_per_component'] for p in parents.values())
        saved=result['model_cycles']['Dy9_minus_La9_minus_Dy8_plus_La8'][kind]
        vector_close(value,saved['induced_Debye']);check(bound==saved['induced_print_bound_Debye_per_component'],'Contrast roundoff differs')
        contrasts[kind]={'induced_Debye':value,'induced_print_bound_Debye_per_component':bound}
    check(result['model_cycles']['all_M_consistency_checks_passed'] and result['model_cycles']['new_targets'] is False and result['model_cycles']['polarizability_measured'] is False,'Interpretation status differs')
    return {'status':'passed_native_actual_review','implementation':pin(__file__),'declaration':pin(HERE/'DECLARATION.json'),'input_result':pin(result_path),'source_review':pin(BRANCH/'independent_review_v1/RECOVERY_REVIEW.json'),'manifest':result['manifest'],'root_release':config['root_release'],'job_id':result['job_id'],'checks_passed':checks,'inventory':counts,'old_energy_replays':38,'maximum_old_energy_error_kcal_mol':max(abs(r['energy_replay_error_kcal_mol']) for r in rows.values() if r['energy_replay_error_kcal_mol'] is not None),'maximum_induction_residual_Debye':max(r['maximum_residual_Debye'] for r in rows.values()),'maximum_induction_iterations':max(r['maximum_iterations'] for r in rows.values()),'parents':parents,'Dy9_minus_La9_minus_Dy8_plus_La8':contrasts,'new74_native_cost':result['cost'],'inherited4_native_cost':prior['cost'],'previous_result':config['previous_result'],'sum_single_thread_native_elapsed_seconds':math.fsum(r['native_elapsed_seconds'] for r in rows.values()),'scheduler_cost':'Root records allocated cost separately; no scheduler command in this review','new_native_calls':0,'limitations':['Printing bounds do not bound induction or physical error.','Fixed-model zero-field dipole redistribution, not polarizability or physical validation.','No QM/model comparison or new acceptance threshold in this review.']}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--result',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    try:
        report=main(a.result)
        report['review_cost']={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'max_local_threads':1}
        (out/'REVIEW.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(json.dumps({'review':pin(out/'REVIEW.json'),'status':report['status'],'checks_passed':checks,'contrasts':report['Dy9_minus_La9_minus_Dy8_plus_La8']}))
    except Exception as exc:
        (out/'FAILURE.json').write_text(json.dumps({'status':'failed_saved_review','implementation':pin(__file__),'error':repr(exc),'checks_passed':checks,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'new_native_calls':0},indent=2)+'\n')
        raise
