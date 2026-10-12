#!/usr/bin/env python3
"""Aquo stage of the finite exact-Linux runtime qualification, no new solvers."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, copy, ctypes, json, os, resource, shutil, tomllib
import orca_parameter_runtime_qualification as base

W=base.W
TRANSITION=['Cu1','Cu2','Mn2','Fe2','Fe3','Co2','Ni2']
CLOSED=['Ca2','Mg2','Sr2','Zn2']
FIELDS=['kcn','shpoly','lgam','refocc']

def changes(a,b,p=''):
    if type(a)!=type(b) and not isinstance(a,(int,float)):return [p]
    if isinstance(a,dict):
        if a.keys()!=b.keys():return [p]
        return [z for k in a for z in changes(a[k],b[k],p+'/'+k)]
    if isinstance(a,list):return [p] if a!=b else []
    return [p] if a!=b else []

def check_La8():
    rows=[base.read(W/'results'/j['id']/'receipt.json') for j in base.read(W/'manifest_La8.json')['jobs']]
    if not all(r['passed'] and r['comparison_to_archived']['accepted_parity'] for r in rows):raise ValueError('La8 stage has not passed')
    return rows

def prepare():
    old=check_La8();source=W/'source'
    for src,dest in [('metal_extension/runtime/data.json','aquo_data.json'),('metal_extension/runtime/sources/author_gfn2.toml','author_gfn2.toml'),('metal_extension/runtime/state_manifest.json','aquo_state_manifest.json'),('metal_extension/runtime/validate_runtime.py','archived_validate_runtime.py')]:
        shutil.copyfile(base.checked_source(src),source/dest)
    data=base.read(source/'aquo_data.json');states=base.read(source/'aquo_state_manifest.json');author=tomllib.loads((source/'author_gfn2.toml').read_text())
    original=W/'results/Ca_vacuum_default/endpoint.xtb.json';native=base.read(original)
    expected_original=base.read(base.checked_source('metal_extension/runtime/sources/historical_native_export.xtb.json'))
    if native!=expected_original:raise ValueError('Linux full export differs from pinned historical schema')
    jobs=[];permutations={}
    for state in CLOSED+TRANSITION:
        s=states[state];symbol=s['metal'];xyzsource=base.checked_source('metal_extension/runtime/inputs/'+state+'_aquo.xyz')
        rows=[r.split() for r in xyzsource.read_text().splitlines()[2:]]
        elems=[r[0] for r in rows]
        if len(rows)!=int(xyzsource.read_text().splitlines()[0]):raise ValueError('atom count')
        nel=sum(sum(author['element'][e]['refocc']) for e in elems)-s['charge']
        nbf=sum(sum(2*'spdfg'.index(x[-1])+1 for x in native['element'][e]['shells']) for e in elems)
        if author['element'][symbol]['shells']!=native['element'][symbol]['shells']:raise ValueError('declared shell mismatch')
        corrected=copy.deepcopy(native)
        if state in TRANSITION:
            for field in FIELDS:corrected['element'][symbol][field]=author['element'][symbol][field]
            diff=changes(native,corrected)
            if not diff or not set(diff)<=set('/element/'+symbol+'/'+f for f in FIELDS):raise ValueError('out-of-scope parameter change')
            permutations[symbol]={'native_shells':native['element'][symbol]['shells'],'author_shells':author['element'][symbol]['shells'],'changes':{f:{'original':native['element'][symbol][f],'corrected':corrected['element'][symbol][f]} for f in FIELDS},'changed_paths':diff,'all_unrelated_fields_equal':True}
        for medium in ['vacuum','water']:
            refs={r['engine']:r for r in data['runs'] if r['state']==state and r['scope']=='aquo' and r['medium']==medium and r['engine'] in ['xtb','tblite']}
            if set(refs)!=set(['xtb','tblite']) or not all(x['accepted'] for x in refs.values()):raise ValueError('essential upstream accepted state missing')
            if any(r['input_sha256']!=base.sha(xyzsource) for r in refs.values()):raise ValueError('upstream geometry binding differs')
            for variant in (['default'] if state in CLOSED else ['default','roundtrip','corrected_tight']):
                jid=f'{state}_aquo_{medium}_{variant}';directory=W/'inputs_aquo'/jid;directory.mkdir(parents=True,exist_ok=False)
                shutil.copyfile(xyzsource,directory/'context.xyz')
                method=' WriteXTBParam true\n ReadXTBParam false' if variant=='default' else ' XTBParamFile "parameters.xtb.json"\n WriteXTBParam false\n ReadXTBParam true'
                tol='1e-12' if variant=='corrected_tight' else '1e-10'
                body='! Native-GFN2-xTB NoAutostart'+(' ALPB(water)' if medium=='water' else '')+'\n%maxcore 1000\n%method\n'+method+'\nend\n%scf\n MaxIter 500\n TolE '+tol+'\n SmearTemp 300\n UseXTBMixer true\nend\n%pal\n nprocs 1\nend\n* xyzfile '+str(s['charge'])+' '+str(s['multiplicity'])+' context.xyz\n'
                (directory/'endpoint.inp').write_text(body)
                archived=next(r for r in data['runs'] if r['state']==state and r['scope']=='aquo' and r['medium']==medium and r['engine']=='native' and r['variant']==variant)
                if base.sha(directory/'endpoint.inp')!=archived['outputs']['endpoint.inp']['sha256']:raise ValueError('native input differs from actual Mac input '+jid)
                if variant=='roundtrip':shutil.copyfile(original,directory/'parameters.xtb.json')
                elif variant=='corrected_tight':(directory/'parameters.xtb.json').write_text(json.dumps(corrected,indent=2,allow_nan=False)+'\n')
                jobs.append({'id':jid,'stage':'aquo50','state':state,'metal':symbol,'medium':medium,'variant':variant,'input_dir':str(directory.relative_to(W)),'expected':{'charge':s['charge'],'multiplicity':s['multiplicity'],'valence_electrons':nel,'atoms':len(rows),'basis_functions':nbf,'smear_temperature_K':300,'max_iterations':500,'tol_energy_hartree':float(tol)},'files':{p.name:base.pin(p) for p in directory.iterdir()},'source_geometry':base.pin(xyzsource),'archived_native_id':archived['id'],'parameter_source':base.pin(original),'expected_counts_derivation':'electron count from author neutral refocc minus formal charge; basis functions from declared native angular shells','upstream_accepted_ids':{k:v['id'] for k,v in refs.items()}})
    if len(jobs)!=50:raise ValueError('finite scope differs')
    manifest={'schema':'native_aquo_qualification_manifest.v1','created_utc':base.now(),'jobs':jobs,'runtime_binding':base.pin(W/'RUNTIME_BINDING.json'),'helper':base.pin(__file__),'base_helper':base.pin(base.__file__),'archived_runner':base.pin(source/'archived_validate_runtime.py'),'author':base.pin(source/'author_gfn2.toml'),'upstream_data':base.pin(source/'aquo_data.json'),'state_manifest':base.pin(source/'aquo_state_manifest.json'),'source_default_export':base.pin(original),'targeted_parameter_changes':permutations,'policy':'default and exact unchanged roundtrip1e-10; corrected_tight1e-12; all printed convergence gates unchanged','threads':{'workers':4,'numerical_threads_per_process':1,'aggregate_numerical_threads':4},'maximum_native_endpoints':50,'first_scout':'Cu1 six aquo cells','no_new_upstream_calls':True}
    base.write(W/'manifest_aquo50.json',manifest)
    print(json.dumps({'prepared_aquo':50,'source_native_input_hashes_all_match':True,'state_counts_from_supplied_schema':True}))

def launch(stage):
    base.verify_binding();check_La8();libc=ctypes.CDLL(None,use_errno=True)
    if libc.prctl(36,1,0,0,0)!=0:raise OSError('subreaper unavailable')
    mp=W/'manifest_aquo50.json';m=base.read(mp)
    if stage=='Cu1':jobs=[j for j in m['jobs'] if j['state']=='Cu1']
    else:
        if not base.read(W/'SCOUT_Cu1.json')['corrected_and_roundtrip_pass']:raise ValueError('Cu1 scout not admitted')
        jobs=[j for j in m['jobs'] if j['state']!='Cu1']
    usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    start=__import__('time').monotonic();results=[]
    base.status('running',[],stage=stage,runner_pid=os.getpid(),requested_jobs=[j['id'] for j in jobs])
    for group in [[j for j in jobs if j['variant']=='default'],[j for j in jobs if j['variant']!='default']]:
        if not group:continue
        if group[0]['variant']!='default':
            for j in group:
                exported=W/'results'/f'{j["state"]}_aquo_{j["medium"]}_default'/'endpoint.xtb.json'
                if base.sha(exported)!=m['source_default_export']['sha256']:raise ValueError('same-state default full export changed; imports not launched')
        with ThreadPoolExecutor(max_workers=min(4,len(group))) as ex:results.extend(ex.map(lambda j:base.run_job(j,mp),group))
    end=resource.getrusage(resource.RUSAGE_CHILDREN)
    base.write(W/f'EXECUTION_{stage}.json',{'stage':stage,'receipts':[base.pin(W/'results'/r['job']['id']/'receipt.json') for r in results],'completed':len(results),'accepted':sum(r['passed'] for r in results),'elapsed_wall_seconds':__import__('time').monotonic()-start,'summed_endpoint_wall_seconds':sum(r['process']['elapsed_seconds'] for r in results),'reaped_child_user_seconds':end.ru_utime-usage.ru_utime,'reaped_child_system_seconds':end.ru_stime-usage.ru_stime})
    base.status('stage_complete',[],stage=stage,completed=len(results),accepted=sum(r['passed'] for r in results))
    summarize(stage)

def summarize(stage='current'):
    rows=[base.read(p) for p in sorted((W/'results').glob('*/receipt.json'))];by={r['job']['id']:r for r in rows}
    pairs=[]
    for state in TRANSITION:
        for medium in ['vacuum','water']:
            ids={v:f'{state}_aquo_{medium}_{v}' for v in ['default','roundtrip','corrected_tight']}
            if not all(x in by for x in ids.values()):continue
            default=by[ids['default']]['output'];rt=by[ids['roundtrip']]['output'];corr=by[ids['corrected_tight']]
            pairs.append({'state':state,'medium':medium,'roundtrip_vs_default':base.compare(rt,default,1e-8,1e-7),'corrected_vs_xtb':corr['comparison_to_archived_engines']['xtb'],'corrected_vs_tblite':corr['comparison_to_archived_engines']['tblite'],'receipts':{k:base.pin(W/'results'/v/'receipt.json') for k,v in ids.items()}})
    cu=[x for x in pairs if x['state']=='Cu1']
    if len(cu)==2 and not (W/'SCOUT_Cu1.json').exists():base.write(W/'SCOUT_Cu1.json',{'pairs':cu,'corrected_and_roundtrip_pass':all(x[k]['accepted_parity'] for x in cu for k in ['roundtrip_vs_default','corrected_vs_xtb','corrected_vs_tblite'])})
    libs=set();exes=set()
    for r in rows:
        libs.update(r.get('observed_library_paths',[]));exes.update(x['exe'] for x in r.get('observed_processes',[]))
    runtime={'source_binding':base.pin(W/'RUNTIME_BINDING.json'),'observed_files':[base.pin(p) for p in sorted(libs|exes) if Path(p).is_file()],'observation':'actual process executables and mapped libraries sampled during owned native endpoints'}
    base.write(W/'OBSERVED_RUNTIME.json',runtime,replace=True)
    result={'schema':'native_runtime_qualification_results.v1','updated_utc':base.now(),'completed_native_endpoints':len(rows),'accepted_native_endpoints':sum(r['passed'] for r in rows),'failed_or_unaccepted':[r['job']['id'] for r in rows if not r['passed']],
            'La8':{'complete':sum(r['job']['stage']=='La8' for r in rows),'accepted_and_archived_parity':sum(r['job']['stage']=='La8' and r['passed'] and r['comparison_to_archived']['accepted_parity'] for r in rows)},'aquo_pairs':pairs,
            'closed_shell_upstream_comparisons':[{'job':r['job']['id'],'accepted':r['passed'],'comparisons':r['comparison_to_archived_engines']} for r in rows if r['job'].get('state') in CLOSED],
            'endpoint_receipts':[base.pin(W/'results'/r['job']['id']/'receipt.json') for r in rows],'runtime_binding':base.pin(W/'RUNTIME_BINDING.json'),'observed_runtime':base.pin(W/'OBSERVED_RUNTIME.json'),'new_xtb_or_tblite_calls':0,'production_changed':False,'all_endpoint_summed_wall_seconds':sum(r['process']['elapsed_seconds'] for r in rows),'max_observed_threads_per_process':max((r['max_process_threads'] for r in rows),default=0),'numerical_gates':'all original printed SCF density/energy criteria, state, charge sum, process and parameter identity; raw and accepted comparisons remain separate','scope_limit':'finite implementation controls; physical spin energetics and protein affinities unqualified'}
    base.write(W/'RESULTS.json',result,replace=True);base.write(base.REPORT/'RESULTS.json',result,replace=True)
    print(json.dumps({'completed':len(rows),'accepted':result['accepted_native_endpoints'],'pairs':len(pairs),'Cu1_gate':base.read(W/'SCOUT_Cu1.json')['corrected_and_roundtrip_pass'] if (W/'SCOUT_Cu1.json').exists() else None}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('op',choices=['prepare','Cu1','remaining44','summarize']);a=p.parse_args()
    if a.op=='prepare':prepare()
    elif a.op=='summarize':summarize()
    else:launch(a.op)
