#!/usr/bin/env python3
"""Finite zero-field moment counterpart; unchanged consumed water mechanics."""
import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shlex
import shutil
import threading
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
PROJECT='lanm_selectivity_free_energy_20261021'
W=ROOT/'workspaces'/PROJECT
D=ROOT/'diagnostics'/PROJECT
GIDS=tuple(f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy'))
POINT=[3.919394849221829,.17662704152555941,3.640166322240097,.15457541153280463,.15863887317162664,.202479823314524]
MAX_CALLS=78
DEBYE=4.803206802
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
NUMBER=r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,value):
    p=Path(path);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(t,p)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def fingerprint(cfg):return hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k not in ('root_release','preparation_result')},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def default_config():
    prior=pin(W/'reference_chemistry_v1/matched_pair_decomposition_v1/model/run_v1/execution/RESULT.json');old=read(prior);src=W/'explicit_model_v1/tinker_source/source'
    return {'protocol':'zero_field_model_redistribution_v1','implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),
        'pair_adapter':old['implementation'],'pair_result':prior,'data':old['config']['data'],'panel':old['panel'],
        'candidate_result':old['candidate_result'],'selected_review':old['selected_review'],
        'resolver':old['config']['resolver'],'native_helper':old['config']['native_helper'],'analyze':old['config']['analyze'],
        'induction_helper':pin(D/'explicit_model_v1/collective_response_design_v1/RUN.py'),
        'physics_sources':{n:pin(src/n) for n in ('analyze.f','moments.f','induce.f','kmpole.f','rotpole.f','optinit.f','units.f','katom.f')},
        'maximum_native_programs':MAX_CALLS,'workers':4,'root_release':None,'preparation_result':None}

def checked_config(cfg,threads=False):
    require(cfg['protocol']=='zero_field_model_redistribution_v1' and cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Source/specification differs')
    require(cfg['maximum_native_programs']==MAX_CALLS and cfg['workers']==4,'Finite scope differs')
    if threads:require(all(os.environ.get(k)=='1' for k in THREADS),'One-thread controls required')
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
    for v in cfg['physics_sources'].values():verify(v)
    require(cfg['candidate_result']['sha256']=='5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832' and cfg['selected_review']['sha256']=='fd7ea77167c833c3990d7624ff136042ce03f3ccf9dc4169ba7cff2aa51bb6a2','Frozen selected authority differs')

def render(key,prm):
    forbidden={'external-field','exfld-freq','exchpol','exchange-polar','a-axis','b-axis','c-axis','ewald','solvate','solvation','active','inactive','save-uinduce','save-induced','polar-print','save-force','mass','charge-flux','chgflux'}
    for text in (key,prm):
        for line in text.splitlines():
            a=line.split('#',1)[0].split()
            if a:require(a[0].lower() not in forbidden and not a[0].lower().startswith('restrain'),'Nonbaseline term/reporting/mass override present')
    require(re.findall(r'(?im)^\s*parameters\s+(\S+)',key)==['amoeba09.prm'],'Parameter filename differs')
    return key+'\n# Reporting only: zero-field redistribution counterpart\npolar-print\nsave-uinduce\n'

def multipoles(prm,key):
    records={}
    for text in (prm,key):
        lines=text.splitlines()
        for i,line in enumerate(lines):
            a=line.split()
            if a and a[0].lower()=='multipole' and int(a[1]) in (36,37,500,501):
                a=shlex.split(line,comments=True)
                tails=[shlex.split(s,comments=True) for s in lines[i+1:i+5]]
                require([len(x) for x in tails]==[3,1,2,3],'Multipole tail shape differs')
                records[int(a[1])]={'definition':a,'components':[float(a[-1])]+[float(x) for row in tails for x in row]}
    require(records[36]['definition'][1:4]==['36','-37','-37'] and records[37]['definition'][1:4]==['37','36','37'],'Water frame definition differs')
    return records

def extract(parent,indices,cid):
    require(indices==sorted(set(indices)) and all(0<=i<len(parent) for i in indices),'Invalid atom mask')
    mapping={a:j+1 for j,a in enumerate(indices)};rows=[]
    for j,a in enumerate(indices):
        r=parent[a];require(all(int(b)-1 in mapping for b in r[6:]),'Atom mask cuts a bond/frame')
        rows.append([str(j+1),*r[1:6],*[str(mapping[int(b)-1]) for b in r[6:]]])
    return str(len(rows))+'  zero-field '+cid+'\n'+'\n'.join(' '.join(r) for r in rows)+'\n'

def frame_signature(parent,a,records):
    typ=int(parent[a][5]);n=len(parent)
    if a==0:
        require(typ in (500,501) and parent[a][6:]==[] and records[typ]['components']==[3.]+[0.]*9,'Metal local frame differs');axis=[];kind='None'
    else:
        o=1+3*((a-1)//3);h1=o+1;h2=o+2;require(h2<n,'Incomplete water')
        if a==o:require(typ==36,'Water oxygen type differs');axis=[h1,h2];kind='Bisector'
        else:require(typ==37,'Water H type differs');axis=[o,h2 if a==h1 else h1];kind='Z-then-X'
    return {'type':typ,'permanent_record':records[typ],'frame_kind':kind,'axis_parent_indices':axis,'axis_positions_A':[[float(x) for x in parent[j][2:5]] for j in axis],
        'parent_index':a,'position_A':[float(x) for x in parent[a][2:5]],'bond_parent_indices':[int(x)-1 for x in parent[a][6:]]}

def terms(gid,kind):
    n=int(gid.split('_')[0][4:]);metal=gid.rsplit('_',1)[1];ion='ion_'+metal
    if kind=='C':return {gid:1,ion:-1,gid+'__shell':-1}
    if kind=='P':return {**{gid+f'__dimer_{i:02d}':1 for i in range(1,n+1)},ion:-n,**{gid+f'__water_{i:02d}':-1 for i in range(1,n+1)}}
    require(kind=='B','Unknown cycle')
    return {gid:1,**{gid+f'__dimer_{i:02d}':-1 for i in range(1,n+1)},ion:n-1,gid+'__shell':-1,**{gid+f'__water_{i:02d}':1 for i in range(1,n+1)}}

def cancellation(rows):
    proof={}
    for gid in GIDS:
        original=rows[gid];n=original['atoms'];proof[gid]={}
        for kind in ('C','P','B'):
            coefficients=terms(gid,kind);counts=[0]*n
            for cid,c in coefficients.items():
                row=rows[cid]
                for a,sig in zip(row['parent_atom_indices'],row['frame_signatures']):
                    require(sig==original['frame_signatures'][a],'Permanent frame/atom signature changed')
                    counts[a]+=c
            require(counts==[0]*n,'Permanent atomic coefficient cancellation failed')
            proof[gid][kind]={'terms':coefficients,'permanent_atom_coefficients':counts,'permanent_cancellation_exact':True,
                'induced_print_bound_Debye_per_component':math.fsum(abs(c)*rows[cid]['atoms']*0.5e-6 for cid,c in coefficients.items()),
                'M_print_bound_Debye_per_component':math.fsum(abs(c)*0.0005 for c in coefficients.values()),
                'COM_arithmetic_allowance_Debye_per_component':math.fsum(abs(c)*2e-9 for c in coefficients.values())}
    return proof

def prepare(cfg,out):
    checked_config(cfg,True);pair=module(cfg['pair_adapter'],'redistribution_pair');old=read(cfg['pair_result']);oldmanifest=read(old['manifest'])
    pair.admit_result(old['config'],oldmanifest,old)
    require(old['candidate_result']==cfg['candidate_result'] and old['selected_review']==cfg['selected_review'] and old['panel']==cfg['panel'],'Prior model lineage differs')
    data=read(cfg['data']);require(data['baseline_parameters']==POINT,'Selected parameters changed');resolver=module(cfg['resolver'],'redistribution_resolver');rows={}
    for gid in GIDS:
        n=int(gid.split('_')[0][4:]);metal=gid.rsplit('_',1)[1];source=oldmanifest['saved_cluster_models'][gid]['input'];parent=pair.typed_rows(verify(source['xyz']).read_text())
        key=verify(source['key']).read_text();prm=verify(source['parameters']).read_text();rendered=render(key,prm);p=resolver.parse_parameters(prm,key);mp=multipoles(prm,key)
        require(len(parent)==1+3*n and [[float(x) for x in r[2:5]] for r in parent]==data['rows'][gid]['geometry']['positions_A'],'Parent coordinate mismatch')
        masks=[(gid,'cluster',list(range(len(parent))),oldmanifest['saved_cluster_models'][gid]['energy_kcal_mol'],oldmanifest['saved_cluster_models'][gid]['native']['analyze_output'])]
        for i in range(1,n+1):
            oldrow=old['rows'][f'pbe0__aquo{n}__{metal}__water_{i:02d}'];idx=[0,3*i-2,3*i-1,3*i]
            require(pair.typed_rows(verify(oldrow['input']['xyz']).read_text())==pair.typed_rows(extract(parent,idx,'check')),'Admitted literal dimer changed')
            masks.append((gid+f'__dimer_{i:02d}','dimer',idx,oldrow['energy_kcal_mol'],oldrow['native']['analyze_output']))
            masks.append((gid+f'__water_{i:02d}','water',idx[1:],None,None))
        masks.append((gid+'__shell','shell',list(range(1,len(parent))),None,None))
        if n==8:masks.append(('ion_'+metal,'ion',[0],0.0,None))
        else:
            ion=rows['ion_'+metal];require(ion['frame_signatures']==[frame_signature(parent,0,mp)] and all(source[k]['sha256']==ion['source_input'][k]['sha256'] for k in ('key','parameters')),'Bare-ion reuse identity differs')
        for cid,kind,idx,energy,authority in masks:
            folder=out/'inputs'/cid;folder.mkdir(parents=True);text=extract(parent,idx,cid);(folder/'system.xyz').write_text(text);(folder/'system.key').write_text(rendered);shutil.copyfile(verify(source['parameters']),folder/'amoeba09.prm')
            rr=pair.typed_rows(text);types=[int(r[5]) for r in rr];require(all(p['polarize'][t]['alpha_A3']>0 and math.isfinite(p['polarize'][t]['alpha_A3']) for t in types),'Nonpositive/nonfinite site alpha')
            masses=[p['atom'][t]['mass'] for t in types];require(all(math.isfinite(x) and x>0 for x in masses),'Invalid mass');positions=[[float(x) for x in r[2:5]] for r in rr]
            com=[math.fsum(m*x[j] for m,x in zip(masses,positions))/math.fsum(masses) for j in range(3)]
            frames=[frame_signature(parent,a,mp) for a in idx];require(all(set(s['axis_parent_indices'])<=set(idx) for s in frames),'Local frame removed')
            rows[cid]={'id':cid,'kind':kind,'metal':metal,'parent_geometry_id':gid,'atoms':len(idx),'parent_atom_indices':idx,'source_input':source,
                'input':{k:pin(folder/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))},
                'types':types,'names':[r[1] for r in rr],'masses_amu':masses,'COM_A':com,'charge':3 if 0 in idx else 0,
                'frame_signatures':frames,'energy_replay_kcal_mol':energy,'energy_replay_authority':authority}
    require(len(rows)==MAX_CALLS and sum(r['kind']=='dimer' for r in rows.values())==34 and sum(r['energy_replay_authority'] is not None for r in rows.values())==38,'78 cells/38 saved energies required')
    manifest={'status':'prepared','implementation':cfg['implementation'],'config_fingerprint':fingerprint(cfg),'candidate_result':cfg['candidate_result'],'selected_review':cfg['selected_review'],
        'data':cfg['data'],'pair_result':cfg['pair_result'],'panel':cfg['panel'],'rows':rows,'permanent_cancellation':cancellation(rows),
        'selected_parameters':POINT,'maximum_native_programs':MAX_CALLS,'workers':4,'new_native_calls':0,'field':False,'exchpol':False,'forces':False,'fit':False}
    save(out/'MANIFEST.json',manifest);return manifest

def bound_inputs(cfg,m):
    checked_config(cfg);require(m['implementation']==cfg['implementation'] and m['config_fingerprint']==fingerprint(cfg) and m['candidate_result']==cfg['candidate_result'] and m['pair_result']==cfg['pair_result'],'Prepared provenance differs')
    require(len(m['rows'])==MAX_CALLS and m['selected_parameters']==POINT and m['permanent_cancellation']==cancellation(m['rows']),'Prepared inventory/cancellation differs')
    pair=module(cfg['pair_adapter'],'redistribution_bound_pair');resolver=module(cfg['resolver'],'redistribution_bound_resolver')
    for cid,row in m['rows'].items():
        src=row['source_input'];parent=pair.typed_rows(verify(src['xyz']).read_text());key=verify(src['key']).read_text();prm=verify(src['parameters']).read_text()
        require(verify(row['input']['xyz']).read_text()==extract(parent,row['parent_atom_indices'],cid) and verify(row['input']['key']).read_text()==render(key,prm) and verify(row['input']['parameters']).read_text()==prm,'Prepared physical input changed')
        mp=multipoles(prm,key);require(row['frame_signatures']==[frame_signature(parent,a,mp) for a in row['parent_atom_indices']],'Prepared local frames differ')
        p=resolver.parse_parameters(prm,key);require(row['masses_amu']==[p['atom'][t]['mass'] for t in row['types']],'Prepared masses differ')
    return m

def parse_energy(path):
    """Same native energy grammar as pair.parse_analyze; isolated terms may be absent."""
    text=Path(path).read_text();require('Tinker' in text and 'Version 26.2' in text,'Native version missing')
    require(not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b',text),'Native failure/nonfinite marker')
    values=re.findall(r'Total Potential Energy\s*:\s*('+NUMBER+')',text);require(len(values)==1 and text.count('Energy Component Breakdown :')==1,'Energy denominator differs');energy=float(values[0].replace('D','E'))
    components={}
    for line in text.split('Energy Component Breakdown :',1)[1].splitlines():
        a=re.fullmatch(r'\s*([A-Za-z][A-Za-z /()\-]+?)\s+('+NUMBER+r')\s+(\d+)\s*',line)
        if a:
            name=a[1].strip();require(name not in components,'Duplicate energy component');components[name]={'energy_kcal_mol':float(a[2].replace('D','E')),'interactions':int(a[3])}
    require(math.isfinite(energy) and all(math.isfinite(v['energy_kcal_mol']) for v in components.values()),'Nonfinite energy');closure=energy-math.fsum(v['energy_kcal_mol'] for v in components.values());require(abs(closure)<=1e-6,'Energy component closure failed')
    return {'energy_kcal_mol':energy,'components':components,'component_closure_kcal_mol':closure}

def parse_moments(output,induced,row,induction):
    text=Path(output).read_text();energy=parse_energy(output);convergence=induction(output);require(convergence['solutions']==2,'EM must contain exactly two induction solutions')
    charges=re.findall(r'Total Electric Charge\s*:\s*(\S+)\s+Electrons',text);dipoles=re.findall(r'Dipole X,Y,Z-Components\s*:\s*([^\n]+)',text)
    require(len(charges)==len(dipoles)==1,'Total moment denominator differs');charge=float(charges[0]);mu=[float(x) for x in dipoles[0].split()];require(len(mu)==3 and all(math.isfinite(x) for x in mu+[charge]) and abs(charge-row['charge'])<=1e-5,'Moment charge/nonfinite mismatch')
    lines=Path(induced).read_text().splitlines();require(int(lines[0].split()[0])==row['atoms'],'Induced header count differs');rr=[s.split() for s in lines[1:] if s.strip()]
    require(len(rr)==row['atoms'] and all(len(r)==5 for r in rr) and [int(r[0]) for r in rr]==list(range(1,row['atoms']+1)) and [r[1] for r in rr]==row['names'],'Incomplete/duplicate induced site inventory')
    vectors=[[float(x) for x in r[2:]] for r in rr];require(all(math.isfinite(x) for r in vectors for x in r),'Nonfinite induced component');summed=[math.fsum(v[j] for v in vectors) for j in range(3)];lab=[mu[j]+DEBYE*row['charge']*row['COM_A'][j] for j in range(3)]
    if row['energy_replay_kcal_mol'] is not None:require(abs(energy['energy_kcal_mol']-row['energy_replay_kcal_mol'])<=1e-4,'Saved native energy replay failed')
    if row['kind']=='ion':require(abs(energy['energy_kcal_mol'])<=1e-8 and max(map(abs,summed))<=0.5e-6,'Isolated ion zero-response proof failed')
    return dict(energy,induction=convergence,charge_e=charge,induced_Debye=vectors,induced_sum_Debye=summed,total_COM_Debye=mu,total_lab_Debye=lab,
        induced_print_bound_Debye_per_component=row['atoms']*0.5e-6,energy_replay_difference_kcal_mol=None if row['energy_replay_kcal_mol'] is None else energy['energy_kcal_mol']-row['energy_replay_kcal_mol'],
        analyze_output=pin(output),induced_output=pin(induced))

def cycles(manifest,rows):
    require(set(rows)==set(manifest['rows']) and all(r['status']=='complete' for r in rows.values()),'Incomplete 78-cell moment denominator');parents={}
    for gid in GIDS:
        parents[gid]={}
        for kind,proof in manifest['permanent_cancellation'][gid].items():
            weights=proof['terms'];ind=[math.fsum(c*rows[cid]['moments']['induced_sum_Debye'][j] for cid,c in weights.items()) for j in range(3)];total=[math.fsum(c*rows[cid]['moments']['total_lab_Debye'][j] for cid,c in weights.items()) for j in range(3)]
            delta=[a-b for a,b in zip(total,ind)];tol=proof['induced_print_bound_Debye_per_component']+proof['M_print_bound_Debye_per_component']+proof['COM_arithmetic_allowance_Debye_per_component']
            parents[gid][kind]=dict(proof,induced_Debye=ind,corrected_M_Debye=total,M_minus_induced_Debye=delta,consistency_tolerance_Debye_per_component=tol,consistent=all(abs(x)<=tol for x in delta))
        require(max(abs(parents[gid]['C']['induced_Debye'][j]-parents[gid]['P']['induced_Debye'][j]-parents[gid]['B']['induced_Debye'][j]) for j in range(3))<=1e-12,'C/P/B arithmetic closure failed')
    signs={GIDS[0]:1,GIDS[1]:-1,GIDS[2]:-1,GIDS[3]:1};contrasts={}
    for kind in ('C','P','B'):
        contrasts[kind]={'induced_Debye':[math.fsum(signs[g]*parents[g][kind]['induced_Debye'][j] for g in GIDS) for j in range(3)],
            'induced_print_bound_Debye_per_component':math.fsum(parents[g][kind]['induced_print_bound_Debye_per_component'] for g in GIDS)}
    return {'units':'Debye','parents':parents,'Dy9_minus_La9_minus_Dy8_plus_La8':contrasts,'all_M_consistency_checks_passed':all(v['consistent'] for p in parents.values() for v in p.values()),'new_targets':False,'polarizability_measured':False}

class Counter:
    def __init__(self):self.lock=threading.RLock();self.stop=threading.Event();self.count=0
    def reserve(self):
        with self.lock:
            require(not self.stop.is_set() and self.count<MAX_CALLS,'Global stop or 78-program ceiling');self.count+=1

def release_gate(cfg,prep,runtime=True):
    require(cfg['root_release'] is not None,'Root release required');r=read(cfg['root_release'])
    require(r['status']=='approved_zero_field_model_redistribution' and r['implementation']==cfg['implementation'] and r['config_fingerprint']==fingerprint(cfg) and r['preparation_result']==cfg['preparation_result'],'Release source/preparation binding differs')
    require(r['maximum_native_programs']==MAX_CALLS and r['workers']==4 and r['finite_field_authorized'] is False and r['fit_authorized'] is False,'Release finite scope differs')
    require(prep['status']=='prepared' and prep['new_native_calls']==0,'Preparation unavailable')
    if runtime:require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4,'One four-CPU Slurm task required')

def run(cfg,out):
    start=time.monotonic();cpu=time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);checked_config(cfg,True);prep=read(cfg['preparation_result']);m=bound_inputs(cfg,read(prep['manifest']));release_gate(cfg,prep)
    native=module(cfg['native_helper'],'redistribution_native');helper=module(cfg['induction_helper'],'redistribution_induction');count=Counter()
    result={'status':'running','implementation':cfg['implementation'],'config':cfg,'manifest':prep['manifest'],'candidate_result':cfg['candidate_result'],'selected_review':cfg['selected_review'],
        'job_id':os.environ['SLURM_JOB_ID'],'rows':{cid:{'id':cid,'status':'not_attempted'} for cid in m['rows']},'native_program_invocations':0,'global_stop':False,'error':None,'model_cycles':None}
    def persist():
        with count.lock:result.update(native_program_invocations=count.count,global_stop=count.stop.is_set());save(out/'RESULT.json',result)
    def worker(cid):
        if count.stop.is_set():return
        row=m['rows'][cid];folder=out/'cells'/cid;record={'id':cid,'status':'started','input':row['input'],'receipts':[]}
        try:
            folder.mkdir(parents=True)
            for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['input'][k]),folder/n)
            record['runtime_input']={k:pin(folder/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            with count.lock:count.reserve();save(out/'CALL_COUNT.json',{'native_program_invocations':count.count,'maximum':MAX_CALLS});save(folder/'RESULT.json',record)
            native.invoke([cfg['analyze']['path'],'system.xyz','EM'],folder,folder/'analyze.out',record['receipts'])
            require(sorted(p.name for p in folder.glob('system.*u'))==['system.001u'],'Induced-file denominator differs')
            parsed=parse_moments(folder/'analyze.out',folder/'system.001u',row,helper.induction)
            record.update(status='complete',moments=parsed)
        except Exception as exc:
            count.stop.set();record.update(status='failed_or_unavailable',error=repr(exc))
            with count.lock:result['error']=repr(exc)
        finally:
            with count.lock:
                if (folder/'analyze.out').exists():record['actual_output']=pin(folder/'analyze.out')
                result['rows'][cid]=record
                if folder.exists():save(folder/'RESULT.json',record)
                persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(worker,m['rows']))
        require(not count.stop.is_set() and count.count==MAX_CALLS,'Incomplete native execution');result['model_cycles']=cycles(m,result['rows']);require(result['model_cycles']['all_M_consistency_checks_passed'],'Permanent/origin consistency failed');result['status']='zero_field_model_redistribution_complete'
    except Exception as exc:count.stop.set();result.update(status='failed_or_unavailable',error=repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-start,'parent_CPU_seconds':time.process_time()-cpu,'aggregate_child_user_CPU_seconds':ch.ru_utime-ch0.ru_utime,'aggregate_child_system_CPU_seconds':ch.ru_stime-ch0.ru_stime,'allocated_CPUs':4,'per_receipt_child_CPU_deltas_not_summed':True};persist()
    return result

def admit_result(cfg,manifest,result):
    """Saved-only raw re-admission; no molecular calls or reference comparison."""
    bound_inputs(cfg,manifest);release_gate(cfg,read(cfg['preparation_result']),False);helper=module(cfg['induction_helper'],'redistribution_read_induction')
    require(result['status']=='zero_field_model_redistribution_complete' and result['error'] is None and result['global_stop'] is False and result['native_program_invocations']==MAX_CALLS,'Terminal run not admitted')
    require(result['config']==cfg and result['implementation']==cfg['implementation'] and read(result['manifest'])==manifest and result['candidate_result']==cfg['candidate_result'] and result['selected_review']==cfg['selected_review'],'Runtime lineage differs')
    require(set(result['rows'])==set(manifest['rows']),'Runtime denominator differs')
    for cid,row in result['rows'].items():
        frozen=manifest['rows'][cid];require(row['status']=='complete' and row['input']==frozen['input'],'Native row admission differs')
        for k in ('xyz','key','parameters'):require(verify(row['runtime_input'][k]).read_bytes()==verify(frozen['input'][k]).read_bytes(),'Actual native input differs')
        native=row['moments'];parsed=parse_moments(verify(native['analyze_output']),verify(native['induced_output']),frozen,helper.induction);require(parsed==native,'Raw moment replay differs')
        require(len(row['receipts'])==1,'Native receipt denominator differs');r=row['receipts'][0];folder=Path(row['runtime_input']['xyz']['path']).parent
        require(r['argv']==[cfg['analyze']['path'],'system.xyz','EM'] and r['cwd']==str(folder) and r['returncode']==0 and r['output']==native['analyze_output'] and all(r['threads'][k]=='1' for k in THREADS),'Native invocation receipt differs')
    replay=cycles(manifest,result['rows']);require(replay==result['model_cycles'] and replay['all_M_consistency_checks_passed'],'Saved cycles differ')
    return {'rows':result['rows'],'model_cycles':replay,'candidate_result':cfg['candidate_result'],'new_native_calls':0}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('config','prepare','run'));p.add_argument('--config');p.add_argument('--output',required=True);a=p.parse_args()
    if a.stage=='config':save(a.output,default_config());raise SystemExit(0)
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);cfg=json.loads(Path(a.config).read_text());start=time.monotonic();cpu=time.process_time()
    try:
        if a.stage=='prepare':
            require(cfg['root_release'] is None and cfg['preparation_result'] is None,'Preparation must remain unreleased');prepare(cfg,out)
            r={'status':'prepared','implementation':pin(__file__),'config':pin(a.config),'manifest':pin(out/'MANIFEST.json'),'new_native_calls':0,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu};save(out/'RESULT.json',r)
        else:r=run(cfg,out)
    except Exception as exc:r={'status':'precall_or_preparation_failed','implementation':pin(__file__),'error':repr(exc)};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'output':str(out)}));raise SystemExit(0 if r['status'] in ('prepared','zero_field_model_redistribution_complete') else 1)
