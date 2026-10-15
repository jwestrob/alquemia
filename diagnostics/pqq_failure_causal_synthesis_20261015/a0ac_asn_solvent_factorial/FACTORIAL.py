#!/usr/bin/env python3
"""Eight new native charge/geometry crosses; no SCC, ML, DFT or new geometry."""
import argparse
import copy
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import socket
import subprocess
import sys
import time
import traceback

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_solvent_factorial'
SOLVENT=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_solvent_v1'
PROTOCOL='A0AC_sample4_saved_geometry_vacuum_charge_native_ALPB_factorial_v1'
GEOMETRIES=('adaptive_La','amide_flip');METALS=('Ca','La');TERMS=('Gelec','Ghb','Gsasa','Gshift')

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def require(value,reason):
    if not value:raise ValueError(reason)
def load(name,p):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(output):
    """Metadata only; both sample4 sources and both metals are mandatory."""
    cp=SOLVENT/'CONFIG.json';rp=SOLVENT/'execution/FINAL_COLLECTION.json';ep=SOLVENT/'execution/EXECUTION.json'
    c=read(cp);result=read(rp);ex=read(ep);require(result['config']==pin(cp)==ex['config'] and ex['status']=='complete' and result['available_new_cells']==20,'Completed parent stage required')
    cases=[]
    for conditioning in METALS:
        source,=[r for r in result['rows'] if r['conditioning']==conditioning and r['sample']==4]
        for metal in METALS:
            inputs={}
            for geometry in GEOMETRIES:
                evidence=source['evidence'][metal][geometry];native=read(evidence['native_result']);t=native['input']
                require(evidence['status']==native['status']=='available','Required accepted diagonal missing')
                original={k:t[k] for k in ('xyz','charge','multiplicity','atom_count','valence_electrons')}
                original.update(MACE=evidence['MACE'],MACE_eV=evidence['MACE_eV'],native_frozen_ALPB_result=evidence['native_result'],native_frozen_ALPB_hartree=evidence['G_qvac_hartree'])
                inputs[geometry]={'original':original,**{k:t[k] for k in ('maps','state_signature','source_preparation')},
                    'native_input_result':evidence['native_result'],'vacuum_charges':t['endpoints']['vacuum']['charges'],'saved_diagonal':native['vacuum']['result']}
            cases.append({'case_id':conditioning+'4_'+metal,'source_case_id':source['source_case_id'],'conditioning':conditioning,'sample':4,'metal':metal,'inputs':inputs})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name in ('FACTORIAL.py','RUN.sbatch'):shutil.copyfile(HERE/name,out/name)
    manifest={'protocol':PROTOCOL,'implementation':pin(out/'FACTORIAL.py'),'wrapper':pin(out/'RUN.sbatch'),'plan':pin(HERE/'PLAN.md'),
        'interpretation':pin(HERE.parent/'a0ac_asn_solvent/INTERPRETATION.md'),'parent_config':pin(cp),'parent_collection':pin(rp),'parent_execution':pin(ep),
        'python':c['numerical_python'],'cases':cases,'geometries':list(GEOMETRIES),'declared_sources':2,'declared_metal_pairs':4,
        'declared_factorial_cells':16,'reused_diagonals':8,'maximum_new_native_calls':8,'new_diagonal_calls':0,'CPUs':1,
        'native_manifest':c['native_manifest'],'native_build':c['native_build'],'native_binary':c['native_binary'],'native_adapter':c['native_adapter'],
        'factorial':c['factorial'],'component':c['component'],'solvent_adapter':c['solvent_adapter'],'runtime':c['runtime'],
        'new_SCC_calls':0,'new_MACE_calls':0,'new_DFT_calls':0,'new_optimizations':0,'classification':None,'output':str(out/'execution')}
    write(out/'manifest.json',manifest);print(json.dumps(pin(out/'manifest.json')))

def allocation():
    job=os.environ.get('SLURM_JOB_ID');require(job and job.isdigit(),'Actual Slurm required before scientific imports')
    raw=subprocess.check_output(['scontrol','show','job','-o',job],text=True);fields=dict(x.split('=',1) for x in raw.split() if '=' in x)
    require(fields['JobState']=='RUNNING' and fields['NumCPUs']=='1' and fields['NumNodes']=='1' and os.environ.get('SLURM_CPUS_ON_NODE')=='1','One actual allocated CPU required')
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):require(os.environ.get(k)=='1','Unconstrained numeric threads: '+k)
    return {'job_id':job,'hostname':socket.gethostname(),'allocated_CPUs':1,'scontrol':raw}

def values(text,label,length):
    lines=re.findall(r'^'+label+r'\s+(.+)$',text,re.M);require(len(lines)==1,'Missing unique '+label)
    out=list(map(float,lines[0].split()));require(len(out)==length and all(math.isfinite(x) for x in out),'Invalid '+label);return out

def diagonal(m,native,group_input):
    """Reparse the actual completed native command; no duplicate evaluation."""
    rr=read(group_input['saved_diagonal']);record=read(group_input['native_input_result']);t=record['input'];ex=read(rr['execution']);text=verify(ex['stdout']).read_text()
    for k in ('xyz','charges','cavity','born_matrix_with_HB'):verify(rr[k])
    for k in ('stdout','stderr','resource_receipt'):verify(ex[k])
    require(ex['returncode']==0 and ex['command']==[m['native_binary']['path'],rr['xyz']['path'],rr['charges']['path'],str(t['atom_count'])],'Archived diagonal executable/input identity differs')
    require(rr['xyz']['sha256']==t['xyz']['sha256'] and rr['charges']['sha256']==group_input['vacuum_charges']['sha256'],'Archived diagonal XYZ/charge identity differs')
    require(rr['charge_source']=='vacuum' and 'internal GFN2-xTB/ALPB' in text,'Archived diagonal model differs')
    parsed=dict(zip(TERMS,values(text,'COMPONENTS',4)));bounds=dict(zip(('Gelec','Ghb'),values(text,'ROUNDING_BOUNDS',2)))
    require(parsed==rr['components_Eh'] and bounds==rr['rounding_bounds_Eh'],'Actual archived diagonal stdout differs from stored result')
    require(abs(math.fsum(parsed.values())-rr['G_Eh'])<=native['tolerances']['native_component_Eh'],'Archived diagonal component closure failed')
    charge_sum,=values(text,'CHARGE_SUM',1);require(abs(charge_sum-t['charge'])<=t['atom_count']*native['tolerances']['charge_rounding_e']+1e-12,'Archived printed charge sum differs')
    return rr|{'result':group_input['saved_diagonal'],'status':'available','reused':True,'diagonal_qualification':'Exact actual command/XYZ/charges/binary; raw four terms and bounds match saved replay','printed_charge_sum_e':charge_sum}

def precision(native,row):
    allowance=native['tolerances']['native_component_Eh'];limit=native['tolerances']['maximum_informative_bound_kcal_mol'];ha=native['HA_TO_KCAL']
    return [{'component':k,'rounding_bound_Eh':row['rounding_bounds_Eh'].get(k,0.0),'tolerance_Eh':row['rounding_bounds_Eh'].get(k,0.0)+allowance,
        'informative':(row['rounding_bounds_Eh'].get(k,0.0)+allowance)*ha<=limit} for k in TERMS]

def preflight(m,n,a,r):
    require(m['protocol']==PROTOCOL and m['implementation']==pin(__file__) and m['geometries']==list(GEOMETRIES),'Wrong frozen implementation/protocol')
    require((m['declared_factorial_cells'],m['reused_diagonals'],m['maximum_new_native_calls'],m['new_diagonal_calls'],len(m['cases']))==(16,8,8,0,4),'Finite scope changed')
    for k in ('wrapper','plan','interpretation','parent_collection','parent_execution'):verify(m[k])
    require(pin(sys.executable)==m['python'],'Qualified numerical interpreter differs');parent=read(m['parent_config']);native=read(m['native_manifest']);build=read(m['native_build'])
    for k in ('native_manifest','native_build','native_binary','native_adapter','factorial','component','solvent_adapter','runtime'):require(m[k]==parent[k],'Parent scientific helper differs: '+k)
    require(native['implementation']['NATIVE.py']==m['native_adapter'] and build['manifest']==m['native_manifest'] and build['binary']==m['native_binary'],'Build/implementation identity differs')
    for p in list(native['implementation'].values())+native['build']['pins']+build['linked_libraries']:verify(p)
    verify(m['native_binary']);verify(parent['corrected_ORCA_parameter_provenance'])
    require(read(parent['ORCA_correction_policy'])['candidate_sha256']==parent['corrected_ORCA_parameter_provenance']['sha256'],'Parameter correction lineage differs')
    expected={f'{c}4_{z}' for c in METALS for z in METALS};require({g['case_id'] for g in m['cases']}==expected,'Both sample4 sources and both metals required')
    groups={}
    for group in m['cases']:
        inp=group['inputs'];require(set(inp)==set(GEOMETRIES),'Missing geometry/charge owner');evidence={};atoms={};states={}
        for geometry in GEOMETRIES:
            x=inp[geometry];record=read(x['native_input_result']);t=record['input']
            require(record['manifest']==m['native_manifest'] and x['original']['native_frozen_ALPB_result']==x['native_input_result'],'Native diagonal lineage differs')
            require(t['metal']==group['metal'] and t['atom_count']==204 and t['charge']=={'Ca':-4,'La':-3}[group['metal']] and t['multiplicity']==1 and t['valence_electrons']==636,'204-atom paired state differs')
            require(t['endpoints']['vacuum']['charges']==x['vacuum_charges'],'Charge-owner pin differs')
            admitted=a.reusable(parent,x,n,r);require(admitted['status']=='available','Original component qualification unavailable')
            v=diagonal(m,native,x);p=precision(native,v);require(all(q['informative'] for q in p),'Diagonal stored-charge precision failed')
            v['precision']=p;v['historical_absolute_point_accepted']=admitted['historical_absolute_point_accepted'];evidence[geometry]=v
            atoms[geometry]=[line.split()[0] for line in verify(t['xyz']).read_text().splitlines()[2:]];states[geometry]=t
        t0,t1=[states[q] for q in GEOMETRIES]
        for k in ('metal','atom_count','charge','multiplicity','valence_electrons','maps','state_signature','source_preparation'):require(t0[k]==t1[k],'Crossed atom/state mapping differs: '+k)
        require(atoms[GEOMETRIES[0]]==atoms[GEOMETRIES[1]],'Crossed atom order differs')
        groups[group['case_id']]={'status':'passed','diagonals':evidence,'exact_atom_mapping':'Same physical maps, source preparation, state and ordered204 element sites; charges remain at their original indices'}
    return {'status':'passed','groups':groups,'reused_diagonals':8,'new_native_calls':0,'diagonal_qualification':'Eight actual saved binary executions verified and reparsed; no synthetic diagonal or duplicate native call'}

def combination(matrix,coefficients,native):
    """Conservative stored-charge and unchanged native-allowance propagation."""
    allowance=native['tolerances']['native_component_Eh'];ha=native['HA_TO_KCAL'];ceiling=native['tolerances']['maximum_informative_bound_kcal_mol']
    terms={k:math.fsum(coef*matrix[key]['components_Eh'][k] for key,coef in coefficients.items()) for k in TERMS}
    bounds={k:math.fsum(abs(coef)*(matrix[key]['rounding_bounds_Eh'].get(k,0.0)+allowance) for key,coef in coefficients.items()) for k in TERMS}
    total=math.fsum(terms.values());bound=math.fsum(bounds.values())
    return {'coefficients':coefficients,'components_Eh':terms,'total_Eh':total,'components_kcal_mol':{k:v*ha for k,v in terms.items()},'total_kcal_mol':total*ha,
        'conservative_bounds_Eh':bounds,'total_conservative_bound_Eh':bound,'component_informative_flags':{k:v*ha<=ceiling for k,v in bounds.items()},
        'total_informative':bound*ha<=ceiling,'informative_ceiling_kcal_mol':ceiling,'precision_scope':'Conservative independent stored-charge rounding plus existing native allowance; false aggregate flags are retained, not hidden or relaxed'}

def collect(m,pf,crossed,native):
    rows=[]
    for group in m['cases']:
        gid=group['case_id'];d=pf['groups'][gid]['diagonals'];matrix={'00':d[GEOMETRIES[0]],'11':d[GEOMETRIES[1]],
            **{key:crossed.get(gid+'__'+key,{'status':'unavailable','reason':'not executed'}) for key in ('01','10')}}
        row={k:group[k] for k in ('case_id','source_case_id','conditioning','sample','metal')};row.update(matrix=matrix,status='unavailable',decomposition=None)
        if all(v['status']=='available' for v in matrix.values()):
            pieces={name:combination(matrix,coefs,native) for name,coefs in (
                ('diagonal_change',{'11':1,'00':-1}),('geometry_at_q0',{'10':1,'00':-1}),('charge_at_x1',{'11':1,'10':-1}),
                ('charge_at_x0',{'01':1,'00':-1}),('geometry_at_q1',{'11':1,'01':-1}),('interaction',{'11':1,'10':-1,'01':-1,'00':1}))}
            closure={}
            for label,names in (('geometry_then_charge',('geometry_at_q0','charge_at_x1')),('charge_then_geometry',('charge_at_x0','geometry_at_q1'))):
                closure[label]={k:math.fsum(pieces[p]['components_Eh'][k] for p in names)-pieces['diagonal_change']['components_Eh'][k] for k in TERMS}
            closure['interaction_geometry']={k:(pieces['geometry_at_q1']['components_Eh'][k]-pieces['geometry_at_q0']['components_Eh'][k])-pieces['interaction']['components_Eh'][k] for k in TERMS}
            closure['interaction_charge']={k:(pieces['charge_at_x1']['components_Eh'][k]-pieces['charge_at_x0']['components_Eh'][k])-pieces['interaction']['components_Eh'][k] for k in TERMS}
            require(all(abs(v)<=native['tolerances']['native_component_Eh'] for path in closure.values() for v in path.values()),'Factorial path closure failed')
            row.update(status='available',decomposition=pieces,closure_residuals_Eh=closure)
        rows.append(row)
    return {'protocol':PROTOCOL,'terminal':True,'declared_sources':2,'declared_metal_pairs':4,'declared_factorial_cells':16,'reused_diagonal_cells':8,
        'available_new_offdiagonal_cells':sum(v['status']=='available' for v in crossed.values()),'available_metal_factorials':sum(row['status']=='available' for row in rows),
        'axis_definition':{'geometry0':'adaptive_La','geometry1':'amide_flip','charge0':'actual adaptive_La vacuum GFN2 charges','charge1':'actual amide_flip vacuum GFN2 charges','cell_ij':'geometry i, charge vector j'},
        'rows':rows,'new_SCC_calls':0,'new_MACE_calls':0,'new_DFT_calls':0,'new_optimizations':0,'classification':None,'new_calibration':False,
        'interpretation':'Counterfactual fixed-distribution diagnostic: geometry includes Born/reaction geometry and HB/SASA weights, not pure cavity volume. Both paths and interaction retained; no unique attribution or physical-state validation.'}

def execute(path,sha):
    placement=allocation();mp={'path':str(Path(path).resolve()),'sha256':sha};m=read(mp);out=Path(m['output']);start=time.monotonic();cpu=time.process_time();pf=None;crossed={};old_handlers={}
    receipt={'status':'failed','manifest':mp,'placement':placement,'new_SCC_calls':0,'new_MACE_calls':0,'new_DFT_calls':0,'new_optimizations':0}
    with (Path(path).parent/'RUN.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);out.mkdir(exist_ok=False)
        try:
            sys.path.insert(0,str(ROOT/'scripts'));n=load('_asn_factorial_native',m['native_adapter']);f=load('_asn_factorial_environment',m['factorial']);r=load('_asn_factorial_component',m['component']);a=load('_asn_factorial_saved_admission',m['solvent_adapter'])
            native=read(m['native_manifest']);pf=preflight(m,n,a,r);write(out/'PREFLIGHT.json',pf);runtime=read(m['runtime']);env=f.clean_env(runtime['xtb_parameter_file']['path'])
            def stop(signum,frame):n.STOP.set()
            for sig in (signal.SIGTERM,signal.SIGINT):old_handlers[sig]=signal.signal(sig,stop)
            (out/'tasks').mkdir()
            for group in m['cases']:
                for key in ('01','10'):
                    gid=group['case_id'];xowner,qowner=[GEOMETRIES[int(i)] for i in key];x=group['inputs'][xowner];q=group['inputs'][qowner];tid=gid+'__'+key;directory=out/'tasks'/tid
                    row={'status':'unavailable','geometry_owner':xowner,'charge_owner':qowner,'geometry_native_input':x['native_input_result'],'charge_native_input':q['native_input_result'],'xyz_input':x['original']['xyz'],'charge_input':q['vacuum_charges'],'reason':None,'reused':False}
                    try:
                        require(not n.STOP.is_set(),'Interrupted before new native call')
                        task={'xyz':x['original']['xyz'],'atom_count':204,'endpoints':{'vacuum':{'charges':q['vacuum_charges']}}}
                        value=n.evaluate(native,task,'vacuum',m['native_binary'],directory,env);row.update(value)
                        same=pf['groups'][gid]['diagonals'][xowner]
                        require(all(value[k]['sha256']==same[k]['sha256'] for k in ('cavity','born_matrix_with_HB')),'Crossed charges changed geometry-specific native cavity')
                        for k in ('Gsasa','Gshift'):require(abs(value['components_Eh'][k]-same['components_Eh'][k])<=native['tolerances']['native_component_Eh'],'Charge-independent native term changed: '+k)
                        ex=read(value['execution']);text=verify(ex['stdout']).read_text();charge_sum,=values(text,'CHARGE_SUM',1)
                        require(abs(charge_sum-q['original']['charge'])<=204*native['tolerances']['charge_rounding_e']+1e-12,'Crossed native charge sum failed')
                        p=precision(native,value);require(all(v['informative'] for v in p),'Crossed stored-charge precision failed')
                        row.update(status='available',precision=p,printed_charge_sum_e=charge_sum)
                    except Exception as e:row.update(reason=type(e).__name__+': '+str(e),traceback=traceback.format_exc())
                    crossed[tid]=row;write(out/(tid+'_RESULT.json'),row)
                    if n.STOP.is_set():raise RuntimeError('Interrupted; no further native calls')
            receipt['status']='complete' if len(crossed)==8 and all(v['status']=='available' for v in crossed.values()) else 'partial'
        except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc());raise
        finally:
            try:
                write(out/'PARTIAL_RESULTS.json',{'manifest':mp,'rows':crossed})
                for sig,handler in old_handlers.items():signal.signal(sig,handler)
                if pf:write(out/'FINAL_COLLECTION.json',{'manifest':mp,**collect(m,pf,crossed,native)})
            except BaseException as e:receipt.update(status='failed',collection_failure={'reason':str(e),'traceback':traceback.format_exc()});raise
            finally:
                receipt.update(elapsed_seconds=time.monotonic()-start,parent_process_CPU_seconds=time.process_time()-cpu)
                try:
                    receipt.update(actual_native_execution_receipts=[pin(p) for p in sorted(out.glob('tasks/*/native_EXECUTION.json'))],
                        new_native_calls=sum(1 for _ in out.glob('tasks/*/native_EXECUTION.json')),new_diagonal_calls=0,reused_diagonals=8 if pf else 0)
                except BaseException as e:receipt.update(status='failed',resource_reporting_failure=str(e))
                write(out/'EXECUTION.json',receipt)
    if receipt['status']=='failed':raise RuntimeError('Execution/reporting failed; inspect EXECUTION.json')

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('prepare');a.add_argument('--output',required=True)
    a=s.add_parser('execute');a.add_argument('--manifest',required=True);a.add_argument('--manifest-sha',required=True)
    args=p.parse_args()
    if args.action=='prepare':prepare(args.output)
    else:execute(args.manifest,args.manifest_sha)
