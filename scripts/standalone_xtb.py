"""Finite installed-xTB backend check on pinned PQQ pools and donor derivatives."""
from __future__ import annotations
import argparse
import concurrent.futures
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,BOHR_TO_A,read_json,record,verify,write_new,xyz
from affordable_response import read_engrad
from mace_site_kinematics import Kinematics
from native_pool_continuation import CASES,CANDIDATES
from native_solvent_force import OFFSETS,derivative,derivative_check

PROTOCOL='standalone_xtb671_PQQ_pool_and_force_accuracy_v1'
XTB=Path('/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb')
ACCURACIES=('0.2','0.02')
SETTINGS={'accuracies':list(ACCURACIES),'reported_accuracy':'0.02','gfn':2,'temperature_K':300,
          'maximum_iterations':500,'restart':False,'solvent':'water','reference_state':'gsolv',
          'surface_grid':'normal','surface_points':230,'ionic_strength_M':0,
          'workers':8,'threads_per_task':8,'cell_tolerance_kcal_mol':.1,'contrast_tolerance_kcal_mol':.2,
          'derivative_absolute_tolerance_kcal_mol_radian':.2,'derivative_relative_tolerance':.05,
          'maximum_calls':224,'q0_analytic_gradients':16}

def control():return '$scc\n maxiterations=500\n$gbsa\n grid=normal\n ion_st=0\n$write\n gbsa=true\n$end\n'
def command(t,executable):
    a=[str(executable),'core.xyz','--gfn','2','--chrg',str(t['charge']),'--uhf','0','--etemp','300',
       '--iterations','500','--acc',t['accuracy'],'--norestart','--parallel','8','--input','xcontrol','--json']
    if t['medium']=='alpb':a+=['--alpb','water','gsolv']
    return a+(['--grad'] if t['gradient_requested'] else ['--sp'])

def prepare(pool_manifest,force_design,agreement,output):
    pm=read_json(pool_manifest);fd=read_json(force_design)
    if len(pm['tasks'])!=80 or len(fd['cells'])!=32:raise InvalidArtifact('source population differs')
    expected={(c,q,z,s) for c in CASES for q in CANDIDATES for z in ('Ca','La') for s in ('vacuum','alpb')}
    if {(t['case_id'],t['candidate'],t['metal'],t['medium']) for t in pm['tasks']}!=expected:raise InvalidArtifact('pool cells differ')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir()
    pins={}
    for p in Path(__file__).parent.glob('*.py'):
        q=impl/p.name;shutil.copyfile(p,q);pins[p.name]=record(q)
    params=out/'parameters';params.mkdir();src=XTB.parent.parent/'share/xtb/param_gfn2-xtb.txt'
    shutil.copyfile(src,params/src.name)
    allcells=[{**t,'kind':'pool','cell_id':t['task_id'],'offset_radian':0. if t['candidate']=='origin' else None} for t in pm['tasks']]
    allcells +=[{**c,'kind':'displacement','candidate':'displacement','gradient_requested':False} for c in fd['cells']]
    tasks=[]
    for acc in ACCURACIES:
        for c in allcells:
            tid='acc'+acc.replace('.','p')+'__'+c['cell_id'];d=out/'tasks'/tid;d.mkdir(parents=True)
            shutil.copyfile(verify(c['xyz']),d/'core.xyz');(d/'xcontrol').write_text(control())
            t={k:c[k] for k in ('kind','cell_id','case_id','candidate','metal','medium','charge','multiplicity','gradient_requested','offset_radian')}
            t.update(task_id=tid,accuracy=acc,xyz=record(d/'core.xyz'),source_xyz=c['xyz'],input=record(d/'xcontrol'),directory=str(d))
            tasks.append(t)
    m={'protocol_id':PROTOCOL,'settings':SETTINGS,'tasks':tasks,'pool_manifest':record(pool_manifest),'force_design':record(force_design),
       'agreement':record(agreement),'implementation':pins,'executable':record(XTB),'parameter':record(params/src.name),
       'installed_parameter':record(src),'package':record(XTB.parent.parent/'conda-meta/xtb-6.7.1-h8876d29_4.json'),
       'inputs':pm['inputs'],'reference':pm['reference'],'inventory':pm['inventory'],'new_MACE_DFT_calls':0,'production_changed':False}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v

def validate(manifest,fresh=False):
    m=read_json(manifest);pm=read_json(verify(m['pool_manifest']));fd=read_json(verify(m['force_design']))
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or len(m['tasks'])!=224:raise InvalidArtifact('fixed backend/scope changed')
    for k in ('agreement','executable','parameter','installed_parameter','package','reference','inputs','inventory'):verify(m[k])
    if verify(m['parameter']).read_bytes()!=verify(m['installed_parameter']).read_bytes():raise InvalidArtifact('parameters changed')
    if {p.name for p in verify(m['parameter']).parent.iterdir()}!={'param_gfn2-xtb.txt'}:raise InvalidArtifact('parameter/rc override')
    for p in m['implementation'].values():verify(p)
    sources={('pool',c['task_id']):c for c in pm['tasks']}|{('displacement',c['cell_id']):c for c in fd['cells']}
    expected={(a,k,cid) for a in ACCURACIES for k,cid in sources}
    if {(t['accuracy'],t['kind'],t['cell_id']) for t in m['tasks']}!=expected:raise InvalidArtifact('case/accuracy membership changed')
    for t in m['tasks']:
        s=sources[t['kind'],t['cell_id']]
        if any(t[k]!=s[k] for k in ('case_id','metal','medium','charge','multiplicity')) or t['multiplicity']!=1:raise InvalidArtifact('state differs')
        grad=t['kind']=='pool' and s['gradient_requested']
        if t['gradient_requested']!=grad:raise InvalidArtifact('analytic center set differs')
        if verify(t['xyz']).read_bytes()!=verify(s['xyz']).read_bytes() or verify(t['input']).read_text()!=control():raise InvalidArtifact('coordinates/control differ')
        d=Path(t['directory']);d.relative_to(Path(manifest).resolve().parent/'tasks')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','xcontrol'}:raise InvalidArtifact('nonfresh task')
    if sum(t['gradient_requested'] for t in m['tasks'])!=16:raise InvalidArtifact('analytic count differs')
    return {'status':'validated','manifest':record(manifest),'tasks':224,'pool_cells':160,'displaced_cells':64,'gradient_calls_within_total':16,'new_calls_in_preflight':0}

def execute_cell(manifest,m,t):
    d=Path(t['directory']);rp=d/'execution.json'
    if rp.exists():raise InvalidArtifact('existing attempt must not rerun')
    env=os.environ.copy();env.update(XTBPATH=str(verify(m['parameter']).parent),OMP_NUM_THREADS='8,1',OMP_MAX_ACTIVE_LEVELS='1',
         OMP_STACKSIZE='1G',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    argv=command(t,verify(m['executable']));begin=time.monotonic();code=None;error=None
    try:
        with (d/'endpoint.out').open('x') as out:code=subprocess.run(argv,cwd=d,env=env,stdout=out,stderr=subprocess.STDOUT,check=False).returncode
    except Exception as exc:error=str(exc)
    receipt={'task_id':t['task_id'],'manifest':record(manifest),'command':argv,'returncode':code,'launch_error':error,
       'wall_seconds':time.monotonic()-begin,'threads':8,'model_call_started':code is not None,'executable':m['executable'],
       'parameter':m['parameter'],'job_id':os.environ.get('SLURM_JOB_ID'),'environment':{k:env[k] for k in ('XTBPATH','OMP_NUM_THREADS','OMP_MAX_ACTIVE_LEVELS','OMP_STACKSIZE','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')},
       'artifacts':{p.name:record(p) for p in d.iterdir() if p.is_file() and p.name not in ('core.xyz','xcontrol')}}
    write_new(rp,receipt);print(json.dumps({'task':t['task_id'],'returncode':code,'wall_seconds':receipt['wall_seconds']}),flush=True)
    return record(rp)

def execute(manifest):
    validate(manifest,True);m=read_json(manifest);root=Path(manifest).parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE',0))!=64:raise InvalidArtifact('64CPU allocation required')
    with (root/'execute.lock').open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);begin=time.monotonic();receipts=[];errors=[]
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
            futures={ex.submit(execute_cell,manifest,m,t):t for t in m['tasks']}
            for future in concurrent.futures.as_completed(futures):
                try:receipts.append(future.result())
                except Exception as e:errors.append({'task_id':futures[future]['task_id'],'error':str(e)})
        elapsed=time.monotonic()-begin
        result={'manifest':record(manifest),'receipts':receipts,'errors':errors,'wall_seconds':elapsed,'allocated_core_seconds':elapsed*64,
                'allocated_cpus':64,'job_id':os.environ['SLURM_JOB_ID'],'GPU_seconds':0}
        write_new(root/'EXECUTION.json',result)
    return {'receipts':len(receipts),'errors':len(errors)}

def parse_task(manifest,t):
    d=Path(t['directory']);row={k:t[k] for k in ('task_id','kind','cell_id','accuracy','case_id','candidate','metal','medium','offset_radian')}
    row.update(status='unavailable',energy_hartree=None,gradient=None,reason=None)
    try:
        rec=read_json(d/'execution.json');text=verify(rec['artifacts']['endpoint.out']).read_text()
        if rec['manifest']!=record(manifest) or rec['task_id']!=t['task_id'] or rec['returncode']!=0 or 'normal termination of xtb' not in text:raise InvalidArtifact('execution/normal termination failed')
        if re.search(r'failed to converge|NOT CONVERGED|abnormal termination',text,re.I):raise InvalidArtifact('SCC not converged')
        energies=re.findall(r'TOTAL ENERGY\s+([-+\d.Ee]+)\s+Eh',text)
        if not energies:raise InvalidArtifact('total energy absent')
        energy=float(energies[-1]);atoms=xyz(verify(t['xyz']))
        charge=re.findall(r'net charge\s+([-+\d.]+)',text);spin=re.findall(r'unpaired electrons\s+([-+\d.]+)',text)
        if not charge or float(charge[-1])!=t['charge'] or not spin or float(spin[-1])!=0:raise InvalidArtifact('printed charge/spin differs')
        data=read_json(verify(rec['artifacts']['xtbout.json']))
        if data['method']!='GFN2-xTB' or not data['xtb version'].startswith('6.7.1') or data['number of unpaired electrons']!=0 or data['number of electrons']%2:
            raise InvalidArtifact('actual backend/electronic state differs')
        if abs(sum(data['partial charges'])-t['charge'])>1e-4 or abs(data['total energy']-energy)>1e-7:
            raise InvalidArtifact('charge closure or JSON/printed energy differs')
        required=[r'max\. iterations\s+500',r'restarted\?\s+false',r'electronic temp\.\s+300\.',r'PC potential\s+false',r'convergence criteria satisfied after']
        if any(not re.search(pattern,text) for pattern in required):raise InvalidArtifact('actual SCC settings/convergence differ')
        acc=re.findall(r':  accuracy\s+([\d.Ee+-]+)',text)
        if not acc or float(acc[-1])!=float(t['accuracy']):raise InvalidArtifact('actual accuracy differs')
        if t['medium']=='alpb' and any(term not in text for term in ['Solvation model:               ALPB','Reference state                gsolv','Grid points                               230','Ion screening                  false']):
            raise InvalidArtifact('actual ALPB reference/surface differs')
        if t['medium']=='vacuum' and 'Solvation model:               ALPB' in text:raise InvalidArtifact('unexpected solvent')
        cycles=re.findall(r'convergence criteria satisfied after\s+(\d+) iterations',text)
        row.update(scc_iterations=int(cycles[-1]),charge_sum_e=sum(data['partial charges']),effective_electrons=data['number of electrons'])
        row.update(energy_hartree=energy,receipt=record(d/'execution.json'),output=rec['artifacts']['endpoint.out'],xtbout=rec['artifacts']['xtbout.json'],
                   printed_charge=float(charge[-1]),printed_spin=float(spin[-1]),raw_json=data,
                   output_summary_lines=[line.strip() for line in text.splitlines() if any(k in line for k in ['accuracy','SCC convergence','wavefunction','Electronic Temp','electronic temp','Gsolv','Gelec','Gsasa','Ghb','Gshift','parameter','iterations','spin','charge','Solvent','solvent','grid'])])
        if t['gradient_requested']:
            import gemmi
            ep=next((d/n for n in ('xtb-orca.engrad','core.engrad') if (d/n).exists()),None)
            if ep is None or re.search('numerical gradient',text,re.I):raise InvalidArtifact('analytic gradient absent/unsupported')
            g=read_engrad(ep)
            if g['atom_count']!=len(atoms) or abs(g['energy_Ha']-energy)>1e-8 or not np.array_equal(g['atomic_numbers'],[gemmi.Element(a[0]).atomic_number for a in atoms]):raise InvalidArtifact('gradient energy/order differs')
            if not np.allclose(g['coordinates_bohr']*BOHR_TO_A,[a[1:] for a in atoms],atol=1e-6,rtol=0):raise InvalidArtifact('gradient coordinates differ')
            row['gradient']={'engrad':record(ep),'kcal_mol_A':(g['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A).tolist(),'analytic_driver':'xtb --grad, version6.7.1 writeResultsOrca','quantity':'gradient_not_force'}
        row['status']='complete'
    except Exception as exc:row['reason']=str(exc)
    return row

def collect(manifest,output):
    validate(manifest);m=read_json(manifest);rows=[parse_task(manifest,t) for t in m['tasks']]
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'rows':rows,'denominator':224,
        'complete':sum(r['status']=='complete' for r in rows),'analysis_implementation':record(__file__),'new_calls_in_collection':0}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}

def compare(collection,output):
    from nikasha_pool import choose_rows
    from native_pool_continuation import matrices
    from accommodation_folds_compare import decision
    r=read_json(collection);m=read_json(verify(r['manifest']));fd=read_json(verify(m['force_design']));pm=read_json(verify(m['pool_manifest']))
    by={(x['accuracy'],x['kind'],x['case_id'],x['candidate'],x['metal'],x['medium'],x['offset_radian']):x for x in r['rows']}
    native=matrices(verify(m['inventory']),Path(m['pool_manifest']['path']).parent/'collection.json')
    ref=read_json(verify(m['reference']));inputs=read_json(verify(m['inputs']));labels={x['case_id']:x['known_class'] for x in inputs['cases']}
    bands={'released_transfer_only':ref['old_frozen_bands'],'native_adaptive_transfer_only':ref['variants']['operational']['bands']}
    pools=[]
    for cid in CASES:
        levels={};sensitivity=[]
        for acc in ACCURACIES:
            matrix={z:{} for z in ('Ca','La')}
            for z in matrix:
                for q in CANDIDATES:
                    rows={s:by[acc,'pool',cid,q,z,s,0. if q=='origin' else None] for s in ('vacuum','alpb')}
                    good=all(x['status']=='complete' for x in rows.values())
                    matrix[z][q]={'status':'complete' if good else 'unavailable','components':{'MACE_eV':native[cid]['old']['matrix'][z][q]['components']['MACE_eV'],
                        'GFN2_vacuum_hartree':rows['vacuum']['energy_hartree'],'GFN2_ALPB_hartree':rows['alpb']['energy_hartree']} if good else None}
            p=choose_rows(matrix,CANDIDATES)
            value=p['operational']['composite_R_model_kcal_mol'] if p['operational'] else None
            levels[acc]={'matrix':matrix,'pool':p,'incompatible_native_band_transfer':{k:decision(value,b) for k,b in bands.items()}}
        deltas=[]
        for q in CANDIDATES:
            vals={}
            for z in ('Ca','La'):
                for s in ('vacuum','alpb'):
                    aa=[by[acc,'pool',cid,q,z,s,0. if q=='origin' else None] for acc in ACCURACIES]
                    delta=(aa[1]['energy_hartree']-aa[0]['energy_hartree'])*HA_TO_KCAL if all(x['status']=='complete' for x in aa) else None
                    vals[z+'_'+s]=delta;deltas.append(delta)
            contrast=vals['Ca_alpb']-vals['Ca_vacuum']-vals['La_alpb']+vals['La_vacuum'] if all(x is not None for x in vals.values()) else None
            sensitivity.append({'candidate':q,'cell_deltas_kcal_mol':vals,'contrast_delta_kcal_mol':contrast,'pass':contrast is not None and abs(contrast)<=.2})
        changes={mode:(levels['0.02']['pool'][mode]['composite_R_model_kcal_mol']-levels['0.2']['pool'][mode]['composite_R_model_kcal_mol']) if all(levels[a]['pool'][mode] for a in ACCURACIES) else None for mode in ('mathematical','operational')}
        good=all(x is not None and abs(x)<=.1 for x in deltas) and all(x['pass'] for x in sensitivity) and all(v is not None and abs(v)<=.2 for v in changes.values())
        pools.append({'case_id':cid,'known_class_for_report_only':labels[cid],'levels':levels,'native_primary':native[cid]['old']['pool'],
            'native_continued':native[cid]['pool'],'sensitivity':sensitivity,'pool_change_kcal_mol':changes,'energy_qualification_pass':good,
            'qualified_R':levels['0.02']['pool']['operational']['composite_R_model_kcal_mol'] if good else None})
    derivatives=[]
    for cid in ('4MAE','q88jh5-pqq-la_model'):
        for acc in ACCURACIES:
            values={};cell_deltas=[];contrast_changes=[]
            for z in ('Ca','La'):
                meta=fd['maps'][cid+'__'+z];kin=Kinematics(read_json(verify(meta['physical_mapping']))['context']);jac=kin.evaluate(np.zeros(len(kin.modes)))[3][meta['mode_index']]
                for s in ('vacuum','alpb'):
                    center=by[acc,'pool',cid,'origin',z,s,0.];ds={off:by[acc,'displacement',cid,'displacement',z,s,off] for off in OFFSETS}
                    if center['status']!='complete' or not center['gradient'] or any(x['status']!='complete' for x in ds.values()):values[z+'_'+s]=None;continue
                    values[z+'_'+s]={**derivative({off:x['energy_hartree'] for off,x in ds.items()}),'analytic':float(np.sum(jac*np.array(center['gradient']['kcal_mol_A'])))}
            for z in ('Ca','La'):
                a,b=(values[z+'_'+s] for s in ('alpb','vacuum'));values[z+'_solvent']={k:a[k]-b[k] for k in ('h','half','analytic')} if a and b else None
            a,b=(values[z+'_solvent'] for z in ('Ca','La'));values['Ca_minus_La_solvent']={k:a[k]-b[k] for k in ('h','half','analytic')} if a and b else None
            speed=fd['maps'][cid+'__Ca']['RMS_moving_heavy_speed_A_radian'];checks=[]
            for name,v in values.items():checks.append({'quantity':name,**derivative_check(v['h'],v['half'],v['analytic'],speed)} if v else {'quantity':name,'pass':False,'status':'unavailable'})
            for off in (*OFFSETS,0.):
                dd={}
                for z in ('Ca','La'):
                    for s in ('vacuum','alpb'):
                        aa=[by[level,'pool' if off==0 else 'displacement',cid,'origin' if off==0 else 'displacement',z,s,off] for level in ACCURACIES]
                        delta=(aa[1]['energy_hartree']-aa[0]['energy_hartree'])*HA_TO_KCAL if all(x['status']=='complete' for x in aa) else None
                        dd[z+'_'+s]=delta;cell_deltas.append(delta)
                contrast=dd['Ca_alpb']-dd['Ca_vacuum']-dd['La_alpb']+dd['La_vacuum'] if all(x is not None for x in dd.values()) else None
                contrast_changes.append({'offset_radian':off,'delta_kcal_mol':contrast,'pass':contrast is not None and abs(contrast)<=.2})
            settling=all(x is not None and abs(x)<=.1 for x in cell_deltas) and all(x['pass'] for x in contrast_changes)
            derivatives.append({'case_id':cid,'accuracy':acc,'checks':checks,'energy_settling_pass':settling,'contrast_changes':contrast_changes,
                'all_derivatives_pass':all(x['pass'] for x in checks),'qualified':settling and all(x['pass'] for x in checks)})
    result={'protocol_id':PROTOCOL,'collection':record(collection),'analysis_implementation':record(__file__),'pool_cases':pools,
        'derivative_cases':derivatives,'reported_accuracy':'0.02','own_reference':None,'new_calls_in_comparison':0,'production_changed':False,
        'settings':SETTINGS,'case_denominator':4,'gradient_quantities_per_accuracy':14,'derivative_source_maps':fd['maps']}
    write_new(output,result);return {'energy_qualified_pools':sum(x['energy_qualification_pass'] for x in pools),'derivative_checks_tight':sum(c['pass'] for d in derivatives if d['accuracy']=='0.02' for c in d['checks'])}

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op in ('prepare','validate','execute','collect','compare'):
        a=sub.add_parser(op)
        names={'prepare':('pool-manifest','force-design','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'compare':('collection','output')}[op]
        for n in names:a.add_argument('--'+n,required=True,type=Path)
        if op=='validate':a.add_argument('--fresh',action='store_true')
    args=vars(p.parse_args());fn=globals()[args.pop('op')];print(json.dumps(fn(**args),indent=2))
if __name__=='__main__':main()
