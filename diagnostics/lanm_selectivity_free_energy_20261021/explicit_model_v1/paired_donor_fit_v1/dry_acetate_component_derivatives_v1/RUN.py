#!/usr/bin/env python3
"""Fixed16 consumed dry-acetate model component derivatives; no fitting."""
import argparse
import copy
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import socket
import time

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
POINT=[3.616264519265877,.10825158943928785,3.3717073252558194,.1962751789555681,3.2616716485427766,.08618293783769042,3.228372064914772,.19287981074489938]
PROTOCOL='dry_acetate_component_derivatives_v1'
METALS=('La','Dy');DIRECTIONS=('radial','rotation');STEPS=(.01,.005);RESOLUTION=.1
COMPONENTS=('Bond Stretching','Angle Bending','Stretch-Bend','Out-of-Plane Bend','Torsional Angle','Van der Waals','Atomic Multipoles','Polarization')

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed evidence: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,x):
    p=Path(path);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def module(p,name):
    sp=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def dec(s):return Decimal(s.replace('D','E').replace('d','e'))
def halfunit(s):return Decimal(5).scaleb(dec(s).as_tuple().exponent-1)

def parse_energy(path,helper):
    parsed=helper.parse_analyze(path);text=Path(path).read_text();number=helper.NUMBER
    total=re.findall(r'Total Potential Energy\s*:\s*('+number+')',text)
    tokens={};counts={}
    for line in text.split('Energy Component Breakdown :',1)[1].splitlines():
        m=re.fullmatch(r'\s*([A-Za-z][A-Za-z /()\-]+?)\s+('+number+r')\s+(\d+)\s*',line)
        if m:tokens[m[1].strip()]=m[2];counts[m[1].strip()]=int(m[3])
    require(set(tokens)==set(COMPONENTS) and len(total)==1,'Exact eight components required')
    tokens['Total']=total[0]
    with localcontext() as ctx:
        ctx.prec=50;remainder=dec(total[0])-sum((dec(tokens[k]) for k in COMPONENTS),Decimal(0));bound=sum((halfunit(s) for s in tokens.values()),Decimal(0))
        require(abs(remainder)<=bound,'Native component closure exceeds printed precision')
    return {'native':parsed,'tokens':tokens,'interactions':counts,'half_print_units_kcal_mol':{k:str(halfunit(s)) for k,s in tokens.items()},'component_sum_remainder_kcal_mol':str(remainder),'component_sum_rounding_bound_kcal_mol':str(bound)}

def context():
    cfg=read(pin(HERE/'INPUTS.json'));require(cfg['protocol']==PROTOCOL,'Protocol differs')
    require(cfg['specification']==pin(HERE/'SPEC.md'),'Specification differs')
    for key,pin_value in cfg.items():
        if key!='selected_fit_result' and isinstance(pin_value,dict) and set(pin_value)=={'path','sha256'}:verify(pin_value)
    # The exact completed fit pin is vouched by the immutable compact actual review; do not reread its trajectory.
    ref=module(cfg['reference_adapter'],'component_reference_geometry');p=ref.load(verify(cfg['reference_panel']))
    helper=module(cfg['energy_helper'],'component_energy_parser');native=module(cfg['native_helper'],'component_native_helper')
    evmod=module(cfg['evaluator'],'component_selected_renderer');ecfg=read(cfg['evaluator_config']);data=read(cfg['training_result']);review=read(cfg['selected_review'])
    require(review['status']=='passed' and review['error'] is None and review['result']==cfg['selected_fit_result'],'Exact admitted terminal candidate required')
    require([v['selected'] for v in review['parameter_bounds']]==POINT and len(review['selected_sources'])==28,'Selected eight-parameter point differs')
    require(ecfg['native_helper']==cfg['native_helper'] and ecfg['executables']['analyze']==cfg['analyze'] and ecfg['training_result']==cfg['training_result'] and ecfg['implementation']==cfg['evaluator'],'Original native/model joins differ')
    require(p['steps_A']==list(STEPS) and p['derivative_resolution_kcal_mol_A']==RESOLUTION and len(p['cells'])==16,'Shared finite scope differs')
    ev=evmod.Evaluator.__new__(evmod.Evaluator);ev.config=ecfg;ev.data=data;ev.family={'parameters':review['parameter_bounds']};ev.resolver=module(ecfg['resolver'],'component_resolver');ev._base_text={};ev._base_parameters={};central={}
    for metal in METALS:
        gid='matched_dry_acetate_mono_central_'+metal;row=data['rows'][gid];s=review['selected_sources'][gid];a=read(s['authority']);inp=s['effective_input']
        for q in inp.values():verify(q)
        old=ecfg['baseline_inputs'][gid];prm=verify(old['parameters']).read_text();key=verify(old['key']).read_text();ev._base_text[gid]=(prm,key);ev._base_parameters[gid]=ev.resolver.parse_parameters(prm,key)
        token,identity,rendered,render=ev.cache_identity(POINT,gid)
        require(a['status']=='complete' and a['cell_id']==gid and a['job_id']=='17621' and not s['reused_baseline'],'Original completed center required')
        require(a['cache_key']==token and a['cache_identity']==identity and a['render']==render and a['input']==inp,'Selected key/authority identity differs')
        require(verify(inp['key']).read_text()==rendered and all(inp[k]['sha256']==row['input'][k]['sha256'] for k in ('xyz','parameters')),'Selected physical bytes differ')
        require(row['geometry']==p['central'][metal]['original_cell']['geometry'],'Reference/model central geometry differs')
        n=native.parse(verify(s['native_outputs']['analyze_output']),verify(s['native_outputs']['gradient_output']),dict(row['input'],derivative_test=False))
        require(n==a['native'] and n['total_energy_kcal_mol']==a['energy'] and n['force_kcal_mol_A']==a['forces'],'Saved raw central E/F differs')
        energy=parse_energy(verify(s['native_outputs']['analyze_output']),helper);typed=helper.typed_rows(verify(inp['xyz']).read_text());g=row['geometry']
        require(len(typed)==8 and [[float(v) for v in r[2:5]] for r in typed]==g['positions_A'] and 36 not in [int(r[5]) for r in typed],'Exact dry central coordinates/types required')
        grad_text=verify(s['native_outputs']['gradient_output']).read_text();grad_tokens=re.findall(r'^\s*Anlyt\s+\d+\s+('+helper.NUMBER+')\s+('+helper.NUMBER+')\s+('+helper.NUMBER+')',grad_text,re.M)
        require(len(grad_tokens)==8,'Central gradient tokens missing');f=ref.frame(g);rounding={}
        for direction in DIRECTIONS:
            tangent=[f['radial_unit'] if direction=='radial' else ref.b.mul(ref.b.cross(f['rotation_unit'],ref.b.sub(g['positions_A'][i],f['centroid_A'])),1/f['rotation_length_A']) for i in f['ligand_indices']]
            rounding[direction]=sum(abs(tangent[j][k])*float(halfunit(grad_tokens[i][k])) for j,i in enumerate(f['ligand_indices']) for k in range(3))
        central[metal]={'geometry_id':gid,'source':s,'selected_render':render,'analytic':ref.projection(a['forces'],g),'analytic_gradient_print_bound_kcal_mol_A':rounding,'energy':energy,'typed_rows':typed,'geometry':g}
    return cfg,ref,p,helper,native,central

def render_xyz(cell,center,helper):
    raw=verify(cell['source_geometry']).read_text().splitlines();ordinary=[line.split() for line in raw[2:] if line.strip()]
    require(int(raw[0])==len(ordinary)==8 and [r[0] for r in ordinary]==cell['symbols'],'Shared ordinary XYZ malformed')
    require([[float(x) for x in r[1:4]] for r in ordinary]==cell['positions_A'],'Shared coordinate tokens differ')
    rows=copy.deepcopy(center['typed_rows'])
    for r,x in zip(rows,ordinary):r[2:5]=x[1:4]
    text='8  '+cell['cell_id']+' frozen donor17621\n'+'\n'.join(' '.join(r) for r in rows)+'\n'
    actual=helper.typed_rows(text)
    require([[r[0],r[1],*r[5:]] for r in actual]==[[r[0],r[1],*r[5:]] for r in center['typed_rows']],'Types/bonds/atom order changed')
    return text

def prepare(out):
    tick=time.monotonic();cpu=time.process_time();cfg,ref,p,helper,native,central=context();out.mkdir(parents=True,exist_ok=False);rows={}
    for c in p['cells']:
        cid=c['cell_id'];directory=out/'cells'/cid;directory.mkdir(parents=True);(directory/'system.xyz').write_text(render_xyz(c,central[c['metal']],helper))
        inp=central[c['metal']]['source']['effective_input'];shutil.copyfile(verify(inp['key']),directory/'system.key');shutil.copyfile(verify(inp['parameters']),directory/'amoeba09.prm')
        rows[cid]={'cell_id':cid,'metal':c['metal'],'direction':c['direction'],'shift_label':c['shift_label'],'q_A':c['q_A'],'source_geometry':c['source_geometry'],'input':{k:pin(directory/f) for k,f in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}}
    manifest={'protocol':PROTOCOL,'implementation':pin(__file__),'inputs':pin(HERE/'INPUTS.json'),'specification':cfg['specification'],'reference_panel':cfg['reference_panel'],'selected_review':cfg['selected_review'],'selected_fit_result':cfg['selected_fit_result'],'selected_parameters':POINT,'analyze':cfg['analyze'],'central':central,'rows':rows,'maximum_native_calls':16,'workers':1,'new_native_calls':0}
    save(out/'MANIFEST.json',manifest);load_manifest(out/'MANIFEST.json');result={'status':'dry_component_preparation_passed','implementation':pin(__file__),'manifest':pin(out/'MANIFEST.json'),'cells':16,'saved_central_raw_replays':2,'new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'process_CPU_seconds':time.process_time()-cpu}}
    save(out/'RESULT.json',result);print(json.dumps(result,indent=2))

def load_manifest(path):
    cfg,ref,p,helper,native,central=context();m=json.loads(Path(path).read_text())
    require(m['protocol']==PROTOCOL and m['implementation']==pin(__file__) and m['inputs']==pin(HERE/'INPUTS.json') and m['specification']==cfg['specification'],'Prepared source differs')
    require(m['reference_panel']==cfg['reference_panel'] and m['selected_review']==cfg['selected_review'] and m['selected_fit_result']==cfg['selected_fit_result'] and m['selected_parameters']==POINT and m['analyze']==cfg['analyze'],'Prepared lineage differs')
    require(m['central']==central and set(m['rows'])==set(ref.ALL) and m['maximum_native_calls']==16 and m['workers']==1,'Prepared central/scope differs')
    for c in p['cells']:
        r=m['rows'][c['cell_id']];require(all(r[k]==c[k] for k in ('cell_id','metal','direction','shift_label','q_A','source_geometry')),'Cell mapping differs')
        require(verify(r['input']['xyz']).read_text()==render_xyz(c,central[c['metal']],helper),'Displaced typed XYZ differs')
        for k in ('key','parameters'):require(r['input'][k]['sha256']==central[c['metal']]['source']['effective_input'][k]['sha256'],'Frozen key/parameters changed');verify(r['input'][k])
    return m,helper,native

def compare(m,rows):
    require(set(rows)==set(m['rows']) and len(rows)==16 and all(r['status']=='complete' for r in rows.values()),'All16 complete required; no partial comparison')
    output={};checks={}
    with localcontext() as ctx:
        ctx.prec=50
        for direction in DIRECTIONS:
            output[direction]={}
            for kind in ('La','Dy','paired'):
                coeff={'La':1} if kind=='La' else {'Dy':1} if kind=='Dy' else {'La':-1,'Dy':1};values={}
                analytic=sum(c*m['central'][metal]['analytic'][direction] for metal,c in coeff.items());analytic_bound=sum(abs(c)*m['central'][metal]['analytic_gradient_print_bound_kcal_mol_A'][direction] for metal,c in coeff.items())
                for step in STEPS:
                    h=Decimal(str(step));derivs={};bounds={}
                    for term in (*COMPONENTS,'Total'):
                        value=Decimal(0);bound=Decimal(0)
                        for metal,c in coeff.items():
                            samples={sign:next(rows[cid]['parsed'] for cid,r in m['rows'].items() if r['metal']==metal and r['direction']==direction and r['q_A']==sign*step) for sign in (-1,1)}
                            value+=c*(dec(samples[1]['tokens'][term])-dec(samples[-1]['tokens'][term]))/(2*h)
                            bound+=abs(c)*(halfunit(samples[1]['tokens'][term])+halfunit(samples[-1]['tokens'][term]))/(2*h)
                        derivs[term]=value;bounds[term]=bound
                    closure=derivs['Total']-sum((derivs[k] for k in COMPONENTS),Decimal(0));closure_bound=sum(bounds.values(),Decimal(0));require(abs(closure)<=closure_bound,'Finite-difference component closure exceeds print bound')
                    error=float(derivs['Total'])-analytic;label=f'{direction}/{kind}/h={step}/analytic';checks[label]=abs(error)<=RESOLUTION
                    values[str(step)]={'derivatives_kcal_mol_A':{k:float(v) for k,v in derivs.items()},'derivative_decimal_tokens':{k:str(v) for k,v in derivs.items()},'generalized_forces_kcal_mol_A':{k:-float(v) for k,v in derivs.items()},'print_bounds_kcal_mol_A':{k:float(v) for k,v in bounds.items()},'component_sum_remainder_kcal_mol_A':float(closure),'component_sum_bound_kcal_mol_A':float(closure_bound),'total_minus_analytic_kcal_mol_A':error,'total_analytic_difference_print_bound_kcal_mol_A':float(bounds['Total'])+analytic_bound,'analytic_resolution_passed':checks[label]}
                delta={k:values['0.005']['derivatives_kcal_mol_A'][k]-values['0.01']['derivatives_kcal_mol_A'][k] for k in (*COMPONENTS,'Total')};label=f'{direction}/{kind}/two_step';checks[label]=abs(delta['Total'])<=RESOLUTION
                output[direction][kind]={'central_analytic_kcal_mol_A':analytic,'central_gradient_print_bound_kcal_mol_A':analytic_bound,'steps':values,'smaller_minus_larger_step_kcal_mol_A':delta,'two_step_total_resolution_passed':checks[label],'smaller_step_descriptively_resolved':checks[label]}
    require(len(checks)==18,'Original numerical denominator differs')
    return {'directions':output,'numerical_checks':checks,'numerical_checks_passed':all(checks.values()),'numerical_check_count':18,'physical_accuracy_qualified':False,'component_physical_attribution_qualified':False}

def release_gate(m,path,out):
    r=json.loads(Path(path).read_text());require(r['status']=='released_dry_acetate_component_derivatives','Root release required')
    require(r['implementation']==pin(__file__) and r['manifest']==pin(m['_path']) and r['analyze']==m['analyze'] and r['selected_review']==m['selected_review'],'Release source/candidate differs')
    require(r['maximum_native_calls']==16 and r['workers']==1 and r['cell_ids']==list(m['rows']) and Path(r['output_directory']).resolve()==out.resolve(),'Release scope/path differs')
    review=read(r['independent_source_review']);require(review['status']=='passed_dry_acetate_component_source_review' and review['implementation']==pin(__file__) and review['manifest']==pin(m['_path']),'Source review differs')
    return r

def run(path,release,out):
    m,helper,native=load_manifest(path);m['_path']=str(Path(path).resolve());release_gate(m,release,out)
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Native work needs one Slurm CPU')
    require(all(os.environ.get(k)=='1' for k in THREADS),'Every numerical thread control must be1');out.mkdir(parents=True,exist_ok=False)
    tick=time.monotonic();cpu=time.process_time();children=resource.getrusage(resource.RUSAGE_CHILDREN)
    result={'status':'running','protocol':PROTOCOL,'implementation':pin(__file__),'manifest':pin(path),'release':pin(release),'analyze':m['analyze'],'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'rows':{cid:{'status':'unattempted','cell_id':cid} for cid in m['rows']},'new_native_calls':0,'error':None,'comparison':None}
    save(out/'START.json',{k:v for k,v in result.items() if k not in ('rows','comparison')});save(out/'RESULT.json',result)
    try:
        for cid,row in m['rows'].items():
            require(result['new_native_calls']<16,'Native call cap');cellout=out/'cells'/cid;cellout.mkdir(parents=True);rec={'cell_id':cid,'status':'started','receipts':[]};result['rows'][cid]=rec
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['input'][k]),cellout/name)
            rec['input']={k:pin(cellout/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))};save(cellout/'START.json',{'cell_id':cid,'input':rec['input'],'manifest':pin(path),'release':pin(release),'analyze':m['analyze']})
            result['new_native_calls']+=1;save(out/'RESULT.json',result)
            try:
                native.invoke([str(verify(m['analyze'])),'system.xyz','E'],cellout,cellout/'analyze.out',rec['receipts']);rec['parsed']=parse_energy(cellout/'analyze.out',helper)
                require(rec['parsed']['interactions']==m['central'][row['metal']]['energy']['interactions'],'Native component interaction inventory changed');rec['status']='complete'
            except Exception as exc:rec.update(status='failed',error=repr(exc));raise
            finally:save(cellout/'RESULT.json',rec);save(out/'RESULT.json',result)
        result['comparison']=compare(m,result['rows']);result['status']='dry_acetate_component_derivatives_complete' if result['comparison']['numerical_checks_passed'] else 'dry_acetate_component_derivatives_failed_numerical_checks'
    except Exception as exc:result.update(status='dry_acetate_component_derivatives_failed',error=repr(exc))
    finally:
        after=resource.getrusage(resource.RUSAGE_CHILDREN);wall=time.monotonic()-tick;result['cost']={'wall_seconds':wall,'process_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':after.ru_utime-children.ru_utime,'child_system_CPU_seconds':after.ru_stime-children.ru_stime,'native_elapsed_seconds':sum(r.get('elapsed_seconds',0) for v in result['rows'].values() for r in v.get('receipts',[])),'allocated_CPU_seconds':wall,'allocated_CPUs':1};result['complete_cells']=sum(r['status']=='complete' for r in result['rows'].values());save(out/'RESULT.json',result)
    return 0 if result['status']=='dry_acetate_component_derivatives_complete' else 1

def admit_result(manifest_path,result_path):
    m,helper,native=load_manifest(manifest_path);r=json.loads(Path(result_path).read_text());m['_path']=str(Path(manifest_path).resolve());release_gate(m,verify(r['release']),Path(result_path).parent)
    require(r['implementation']==pin(__file__) and r['manifest']==pin(manifest_path) and r['analyze']==m['analyze'] and r['error'] is None and r['complete_cells']==r['new_native_calls']==16,'Full16 lineage/admission required')
    require(r['status'] in ('dry_acetate_component_derivatives_complete','dry_acetate_component_derivatives_failed_numerical_checks') and set(r['rows'])==set(m['rows']),'Exact complete16 required')
    for cid,row in m['rows'].items():
        rec=r['rows'][cid];cellout=Path(result_path).parent/'cells'/cid;require(rec==json.loads((cellout/'RESULT.json').read_text()) and rec['status']=='complete','Per-cell persistence differs')
        for k in ('xyz','key','parameters'):verify(rec['input'][k]);require(rec['input'][k]['sha256']==row['input'][k]['sha256'],'Executed physical input differs')
        require(len(rec['receipts'])==1,'One analyze per cell required');q=rec['receipts'][0];require(q['returncode']==0 and q['argv']==[m['analyze']['path'],'system.xyz','E'] and Path(q['cwd'])==cellout and all(v=='1' for v in q['threads'].values()),'Native invocation differs')
        require(verify(q['output'])==cellout/'analyze.out' and parse_energy(cellout/'analyze.out',helper)==rec['parsed'],'Actual raw energy replay differs')
        require(rec['parsed']['interactions']==m['central'][row['metal']]['energy']['interactions'],'Component inventory changed')
    require(compare(m,r['rows'])==r['comparison'],'Saved component derivatives differ');return r

if __name__=='__main__':
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='mode',required=True);p=sub.add_parser('prepare');p.add_argument('--out',type=Path,required=True);r=sub.add_parser('run');r.add_argument('--manifest',type=Path,required=True);r.add_argument('--release',type=Path,required=True);r.add_argument('--out',type=Path,required=True);a=sub.add_parser('admit');a.add_argument('--manifest',type=Path,required=True);a.add_argument('--result',type=Path,required=True);args=ap.parse_args()
    if args.mode=='prepare':prepare(args.out)
    elif args.mode=='run':raise SystemExit(run(args.manifest,args.release,args.out))
    else:print(json.dumps({'status':admit_result(args.manifest,args.result)['status'],'new_native_calls':0}))
