#!/usr/bin/env python3
"""Eight consumed rotation SPs; canonical donor deck/parser, staged one/seven."""
import argparse
import copy
from decimal import Decimal
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('qualified_donor_cc',HERE.parent/'donor_cp_v1/REFERENCE.py')
donor=importlib.util.module_from_spec(spec);spec.loader.exec_module(donor)
b=donor.b;WS=b.ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
PROTOCOL='dry_acetate_correlated_rotation_v1';METALS=('La','Dy')
SHIFTS={'minus_010':-.01,'minus_005':-.005,'plus_005':.005,'plus_010':.01}
ALL=tuple(f'cc__{m}__rotation__{s}__complex' for m in METALS for s in SHIFTS)
SCOUT=ALL[0];REMAINING=ALL[1:];RESOLUTION=.1
PARENT=WS/'reference_chemistry_v1/dry_acetate_force_qualification_v1/preparation_v1/PANEL.json'
PARENT_SHA='aee6e6fc5c78ccf11e3aa86bb5271aa284bc6007a1460edecfadba3c38703ecf'
CANONICAL=WS/'reference_chemistry_v1/donor_cp_v1/prepared_v1/PANEL.json'
CANONICAL_SHA='1a4bca36f2e5f0d0e223b235b0915708a455627444bcb355fc7cbcba2930a701'
DONOR_SHA='f4baaeddce09cb26551e845f05e7f3c1bfff604a84fc02c666908ada2420dea2'
EVIDENCE=HERE.parent/'donor_cp_v1/FINAL_EVIDENCE.json'
EVIDENCE_SHA='3d442f16c25af5c488047899e0edfc9502ce1200d8ef7853720eef9562b04821'
INVENTORY=WS/'reference_chemistry_v1/dry_donor_response_design_v1/SOURCE_INVENTORY.json'
INVENTORY_SHA='6d169c6d1361d2ad4e2b20871c11dd5cdb633149a4a887d4775d59959258bed5'

def read(path):return json.loads(Path(path).read_text())
def verify_tree(value):
    if isinstance(value,dict):
        if set(value)=={'path','sha256'}:b.verify(value)
        else:
            for v in value.values():verify_tree(v)
    elif isinstance(value,list):
        for v in value:verify_tree(v)

def source_gates():
    b.check(b.pin(donor.__file__)['sha256']==DONOR_SHA,'Canonical donor source changed')
    for path,sha in [(PARENT,PARENT_SHA),(CANONICAL,CANONICAL_SHA),(EVIDENCE,EVIDENCE_SHA),(INVENTORY,INVENTORY_SHA)]:
        b.check(b.pin(path)['sha256']==sha,'Frozen authority changed: '+str(path))
    old=read(CANONICAL)
    for key,mod in [('base_helper',b),('canonical_helper',donor.cp),('dft_helper',donor.dft)]:
        b.check(b.pin(mod.__file__)==old[key],'Canonical helper changed: '+key)
    return old

def deck(panel,cell,ranks,maxcore):return donor.deck(panel,cell,ranks,maxcore)

def parse(text,cell,returncode):
    result=donor.parse(text,cell,returncode);native=text.split('****END OF INPUT****')[-1]
    tokens={}
    for name,pattern in [('CCSD_T',r'FINAL SINGLE POINT ENERGY\s+'),('CCSD',r'^E\(CCSD\)\s*\.{2,}\s*'),('triples',r'^Triples Correction \(T\)\s*\.{2,}\s*')]:
        matches=re.findall(pattern+'('+donor.F+')',native,re.M)
        b.check(len(matches)==1 and math.isfinite(donor.cp.number(matches[0])),'Missing/ambiguous finite '+name+' energy')
        tokens[name]=matches[0].replace('D','E').replace('d','e')
    result['energy_decimal_tokens_hartree']=tokens
    result['checks']['CCSD_plus_triples_closure']=abs(float(Decimal(tokens['CCSD_T'])-Decimal(tokens['CCSD'])-Decimal(tokens['triples'])))<1e-8
    result['status']='complete' if all(result['checks'].values()) else 'failed_admission'
    return result

def prepare(out):
    tick=time.monotonic();cpu=time.process_time();old=source_gates();parent=read(PARENT);inv=read(INVENTORY);ev=read(EVIDENCE)
    b.check(parent['derivative_resolution_kcal_mol_A']==RESOLUTION and parent['steps_A']==[.01,.005] and parent['rotation_length_A']==1.,'Original direction/resolution differs')
    out.mkdir(parents=True,exist_ok=False);(out/'inputs').mkdir()
    authorities={k:inv['authorities'][k] for k in ['reference_derivative_result','reference_actual_review','model_manifest','model_derivative_result','selected_donor_review']}
    verify_tree(authorities);cells=[]
    for m in METALS:
        b.check(old['assets'][m]['basis']['sha256']==parent['assets'][m]['basis']['sha256'],'Canonical/PBE0 orbital or ECP asset differs')
        for s,q in SHIFTS.items():
            qid=f'pbe0__{m}__rotation__{s}__complex';original=next(c for c in parent['cells'] if c['cell_id']==qid)
            c=copy.deepcopy(original);c.pop('template_input');c.update(cell_id=qid.replace('pbe0__','cc__',1),method='cc',d4=False,donor='acetate_mono',parent_cell_id=qid,parent_cell=original)
            b.check(c['q_A']==q and inv['rows'][qid]['model_input']['xyz'],'Original rotation join differs');cells.append(c)
    costs=[r for r in ev['native_rows'] if r['cell_id'] in ('cc__acetate_bi__La__complex','cc__acetate_bi__Dy__complex')]
    b.check(len(costs)==2,'Two measured canonical size controls required')
    panel={'protocol':PROTOCOL,'adapter':b.pin(__file__),'specification':b.pin(HERE/'SPEC.md'),'parent_panel':b.pin(PARENT),'canonical_panel':b.pin(CANONICAL),'source_inventory':b.pin(INVENTORY),'method_evidence':b.pin(EVIDENCE),'canonical_review':ev['strict_review'],'helpers':{'donor':b.pin(donor.__file__),'base':b.pin(b.__file__),'canonical':b.pin(donor.cp.__file__),'dft':b.pin(donor.dft.__file__)},'orca':old['orca'],'assets':old['assets'],'cells':cells,'comparison_authorities':authorities,'cost_controls':costs,'estimated_allocated_CPU_seconds':4*sum(r['allocated_cpu_seconds'] for r in costs),'scout_cell_id':SCOUT,'remaining_cell_ids':list(REMAINING),'maximum_native_calls':8,'ranks':4,'nodes':1,'steps_A':[.01,.005],'derivative_resolution_kcal_mol_A':RESOLUTION,'physical_state_qualified':False,'method_accuracy_qualified':False,'new_native_calls':0}
    for c in cells:
        filename=out/'inputs'/(c['cell_id']+'.inp');filename.write_text(deck(panel,c,4,3000));c['template_input']=b.pin(filename)
    pp=b.write(out/'PANEL.json',panel);check=validate(pp['path']);check['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
    pf=b.write(out/'PREFLIGHT.json',check);print(json.dumps({'panel':pp,'preflight':pf,'status':check['status'],'new_native_calls':0},indent=2))

def load(path):
    source_gates();p=read(path);verify_tree(p)
    b.check(p['protocol']==PROTOCOL and p['adapter']==b.pin(__file__) and p['specification']==b.pin(HERE/'SPEC.md'),'Source/specification differs')
    b.check(p['parent_panel']==b.pin(PARENT) and p['canonical_panel']==b.pin(CANONICAL) and p['source_inventory']==b.pin(INVENTORY),'Parent bindings differ')
    b.check(tuple(c['cell_id'] for c in p['cells'])==ALL and p['maximum_native_calls']==8 and p['scout_cell_id']==SCOUT and p['remaining_cell_ids']==list(REMAINING),'Exact eight/stages differ')
    b.check(p['ranks']==4 and p['nodes']==1 and p['steps_A']==[.01,.005] and p['derivative_resolution_kcal_mol_A']==RESOLUTION,'Resource/numerical scope differs')
    p['_path']=str(Path(path).resolve());return p

def validate(path):
    p=load(path);parent=read(PARENT);old=read(CANONICAL);records=[]
    for c in p['cells']:
        original=next(x for x in parent['cells'] if x['cell_id']==c['parent_cell_id']);b.check(c['parent_cell']==original,'Exact original cell differs')
        for k,v in original.items():
            if k not in ('cell_id','method','d4','template_input'):b.check(c[k]==v,'Changed geometry/state field '+k)
        symbols,x=b.xyz_read(b.verify(c['source_geometry']));b.check(symbols==c['symbols'] and x==c['positions_A'],'Literal XYZ changed')
        b.check(c['method']=='cc' and not c['d4'] and c['role']=='complex' and c['charge']==2 and c['multiplicity']==1 and c['expected_electrons']==c['expected_correlated_electrons']==40 and c['expected_basis_functions']==189 and not any(c['ghost_mask']),'Chemistry differs')
        raw=deck(p,c,4,3000);b.check(raw==b.verify(c['template_input']).read_text(),'Prepared deck differs')
        peer=next(x for x in old['cells'] if x['cell_id']==f'cc__acetate_bi__{c["metal"]}__complex')
        b.check(raw.split('* xyz ',1)[0]==deck(old,peer,4,3000).split('* xyz ',1)[0],'Original canonical recipe changed')
        lines=raw.split('* xyz 2 1\n',1)[1].splitlines()[:-1]
        b.check([s.split()[0] for s in lines]==symbols and [[float(v) for v in s.split()[1:]] for s in lines]==x,'Rendered coordinates differ')
        b.check(len(lines)==8 and all(k not in raw for k in ('EnGrad','RIJCOSX','D4','DLPNO','DelECP')),'Unexpected scientific directive')
        records.append({'cell_id':c['cell_id'],'exact_original_xyz':True,'unchanged_canonical_recipe_prefix':True,'real_centers':8,'electrons':40,'orbital_functions':189})
    for s in SHIFTS:
        pair=[next(c for c in p['cells'] if c['cell_id']==f'cc__{m}__rotation__{s}__complex') for m in METALS]
        b.check(pair[0]['positions_A']==pair[1]['positions_A'],'Matched metal coordinates differ')
    return {'status':'passed_static_preflight','implementation':b.pin(__file__),'panel':b.pin(path),'records':records,'new_native_calls':0}

def release_gate(p,release_pin,cell_id):
    r=read(b.verify(release_pin));stage=r['stage'];ids=(SCOUT,) if stage=='scout' else REMAINING
    pp=b.pin(p['_path']);b.check(stage in ('scout','remaining') and r['status']=='released_dry_acetate_correlated_rotation_'+stage,'Unreleased stage')
    b.check(r['implementation']==b.pin(__file__) and r['panel']==pp and r['cell_ids']==list(ids) and cell_id in ids,'Release/source/inventory differs')
    b.check(r['maximum_total_native_calls']==8 and r['maximum_stage_calls']==len(ids) and r['ranks']==4 and r['nodes']==1,'Release bounds differ')
    root=Path(r['execution_root']).resolve();expected={cid:str(root/'tasks'/cid/'execution') for cid in ids}
    b.check(r['task_output_directories']==expected,'One immutable directory per endpoint required')
    review=read(b.verify(r['source_review']));b.check(review['status']=='passed_dry_acetate_correlated_rotation_source_review' and review['implementation']==b.pin(__file__) and review['panel']==pp,'Source review join differs')
    if stage=='remaining':
        scout=admit(p,r['scout_result']['path']);sr=read(b.verify(r['scout_review']))
        b.check(b.pin(r['scout_result']['path'])==r['scout_result'] and scout['cell_id']==SCOUT,'Wrong scout')
        b.check(sr['status']=='passed_dry_acetate_correlated_rotation_scout_review' and sr['result']==r['scout_result'] and sr['implementation']==b.pin(__file__) and sr['panel']==pp and sr['cost_reviewed'] is True,'Actual scout/cost review missing')
        first=read(b.verify(scout['release']));b.check(first['execution_root']==r['execution_root'],'Continuation moved execution root')
    return r

def run(path,cell_id,release_path,out):
    p=load(path);rp=b.pin(release_path);r=release_gate(p,rp,cell_id);c=next(c for c in p['cells'] if c['cell_id']==cell_id)
    b.check(str(out.resolve())==r['task_output_directories'][cell_id],'Output differs from release')
    b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Native calls require Slurm')
    ranks=int(os.environ.get('SLURM_NTASKS','0'));nodes=int(os.environ.get('SLURM_JOB_NUM_NODES',os.environ.get('SLURM_NNODES','0')))
    b.check(ranks==4 and nodes==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1,'Four single-thread ranks on one node required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or ranks*int(os.environ.get('SLURM_MEM_PER_CPU','0'));b.check(mem>0,'Actual memory grant missing');maxcore=math.floor(.75*mem/ranks)
    out.mkdir(parents=True,exist_ok=False);(out/'endpoint.inp').write_text(deck(p,c,ranks,maxcore))
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):env[key]='1'
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE');env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
    cmd=[p['orca']['path'],'endpoint.inp','--bind-to none'];pp=b.pin(path)
    start={'adapter':b.pin(__file__),'cell':c,'panel':pp,'release':rp,'input':b.pin(out/'endpoint.inp'),'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'ranks':ranks,'nodes':nodes,'memory_MiB':mem,'maxcore_MiB':maxcore,'command':cmd,'utc':b.now()}
    sp=b.write(out/'START.json',start);tick=time.monotonic();returncode=None
    try:
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:proc=subprocess.run(cmd,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False)
        returncode=proc.returncode;result=parse((out/'endpoint.out').read_text(errors='replace'),c,returncode)
    except Exception as exc:result={'status':'failed_admission','error':repr(exc)}
    elapsed=time.monotonic()-tick
    result.update(cell_id=cell_id,cell=c,panel=pp,adapter=b.pin(__file__),release=rp,start=sp,input=start['input'],output=b.pin(out/'endpoint.out'),stderr=b.pin(out/'endpoint.err'),returncode=returncode,job_id=start['job_id'],elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*ranks,new_native_calls=1,finished_utc=b.now())
    b.write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'cell_id':cell_id,'seconds':elapsed}));return 0 if result['status']=='complete' else 1

def admit(p,filename):
    r=read(filename);cid=r['cell_id'];c=next(c for c in p['cells'] if c['cell_id']==cid);pp=b.pin(p['_path'])
    b.check(r['status']=='complete' and r['adapter']==b.pin(__file__) and r['panel']==pp and r['cell']==c and r['new_native_calls']==1,'Native result identity/status differs')
    release=release_gate(p,r['release'],cid);directory=Path(filename).resolve().parent
    b.check(str(directory)==release['task_output_directories'][cid],'Native directory differs')
    for key,name in [('input','endpoint.inp'),('output','endpoint.out'),('stderr','endpoint.err'),('start','START.json')]:
        b.check(r[key]==b.pin(directory/name),'Native artifact join differs: '+key)
    s=read(b.verify(r['start']));b.check(s['adapter']==r['adapter'] and s['cell']==c and s['panel']==pp and s['release']==r['release'] and s['input']==r['input'] and s['job_id']==r['job_id'],'START identity differs')
    b.check(s['ranks']==4 and s['nodes']==1 and s['maxcore_MiB']==math.floor(.75*s['memory_MiB']/4) and s['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'],'Native resources/command differ')
    b.check(b.verify(r['input']).read_text()==deck(p,c,4,s['maxcore_MiB']),'Executed input differs')
    replay=parse(b.verify(r['output']).read_text(),c,r['returncode']);b.check(replay['status']=='complete','Raw canonical admission failed')
    for k,v in replay.items():b.check(r[k]==v,'Stored/raw parse differs: '+k)
    b.check(math.isfinite(r['elapsed_seconds']) and r['elapsed_seconds']>=0 and r['allocated_cpu_seconds']==4*r['elapsed_seconds'],'Native cost differs')
    return r

def half_token(token):return Decimal('0.5')*Decimal(10)**Decimal(token).as_tuple().exponent

def score(rows):
    groups={m:None for m in (*METALS,'paired')};checks={m:{'passed':None,'difference_kcal_mol_A':None,'resolution_kcal_mol_A':RESOLUTION} for m in groups}
    if set(rows)!=set(ALL):return {'groups':groups,'checks':checks,'declared_checks':3,'all_checks_passed':False}
    b.check(all(Decimal(token).is_finite() for row in rows.values() for token in row['energy_decimal_tokens_hartree'].values()),'Nonfinite score input')
    ha=Decimal(str(b.HA))
    for m in METALS:
        steps={}
        for q,label in [(.01,'010'),(.005,'005')]:
            plus=rows[f'cc__{m}__rotation__plus_{label}__complex']['energy_decimal_tokens_hartree'];minus=rows[f'cc__{m}__rotation__minus_{label}__complex']['energy_decimal_tokens_hartree'];den=2*Decimal(str(q))
            steps[str(q)]={k:{'derivative_kcal_mol_A':float((Decimal(plus[k])-Decimal(minus[k]))*ha/den),'print_bound_kcal_mol_A':float((half_token(plus[k])+half_token(minus[k]))*ha/den)} for k in ('CCSD_T','CCSD','triples')}
        groups[m]={'steps':steps}
    groups['paired']={'steps':{q:{k:{'derivative_kcal_mol_A':groups['Dy']['steps'][q][k]['derivative_kcal_mol_A']-groups['La']['steps'][q][k]['derivative_kcal_mol_A'],'print_bound_kcal_mol_A':groups['Dy']['steps'][q][k]['print_bound_kcal_mol_A']+groups['La']['steps'][q][k]['print_bound_kcal_mol_A']} for k in ('CCSD_T','CCSD','triples')} for q in ('0.01','0.005')}}
    for m,g in groups.items():
        delta=g['steps']['0.005']['CCSD_T']['derivative_kcal_mol_A']-g['steps']['0.01']['CCSD_T']['derivative_kcal_mol_A'];checks[m].update(difference_kcal_mol_A=delta,passed=abs(delta)<=RESOLUTION)
        g['smaller_minus_larger_step_kcal_mol_A']=delta
    return {'groups':groups,'checks':checks,'declared_checks':3,'all_checks_passed':all(v['passed'] is True for v in checks.values())}

def collect(path,inventory_path,out):
    tick=time.monotonic();cpu=time.process_time();p=load(path);inventory=read(inventory_path);b.check(set(inventory)==set(ALL),'Exact eight-slot inventory required')
    rows={};slots={cid:{'status':'unavailable'} for cid in ALL};costs=[]
    for cid,filename in inventory.items():
        if filename is None or not Path(filename).is_file():continue
        try:
            raw=read(filename)
            if raw.get('new_native_calls')==1:costs.append({'cell_id':cid,'job_id':raw.get('job_id'),'elapsed_seconds':raw.get('elapsed_seconds'),'allocated_cpu_seconds':raw.get('allocated_cpu_seconds'),'reported_status':raw.get('status')})
            admitted=admit(p,filename);b.check(admitted['cell_id']==cid,'Result inventory mismatch');rows[cid]=admitted;slots[cid]={'status':'complete','result':b.pin(filename)}
        except Exception as exc:slots[cid]={'status':'unavailable','result':b.pin(filename),'error':repr(exc)}
    scored=score(rows);qm=read(b.verify(p['comparison_authorities']['reference_derivative_result']));model=read(b.verify(p['comparison_authorities']['model_derivative_result']))
    b.check(qm['admitted_cells']==16 and qm['stage_numerical_checks_passed'] is True and model['complete_cells']==16 and model['error'] is None,'Saved comparison unavailable')
    comparison={}
    for m in (*METALS,'paired'):
        comparison[m]={}
        for q in ('0.01','0.005'):
            cc=None if scored['groups'][m] is None else scored['groups'][m]['steps'][q]['CCSD_T']['derivative_kcal_mol_A']
            ref=qm['scoring']['directions']['rotation'][m]['derivatives_kcal_mol_A'][q];mm=model['comparison']['directions']['rotation'][m]['steps'][q]['derivatives_kcal_mol_A']
            comparison[m][q]={'canonical_CCSD_T':cc,'PBE0_electronic':ref['electronic'],'PBE0_D4':ref['PBE0_D4'],'D4_only':ref['D4'],'fixed_model':mm['Total'],'fixed_model_components':mm,'model_minus_CC':None if cc is None else mm['Total']-cc,'PBE0_electronic_minus_CC':None if cc is None else ref['electronic']-cc,'PBE0_D4_minus_CC':None if cc is None else ref['PBE0_D4']-cc,'units':'kcal/mol/Angstrom','cross_method_accuracy_gate':None}
    status='canonical_rotation_complete' if len(rows)==8 else 'canonical_rotation_incomplete'
    result={'status':status,'implementation':b.pin(__file__),'panel':b.pin(path),'inventory':b.pin(inventory_path),'slots':slots,'admitted_cells':len(rows),'declared_cells':8,'scoring':scored,'comparison':comparison,'comparison_authorities':p['comparison_authorities'],'native_admissions':rows,'native_costs':costs,'native_allocated_CPU_seconds':sum(x['allocated_cpu_seconds'] for x in costs if isinstance(x['allocated_cpu_seconds'],(int,float)) and math.isfinite(x['allocated_cpu_seconds'])),'new_native_calls':0,'physical_state_qualified':False,'method_accuracy_qualified':False,'CP_force_qualified':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    out.mkdir(parents=True,exist_ok=False);rp=b.write(out/'RESULT.json',result)
    lines=['# Consumed dry-acetate canonical rotation','',f'{len(rows)}/8 admitted; all three two-step numerical checks pass: {scored["all_checks_passed"]}. No new native calls in collection.','','| Comparison | q (Å) | CCSD(T) | PBE0 electronic | PBE0-D4 | Fixed model |','|---|---:|---:|---:|---:|---:|']
    for m,steps in comparison.items():
        for q,v in steps.items():lines.append('| '+m+' | '+q+' | '+('unavailable' if v['canonical_CCSD_T'] is None else f'{v["canonical_CCSD_T"]:.9f}')+f' | {v["PBE0_electronic"]:.9f} | {v["PBE0_D4"]:.9f} | {v["fixed_model"]:.9f} |')
    lines+=['','Slopes are kcal/mol/Å, q=1Å×rotation angle. The0.1 resolution applies only to CC two-step consistency. Cross-method differences are descriptive. D4 is separate and never added to CC. Frozen-f/state/basis/diffuse/single-reference limitations and raw SCF/CC signals remain in RESULT; this is not a CP-corrected force or physical-truth claim.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'result':rp,'status':status,'numerical_checks_passed':scored['all_checks_passed']},indent=2))

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--cell',required=True);a.add_argument('--release',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.action=='prepare':prepare(a.out)
    elif a.action=='run':return run(a.panel,a.cell,a.release,a.out)
    else:collect(a.panel,a.inventory,a.out)

if __name__=='__main__':sys.exit(main())
