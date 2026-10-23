#!/usr/bin/env python3
"""Finite saved-operator response diagnostic. No native, optimization or selection."""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
SCOUT='aquo9_central_La__dimer_01'

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):require(pin(q['path'])==q,'Changed pin: '+q['path']);return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def module(q,name):
    s=importlib.util.spec_from_file_location(name,verify(q));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def arithmetic(n,*values):return 128*n*np.finfo(float).eps*max([1.]+[float(np.max(abs(np.asarray(v)))) for v in values])
def at(value,path):
    for k in path:value=value[k]
    return value
def geometry(types,positions):
    xyz=np.asarray(positions,float);metal=[i for i,t in enumerate(types) if t in (500,501)];oxygen=[i for i,t in enumerate(types) if t==36]
    require(len(metal)==1 and set(types)<={36,37,500,501},'Only one metal and intact water types')
    mi=metal[0];drive=np.zeros_like(xyz);factors=[]
    for i in oxygen:
        v=xyz[i]-xyz[mi];r=float(np.linalg.norm(v));require(r>0 and math.isfinite(r),'Invalid metal/O separation')
        t=2*r;f=((1+t+t*t/3)*math.exp(-t)/((1+5+25/3)*math.exp(-5)))**2
        factors.append(f);drive[i]=f*3*v/(r**3)
    s=math.fsum(factors);return {'metal_index':mi,'factors':factors,'s':s,'g':2*s/(1+s),'drive_per_a':drive.reshape(-1).tolist()}

def prepare(cfg,out):
    tick,cpu=time.monotonic(),time.process_time()
    require(cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md'),'Source/specification binding')
    inv=read(cfg['inventory']);target=read(cfg['targets']);field=module(cfg['field_helper'],'screen_field');moment=module(cfg['moment_helper'],'screen_moment')
    require(inv['artifacts']['target_inventory']==cfg['targets'] and len(inv['rows'])==54 and len(inv['cached_metal_absent'])==46,'Frozen target/fixture denominator')
    loaded={}
    def artifact(name):
        if name not in loaded:loaded[name]=read(inv['artifacts'][name])
        return loaded[name]
    def coordinates(item,frame_source=None):
        q={k:inv['artifacts'][v] for k,v in item['input'].items()};rr=[s.split() for s in verify(q['xyz']).read_text().splitlines()[1:] if s.strip()]
        require([int(s[0]) for s in rr]==list(range(1,len(rr)+1)),'Literal atom order differs')
        mp=moment.multipoles(verify(q['parameters']).read_text(),verify(q['key']).read_text())
        if frame_source is None:frames=[moment.frame_signature(rr,i,mp) for i in range(len(rr))]
        else:
            frames=frame_source['frame_signatures']
            require(all(verify(q[k]).read_bytes()==verify(frame_source['input'][k]).read_bytes() for k in q),'Static frame/input authority differs')
            require([f['type'] for f in frames]==[int(s[5]) for s in rr] and [f['position_A'] for f in frames]==[[float(x) for x in s[2:5]] for s in rr],'Static literal frame map differs')
            require(all(f['permanent_record']==mp[f['type']] for f in frames),'Static permanent multipole authority differs')
        p=field.permanent_moment({'frame_signatures':frames})
        return {'input':q,'atoms':len(rr),'types':[int(s[5]) for s in rr],'positions_A':[[float(x) for x in s[2:5]] for s in rr],
                'permanent_lab_Debye':p,'frame_signatures':frames}
    rows={}
    for cid,item in inv['rows'].items():
        row=coordinates(item);require(row['atoms']==item['atoms'],'Original fixture size differs')
        row.update(metal=item['metal'],geometry=geometry(row['types'],row['positions_A']),baseline_energy_kcal_mol=item['cached'][0]['energy_kcal_mol'])
        rows[cid]=row
    static={}
    for cid,item in inv['cached_metal_absent'].items():
        original=artifact('moment_manifest')['rows'].get(cid) or artifact('trimer_manifest')['model_rows'].get(cid)
        require(original is not None,'Static parent frame authority unavailable')
        row=coordinates(item,original);require(not set(row['types'])&{500,501},'Metal-free cache contains metal')
        src=item['source'];value=at(artifact(src['artifact']),src['json_path'])
        for k in ('analyze_output','induced_output'):verify(value[k])
        mu=np.asarray(row['permanent_lab_Debye'])+np.sum(value['induced_Debye'],axis=0)
        require(np.max(abs(mu-np.asarray(value['total_lab_Debye'])))<=.0005+row['atoms']*.5e-6+2e-9,'Static permanent/common-origin moment replay')
        static[cid]={'energy_kcal_mol':value['energy_kcal_mol'],'mu_Debye':mu.tolist(),'mu_print_bound':row['atoms']*.5e-6+2e-9,'source':src,'input':row['input'],'permanent_lab_Debye':row['permanent_lab_Debye']}
    gmax=max(r['geometry']['g'] for r in rows.values());limit=math.log(2)/gmax
    require(gmax>0 and math.isfinite(limit),'Geometry-only global L')
    grid=[{'id':f'a{ia}_l{il}','a':a,'l':l,'a_index':ia,'l_index':il} for ia,a in enumerate((-.1,0.,.1)) for il,l in enumerate((-limit,0.,limit))]
    m={'status':'prepared_saved_screen_unbound','implementation':pin(__file__),'configuration':cfg,'specification':cfg['specification'],'inventory':cfg['inventory'],'targets':cfg['targets'],
       'rows':rows,'static':static,'geometry_gmax':gmax,'L':limit,'grid_per_metal':grid,'metals':['La','Dy'],'field_aliases':inv['field_aliases'],
       'units':artifact('field_manifest')['conversion'],
       'counts':{'fixtures':54,'static':46,'grid_per_metal':9,'matrix_states':486,'paired_compositions':81,'force_checks_unevaluated':15},
       'new_native_calls':0,'screen_executed':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    save(out/'MANIFEST.json',m)
    save(out/'CONFIG_RUN_TEMPLATE.json',{'implementation':pin(__file__),'specification':cfg['specification'],'manifest':pin(out/'MANIFEST.json'),'full_result':None,'full_review':None,'threads':1})
    return m

def kernel(A,ed,ep,alpha,geo,a,l,units):
    """Stationary dual energy, conjugate induced moment and uniform-z column."""
    n=len(alpha)//3;mi=geo['metal_index'];Ht=np.zeros_like(A);sigma=l*geo['g']
    require(abs(sigma)<=math.log(2)+2e-15,'Prospective local response domain')
    idx=slice(3*mi,3*mi+3);Ht[idx,idx]=np.eye(3)*math.expm1(sigma)/alpha[3*mi]
    At=A+Ht;norm=np.sqrt(alpha)[:,None]*At*np.sqrt(alpha)[None,:]
    require(np.isfinite(norm).all(),'Nonfinite candidate matrix')
    eigen=np.linalg.eigvalsh(norm);spectrum={'normalized_min':float(eigen[0]),'normalized_max':float(eigen[-1]),'normalized_condition':float(eigen[-1]/eigen[0]) if eigen[0]>0 else None}
    if not eigen[0]>1e-6*eigen[-1]>0:return {'status':'unstable','spectrum':spectrum,'sigma':sigma}
    C=np.linalg.cholesky(At);h=a*np.asarray(geo['drive_per_a']);d=ed+h;p=ep+h
    Wz=np.tile([0.,0.,1.],n)
    sol=np.linalg.solve(C.T,np.linalg.solve(C,np.column_stack((d,p,Wz))))
    u,v,x=sol.T;K=units['electric'];B=units['bohr_A'];H=units['hartree_kcal_mol'];D=units['debye_per_eA']
    residual=float(np.max(abs(At@sol-np.column_stack((d,p,Wz)))));bound=arithmetic(3*n,At,sol,d,p)*max(1.,float(np.max(abs(At))),float(np.max(abs(sol))))
    require(residual<=bound,'Candidate dense solve arithmetic failed')
    return {'status':'complete','polarization_energy':-.5*K*float(u@p),'induced_conjugate_Debye':(D*(u+v).reshape(n,3).sum(axis=0)/2).tolist(),
            'chi_column_a0_cubed':(H/(K*B*B)*x.reshape(n,3).sum(axis=0)).tolist(),'u':u.tolist(),'v':v.tolist(),
            'spectrum':spectrum,'sigma':sigma,'solve_residual_max':residual,'solve_arithmetic_bound':bound}

def matrix_context(cfg,m):
    require(cfg['implementation']==pin(__file__) and cfg['specification']==pin(HERE/'SPEC.md')==m['specification'] and m['implementation']==cfg['implementation'],'Runtime source/specification')
    require(cfg['threads']==1 and all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
    r=read(cfg['full_result']);review=read(cfg['full_review']);pcfg=m['configuration'];scout=module(pcfg['scout_helper'],'screen_scout')
    require(review['status']=='passed_baseline_matrix_inventory_review' and review['result']==cfg['full_result'] and review['implementation']==r['implementation']==pcfg['runner'],'Full actual54 independent admission')
    require(read(review['config'])==r['config'] and review['source_manifest']==r['config']['export_source_manifest'] and review['manifest']==r['manifest'],'Full actual provenance joins')
    require(r['status']=='baseline_matrix_inventory_complete' and r['complete_fixtures']==54 and r['admitted_snapshots']==108 and r['native_program_invocations']==53 and r['error'] is None and r['global_stop'] is False,'Full54 terminal denominator')
    rm=read(r['manifest']);require(rm['inventory']==m['inventory'] and rm['units']==m['units'] and set(r['rows'])==set(m['rows']),'Exact54 inventory/units')
    matrices={};baseline={};replays={}
    for cid,row in m['rows'].items():
        rec=r['rows'][cid]
        if cid==SCOUT:
            require(rec['inherited'] and rec['admission']['status']=='admitted_response_matrix_scout','Inherited scout admission')
            native=rec['admission']['on_row'];qualification=rec['admission']['qualification']
        else:native=rec;qualification=rec['numerical']
        require(rec['status']=='complete' and qualification['passed'],'Incomplete fixture qualification')
        rr=rm['rows'][cid];require(rr['source_input']==row['input'] and rr['types']==row['types'] and rr['positions_A']==row['positions_A'],'Actual literal input mapping')
        parsed=[scout.parse_export(q,rr,m['units']) for q in native['exports']]
        for p,saved in zip(parsed,qualification['exports']):
            require(all(p[k]==saved[k] for k in ('source','sequence','sites','vectors','A','checks','peek')),'Raw qualified matrix changed')
        p=parsed[0];A=np.asarray(p['A']);ed,ep,u,v=np.asarray(p['vectors']).T;alpha=np.repeat([s['alpha_A3'] for s in p['sites']],3)
        require([s['atom'] for s in p['sites']]==list(range(1,row['atoms']+1)),'Literal site order')
        require(np.max(abs(A-A.T))<=scout.arithmetic(3*row['atoms'],A) and np.array_equal(ed,ep),'Baseline symmetric intact-water dual convention')
        zero=kernel(A,ed,ep,alpha,row['geometry'],0.,0.,m['units']);require(zero['status']=='complete','Baseline matrix unstable')
        mom=native['moments'];bound_u=qualification['exports'][0]['numerical_bounds']['native_vs_dense_u'];bound_v=qualification['exports'][0]['numerical_bounds']['native_vs_dense_v']
        D=m['units']['debye_per_eA'];n=row['atoms'];mu=np.asarray(row['permanent_lab_Debye'])+np.asarray(zero['induced_conjugate_Debye'])
        mbound=n*.5e-6+D*math.sqrt(n)*(bound_u+bound_v)/2+2e-9
        require(np.max(abs(mu-np.asarray(mom['total_lab_Debye'])))<=.0005+mbound,'Baseline total laboratory moment replay')
        pol=mom['components'].get('Polarization',{'energy_kcal_mol':0.})['energy_kcal_mol']
        ebound=5e-9+.5*m['units']['electric']*np.linalg.norm(ep)*bound_u+arithmetic(3*n,pol)
        require(abs(zero['polarization_energy']-pol)<=ebound,'Unanchored baseline polarization replay')
        require(abs(mom['energy_kcal_mol']-row['baseline_energy_kcal_mol'])<=1e-8,'Baseline original total energy replay')
        baseline[cid]={'energy_kcal_mol':mom['energy_kcal_mol'],'mu_Debye':mu.tolist(),'mu_print_bound':mbound,'chi_column_a0_cubed':zero['chi_column_a0_cubed']}
        matrices[cid]={'A':A,'ed':ed,'ep':ep,'alpha':alpha,'zero':zero}
        replays[cid]={'energy_difference_kcal_mol':mom['energy_kcal_mol']-row['baseline_energy_kcal_mol'],'polarization_difference_kcal_mol':zero['polarization_energy']-pol,'polarization_bound':ebound,'total_moment_print_bound_Debye':.0005+mbound,'raw_exports':native['exports'],'spectrum':zero['spectrum']}
    return r,review,matrices,baseline,replays

def combine(values,coefficients,key):
    return sum((float(c)*np.asarray(values[cid][key],float) for cid,c in coefficients.items()),np.asarray(0.))
def response(value,reference,baseline):
    x=np.asarray(value,float);q=np.asarray(reference,float);b=np.asarray(baseline,float)
    return {'model':x.tolist(),'reference':q.tolist(),'model_minus_reference':(x-q).tolist(),'change_from_baseline':(x-b).tolist()}

def observables(values,target,aliases):
    pairs={};parents={};trimers={};field={}
    group=target['groups']['pair_and_collective_dipoles']
    for row in group['pair_cycles']:
        value=combine(values,row['model_terms'],'mu_Debye')
        pairs[row['id']]=response(value,row['reference_pair_mu_Debye'],row['selected17676_pair_mu_Debye'])
    for parent,targets in group['parents'].items():
        metal=parent.rsplit('_',1)[1];terms={parent:1,'ion_'+metal:-1,parent+'__shell':-1}
        C=combine(values,terms,'mu_Debye');P=sum((np.asarray(pairs[q['id']]['model']) for q in group['pair_cycles'] if q['parent']==parent),np.zeros(3))
        parents[parent]={k:response(v,targets[k]['reference_mu_Debye'],targets[k]['selected17676_model_mu_Debye']) for k,v in [('C',C),('P',P),('B',C-P)]}
    for cid,row in target['groups']['trimer_energy_and_dipoles']['cycles'].items():
        trimers[cid]={'energy_kcal_mol':response(combine(values,row['model_terms'],'energy_kcal_mol'),row['reference_electronic_E3_kcal_mol'],row['model_E3_kcal_mol']),
                      'mu_Debye':response(combine(values,row['model_terms'],'mu_Debye'),row['reference_mu3_Debye'],row['model_mu3_Debye'])}
    for state,item in aliases.items():field[state]=values[item['fixture']]['chi_column_a0_cubed'][2]
    for metal in ('La','Dy'):field[metal]=field[metal+'_trimer']-field[metal+'_dimer']
    field['Dy_minus_La']=field['Dy']-field['La']
    fields=[dict(row,candidate_a0_cubed=field[row['quantity']],candidate_minus_reference_a0_cubed=field[row['quantity']]-row['reference_a0_cubed'],change_from_baseline_a0_cubed=field[row['quantity']]-row['model_a0_cubed']) for row in target['groups']['conditional_susceptibility']['reported_rows']]
    energy=[]
    for row in target['groups']['original_energy_force']['energy_checks']:
        v=float(combine(values,row['coefficients'],'energy_kcal_mol'));err=v-row['reference_kcal_mol']
        energy.append(dict(row,candidate_kcal_mol=v,error_kcal_mol=err,passed=abs(err)<=row['absolute_error_limit_kcal_mol']))
    return {'original_energy_checks':energy,'pair_dipoles':pairs,'CPB_dipoles':parents,'trimer_E3_mu3':trimers,'field_rows':fields,
            'force_checks':{'status':'unevaluated','count':15},'physical_response_acceptance':None,'aggregate_loss':None,'selection':None}

def zero_replay(observed,target,baseline):
    for row in observed['original_energy_checks']:
        require(abs(row['candidate_kcal_mol']-row['selected17676_model_kcal_mol'])<=sum(abs(c) for c in row['coefficients'].values())*1e-8+2e-10,'Zero original energy cycle')
    group=target['groups']['pair_and_collective_dipoles']
    for row in group['pair_cycles']:
        b=sum(abs(c)*baseline[cid]['mu_print_bound'] for cid,c in row['model_terms'].items())
        require(np.max(abs(observed['pair_dipoles'][row['id']]['change_from_baseline']))<=b+2e-9,'Zero pair dipole cycle')
    for parent,row in observed['CPB_dipoles'].items():
        for k,v in row.items():require(np.max(abs(v['change_from_baseline']))<=group['parents'][parent][k]['model_print_bound_Debye_per_component']+sum(x['mu_print_bound'] for cid,x in baseline.items() if cid.startswith(parent))+2e-9,'Zero CPB cycle')
    for cid,row in target['groups']['trimer_energy_and_dipoles']['cycles'].items():
        o=observed['trimer_E3_mu3'][cid];b=sum(abs(c)*baseline[x]['mu_print_bound'] for x,c in row['model_terms'].items())
        require(abs(o['energy_kcal_mol']['change_from_baseline'])<=sum(abs(c) for c in row['model_terms'].values())*1e-8+2e-10,'Zero E3 energy')
        require(np.max(abs(o['mu_Debye']['change_from_baseline']))<=b+2e-9,'Zero mu3')
    require(len(observed['original_energy_checks'])==10 and len(observed['pair_dipoles'])==34 and len(observed['CPB_dipoles'])==4 and len(observed['trimer_E3_mu3'])==8 and len(observed['field_rows'])==28,'Full observable denominator')

def run(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();m=read(cfg['manifest']);target=read(m['targets']);read(m['inventory'])
    r,review,matrices,base,replays=matrix_context(cfg,m);base.update(m['static'])
    zero=observables(base,target,m['field_aliases']);zero_replay(zero,target,base)
    # All four full-precision columns must join the independent full54 field replay.
    fields=review['field_alias_consistency']
    for state,alias in m['field_aliases'].items():
        q=fields[state];require(q['passed'],'Independent baseline field qualification unavailable')
        a=np.asarray(base[alias['fixture']]['chi_column_a0_cubed']);b=np.asarray(q['susceptibility']['predicted_column_a0_cubed'])
        require(np.max(abs(a-b))<=arithmetic(3*m['rows'][alias['fixture']]['atoms'],a,b),'Baseline four-state susceptibility replay')
    profiles={metal:{} for metal in m['metals']}
    for metal in m['metals']:
        for point in m['grid_per_metal']:
            cells={};failed=[]
            for cid,row in m['rows'].items():
                if row['metal']!=metal:continue
                x=matrices[cid]
                v=x['zero'] if point['a']==point['l']==0 else kernel(x['A'],x['ed'],x['ep'],x['alpha'],row['geometry'],point['a'],point['l'],m['units'])
                if v['status']!='complete':failed.append(cid);cells[cid]=v;continue
                delta=v['polarization_energy']-x['zero']['polarization_energy']
                cells[cid]={'status':'complete','energy_kcal_mol':base[cid]['energy_kcal_mol']+delta,'energy_change_kcal_mol':delta,
                            'mu_Debye':(np.asarray(row['permanent_lab_Debye'])+v['induced_conjugate_Debye']).tolist(),'mu_change_Debye':(np.asarray(v['induced_conjugate_Debye'])-x['zero']['induced_conjugate_Debye']).tolist(),
                            'chi_column_a0_cubed':v['chi_column_a0_cubed'],'spectrum':v['spectrum'],'sigma':v['sigma'],'solve_residual_max':v['solve_residual_max']}
            require(len(cells)==27,'Metal-local denominator')
            profiles[metal][point['id']]={'parameters':point,'status':'unstable' if failed else 'complete','failed_fixtures':failed,'cells':cells}
    paired={}
    for la,dy in itertools.product(m['grid_per_metal'],repeat=2):
        a=profiles['La'][la['id']];b=profiles['Dy'][dy['id']];cid=la['id']+'__'+dy['id']
        if a['status']!='complete' or b['status']!='complete':paired[cid]={'status':'unavailable_unstable','La':la,'Dy':dy};continue
        values=dict(m['static']);values.update(a['cells']);values.update(b['cells'])
        paired[cid]={'status':'complete','La':la,'Dy':dy,'observables':observables(values,target,m['field_aliases'])}
    require(len(paired)==81 and all(len(v)==9 for v in profiles.values()),'Finite grid denominator')
    result={'status':'finite_matrix_screen_complete','implementation':pin(__file__),'configuration':cfg,'manifest':cfg['manifest'],'full_result':cfg['full_result'],'full_review':cfg['full_review'],
            'baseline_replay':{'passed':True,'fixtures':replays,'observables':zero,'field_alias_qualification':fields},'profiles':profiles,'paired_combinations':paired,'counts':m['counts'],
            'new_native_calls':0,'fit_performed':False,'aggregate_loss':None,'selected_point':None,'force_checks_evaluated':0,
            'limitations':['Nine tested points per metal do not exhaust the continuous family/domain.','Known analytic Dy same-sign stiffness obstruction remains; no exact repair claimed.','Fixed nonpolarization terms may contain prior-fit compensation; energy tradeoffs do not exclude a future joint energy/response refit.','No coordinate derivatives, force acceptance, reference promotion or production qualification.'],
            'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1},'error':None}
    save(out/'RESULT.json',result)
    lines=['# Finite saved-matrix screen','',f"Retained 9 points per metal, 81 paired combinations; L={m['L']:.12g}, geometry-only gmax={m['geometry_gmax']:.12g}.",'','No loss, ranking, selected candidate or candidate force evaluation.','', '| La point | Dy point | Original energy checks passed /10 | W9 coordination error kcal/mol |','|---|---|---:|---:|']
    for cid,v in paired.items():
        if v['status']!='complete':lines.append(f"| {v['La']['id']} | {v['Dy']['id']} | unavailable: unstable | — |");continue
        e=v['observables']['original_energy_checks'];w=next(x for x in e if x['id']=='exchange__aquo9_central')
        lines.append(f"| {v['La']['id']} | {v['Dy']['id']} | {sum(x['passed'] for x in e)} | {w['error_kcal_mol']:.8f} |")
    lines+=['','RESULT retains every fixture energy/moment change,34 pair cycles,4 C/P/B sets,8 E3/mu3 cycles and28 susceptibility rows for each available combination.','Only tested discrete points are assessed. Existing numerical/reference limitations remain; all15 candidate force checks are unevaluated.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');return result

def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('prepare','run'));p.add_argument('--config',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'One local numerical thread');a.out.mkdir(parents=True,exist_ok=False)
    try:r=prepare(json.loads(a.config.read_text()),a.out) if a.mode=='prepare' else run(json.loads(a.config.read_text()),a.out)
    except Exception as e:save(a.out/'FAILURE.json',{'status':'failed','error':repr(e),'new_native_calls':0});raise
    print(json.dumps({'status':r['status'],'output':str(a.out.resolve())}))

if __name__=='__main__':main()
