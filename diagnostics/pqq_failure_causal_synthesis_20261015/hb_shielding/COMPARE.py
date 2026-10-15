#!/usr/bin/env python3
"""Saved-only native component sensitivities; numerical execution requires Slurm."""
import argparse,copy,hashlib,importlib.util,json,math,os,shutil,sys,time,traceback
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/hb_shielding'
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015/hb_shielding_v1'
ORDER=('origin','adaptive_Ca','adaptive_La');METALS=('Ca','La');POLICIES=('mathematical','operational')
TERMS={'Gelec':('Gelec',),'Ghb':('Ghb',),'Gelec+Ghb':('Gelec','Ghb')};TOL=1e-9;CEILING=1e-6
def require(v,why):
    if not v:raise ValueError(why)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def check(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(p):
    s=importlib.util.spec_from_file_location('_saved_HB_chooser',check(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def bind():
    declaration=read(WORK/'PLAN_DECLARATION.json');require(declaration['plan']==pin(HERE/'PLAN.md'),'Declared plan changed')
    summary=pin(ROOT/'workspaces/native_cavity_accommodation_20261014/summary_v1/execution/RESULT.json');prior=read(summary)
    old,new=[read(prior['inputs'][k]) for k in ('origin_collection','new_collection')]
    ip=pin(ROOT/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json');inv=read(ip)
    chooser=read(ROOT/'workspaces/exaa_metal_position_20261014/full_pool_check_v1/CONFIG.json')['scorer'];rows=[]
    for sample in (0,3):
        src,=[r for r in inv['rows'] if r['case_id'].endswith('sample-'+str(sample))]
        for q in ORDER:
            for z in METALS:
                key=f'sample{sample}_{z}' if q=='origin' else f'sample{sample}_{q}_{z}'
                rows.append({'id':key,'sample':sample,'candidate':q,'metal':z,'source_case_id':src['case_id'],
                    'state_signature':src['state_signature'],'maps':src['maps'],'source_preparation':src['source_preparation'],
                    'cell':src['matrix'][z][q]})
    for n in ('COMPARE.py','RUN.sbatch','PLAN.md'):shutil.copy2(HERE/n,WORK/n)
    c={'protocol':'saved_native_HB_shielding_common3_v1','implementation':pin(WORK/'COMPARE.py'),'wrapper':pin(WORK/'RUN.sbatch'),
       'plan':pin(WORK/'PLAN.md'),'declaration':pin(WORK/'PLAN_DECLARATION.json'),'inventory':ip,'rows':rows,'chooser':chooser,
       'origin_collection':prior['inputs']['origin_collection'],'nonorigin_collection':prior['inputs']['new_collection'],
       'origin_preparation':old['preparation'],'nonorigin_preparation':new['preparation'],'baseline':prior['inputs']['baseline'],
       'prior_Gelec_summary':summary,'origin_execution':pin(Path(prior['inputs']['origin_collection']['path']).parent/'EXECUTION.json'),
       'nonorigin_execution':pin(Path(prior['inputs']['new_collection']['path']).parent/'EXECUTION.json'),
       'term_sets':TERMS,'order':ORDER,'numerical_allowance_Eh':TOL,'informative_ceiling_Eh':CEILING,'new_molecular_calls':0,
       'python':pin(sys.executable),'output':str(WORK/'execution')}
    write(WORK/'CONFIG.json',c);print(json.dumps(pin(WORK/'CONFIG.json')))
def component(v,key,terms):
    vals=[v[key][t] for t in terms];require(all(math.isfinite(x) for x in vals),'Nonfinite native component');return math.fsum(vals)
def fact(value,bound):return {'value_Eh':value,'charge_bound_Eh':bound,'precision_informative':0<=bound<CEILING}
def available_update(u):
    require(u['status']=='available' and u['checks'] and all(x['passed'] for x in u['checks']),'Saved native update not admitted')
    actual=read(u['result']);require(actual=={k:v for k,v in u.items() if k!='result'},'Native collection/result mismatch')
    for k in ('xyz','charges','masks','mapping'):check(u['input'][k])
    read(u['execution']);check(u['cavity'])
def controls(origin):
    scout=next((r for r in origin['rows'] if r['id']=='sample0_Ca'),None);result={}
    for label,terms in TERMS.items():
        checks=[]
        try:
            require(scout is not None,'Scout row unavailable');u=scout['updates']['union'];available_update(u)
            for name in ('translation','rotation','permutation'):
                v=scout['updates'][name];available_update(v)
                d=component(v,'components_Eh',terms)-component(u,'components_Eh',terms)
                checks.append({'control':name,'residual_Eh':d,'allowance_Eh':TOL*len(terms),'passed':abs(d)<=TOL*len(terms)})
            result[label]={'status':'passed' if all(x['passed'] for x in checks) else 'failed','checks':checks,
                'scope':'sample0 Ca origin shared implementation controls; HB rigid pass was not asserted by the original run'}
        except Exception as e:result[label]={'status':'unavailable','reason':str(e),'checks':checks}
    return result
def endpoint(row,record,prepared):
    cell=row['cell'];require(record['status']=='available' and prepared['status']=='prepared','Saved endpoint unavailable')
    require(prepared['source_case_id']==row['source_case_id'] and (prepared['sample'],prepared['metal'])==(row['sample'],row['metal']),'Source join differs')
    for k in ('state_signature','source_preparation'):require(prepared[k]==row[k],'Source state/preparation changed')
    require(record['archived_local']==prepared['archived_native_result']==cell['native_frozen_ALPB_result'],'Local native lineage differs')
    require(prepared['core_atom_count']==cell['atom_count']==192 and prepared['cap_count']==14,'Physical context inventory differs')
    require(prepared['core_charge']==cell['charge'] and prepared['multiplicity']==cell['multiplicity']==1,'Endpoint charge/spin differs')
    require(prepared['variants']['local']['xyz']['sha256']==cell['xyz']['sha256'],'Local coordinate pin differs')
    require(prepared['saved_solvent']['component_v2_acceptance']['accepted'],'Saved component-v2 qualification failed')
    archived=read(record['archived_local']);require(archived['status']=='available' and archived['G_qvac_hartree']==cell['native_frozen_ALPB_hartree'],'Local native value unavailable')
    require(archived['input']['xyz']['sha256']==cell['xyz']['sha256'] and archived['input']['state_signature']==row['state_signature'],'Native input differs')
    mace=read(cell['MACE']);request=read(cell['MACE_request'])
    require(mace['energy_eV']==cell['MACE_eV'] and mace['request']==cell['MACE_request'] and mace['model']==cell['MACE_model'],'Saved MACE receipt changed')
    require(request['xyz']==cell['xyz'] and (request['charge'],request['multiplicity'])==(cell['charge'],cell['multiplicity']),'MACE coordinate/state changed')
    check(cell['xyz']);check(row['state_signature'])
    for p in row['maps'].values():check(p)
    u,removed=[record['updates'][k] for k in ('union','cap_removed')]
    for name,v in (('union',u),('cap_removed',removed)):
        available_update(v);require(v['input']==prepared['variants'][name],'Native variant differs from preparation')
    local=archived['vacuum'];q=[float(x) for x in check(local['charges']).read_text().split()]
    q0=[float(x) for x in check(prepared['variants']['local']['charges']).read_text().split()];qu=[float(x) for x in check(u['input']['charges']).read_text().split()]
    require(len(q)==192 and q==q0==qu[:192] and all(x==0 for x in qu[192:]),'Fixed core/exterior charges differ')
    r={**{k:row[k] for k in ('id','sample','candidate','metal','source_case_id','state_signature','maps','source_preparation')},
       'status':'available','local_native':record['archived_local'],'MACE':cell['MACE'],'union_native':u['result'],'cap_removed_native':removed['result'],
       'baseline':{'MACE_eV':cell['MACE_eV'],'solvent_hartree':cell['native_frozen_ALPB_hartree']},'terms':{},
       'local_cap_zero_status':'available' if 'local' in record['updates'] else 'unavailable_not_saved',
       'cap_zero_repeat_status':'unavailable_not_separately_repeated','source_charge_precision':archived['vacuum_precision']}
    local_update=record['updates'].get('local')
    if local_update:available_update(local_update)
    for label,terms in TERMS.items():
        lv=component(local,'components_Eh',terms);uv=component(u,'components_Eh',terms)
        lb=component(local,'rounding_bounds_Eh',terms);ub=component(u,'rounding_bounds_Eh',terms)
        repeat=component(u,'repeat_components_Eh',terms)-uv
        t={'local_Eh':lv,'union_Eh':uv,'transfer':fact(uv-lv,ub+lb),'repeat':{'residual_Eh':repeat,'passed':abs(repeat)<=TOL*len(terms),'allowance_Eh':TOL*len(terms)},
           'union_cap_zero_Eh':component(u,'cap_zero_components_Eh',terms),'local_cap_zero_Eh':component(local_update,'cap_zero_components_Eh',terms) if local_update else None,
           'cap_volume':fact(component(u,'cap_zero_components_Eh',terms)-component(removed,'components_Eh',terms),component(u,'cap_zero_rounding_bounds_Eh',terms)+component(removed,'rounding_bounds_Eh',terms)),
           'charged_cap_part_in_union':fact(uv-component(u,'cap_zero_components_Eh',terms),ub+component(u,'cap_zero_rounding_bounds_Eh',terms)),
           'cap_removed_repeat':{'residual_Eh':component(removed,'repeat_components_Eh',terms)-component(removed,'components_Eh',terms),'passed':abs(component(removed,'repeat_components_Eh',terms)-component(removed,'components_Eh',terms))<=TOL*len(terms)},
           'baseline_total_charge_bound_Eh':component(local,'rounding_bounds_Eh',('Gelec','Ghb'))}
        if label=='Gelec':require(abs(t['transfer']['value_Eh']-record['transfer_Eh']['Gelec'])<=TOL,'Prior Gelec endpoint replay differs')
        r['terms'][label]=t
    return r
def pool_summary(selector,records,rigid,ha):
    by={(r['metal'],r['candidate']):r for r in records};base={z:{q:by[z,q].get('baseline') for q in ORDER} for z in METALS}
    baseline=selector.choose(base);panels={}
    for label in TERMS:
        matrix={z:{q:({**base[z][q],'solvent_hartree':base[z][q]['solvent_hartree']+by[z,q]['terms'][label]['transfer']['value_Eh']} if base[z][q] else None) for q in ORDER} for z in METALS}
        p=selector.choose(matrix);panel={'status':p['status'],'raw_policies':{},'admitted_policies':{},'admission_reasons':{},'same_candidate_paths':{},'rigid_control':rigid[label]}
        for q in ORDER:
            if all(base[z][q] and base[z]['origin'] for z in METALS):
                v={z:by[z,q]['terms'][label]['transfer']['value_Eh']-by[z,'origin']['terms'][label]['transfer']['value_Eh'] for z in METALS}
                bound=math.fsum(by[z,k]['terms'][label]['transfer']['charge_bound_Eh'] for z in METALS for k in (q,'origin'))
                panel['same_candidate_paths'][q]={'per_metal_change_kcal_mol':{z:v[z]*ha for z in METALS},'paired_change_kcal_mol':(v['Ca']-v['La'])*ha,'charge_bound_Eh':bound,'precision_informative':bound<CEILING}
            else:panel['same_candidate_paths'][q]=None
        for policy in POLICIES:
            value=None
            if p[policy] and baseline[policy]:
                metals={};bound=0.;selection_robust=True
                for z in METALS:
                    old,new=baseline[policy][z],p[policy][z];origin=by[z,'origin']['terms'][label]['transfer']['value_Eh']
                    oldb=max(by[z,q]['terms'][label]['baseline_total_charge_bound_Eh'] for q in ORDER)
                    newb=max(by[z,q]['terms'][label]['baseline_total_charge_bound_Eh']+by[z,q]['terms'][label]['transfer']['charge_bound_Eh'] for q in ORDER)
                    wb=oldb+newb+by[z,'origin']['terms'][label]['transfer']['charge_bound_Eh'];bound+=wb
                    for chosen,err in ((old,oldb),(new,newb)):
                        works=chosen['all_works_kcal_mol'];ordered=sorted(ORDER,key=works.get)
                        if works[ordered[1]]-works[ordered[0]]<=2*err*ha:selection_robust=False
                        if policy=='operational' and abs(works[ordered[0]]+.1)<=2*err*ha:selection_robust=False
                    dw=new['work_kcal_mol']-old['work_kcal_mol']
                    metals[z]={'old_candidate':old['candidate'],'candidate':new['candidate'],'old_work_kcal_mol':old['work_kcal_mol'],'work_kcal_mol':new['work_kcal_mol'],
                        'all_works_kcal_mol':new['all_works_kcal_mol'],'work_change_kcal_mol':dw,'origin_transfer_kcal_mol':origin*ha,'selected_energy_change_kcal_mol':dw+origin*ha,'work_change_charge_bound_Eh':wb}
                value={'metals':metals,'paired_work_change_kcal_mol':metals['Ca']['work_change_kcal_mol']-metals['La']['work_change_kcal_mol'],
                    'paired_pool_change_kcal_mol':metals['Ca']['selected_energy_change_kcal_mol']-metals['La']['selected_energy_change_kcal_mol'],
                    'work_change_charge_bound_Eh':bound,'precision_informative':bound<CEILING,'selection_robust_to_recorded_bounds':selection_robust}
            panel['raw_policies'][policy]=value
            admitted=value is not None and rigid[label]['status']=='passed' and value['precision_informative'] and value['selection_robust_to_recorded_bounds'] and all(r['terms'][label]['repeat']['passed'] for r in records)
            panel['admitted_policies'][policy]=value if admitted else None
            reasons=[]
            if value is None:reasons.append('Required complete pool unavailable')
            if rigid[label]['status']!='passed':reasons.append('Rigid controls '+rigid[label]['status'])
            if value and not value['precision_informative']:reasons.append('Aggregate charge bound exceeds unchanged informative ceiling')
            if value and not value['selection_robust_to_recorded_bounds']:reasons.append('Discrete selection ambiguous within recorded charge bounds')
            if any(r.get('terms',{}).get(label,{}).get('repeat',{}).get('passed') is False for r in records):reasons.append('Native repeat failed')
            panel['admission_reasons'][policy]=reasons
        panels[label]=panel
    return baseline,panels
def run(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'All numerical analysis requires Slurm')
    cp={'path':str(Path(path).resolve()),'sha256':sha};c=read(cp);require(pin(__file__)==c['implementation'],'Wrong comparator');out=Path(c['output']);out.mkdir(exist_ok=False)
    began=time.monotonic();cpu=time.process_time();receipt={'status':'failed','config':cp,'job_id':os.environ['SLURM_JOB_ID'],'new_molecular_calls':0};rows=[]
    try:
        check(c['plan']);check(c['wrapper']);read(c['declaration']);require(c['term_sets']=={k:list(v) for k,v in TERMS.items()} and c['order']==list(ORDER),'Frozen policy differs')
        inventory=read(c['inventory'])
        for row in c['rows']:
            source,=[r for r in inventory['rows'] if r['case_id']==row['source_case_id']]
            require(row['cell']==source['matrix'][row['metal']][row['candidate']] and all(row[k]==source[k] for k in ('state_signature','maps','source_preparation')),'Frozen source/cell inventory differs')
        selector=load(c['chooser']);require(selector.ORDER==ORDER,'Historical order differs');ha=selector.HA
        old,new=[read(c[k]) for k in ('origin_collection','nonorigin_collection')];require(old['status']==new['status']=='complete' and new['origin_collection']==c['origin_collection'],'Saved native matrix incomplete')
        require(all(read(c[k])['status']=='complete' for k in ('origin_execution','nonorigin_execution')),'Saved execution incomplete')
        ops,nps=[{r['id']:r for r in read(c[k])['rows']} for k in ('origin_preparation','nonorigin_preparation')]
        obs,nbs=[{r['id']:r for r in v['rows']} for v in (old,new)];rigid=controls(old)
        for row in c['rows']:
            try:rows.append(endpoint(row,(obs if row['candidate']=='origin' else nbs)[row['id']],(ops if row['candidate']=='origin' else nps)[row['id']]))
            except Exception as e:rows.append({**{k:row[k] for k in ('id','sample','candidate','metal','source_case_id')},'status':'unavailable','reason':str(e),'traceback':traceback.format_exc()})
        require(len(rows)==12 and len({r['id'] for r in rows})==12,'Endpoint denominator changed')
        baseline=read(c['baseline']);prior=read(c['prior_Gelec_summary']);sources=[]
        for sample in (0,3):
            rs=[r for r in rows if r['sample']==sample];before,panels=pool_summary(selector,rs,rigid,ha);saved=next(r for r in baseline['rows'] if r['sample']==sample)
            if before['status']=='available':
                for pol in POLICIES:
                    for z in METALS:
                        a,b=before[pol][z],saved['pools']['historical3'][pol][z];require(a['candidate']==b['candidate'] and abs(a['work_kcal_mol']-b['work_kcal_mol'])<=1e-8,'Old complete-pool replay differs')
                oldp=next(r for r in prior['rows'] if r['sample']==sample);actual=panels['Gelec']['raw_policies']['mathematical']
                require(abs(actual['paired_pool_change_kcal_mol']-oldp['paired_pool_change_kcal_mol'])<=1e-8 and abs(actual['paired_work_change_kcal_mol']-oldp['paired_accommodation_change_kcal_mol'])<=1e-8,'Previous Gelec-only summary not reproduced')
            sources.append({'sample':sample,'source_case_id':saved['source_case_id'],'complete_twelve_matrix_half':before['status']=='available','panels':panels})
        contrasts={}
        for label in TERMS:
            contrasts[label]={}
            for pol in POLICIES:
                a,b=[s['panels'][label]['raw_policies'][pol] for s in sources]
                v=None
                if a and b:
                    bound=a['work_change_charge_bound_Eh']+b['work_change_charge_bound_Eh']
                    v={'pool_change_kcal_mol':b['paired_pool_change_kcal_mol']-a['paired_pool_change_kcal_mol'],'accommodation_change_kcal_mol':b['paired_work_change_kcal_mol']-a['paired_work_change_kcal_mol'],
                        'work_change_charge_bound_Eh':bound,'precision_informative':bound<CEILING,'qualified':bound<CEILING and all(s['panels'][label]['admitted_policies'][pol] is not None for s in sources)}
                contrasts[label][pol]=v
        result={'status':'complete' if all(r['status']=='available' for r in rows) else 'partial','config':cp,'rows':rows,'sources':sources,'source3_minus_source0':contrasts,
            'rigid_controls':rigid,'classification':None,'new_molecular_calls':0,'interpretation':'Raw component sensitivities; failed precision/rigid/selection gates leave admitted comparisons unavailable. Missing local cap-zero controls remain explicit.'}
        write(out/'RESULT.json',result);receipt.update(status=result['status'],result=pin(out/'RESULT.json'))
    except BaseException as e:
        receipt.update(reason=str(e),traceback=traceback.format_exc());done={r['id']:r for r in rows}
        rows=[done.get(r['id'],{**{k:r[k] for k in ('id','sample','candidate','metal','source_case_id')},'status':'unavailable','reason':'Comparator stopped before endpoint admission: '+str(e)}) for r in c['rows']]
        write(out/'PARTIAL_ENDPOINTS.json',rows)
    finally:receipt.update(wall_seconds=time.monotonic()-began,process_CPU_seconds=time.process_time()-cpu);write(out/'EXECUTION.json',receipt)
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='op',required=True);s.add_parser('bind');r=s.add_parser('run');r.add_argument('--config',required=True);r.add_argument('--config-sha',required=True);a=p.parse_args()
    if a.op=='bind':bind()
    else:sys.exit(run(a.config,a.config_sha))
