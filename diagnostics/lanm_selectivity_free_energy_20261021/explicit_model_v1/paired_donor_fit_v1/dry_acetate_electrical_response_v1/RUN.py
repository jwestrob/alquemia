#!/usr/bin/env python3
"""Forty fixed original-model field energies on ten consumed dry-acetate geometries."""
import argparse
import copy
from decimal import Decimal,localcontext
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
PROTOCOL='dry_acetate_electrical_response_v1'
FIELDS={'minus_h':-.001,'minus_half':-.0005,'plus_half':.0005,'plus_h':.001}
QSTEPS=(.01,.005)

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed evidence: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def module(p,name):
    sp=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def save(path,x):
    path=Path(path);tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)
def dot(x,y):return math.fsum(a*b for a,b in zip(x,y))
def dec(s):return Decimal(s.replace('D','E').replace('d','e'))
def half(s):return Decimal(5).scaleb(dec(s).as_tuple().exponent-1)

def render(key,prm,field,direction,units):
    require(field in FIELDS.values() and math.isfinite(field) and len(direction)==3 and all(math.isfinite(v) for v in direction) and abs(dot(direction,direction)-1)<1e-12,'Undeclared field')
    forbidden={'external-field','exfld-freq','polar-local-response','polar-pair-response','exchpol','exchange-polar','a-axis','b-axis','c-axis','ewald','solvate','solvation','active','inactive','charge-flux','chgflux','save-uinduce','save-induced','polar-print'}
    for text in (prm,key):
        for line in text.splitlines():
            a=line.split('#',1)[0].split()
            if a:require(a[0].lower() not in forbidden,'Preexisting field/modified route or reporting keyword')
    require(re.findall(r'(?im)^\s*parameters\s+(\S+)',key)==['amoeba09.prm'],'Original parameter basename differs')
    vector=[field*units.FIELD_FACTOR*x for x in direction]
    return key+'\n# Fixed dry-acetate field-energy discriminator; original donor17621\nEXTERNAL-FIELD '+' '.join(format(v,'.17g') for v in vector)+'\n'

def charge(prm,key,typed):
    types=[int(r[5]) for r in typed];records={}
    for text in (prm,key):
        in_this={}
        for line in text.splitlines():
            a=line.split('#',1)[0].split()
            if a and a[0].lower()=='multipole' and int(a[1]) in types:
                require(int(a[1]) not in in_this,'Ambiguous represented multipole charge');in_this[int(a[1])]=dec(a[-1])
        records.update(in_this)
    require(set(types)<=set(records),'Represented charge absent');total=sum((records[t] for t in types),Decimal(0));require(total==2,'Common model charge must equal+2')
    return {'charge_e':str(total),'site_monopoles_e':[str(records[t]) for t in types],'origin_A':[0.,0.,0.]}

def quantum_replay(row,helper):
    q=read(row['quantum_result']);verify(row['quantum_output']);require(q['output']==row['quantum_output'] and q['status']=='complete','Saved quantum authority differs')
    prop=verify(row['property_file']).read_text();geo=helper.section(prop,'Geometry');require('&Units "Bohr"' in geo,'Quantum coordinate units differ')
    values=re.findall(r'^\s*([A-Za-z]+)\s+('+helper.N+r')\s+('+helper.N+r')\s+('+helper.N+r')\s*$',geo,re.M)
    require(len(values)==8,'Quantum atom inventory differs');xyz=[[float(v) for v in r[1:]] for r in values];z=[11,6,6,8,8,1,1,1]
    dip=helper.section(prop,'SCF_Dipole_Moment');total=helper.array(dip,'dipoleTotal',3);nuc=helper.array(dip,'dipoleNucContrib',3);elec=helper.array(dip,'dipoleElecContrib',3)
    require(max(abs(total[k]-nuc[k]-elec[k]) for k in range(3))<1e-11,'Quantum moment components differ')
    origin=[(math.fsum(z[i]*xyz[i][k] for i in range(8))-nuc[k])/42 for k in range(3)];lab=[total[k]+2*origin[k] for k in range(3)]
    require(row['charge']==2 and row['effective_nuclear_charge']==42 and row['electrons']==40 and max(abs(origin[k]-row['precise_origin_Bohr'][k]) for k in range(3))<1e-12 and max(abs(lab[k]-row['total_dipole_lab_au'][k]) for k in range(3))<1e-12,'Quantum lab-origin replay differs')
    return {'total_lab_au':lab,'origin_Bohr':origin,'charge_e':2,'property_file':row['property_file'],'roundoff_bound_au_per_component':row['printed_origin_total_bound_au'],'origin_rule':'mu(lab0)=mu(native)+2*origin; nuclear42/electronic40'}

def context():
    cfg=read(pin(HERE/'INPUTS.json'));require(cfg['protocol']==PROTOCOL and cfg['specification']==pin(HERE/'SPEC.md'),'Scope/source differs')
    for value in cfg.values():
        if isinstance(value,dict) and set(value)=={'path','sha256'}:verify(value)
    for value in cfg['native_source'].values():verify(value)
    comp=module(cfg['component_adapter'],'field_component_helper');m,energy,native=comp.load_manifest(verify(cfg['component_manifest']));old=read(cfg['component_result']);review=read(cfg['component_actual_review'])
    require(review['status']=='passed_dry_acetate_component_actual_review' and review['result']==cfg['component_result'],'Admitted model component result required')
    require(old['implementation']==cfg['component_adapter'] and old['manifest']==cfg['component_manifest'] and old['status']=='dry_acetate_component_derivatives_complete' and old['complete_cells']==old['new_native_calls']==16 and old['error'] is None,'Complete model zeros unavailable')
    inv=read(cfg['shared_inventory']);require(inv['status']=='saved_response_evidence_and_design_inventory' and inv['authorities']['model_manifest']==cfg['component_manifest'] and inv['authorities']['model_derivative_result']==cfg['component_result'] and inv['authorities']['original_analyze']==m['analyze'],'Shared inventory/model lineage differs')
    units=module(cfg['field_helper'],'qualified_original_field_units');fr=read(cfg['field_actual_review']);require(fr['status']=='passed_independent_actual_model_susceptibility_review' and fr['implementation']==cfg['field_helper'],'Prior field qualification differs')
    prior=read(fr['input_result']);require(prior['controls_passed'] and read(prior['manifest'])['sources']['analyze']==m['analyze'],'Same original binary must have passed field sign/scale/zero controls')
    require(units.FIELDS==FIELDS and units.BOHR==.529177210544 and units.HARTREE==627.509474063 and units.ELECTRIC==332.0637133 and units.ELEFIELD==1439.96455,'Qualified field constants differ')
    qhelper=module(inv['authorities']['property_origin_helper'],'field_quantum_property');direction=inv['geometry_defined_field_direction'];require(len(inv['rows'])==10 and abs(dot(direction,direction)-1)<1e-12,'Ten fixed geometries/direction required')
    ref=module(read(m['inputs'])['reference_adapter'],'field_fixed_geometry');g=m['central']['La']['geometry'];frame=ref.frame(g);t=ref.b.cross(frame['rotation_unit'],ref.b.sub(g['positions_A'][3],frame['centroid_A']));expected=ref.b.unit(t);require(max(abs(x-y) for x,y in zip(direction,expected))<1e-12,'Geometry-only fixed tangent differs')
    refpanel=read(m['reference_panel']);qcollection=read(inv['authorities']['reference_derivative_result']);rows={}
    for gid,shared in inv['rows'].items():
        if gid in ('La_central','Dy_central'):
            metal=gid[:2];c=m['central'][metal];inp=c['source']['effective_input'];authority=c['source']['authority'];parsed=c['energy'];q=0.;geometry=c['geometry'];require(shared['quantum_result']==refpanel['central'][metal]['saved_result'],'Center quantum mapping differs');expected_q_geometry=geometry['xyz']
        else:
            require(gid in m['rows'] and m['rows'][gid]['direction']=='rotation','Only eight existing rotation endpoints allowed');c=m['rows'][gid];metal=c['metal'];q=c['q_A'];inp=c['input'];rec=old['rows'][gid];authority=pin(Path(rec['receipts'][0]['output']['path']).parent/'RESULT.json');require(read(authority)==rec and rec['status']=='complete','Displacement zero authority differs');parsed=comp.parse_energy(verify(rec['receipts'][0]['output']),energy);require(parsed==rec['parsed'],'Raw zero energy differs');geometry=next(x for x in refpanel['cells'] if x['cell_id']==gid);expected_q_geometry=geometry['source_geometry'];require(shared['quantum_result']==qcollection['slots'][gid]['result'],'Rotation quantum mapping differs')
        require(shared['model_input']==inp and shared['model_zero_energy_authority']==authority,'Shared exact model zeros/input differ')
        typed=energy.typed_rows(verify(inp['xyz']).read_text());require(len(typed)==8 and [[float(x) for x in r[2:5]] for r in typed]==geometry['positions_A'] and typed[0][2:5]==['0.000000000000','0.000000000000','0.000000000000'],'Literal lab coordinates/origin differ')
        qraw=read(shared['quantum_result']);verify(qraw['geometry']);require(qraw['geometry']['sha256']==expected_q_geometry['sha256'],'Quantum/model physical geometry differs');qrec=quantum_replay(shared,qhelper);charges=charge(verify(inp['parameters']).read_text(),verify(inp['key']).read_text(),typed)
        rows[gid]={'geometry_id':gid,'metal':metal,'q_A':q,'input':inp,'zero_authority':authority,'zero_energy':parsed,'charge_origin':charges,'quantum':qrec,'quantum_projected_au':dot(qrec['total_lab_au'],direction),'quantum_projected_print_bound_au':sum(abs(v) for v in direction)*qrec['roundoff_bound_au_per_component'],'expected_field_interactions':{k:v+(8 if k in ('Atomic Multipoles','Polarization') else 0) for k,v in parsed['interactions'].items()}}
    require({(r['metal'],r['q_A']) for r in rows.values()}=={(m,q) for m in ('La','Dy') for q in (0.,-.01,-.005,.005,.01)},'Ten geometry denominator differs')
    return cfg,comp,m,energy,native,units,direction,rows


def prepare(out):
    tick=time.monotonic();cpu=time.process_time();cfg,comp,old,energy,native,units,direction,rows=context();out.mkdir(parents=True,exist_ok=False);jobs={}
    for gid,row in rows.items():
        for label,field in FIELDS.items():
            jid=gid+'__'+label;folder=out/'cells'/jid;folder.mkdir(parents=True);src=row['input']
            shutil.copyfile(verify(src['xyz']),folder/'system.xyz');shutil.copyfile(verify(src['parameters']),folder/'amoeba09.prm');(folder/'system.key').write_text(render(verify(src['key']).read_text(),verify(src['parameters']).read_text(),field,direction,units))
            jobs[jid]={'id':jid,'geometry_id':gid,'field_label':label,'signed_field_au':field,'field_au':[field*v for v in direction],'field_MV_cm':[field*units.FIELD_FACTOR*v for v in direction],'input':{k:pin(folder/f) for k,f in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}}
    m={'protocol':PROTOCOL,'implementation':pin(__file__),'specification':cfg['specification'],'inputs':pin(HERE/'INPUTS.json'),'shared_specification':cfg['shared_specification'],'shared_inventory':cfg['shared_inventory'],'analyze':old['analyze'],'selected_review':old['selected_review'],'selected_fit_result':old['selected_fit_result'],'direction':direction,'rows':rows,'jobs':jobs,'maximum_native_calls':40,'workers':1,'conversion':{k:getattr(units,k) for k in ('BOHR','HARTREE','ELECTRIC','ELEFIELD','DEBYE','FIELD_FACTOR')},'new_native_calls':0}
    save(out/'MANIFEST.json',m);load_manifest(out/'MANIFEST.json');r={'status':'dry_electrical_preparation_passed','implementation':pin(__file__),'manifest':pin(out/'MANIFEST.json'),'cells':40,'zero_energies_reused':10,'quantum_dipoles_replayed':10,'new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'process_CPU_seconds':time.process_time()-cpu}};save(out/'RESULT.json',r);print(json.dumps(r,indent=2))

def load_manifest(path):
    cfg,comp,old,energy,native,units,direction,rows=context();m=json.loads(Path(path).read_text());require(m['protocol']==PROTOCOL and m['implementation']==pin(__file__) and m['inputs']==pin(HERE/'INPUTS.json') and m['specification']==cfg['specification'] and m['shared_specification']==cfg['shared_specification'] and m['shared_inventory']==cfg['shared_inventory'],'Prepared source/evidence differs')
    require(m['rows']==rows and m['direction']==direction and m['analyze']==old['analyze'] and m['selected_review']==old['selected_review'] and m['selected_fit_result']==old['selected_fit_result'] and m['maximum_native_calls']==40 and m['workers']==1,'Prepared model/scope differs')
    require(m['conversion']=={k:getattr(units,k) for k in ('BOHR','HARTREE','ELECTRIC','ELEFIELD','DEBYE','FIELD_FACTOR')} and len(m['jobs'])==40,'Conversion/inventory differs')
    require(set(m['jobs'])=={g+'__'+f for g in rows for f in FIELDS} and list(m['jobs'])[:8]==[g+'__'+f for g in ('La_central','Dy_central') for f in FIELDS],'Exact forty field slots/central-first ordering required')
    for jid,job in m['jobs'].items():
        row=rows[job['geometry_id']];f=FIELDS[job['field_label']];src=row['input'];require(jid==job['geometry_id']+'__'+job['field_label'] and job['signed_field_au']==f and job['field_au']==[f*x for x in direction] and job['field_MV_cm']==[f*units.FIELD_FACTOR*x for x in direction],'Rendered field vector differs')
        for key in ('xyz','parameters'):verify(job['input'][key]);require(job['input'][key]['sha256']==src[key]['sha256'],'Geometry/parameters changed')
        require(verify(job['input']['key']).read_text()==render(verify(src['key']).read_text(),verify(src['parameters']).read_text(),f,direction,units),'Field-only key delta differs')
    return m,comp,energy,native

def summarize(m,records):
    require(set(records)==set(m['jobs']) and len(records)==40 and all(r['status']=='complete' for r in records.values()),'Complete40 required; no partial response')
    hconv=dec(str(m['conversion']['HARTREE']));states={};checks={}
    with localcontext() as ctx:
        ctx.prec=50
        for gid,row in m['rows'].items():
            values={};zero=row['zero_energy']['tokens']['Total'];e0=dec(zero)
            for tag,hf in (('h',.001),('half',.0005)):
                h=dec(str(hf));plus=records[gid+'__plus_'+tag]['parsed']['tokens']['Total'];minus=records[gid+'__minus_'+tag]['parsed']['tokens']['Total'];ep,en=dec(plus),dec(minus)
                mu=-(ep-en)/(2*h*hconv);chi=-(ep-2*e0+en)/(h*h*hconv);mub=(half(plus)+half(minus))/(2*h*hconv);chib=(half(plus)+2*half(zero)+half(minus))/(h*h*hconv)
                values[tag]={'field_au':hf,'mu_E_au':float(mu),'mu_E_decimal':str(mu),'mu_print_bound_au':float(mub),'chi_E_a0_cubed':float(chi),'chi_E_decimal':str(chi),'chi_print_bound_a0_cubed':float(chib),'QM_projected_au':row['quantum_projected_au'],'model_minus_QM_mu_au':float(mu)-row['quantum_projected_au'],'model_minus_QM_print_bound_au':float(mub)+row['quantum_projected_print_bound_au']}
            dm=values['half']['mu_E_au']-values['h']['mu_E_au'];dc=values['half']['chi_E_a0_cubed']-values['h']['chi_E_a0_cubed'];checks[gid+'/mu_field_step']=abs(dm)<=1e-4;checks[gid+'/chi_field_step']=abs(dc)<=.5
            states[gid]={'metal':row['metal'],'q_A':row['q_A'],'steps':values,'mu_half_minus_h_au':dm,'chi_half_minus_h_a0_cubed':dc,'mu_step_print_bound_au':sum(v['mu_print_bound_au'] for v in values.values()),'chi_step_print_bound_a0_cubed':sum(v['chi_print_bound_a0_cubed'] for v in values.values())}
        def state(metal,q):return next(v for v in states.values() if v['metal']==metal and v['q_A']==q)
        paired={}
        for q in (0.,-.01,-.005,.005,.01):
            la,dy=state('La',q),state('Dy',q);steps={}
            for tag in ('h','half'):
                a,b=la['steps'][tag],dy['steps'][tag];steps[tag]={'mu_E_au':b['mu_E_au']-a['mu_E_au'],'mu_print_bound_au':a['mu_print_bound_au']+b['mu_print_bound_au'],'QM_projected_au':b['QM_projected_au']-a['QM_projected_au'],'chi_E_a0_cubed':b['chi_E_a0_cubed']-a['chi_E_a0_cubed'],'chi_print_bound_a0_cubed':a['chi_print_bound_a0_cubed']+b['chi_print_bound_a0_cubed'],'model_minus_QM_print_bound_au':a['model_minus_QM_print_bound_au']+b['model_minus_QM_print_bound_au']};steps[tag]['model_minus_QM_mu_au']=steps[tag]['mu_E_au']-steps[tag]['QM_projected_au']
            delta=steps['half']['chi_E_a0_cubed']-steps['h']['chi_E_a0_cubed'];checks[f'paired/q={q}/chi_field_step']=abs(delta)<=1.;paired[str(q)]={'steps':steps,'chi_half_minus_h_a0_cubed':delta,'mu_half_minus_h_au':steps['half']['mu_E_au']-steps['h']['mu_E_au'],'implied_mu_step_envelope_au':2e-4,'mu_step_print_bound_au':sum(v['mu_print_bound_au'] for v in steps.values()),'chi_step_print_bound_a0_cubed':sum(v['chi_print_bound_a0_cubed'] for v in steps.values())}
        derivatives={}
        for kind in ('La','Dy','paired'):
            coeff={'La':1} if kind=='La' else {'Dy':1} if kind=='Dy' else {'La':-1,'Dy':1};derivatives[kind]={}
            for q in QSTEPS:
                values={}
                for tag in ('h','half'):
                    value=0.;bound=0.;qm=0.;qmb=0.
                    for metal,c in coeff.items():
                        p,n=state(metal,q)['steps'][tag],state(metal,-q)['steps'][tag];value+=c*(p['mu_E_au']-n['mu_E_au'])/(2*q);bound+=abs(c)*(p['mu_print_bound_au']+n['mu_print_bound_au'])/(2*q);qm+=c*(p['QM_projected_au']-n['QM_projected_au'])/(2*q)
                        qmb+=abs(c)*sum(r['quantum_projected_print_bound_au'] for r in m['rows'].values() if r['metal']==metal and abs(r['q_A'])==q)/(2*q)
                    values[tag]={'dmu_dq_ea0_per_A':value,'print_bound_ea0_per_A':bound,'QM_dmu_dq_ea0_per_A':qm,'QM_print_bound_ea0_per_A':qmb,'model_minus_QM_ea0_per_A':value-qm,'model_minus_QM_print_bound_ea0_per_A':bound+qmb}
                delta=values['half']['dmu_dq_ea0_per_A']-values['h']['dmu_dq_ea0_per_A'];envelope=1e-4*sum(abs(c) for c in coeff.values())/q;checks[f'{kind}/q={q}/dmu_field_step']=abs(delta)<=envelope
                derivatives[kind][str(q)]={'field_steps':values,'half_minus_h_ea0_per_A':delta,'inherited_field_step_envelope_ea0_per_A':envelope,'field_step_print_bound_ea0_per_A':sum(v['print_bound_ea0_per_A'] for v in values.values())}
            derivatives[kind]['smaller_minus_larger_q_descriptive']={tag:{k:derivatives[kind]['0.005']['field_steps'][tag][k]-derivatives[kind]['0.01']['field_steps'][tag][k] for k in ('dmu_dq_ea0_per_A','QM_dmu_dq_ea0_per_A','model_minus_QM_ea0_per_A')} for tag in ('h','half')}
    require(len(checks)==31,'Declared numerical denominator differs')
    return {'states':states,'paired_Dy_minus_La':paired,'rotation_dipole_derivatives':derivatives,'numerical_checks':checks,'numerical_check_count':31,'numerical_checks_passed':all(checks.values()),'physical_accuracy_qualified':False,'unique_charge_transfer_or_parameter_attribution':False,'new_quantum_calls':0}

def central_gate(m,records):
    expected=[g+'__'+f for g in ('La_central','Dy_central') for f in FIELDS]
    require(all(g in records and records[g]['status']=='complete' for g in expected),'Both central quartets required')
    checks={};values={}
    with localcontext() as ctx:
        ctx.prec=50;K=dec(str(m['conversion']['HARTREE']))
        for gid in ('La_central','Dy_central'):
            zero=m['rows'][gid]['zero_energy']['tokens']['Total'];values[gid]={}
            for tag,h in (('h',Decimal('.001')),('half',Decimal('.0005'))):
                p=records[gid+'__plus_'+tag]['parsed']['tokens']['Total'];n=records[gid+'__minus_'+tag]['parsed']['tokens']['Total']
                values[gid][tag]={'mu':float(-(dec(p)-dec(n))/(2*h*K)),'chi':float(-(dec(p)+dec(n)-2*dec(zero))/(h*h*K)),'mu_print_bound':float((half(p)+half(n))/(2*h*K)),'chi_print_bound':float((half(p)+half(n)+2*half(zero))/(h*h*K))}
            checks[gid+'/mu']=abs(values[gid]['h']['mu']-values[gid]['half']['mu'])<=1e-4;checks[gid+'/chi']=abs(values[gid]['h']['chi']-values[gid]['half']['chi'])<=.5
        checks['paired_chi']=abs((values['Dy_central']['h']['chi']-values['La_central']['h']['chi'])-(values['Dy_central']['half']['chi']-values['La_central']['half']['chi']))<=1
    return {'passed':all(checks.values()),'checks':checks,'values':values,'native_calls':8}

def release_gate(m,manifest,path,out):
    r=json.loads(Path(path).read_text());require(r['status']=='released_dry_acetate_electrical_response','Root release required')
    require(r['implementation']==pin(__file__) and r['manifest']==pin(manifest) and r['analyze']==m['analyze'] and r['shared_specification']==m['shared_specification'] and r['maximum_native_calls']==40 and r['workers']==1 and r['cell_ids']==list(m['jobs']) and Path(r['output_directory']).resolve()==out.resolve(),'Root release binding/scope differs')
    review=read(r['independent_source_review']);require(review['status']=='passed_dry_acetate_electrical_source_review' and review['implementation']==pin(__file__) and review['manifest']==pin(manifest),'Source review differs');return r

def run(path,release,out):
    m,comp,energy,native=load_manifest(path);release_gate(m,path,release,out);require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'One Slurm CPU required');require(all(os.environ.get(k)=='1' for k in THREADS),'One thread required');out.mkdir(parents=True,exist_ok=False)
    tick=time.monotonic();cpu=time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN);result={'protocol':PROTOCOL,'status':'running','implementation':pin(__file__),'manifest':pin(path),'release':pin(release),'analyze':m['analyze'],'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'rows':{jid:{'status':'unattempted','cell_id':jid} for jid in m['jobs']},'new_native_calls':0,'error':None,'summary':None,'central_gate':None};save(out/'START.json',{k:v for k,v in result.items() if k not in ('rows','summary')});save(out/'RESULT.json',result)
    try:
        for jid,job in m['jobs'].items():
            require(result['new_native_calls']<40,'Forty-call cap');require(result['new_native_calls']<8 or result['central_gate']['passed'],'Central8 numerical admission required');folder=out/'cells'/jid;folder.mkdir(parents=True);rec={'cell_id':jid,'status':'started','receipts':[]};result['rows'][jid]=rec
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(job['input'][k]),folder/name)
            rec['input']={k:pin(folder/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))};save(folder/'START.json',{'job':job,'manifest':pin(path),'release':pin(release),'input':rec['input'],'analyze':m['analyze']});result['new_native_calls']+=1;save(out/'RESULT.json',result)
            try:
                native.invoke([str(verify(m['analyze'])),'system.xyz','E'],folder,folder/'analyze.out',rec['receipts']);rec['parsed']=comp.parse_energy(folder/'analyze.out',energy);require(rec['parsed']['interactions']==m['rows'][job['geometry_id']]['expected_field_interactions'],'Field-enabled component interaction counts differ');rec['status']='complete'
            except Exception as exc:rec.update(status='failed',error=repr(exc));raise
            finally:save(folder/'RESULT.json',rec);save(out/'RESULT.json',result)
            if result['new_native_calls']==8:
                result['central_gate']=central_gate(m,result['rows']);save(out/'RESULT.json',result);require(result['central_gate']['passed'],'Central field-step qualification failed; remaining32 stopped')
        result['summary']=summarize(m,result['rows']);result['status']='dry_acetate_electrical_response_complete' if result['summary']['numerical_checks_passed'] else 'dry_acetate_electrical_response_failed_numerical_checks'
    except Exception as exc:result.update(status='dry_acetate_electrical_response_failed',error=repr(exc))
    finally:
        end=resource.getrusage(resource.RUSAGE_CHILDREN);wall=time.monotonic()-tick;result['complete_cells']=sum(r['status']=='complete' for r in result['rows'].values());result['cost']={'wall_seconds':wall,'process_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':end.ru_utime-child.ru_utime,'child_system_CPU_seconds':end.ru_stime-child.ru_stime,'native_elapsed_seconds':sum(r.get('elapsed_seconds',0) for v in result['rows'].values() for r in v.get('receipts',[])),'allocated_CPUs':1,'allocated_CPU_seconds':wall};save(out/'RESULT.json',result)
    return 0 if result['status']=='dry_acetate_electrical_response_complete' else 1

def admit_result(manifest_path,result_path):
    m,comp,energy,native=load_manifest(manifest_path);r=json.loads(Path(result_path).read_text());release_gate(m,manifest_path,verify(r['release']),Path(result_path).parent)
    require(r['implementation']==pin(__file__) and r['manifest']==pin(manifest_path) and r['analyze']==m['analyze'] and r['error'] is None and r['new_native_calls']==r['complete_cells']==40 and r['status'] in ('dry_acetate_electrical_response_complete','dry_acetate_electrical_response_failed_numerical_checks'),'Complete40 authority required')
    require(set(r['rows'])==set(m['jobs']),'Forty-slot denominator differs')
    for jid,job in m['jobs'].items():
        rec=r['rows'][jid];folder=Path(result_path).parent/'cells'/jid;require(rec==json.loads((folder/'RESULT.json').read_text()) and rec['status']=='complete','Stored cell differs')
        start=json.loads((folder/'START.json').read_text());require(start=={'job':job,'manifest':pin(manifest_path),'release':r['release'],'input':rec['input'],'analyze':m['analyze']},'Native START differs')
        for k in ('xyz','key','parameters'):verify(rec['input'][k]);require(rec['input'][k]['sha256']==job['input'][k]['sha256'],'Actual input differs')
        require(len(rec['receipts'])==1,'One program per cell required');q=rec['receipts'][0];require(q['returncode']==0 and q['argv']==[m['analyze']['path'],'system.xyz','E'] and Path(q['cwd'])==folder and all(v=='1' for v in q['threads'].values()),'Native call differs')
        require(verify(q['output'])==folder/'analyze.out' and comp.parse_energy(folder/'analyze.out',energy)==rec['parsed'],'Actual raw energy replay differs');require(rec['parsed']['interactions']==m['rows'][job['geometry_id']]['expected_field_interactions'],'Field term counts differ')
    require(central_gate(m,r['rows'])==r['central_gate'] and r['central_gate']['passed'],'Saved central gate differs');require(summarize(m,r['rows'])==r['summary'],'Saved response arithmetic differs');return r

if __name__=='__main__':
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='mode',required=True);p=sub.add_parser('prepare');p.add_argument('--out',type=Path,required=True);r=sub.add_parser('run');r.add_argument('--manifest',type=Path,required=True);r.add_argument('--release',type=Path,required=True);r.add_argument('--out',type=Path,required=True);a=sub.add_parser('admit');a.add_argument('--manifest',type=Path,required=True);a.add_argument('--result',type=Path,required=True);args=ap.parse_args()
    if args.mode=='prepare':prepare(args.out)
    elif args.mode=='run':raise SystemExit(run(args.manifest,args.release,args.out))
    else:print(json.dumps({'status':admit_result(args.manifest,args.result)['status'],'new_native_calls':0}))
