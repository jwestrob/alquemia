#!/usr/bin/env python3
"""Ten fixed dry donor SP fields; qualified native/parser/estimator reuse, staged 5+5."""
import argparse, copy, importlib.util, json, math, os, re, shutil, socket, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def module(path,name):
    s=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
old=module(HERE.parent/'conditional_susceptibility_v1/REFERENCE.py','dry_field_old')
dry=module(HERE.parent/'dry_acetate_force_qualification_v1/REFERENCE.py','dry_field_native')
b=dry.b;read=old.read;save=old.save;canonical=old.canonical
prop=old.prop;moment=old.moment;token_roundoff=old.token_roundoff
PROTOCOL='dry_donor_susceptibility_v1';METALS=('La','Dy');THREADS=old.THREADS
WS=b.ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
INVENTORY=WS/'reference_chemistry_v1/dry_donor_susceptibility_v1/SOURCE_INVENTORY.json'
INVENTORY_SHA='cfa23743f3431468faa4e92d13d53c64f518a19d8e233edad4dc927f13a12221'
MODEL_REVIEW=WS/'explicit_model_v1/paired_donor_fit_v1/dry_acetate_electrical_response_v1/independent_actual_review_v1/REVIEW.json'
MODEL_REVIEW_SHA='f21ff157ae99273ef7df326dca64a5647852e79492ed4c6839cd7c0e6875c122'
FIELDS=dict(zero=0.,**old.FIELDS)
HELPERS=dict(old_field=Path(old.__file__),dry_native=Path(dry.__file__),**old.HELPERS,
             dry_amide=Path(dry.h.__file__),dry_water=Path(dry.w.__file__))


def default_config():
    b.check(b.pin(INVENTORY)['sha256']==INVENTORY_SHA,'Frozen design inventory changed')
    inv=read(b.pin(INVENTORY));b.check(inv['specification']==b.pin(HERE/'SPEC.md'),'Design specification changed')
    b.check(b.pin(MODEL_REVIEW)['sha256']==MODEL_REVIEW_SHA,'Model actual review changed')
    return dict(protocol=PROTOCOL,implementation=b.pin(__file__),specification=inv['specification'],source_inventory=b.pin(INVENTORY),
        model_actual_review=b.pin(MODEL_REVIEW),orca=inv['centers']['La']['original_orca'],run_wrapper=b.pin(HERE/'RUN.sbatch'),
        **{k:b.pin(v) for k,v in HELPERS.items()})


def config_check(cfg):
    b.check(cfg==default_config(),'Configuration/source changed');dry.source_gates()
    for v in cfg.values():
        if isinstance(v,dict) and set(v)=={'path','sha256'}:b.verify(v)
    inv=read(cfg['source_inventory'])
    b.check(inv['old_conditional_field_adapter']==cfg['old_field'] and inv['dry_native_parser']==cfg['dry_native'],'Design/helper joins differ')


def field_line(c):return ' EField '+','.join(f'{x:.17g}' for x in c['field_vector_au'])
def field_block(c):return field_line(c)+'\n EFieldOrigin 0.0,0.0,0.0\n'
PROPERTY_BLOCK='%elprop\n Dipole true\n Origin 0.0,0.0,0.0\nend\n'

def base_deck(p,c,ranks,maxcore):
    original=p['centers'][c['metal']]['zero_EnGrad_input']
    text=old.base_deck({'input':original},ranks,maxcore)
    b.check(text.count(' EnGrad')==1 and text.count('* xyzfile 2 1 geometry.xyz')==1,'Original central recipe differs')
    return text.replace(' EnGrad','',1)


def deck(p,c,ranks,maxcore):
    text=base_deck(p,c,ranks,maxcore)
    b.check(text.count('%scf\n')==1 and 'EField' not in text and '%elprop' not in text,'Original field/property recipe differs')
    return text.replace('%scf\n','%scf\n'+field_block(c),1).replace('* xyzfile ',PROPERTY_BLOCK+'* xyzfile ',1)


def cells(inv):
    out=[]
    for row in inv['cells']:
        c=dry.h.metadata(inv['centers'][row['metal']]['geometry']);c.update(copy.deepcopy(row))
        c['expected_effective_nuclear_charges']=inv['centers'][row['metal']]['expected_effective_nuclear_charges']
        c['expected_aux_J_basis_functions']=inv['centers'][row['metal']]['expected_aux_J_functions']
        out.append(c)
    return out


def old_center(center):
    r=read(center['zero_EnGrad_result']);start=read(center['zero_EnGrad_start']);oc=center['original_cell']
    b.check(r['status']=='complete' and r['returncode']==0 and r['input']==center['zero_EnGrad_input'] and r['output']==center['zero_EnGrad_output'],'Old central result differs')
    b.check(start['cell']==oc and start['panel']==center['original_panel'] and start['input']==r['input'],'Old central runtime differs')
    parent=read(center['original_panel']);b.check(parent['orca']==center['original_orca'] and parent['assets'][oc['metal']]==center['assets'],'Old binary/assets differ')
    b.check(b.verify(r['input']).read_text()==b.deck(parent,oc,start['ranks'],start['maxcore_MiB']),'Old central input renderer differs')
    q=dry.parse(b.verify(r['output']).read_text(),dry.h.metadata(center['geometry']),0)
    b.check(q['status']=='complete' and q['energy_hartree']==r['energy_hartree'],'Old native replay differs')
    mu=moment.quantum_moment(r['output'],2,prop)
    b.check(mu['property_file']==center['zero_EnGrad_property'] and mu['total_dipole_au']==center['old_total_dipole_lab_au'],'Old common-origin dipole differs')
    return dict(native=q,moment=mu,result=center['zero_EnGrad_result'],context_only=True,not_SP_zero=True)


def prepare(config,out):
    tick,cpu=time.monotonic(),time.process_time();cfg=json.loads(Path(config).read_text());config_check(cfg);inv=read(cfg['source_inventory'])
    dry.w.recursive_verify(inv);zeros={m:old_center(inv['centers'][m]) for m in METALS}
    b.check(inv['centers']['La']['geometry']['positions_A']==inv['centers']['Dy']['geometry']['positions_A'],'Paired coordinates differ')
    out.mkdir(parents=True,exist_ok=False);(out/'inputs').mkdir()
    p=dict(protocol=PROTOCOL,adapter=cfg['implementation'],configuration=b.pin(config),specification=cfg['specification'],
        source_inventory=cfg['source_inventory'],orca=cfg['orca'],centers=inv['centers'],old_EnGrad_context=zeros,
        field_direction_lab=inv['field_direction_lab'],cells=cells(inv),scout_cell_ids=inv['first_stage_ids'],remaining_cell_ids=inv['remaining_stage_ids'],
        maximum_native_calls=10,ranks=4,nodes=1,run_wrapper=cfg['run_wrapper'],model_actual_review=cfg['model_actual_review'],**{k:cfg[k] for k in HELPERS})
    for c in p['cells']:
        f=out/'inputs'/(c['cell_id']+'.inp');f.write_text(deck(p,c,4,3000));c['template_input']=b.pin(f)
    pp=save(out/'PANEL.json',p);v=validate(pp['path']);v.update(old_central_raw_replays=2,old_SP_zeros_reused=0,
        cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    pf=save(out/'PREFLIGHT.json',v);print(json.dumps(dict(panel=pp,preflight=pf,status=v['status'],new_native_calls=0)))


def load(path):
    p=json.loads(Path(path).read_text());cfg=read(p['configuration']);config_check(cfg);inv=read(cfg['source_inventory'])
    for k in ('specification','source_inventory','orca','run_wrapper','model_actual_review',*HELPERS):b.check(p[k]==cfg[k],'Panel binding differs: '+k)
    b.check(p['protocol']==PROTOCOL and p['adapter']==cfg['implementation'] and p['centers']==inv['centers'],'Panel identity/centers differ')
    b.check(p['maximum_native_calls']==len(p['cells'])==10 and p['ranks']==4 and p['nodes']==1 and
        p['scout_cell_ids']==inv['first_stage_ids'] and p['remaining_cell_ids']==inv['remaining_stage_ids'] and
        p['field_direction_lab']==inv['field_direction_lab'],'Finite inventory differs')
    for a,e in zip(p['cells'],cells(inv)):
        b.check({k:v for k,v in a.items() if k!='template_input'}==e,'Exact cell changed')
        b.check(b.verify(a['template_input']).read_text()==deck(p,a,4,3000),'Exact template changed')
    return p


def validate(path):
    p=load(path);records={};b.check(abs(math.fsum(x*x for x in p['field_direction_lab'])-1)<1e-14,'Field direction not unit')
    for c in p['cells']:
        text=b.verify(c['template_input']).read_text();stripped=text.replace(field_block(c),'',1).replace(PROPERTY_BLOCK,'',1)
        b.check(stripped==base_deck(p,c,4,3000),'More than declared field/property/EnGrad delta')
        symbols,x=b.xyz_read(b.verify(c['source_geometry']))
        b.check(x==c['positions_A'] and symbols==c['symbols'] and x[0]==[0.,0.,0.],'Literal geometry differs')
        b.check(c['charge']==2 and c['multiplicity']==1 and c['expected_electrons']==40 and c['expected_basis_functions']==189 and
            c['expected_aux_J_basis_functions']==(580 if c['metal']=='La' else 598) and not any(c['ghost_mask']) and c['d4'],'Native state/basis/D4 differs')
        b.check(c['field_au']==FIELDS[c['field_tag']] and c['field_vector_au']==[c['field_au']*z for z in p['field_direction_lab']],'Exact field differs')
        records[c['cell_id']]=dict(only_declared_delta=True,charge=2,electrons=40,AO=189,aux_J=c['expected_aux_J_basis_functions'],new_SP_zero=c['field_tag']=='zero')
    return dict(status='passed_static_preflight',panel=b.pin(path),declared=10,scout_cell_ids=p['scout_cell_ids'],records=records,new_native_calls=0)


def parse(text,c,returncode):
    q=dry.parse(text,c,returncode)
    echo=[re.sub(r'^\s*\|\s*\d+> ?','',x) for x in text.split('****END OF INPUT****')[0].splitlines() if re.match(r'^\s*\|\s*\d+>',x)]
    q['checks']['declared_field_echo']=echo.count(field_line(c))==1 and sum(bool(re.search(r'\bEField\b',x)) for x in echo)==1
    q['checks']['declared_origins_echo']=echo.count(' EFieldOrigin 0.0,0.0,0.0')==1 and echo.count(' Origin 0.0,0.0,0.0')==1
    q['status']='complete' if all(q['checks'].values()) else 'failed_admission';return q


def one(pattern,text,label):
    hits=re.findall(pattern,text,re.M);b.check(len(hits)==1,'Missing/duplicate '+label);return hits[0]
def array_tokens(block,name,count):
    parts=block.split('&'+name+' ');b.check(len(parts)==2,'Missing/duplicate property array '+name)
    tail=parts[1].split('\n   &',1)[0];values=re.findall(r'^\s*(\d+)\s+('+moment.NUMBER+r')\s*$',tail,re.M)
    b.check([int(i) for i,_ in values]==list(range(count)),'Property array indices differ: '+name)
    tokens=[t for _,t in values];b.check(all(math.isfinite(float(t)) for t in tokens),'Nonfinite property array');return tokens

def field_properties(output,property_path,c,q):
    text=Path(output).read_text();raw=Path(property_path).read_text();b.check('ORCA 6.1.1' in raw,'Property version differs')
    origin=one(r'^Position of electric origin\s*\.\.\.\s*('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s*$',text,'electric origin')
    b.check(list(map(float,origin))==[0.,0.,0.],'Native electric origin is not lab zero')
    geo=prop.section(raw,'Geometry');b.check('&Units "Bohr"' in geo,'Property positions not Bohr')
    rr=re.findall(r'^\s*([A-Za-z]+)(:?)\s+('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s*$',geo.split('&Units "Bohr"]',1)[1],re.M)
    b.check(len(rr)==8 and [r[0] for r in rr]==c['symbols'] and not any(r[1] for r in rr),'Property eight-real-center inventory differs')
    xyz=[[float(v) for v in r[2:]] for r in rr]
    b.check(max(abs(x*moment.BOHR-y) for row,z in zip(xyz,c['positions_A']) for x,y in zip(row,z))<=1e-6,'Property geometry differs')
    info=prop.section(raw,'Calculation_Info')
    for name,value in dict(Charge=2,Mult=1,NumOfAtoms=8,NumOfElectrons=40,NumOfBasisFuncts=189,NumOfAuxJBasisFuncts=c['expected_aux_J_basis_functions']).items():
        b.check(old.scalar(info,name)==value,'Property state/basis differs: '+name)
    block=prop.section(raw,'SCF_Dipole_Moment');values={};bounds={};tokens={}
    for key,name,label in [('total','dipoleTotal','Total Dipole Moment'),('nuclear','dipoleNucContrib','Nuclear contribution'),('electronic','dipoleElecContrib','Electronic contribution')]:
        ts=array_tokens(block,name,3);values[key]=[float(t) for t in ts];bounds[key]=[token_roundoff(t) for t in ts];tokens[key]=ts
        printed=one(r'^\s*'+re.escape(label)+r'\s*:\s*('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s*$',text,label)
        b.check(all(abs(values[key][j]-float(t))<=bounds[key][j]+token_roundoff(t) for j,t in enumerate(printed)),'Property/native dipole token closure differs')
    z=c['expected_effective_nuclear_charges'];b.check(z==[11,6,6,8,8,1,1,1] and sum(z)-40==2,'Effective nuclear charges differ')
    nuclear=[math.fsum(a*r[j] for a,r in zip(z,xyz)) for j in range(3)]
    nb=[math.fsum(abs(a)*token_roundoff(r[j+2]) for a,r in zip(z,rr))+8*math.ulp(nuclear[j]) for j in range(3)]
    b.check(all(abs(nuclear[j]-values['nuclear'][j])<=nb[j]+bounds['nuclear'][j] for j in range(3)),'Effective-ECP nuclear origin closure differs')
    b.check(all(abs(values['nuclear'][j]+values['electronic'][j]-values['total'][j])<=sum(bounds[k][j] for k in bounds)+4*math.ulp(abs(values['nuclear'][j])+abs(values['electronic'][j])) for j in range(3)),'Dipole component closure differs')
    sp=prop.section(raw,'Single_Point_Data');b.check('&Converged [&Type "Boolean"] true' in sp,'Property convergence absent')
    et=one(r'&FinalEnergy\s+\[&Type "Double"\]\s+('+moment.NUMBER+r')',sp,'FinalEnergy')
    st=array_tokens(prop.section(raw,'SCF_Energy'),'totalEnergy',1)[0]
    nt=one(r'^FINAL SINGLE POINT ENERGY\s+('+moment.NUMBER+r')\s*$',text,'native final energy')
    ns=one(r'^Total Energy\s*:\s*('+moment.NUMBER+r')\s+Eh',text,'native electronic energy')
    dt=one(r'^Dispersion correction\s+('+moment.NUMBER+r')\s*$',text,'native D4')
    e,s,d=float(et),float(st),float(dt);eb,sb,db=map(token_roundoff,(et,st,dt))
    b.check(float(nt)==q['energy_hartree'] and float(ns)==q['electronic_energy_hartree'] and d==q['D4_energy_hartree'],'Parsed energy stages differ')
    b.check(abs(e-float(nt))<=eb+token_roundoff(nt) and abs(s-float(ns))<=sb+token_roundoff(ns),'Property/native energy token closure differs')
    closure=e-s-d;cb=eb+sb+db+4*math.ulp(abs(e)+abs(s))
    b.check(abs(closure)<=cb,'Final = electronic + D4 token closure failed')
    return dict(property_file=b.pin(property_path),dipole_au=values['total'],dipole_components_au=values,positions_Bohr=xyz,origin_Bohr=[0.,0.,0.],
        energy_hartree=s,SCF_energy_hartree=s,final_energy_hartree=e,D4_energy_hartree=d,energy_stage_closure_hartree=closure,energy_stage_closure_bound_hartree=cb,
        raw_tokens=dict(electronic=st,final=et,native_electronic=ns,native_final=nt,D4=dt,dipoles=tokens),
        precision=dict(energy_hartree=sb,final_energy_hartree=eb,D4_energy_hartree=db,moment_components_au=bounds['total'],
            meaning='Raw decimal plus binary bounds only; no SCF or physical accuracy bound'))


def view(props,t,final=False):
    # Slot 2 is the inherited estimator's scalar API; the actual field is along t, not lab z.
    mu=math.fsum(a*x for a,x in zip(t,props['dipole_au']));pb=math.fsum(abs(a)*x for a,x in zip(t,props['precision']['moment_components_au']))+4*math.ulp(sum(abs(a*x) for a,x in zip(t,props['dipole_au'])))
    key='final_energy_hartree' if final else 'energy_hartree'
    return dict(energy_hartree=props[key],dipole_au=[0.,0.,mu],precision=dict(energy_hartree=props['precision'][key],moment_z_au=pb))


def score(rows,p,stage):
    metals=['La'] if stage=='scout' else list(METALS);per={};checks={};d4checks={};context={}
    for m in metals:
        cc=[c for c in p['cells'] if c['metal']==m];needed=[c['cell_id'] for c in cc]
        if not all(rows.get(cid,{}).get('status')=='complete' for cid in needed):
            per[m]=None
            for k in ('energy_vs_moment_h','energy_vs_moment_half','chi_energy_step','chi_moment_step','odd_slope_h','odd_slope_half'):checks[m+'__'+k]=False
            d4checks[m]=False;continue
        by={c['field_tag']:rows[c['cell_id']]['properties'] for c in cc};vv={k:view(v,p['field_direction_lab']) for k,v in by.items()}
        est=old.estimate(vv,vv['zero']);per[m]=old.consistency(est,.5,odd=True);checks.update({m+'__'+k:v for k,v in per[m]['checks'].items()})
        d0=by['zero']['D4_energy_hartree'];d0b=by['zero']['precision']['D4_energy_hartree']
        d4checks[m]=all(abs(v['D4_energy_hartree']-d0)<=v['precision']['D4_energy_hartree']+d0b for v in by.values())
        tv={k:view(v,p['field_direction_lab'],True) for k,v in by.items()};total=old.estimate(tv,tv['zero'])
        per[m]['final_total_estimates']=total;per[m]['total_minus_electronic_curvature']={s:total[s]['chi_energy']-est[s]['chi_energy'] for s in ('h','half')}
        per[m]['D4_values_hartree']={k:v['D4_energy_hartree'] for k,v in by.items()}
        per[m]['D4_field_independence_passed']=d4checks[m]
        oc=p['old_EnGrad_context'][m];context[m]=dict(new_SP_minus_old_EnGrad_total_hartree=by['zero']['final_energy_hartree']-oc['native']['energy_hartree'],
            new_SP_minus_old_EnGrad_dipole_au=[x-y for x,y in zip(by['zero']['dipole_au'],oc['moment']['total_dipole_au'])],descriptive_only=True)
    paired=None
    if stage=='full':
        if all(per[m] is not None for m in METALS):paired=old.consistency(old.combine({m:per[m]['estimates'] for m in METALS},{'Dy':1,'La':-1}),1.)
        for k in ('energy_vs_moment_h','energy_vs_moment_half','chi_energy_step','chi_moment_step'):checks['Dy_minus_La__'+k]=False if paired is None else paired['checks'][k]
    return dict(per_metal=per,Dy_minus_La=paired,checks=checks,declared_checks=6 if stage=='scout' else 16,
        numerical_checks_passed=all(checks.values()),D4_field_independence=d4checks,D4_field_independence_passed=all(d4checks.values()),old_EnGrad_context_differences=context,
        passed=all(checks.values()) and all(d4checks.values()),field_direction_lab=p['field_direction_lab'],physical_accuracy_qualified=False)


def release_contract(p,path,r):
    b.check(r['status']=='approved_dry_donor_susceptibility_reference' and r['implementation']==p['adapter'] and r['panel']==b.pin(path) and r['specification']==p['specification'],'Root source release differs')
    ids=[c['cell_id'] for c in p['cells']]
    b.check(r['cell_ids']==ids and r['maximum_native_calls']==10 and r['ranks']==4 and r['nodes']==1,'Released finite scope differs')
    prep=read(r['preparation']);review=read(r['independent_review'])
    b.check(prep['status']=='passed_static_preflight' and prep['panel']==b.pin(path),'Preparation unavailable')
    b.check(review['status']=='passed_dry_donor_susceptibility_source_review' and review['implementation']==p['adapter'] and review['panel']==b.pin(path),'Independent source/preparation review differs')
    b.check(r['stage'] in ('scout','remaining'),'Unknown release stage')
    expected=p['scout_cell_ids'] if r['stage']=='scout' else p['remaining_cell_ids']
    b.check(r['released_cell_ids']==expected and len(expected)==5,'Exact five-cell stage differs')
    if r['stage']=='remaining':
        s=read(r['scout_collection']);sr=read(r['scout_review'])
        b.check(s['status']=='dry_donor_susceptibility_scout_complete' and s['panel']==b.pin(path) and s['complete']==5 and
            s['cell_ids']==p['scout_cell_ids'] and s['numerical_checks_passed'] is True,'Actual La five-cell block has not qualified')
        b.check(sr['status']=='passed_dry_donor_susceptibility_scout_review' and sr['collection']==r['scout_collection'] and sr['panel']==b.pin(path) and sr['implementation']==p['adapter'],'Actual scout/cost review differs')
        b.check(score(s['rows'],p,'scout')==s['response'],'Scout numerical replay differs')
        for cid in p['scout_cell_ids']:
            row=s['rows'][cid];b.verify(row['result']);b.verify(row['output']);b.verify(row['properties']['property_file'])
    b.check(Path(r['execution_root']).is_absolute(),'Absolute execution root required');return r


claim=old.claim;stop=old.stop


def run(path,release,cid,out):
    p=load(path); r=release_contract(p,path,json.loads(Path(release).read_text()))
    b.check(cid in r['released_cell_ids'],'Cell not released'); c=next(x for x in p['cells'] if x['cell_id']==cid)
    b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm is required')
    b.check(int(os.environ.get('SLURM_NTASKS','0'))==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Four single-thread MPI ranks required')
    b.check(int(os.environ.get('SLURM_JOB_NUM_NODES','0'))==1,'One node required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or 4*int(os.environ.get('SLURM_MEM_PER_CPU','0'))
    b.check(mem>0,'Actual memory grant required'); maxcore=math.floor(.75*mem/4); b.check(maxcore>0,'Invalid MaxCore')
    root=Path(r['execution_root']); claim(root,out,cid,[x['cell_id'] for x in p['cells']])
    tick=time.monotonic(); code=None; parsed=None; properties=None; error=None; start=None; attempted=0
    try:
        shutil.copyfile(b.verify(c['source_geometry']),out/'geometry.xyz'); (out/'endpoint.inp').write_text(deck(p,c,4,maxcore)); env=os.environ.copy()
        for key in THREADS:env[key]='1'
        env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE'); env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
        command=[p['orca']['path'],'endpoint.inp','--bind-to none']
        start=save(out/'START.json',dict(cell=c,panel=b.pin(path),adapter=p['adapter'],release=b.pin(release),input=b.pin(out/'endpoint.inp'),
            geometry=b.pin(out/'geometry.xyz'),command=command,ranks=4,nodes=1,memory_MiB=mem,maxcore_MiB=maxcore,job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),
            maximum_cell_native_programs=1,threads={k:env[k] for k in THREADS},utc=b.now()))
        with (out/'endpoint.out').open('x') as stdout,(out/'endpoint.err').open('x') as stderr:
            attempted=1; code=subprocess.run(command,cwd=out,env=env,stdout=stdout,stderr=stderr,check=False).returncode
        parsed=parse((out/'endpoint.out').read_text(errors='replace'),c,code)
        b.check(parsed['status']=='complete','Native admission failed')
        properties=field_properties(out/'endpoint.out',out/'endpoint.property.txt',c,parsed)
    except Exception as exc:
        error=repr(exc); stop(root,cid,error)
    elapsed=time.monotonic()-tick
    result=dict(status='complete' if error is None else 'failed_admission_or_execution',adapter=p['adapter'],panel=b.pin(path),
                cell_id=cid,cell=c,release=b.pin(release),start=start,returncode=code,parsed=parsed,properties=properties,
                native_program_invocations=attempted,elapsed_seconds=elapsed,allocated_elapsed_CPU_seconds_proxy=4*elapsed,
                job_id=os.environ['SLURM_JOB_ID'],error=error)
    for name,suffix in [('input','inp'),('output','out'),('stderr','err')]:
        f=out/('endpoint.'+suffix); result[name]=b.pin(f) if f.exists() else None
    save(out/'RESULT.json',result); print(json.dumps({'cell_id':cid,'status':result['status']})); return 0 if error is None else 1


def admit_result(path,c,p,panel):
    path=Path(path); r=json.loads(path.read_text()); b.check(r['status']=='complete' and r['error'] is None and r['returncode']==0,'Native field unavailable')
    b.check(r['cell']==c and r['cell_id']==c['cell_id'] and r['adapter']==p['adapter'] and r['panel']==b.pin(panel),'Native field lineage differs')
    release=read(r['release']); release_contract(p,panel,release)
    b.check(c['cell_id'] in release['released_cell_ids'] and path.resolve()==Path(release['execution_root']).resolve()/'tasks'/c['cell_id']/'execution/RESULT.json','Native attempt path differs')
    start=read(r['start']); b.check(r['start']==b.pin(path.parent/'START.json') and start['cell']==c and start['panel']==r['panel'] and start['release']==r['release'] and start['input']==r['input'],'Native START differs')
    b.check(start['ranks']==4 and start['nodes']==1 and start['memory_MiB']>0 and start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4),'Actual resources differ')
    b.check(start['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'] and start['maximum_cell_native_programs']==r['native_program_invocations']==1,'Native invocation differs')
    b.check(all(start['threads'][k]=='1' for k in THREADS),'Native numerical threading differs')
    b.check(b.verify(start['geometry']).read_bytes()==b.verify(c['source_geometry']).read_bytes(),'Executed geometry differs')
    b.check(b.verify(r['input']).read_text()==deck(p,c,4,start['maxcore_MiB']),'Executed deck differs'); b.verify(r['stderr'])
    q=parse(b.verify(r['output']).read_text(),c,0); b.check(canonical(q)==r['parsed'] and q['status']=='complete','Native raw replay differs')
    props=field_properties(r['output']['path'],b.verify(r['properties']['property_file']),c,q)
    b.check(props==r['properties'],'Native full-precision property replay differs')
    return dict(r,result=b.pin(path))



def compare_model(response,p):
    inv=read(p['source_inventory']);review=read(p['model_actual_review']);result=read(inv['model40_result'])
    b.check(review['result']==inv['model40_result'] and review['status'].startswith('passed_'),'Model actual review/result differs')
    # Producer arithmetic already independently admitted; preserve its exact center estimates and raw authority.
    values={m:p['centers'][m]['model_central_field_estimates'] for m in METALS}
    return dict(model_result=inv['model40_result'],model_review=p['model_actual_review'],model_center_estimates=values,
        quantum_estimates={m:response['per_metal'][m]['estimates'] for m in METALS},
        meaning='Signed model minus quantum susceptibility at both steps; no cross-method accuracy threshold',
        model_minus_quantum={m:{s:values[m][s]['chi_E_a0_cubed']-response['per_metal'][m]['estimates'][s]['chi_energy'] for s in ('h','half')} for m in METALS})


def collect(panel,inventory,out,stage):
    tick,cpu=time.monotonic(),time.process_time();p=load(panel);inv=json.loads(Path(inventory).read_text())
    ids=p['scout_cell_ids'] if stage=='scout' else [c['cell_id'] for c in p['cells']]
    b.check(set(inv)==set(ids),'Exact field inventory denominator required');out.mkdir(parents=True,exist_ok=False)
    rows={cid:dict(status='unavailable',result=None,error='Missing output') for cid in ids}
    for cid in ids:
        source=Path(inv[cid]['path'] if isinstance(inv[cid],dict) else inv[cid])
        try:
            if isinstance(inv[cid],dict):b.verify(inv[cid])
            c=next(c for c in p['cells'] if c['cell_id']==cid);rows[cid]=admit_result(source,c,p,panel)
        except Exception as exc:rows[cid]=dict(status='unavailable',result=b.pin(source) if source.is_file() else None,error=repr(exc))
    error=None;response=None;comparison=None;complete=sum(r['status']=='complete' for r in rows.values())
    try:
        response=score(rows,p,stage)
        if stage=='full' and complete==10:comparison=compare_model(response,p)
    except Exception as exc:error=repr(exc)
    passed=error is None and complete==len(ids) and response['passed']
    r=dict(status='dry_donor_susceptibility_'+('scout' if stage=='scout' else 'reference')+('_complete' if passed else '_unqualified'),
        implementation=b.pin(__file__),panel=b.pin(panel),specification=p['specification'],inventory=b.pin(inventory),stage=stage,cell_ids=ids,
        rows=rows,complete=complete,declared=len(ids),response=response,model_comparison=comparison,numerical_checks_passed=passed,new_native_calls=0,error=error,
        native_programs_in_admitted_rows=sum(r.get('native_program_invocations',0) for r in rows.values()),
        admitted_native_allocated_CPU_seconds_proxy=sum(r.get('allocated_elapsed_CPU_seconds_proxy',0) for r in rows.values()),
        cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,local_threads=1))
    rp=save(out/'RESULT.json',r)
    (out/'REPORT.md').write_text(f'# Dry donor susceptibility\n\nAdmitted {complete}/{len(ids)} native SPs. Numerical consistency: {passed}.\n\n'
        'RESULT retains both steps/estimators, matching SP zeros, electronic/final/D4 stages, all native SCF criteria and COSX/XC metadata, and failed slots. '
        'The older EnGrad zeros are contextual only. Printing bounds and numerical resolutions are not physical error bars. '
        'Failed consistency leaves nonlinearity and numerical error unresolved. No unique permanent/induced/charge-transfer attribution follows.\n')
    print(json.dumps(dict(result=rp,status=r['status'])));return 0 if passed else 1


def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('config');a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('prepare');a.add_argument('--config',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('run');a.add_argument('--panel',required=True);a.add_argument('--release',required=True);a.add_argument('--cell',required=True);a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect');a.add_argument('--panel',required=True);a.add_argument('--inventory',required=True);a.add_argument('--stage',choices=('scout','full'),required=True);a.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.action=='config':b.check(not a.out.exists(),'Preserve config');save(a.out,default_config())
    elif a.action=='prepare':prepare(a.config,a.out)
    elif a.action=='run':return run(a.panel,a.release,a.cell,a.out)
    else:return collect(a.panel,a.inventory,a.out,a.stage)
    return 0
if __name__=='__main__':sys.exit(main())
