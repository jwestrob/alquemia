#!/usr/bin/env python3
"""Sixteen declared fields, compatible SP zeros, and a four-field scout gate."""
import argparse
import copy
import fcntl
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

HERE = Path(__file__).resolve().parent


def module(path, name):
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


common = module(HERE.parent/'common_basis_decomposition_v1/REFERENCE.py', 'field_common')
b = common.b; read = common.read
ROOT = b.ROOT
TRIMER = HERE.parent.parent/'explicit_model_v1/collective_mechanics_options_v1/collection_v1/TRIMER_COMPARE_v2.py'
PROPERTY = HERE.parent/'saved_electronic_redistribution_v1/independent_actual_review_v1/REVIEW_v3.py'
moment = module(TRIMER, 'field_old_moments'); prop = module(PROPERTY, 'field_property_arrays')
PROTOCOL = 'conditional_susceptibility_reference_v1'
FIELDS = {'minus_h': -.001, 'minus_half': -.0005, 'plus_half': .0005, 'plus_h': .001}
BASES = [f'pbe0__aquo9__{m}__{tail}' for m in ('La', 'Dy') for tail in
         ('water_01__metal_water', 'far__metal_two_water')]
SCOUT_BASE = BASES[0]
SCOUTS = [SCOUT_BASE+'__Fz_'+tag for tag in FIELDS]
HELPERS = dict(common_adapter=Path(common.__file__), property_helper=PROPERTY,
               old_moment_helper=TRIMER, **common.HELPERS)
THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
           'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')


def canonical(value):
    return json.loads(json.dumps(value, allow_nan=False))


def token_roundoff(token):
    """Printed decimal half-unit plus conservative binary conversion allowance."""
    mantissa, *exponent = re.split('[EeDd]', token)
    places = len(mantissa.split('.', 1)[1]) if '.' in mantissa else 0
    value = float(token.replace('D','E').replace('d','e'))
    b.check(math.isfinite(value), 'Nonfinite precision token')
    return .5*10.**((int(exponent[0]) if exponent else 0)-places)+2*math.ulp(value)


def save(path, value):
    path = Path(path); tmp = path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n'); os.replace(tmp, path)
    return b.pin(path)


def default_config():
    inv = read(b.pin(HERE/'SOURCE_INVENTORY.json'))
    trimer = read(inv['trimer_geometry_manifest'])
    parent = read(trimer['sources']['common_panel'])
    return dict(protocol=PROTOCOL, implementation=b.pin(__file__), specification=b.pin(HERE/'SPEC.md'),
                source_inventory=b.pin(HERE/'SOURCE_INVENTORY.json'), common_panel=trimer['sources']['common_panel'],
                orca=parent['orca'], run_wrapper=b.pin(HERE.parent.parent/'conditional_susceptibility_v1/REFERENCE.sbatch'),
                **{key:b.pin(path) for key,path in HELPERS.items()})


def config_check(cfg):
    b.check(cfg == default_config(), 'Exact source/configuration differs')
    for value in cfg.values():
        if isinstance(value, dict) and set(value) == {'path','sha256'}: b.verify(value)


def base_deck(row, ranks, maxcore):
    text = b.verify(row['input']).read_text()
    text, n = re.subn(r'^%maxcore \d+\s*$', '%maxcore '+str(maxcore), text, flags=re.M)
    b.check(n == 1, 'Resource MaxCore denominator differs')
    text, n = re.subn(r'^%pal nprocs \d+ end\s*$', f'%pal nprocs {ranks} end', text, flags=re.M)
    b.check(n == 1, 'Resource ranks denominator differs'); return text


def deck(p, c, ranks, maxcore):
    text = base_deck(p['zero_field_records'][c['zero_id']], ranks, maxcore)
    b.check(text.count('%scf\n') == 1 and '%elprop' not in text and 'EField' not in text,
            'Original field/property recipe is not the declared zero')
    text = text.replace('%scf\n', f'%scf\n EField 0.0,0.0,{c["field_au"]:.7f}\n EFieldOrigin 0.0,0.0,0.0\n', 1)
    return text.replace('* xyz ', '%elprop\n Dipole true\n Origin 0.0,0.0,0.0\nend\n* xyz ', 1)


def cells(records):
    out = []
    for gid in BASES:
        for tag, field in FIELDS.items():
            c = copy.deepcopy(records[gid]['cell']); c.pop('template_input', None)
            c.update(cell_id=gid+'__Fz_'+tag, zero_id=gid, field_tag=tag, field_au=field)
            out.append(c)
    return out


def zero_replay(record, parent):
    r = read(record['result']); c = record['cell']; old = r.get('parsed', r)
    b.check(r['status'] == 'complete' and r['returncode'] == 0 and r['cell'] == c,
            'Compatible zero result unavailable or changed')
    b.check(r['input'] == record['input'] and r['output'] == record['output'], 'Zero input/output join differs')
    b.verify(r['adapter']); start = read(r['start'])
    b.check(start['input'] == r['input'] and start['ranks'] == 4, 'Zero START join differs')
    b.check(base_deck(record, 4, 3000) == common.deck(parent, c, 4, 3000), 'Zero recipe differs from original renderer')
    q = common.parse(b.verify(r['output']).read_text(), c, 0)
    b.check(q['status'] == 'complete' and all(canonical(old[k]) == canonical(v) for k,v in q.items()), 'Zero raw admission differs')
    mu = moment.quantum_moment(r['output'], 3, prop)
    b.check(mu['property_file'] == record['property_file'] and mu['total_dipole_au'] == record['dipole_at_metal_origin_au'], 'Zero full-precision moment differs')
    b.check(q['electronic_energy_hartree'] == record['E0_electronic_hartree'] and q['D4_energy_hartree'] == 0., 'Zero energy/D4 differs')
    et=re.findall(r'^Total Energy\s*:\s*('+moment.NUMBER+r')\s+Eh',b.verify(r['output']).read_text(),re.M)
    b.check(len(et)==1,'Zero energy precision token denominator differs')
    return dict(result=record['result'], energy_hartree=q['electronic_energy_hartree'],
                dipole_au=mu['total_dipole_au'], native=q, moment=mu,
                precision=dict(energy_hartree=token_roundoff(et[0]),
                    moment_z_au=mu['print_bounds_au']['total_dipole_au'],
                    meaning='Decimal/binary bounds only; zero moment conservatively retains original printed-origin conversion bound; excludes SCF/method error'))


def prepare(config, out):
    tick, cpu = time.monotonic(), time.process_time(); cfg = json.loads(Path(config).read_text()); config_check(cfg)
    inv = read(cfg['source_inventory']); parent = read(cfg['common_panel']); records = inv['zero_field_records']
    b.check(set(records) == set(BASES) and inv['specification'] == cfg['specification'], 'Four-system inventory differs')
    zeros = {gid:zero_replay(records[gid], parent) for gid in BASES}
    for m in ('La','Dy'):
        one, two = [records[x]['cell'] for x in BASES if records[x]['metal'] == m]
        b.check(one['positions_A'] == two['positions_A'] and one['symbols'] == two['symbols'], 'Parent frames differ')
        b.check(one['expected_basis_functions_by_center'] == two['expected_basis_functions_by_center'] and
                one['expected_aux_J_basis_functions'] == two['expected_aux_J_basis_functions'], 'Parent bases differ')
        b.check(one['real_parent_atom_indices'] == [0,1,2,3] and two['real_parent_atom_indices'] == [0,1,2,3,7,8,9], 'Fixed far selection differs')
    out.mkdir(parents=True, exist_ok=False); (out/'inputs').mkdir()
    p = dict(protocol=PROTOCOL, adapter=cfg['implementation'], configuration=b.pin(config), specification=cfg['specification'],
             source_inventory=cfg['source_inventory'], common_panel=cfg['common_panel'], orca=cfg['orca'],
             zero_field_records=records, zero_replays=zeros, cells=cells(records), maximum_native_calls=16,
             ranks=4, scout_cell_ids=SCOUTS, run_wrapper=cfg['run_wrapper'], **{k:cfg[k] for k in HELPERS})
    for c in p['cells']:
        f=out/'inputs'/(c['cell_id']+'.inp'); f.write_text(deck(p,c,4,3000)); c['template_input']=b.pin(f)
    pp=save(out/'PANEL.json',p); v=validate(pp['path']); v.update(wall_seconds=time.monotonic()-tick,
        parent_cpu_seconds=time.process_time()-cpu, local_threads=1, zero_raw_replays=4, zero_full_precision_moment_replays=4)
    vp=save(out/'PREFLIGHT.json',v); print(json.dumps(dict(panel=pp,preflight=vp,status=v['status'],new_native_calls=0)))


def load(path):
    p=json.loads(Path(path).read_text()); cfg=read(p['configuration']); config_check(cfg)
    b.check(p['protocol']==PROTOCOL and p['adapter']==cfg['implementation'] and p['specification']==cfg['specification'] and
            p['source_inventory']==cfg['source_inventory'] and p['common_panel']==cfg['common_panel'] and p['orca']==cfg['orca'] and
            p['run_wrapper']==cfg['run_wrapper'], 'Panel source/configuration differs')
    for k in HELPERS: b.check(p[k]==cfg[k], 'Executing helper differs: '+k)
    b.check(p['zero_field_records']==read(cfg['source_inventory'])['zero_field_records'] and
            set(p['zero_replays'])==set(BASES) and p['maximum_native_calls']==len(p['cells'])==16 and
            p['ranks']==4 and p['scout_cell_ids']==SCOUTS, 'Finite inventory differs')
    for actual, expected in zip(p['cells'], cells(p['zero_field_records'])):
        b.check({k:v for k,v in actual.items() if k!='template_input'}==expected, 'Field definition differs')
        b.check(b.verify(actual['template_input']).read_text()==deck(p,actual,4,3000), 'Field template differs')
    return p


def validate(path):
    p=load(path); checks={}
    for c in p['cells']:
        text=b.verify(c['template_input']).read_text(); original=base_deck(p['zero_field_records'][c['zero_id']],4,3000)
        removed=text.replace(f' EField 0.0,0.0,{c["field_au"]:.7f}\n EFieldOrigin 0.0,0.0,0.0\n','',1).replace('%elprop\n Dipole true\n Origin 0.0,0.0,0.0\nend\n','',1)
        b.check(removed==original,'More than declared field/property/resource delta')
        b.check(c['charge']==3 and c['multiplicity']==1 and c['expected_basis_functions']==434 and
                c['expected_electrons'] in (18,28) and c['d4'] is False and 'EnGrad' not in text,'Reference scope differs')
        checks[c['cell_id']]=dict(sole_declared_delta=True,field_au=c['field_au'],zero_result=p['zero_field_records'][c['zero_id']]['result'])
    return dict(status='passed_static_preflight',panel=b.pin(path),specification=p['specification'],declared=16,
                zero_reuse=4,scout_cell_ids=SCOUTS,checks=checks,new_native_calls=0,native_field_admission_pending=True)


def scalar(block, name):
    hits=re.findall(r'&'+re.escape(name)+r'\s+\[&Type "(?:Double|Integer)"[^\]]*\]\s+('+moment.NUMBER+r')',block)
    b.check(len(hits)==1,'Missing/duplicate property scalar '+name); value=float(hits[0])
    b.check(math.isfinite(value),'Nonfinite property scalar'); return value


def field_properties(output, property_path, c, q):
    text=Path(output).read_text(); raw=Path(property_path).read_text()
    b.check('ORCA 6.1.1' in raw,'Property version differs')
    origins=re.findall(r'^Position of electric origin\s*\.\.\.\s*('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s*$',text,re.M)
    b.check(len(origins)==1 and list(map(float,origins[0]))==[0.,0.,0.],'Native property origin differs')
    geometry=prop.section(raw,'Geometry'); b.check('&Units "Bohr"' in geometry,'Property coordinate units differ')
    tail=geometry.split('&Units "Bohr"]',1)[1]; rows=re.findall(r'^\s*([A-Za-z]+)(:?)\s+('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s+('+moment.NUMBER+r')\s*$',tail,re.M)
    b.check(len(rows)==len(c['symbols'])==28 and [v[0] for v in rows]==c['symbols'] and
            [bool(v[1]) for v in rows]==c['ghost_mask'],'Property geometry/ghost inventory differs')
    positions=[list(map(float,v[2:])) for v in rows]
    b.check(max(abs(x*moment.BOHR-y) for a,z in zip(positions,c['positions_A']) for x,y in zip(a,z))<=1e-6,'Property positions differ')
    info=prop.section(raw,'Calculation_Info')
    for key,expected in {'Charge':3,'Mult':1,'NumOfAtoms':28,'NumOfElectrons':c['expected_electrons'],
                         'NumOfBasisFuncts':434,'NumOfAuxJBasisFuncts':c['expected_aux_J_basis_functions']}.items():
        b.check(scalar(info,key)==expected,'Property metadata differs: '+key)
    block=prop.section(raw,'SCF_Dipole_Moment'); values={}
    for key,name,label in [('total','dipoleTotal','Total Dipole Moment'),('nuclear','dipoleNucContrib','Nuclear contribution'),('electronic','dipoleElecContrib','Electronic contribution')]:
        values[key]=prop.array(block,name,3)
        b.check(moment.difference(values[key],moment.vector(text,label))<=5e-10,'Full-precision dipole does not round to output')
    expected=[math.fsum(z*x[j] for z,x in zip(c['expected_effective_nuclear_charges'],positions)) for j in range(3)]
    b.check(moment.difference(values['nuclear'],expected)<=1e-8,'Effective-ECP nuclear dipole/lab origin differs')
    b.check(moment.difference(values['total'],[a+z for a,z in zip(values['nuclear'],values['electronic'])])<=5e-9,'Dipole component closure differs')
    ep=scalar(prop.section(raw,'Single_Point_Data'),'FinalEnergy')
    es=prop.array(prop.section(raw,'SCF_Energy'),'totalEnergy',1)[0]
    b.check(abs(ep-q['energy_hartree'])<=1e-10 and abs(es-ep)<=1e-10,'Property/output energy stages differ')
    b.check('&Converged [&Type "Boolean"] true' in prop.section(raw,'Single_Point_Data'),'Property convergence absent')
    et=re.findall(r'&FinalEnergy\s+\[&Type "Double"\]\s+('+moment.NUMBER+r')',prop.section(raw,'Single_Point_Data'))
    tail=block.split('&dipoleTotal ',1)[1].split('\n   &',1)[0]
    mt=re.findall(r'^\s*(\d+)\s+('+moment.NUMBER+r')\s*$',tail,re.M)
    b.check(len(et)==1 and [int(i) for i,_ in mt]==[0,1,2],'Property precision token denominator differs')
    return dict(property_file=b.pin(property_path),dipole_au=values['total'],dipole_components_au=values,
                origin_Bohr=[0.,0.,0.],energy_hartree=ep,SCF_energy_hartree=es,positions_Bohr=positions,
                precision=dict(energy_hartree=token_roundoff(et[0]),moment_z_au=token_roundoff(mt[2][1]),
                    meaning='Full-precision property decimal/binary bounds only; excludes SCF/method error'))


def parse(text,c,returncode):
    q=common.parse(text,c,returncode)
    echo=[re.sub(r'^\s*\|\s*\d+> ?','',x) for x in text.split('****END OF INPUT****')[0].splitlines() if re.match(r'^\s*\|\s*\d+>',x)]
    expected=f' EField 0.0,0.0,{c["field_au"]:.7f}'
    q['checks']['declared_field_echo']=echo.count(expected)==1 and sum(bool(re.search(r'\bEField\b',x)) for x in echo)==1
    q['checks']['declared_origins_echo']=echo.count(' EFieldOrigin 0.0,0.0,0.0')==1 and echo.count(' Origin 0.0,0.0,0.0')==1
    q['status']='complete' if all(q['checks'].values()) else 'failed_admission'; return q


def estimate(rows,zero):
    out={}
    for label,s,minus,plus in [('h',.001,'minus_h','plus_h'),('half',.0005,'minus_half','plus_half')]:
        a,z=rows[minus],rows[plus]; e0=zero['energy_hartree']; mu0=zero['dipole_au'][2]
        vals=[a['energy_hartree'],z['energy_hartree'],e0,a['dipole_au'][2],z['dipole_au'][2],mu0]
        b.check(all(math.isfinite(x) for x in vals),'Nonfinite response data')
        out[label]=dict(chi_energy=-(z['energy_hartree']+a['energy_hartree']-2*e0)/s**2,
                        chi_moment=(z['dipole_au'][2]-a['dipole_au'][2])/(2*s),
                        odd_energy_moment=-(z['energy_hartree']-a['energy_hartree'])/(2*s),
                        zero_moment=mu0,even_moment_numerator=z['dipole_au'][2]+a['dipole_au'][2]-2*mu0)
        ae,ze,e0b=[x['precision']['energy_hartree'] for x in (a,z,zero)]
        am,zm,m0b=[x['precision']['moment_z_au'] for x in (a,z,zero)]
        out[label]['roundoff_bounds']=dict(chi_energy=(ae+ze+2*e0b)/s**2,
            chi_moment=(am+zm)/(2*s),odd_energy_moment=(ae+ze)/(2*s),zero_moment=m0b,
            odd_minus_zero=(ae+ze)/(2*s)+m0b,even_moment_numerator=am+zm+2*m0b)
    return out


def consistency(est,limit,odd=False):
    checks={f'energy_vs_moment_{s}':abs(est[s]['chi_energy']-est[s]['chi_moment'])<=limit for s in ('h','half')}
    checks.update({f'{key}_step':abs(est['h'][key]-est['half'][key])<=limit for key in ('chi_energy','chi_moment')})
    if odd: checks.update({f'odd_slope_{s}':abs(est[s]['odd_energy_moment']-est[s]['zero_moment'])<=1e-4 for s in ('h','half')})
    return dict(estimates=est,checks=checks,passed=all(checks.values()),limit_a0_cubed=limit)


def combine(estimates,coefficients):
    out={}
    for s in ('h','half'):
        first=next(iter(estimates.values()))[s]
        out[s]={key:math.fsum(coefficients[gid]*estimates[gid][s][key] for gid in coefficients)
                for key in first if key!='roundoff_bounds'}
        out[s]['roundoff_bounds']={key:math.fsum(abs(coefficients[gid])*estimates[gid][s]['roundoff_bounds'][key]
                                               for gid in coefficients) for key in first['roundoff_bounds']}
    return out


def score(rows,p,stage):
    wanted=[SCOUT_BASE] if stage=='scout' else BASES; states={gid:None for gid in wanted}; errors={}
    for gid in wanted:
        try:
            rr={tag:rows[gid+'__Fz_'+tag]['properties'] for tag in FIELDS}
            states[gid]=consistency(estimate(rr,p['zero_replays'][gid]),.5,True)
        except (KeyError,TypeError,ValueError) as exc: errors[gid]=repr(exc)
    increments={m:None for m in ('La','Dy')}; paired=None
    if stage=='full' and all(states.values()):
        estimates={gid:r['estimates'] for gid,r in states.items()}
        for m in ('La','Dy'):
            one,two=[gid for gid in BASES if p['zero_field_records'][gid]['metal']==m]
            increments[m]=consistency(combine(estimates,{two:1,one:-1}),1.)
        paired=consistency(combine({m:r['estimates'] for m,r in increments.items()},{'Dy':1,'La':-1}),1.)
    passed=all(r is not None and r['passed'] for r in states.values())
    if stage=='full':passed=passed and all(r is not None and r['passed'] for r in increments.values()) and paired is not None and paired['passed']
    return dict(per_state=states,increments=increments,paired=paired,errors=errors,numerical_checks_passed=bool(passed),
                interpretation='Conditional susceptibility increment, not alpha3; failed consistency cause unresolved')


def release_contract(p,path,r):
    b.check(r['status']=='approved_conditional_susceptibility_reference' and r['implementation']==p['adapter'] and
            r['panel']==b.pin(path) and r['specification']==p['specification'],'Root source release differs')
    ids=[c['cell_id'] for c in p['cells']]
    b.check(r['cell_ids']==ids and r['maximum_native_calls']==16 and r['ranks']==4,'Released finite scope differs')
    prep=read(r['preparation']); review=read(r['independent_review'])
    b.check(prep['status']=='passed_static_preflight' and prep['panel']==b.pin(path),'Preparation unavailable')
    b.check(review['status']=='passed_conditional_susceptibility_reference_review' and review['implementation']==p['adapter'] and review['panel']==b.pin(path),'Independent source/preparation review differs')
    b.check(r['stage'] in ('scout','remaining'),'Unknown release stage')
    expected=SCOUTS if r['stage']=='scout' else [x for x in ids if x not in SCOUTS]
    b.check(r['released_cell_ids']==expected,'Exact stage inventory differs')
    if r['stage']=='remaining':
        s=read(r['scout_collection']); sr=read(r['scout_review'])
        b.check(s['status']=='conditional_susceptibility_scout_complete' and s['panel']==b.pin(path) and
                s['complete']==4 and s['cell_ids']==SCOUTS and s['numerical_checks_passed'] is True,'Actual four-field scout has not qualified')
        b.check(sr['status']=='passed_conditional_susceptibility_scout_review' and sr['collection']==r['scout_collection'] and
                sr['panel']==b.pin(path) and sr['implementation']==p['adapter'],'Actual scout review differs')
        b.check(score(s['rows'],p,'scout')==s['response'],'Saved scout numerical replay differs')
        for cid in SCOUTS:
            row=s['rows'][cid]; b.verify(row['result']); b.verify(row['output']); b.verify(row['properties']['property_file'])
    b.check(Path(r['execution_root']).is_absolute(),'Absolute execution root required'); return r


def claim(root,out,cid,ids):
    b.check(out.resolve()==root.resolve()/'tasks'/cid/'execution','Exact one-attempt output differs')
    root.mkdir(parents=True,exist_ok=True)
    with (root/'.native.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX); b.check(not (root/'STOP.json').exists(),'Reference stop is set')
        lp=root/'NATIVE_LEDGER.json'; ledger=json.loads(lp.read_text()) if lp.exists() else []
        b.check(len(ledger)==len(set(ledger)) and set(ledger)<=set(ids) and len(ledger)<16 and cid not in ledger,'Native cap or duplicate attempt')
        out.mkdir(parents=True,exist_ok=False); save(lp,ledger+[cid])


def stop(root,cid,error):
    with (root/'.native.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if not (root/'STOP.json').exists():save(root/'STOP.json',dict(cell_id=cid,error=error,utc=b.now()))


def run(path,release,cid,out):
    p=load(path); r=release_contract(p,path,json.loads(Path(release).read_text()))
    b.check(cid in r['released_cell_ids'],'Cell not released'); c=next(x for x in p['cells'] if x['cell_id']==cid)
    b.check(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm is required')
    b.check(int(os.environ.get('SLURM_NTASKS','0'))==4 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Four single-thread MPI ranks required')
    mem=int(os.environ.get('SLURM_MEM_PER_NODE','0')) or 4*int(os.environ.get('SLURM_MEM_PER_CPU','0'))
    b.check(mem>0,'Actual memory grant required'); maxcore=math.floor(.75*mem/4); b.check(maxcore>0,'Invalid MaxCore')
    root=Path(r['execution_root']); claim(root,out,cid,[x['cell_id'] for x in p['cells']])
    tick=time.monotonic(); code=None; parsed=None; properties=None; error=None; start=None; attempted=0
    try:
        (out/'endpoint.inp').write_text(deck(p,c,4,maxcore)); env=os.environ.copy()
        for key in THREADS:env[key]='1'
        env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE'); env['PATH']=str(b.verify(p['orca']).parent)+os.pathsep+env.get('PATH','')
        command=[p['orca']['path'],'endpoint.inp','--bind-to none']
        start=save(out/'START.json',dict(cell=c,panel=b.pin(path),adapter=p['adapter'],release=b.pin(release),input=b.pin(out/'endpoint.inp'),
            command=command,ranks=4,memory_MiB=mem,maxcore_MiB=maxcore,job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),
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
    b.check(start['ranks']==4 and start['memory_MiB']>0 and start['maxcore_MiB']==math.floor(.75*start['memory_MiB']/4),'Actual resources differ')
    b.check(start['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'] and start['maximum_cell_native_programs']==r['native_program_invocations']==1,'Native invocation differs')
    b.check(all(start['threads'][k]=='1' for k in THREADS),'Native numerical threading differs')
    b.check(b.verify(r['input']).read_text()==deck(p,c,4,start['maxcore_MiB']),'Executed deck differs'); b.verify(r['stderr'])
    q=parse(b.verify(r['output']).read_text(),c,0); b.check(canonical(q)==r['parsed'] and q['status']=='complete','Native raw replay differs')
    props=field_properties(r['output']['path'],b.verify(r['properties']['property_file']),c,q)
    b.check(props==r['properties'],'Native full-precision property replay differs')
    return dict(r,result=b.pin(path))


def collect(panel,inventory,out,stage):
    tick,cpu=time.monotonic(),time.process_time(); p=load(panel); inv=json.loads(Path(inventory).read_text())
    ids=SCOUTS if stage=='scout' else [c['cell_id'] for c in p['cells']]
    b.check(set(inv)==set(ids),'Exact field inventory denominator required'); out.mkdir(parents=True,exist_ok=False)
    rows={cid:dict(status='unavailable',result=None,error='Missing output') for cid in ids}
    for cid in ids:
        source=Path(inv[cid]['path'] if isinstance(inv[cid],dict) else inv[cid])
        try:
            if isinstance(inv[cid],dict):b.verify(inv[cid])
            c=next(c for c in p['cells'] if c['cell_id']==cid); row=admit_result(source,c,p,panel); rows[cid]=row
        except Exception as exc:
            rows[cid]=dict(status='unavailable',result=b.pin(source) if source.is_file() else None,error=repr(exc))
    response=score(rows,p,stage); complete=sum(r['status']=='complete' for r in rows.values()); passed=complete==len(ids) and response['numerical_checks_passed']
    r=dict(status='conditional_susceptibility_'+('scout' if stage=='scout' else 'reference')+('_complete' if passed else '_unqualified'),
           implementation=b.pin(__file__),panel=b.pin(panel),specification=p['specification'],inventory=b.pin(inventory),stage=stage,cell_ids=ids,
           rows=rows,complete=complete,declared=len(ids),response=response,numerical_checks_passed=passed,new_native_calls=0,
           native_programs_in_admitted_rows=sum(r.get('native_program_invocations',0) for r in rows.values()),
           admitted_native_allocated_CPU_seconds_proxy=sum(r.get('allocated_elapsed_CPU_seconds_proxy',0) for r in rows.values()),
           cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,local_threads=1))
    rp=save(out/'RESULT.json',r)
    (out/'REPORT.md').write_text(f'# Conditional reference susceptibility\n\nNative fields admitted: {complete}/{len(ids)}. Numerical consistency: {passed}.\n\n'
        'All declared slots and failed checks are retained in RESULT.json. Consistency limits are not physical error bars. '
        'A failed step/energy–dipole check leaves finite-step nonlinearity and numerical error unresolved. '
        'The conditional increment is not alpha3 or a unique mechanism. Original SCF/state/basis limitations remain.\n')
    print(json.dumps(dict(result=rp,status=r['status']))); return 0 if passed else 1


def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('config'); a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('prepare'); a.add_argument('--config',required=True); a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('run'); a.add_argument('--panel',required=True); a.add_argument('--release',required=True); a.add_argument('--cell',required=True); a.add_argument('--out',type=Path,required=True)
    a=sub.add_parser('collect'); a.add_argument('--panel',required=True); a.add_argument('--inventory',required=True); a.add_argument('--stage',choices=('scout','full'),required=True); a.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.action=='config':b.check(not a.out.exists(),'Preserve existing config');save(a.out,default_config())
    elif a.action=='prepare':prepare(a.config,a.out)
    elif a.action=='run':return run(a.panel,a.release,a.cell,a.out)
    else:return collect(a.panel,a.inventory,a.out,a.stage)
    return 0


if __name__=='__main__':sys.exit(main())
