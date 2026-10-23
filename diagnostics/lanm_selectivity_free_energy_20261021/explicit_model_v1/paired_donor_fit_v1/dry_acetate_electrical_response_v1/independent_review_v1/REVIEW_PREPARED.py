#!/usr/bin/env python3
"""Literal source/input/property joins only; no adapter or control-suite execution."""
import hashlib,importlib.util,json,math,re,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;BASE=HERE.parent
ROOT=next(p for p in HERE.parents if (p/'diagnostics').is_dir() and (p/'workspaces').is_dir())
OUT=ROOT/'workspaces'/BASE.relative_to(ROOT/'diagnostics')/'independent_review_v1'
checks={};tick,cpu=time.monotonic(),time.process_time()
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def check(k,v):checks[k]=bool(v);assert v,k
def verify(p):check('pin/'+p['path'],pin(p['path'])==p);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
ready=json.loads((BASE/'READY.json').read_text())
for k,v in ready.items():
    if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
cfg=read(ready['inputs']);m=read(ready['manifest']);inv=read(ready['shared_inventory'])
for k,v in cfg.items():
    if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
for v in cfg['native_source'].values():verify(v)
check('source_preparation_join',m['implementation']==ready['implementation'] and m['inputs']==ready['inputs'] and m['shared_specification']==ready['shared_specification'] and m['shared_inventory']==ready['shared_inventory'])
fields={'minus_h':-.001,'minus_half':-.0005,'plus_half':.0005,'plus_h':.001};direction=inv['geometry_defined_field_direction'];u=m['conversion']
check('native_conversion',u['FIELD_FACTOR']==u['HARTREE']*u['ELEFIELD']/(u['ELECTRIC']*u['BOHR']) and u['FIELD_FACTOR']==5142.206758812493)
check('fixed_direction',direction==m['direction'] and abs(sum(x*x for x in direction)-1)<1e-15)
check('exact_finite_inventory',len(m['rows'])==10 and set(m['jobs'])=={g+'__'+f for g in m['rows'] for f in fields} and len(m['jobs'])==m['maximum_native_calls']==40 and m['workers']==1)
check('central_eight_first',list(m['jobs'])[:8]==[g+'__'+f for g in ('La_central','Dy_central') for f in fields])
component=read(cfg['component_manifest']);actual=read(cfg['component_result']);actualreview=read(cfg['component_actual_review']);check('actual_zero_authority',actualreview['result']==cfg['component_result'] and actualreview['numerical_checks_passed'])
qpanel=read(inv['authorities']['reference_panel']);qcollection=read(inv['authorities']['reference_derivative_result'])
helperpin=inv['authorities']['property_origin_helper'];sp=importlib.util.spec_from_file_location('frozen_property_format',verify(helperpin));helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)
for gid,row in m['rows'].items():
    shared=inv['rows'][gid];check(gid+'/shared_input',row['input']==shared['model_input'] and row['zero_authority']==shared['model_zero_energy_authority'])
    if row['q_A']==0:
        c=component['central'][row['metal']];check(gid+'/central_zero',row['input']==c['source']['effective_input'] and row['zero_energy']==c['energy'] and shared['quantum_result']==qpanel['central'][row['metal']]['saved_result']);qgeom=c['geometry']['xyz']
    else:
        c=component['rows'][gid];check(gid+'/rotation_zero',row['input']==c['input'] and row['zero_energy']==actual['rows'][gid]['parsed'] and shared['quantum_result']==qcollection['slots'][gid]['result']);qgeom=c['source_geometry']
    typed=[x.split() for x in verify(row['input']['xyz']).read_text().splitlines()[1:] if x.strip()]
    qr=read(shared['quantum_result']);ordinary=[x.split() for x in verify(qr['geometry']).read_text().splitlines()[2:] if x.strip()]
    check(gid+'/quantum_geometry_pin',qr['geometry']['sha256']==qgeom['sha256'])
    check(gid+'/literal_model_quantum_atoms',len(typed)==len(ordinary)==8 and [x[1] for x in typed]==[x[0] for x in ordinary] and [[float(x) for x in t[2:5]] for t in typed]==[[float(x) for x in t[1:4]] for t in ordinary])
    check(gid+'/origin_charge',row['charge_origin']['origin_A']==[0.,0.,0.] and float(row['charge_origin']['charge_e'])==2 and abs(sum(float(x) for x in row['charge_origin']['site_monopoles_e'])-2)<1e-14)
    raw=verify(shared['property_file']).read_text();geo=helper.section(raw,'Geometry');check(gid+'/property_units','&Units "Bohr"' in geo)
    atoms=re.findall(r'^\s*([A-Za-z]+)\s+('+helper.N+r')\s+('+helper.N+r')\s+('+helper.N+r')\s*$',geo,re.M);check(gid+'/property_symbols',[a[0] for a in atoms]==[a[0] for a in ordinary])
    xyz=[[float(x) for x in a[1:]] for a in atoms];dip=helper.section(raw,'SCF_Dipole_Moment');total=helper.array(dip,'dipoleTotal',3);nuc=helper.array(dip,'dipoleNucContrib',3);elec=helper.array(dip,'dipoleElecContrib',3);z=[11,6,6,8,8,1,1,1]
    origin=[(math.fsum(z[i]*xyz[i][k] for i in range(8))-nuc[k])/42 for k in range(3)];lab=[total[k]+2*origin[k] for k in range(3)]
    check(gid+'/total_component_sum',max(abs(total[k]-nuc[k]-elec[k]) for k in range(3))<1e-11)
    check(gid+'/full_precision_origin',max(abs(origin[k]-row['quantum']['origin_Bohr'][k]) for k in range(3))<1e-12 and max(abs(lab[k]-row['quantum']['total_lab_au'][k]) for k in range(3))<1e-12)
    check(gid+'/fixed_projection',abs(math.fsum(a*b for a,b in zip(lab,direction))-row['quantum_projected_au'])<1e-12)
    check(gid+'/native_counters',all(v==row['zero_energy']['interactions'][k]+(8 if k in ('Atomic Multipoles','Polarization') else 0) for k,v in row['expected_field_interactions'].items()))
for jid,job in m['jobs'].items():
    row=m['rows'][job['geometry_id']];f=fields[job['field_label']]
    check(jid+'/field_vector',job['signed_field_au']==f and job['field_MV_cm']==[f*u['FIELD_FACTOR']*x for x in direction] and job['field_au']==[f*x for x in direction])
    for k in ('xyz','parameters'):
        verify(job['input'][k]);check(jid+'/'+k+'_unchanged',job['input'][k]['sha256']==row['input'][k]['sha256'])
    before=verify(row['input']['key']).read_text();after=verify(job['input']['key']).read_text();record='EXTERNAL-FIELD '+' '.join(format(f*u['FIELD_FACTOR']*x,'.17g') for x in direction)
    check(jid+'/field_only_key_delta',after==before+'\n# Fixed dry-acetate field-energy discriminator; original donor17621\n'+record+'\n')
    check(jid+'/one_field_record',len(re.findall(r'(?im)^\s*EXTERNAL-FIELD\s+',after))==1)
controls=read(ready['saved_controls']);prep=read(ready['preparation_result'])
check('author39_controls_join',controls['implementation']==ready['implementation'] and controls['manifest']==ready['manifest'] and controls['check_count']==39 and all(controls['checks'].values()) and controls['actual_new_native_calls']==0)
check('preparation40_10zeros_join',prep['implementation']==ready['implementation'] and prep['manifest']==ready['manifest'] and prep['cells']==40 and prep['zero_energies_reused']==prep['quantum_dipoles_replayed']==10 and prep['new_native_calls']==0)
check('same_previously_field_qualified_binary',m['analyze']==cfg['native_source']['analyze'] and m['analyze']['sha256']=='f67420e79f8c1793bf491569f77ab4354219926fe5b4c49a2311e64f217f7a28')
OUT.mkdir(parents=True,exist_ok=True);result=dict(status='passed_independent_field_preparation_checks',implementation=ready['implementation'],manifest=ready['manifest'],reviewer_source=pin(__file__),checks=checks,check_count=len(checks),new_native_calls=0,author_controls_rerun=False,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu));(OUT/'CHECKS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(checks=len(checks),result=pin(OUT/'CHECKS.json')),indent=2))
