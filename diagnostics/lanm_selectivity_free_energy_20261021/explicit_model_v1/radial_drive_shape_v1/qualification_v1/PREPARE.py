#!/usr/bin/env python3
"""Literal49-input preparation at admitted feca; no molecular or optimizer call."""
import argparse,copy,hashlib,importlib.util,json,math,shutil,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
FILES=(('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))
DIMERS=tuple('water_central_'+m for m in ('La','Dy'))
CENTERS=tuple('aquo9_central_'+m for m in ('La','Dy'))
TRIMERS=tuple('aquo9_central_'+m+'__far__metal_two_water' for m in ('La','Dy'))
METAL=DIMERS+CENTERS+TRIMERS
BYPASS=('ion_La','ion_Dy','aquo9_central_La__near__two_water')
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):require(pin(q['path'])==q,'Changed pin '+q['path']);return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def inventory():
    jobs=[]
    def add(stage,g,program,**extra):
        job=dict(stage=stage,cell=g,program=program,fd=False,profile='active' if stage not in ('absent','zero') else stage,**extra)
        job['id']='__'.join([stage,g,program]+([job['kind'],str(job['h']),str(job['sign'])] if 'q' in job else [str(job['sign'])] if 'field_au' in job else []));jobs.append(job)
    for g in DIMERS:
        add('absent',g,'analyze');add('absent',g,'testgrad')
    for g in METAL:add('zero',g,'analyze')
    for g in DIMERS+CENTERS:add('zero',g,'testgrad')
    for g in METAL:
        add('active',g,'analyze');add('active',g,'testgrad');jobs[-1]['fd']=g in DIMERS+TRIMERS
    for kind,cells in (('metal_x',CENTERS),('water_rotation_x',DIMERS)):
        for g in cells:
            for h in (1e-4,5e-5):
                for sign in (-1,1):add('directions',g,'analyze',kind=kind,h=h,sign=sign,q=sign*h)
    for g in DIMERS:
        for sign in (-1,1):add('fields',g,'analyze',field_au=sign*.001,sign=sign)
    for g in BYPASS:add('bypass',g,'analyze')
    require(len(jobs)==49 and len({j['id'] for j in jobs})==49,'49 unique jobs')
    return jobs

def path_geometry(row,kind,q):
    xyz=np.asarray(row['positions_A'],float);out=xyz.copy();tangent=np.zeros_like(xyz);metal=next(i for i,t in enumerate(row['types']) if t in (500,501))
    if kind=='metal_x':out[metal,0]+=q;tangent[metal,0]=1.
    else:
        require(kind=='water_rotation_x' and len(xyz)==4,'Fixed dimer rotation')
        oxy=next(i for i,t in enumerate(row['types']) if t==36);axis=np.array([1.,0.,0.]);ell=float(np.linalg.norm(np.cross(axis,xyz[oxy]-xyz[metal])));require(ell>0,'Nonzero tangent scale')
        theta=q/ell;c=math.cos(theta);s=math.sin(theta);rot=np.array([[1.,0.,0.],[0.,c,-s],[0.,s,c]])
        for i in range(len(xyz)):
            if i!=metal:out[i]=xyz[metal]+rot@(xyz[i]-xyz[metal]);tangent[i]=np.cross(axis,xyz[i]-xyz[metal])/ell
    return out.tolist(),tangent.tolist()

def render(row,job,units):
    key=row['base_key']
    if job['profile']!='absent':key+='POLAR-RADIAL-DRIVE '+str(row['metal_type'])+' 36 '+format(0. if job['profile']=='zero' else row['e'],'.17g')+'\n'
    if job['program']=='analyze' and job['stage']!='directions':key+='RESPONSE-MATRIX-EXPORT\n'
    if 'field_au' in job:key+='EXTERNAL-FIELD 0 0 '+format(job['field_au']*units['MV_cm_per_field_au'],'.17g')+'\n'
    return key

def prepare(config,out):
    t=time.monotonic();cpu=time.process_time();s=config['sources'];[verify(p) for p in s.values()]
    inp=read(s['radial_inputs'])['sources'];cfg=read(inp['parent_evaluator_configuration']);ev=module(inp['inherited_evaluator'],'radial_old_evaluator');resolver=module(cfg['resolver'],'radial_resolver')
    parent=read(inp['inherited_manifest']);qual=read(inp['parent_qualification_manifest']);selected=read(inp['selected']);domain=read(inp['inherited_domain']);values=selected['parameters'];require(selected['status']=='complete' and len(values)==12,'Admitted final feca')
    old=module(inp['parent_qualification_runner'],'radial_old_qualification');h=old.helpers(qual['sources']);rows={}
    for g in METAL+BYPASS:
        metal='Dy' if g.endswith('_Dy') or '_Dy__' in g else 'La';mt=500 if metal=='La' else 501
        if g in METAL:
            row=copy.deepcopy(parent['rows'][g]);cell=selected['cells'][g];require(cell['status']=='complete','Saved final cell incomplete')
            row['baseline_final_cell']=cell;key,meta=ev.render(values,row,domain,resolver,exports=False);oldkey=verify(cell['program_inputs']['analyze']['key']).read_text()
            require(old.scientific_key(key)==old.scientific_key(oldkey),'Literal feca scientific key')
            require(all(verify(row['source_input'][k]).read_bytes()==verify(cell['program_inputs']['analyze'][k]).read_bytes() for k in ('xyz','parameters')),'Literal feca coordinates/parameters')
            row['input']=row['source_input'];row['physical']=meta['physical_values_a_b_c_d_R_epsilon']
        else:
            row=copy.deepcopy(qual['rows'][g]);row['source_input']=copy.deepcopy(row['input']);row['metal']=metal
            physical=ev.effective(values,metal,domain);key=old.clean_key(verify(row['input']['key']).read_text());key+='POLAR-LOCAL-RESPONSE '+str(mt)+' 36 '+' '.join(format(v,'.17g') for v in physical[:4])+'\n';row['physical']=physical
            for kw in ('polar-print','save-uinduce'):
                if not any(tok[0].lower()==kw for tok in old.tokens(key)):key+=kw+'\n'
            row['permanent_lab_Debye']=h['field'].permanent_moment(row)
        row.update(base_key=key,metal_type=mt,e=.1 if metal=='La' else -.1,physical_POLPAIR=domain['fixed_native_POLPAIR'][metal])
        require(not any(tok[0].lower() in ('polar-radial-drive','response-matrix-export','external-field') for tok in old.tokens(key)),'Clean base reporting/shape key')
        header,atoms,xyz=h['old_harness'].xyz(verify(row['input']['xyz']));require(len(atoms)==row['atoms'] and [int(a[5]) for a in atoms]==row['types'] and xyz==row['positions_A'],'Exact consumed coordinates/types')
        rows[g]=row
    jobs=inventory();require(sum(6*rows[j['cell']]['atoms'] for j in jobs if j['fd'])==132,'132 FD denominator')
    for job in jobs:
        row=rows[job['cell']];folder=out/'inputs'/job['id'];folder.mkdir(parents=True)
        for k,name in FILES:
            if k!='key':shutil.copyfile(verify(row['input'][k]),folder/name)
        (folder/'system.key').write_text(render(row,job,parent['units']))
        if 'q' in job:
            header,atoms,xyz=h['old_harness'].xyz(folder/'system.xyz');shift,tangent=path_geometry(row,job['kind'],job['q']);(folder/'system.xyz').write_text('\n'.join([header]+[' '.join(a[:2]+[format(v,'.17g') for v in x]+a[5:]) for a,x in zip(atoms,shift)])+'\n');job['tangent']=tangent;job['runtime_positions_A']=shift
        job['input']={k:pin(folder/name) for k,name in FILES}
    manifest={'status':'prepared_radial_qualification_no_calls','implementation':pin(__file__),'configuration':pin(config['configuration_path']),'sources':s,'radial_inputs':inp,'helper_sources':qual['sources'],'evaluator_config':cfg,'rows':rows,'programs':jobs,'units':parent['units'],'maximum_native_programs':49,'maximum_internal_FD_energies':132,'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1},'native_calls':0,'solver_calls':0}
    save(out/'MANIFEST.json',manifest);return manifest
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True);args=ap.parse_args();cfg=json.loads(Path(args.config).read_text());require(cfg['implementation']==pin(__file__),'Preparation source');cfg['configuration_path']=args.config;out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    try:r=prepare(cfg,out);print(json.dumps({'status':r['status'],'count':len(r['programs']),'cost':r['cost']}))
    except Exception as e:save(out/'FAILURE.json',{'error':repr(e),'implementation':pin(__file__),'native_calls':0});raise
