"""Literal sites/model/settings; no operator or molecular calls."""
import math
from COMMON import *

def sites(inputs):
    xyz=verify(inputs['xyz']).read_text().splitlines();key=verify(inputs['key']).read_text();prm=verify(inputs['parameters']).read_text()
    require(len(xyz)==11960 and int(xyz[0].split()[0])==11958,'literal atom count')
    cell=list(map(float,xyz[1].split()));require(len(cell)==6 and cell[3:]==[90.]*3 and min(cell[:3])>=44,'orthogonal declared cell')
    alpha={}
    for text in (prm,key):
        for line in text.splitlines():
            v=line.split()
            if v and v[0].lower()=='polarize':alpha[int(v[1])]=float(v[2])
    output=[]
    for i,line in enumerate(xyz[2:],1):
        v=line.split();t=int(v[5]);r=list(map(float,v[2:5]));a=alpha[t]
        require(int(v[0])==i and t in (900,402,403,405,414) and a>=1e-8 and all(math.isfinite(x) for x in [a,*r]),'native input site')
        output.append({'atom':i,'type':t,'xyz':r,'alpha':a})
    return output,cell

def physical_metadata(row,m):
    expected,cell=sites(row['inputs']);meta=read(m['physics_metadata'])
    require(meta['all_alpha_positive'] is True and meta['metal_site_index']==0 and meta['sites']==11958 and meta['water_count']==3980,'frozen active map')
    require(meta['mutual_scales']=={'mutual-'+s+'-scale':1. for s in ('11','12','13','14')},'all-unity mutual scaling')
    coefficients,radial=response_records([verify(row['inputs'][k]).read_text() for k in ('parameters','key')])
    scales={};pairs={};polar={}
    for key in ('parameters','key'):
        seen=set()
        for line in verify(row['inputs'][key]).read_text().splitlines():
            v=line.split('#',1)[0].split()
            if not v:continue
            tag=v[0].lower()
            if tag in ('polar-local-response','polar-radial-drive'):continue
            elif tag=='polpair':
                types=tuple(sorted(map(int,v[1:3])));require((tag,types) not in seen,'duplicate pair');seen.add((tag,types));pairs[types]=float(v[3])
            elif tag.startswith('mutual-') and tag.endswith('-scale'):
                require(tag not in seen,'duplicate scale');seen.add(tag);scales[tag]=float(v[1])
            elif tag=='polarize' and int(v[1]) in (402,403,405,414,900,901):
                ty=v[1];require((tag,ty) not in seen,'duplicate polarize');seen.add((tag,ty));polar[ty]={'alpha':float(v[2]),'thole':float(v[3])}
            elif tag=='polar-donor-response':raise ValueError('donor chemistry outside fixed panel')
    require(polar==meta['type_parameters'] and scales==meta['mutual_scales'],'exact intrinsic alpha/damping/mutual scaling')
    require(pairs=={tuple(x['types']):x['gamma'] for x in meta['pair_overrides']},'exact pair damping')
    require(coefficients==m['model']['coefficients'],'Both immutable original coefficient slots')
    require(radial==m['model']['radial'],'Both explicit immutable radial slots')
    I.coefficients(coefficients);I.radial(radial)
    require(cell[:3]==row['cell'],'Literal declared cell')
    settings=key_settings(verify(row['inputs']['key']).read_text())
    require(settings==m['model']['settings'],'Literal key settings')
    require(expected[0]['type']==900 and sum(x['type']==900 for x in expected)==1 and sum(x['type']==402 for x in expected)==3980,'La/oxygen active map')
    return expected,cell,meta,coefficients,radial,pairs

def key_settings(text):
    wanted={'ewald':[],'ewald-cutoff':[12.],'vdw-cutoff':[9.],'vdw-taper':[8.1],
            'ewald-alpha':[.4],'pewald-alpha':[.4],'pme-grid':[64.,64.,64.],
            'pme-order':[5.],'ppme-order':[5.],'fft-package':['FFTW'],
            'polarization':['MUTUAL'],'polar-iter':[100.],'polar-eps':[1e-8],
            'pme-norm-export':[1.]}
    seen={}
    forbidden={'ewald-boundary','external-field','dewald','solvate','vdw-correction',
               'response-matrix-export','polar-donor-response','a-axis','b-axis','c-axis','alpha','beta','gamma'}
    for line in text.splitlines():
        v=line.split('#',1)[0].split()
        if not v:continue
        tag=v[0].lower();require(tag not in forbidden,'Unsupported key setting '+tag)
        if tag not in wanted:continue
        require(tag not in seen,'Duplicate setting '+tag)
        value=[x.upper() for x in v[1:]] if tag in ('fft-package','polarization') else list(map(float,v[1:]))
        require(value==wanted[tag],'Exact key setting '+tag);seen[tag]=value
    require(set(seen)==set(wanted),'All native settings explicit in key')
    return dict(I.SETTINGS)

def response_records(texts):
    coefficients={};radial={}
    for text in texts:
        for line in text.splitlines():
            v=line.split('#',1)[0].split()
            if not v:continue
            tag=v[0].lower()
            if tag not in ('polar-local-response','polar-radial-drive'):continue
            size=7 if tag=='polar-local-response' else 4
            require(len(v)==size,'Exact response arity')
            metal,water=map(int,v[1:3]);require(metal in (900,901) and water==402,'Response roles')
            target=coefficients if tag=='polar-local-response' else radial
            require(str(metal) not in target,'Duplicate response across input files')
            value=list(map(float,v[3:]));require(all(math.isfinite(x) for x in value),'Finite response')
            target[str(metal)]=value if tag=='polar-local-response' else value[0]
    I.coefficients(coefficients);I.radial(radial)
    return coefficients,radial
