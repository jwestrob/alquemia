"""Strict initialized-PME snapshot reader and ordinary-double norm replay."""
import itertools
import math
import numpy as np


def read(stream):
    def row(tag, count):
        fields = stream.readline().split()
        if len(fields) != count+1 or fields[0] != tag:
            raise ValueError('Missing/malformed/out-of-order '+tag)
        return fields[1:]
    def floats(values):
        out = [float(x) for x in values]
        if not all(math.isfinite(x) for x in out):
            raise ValueError('Nonfinite snapshot')
        return out
    row('PME_NORM_EXPORT_V1', 0)
    sequence,n1,n2,n3,p,active,npole,M = map(int,row('GRID',8))
    shape=(n1,n2,n3)
    if sequence != 1 or min(shape)<p or not 2<=p<=8 or math.prod(shape)!=M or M>1048576:
        raise ValueError('Unsupported grid/sequence/stencil')
    if not 1<=active<=20000 or npole<active:
        raise ValueError('Unsupported active-site count')
    fft=row('FFT',1)[0]
    if fft not in ('FFTW','FFTPACK'):
        raise ValueError('Unknown FFT backend')
    chunks=list(map(int,row('CHUNKS',9)))
    nc,ng=chunks[:3],chunks[3:6]
    nl,nr,offset=chunks[6:]
    if (any(x<=0 for x in nc+ng) or tuple(x*y for x,y in zip(nc,ng))!=shape
            or nl+nr+1!=p or offset-nl!=2 or min(ng)<max(nl,nr)):
        raise ValueError('Unsupported chunk stencil')
    ewald,volume=floats(row('EWALD_VOLUME',2))
    if ewald<1e-6 or volume<=0:
        raise ValueError('Invalid initialized state')
    trans=[]
    for axis in range(1,4):
        rr=row('TRANS',4)
        if int(rr[0])!=axis: raise ValueError('Transform ordering')
        trans.append(floats(rr[1:]))
    q_stats=floats(row('Q_STATS',5)); native_norms=floats(row('NORMS',5))
    sites=[];seen=set()
    for _ in range(active):
        rr=row('SITE',9)
        atom,atype,*igrid=map(int,rr[:5]);alpha,*xyz=floats(rr[5:])
        if atom<=0 or atom in seen or alpha<1e-8:
            raise ValueError('Invalid/duplicate active site or PCG floor')
        seen.add(atom)
        theta=np.empty((3,p,2))
        for slot in range(1,p+1):
            for axis in range(1,4):
                tr=row('THETA',5)
                if list(map(int,tr[:3])) != [atom,axis,slot]:
                    raise ValueError('Incomplete/duplicate spline record')
                theta[axis-1,slot-1]=floats(tr[3:])
        sites.append(dict(atom=atom,type=atype,igrid=igrid,alpha=alpha,xyz=xyz,theta=theta))
    qfac=np.empty(shape)
    for z in range(n3):
        for y in range(n2):
            for x in range(n1):
                qr=row('QFAC',4)
                if list(map(int,qr[:3])) != [x+1,y+1,z+1]:
                    raise ValueError('Incomplete/duplicate multiplier record')
                qfac[x,y,z]=floats(qr[3:])[0]
    if list(map(int,row('END',4))) != [active,3*p*active,M,active*p**3]:
        raise ValueError('Incomplete snapshot counters')
    if stream.read().strip(): raise ValueError('Trailing snapshot data')
    return dict(shape=shape,p=p,sites=sites,trans=np.asarray(trans),qfac=qfac,
                native_norms=native_norms,q_stats=q_stats,sequence=sequence,
                fft=fft,chunks=chunks,ewald=ewald,volume=volume)


def entries(snapshot):
    shape,p=snapshot['shape'],snapshot['p']
    if p>min(shape): raise ValueError('Aliased stencil unsupported')
    A=np.asarray(snapshot['trans'])
    for site_index,site in enumerate(snapshot['sites']):
        theta=site['theta'];alpha=site['alpha']
        for j3 in range(p):
            for j2 in range(p):
                for j1 in range(p):
                    slots=(j1,j2,j3)
                    # Native 1-based grid = modulo(igrid+slot1,nfft)+1.
                    grid=tuple((site['igrid'][b]+slots[b]+1)%shape[b] for b in range(3))
                    t,u,v=[theta[b,slots[b]] for b in range(3)]
                    dv=np.array([t[1]*u[0]*v[0],t[0]*u[1]*v[0],t[0]*u[0]*v[1]])
                    yield grid,site_index,math.sqrt(alpha)*(dv@A)


def replay(snapshot):
    shape=snapshot['shape'];M=math.prod(shape)
    rows=np.zeros(shape);columns=np.zeros((len(snapshot['sites']),3));frob=0.
    for grid,i,b in entries(snapshot):
        rows[grid]+=np.abs(b).sum();columns[i]+=np.abs(b);frob+=float(b@b)
    c1=float(columns.max());cinf=float(rows.max())
    q=snapshot['qfac'];neg=q[np.ix_(*[(-np.arange(n))%n for n in shape])]
    qabs=float(np.abs(q).max());even=float(np.abs(q-neg).max())
    scale=M*min(c1*cinf,frob)
    result={'c1':c1,'cinf':cinf,'cF':frob,'U_R_nominal':scale*qabs,
            'evenness_perturbation_nominal':scale*even/2,
            'q_min':float(q.min()),'q_max':float(q.max()),'q_abs_max':qabs,
            'q_zero':float(q[0,0,0]),'q_conjugate_difference_max':even,
            'q_negative_count':int((q<0).sum()),'M':M,
            'arithmetic':'ordinary_double_not_outward','PSD_admitted':False,
            'full_operator_stability_admitted':False}
    if not all(math.isfinite(x) for x in result.values() if isinstance(x,float)):
        raise ValueError('Nonfinite replay')
    return result
