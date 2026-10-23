"""Saved Stage1 parity only: sparse actions and complex PME replay, no solver/native."""
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.sparse import coo_array,csc_array
from PREPARE import pin,verify,load,probes


class Rows:
    def __init__(self,path):self.stream=Path(path).open()
    def row(self,tag,n):
        fields=self.stream.readline().split()
        if len(fields)!=n+1 or fields[0]!=tag:raise ValueError('Malformed/missing '+tag)
        return fields[1:]
    def literal(self,text):
        if self.stream.readline().rstrip('\n')!=text:raise ValueError('Literal header '+text)
    def end(self):
        if self.stream.read().strip():raise ValueError('Unexpected trailing output')
        self.stream.close()


def floats(fields):
    array=np.asarray(list(map(float,fields)))
    if not np.isfinite(array).all():raise ValueError('Nonfinite output')
    return array


def compare(observed,predicted,label,checks):
    observed=np.asarray(observed);predicted=np.asarray(predicted)
    if observed.shape!=predicted.shape or not np.isfinite(observed).all() or not np.isfinite(predicted).all():
        raise ValueError('Shape/nonfinite '+label)
    error=float(np.max(np.abs(observed-predicted)))
    scale=max(1.,float(np.max(np.abs(predicted))));bound=1e-11*scale
    item=dict(label=label,max_abs_error=error,scale=scale,bound=bound,passed=error<=bound)
    checks.append(item)
    if not item['passed']:raise ValueError('Parity failed '+label)
    return item


def compact_response(sites,coefficients,gamma,pdamp):
    """Independent ordered scalar reconstruction; native damping values read back."""
    xyz=np.array([s['xyz'] for s in sites]);alpha=np.array([s['alpha'] for s in sites])
    oxygen=np.array([i for i,s in enumerate(sites) if s['type']==402]);d=xyz[oxygen]-xyz[0]
    # At exact half-cell ties compact support is zero; np.rint tie choice cannot affect h/H.
    d-=49.324*np.rint(d/49.324);r=np.linalg.norm(d,axis=1)
    weights=np.zeros(len(r));phi=np.zeros(len(r))
    factor=1/((1+5+25/3)*math.exp(-5))**2
    for j,rr in enumerate(r):
        if rr>=21:continue
        t=2*rr;p=1+t+t*t/3;weight=factor*(p*math.exp(-t))**2
        if rr>20:
            u=rr-20;weight*=(1-u)**3*(1+3*u+6*u*u)
        weights[j]=weight
        kappa=gamma[oxygen[j]]/(pdamp[0]*pdamp[oxygen[j]])**3;z=kappa*rr**3
        if z<1e-4:phi[j]=kappa*(1-z/2+z*z/6-z**3/24+z**4/120-z**5/720)
        else:phi[j]=-math.expm1(-z)/rr**3
    prefix=np.zeros(len(r)+1);suffix=np.zeros(len(r)+1);cross=0.
    for j,f in enumerate(weights):cross+=2*f*prefix[j];prefix[j+1]=prefix[j]+f
    for j in range(len(r)-1,-1,-1):suffix[j]=suffix[j+1]+weights[j]
    other=prefix[:-1]+suffix[1:];eta_i=other/(1+other)
    s=prefix[-1];eta=cross/(1+s*s);g=2*s/(1+s)
    a,b,c,dcoef=coefficients
    beta=np.zeros(len(sites));beta[0]=math.expm1(g*(c+dcoef*eta))/alpha[0]
    h=np.zeros((len(sites),3));h[oxygen]=(3*weights*(a+b*eta_i)*phi)[:,None]*d
    return h,beta


def read_current(path,snapshot,expected,coefficients):
    rows=Rows(path);rows.literal('PERIODIC_PARITY_STATE_V1')
    rows.literal('UNITS xyz=angstrom alpha=angstrom^3 beta=angstrom^-3 h=electron/angstrom^2')
    rows.literal('RHS LABELS raw_Ed raw_Ep; add_shared_h_to_both; no_solved_dipoles')
    grid=rows.row('GRID',6)
    if list(map(int,grid[:4]))!=[64,64,64,5] or list(map(float,grid[4:]))!=[snapshot['ewald'],snapshot['volume']]:
        raise ValueError('Current grid mismatch')
    rows.literal('MODE FFTW TINFOIL MUTUAL ALL_ACTIVE SIMPLE_LOOP')
    for axis in range(3):
        row=rows.row('TRANS',4)
        if int(row[0])!=axis+1 or not np.array_equal(floats(row[1:]),snapshot['trans'][axis]):raise ValueError('Transform')
    for metal in (900,901):
        row=rows.row('COEFFICIENT',5)
        if int(row[0])!=metal or not np.array_equal(floats(row[1:]),coefficients[str(metal)]):raise ValueError('Coefficients')
    n=len(expected);h=np.empty((n,3));beta=np.empty(n);direct=np.empty((n,6));pdamp=np.empty(n);thole=np.empty(n)
    gamma=np.zeros(n)
    for i,(site,raw) in enumerate(zip(snapshot['sites'],expected)):
        row=rows.row('SITE',9);values=floats(row[3:])
        if list(map(int,row[:3]))!=[i+1,site['type'],i+1] or site['atom']!=i+1:raise ValueError('Active map')
        if not np.array_equal(values[:4],raw['xyz']+[raw['alpha']]) or values[3]!=site['alpha'] or site['xyz']!=raw['xyz']:
            raise ValueError('Exact current alpha/coordinates')
        pdamp[i],thole[i]=values[4:]
        row=rows.row('RESPONSE',5)
        if int(row[0])!=i+1:raise ValueError('Response order')
        values=floats(row[1:]);beta[i]=values[0];h[i]=values[1:]
        row=rows.row('DIRECT',7)
        if int(row[0])!=i+1:raise ValueError('Direct order')
        direct[i]=floats(row[1:])
        if site['type']==402:
            row=rows.row('PAIR_GAMMA',2)
            if int(row[0])!=i+1:raise ValueError('Pair damping order')
            gamma[i]=floats(row[1:])[0]
        if list(map(int,rows.row('IGRID',4)))!=[i+1,*site['igrid']]:raise ValueError('Current grid map')
        for slot in range(5):
            for axis in range(3):
                row=rows.row('THETA',5)
                if list(map(int,row[:3]))!=[i+1,axis+1,slot+1] or not np.array_equal(floats(row[3:]),site['theta'][axis,slot]):
                    raise ValueError('Current native spline mismatch')
    for z in range(64):
        for y in range(64):
            for x in range(64):
                row=rows.row('QFAC',4)
                if list(map(int,row[:3]))!=[x+1,y+1,z+1] or floats(row[3:])[0]!=snapshot['qfac'][x,y,z]:
                    raise ValueError('Current native qfac mismatch')
    if list(map(int,rows.row('END',3)))!=[n,15*n,64**3]:raise ValueError('Current state denominator')
    rows.end()
    return h,beta,direct,pdamp,thole,gamma


def scatter_matrix(snapshot,reader):
    # At most 11958*125 stencil entries, each a three-component row. No dense A.
    count=len(snapshot['sites'])*snapshot['p']**3*3
    row=np.empty(count,dtype=np.int64);col=np.empty(count,dtype=np.int64);data=np.empty(count)
    k=0
    for grid,site,weights in reader.entries(snapshot):
        flat=np.ravel_multi_index(grid,snapshot['shape'])
        row[k:k+3]=flat;col[k:k+3]=np.arange(3*site,3*site+3);data[k:k+3]=weights;k+=3
    if k!=count:raise ValueError('Scatter denominator')
    return coo_array((data,(row,col)),shape=(math.prod(snapshot['shape']),3*len(snapshot['sites']))).tocsr()


def complex_reciprocal(B,snapshot,z):
    grid=np.asarray(B@z).reshape(snapshot['shape'])
    back=np.fft.ifftn(np.fft.fftn(grid)*snapshot['qfac'])*math.prod(snapshot['shape'])
    return np.asarray(B.T@back.ravel())


def real_action(lower,z):
    return np.asarray(lower@z+lower.T@z-lower.diagonal()*z)


def replay(manifest_path,folder,checks):
    tick,cpu=time.monotonic(),time.process_time();m=json.loads(Path(manifest_path).read_text());folder=Path(folder)
    static=json.loads(verify(m['static_input']).read_text());reader=load(m['snapshot_reader'],'parity_snapshot')
    for name,item in m['raw_inputs'].items():
        if pin(folder/name)['sha256']!=item['sha256']:raise ValueError('Raw input changed '+name)
    actual_snapshot=pin(folder/'pme_norm_export_001.dat')
    if actual_snapshot['sha256']!=static['pme_export']['sha256']:
        raise ValueError('PME state differs from admitted exact static geometry')
    with Path(actual_snapshot['path']).open() as stream:snapshot=reader.read(stream)
    sites_helper=load(m['input_site_helper'],'parity_expected_site_reader')
    expected,_=sites_helper.input_sites((folder/'system.xyz').read_text(),(folder/'amoebabio09.prm').read_text(),
                                      (folder/'system.key').read_text())
    h,beta,direct,pdamp,thole,gamma=read_current(folder/'parity_state.dat',snapshot,expected,m['coefficients'])
    # Alpha/Thole/POLPAIR are frozen original model parameters, independent of h/H.
    polar={};pair={}
    for name in ('amoebabio09.prm','system.key'):
        for line in (folder/name).read_text().splitlines():
            f=line.split()
            if f and f[0].lower()=='polarize':polar[int(f[1])]=(float(f[2]),float(f[3]))
            if f and f[0].lower()=='polpair':pair[tuple(sorted(map(int,f[1:3])))]=float(f[3])
    for i,s in enumerate(expected):
        a,t=polar[s['type']]
        if thole[i]!=t:raise ValueError('Literal native Thole mismatch')
        compare(np.array([pdamp[i]]),np.array([a**(1/6)]),'pdamp '+str(i+1),checks)
        if s['type']==402 and gamma[i]!=pair[(402,900)]:raise ValueError('Native POLPAIR differs')
    expected_h,expected_beta=compact_response(expected,m['coefficients']['900'],gamma,pdamp)
    compare(h,expected_h,'native compact h',checks);compare(beta,expected_beta,'native compact H',checks)
    compare(np.array([beta[0]*expected[0]['alpha']]),np.array([static['metal_stiffness']['normalized_beta']]),
            'normalized H versus job20053',checks)
    # Exact static lower CSC; no refactorization, inverse, or eigensolve.
    arrays={key:np.load(verify(item),mmap_mode='r',allow_pickle=False) for key,item in static['arrays'].items()}
    lower=csc_array((arrays['data'],arrays['indices'],arrays['indptr']),shape=tuple(static['shape']))
    B=scatter_matrix(snapshot,reader);ps=probes(m['probe_source'],static['shape'][0])
    rows=Rows(folder/'parity_apply.dat');rows.literal('PERIODIC_PARITY_APPLY_V1')
    rows.literal('CHANNELS real=uind imaginary=uinp NORMALIZED_M=I-C*(T-H)*C RECIPROCAL=positive_R')
    applications=[];n=len(expected)
    for mode in (1,2,3):
        for j in range(10):
            real=ps[j][1] if mode!=2 else np.zeros(3*n)
            imaginary=ps[j][1] if mode==2 else ps[(j+1)%10][1] if mode==3 else np.zeros(3*n)
            z=real+1j*imaginary;reciprocal=complex_reciprocal(B,snapshot,z);s=real_action(lower,z)
            full=np.empty(3*n,dtype=complex);recip=np.empty(3*n,dtype=complex)
            for i in range(n):
                row=rows.row('APPLY',15)
                if list(map(int,row[:3]))!=[mode,j+1,i+1]:raise ValueError('Action order')
                v=floats(row[3:]);full[3*i:3*i+3]=v[:3]+1j*v[3:6];recip[3*i:3*i+3]=v[6:9]+1j*v[9:12]
            own=[]
            for channel in ('real','imag'):
                for observed,predicted,label in ((recip,reciprocal,'reciprocal'),(full,s+reciprocal,'complete'),
                                                  (full-recip,s,'real+self+H')):
                    item=compare(getattr(observed,channel),getattr(predicted,channel),
                                 f'{mode}/{j+1}/{channel}/{label}',checks);own.append(item)
            applications.append(dict(mode=mode,probe=j+1,label=ps[j][0],checks=own))
    if rows.row('END',2)!=['30','30']:raise ValueError('Action denominator')
    rows.end()
    expected_receipt='\n'.join(['PERIODIC_PARITY_RECEIPT_V1',
        'COUNTS full=30 separate_reciprocal=30 total_reciprocal=60 direct_initializations=1',
        'COUNTS induced_solves=0 energies=0 gradients=0 factorizations=0',
        'NORMAL_GUARDS_UNCHANGED TRUE','STATE_RESTORED TRUE','END',''])
    if (folder/'parity_receipt.dat').read_text()!=expected_receipt:raise ValueError('Native terminal receipt')
    return dict(status='passed_finite_complete_periodic_operator_parity',implementation=pin(__file__),
                preparation=pin(manifest_path),candidate=m['candidate'],static_input=m['static_input'],
                applications=applications,checks=checks,counts=m['counts'],
                direct_channel_difference_max=float(np.max(np.abs(direct[:,:3]-direct[:,3:]))),
                normalized_beta=float(beta[0]*expected[0]['alpha']),
                actual_snapshot=actual_snapshot,raw={name:pin(folder/name) for name in
                    ('parity_state.dat','parity_apply.dat','parity_receipt.dat')},
                e_app=None,uniform_application_bound=False,energy_force_qualified=False,physical_accuracy_qualified=False,
                new_native_calls=0,new_factorizations=0,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))


if __name__=='__main__':
    import sys
    manifest,folder,out=sys.argv[1:]
    destination=Path(out);destination.mkdir(parents=True,exist_ok=False)
    checks=[]
    try:
        result=replay(manifest,folder,checks)
    except Exception as error:
        result=dict(status='failed_finite_complete_periodic_operator_parity',error=repr(error),
                    implementation=pin(__file__),preparation=pin(manifest),checks=checks,
                    new_native_calls=0,new_factorizations=0)
    (destination/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
    sys.exit(0 if result['status'].startswith('passed_') else 1)
