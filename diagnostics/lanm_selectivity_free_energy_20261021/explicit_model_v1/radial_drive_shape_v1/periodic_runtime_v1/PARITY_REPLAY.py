"""Same saved sparse/complex action arithmetic, fresh V2 state/assembly bindings."""
import ast
import math
import time
import numpy as np
from scipy.sparse import coo_array,csc_array
from COMMON import *

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

def probes(source_item,n):
    # Import no factorization API: compile only the identical frozen pure helper.
    tree=ast.parse(verify(source_item).read_text())
    node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='ten_probes')
    namespace={'np':np}
    exec(compile(ast.Module(body=[node],type_ignores=[]),source_item['path'],'exec'),namespace)
    return namespace['ten_probes'](n,(0,1,2))

def replay(cfg,m,cert,folder):
    tick,cpu=time.monotonic(),time.process_time();folder=Path(folder);checks=[]
    static=read(cert['assembly']);initializer=read(cert['initializer'])
    require(static['status']=='stage2_assembly_complete' and static['initializer']==cert['initializer'],'Fresh assembly provenance')
    require(static['cutoff_A']==12. and static['ewald_inverse_A']==.4,'Actual cutoff12 assembly')
    I.validate_identity(static['identity'],cert['identity'])
    reader=load(m['helpers']['snapshot'],'v2_parity_snapshot')
    row=m['rows']['central']
    for key,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoebabio09.prm')]:
        require(pin(folder/name)['sha256']==row['inputs'][key]['sha256'],'Raw input bytes')
    require(pin(folder/'expected_state.dat')['sha256']==cert['expected_state']['sha256'],'Actual expected state')
    require(pin(folder/'parity_state.dat')['sha256']==cert['expected_state']['sha256'],'Exact same initialized V2 state')
    require(pin(folder/'probes.dat')['sha256']==m['probe_data']['sha256'],'Frozen probe input')
    actual_snapshot=pin(folder/'pme_norm_export_001.dat')
    require(actual_snapshot['sha256']==cert['pme_export']['sha256']==initializer['pme_export']['sha256'],'Same initialized PME')
    with Path(actual_snapshot['path']).open() as stream:snapshot=reader.read(stream)
    metadata=load(m['execution_sources']['METADATA.py'],'v2_parity_metadata')
    expected,cell=metadata.sites(row['inputs'])
    state=identity_state(m);model={'cell':cell[:3],**{k:m['model'][k] for k in ('settings','coefficients','radial')}}
    values=state.read_current(folder/'parity_state.dat',snapshot,expected,model)
    h,beta,direct,pdamp,thole,gamma=(values[k] for k in ('h','beta','direct','pdamp','thole','gamma'))
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
    expected_h,expected_beta=state.compact_response(expected,model['coefficients']['900'],gamma,pdamp,cell[:3],model['radial']['900'])
    compare(h,expected_h,'native compact h',checks);compare(beta,expected_beta,'native compact H',checks)
    compare(np.array([beta[0]*expected[0]['alpha']]),np.array([static['metal_stiffness']['normalized_beta']]),
            'normalized H versus fresh cutoff12 assembly',checks)
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
                manifest=cfg['manifest'],candidate=m['model']['candidate'],static_input=cert['assembly'],identity=cert['identity'],campaign_identity=policy.campaign_identity(m,cfg['build']),
                applications=applications,checks=checks,counts={'complete_operator':30,'separate_reciprocal':30,'total_reciprocal':60,'direct_initializations':1,'induced_solves':0,'energies':0,'gradients':0},
                direct_channel_difference_max=float(np.max(np.abs(direct[:,:3]-direct[:,3:]))),
                normalized_beta=float(beta[0]*expected[0]['alpha']),
                actual_snapshot=actual_snapshot,raw={name:pin(folder/name) for name in
                    ('parity_state.dat','parity_apply.dat','parity_receipt.dat')},
                e_app=None,uniform_application_bound=False,energy_force_qualified=False,physical_accuracy_qualified=False,
                new_native_calls=0,new_factorizations=0,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))

