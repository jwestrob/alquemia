"""Stage1 pure readback, with explicit declared cell for compact response."""
import math
from pathlib import Path
import numpy as np

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

def compact_response(sites,coefficients,gamma,pdamp,cell):
    """Independent ordered scalar reconstruction; native damping values read back."""
    xyz=np.array([s['xyz'] for s in sites]);alpha=np.array([s['alpha'] for s in sites])
    oxygen=np.array([i for i,s in enumerate(sites) if s['type']==402]);d=xyz[oxygen]-xyz[0]
    # At exact half-cell ties compact support is zero; np.rint tie choice cannot affect h/H.
    cell=np.asarray(cell,dtype=float)
    if cell.shape!=(3,) or not np.isfinite(cell).all() or min(cell)<44:raise ValueError("Declared orthogonal cell")
    d-=cell*np.rint(d/cell);r=np.linalg.norm(d,axis=1)
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
