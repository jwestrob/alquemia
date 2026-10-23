"""Strict V2 state reader and cell-aware radial reconstruction; no FFT/solve."""
import math
from pathlib import Path
import numpy as np
import IDENTITY as I

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

def _legacy_compact_response(sites,coefficients,gamma,pdamp,cell):
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

def compact_response(sites,coefficients,gamma,pdamp,cell,e):
    """Old arithmetic unchanged at exact e=0; radial term alters h only."""
    I.cell(list(cell));I.require(abs(I.finite(e))<=1.,'Radial domain')
    I.require(len(sites)>1 and sites[0]['type'] in (900,901),'Metal at first site')
    I.require(sum(s['type'] in (900,901) for s in sites)==1,'One metal')
    I.require(all(I.finite(s['alpha'])>0 for s in sites),'Positive alpha')
    for s in sites:
        I.require(len(s['xyz'])==3,'XYZ shape')
        for x in s['xyz']:I.finite(x)
    oxygen=[i for i,s in enumerate(sites) if s['type']==402]
    I.require(0<len(oxygen)<=3980,'Water count domain')
    I.require(len(gamma)==len(sites) and len(pdamp)==len(sites),'Damping shapes')
    I.require(I.finite(pdamp[0])>0,'Metal damping')
    for i in oxygen:
        I.require(I.finite(pdamp[i])>0 and I.finite(gamma[i])>0,'Water damping')
    I.require(len(coefficients)==4,'Original coefficient shape')
    for x in coefficients:I.finite(x)
    h,beta=_legacy_compact_response(sites,coefficients,gamma,pdamp,cell)
    I.require(np.isfinite(h).all() and np.isfinite(beta).all(),'Finite old response')
    if e==0.:return h,beta
    cell=np.asarray(cell);xyz=np.asarray([s['xyz'] for s in sites]);factor=1/((1+5+25/3)*math.exp(-5))**2
    for i in oxygen:
        R=xyz[i]-xyz[0];R-=cell*np.rint(R/cell);r=float(np.linalg.norm(R))
        if r>=21.:continue
        t=2*r;p=1+t+t*t/3;F=factor*(p*math.exp(-t))**2
        if r>20.:
            u=r-20.;F*=(1-u)**3*(1+3*u+6*u*u)
        kappa=gamma[i]/(pdamp[0]*pdamp[i])**3;z=kappa*r**3
        if z<1e-4:phi=kappa*(1-z/2+z*z/6-z**3/24+z**4/120-z**5/720)
        else:phi=-math.expm1(-z)/r**3
        psi=(r*r-6.25)/6.25
        h[i]+=3*e*F*phi*psi*R
    I.require(np.isfinite(h).all() and np.isfinite(beta).all(),'Finite reconstructed response')
    return h,beta

def arithmetic_join(actual,predicted,label):
    """Only inferred cell algebra uses 128 eps; literal/native joins remain exact."""
    actual=floats(np.asarray(actual).ravel());predicted=floats(np.asarray(predicted).ravel())
    I.require(actual.shape==predicted.shape,'Cell algebra shape')
    bound=128*np.finfo(float).eps*max(1.,float(np.max(np.abs(predicted))))
    I.require(float(np.max(np.abs(actual-predicted)))<=bound,'Cell algebra '+label)

def read_header(rows,expected):
    I.cell(expected['cell']);I.settings(expected['settings']);I.coefficients(expected['coefficients']);I.radial(expected['radial'])
    rows.literal(I.SCHEMA)
    rows.literal('UNITS xyz=angstrom alpha=angstrom^3 beta=angstrom^-3 h=electron/angstrom^2')
    rows.literal('RHS LABELS raw_Ed raw_Ep; add_shared_h_to_both; no_solved_dipoles')
    grid=rows.row('GRID',6);g=list(map(int,grid[:4]));ewald,volume=floats(grid[4:])
    I.require(g==[64,64,64,5] and ewald==.4,'Fixed initialized grid')
    rows.literal('MODE FFTW TINFOIL MUTUAL ALL_ACTIVE SIMPLE_LOOP')
    cell=floats(rows.row('CELL',3)).tolist();I.cell(cell)
    cuts=floats(rows.row('CUTOFFS',3)).tolist();splits=floats(rows.row('EWALD_SPLIT',2)).tolist()
    ind=rows.row('INDUCTION',2);eps=float(floats(ind[:1])[0]);iters=int(ind[1])
    settings={'ewaldcut_A':cuts[0],'vdwcut_A':cuts[1],'vdwtaper_A':cuts[2],
              'aeewald_inverse_A':splits[0],'apewald_inverse_A':splits[1],
              'grid':g[:3],'order':g[3],'boundary':'TINFOIL','poleps':eps,'politer':iters}
    I.settings(settings);I.require(settings==expected['settings'] and cell==expected['cell'],'Literal expected cell/settings')
    transform=[]
    for axis in range(3):
        row=rows.row('TRANS',4);I.require(int(row[0])==axis+1,'Transform order');transform.append(floats(row[1:]).tolist())
    arithmetic_join(volume,cell[0]*cell[1]*cell[2],'volume')
    arithmetic_join(transform,np.diag(np.asarray(g[:3])/cell),'orthogonal transform')
    coefficients={};radial={}
    for metal in (900,901):
        row=rows.row('COEFFICIENT',5);I.require(int(row[0])==metal,'Coefficient order');coefficients[str(metal)]=floats(row[1:]).tolist()
    for metal in (900,901):
        row=rows.row('RADIAL_COEFFICIENT',2);I.require(int(row[0])==metal,'Both ordered radial slots');radial[str(metal)]=float(floats(row[1:])[0])
    I.radial(radial);I.require(coefficients==expected['coefficients'] and radial==expected['radial'],'Exact model coefficients')
    return {'schema':I.SCHEMA,'cell':cell,'settings':settings,'coefficients':coefficients,'radial':radial,
            'ewald':float(ewald),'volume':float(volume),'trans':transform}

def _read_body(rows,snapshot,expected,grid):
    """Private grammar helper. Public reader fixes grid64/order5; tiny controls only."""
    n=len(expected);I.require(n>0 and len(snapshot['sites'])==n,'Site denominator')
    h=np.empty((n,3));beta=np.empty(n);direct=np.empty((n,6));pdamp=np.empty(n);thole=np.empty(n);gamma=np.zeros(n)
    for i,(site,raw) in enumerate(zip(snapshot['sites'],expected)):
        row=rows.row('SITE',9);values=floats(row[3:])
        I.require(list(map(int,row[:3]))==[i+1,site['type'],i+1] and site['atom']==i+1 and raw['type']==site['type'],'Exact active map/type')
        I.require(np.array_equal(values[:4],raw['xyz']+[raw['alpha']]) and values[3]==site['alpha'] and site['xyz']==raw['xyz'],'Exact alpha/coordinates')
        I.require(values[3]>0 and values[4]>0 and values[5]>=0,'Site damping domain')
        pdamp[i],thole[i]=values[4:]
        row=rows.row('RESPONSE',5);I.require(int(row[0])==i+1,'Response order');values=floats(row[1:]);beta[i]=values[0];h[i]=values[1:]
        row=rows.row('DIRECT',7);I.require(int(row[0])==i+1,'Direct order');direct[i]=floats(row[1:])
        if site['type']==402:
            row=rows.row('PAIR_GAMMA',2);I.require(int(row[0])==i+1,'Pair damping order');gamma[i]=floats(row[1:])[0];I.require(gamma[i]>0,'Positive pair damping')
        I.require(list(map(int,rows.row('IGRID',4)))==[i+1,*site['igrid']],'Current grid map')
        for slot in range(5):
            for axis in range(3):
                row=rows.row('THETA',5)
                I.require(list(map(int,row[:3]))==[i+1,axis+1,slot+1] and np.array_equal(floats(row[3:]),site['theta'][axis,slot]),'Native spline identity')
    for z in range(grid[2]):
        for y in range(grid[1]):
            for x in range(grid[0]):
                row=rows.row('QFAC',4)
                I.require(list(map(int,row[:3]))==[x+1,y+1,z+1] and floats(row[3:])[0]==snapshot['qfac'][x,y,z],'Native qfac identity')
    I.require(list(map(int,rows.row('END',4)))==[n,15*n,math.prod(grid),2],'State denominator')
    rows.end()
    return {'h':h,'beta':beta,'direct':direct,'pdamp':pdamp,'thole':thole,'gamma':gamma}

def read_current(path,snapshot,expected,model):
    """Strict V2 only. Snapshot arrays are supplied by separately admitted reader."""
    rows=Rows(path)
    try:
        header=read_header(rows,model)
        I.require(header['ewald']==snapshot['ewald'] and header['volume']==snapshot['volume'] and np.array_equal(header['trans'],snapshot['trans']),'Exact same PME initialization')
        arrays=_read_body(rows,snapshot,expected,(64,64,64));arrays['identity_metadata']=header
        return arrays
    finally:rows.stream.close()

def compare_legacy_zero(v1_path,v2_path,authority,model):
    """Read-only byte comparison; returns receipt, never a projected state file.

    Caller must separately admit the complete V2/V1 grammars. This helper only
    compares raw bodies under explicit same-settings, zero-only authority.
    """
    rows=Rows(v2_path)
    try:header=read_header(rows,model)
    finally:rows.stream.close()
    I.validate_legacy_authority(authority,header)
    I.require(str(Path(v1_path).resolve())==authority['proof_pins']['state']['path'],'Legacy raw state binding')
    old=Path(v1_path).read_bytes();new=Path(v2_path).read_bytes();lines=new.splitlines(keepends=True)
    I.require(old.startswith(b'PERIODIC_PARITY_STATE_V1\n'),'Explicit legacy header')
    # Fixed new records only. Old numeric rendering and unequal RHS bytes survive.
    tags=(b'CELL ',b'CUTOFFS ',b'EWALD_SPLIT ',b'INDUCTION ',b'RADIAL_COEFFICIENT ')
    removed=[line for line in lines if line.startswith(tags)];I.require(len(removed)==6,'Projection header counts')
    kept=[line for line in lines if not line.startswith(tags)];kept[0]=b'PERIODIC_PARITY_STATE_V1\n'
    end=kept[-1];I.require(end.startswith(b'END ') and len(end)==53 and end.endswith(b'\n') and int(end[-13:-1])==2,'Native fixed-width V2 END')
    kept[-1]=end[:-13]+b'\n'
    I.require(b''.join(kept)==old,'Legacy zero raw body parity')
    return {'status':'passed_same_settings_zero_legacy_comparison','legacy_state':authority['proof_pins']['state'],
            'settings':header['settings'],'expected_state_created':False,'certificate_admission':False}
