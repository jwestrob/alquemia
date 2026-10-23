"""Original active-set direction solver; parameter dimension follows z."""
import itertools
import numpy as np

def require(ok,message):
    if not ok:raise ValueError(message)

def directions(J,r,z):
    full,_,rank,singular=np.linalg.lstsq(J,-r,rcond=1e-6)
    lower,upper=np.maximum(-.025,-z),np.minimum(.025,1-z);factor=1.
    for v,lo,hi in zip(full,lower,upper):
        if v>0:factor=min(factor,hi/v)
        elif v<0:factor=min(factor,lo/v)
    scaled=full*factor;bounded=None;loss=None
    for pattern in itertools.product((-1,0,1),repeat=len(z)):
        fixed=[j for j,v in enumerate(pattern) if v];free=[j for j,v in enumerate(pattern) if not v];d=np.zeros(len(z))
        for j in fixed:d[j]=lower[j] if pattern[j]<0 else upper[j]
        if free:d[free]=np.linalg.lstsq(J[:,free],-(r+J[:,fixed]@d[fixed]),rcond=1e-6)[0]
        if np.any(d<lower-1e-12) or np.any(d>upper+1e-12):continue
        d=np.clip(d,lower,upper);value=float(np.linalg.norm(r+J@d)**2)
        if loss is None or value<loss:bounded,loss=d,value
    require(bounded is not None,'No feasible bounded direction')
    return scaled,bounded,{'rank':int(rank),'singular_values':singular.tolist(),
        'unrestricted_GN_not_evaluated':full.tolist(),'scaled_GN':scaled.tolist(),'bounded_LS':bounded.tolist()}
