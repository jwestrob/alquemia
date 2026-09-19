"""Physical rigid-water coordinates and coupled curvature; no fitted springs."""
from __future__ import annotations
import numpy as np
from scipy.spatial.transform import Rotation
from affordable_common import InvalidArtifact

MASSES = {'H': 1.00794, 'O': 15.9994}  # atomic mass units, same in both endpoints


def left_jacobian(v):
    theta = float(np.linalg.norm(v)); x, y, z = v
    skew = np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
    if theta < 1e-4:
        a=.5-theta**2/24+theta**4/720; b=1/6-theta**2/120+theta**4/5040
    else:
        a=(1-np.cos(theta))/theta**2; b=(theta-np.sin(theta))/theta**3
    return np.eye(3)+a*skew+b*(skew@skew)


class WaterCoordinates:
    """COM translations in Å, lab-frame exponential rotations in radians."""
    def __init__(self, atoms, groups):
        self.initial=np.array([a[1:] for a in atoms],dtype=float)
        self.symbols=[a[0] for a in atoms]
        self.groups=[w for w in groups if w['role']=='variable']
        if not self.groups:raise InvalidArtifact('no moving physical water')
        self.centers=[];self.masses=[];used=set()
        for w in self.groups:
            ids=w['indices']
            if used.intersection(ids) or sorted(self.symbols[i] for i in ids)!=['H','H','O']:
                raise InvalidArtifact('overlap or incomplete physical water')
            used.update(ids)
            masses=np.array([MASSES[self.symbols[i]] for i in ids])
            self.masses.append(masses)
            self.centers.append(np.sum(self.initial[ids]*masses[:,None],axis=0)/masses.sum())
        self.centers=np.array(self.centers);self.dimension=6*len(self.groups)
        self.mobile=sorted(used);self.fixed=[i for i in range(len(atoms)) if i not in used]

    def positions(self,q):
        q=np.asarray(q).reshape(-1,6)
        if q.shape!=(len(self.groups),6) or not np.isfinite(q).all():raise InvalidArtifact('coordinate shape/nonfinite')
        coords=self.initial.copy()
        for w,c,v in zip(self.groups,self.centers,q):
            ids=w['indices'];coords[ids]=(self.initial[ids]-c)@Rotation.from_rotvec(v[3:]).as_matrix().T+c+v[:3]
        return coords

    def gradient(self,q,cartesian_gradient):
        q=np.asarray(q).reshape(-1,6);coords=self.positions(q);parts=[]
        for w,c,v in zip(self.groups,self.centers,q):
            ids=w['indices'];g=cartesian_gradient[ids]
            torque=np.sum(np.cross(coords[ids]-(c+v[:3]),g),axis=0)
            parts.append(np.r_[g.sum(axis=0),left_jacobian(v[3:]).T@torque])
        return np.concatenate(parts)

    def jacobian(self):
        j=np.zeros((len(self.initial),3,self.dimension))
        for k,(w,c) in enumerate(zip(self.groups,self.centers)):
            ids=w['indices'];j[ids,:,6*k:6*k+3]=np.eye(3)
            rel=self.initial[ids]-c
            for axis in range(3):j[ids,:,6*k+3+axis]=np.cross(np.eye(3)[axis],rel)
        return j

    def mass_matrix(self):
        g=np.zeros((self.dimension,self.dimension))
        for k,(w,c,m) in enumerate(zip(self.groups,self.centers,self.masses)):
            g[6*k:6*k+3,6*k:6*k+3]=m.sum()*np.eye(3)
            rel=self.initial[w['indices']]-c
            inertia=sum(mi*((v@v)*np.eye(3)-np.outer(v,v)) for mi,v in zip(m,rel))
            g[6*k+3:6*k+6,6*k+3:6*k+6]=inertia
        return g

    def rows(self,q):return [(z,*v) for z,v in zip(self.symbols,self.positions(q))]


def projected_hessian(gradient,dimension,step):
    columns=[]
    for i in range(dimension):
        d=np.zeros(dimension);d[i]=step
        columns.append((gradient(d)-gradient(-d))/(2*step))
    return np.stack(columns,axis=1)


def curvature_summary(coarse,fine,mass):
    norm=float(np.linalg.norm(fine,2))
    asym=float(np.linalg.norm(fine-fine.T,2));residual=float(np.linalg.norm(fine-coarse,2))
    k=(fine+fine.T)/2
    mv,mu=np.linalg.eigh(mass)
    if np.min(mv)<=0:raise InvalidArtifact('singular physical water mass metric')
    inverse=(mu/np.sqrt(mv))@mu.T
    eigenvalues,eigenvectors=np.linalg.eigh(inverse@k@inverse)
    return {'symmetric_curvature':k.tolist(),'mass_matrix':mass.tolist(),
            'mass_inverse_sqrt':inverse.tolist(),'mass_weighted_eigenvalues':eigenvalues.tolist(),
            'physical_eigenvectors':(inverse@eigenvectors).tolist(),
            'negative_mode_count':int(np.sum(eigenvalues<0)),
            'positive_definite':bool(np.min(eigenvalues)>0),
            'asymmetry_norm':asym,'refinement_norm':residual,'curvature_norm':norm,
            'symmetry_tolerance':max(.05,.001*norm),'refinement_tolerance':max(.05,.01*norm),
            'numerical_pass':bool(asym<=max(.05,.001*norm) and residual<=max(.05,.01*norm))}
