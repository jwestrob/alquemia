"""Analytic coupled physical-coordinate maps; no scientific energy engine."""
from __future__ import annotations
import numpy as np
from affordable_common import InvalidArtifact
from affordable_response import cap_jacobians


class Kinematics:
    def __init__(self, data):
        self.data=data; self.positions=np.array(data['positions_A']); self.core=np.array(data['core_positions_A'])
        self.modes=data['modes']; self.links=data['core_links']; self.metal=data['metal_index']; self.heavy=data['heavy_indices']

    @classmethod
    def from_site(cls, site):
        links=[[x.tolist() if isinstance(x,np.ndarray) else x for x in row] for row in site.core_links]
        return cls({'positions_A':site.positions.tolist(),'core_positions_A':site.core_positions.tolist(),
                    'modes':site.modes,'core_links':links,'metal_index':site.metal,
                    'heavy_indices':[i for i,a in enumerate(site.atoms) if a['element']!='H'],
                    'physical_ids':site.ids,'bonds':site.bonds})

    def evaluate(self,q):
        q=np.asarray(q,dtype=float); n=len(self.modes)
        if q.shape!=(n,) or not np.isfinite(q).all():raise InvalidArtifact('invalid physical coordinates')
        p=self.positions.copy(); jac=np.zeros((n,len(p),3))
        for k,(mode,angle) in enumerate(zip(self.modes,q)):
            ids=mode['moving_indices']
            if mode['unit']=='angstrom':
                p[ids]+=angle*np.array(mode['axis']);jac[k,ids]+=mode['axis'];continue
            i,j=mode['axis_indices'];a=p[i].copy(); da=jac[:,i].copy();v=p[j]-a;length=np.linalg.norm(v)
            if length<=1e-10:raise InvalidArtifact('zero rotation axis')
            axis=v/length;dv=jac[:,j]-da;dn=(dv-(dv@axis)[:,None]*axis)/length
            r=p[ids]-a;dr=jac[:,ids]-da[:,None,:];c,s=np.cos(angle),np.sin(angle)
            nr=r@axis;dnr=np.einsum('mi,ni->mn',dn,r)+np.einsum('i,mni->mn',axis,dr)
            value=a+c*r+s*np.cross(axis,r)+(1-c)*nr[:,None]*axis
            derivative=da[:,None]+c*dr+s*(np.cross(dn[:,None],r)+np.cross(axis,dr))
            derivative+=(1-c)*(dn[:,None]*nr[None,:,None]+axis*dnr[:,:,None])
            derivative[k]+=-s*r+c*np.cross(axis,r)+s*nr[:,None]*axis
            p[ids]=value;jac[:,ids]=derivative
        core=self.core.copy();core[0]=p[self.metal];cj=np.zeros((n,len(core),3));cj[:,0]=jac[:,self.metal]
        for item in self.links:
            i,kind,*args=item
            if kind=='source':core[i]=p[args[0]];cj[:,i]=jac[:,args[0]]
            elif kind=='cap':
                x,y,length,original=args;a,b=p[[x,y]];core[i]+=a+length*(b-a)/np.linalg.norm(b-a)-original
                ja,jb=cap_jacobians(a,b,length);cj[:,i]=jac[:,x]@ja.T+jac[:,y]@jb.T
            else:raise InvalidArtifact('unsupported physical core link')
        return p,core,jac,cj

    def positions_only(self,q):
        """Same rotation algebra without constructing unused Jacobian arrays."""
        q=np.asarray(q,dtype=float)
        if q.shape!=(len(self.modes),) or not np.isfinite(q).all():raise InvalidArtifact('invalid physical coordinates')
        p=self.positions.copy()
        for mode,angle in zip(self.modes,q):
            ids=mode['moving_indices']
            if mode['unit']=='angstrom':p[ids]+=angle*np.array(mode['axis']);continue
            i,j=mode['axis_indices'];a=p[i].copy();v=p[j]-a;length=np.linalg.norm(v)
            if length<=1e-10:raise InvalidArtifact('zero rotation axis')
            axis=v/length;r=p[ids]-a;c,s=np.cos(angle),np.sin(angle);nr=r@axis
            p[ids]=a+c*r+s*np.cross(axis,r)+(1-c)*nr[:,None]*axis
        return p

    def displacement(self,q):return float(np.linalg.norm(self.positions_only(q)[self.heavy]-self.positions[self.heavy],axis=1).max())

    def check(self,q):
        p,c,j,cj=self.evaluate(q);errors=[]
        for k in range(len(q)):
            d=np.zeros(len(q));d[k]=1e-6;plus=self.evaluate(q+d);minus=self.evaluate(q-d)
            errors.extend([float(np.max(np.abs((plus[0]-minus[0])/2e-6-j[k]))),float(np.max(np.abs((plus[1]-minus[1])/2e-6-cj[k])))])
        bonds=self.data['bonds'];err=max(abs(np.linalg.norm(p[i]-p[j])-np.linalg.norm(self.positions[i]-self.positions[j])) for i,j in bonds)
        moved=set(i for mode in self.modes for i in mode['moving_indices']);fixed=[i for i in range(len(p)) if i not in moved]
        return {'jacobian_max_error_A_per_unit':max(errors),'bond_max_error_A':float(err),
                'fixed_atoms_unchanged':bool(np.array_equal(p[fixed],self.positions[fixed])),
                'pass':bool(max(errors)<=1e-7 and err<=1e-9 and np.array_equal(p[fixed],self.positions[fixed]))}

    def direction(self,gradient):
        g=np.asarray(gradient);j=self.evaluate(np.zeros(len(self.modes)))[2][:,self.heavy]
        gram=np.einsum('mij,nij->mn',j,j);eig=np.linalg.eigvalsh(gram)
        if eig[0]/eig[-1]<=1e-10:raise InvalidArtifact('unsupported singular physical displacement metric')
        v=-np.linalg.solve(gram,g);speed=float(np.linalg.norm(np.einsum('m,mij->ij',v,j),axis=1).max())
        if speed<=1e-12:return {'status':'zero_projected_gradient','q':np.zeros(len(g)).tolist(),'metric_eigenvalues':eig.tolist()}
        scale=min(.20/speed,.20/max(np.max(np.abs(v[3:])),1e-300));v*=scale
        # Geometry-only shrink precedes every scientific path evaluation.
        low,high=0.,1.
        if max(self.displacement(f*v) for f in (.25,.5,.75,1.))>.20:
            for _ in range(50):
                mid=(low+high)/2
                if max(self.displacement(f*mid*v) for f in (.25,.5,.75,1.))<=.20:low=mid
                else:high=mid
            v*=low
        return {'status':'prepared','q':v.tolist(),'metric_eigenvalues':eig.tolist(),'initial_directional_derivative':float(g@v),
                'maximum_heavy_displacement_A':max(self.displacement(f*v) for f in (.25,.5,.75,1.))}
