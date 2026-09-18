"""Physical donor torsions/crankshafts and source-aware QM cap Jacobians."""
from __future__ import annotations
import copy
import numpy as np
from affordable_common import InvalidArtifact,verify,xyz
from affordable_response import cap_jacobians
from mace_mechanics import physical_id

POLICY='matched_hybrid_physical_donor_coordinates_v1'
CHI={'ASP':(('CA','CB'),('CB','CG')),'ASN':(('CA','CB'),('CB','CG')),
     'GLU':(('CA','CB'),('CB','CG'),('CG','CD')),'GLN':(('CA','CB'),('CB','CG'),('CG','CD'))}


def rotate(points,origin,axis,angle):
    v=np.asarray(points)-origin;axis=np.asarray(axis);axis=axis/np.linalg.norm(axis)
    return origin+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+np.outer(v@axis,axis)*(1-np.cos(angle))


class SiteCoordinates:
    """Coordinate kinematics only: no forcefield energy, spring or relaxation."""
    def __init__(self,repair,whole,core_xyz):
        from affordable_peptide import SourceGraph
        self.repair=repair;self.whole=whole
        self.atoms=copy.deepcopy(whole['physical_atoms']);self.ids=[a['id'] for a in self.atoms]
        if len(set(self.ids))!=len(self.ids):raise InvalidArtifact('duplicate physical atom identity')
        self.lookup={v:i for i,v in enumerate(self.ids)};self.positions=np.array([a['xyz_A'] for a in self.atoms]);self.metal=self.lookup['metal']
        self.core=xyz(core_xyz);self.core_positions=np.array([a[1:] for a in self.core]);self.graph=SourceGraph(verify(repair['source_structure']),verify(repair['topology_definition']))
        self.bonds=[];self.adjacency={i:set() for i in range(len(self.ids))}
        for b in whole['preparation_details']['bonds']:
            i,j=(self.lookup[b[k]] for k in ('atom_a_id','atom_b_id'))
            self.bonds.append((i,j));self.adjacency[i].add(j);self.adjacency[j].add(i)
        self.modes=[{'id':'metal_'+a,'kind':'metal_translation','unit':'angstrom','axis':v,'moving_indices':[self.metal]} for a,v in zip('xyz',np.eye(3).tolist())]
        ledger=sorted(repair['charge_ledger'],key=lambda e:(e['source']['chain'],e['source']['resnum'],e['source'].get('insertion_code',''),e['kind']))
        for entry in ledger:
            if entry['kind']=='water':continue
            source=entry['source'];res=self.graph.locate(source);prefix=f"{res.chain}/{res.resnum}/{source.get('insertion_code','')}"
            if entry['kind']=='sidechain':
                if res.canonical_resname not in CHI:raise InvalidArtifact('unsupported donor torsion chemistry: '+res.canonical_resname)
                for number,(an,bn) in enumerate(CHI[res.canonical_resname],1):
                    i,j=[self.index(self.graph.key(res.key,n)) for n in (an,bn)]
                    if j not in self.adjacency[i]:raise InvalidArtifact('actual side-chain bond absent')
                    seen=set();stack=[j]
                    while stack:
                        k=stack.pop()
                        if k in seen:continue
                        seen.add(k)
                        stack.extend(n for n in self.adjacency[k] if {n,k}!={i,j})
                    if i in seen or any(self.ids[k].rsplit('/',1)[0]!=prefix for k in seen):raise InvalidArtifact('ring/crosslink or nonlocal side-chain component')
                    self.modes.append({'id':prefix+'/chi'+str(number),'kind':'sidechain_torsion','unit':'radian','axis_indices':[i,j],
                                       'moving_indices':sorted(seen),'source':source,'bond_names':[an,bn]})
            elif entry['kind']=='backbone_carbonyl':
                nitrogen=self.graph.amide_next.get(res.key)
                if nitrogen is None or self.graph.residues[nitrogen[:2]].canonical_resname=='PRO':raise InvalidArtifact('unsupported/missing actual next amide')
                n=self.index(nitrogen);hs=[i for i in self.adjacency[n] if self.atoms[i]['element']=='H']
                if len(hs)!=1:raise InvalidArtifact('peptide crankshaft requires exactly one actual N-H')
                c,o=[self.index(self.graph.key(res.key,s)) for s in ('C','O')]
                if n not in self.adjacency[c]:raise InvalidArtifact('actual source peptide bond absent')
                anchors=[self.index(self.graph.key(k,'CA')) for k in (res.key,nitrogen[:2])]
                self.modes.append({'id':prefix+'/peptide_crankshaft','kind':'peptide_crankshaft','unit':'radian','axis_indices':anchors,
                                   'moving_indices':[c,o,n,*hs],'source':source,'amide_nitrogen':self.graph.meta[nitrogen]})
            else:raise InvalidArtifact('unsupported donor ledger kind')
        if len({m['id'] for m in self.modes})!=len(self.modes):raise InvalidArtifact('duplicate physical coordinate')
        self.mapping=repair['atom_graph']['source_to_qm'];self.core_links=[]
        if sorted(row['qm_index'] for row in self.mapping)!=list(range(1,len(self.core))):raise InvalidArtifact('QM source mapping is not bijective')
        for row in self.mapping:
            i=row['qm_index']
            if row['kind']=='source':
                p=self.lookup[physical_id(row['source'])];self.core_links.append((i,'source',p))
                if self.core[i][0]!=self.atoms[p]['element'] or not np.allclose(self.core_positions[i],self.positions[p],atol=5e-10,rtol=0):raise InvalidArtifact('normalized physical/core source mismatch')
            elif row['kind']=='sigma_link_H':
                x,y=[self.lookup[physical_id(row[k])] for k in ('retained','omitted')];length=row['length_A'];a,b=self.positions[[x,y]]
                expected=a+length*(b-a)/np.linalg.norm(b-a)
                if self.core[i][0]!='H' or not np.allclose(self.core_positions[i],expected,atol=5e-10,rtol=0):raise InvalidArtifact('cap/source anchors differ')
                self.core_links.append((i,'cap',x,y,length,expected))
            else:raise InvalidArtifact('unsupported QM source mapping')
        if not np.allclose(self.core_positions[0],self.positions[self.metal],atol=5e-10,rtol=0):raise InvalidArtifact('core/whole metal positions differ')
        self.full_jacobian=self.jacobian_at_center()
        self.core_jacobian=self.map_jacobian(self.full_jacobian)
        heavy=[i for i,a in enumerate(self.atoms) if a['element']!='H']
        for mode,jac in zip(self.modes,self.full_jacobian):
            speed=float(np.linalg.norm(jac[heavy],axis=1).max())
            if speed<=0:raise InvalidArtifact('physical coordinate has no heavy motion')
            mode.update(maximum_heavy_speed_A_per_unit=speed,linear_probe_amplitude=.02/speed,
                        moving_atom_ids=[self.ids[i] for i in mode['moving_indices']])

    def index(self,key):return self.lookup[physical_id(self.graph.meta[key])]

    def jacobian_at_center(self):
        result=np.zeros((len(self.modes),len(self.atoms),3))
        for k,mode in enumerate(self.modes):
            indices=mode['moving_indices']
            if mode['unit']=='angstrom':result[k,indices]=mode['axis']
            else:
                i,j=mode['axis_indices'];origin=self.positions[i];axis=self.positions[j]-origin;axis/=np.linalg.norm(axis)
                result[k,indices]=np.cross(axis,self.positions[indices]-origin)
        return result

    def map_jacobian(self,full):
        result=np.zeros((len(self.modes),len(self.core),3));result[:,0]=full[:,self.metal]
        for item in self.core_links:
            i,kind,*args=item
            if kind=='source':result[:,i]=full[:,args[0]]
            else:
                x,y,length,_=args;ja,jb=cap_jacobians(self.positions[x],self.positions[y],length)
                result[:,i]=full[:,x]@ja.T+full[:,y]@jb.T
        return result

    def full_positions(self,q):
        q=np.asarray(q,dtype=float)
        if q.shape!=(len(self.modes),) or not np.isfinite(q).all():raise InvalidArtifact('invalid physical coordinate vector')
        result=self.positions.copy()
        for mode,amount in zip(self.modes,q):
            if not amount:continue
            indices=mode['moving_indices']
            if mode['unit']=='angstrom':result[indices]+=amount*np.array(mode['axis'])
            else:
                i,j=mode['axis_indices'];result[indices]=rotate(result[indices],result[i],result[j]-result[i],amount)
        return result

    def core_from_full(self,positions):
        result=self.core_positions.copy();result[0]=positions[self.metal]
        for item in self.core_links:
            i,kind,*args=item
            if kind=='source':result[i]=positions[args[0]]
            else:
                x,y,length,original=args;a,b=positions[[x,y]]
                # Preserve the archived cap's <=5e-10 A serialization offset.
                # Both bond anchors contribute analytically; caps are not DOFs.
                result[i]+=a+length*(b-a)/np.linalg.norm(b-a)-original
        return result

    def source_gradient(self,core_gradient):
        g=np.asarray(core_gradient);result=np.zeros_like(self.positions);result[self.metal]=g[0]
        if g.shape!=self.core_positions.shape or not np.isfinite(g).all():raise InvalidArtifact('invalid actual core gradient')
        for item in self.core_links:
            i,kind,*args=item
            if kind=='source':result[args[0]]+=g[i]
            else:
                x,y,length,_=args;ja,jb=cap_jacobians(self.positions[x],self.positions[y],length)
                result[x]+=ja.T@g[i];result[y]+=jb.T@g[i]
        return result

    def checks(self):
        errors_full=[];errors_core=[];bond_errors=[];unmoved=[]
        base_lengths=np.array([np.linalg.norm(self.positions[i]-self.positions[j]) for i,j in self.bonds])
        for k,mode in enumerate(self.modes):
            q=np.zeros(len(self.modes));q[k]=1e-6;p=self.full_positions(q);n=self.full_positions(-q)
            errors_full.append(float(np.max(np.abs((p-n)/2e-6-self.full_jacobian[k]))))
            errors_core.append(float(np.max(np.abs((self.core_from_full(p)-self.core_from_full(n))/2e-6-self.core_jacobian[k]))))
            q[k]=mode['linear_probe_amplitude'];probe=self.full_positions(q)
            lengths=np.array([np.linalg.norm(probe[i]-probe[j]) for i,j in self.bonds]);bond_errors.append(float(np.max(np.abs(lengths-base_lengths))))
            keep=[i for i in range(len(self.ids)) if i not in mode['moving_indices']];unmoved.append(bool(np.array_equal(probe[keep],self.positions[keep])))
        passed=max(errors_full+errors_core)<=1e-7 and max(bond_errors)<=1e-9 and all(unmoved)
        return {'pass':passed,'full_J_max_absolute_error_A_per_unit':max(errors_full),'core_J_max_absolute_error_A_per_unit':max(errors_core),
                'bond_length_max_change_A':max(bond_errors),'exterior_and_waters_unchanged':all(unmoved),'per_mode_full_errors':errors_full,'per_mode_core_errors':errors_core,
                'new_scientific_energy_calls':0}
