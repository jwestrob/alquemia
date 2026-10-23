"""One finite manufactured qualification. No native, optimizer, or cached replay."""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    if os.environ.get(_name) != '1':
        raise RuntimeError(_name+' must be explicitly 1')
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform
import resource
import time
import traceback
import torch
import KERNEL as k

torch.set_num_threads(1)
torch.set_num_interop_threads(1)
D = torch.float64
K = 332.0637133
UNITS = {'electric': K, 'bohr_A': 0.529177210544, 'hartree_kcal_mol': 627.509474063}


def fixture():
    xyz = torch.tensor([
        [0,0,0], [2.5,0,0], [3.05,.73,0], [3.05,-.73,0],
        [0,2.4,.4], [.8,2.8,-.4], [.5,3.3,.2], [.5,4.7,.2],
        [1.2,5.1,.8], [-.2,5.1,.8], [.5,5.1,-.8],
        [-2.4,0,.2], [-3.5,.1,.3], [-4.1,1.2,.4],
        [-3.9,-.8,.2], [-5.,1.3,.4], [-3.6,2.0,.5]], dtype=D)
    z = (57,8,1,1,8,8,6,6,1,1,1,8,6,7,1,1,1)
    ligands = (
        k.Ligand('W','water',(1,2,3),(1,),{1:(2,3)}),
        k.Ligand('A','acetate',(4,5,6,7,8,9,10),(4,5),{4:(6,),5:(6,)}),
        k.Ligand('N','amide',(11,12,13,14,15,16),(11,),{11:(12,13)}))
    return xyz, k.Geometry(0,z,ligands)


def parameters():
    pairs = {kind:k.Pair(v0=.31*(i+1),vr=-.17,va=.13,hr=.007,hf=-.004,krr=.031,kff=-.009)
             for i,kind in enumerate(('water','acetate','amide'))}
    triples = {tuple(sorted((a,b))):k.Triple(.19,-.12,.08,.003,-.011)
               for a in pairs for b in pairs}
    return k.Parameters(pairs,triples)


def provider(geometry):
    def build(xyz, cell=None):
        # Fully declared, translation/rotation covariant manufactured baseline.
        anchor = xyz[geometry.metal_atom] if geometry.metal_atom is not None else xyz[0]
        R = xyz-anchor
        n = len(xyz); dim=3*n; r2=(R*R).sum(dim=1)
        diag=(.7+.03*r2).repeat_interleave(3)
        q=R.reshape(-1)
        A=torch.diag(diag)+.02*torch.outer(q,q)/(1+torch.dot(q,q))
        weights=torch.arange(1,n+1,dtype=D)/n
        ed=(.01*R*weights[:,None]).reshape(-1)
        ep=(.006*R*(1.1-weights[:,None])).reshape(-1)
        other=.13*r2.sum()
        if cell is not None:
            other=other+.001*torch.linalg.det(cell)
            A=A+torch.eye(dim,dtype=D)*(.0001*(cell*cell).sum())
        mu=(R*(.03*weights[:,None])).sum(dim=0)
        return k.Baseline(A,ed,ep,tuple(range(n)),K,other,mu,True,'manufactured_complete_provider_v1')
    return build


def close(a,b,tol=2e-10):
    a,b=torch.as_tensor(a,dtype=D),torch.as_tensor(b,dtype=D)
    err=float(torch.max(torch.abs(a-b))) if a.numel() else 0.
    scale=max(1.,float(torch.max(torch.abs(a))) if a.numel() else 0.,float(torch.max(torch.abs(b))) if b.numel() else 0.)
    if err>tol*scale: raise AssertionError((err, tol*scale))
    return {'max_absolute_error':err,'allowed':tol*scale}


def rejects(call, text):
    try: call()
    except ValueError as error:
        if text not in str(error): raise
        return {'rejection':str(error)}
    raise AssertionError('expected rejection: '+text)


def permute(xyz,g,order):
    inverse={old:new for new,old in enumerate(order)}
    ligands=[]
    for l in g.ligands:
        ligands.append(k.Ligand(l.id,l.kind,tuple(inverse[a] for a in l.atom_indices),
            tuple(inverse[a] for a in l.donor_atoms),
            {inverse[a]:tuple(inverse[b] for b in bs) for a,bs in l.frame_atoms_by_donor.items()}))
    return xyz[order],k.Geometry(None if g.metal_atom is None else inverse[g.metal_atom],
                                tuple(g.atomic_numbers[a] for a in order),tuple(ligands))


def run_tests():
    rows=[]
    def test(name,call):
        try: rows.append({'name':name,'passed':True,'details':call()})
        except Exception as exc: rows.append({'name':name,'passed':False,'error':repr(exc),'traceback':traceback.format_exc()})
    x,g=fixture(); p=parameters(); sites=tuple(range(len(x))); base=provider(g)(x)
    def ev(y=x, gg=g, pp=p): return k.evaluate_provider(y,gg,pp,provider(gg))
    ref=ev()

    def zero():
        out=k.evaluate(x,g,k.Parameters(),base)
        d=out['delta']; assert float(d['V'])==0 and torch.count_nonzero(d['h'])==0 and torch.count_nonzero(d['H'])==0
        u=torch.linalg.solve(base.A0,base.Ed0); v=torch.linalg.solve(base.A0,base.Ep0)
        return {'energy':close(out['energy'],base.U_other0-.5*K*torch.dot(u,base.Ep0)),
                'u':close(out['u'],u),'v':close(out['v'],v)}
    test('zero_coefficients_and_baseline',zero)
    test('translation_invariance',lambda:close(ev(x+torch.tensor([.7,-1.1,2.3],dtype=D))['energy'],ref['energy']))

    def rotation():
        a=.713; Q=torch.tensor([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]],dtype=D)
        out=ev(x@Q.T); T=torch.kron(torch.eye(len(x),dtype=D),Q)
        return {'energy':close(out['energy'],ref['energy']), 'h':close(out['delta']['h'],T@ref['delta']['h']),
                'H':close(out['delta']['H'],T@ref['delta']['H']@T.T),
                'mu':close(out['mu_eA'],Q@ref['mu_eA']), 'chi':close(out['chi_A3'],Q@ref['chi_A3']@Q.T)}
    import math
    test('proper_rotation_covariance',rotation)
    def atom_swap(a,b):
        order=list(range(len(x)));order[a],order[b]=order[b],order[a]
        xx,gg=permute(x,g,order)
        # Compare correction only: the manufactured baseline's label-dependent
        # coefficients are not a chemical symmetry model.
        out=k.correction(xx,gg,sites,p); before=ref['delta']
        perm=[3*i+j for i in order for j in range(3)]
        return {'V':close(out['V'],before['V']),'h':close(out['h'],before['h'][perm]),
                'H':close(out['H'],before['H'][perm][:,perm])}
    test('equivalent_water_H_permutation',lambda:atom_swap(2,3))
    test('equivalent_acetate_O_permutation',lambda:atom_swap(4,5))
    test('ligand_order_invariance',lambda:close(k.correction(x,replace(g,ligands=tuple(reversed(g.ligands))),sites,p)['V'],ref['delta']['V']))

    def chelate():
        gg=replace(g,ligands=(g.ligands[1],))
        out=k.correction(x,gg,sites,k.Parameters(triples=p.triples))
        assert out['counts']['intact_triples']==0 and float(out['V'])==0 and torch.count_nonzero(out['H'])==0
        return out['counts']
    test('chelate_not_split_into_two_ligands',chelate)

    def frame_energy():
        two=replace(g,ligands=(g.ligands[0],g.ligands[2]))
        pp=k.Parameters(triples={('amide','water'):k.Triple(vf=1.)})
        xx=x.clone(); xx[2:4]=2*x[1]-x[2:4]
        a=k.correction(x,two,sites,pp); b=k.correction(xx,two,sites,pp)
        assert abs(float(a['V']-b['V']))>.01 and torch.count_nonzero(a['h'])==0 and torch.count_nonzero(a['H'])==0
        pair=k.Parameters(pairs={'water':k.Pair(va=1.)})
        pa=k.correction(x,g,sites,pair)['V'];pb=k.correction(xx,g,sites,pair)['V']
        assert abs(float(pa-pb))>.01
        return {'triple_frame_energy_change':float(b['V']-a['V']),'pair_directed_frame_change':float(pb-pa)}
    test('independent_directed_pair_and_triple_frame_energy',frame_energy)

    def separation():
        xx=x.clone()
        for ligand in g.ligands: xx[list(ligand.atom_indices)]+=torch.tensor([20.,20.,20.],dtype=D)
        xx.requires_grad_(); out=k.correction(xx,g,sites,p)
        loss=out['V']+out['h'].sum()+out['H'].sum(); grad=torch.autograd.grad(loss,xx,create_graph=True)[0]
        hess=torch.autograd.grad(grad.sum(),xx)[0]
        assert float(loss)==0 and torch.count_nonzero(grad)==0 and torch.count_nonzero(hess)==0
        return {'correction_gradient_hessian_exact_zero':True}
    test('separation_exact_zero_through_second_derivative',separation)
    def cutoff():
        values=[]
        for r in (6.,6.0001,5.9999,5.99995):
            t=torch.tensor(r,dtype=D,requires_grad=True); v=k.compact_weight(t*t,6.)
            first=torch.autograd.grad(v,t,create_graph=True)[0];second=torch.autograd.grad(first,t)[0]
            values.append([r,float(v),float(first),float(second)])
        assert values[0][1:]==[0.,0.,0.] and values[1][1:]==[0.,0.,0.]
        assert abs(values[3][3])<.51*abs(values[2][3])
        return {'values':values}
    test('C2_cutoff_boundary',cutoff)
    def contact():
        xx=x.clone();xx[1]=xx[0];xx.requires_grad_()
        out=k.correction(xx,g,sites,p); loss=out['V']+out['h'].sum()+out['H'].sum()
        grad=torch.autograd.grad(loss,xx,create_graph=True)[0];hess=torch.autograd.grad(grad.sum(),xx)[0]
        assert bool(torch.isfinite(grad).all()) and bool(torch.isfinite(hess).all())
        return {'finite_first_and_second_derivatives':True}
    test('contact_regularization',contact)

    def fragment(metal,ligands):
        gg=replace(g,metal_atom=metal,ligands=ligands)
        out=k.evaluate(x,gg,p,base); d=out['delta']
        assert float(d['V'])==0 and torch.count_nonzero(d['h'])==0 and torch.count_nonzero(d['H'])==0
        return close(out['energy'],k.evaluate(x,gg,k.Parameters(),base)['energy'])
    test('ligand_only_no_invented_metal',lambda:fragment(None,g.ligands))
    test('ion_only_empty_ligand_list',lambda:fragment(0,()))

    def dual():
        out=ref;u,v=out['u'],out['v'];delta=out['delta'];G=out['G'];ed=base.Ed0+delta['h'];ep=base.Ep0+delta['h']
        assert float(torch.max(torch.abs(u-v)))>.001
        L=base.U_other0+delta['V']+.5*K*(u@G@v-u@ep-v@ed)
        return {'stationary_energy':close(L,out['energy']), 'u_equation':close(G@u,ed),'v_equation':close(G@v,ep)}
    test('unequal_dual_stationary_energy',dual)

    def field_derivatives():
        F=torch.tensor([.0003,-.0002,.0001],dtype=D,requires_grad=True)
        out=k.evaluate(x,g,p,base,F)
        mu=-torch.autograd.grad(out['energy'],F,create_graph=True)[0]/K
        chi=torch.stack([torch.autograd.grad(mu[i],F,retain_graph=True)[0] for i in range(3)])
        return {'mu':close(mu,out['mu_eA']),'chi':close(chi,out['chi_A3'])}
    test('field_autodiff_dipole_susceptibility',field_derivatives)
    def field_fd():
        rows=[]
        for step in (1e-5,5e-6):
            f=torch.tensor([step,0.,0.],dtype=D)
            plus=k.evaluate(x,g,p,base,f);minus=k.evaluate(x,g,p,base,-f)
            rows.append({'step':step,'mu':close(-(plus['energy']-minus['energy'])/(2*step*K),ref['mu_eA'][0],1e-6),
                         'chi':close((plus['mu_eA']-minus['mu_eA'])/(2*step),ref['chi_A3'][:,0],1e-6)})
        return rows
    test('field_finite_difference_two_steps',field_fd)

    def coordinate():
        xx=x.clone().requires_grad_();out=ev(xx)
        grad=torch.autograd.grad(out['energy'],xx)[0]
        direction=torch.sin(torch.arange(x.numel(),dtype=D)+.3).reshape(x.shape);direction/=torch.linalg.vector_norm(direction)
        target=(grad*direction).sum(); rows=[]
        for step in (1e-5,5e-6):
            fd=(ev(x+step*direction)['energy']-ev(x-step*direction)['energy'])/(2*step)
            rows.append({'step':step,**close(fd,target,1e-6)})
        assert out['coordinate_derivative_scope']=='total_by_provider_contract'
        return rows
    test('complete_provider_total_coordinate_derivative',coordinate)
    def strain():
        cell=torch.diag(torch.tensor([9.,10.,11.],dtype=D));S=torch.tensor([[.2,.1,0.],[0.,-.3,.05],[0.,0.,.1]],dtype=D)
        def energy(t):
            F=torch.eye(3,dtype=D)+t*S
            return k.evaluate_provider(x@F.T,g,p,provider(g),cell=F@cell)['energy']
        t=torch.tensor(0.,dtype=D,requires_grad=True);value=energy(t);slope=torch.autograd.grad(value,t)[0]
        return [{'step':step,**close((energy(step)-energy(-step))/(2*step),slope,1e-6)} for step in (1e-5,5e-6)]
    test('complete_provider_coordinate_and_cell_strain_derivative',strain)
    def partial():
        xx=x.clone().requires_grad_();frozen=k.evaluate(xx,g,p,base)
        gf=torch.autograd.grad(frozen['energy'],xx)[0]
        yy=x.clone().requires_grad_();full=ev(yy);gt=torch.autograd.grad(full['energy'],yy)[0]
        assert frozen['coordinate_derivative_scope']=='partial_frozen_baseline_not_total_force'
        difference=float(torch.max(torch.abs(gf-gt))); assert difference>.1
        return {'missing_baseline_derivative_changes_gradient':difference,'partial_scope':frozen['coordinate_derivative_scope']}
    test('frozen_operator_gradient_explicitly_incomplete',partial)
    test('incomplete_provider_rejected',lambda:rejects(lambda:k.evaluate(x,g,p,replace(base,coordinate_derivatives_complete=False),baseline_mode='provider'),'total derivative provider'))
    test('nonSPD_full_operator_rejected',lambda:rejects(lambda:k.evaluate(x,g,p,replace(base,A0=-torch.eye(len(base.Ed0),dtype=D))),'not numerically SPD'))
    def symmetry():
        A=base.A0.clone();A[0,1]+=1e-13
        raw=replace(base,A0=A)
        rejected=rejects(lambda:k.evaluate(x,g,p,raw),'exactly symmetric')
        sym,receipt=k.symmetrize_native_operator(A,max_absolute_asymmetry=2e-13,admission_authority='manufactured inherited bound')
        assert torch.equal(sym,sym.T) and receipt['raw_max_absolute_asymmetry']>0
        rejects(lambda:k.symmetrize_native_operator(A,max_absolute_asymmetry=1e-14,admission_authority='manufactured strict bound'),'exceeds inherited')
        return {'implicit_rejected':rejected,'explicit_receipt':receipt}
    test('explicit_bounded_native_symmetry_policy',symmetry)
    def scopes():
        out=k.evaluate(x,g,p,replace(base,U_other0=None,permanent_mu_eA=None))
        assert out['energy_scope']=='polarization_plus_local_correction' and out['moment_scope']=='induced_only'
        return {'energy_scope':out['energy_scope'],'moment_scope':out['moment_scope']}
    test('missing_baseline_channels_not_fabricated',scopes)
    def units():
        F=torch.tensor([.001,0,0],dtype=D); native=k.field_au_to_native(F,UNITS)
        converted=k.response_to_source_au(ref['mu_eA'],ref['chi_A3'],UNITS)
        return {'field':close(native,F*627.509474063/(332.0637133*.529177210544),0),
                'source_chi':close(converted['chi_source_au'],ref['chi_A3']*627.509474063/(332.0637133*.529177210544**2),0)}
    test('literal_native_unit_conversions',units)
    return rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    path=Path(args.out);path.mkdir(parents=True,exist_ok=False)
    start=time.monotonic(); cpu=time.process_time();result={'status':'started','scope':'manufactured_only','native_calls':0,'quantum_calls':0,'fits':0,'cached_operator_replays':0}
    try:
        result['checks']=run_tests();result['status']='passed' if all(r['passed'] for r in result['checks']) else 'failed'
    except Exception as exc:
        result.update(status='failed',error=repr(exc),traceback=traceback.format_exc())
    result['cost']={'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    result['runtime']={'torch':torch.__version__,'python':platform.python_version(),'torch_threads':torch.get_num_threads(),'interop_threads':torch.get_num_interop_threads(),'thread_environment':{n:os.environ[n] for n in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')}}
    result['sources']={name:{'path':str(Path(__file__).with_name(name).resolve()),'sha256':hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()} for name in ('KERNEL.py','TEST_KERNEL.py','SPEC.md')}
    (path/'RESULT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'passed':sum(r['passed'] for r in result.get('checks',[])),'total':len(result.get('checks',[])),'cost':result['cost']}))
    raise SystemExit(0 if result['status']=='passed' else 1)


if __name__=='__main__': main()
