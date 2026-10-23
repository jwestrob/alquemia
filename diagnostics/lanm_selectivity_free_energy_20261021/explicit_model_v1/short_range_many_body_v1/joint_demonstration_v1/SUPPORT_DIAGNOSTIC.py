"""Saved-operator support bounds; callable only by the declared Slurm comparator.

No fitting, parameter sweep, new native call, or physical tolerance. U_H spans
all columns of the four actual stiffness features. U_mu also spans the three
actual drive features. Pure mechanical V never enters either support.
"""
import os,hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np


def orthogonalize(M):
    M=np.asarray(M,dtype=float)
    u,s,vh=np.linalg.svd(M,full_matrices=False)
    cutoff=np.finfo(float).eps*max(M.shape)*(float(s[0]) if len(s) else 0.)
    rank=int(np.count_nonzero(s>cutoff))
    return u[:,:rank],{'shape':list(M.shape),'singular_values':s.tolist(),'rank':rank,
                      'rank_cutoff':cutoff,'rule':'eps(float64)*max(shape)*largest singular value; algebraic numerical rank, not a physical accuracy gate'}


def report(states,bundle):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('support arithmetic belongs to the declared Slurm comparison')
    units=bundle['units'];D=units['debye_per_eA'];rows={}
    for mode in bundle['chi_targets']:
        id=mode['fixture_id'];s=states[id];f=bundle['operators'][id]
        A=s['A'];L=np.linalg.cholesky(A)
        # H matrices may contain redundant columns. Preserve all four actual
        # features and apply one explicitly reported ordinary SVD rank rule.
        UH,Hrank=orthogonalize(np.concatenate(list(s['H']),axis=1))
        Umu,murank=orthogonalize(np.column_stack((UH,s['h'])))
        b=s['B']@np.asarray(mode['direction']);z=np.linalg.solve(L,b)
        W=np.linalg.solve(L,UH);Q,whiterank=orthogonalize(W)
        perpendicular=z-Q@(Q.T@z)
        conv=units['hartree_kcal_mol']/(s['K']*units['bohr_A']**2)
        floor=float(perpendicular@perpendicular)*conv;chi0=float(z@z)*conv
        solved=np.linalg.solve(L.T,np.linalg.solve(L,np.column_stack((s['Ed'],s['Ep'],Umu))))
        baseline_mu=(s['permanent']+.5*s['B'].T@(solved[:,0]+solved[:,1]))*D
        image=s['B'].T@solved[:,2:]*D;P,image_rank=orthogonalize(image)
        donors=sorted({a for lig in f['geometry']['ligands'] for a in lig['donor_atoms']})
        sites=f['native_operator']['site_atom_indices']
        if any(f['geometry']['atomic_numbers'][a]!=8 for a in donors):raise ValueError('declared donor oxygen support')
        oi=[3*sites.index(a)+axis for a in donors for axis in range(3)]
        UO=np.eye(len(A))[:,oi];WO=np.linalg.solve(L,UO);QO,oxygen_rank=orthogonalize(WO)
        poza=z-QO@(QO.T@z);oxygen_floor=float(poza@poza)*conv
        ZO=np.linalg.solve(L.T,WO);PO,oxygen_image_rank=orthogonalize(s['B'].T@ZO*D)
        targets=[]
        for target in bundle['individual_moment_targets']:
            if target['model_id']!=id:continue
            mu=np.asarray(target['target_mu_e_bohr'])*units['bohr_A']*D
            delta=mu-baseline_mu;outside=delta-P@(P.T@delta);outside_O=delta-PO@(PO.T@delta)
            targets.append({'reference_id':target['id'],'method':target['method'],'source':target.get('source'),
                            'origin_A':target['origin_A'],'target_mu_Debye':mu.tolist(),'baseline_mu_Debye':baseline_mu.tolist(),
                            'target_minus_baseline_Debye':delta.tolist(),'perpendicular_component_Debye':outside.tolist(),
                            'irreducible_norm_for_current_support_Debye':float(np.linalg.norm(outside)),
                            'perpendicular_component_full_oxygen_support_Debye':outside_O.tolist(),
                            'irreducible_norm_full_oxygen_support_Debye':float(np.linalg.norm(outside_O)),
                            'qualification':'Static signed full-vector target; no claim that all three components have independent field-derivative qualification.'})
        g=f['geometry'];ligands=[]
        for lig in g['ligands']:
            ligands.append({'kind':lig['kind'],'donor_atoms':lig['donor_atoms'],
                            'frame_atoms_by_donor':lig['frame_atoms_by_donor']})
        rows[mode['id']]={'fixture_id':id,'fixture':bundle['fixture_pins'][id],
                         'atom_count':len(g['atomic_numbers']),'metal_atom':g['metal_atom'],'ligands':ligands,
                         'donor_count':sum(len(lig['donor_atoms']) for lig in g['ligands']),
                         'actual_H_feature_count':len(s['H']),'actual_h_feature_count':s['h'].shape[1],
                         'H_support':Hrank,'mu_support':murank,'whitened_H_support':whiterank,'moment_image':image_rank,
                         'direction':mode['direction'],'baseline_projected_chi_source_a0_cubed':chi0,
                         'optimistic_unbounded_stiffness_floor_source_a0_cubed':floor,
                         'reference_chi_source_a0_cubed':mode['target_a0_cubed'],
                         'floor_minus_reference_source_a0_cubed':floor-mode['target_a0_cubed'],
                         'reference_below_ordinary_double_floor':bool(mode['target_a0_cubed']<floor),
                         'full_oxygen_support':{'donor_atom_indices':donors,'response_component_indices':oi,'support_dimension':len(oi),'whitened_support':oxygen_rank,'moment_image':oxygen_image_rank,'optimistic_unbounded_chi_floor_source_a0_cubed':oxygen_floor,'floor_minus_reference_source_a0_cubed':oxygen_floor-mode['target_a0_cubed'],'reference_below_ordinary_double_floor':bool(mode['target_a0_cubed']<oxygen_floor),'interpretation':'Superset allowing every Cartesian response component on all donor oxygens; includes two O sites for acetate. No fitted coefficient or finite realization asserted.'},
                         'moment_targets':targets}
    return {'status':'completed_support_diagnostic','rows':rows,'new_solver_optimizations':0,'new_native_calls':0,'new_quantum_calls':0,
            'chi_derivation':'For A0=L L^T and U_H spanning H feature columns, floor=||(I-P_span(L^-1 U_H)) L^-1 b||^2. This equals b^T A0^-1 b-q^T(U_H^T A0^-1 U_H)^+q. Apply literal native/source unit conversion.',
            'chi_interpretation':'Optimistic lower limit allowing arbitrary PSD stiffness of unbounded size over U_H, a superset of the declared coefficient family. It is also an infimum bound for any stable supported symmetric update by block elimination. No finite coefficient attaining the limit is asserted.',
            'mu_derivation':'m=(u+v)/2; A0(m-m0)=h-Hm. Therefore delta_mu lies in range(B^T A0^-1 U_mu), U_mu=span(U_H,h features). Orthogonal target residual cannot be removed by any coefficients confined to this support, even unbounded.',
            'near_collinearity':'Actual radial and frame vectors may be nearly collinear from rounded source coordinates. Inspect all singular values/cutoffs. Full oxygen Cartesian support is a more permissive check independent of these near-collinear frame directions; do not infer strict impossibility solely from a suppressed near-null vector.',
            'scope':'Ordinary-double geometric support calculation with reported SVD ranks, not a rigorous interval certificate or a quantum mechanism assignment. Raw reference/method/printing limitations remain. All targets consumed; no model refit or new acceptance threshold.'}


def analyze(sources,baseline_predictions):
    """Called once in COMPARE. Construct 7 features for each of 5 modes only."""
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Slurm comparison required')
    def read_pin(p):
        path=Path(p['path']);data=path.read_bytes()
        if hashlib.sha256(data).hexdigest()!=p['sha256']:raise ValueError('changed source/input '+str(path))
        return data
    def load(name,p):
        read_pin(p);spec=importlib.util.spec_from_file_location(name,p['path']);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);return module
    import torch
    torch.set_num_threads(1)
    if torch.get_num_interop_threads()!=1:torch.set_num_interop_threads(1)
    C=load('support_admitted_core',sources['core']);K=load('support_frozen_kernel',sources['kernel']);D=load('support_admitted_data',sources['data'])
    bundle=D.load(json.loads(read_pin(sources['data_inputs'])),json.loads(read_pin(sources['endpoint_map'])))
    states={};symmetry={};feature_calls=0
    ids=list(dict.fromkeys(r['fixture_id'] for r in bundle['chi_targets']))
    if len(ids)!=5:raise ValueError('exact five measured susceptibility states required')
    for id in ids:
        f=bundle['operators'][id];op=f['native_operator'];q=f['native_qualification']
        if baseline_predictions[id]['fixture']!=bundle['fixture_pins'][id]:raise ValueError('saved baseline fixture identity')
        bound=min(r['numerical_bounds']['arithmetic'] for r in q['snapshots'])
        A,symmetry[id]=K.symmetrize_native_operator(C.T(op['A0_Aminus3']),max_absolute_asymmetry=bound,admission_authority={'fixture':bundle['fixture_pins'][id],'field':'native_qualification.snapshots'})
        xyz=C.T(f['geometry']['xyz_A']);g=C.geometry(K,f['geometry']);sites=tuple(op['site_atom_indices']);features={}
        for j in C.HIDX+C.hIDX:
            coeff=np.zeros(len(C.ORDER));coeff[j]=1.
            features[j]=K.correction(xyz,g,sites,C.parameters(K,coeff));feature_calls+=1
        states[id]={'A':A.numpy(),'Ed':np.asarray(op['Ed_e_per_A2']),'Ep':np.asarray(op['Ep_e_per_A2']),
                    'permanent':np.asarray(op['permanent_mu_eA']),'K':op['K_kcal_A_per_mol_e2'],
                    'B':np.tile(np.eye(3),(len(sites),1)),
                    'H':np.stack([features[j]['H'].numpy() for j in C.HIDX]),
                    'h':np.column_stack([features[j]['h'].numpy() for j in C.hIDX])}
    answer=report(states,bundle)
    answer.update(sources={k:sources[k] for k in ('core','kernel','data','data_inputs','endpoint_map')},
                  symmetry_projection=symmetry,feature_constructions=feature_calls,
                  baseline_operator_states=5,new_response_energy_evaluations=0,
                  saved_baseline_field_readback={id:{k:baseline_predictions[id][k] for k in ('mu_conjugate_Debye','chi_source_au','fixture')} for id in ids})
    return answer
