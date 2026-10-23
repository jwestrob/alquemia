"""Prospective 15-coefficient metal-energy comparison; frozen kernel and shared h/H."""
import math, time
import numpy as np
import torch
from scipy.optimize import minimize, lsq_linear

PARENT_CORE_SHA256='2aab54fe47c6c1f283fbdcda6d75e6303a441e0b2a3e2d2483ce4a18cd3c5fd0'
ORDER=[('water_La','v0','V'),('water','hf','h'),('water','krr','H'),('water','kff','H'),
       ('water_water_La','v0','V'),('water_water_La','va','V'),('water_water','h_neighbor','h'),('water_water','k_cross','H'),
       ('acetate','vr','V'),('acetate','va','V'),('acetate','hf','h'),('acetate','kff','H'),
       ('water_Dy','v0','V'),('water_water_Dy','v0','V'),('water_water_Dy','va','V')]
HIDX=[2,3,7,11];hIDX=[1,6,10];VIDX=[0,4,5,8,9,12,13,14]
SCALES=np.array([100. if s=='V' else .1 if s=='h' else 1. for _,_,s in ORDER])
DONOR='matched_dry_acetate_mono_central_La'
T=lambda x:torch.tensor(x,dtype=torch.float64)

def metal_atomic_number(g):
    atom=g['metal_atom']
    if atom is None:raise ValueError('operator fixture requires its explicit metal atom')
    z=g['atomic_numbers'][atom]
    if z not in (57,66):raise ValueError('only actual La57/Dy66 fixtures are declared')
    return z

def parameters(K,coeff,atomic_number):
    if len(coeff)!=15 or atomic_number not in (57,66):
        raise ValueError('exact15 coefficients and actual metal57/66 required')
    selected='La' if atomic_number==57 else 'Dy'
    pair={};triple={}
    for value,(kind,name,_) in zip(coeff,ORDER):
        if kind in ('water_La','water_Dy','water_water_La','water_water_Dy'):
            base,metal=kind.rsplit('_',1)
            if metal!=selected:continue
            kind=base
        if kind=='water_water':triple[name]=float(value)
        else:pair.setdefault(kind,{})[name]=float(value)
    return K.Parameters(pairs={k:K.Pair(**v) for k,v in pair.items()},triples={('water','water'):K.Triple(**triple)})

def geometry(K,g):
    return K.Geometry(g['metal_atom'],tuple(g['atomic_numbers']),tuple(K.Ligand.from_dict(l) for l in g['ligands']))

def prepare(bundle,K):
    if not bundle['complete_for_joint_demonstration'] or len(bundle['operators'])!=64:
        raise ValueError('exact 54 water + two donor centers + eight La endpoints required')
    states={};receipts={}
    for id,f in bundle['operators'].items():
        op=f['native_operator'];q=f['native_qualification'];raw=T(op['A0_Aminus3'])
        if 'snapshots' in q:
            bound=min(s['numerical_bounds']['arithmetic'] for s in q['snapshots'])
            authority={'fixture':bundle['fixture_pins'][id],'field':'native_qualification.snapshots'}
        else:
            admission=q['symmetry_projection'];bound=admission['inherited_max_absolute_asymmetry'];authority=admission['admission_authority']
        A,receipts[id]=K.symmetrize_native_operator(raw,max_absolute_asymmetry=bound,admission_authority=authority)
        xyz=T(f['geometry']['xyz_A']);g=geometry(K,f['geometry']);sites=tuple(op['site_atom_indices']);d=3*len(sites)
        features=[]
        for j in range(15):
            c=np.zeros(15);c[j]=1.;z=K.correction(xyz,g,sites,parameters(K,c,metal_atomic_number(f['geometry'])))
            features.append((float(z['V']),z['h'].numpy(),z['H'].numpy()))
        a=A.numpy();diag=np.diag(a)
        if not np.all(diag>0):raise ValueError('positive fixed baseline diagonal required')
        states[id]={'A':a,'Ed':np.array(op['Ed_e_per_A2']),'Ep':np.array(op['Ep_e_per_A2']),
                    'permanent':np.array(op['permanent_mu_eA']),'other':float(op['nonpolar_energy_kcal']),
                    'K':float(op['K_kcal_A_per_mol_e2']),'B':np.tile(np.eye(3),(len(sites),1)),
                    'D0':1/np.sqrt(diag),'H':np.stack([features[j][2] for j in HIDX]),
                    'h':np.column_stack([features[j][1] for j in hIDX]),'V':np.array([features[j][0] for j in VIDX])}
        assert states[id]['K']==bundle['units']['electric'] and a.shape==(d,d)
    return states,receipts

def operator(s,H):return s['A']+np.einsum('i,ijk->jk',H,s['H'])
def margins(states,H,floor):
    return {id:float(np.linalg.eigvalsh((s['D0'][:,None]*operator(s,H))*s['D0'][None,:])[0])-floor for id,s in states.items()}

class HProblem:
    """Exact affine operators and analytic chi/eigenvalue derivatives."""
    def __init__(self,states,bundle,settings):
        self.states=states;self.rows=bundle['chi_targets'];self.settings=settings
        self.conv=bundle['units']['hartree_kcal_mol']/(bundle['units']['electric']*bundle['units']['bohr_A']**2)
        self.last=None;self.cache=None;self.evaluations=0;self.singular_trial_objectives=0
    def geometry_cache(self,x):
        if self.last is None or not np.array_equal(x,self.last):
            operators={};constraint=[];jac=[]
            for id,s in self.states.items():
                G=operator(s,x);operators[id]=G;values,vectors=np.linalg.eigh((s['D0'][:,None]*G)*s['D0'][None,:])
                v=s['D0']*vectors[:,0];constraint.append(values[0]-self.settings['normalized_SPD_floor'])
                jac.append(np.einsum('j,ijk,k->i',v,s['H'],v))
            self.last=x.copy();self.cache=(operators,np.array(constraint),np.array(jac));self.evaluations+=1
        return self.cache
    def fun_jac(self,x):
        operators,_,_=self.geometry_cache(x);res=[];deriv=[]
        try:
            for row in self.rows:
                s=self.states[row['fixture_id']];b=s['B']@np.array(row['direction']);v=np.linalg.solve(operators[row['fixture_id']],b)
                res.append((b@v)*self.conv-row['target_a0_cubed'])
                deriv.append(-np.einsum('j,ijk,k->i',v,s['H'],v)*self.conv)
            r=np.array(res);J=np.array(deriv)
            if not np.isfinite(r).all() or not np.isfinite(J).all():raise np.linalg.LinAlgError('nonfinite off-domain objective')
            return float(np.mean(r*r)),2.*J.T@r/len(r)
        except np.linalg.LinAlgError:
            # No operator or data is clamped. Singular trial has no candidate
            # interpretation; the declared SPD constraints still govern endpoints.
            self.singular_trial_objectives+=1
            return 1e30,np.zeros(4)
    def constraint(self,x):return self.geometry_cache(x)[1]
    def constraint_jac(self,x):return self.geometry_cache(x)[2]

def electrical_base(states,H,units):
    """Conditional exact signed moment affine coefficients in normalized h."""
    ans={};D=units['debye_per_eA']
    for id,s in states.items():
        G=operator(s,H);ch=np.linalg.cholesky(G)
        rhs=np.column_stack((s['Ed'],s['Ep'],s['h']*.1,s['B']))
        q=np.linalg.solve(ch.T,np.linalg.solve(ch,rhs));u,v=q[:,0],q[:,1]
        ans[id]={'mu':(s['permanent']+.5*s['B'].T@(u+v))*D,
                 'mu_columns':s['B'].T@q[:,2:5]*D,
                 'chi_A3':s['B'].T@q[:,5:8]}
    return ans

def baseline_with_h(states,H,h,units):
    ans={};D=units['debye_per_eA']
    for id,s in states.items():
        G=operator(s,H);z=s['h']@h;ch=np.linalg.cholesky(G)
        uv=np.linalg.solve(ch.T,np.linalg.solve(ch,np.column_stack((s['Ed']+z,s['Ep']+z))))
        u,v=uv[:,0],uv[:,1]
        ans[id]={'energy':float(s['other']-.5*s['K']*u@(s['Ep']+z)),
                 'mu':(s['permanent']+.5*s['B'].T@(u+v))*D,'V_columns':s['V']*100.}
    return ans

def static_values(bundle):
    return {id:{'energy':float(s['energy_kcal_mol']),'mu':np.array(s['mu_Debye']),
                'mu_columns':np.zeros((3,3)),'V_columns':np.zeros(8)} for id,s in bundle['static'].items()}

def cycle_value(terms,pred,key):
    values=[c*np.asarray(pred[id][key]) for id,c in terms.items()]
    if values[0].ndim==0:return math.fsum(float(v) for v in values)
    return np.array([math.fsum(float(v.flat[j]) for v in values) for j in range(values[0].size)]).reshape(values[0].shape)

def fit_linear(families,n,settings):
    matrices=[];targets=[];rows={}
    for name,(matrix,target) in families.items():
        M=np.asarray(matrix).reshape((-1,n));y=np.asarray(target).reshape(-1)
        assert len(M)==len(y) and len(y)>0
        weight=1/math.sqrt(len(families)*len(y));matrices.append(M*weight);targets.append(y*weight)
        rows[name]={'scalar_rows':len(y),'row_weight':weight}
    A=np.concatenate(matrices);b=np.concatenate(targets);u,sv,vh=np.linalg.svd(A,full_matrices=False)
    cutoff=np.finfo(float).eps*max(A.shape)*(sv[0] if len(sv) else 0.);rank=int(np.count_nonzero(sv>cutoff))
    result=lsq_linear(A,b,bounds=(-np.ones(n),np.ones(n)),method='bvls',tol=settings['linear_tol'],max_iter=settings['linear_max_iter'])
    evidence={'status':int(result.status),'success':bool(result.success),'message':str(result.message),'iterations':int(result.nit),
              'normalized_solution':result.x.tolist(),'cost':float(result.cost),'optimality':float(result.optimality),
              'active_mask':result.active_mask.tolist(),'families':rows,'rank':rank,'columns':n,
              'singular_values':sv.tolist(),'rank_cutoff':float(cutoff),'null_directions_normalized':vh[rank:].tolist(),
              'rank_interpretation':'Numerical rank of this conditional development matrix; no full-tensor identifiability claim.'}
    if not result.success or not np.isfinite(result.x).all() or not np.all(np.abs(result.x)<=1):
        raise LinearStageFailure(evidence)
    return result.x,evidence

class LinearStageFailure(RuntimeError):
    def __init__(self,evidence):super().__init__('bounded linear stage did not complete');self.evidence=evidence

def h_fit(states,bundle,H,settings):
    pred={**static_values(bundle),**electrical_base(states,H,bundle['units'])};fam={}
    for name,cycles in [('pair_mu',bundle['pair_cycles']),('W9_mu3',[c for c in bundle['trimer_cycles'].values() if c['primary_fit']])]:
        M=[];y=[]
        for c in cycles:
            M.extend(cycle_value(c['model_terms'],pred,'mu_columns'))
            y.extend(np.array(c['stored_signed_mu_Debye'])-cycle_value(c['model_terms'],pred,'mu'))
        fam[name]=(M,y)
    target=np.array(bundle['donor']['moment_target_au'])*bundle['units']['bohr_A']*bundle['units']['debye_per_eA']
    fam['donor_mu']=(pred[DONOR]['mu_columns'],target-pred[DONOR]['mu'])
    x,e=fit_linear(fam,3,settings);return x*.1,e

def slope_ids(direction,tag):
    return (f'pbe0__La__{direction}__plus_{tag}__complex',f'pbe0__La__{direction}__minus_{tag}__complex')

def V_fit(states,bundle,H,h,settings):
    pred={**static_values(bundle),**baseline_with_h(states,H,h,bundle['units'])};fam={}
    for name,cycles in [('E2',bundle['pair_cycles']),('W9_E3',[c for c in bundle['trimer_cycles'].values() if c['primary_fit']])]:
        M=[];y=[]
        for c in cycles:
            M.append(cycle_value(c['model_terms'],pred,'V_columns'))
            y.append(c['energy_kcal_mol']-cycle_value(c['model_terms'],pred,'energy'))
        fam[name]=(M,y)
    M=[];y=[]
    for row in bundle['donor']['mechanical_modes']['qualified_two_step']:
        plus,minus=slope_ids(row['direction'],'005')
        M.append((pred[plus]['V_columns']-pred[minus]['V_columns'])/.01)
        y.append(row['reference_derivative']-(pred[plus]['energy']-pred[minus]['energy'])/.01)
    assert len(M)==2;fam['donor_slopes']=(M,y)
    x,e=fit_linear(fam,8,settings);return x*100.,e

def full_evaluate(states,bundle,K,coeff):
    """Final predictions use frozen KERNEL directly, with complete per-geometry caches."""
    pred=static_values(bundle);records={}
    for id,s in states.items():
        f=bundle['operators'][id];g=f['geometry'];op=f['native_operator'];params=parameters(K,coeff,metal_atomic_number(g));base=K.Baseline(T(s['A']),T(s['Ed']),T(s['Ep']),tuple(op['site_atom_indices']),s['K'],T(s['other']),T(s['permanent']),authority=bundle['fixture_pins'][id]['sha256'])
        out=K.evaluate(T(g['xyz_A']),geometry(K,g),params,base);converted=K.response_to_source_au(out['mu_eA'],out['chi_A3'],bundle['units'])
        G=out['G'].numpy();normalized=(s['D0'][:,None]*G)*s['D0'][None,:]
        mu=out['mu_eA'].numpy()*bundle['units']['debye_per_eA']
        pred[id]={'energy':float(out['energy']),'mu':mu,'chi_source_au':converted['chi_source_au'].numpy()}
        records[id]={'energy_kcal_mol':float(out['energy']),'mu_conjugate_Debye':mu.tolist(),
                     'mu_native_physical_u_Debye':((s['permanent']+s['B'].T@out['u'].numpy())*bundle['units']['debye_per_eA']).tolist(),
                     'chi_source_au':converted['chi_source_au'].tolist(),'chi_geometric_bohr3':converted['chi_geometric_bohr3'].tolist(),
                     'minimum_eigenvalue_Aminus3':float(out['eigenvalues'][0]),'minimum_normalized_eigenvalue':float(np.linalg.eigvalsh(normalized)[0]),
                     'solve_residual':float(out['solve_residual']),'DeltaV_kcal_mol':float(out['delta']['V']),
                     'polarization_energy_kcal_mol':float(out['polarization_energy']),'fixture':bundle['fixture_pins'][id],
                     'coordinate_derivative_scope':out['coordinate_derivative_scope']}
    return pred,records

def residual(value,target):
    v=np.asarray(value);t=np.asarray(target);r=v-t
    return {'prediction':v.tolist(),'target':t.tolist(),'error':r.tolist(),'sign_agreement':(np.sign(v)==np.sign(t)).tolist()}

def assess(pred,bundle):
    pairs={};trimers={};chi={};individual={};groups={k:[] for k in ('chi','pair_mu','W9_mu3','donor_mu','E2','E3','slope')}
    for name,cycles,store in [('pair',bundle['pair_cycles'],pairs),('trimer',list(bundle['trimer_cycles'].values()),trimers)]:
        for c in cycles:
            E=cycle_value(c['model_terms'],pred,'energy');mu=cycle_value(c['model_terms'],pred,'mu')
            er=residual(E,c['energy_kcal_mol']);mr=residual(mu,c['stored_signed_mu_Debye'])
            primary=name=='pair' or c['primary_fit'];store[c['id']]={'energy_kcal_mol':er,'mu_Debye':mr,'primary_fit':primary,'model_terms':c['model_terms'],'method':c['method']}
            if primary:
                groups['E2' if name=='pair' else 'E3'].append(er['error'])
                groups['pair_mu' if name=='pair' else 'W9_mu3'].extend(mr['error'])
    for row in bundle['chi_targets']:
        a=np.array(row['direction']);x=float(a@pred[row['fixture_id']]['chi_source_au']@a);r=residual(x,row['target_a0_cubed']);chi[row['id']]={**r,'direction':row['direction'],'units':'source a0^3','all_estimators':row['all_saved_estimators']};groups['chi'].append(r['error'])
    for row in bundle['individual_moment_targets']:
        target=np.array(row['target_mu_e_bohr'])*bundle['units']['bohr_A']*bundle['units']['debye_per_eA'];individual[row['id']]={**residual(pred[row['model_id']]['mu'],target),'units':'Debye','model_id':row['model_id'],'method':row['method'],'used_as_independent_full_state_fit_row':row['id']=='donor/La_acetate'}
    groups['donor_mu']=individual['donor/La_acetate']['error'];slopes={};relative={}
    for row in bundle['donor']['mechanical_modes']['qualified_two_step']:
        direction=row['direction'];slopes[direction]={}
        for tag,step in [('005',.005),('010',.01)]:
            plus,minus=slope_ids(direction,tag);model=(pred[plus]['energy']-pred[minus]['energy'])/(2*step)
            ep=bundle['donor']['endpoint_rows'][plus];em=bundle['donor']['endpoint_rows'][minus]
            target=row['reference_derivative'] if tag=='005' else (ep['target_relative_energy_kcal_mol']-em['target_relative_energy_kcal_mol'])/(2*step)
            r=residual(model,target);slopes[direction][tag]={**r,'step_A':step,'units':'kcal/mol/A','primary_fit':tag=='005'}
            if tag=='005':groups['slope'].append(r['error'])
    for id,row in bundle['donor']['endpoint_rows'].items():
        relative[id]={**residual(pred[id]['energy']-pred[DONOR]['energy'],row['target_relative_energy_kcal_mol']),'units':'kcal/mol','coordinate':row['coordinate'],'method':row['relative_energy_method']}
    mse={k:float(np.mean(np.square(v))) for k,v in groups.items()};families={'chi':mse['chi'],'mu':(mse['pair_mu']+mse['W9_mu3']+mse['donor_mu'])/3,'E2':mse['E2'],'E3':mse['E3'],'slope':mse['slope']}
    contrast={}
    for id,la in trimers.items():
        if '_La__' not in id:continue
        dyid=id.replace('_La__','_Dy__');dy=trimers[dyid]
        contrast[id.replace('_La__','_Dy_minus_La__')]={k:residual(np.array(dy[k]['prediction'])-np.array(la[k]['prediction']),np.array(dy[k]['target'])-np.array(la[k]['target'])) for k in ('energy_kcal_mol','mu_Debye')}
    water_radial={}
    for metal in ('La','Dy'):
        center=pred['water_central_'+metal]['energy'];water_radial[metal]={name:{'energy_kcal_mol':pred['water_'+name+'_'+metal]['energy'],'relative_to_center_kcal_mol':pred['water_'+name+'_'+metal]['energy']-center} for name in ('minus','central','plus')}
    return {'score':sum(families.values())/5,'family_MSE':families,'mu_subfamily_MSE':{k:mse[k] for k in ('pair_mu','W9_mu3','donor_mu')},
            'family_RMS_in_declared_units':{k:math.sqrt(v) for k,v in families.items()},'pairs':pairs,'trimers':trimers,'trimer_Dy_minus_La':contrast,
            'chi':chi,'individual_signed_moments':individual,'donor_slopes':slopes,'donor_relative_energies':relative,
            'water_radial_energy_curves':water_radial,'water_total_force_prediction':'unavailable from frozen operator cache; no claim of force preservation',
            'donor_canonical_rotation_sensitivity':bundle['donor']['canonical_rotation_corroboration_not_mixed_loss'],
            'Dy_donor_report_only_reference':bundle['donor_Dy_report_only'],
            'score_is_physical_accuracy_gate':False,'all_data_consumed':True}

class JointProblem:
    """Analytic coefficient derivatives of the full five-family objective."""
    def __init__(self,states,bundle,settings):
        self.states=states;self.bundle=bundle;self.settings=settings
        self.hp=HProblem(states,bundle,settings);self.last=None;self.cache=None
        self.evaluations=0;self.singular_trial_objectives=0
    def constraint(self,x):return self.hp.constraint(x[HIDX])
    def constraint_jac(self,x):
        J=np.zeros((len(self.states),15));J[:,HIDX]=self.hp.constraint_jac(x[HIDX]);return J
    def fun_jac(self,x):
        if self.last is not None and np.array_equal(x,self.last):return self.cache
        self.evaluations+=1;b=self.bundle;units=b['units'];D=units['debye_per_eA']
        conv=units['hartree_kcal_mol']/(units['electric']*units['bohr_A']**2)
        pred=static_values(b)
        for row in pred.values():row.update(energy_jac=np.zeros(15),mu_jac=np.zeros((3,15)))
        try:
            for id,s in self.states.items():
                G=operator(s,x[HIDX]);h=s['h']@(x[hIDX]*.1)
                ans=np.linalg.solve(G,np.column_stack((s['Ed']+h,s['Ep']+h,s['B'])))
                u,v=ans[:,0],ans[:,1];Z=ans[:,2:5]
                E=s['other']+s['V']@(x[VIDX]*100.)-.5*s['K']*u@(s['Ep']+h)
                dE=np.zeros(15);dE[VIDX]=s['V']*100.;dE[hIDX]=-.5*s['K']*(u+v)@s['h']*.1
                dE[HIDX]=.5*s['K']*np.einsum('j,ijk,k->i',u,s['H'],v)
                mu=(s['permanent']+.5*s['B'].T@(u+v))*D;dm=np.zeros((3,15))
                dm[:,hIDX]=Z.T@s['h']*.1*D
                dm[:,HIDX]=-.5*np.einsum('ja,ijk,k->ai',Z,s['H'],u+v)*D
                pred[id]={'energy':float(E),'energy_jac':dE,'mu':mu,'mu_jac':dm,
                          'chi_source_au':s['B'].T@Z*conv,'Z':Z}
            families={k:([],[]) for k in ('chi','pair_mu','W9_mu3','donor_mu','E2','E3','slope')}
            def add(name,r,J):
                families[name][0].extend(np.asarray(r).reshape(-1));families[name][1].extend(np.asarray(J).reshape((-1,15)))
            for row in b['chi_targets']:
                a=np.array(row['direction']);id=row['fixture_id'];z=pred[id]['Z']@a;J=np.zeros(15)
                J[HIDX]=-np.einsum('j,ijk,k->i',z,self.states[id]['H'],z)*conv
                add('chi',a@pred[id]['chi_source_au']@a-row['target_a0_cubed'],J)
            for name,cycles in [('pair',b['pair_cycles']),('trimer',[c for c in b['trimer_cycles'].values() if c['primary_fit']])]:
                for c in cycles:
                    add('E2' if name=='pair' else 'E3',cycle_value(c['model_terms'],pred,'energy')-c['energy_kcal_mol'],cycle_value(c['model_terms'],pred,'energy_jac'))
                    add('pair_mu' if name=='pair' else 'W9_mu3',cycle_value(c['model_terms'],pred,'mu')-np.array(c['stored_signed_mu_Debye']),cycle_value(c['model_terms'],pred,'mu_jac'))
            target=np.array(b['donor']['moment_target_au'])*units['bohr_A']*D
            add('donor_mu',pred[DONOR]['mu']-target,pred[DONOR]['mu_jac'])
            for row in b['donor']['mechanical_modes']['qualified_two_step']:
                plus,minus=slope_ids(row['direction'],'005')
                add('slope',(pred[plus]['energy']-pred[minus]['energy'])/.01-row['reference_derivative'],(pred[plus]['energy_jac']-pred[minus]['energy_jac'])/.01)
            loss=0.;gradient=np.zeros(15)
            for name,(r,J) in families.items():
                r=np.array(r);J=np.array(J);weight=1./(15. if name in ('pair_mu','W9_mu3','donor_mu') else 5.)
                loss+=weight*float(np.mean(r*r));gradient+=weight*2.*J.T@r/len(r)
            if not math.isfinite(loss) or not np.isfinite(gradient).all():raise np.linalg.LinAlgError('nonfinite off-domain objective')
        except np.linalg.LinAlgError:
            # Undefined inverse trials are never candidates; the exact SPD endpoint
            # constraint remains unchanged. No matrix, target, or margin is padded.
            self.singular_trial_objectives+=1;loss=1e30;gradient=np.zeros(15)
        self.last=x.copy();self.cache=(loss,gradient);return self.cache

def gradient_check(states,bundle,settings):
    problem=JointProblem(states,bundle,settings);x=np.zeros(15);loss,g=problem.fun_jac(x);rows=[]
    for step in settings['gradient_check_steps']:
        for j in range(15):
            a=x.copy();a[j]=step;c=x.copy();c[j]=-step
            fd=(problem.fun_jac(a)[0]-problem.fun_jac(c)[0])/(2*step)
            bound=settings['gradient_check_relative_allowance']*(1+abs(g[j]))
            rows.append({'column':j,'name':'.'.join(ORDER[j][:2]),'normalized_step':step,'analytic':float(g[j]),'central_difference':float(fd),'absolute_error':float(abs(fd-g[j])),'allowance':float(bound),'passed':bool(abs(fd-g[j])<=bound)})
    return {'status':'passed' if all(r['passed'] for r in rows) else 'failed','zero_loss':loss,'checks':rows,
            'objective_evaluations':problem.evaluations,'new_native_calls':0,'new_quantum_calls':0,
            'scope':'Coefficient Jacobian of fixed cached geometries; not coordinate-force validation.'}

def candidate(states,bundle,K,coeff,settings):
    if not np.isfinite(coeff).all() or not np.all(np.abs(coeff)<=SCALES):
        return {'status':'rejected_parameter_box','coefficients':coeff.tolist()}
    mm=margins(states,coeff[HIDX],settings['normalized_SPD_floor'])
    if min(mm.values())<0:return {'status':'rejected_SPD','coefficients':coeff.tolist(),'SPD_margins':mm}
    pred,records=full_evaluate(states,bundle,K,coeff)
    stable=min(r['minimum_normalized_eigenvalue'] for r in records.values())>=settings['normalized_SPD_floor']
    out={'status':'complete_stable_candidate' if stable else 'rejected_final_kernel_SPD',
         'coefficients':coeff.tolist(),'normalized_coefficients':(coeff/SCALES).tolist(),
         'at_exact_box_boundary':(np.abs(coeff)==SCALES).tolist(),'SPD_margins':mm,
         'predictions':records,'assessment':assess(pred,bundle)}
    return out

def solve_start(index,start,states,bundle,K,settings,checkpoint=None):
    tick=time.monotonic();cpu=time.process_time();result={'index':index,'prospective_H_start':list(start),'status':'started','error':None,'candidates':{},'solver_attempts':[]}
    def save(stage):
        result['stage']=stage
        if stage.endswith('_solver_intent'):result['solver_attempts'].append(stage)
        if checkpoint is not None:checkpoint(result)
    try:
        x=np.array(start,dtype=float);shrinks=0
        while min(margins(states,x,settings['normalized_SPD_floor']).values())<0 and shrinks<settings['start_max_halvings']:
            x=x*.5;shrinks+=1
        result.update(actual_H_start=x.tolist(),start_halvings=shrinks)
        if min(margins(states,x,settings['normalized_SPD_floor']).values())<0:
            result['status']='unavailable_initial_SPD';return result
        save('H_solver_intent');problem=HProblem(states,bundle,settings)
        fit=minimize(problem.fun_jac,x,jac=True,method='SLSQP',bounds=[(-1.,1.)]*4,
                     constraints=[{'type':'ineq','fun':problem.constraint,'jac':problem.constraint_jac}],
                     options={'maxiter':settings['H_max_iter'],'ftol':settings['H_ftol'],'disp':False})
        H=fit.x;end=margins(states,H,settings['normalized_SPD_floor'])
        result['H_solver']={'status':int(fit.status),'success':bool(fit.success),'message':str(fit.message),'iterations':int(fit.nit),
                            'objective':float(fit.fun),'solution':H.tolist(),'operator_evaluations':problem.evaluations,
                            'singular_off_domain_objectives':problem.singular_trial_objectives,'SPD_margins':end,
                            'interpretation':'Finite local initialization; unsuccessful raw status retained.'}
        save('H_solver_returned')
        if not np.isfinite(H).all() or not np.all(np.abs(H)<=1) or min(end.values())<0:
            result['status']='rejected_H_endpoint';return result
        save('h_solver_intent');h,result['h_solver']=h_fit(states,bundle,H,settings);save('h_solver_returned')
        save('V_solver_intent');V,result['V_solver']=V_fit(states,bundle,H,h,settings);save('V_solver_returned')
        coeff=np.zeros(15);coeff[HIDX]=H;coeff[hIDX]=h;coeff[VIDX]=V
        result['candidates']['conditional']=candidate(states,bundle,K,coeff,settings);save('conditional_complete')
        if result['candidates']['conditional']['status']!='complete_stable_candidate':
            result['status']='rejected_conditional_candidate';return result
        joint=JointProblem(states,bundle,settings);init=coeff/SCALES
        analytic_loss=joint.fun_jac(init)[0];literal_loss=result['candidates']['conditional']['assessment']['score']
        result['conditional_objective_join']={'analytic_loss':analytic_loss,'kernel_reported_loss':literal_loss}
        if abs(analytic_loss-literal_loss)>1e-9*(1+abs(literal_loss)):raise ValueError('joint objective and frozen kernel report disagree')
        save('joint_solver_intent')
        try:
            fit=minimize(joint.fun_jac,init,jac=True,method='SLSQP',bounds=[(-1.,1.)]*15,
                         constraints=[{'type':'ineq','fun':joint.constraint,'jac':joint.constraint_jac}],
                         options={'maxiter':settings['joint_max_iter'],'ftol':settings['joint_ftol'],'disp':False})
            result['joint_solver']={'status':int(fit.status),'success':bool(fit.success),'message':str(fit.message),'iterations':int(fit.nit),
                    'objective':float(fit.fun),'normalized_solution':fit.x.tolist(),'evaluations':joint.evaluations,
                    'singular_off_domain_objectives':joint.singular_trial_objectives,'interpretation':'One bounded nonconvex refinement; no global-optimum claim.'}
            save('joint_solver_returned');result['candidates']['joint']=candidate(states,bundle,K,fit.x*SCALES,settings)
        except Exception as exc:result['joint_error']=repr(exc)
        eligible=[(c['assessment']['score'],0 if name=='conditional' else 1,name,c) for name,c in result['candidates'].items() if c['status']=='complete_stable_candidate']
        _,_,name,best=min(eligible,key=lambda row:row[:2])
        result.update(status='complete_stable_candidate',selected_stage=name,coefficients=best['coefficients'],assessment=best['assessment'],
                      coefficient_order=['.'.join(x[:2]) for x in ORDER],optimizer_success_required_for_feasible_candidate=False)
    except LinearStageFailure as exc:
        result.update(status='failed_linear_stage',error=str(exc),failed_linear_stage=exc.evidence)
    except Exception as exc:result.update(status='failed',error=repr(exc))
    finally:
        result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
        save('terminal')
    return result
