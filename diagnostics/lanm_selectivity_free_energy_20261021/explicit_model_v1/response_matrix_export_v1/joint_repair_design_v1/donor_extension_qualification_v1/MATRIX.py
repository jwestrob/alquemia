#!/usr/bin/env python3
"""Narrow donor adaptation of admitted scout6541462b. No native execution."""
import math
import re
import numpy as np

def qualify_exports(matrix,export_paths,row,units,new_moments,report=None):
    """Pinned scout arithmetic with unequal donor duals retained as observations."""
    require=matrix.require; parse_export=matrix.parse_export
    induction_values=matrix.induction_values; arithmetic=matrix.arithmetic
    verify=matrix.verify; half_print=matrix.half_print
    report={} if report is None else report;report.update(checks={},exports=[],field_predictions={})
    def gate(name,error,limit):
        passed=math.isfinite(float(error)) and math.isfinite(float(limit)) and 0<=error<=limit
        report['checks'][name]={'error':float(error),'bound':float(limit),'passed':bool(passed)}
        require(passed,'Numerical export qualification: '+name)
    require(len(export_paths)==2,'Exactly two exports')
    parsed=[parse_export(p,row,units) for p in export_paths];report['exports']=parsed
    require([p['sequence'] for p in parsed]==[1,2] and parsed[0]['sites']==parsed[1]['sites'],'Two solve/site identities')
    q=induction_values(new_moments['analyze_output']);n=row['atoms'];dim=3*n;D=units['debye_per_eA'];K=units['electric'];B=units['bohr_A'];H=units['hartree_kcal_mol']
    for k,p in enumerate(parsed):
        A=np.array(p['A']);v=np.array(p['vectors']);ed,ep,u,w=v.T
        alpha=np.repeat([s['alpha_A3'] for s in p['sites']],3);sqrt=np.sqrt(alpha)
        ar=arithmetic(dim,A,ed,ep,u,w);an=math.sqrt(dim)*ar
        sym=float(np.max(abs(A-A.T)));gate(f'{k}_symmetry',sym,ar)
        lam=np.linalg.eigvalsh(A);scaled=np.linalg.eigvalsh(sqrt[:,None]*A*sqrt[None,:])
        require(np.isfinite(lam).all() and lam[0]>0 and scaled[0]>1e-6*scaled[-1]>0,'Full matrix stability failed')
        p['spectrum']={'A':lam.tolist(),'normalized_A':scaled.tolist(),'condition_A':float(lam[-1]/lam[0]),'condition_normalized':float(scaled[-1]/scaled[0])}
        ud=np.linalg.solve(A,ed);wd=np.linalg.solve(A,ep)
        r=A@u-ed;rp=A@w-ep;dr=A@ud-ed;dp=A@wd-ep
        amp=float(np.linalg.norm(np.eye(dim)-p['peek']*A*alpha[None,:],2))
        cg=amp*math.sqrt(n)*(q[k]+.5e-10)/D+an
        gate(f'{k}_Au_Ed',np.linalg.norm(r),cg);gate(f'{k}_Av_Ep',np.linalg.norm(rp),cg)
        bound_u=(np.linalg.norm(r)+np.linalg.norm(dr)+2*an)/lam[0]
        bound_v=(np.linalg.norm(rp)+np.linalg.norm(dp)+2*an)/lam[0]
        gate(f'{k}_native_vs_dense_u',np.linalg.norm(u-ud),bound_u);gate(f'{k}_native_vs_dense_v',np.linalg.norm(w-wd),bound_v)
        p['dual_differences_observed']={'RHS_max':float(np.max(abs(ed-ep))),
            'solution_norm':float(np.linalg.norm(u-w))}
        # Donors retain both physically distinct native direct-field channels.
        for j,name in enumerate(('channel_operator','direct_RHS','native_symmetry','native_operator_action')):gate(f'{k}_{name}',p['checks'][j],ar)
        gate(f'{k}_native_residual_diagnostic_u',abs(p['checks'][4]-np.max(abs(r))),ar)
        gate(f'{k}_native_residual_diagnostic_v',abs(p['checks'][5]-np.max(abs(rp))),ar)
        p['numerical_bounds']=dict(arithmetic=ar,arithmetic_norm=an,post_peek_amplification=amp,CG_post_peek_residual_norm=cg,native_vs_dense_u=bound_u,native_vs_dense_v=bound_v)
        if k==0:
            text=verify(new_moments['analyze_output']).read_text()
            token=re.findall(r'^\s*Polarization\s+([-+\d.EeDd]+)\s+\d+\s*$',text,re.M)
            if token:require(len(token)==1,'Polarization component denominator');printed=float(token[0]);print_error=half_print(token[0])
            else:
                require('Polarization' not in new_moments['components'] and not np.any(ed) and not np.any(ep),'Missing nonzero polarization component')
                printed=0.;print_error=0.
            native_E=-.5*K*float(u@ep);dense_E=-.5*K*float(ud@ep);ea=arithmetic(dim,native_E,printed)
            gate('native_polarization_energy',abs(native_E-printed),print_error+ea)
            gate('dense_polarization_energy',abs(dense_E-printed),print_error+ea+.5*K*np.linalg.norm(ep)*bound_u)
            report['polarization_energy']=dict(native_contraction=native_E,dense_contraction=dense_E,printed=printed,print_half_unit=print_error)
        else:
            mapped=np.empty((n,3))
            for pos,s in enumerate(p['sites']):mapped[s['atom']-1]=D*u[3*pos:3*pos+3]
            gate('native_uind_site_print',np.max(abs(mapped-np.array(new_moments['induced_Debye']))),.5e-6+D*ar)
    a,b=parsed;gate('two_solve_matrix',np.max(abs(np.array(a['A'])-b['A'])),arithmetic(dim,a['A'],b['A']))
    gate('two_solve_RHS',np.max(abs(np.array(a['vectors'])[:,:2]-np.array(b['vectors'])[:,:2])),arithmetic(dim,a['vectors'],b['vectors']))
    report['two_solve_exact']={'matrix':a['A']==b['A'],'RHS':np.array_equal(np.array(a['vectors'])[:,:2],np.array(b['vectors'])[:,:2])}
    report['passed']=all(x['passed'] for x in report['checks'].values());return report
