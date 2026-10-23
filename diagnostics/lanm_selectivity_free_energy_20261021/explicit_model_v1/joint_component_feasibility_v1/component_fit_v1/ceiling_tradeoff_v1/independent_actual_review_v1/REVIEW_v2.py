"""Saved-only tradeoff review; imports independent arithmetic, never a solver."""
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback
from decimal import Decimal, localcontext

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
cfg = json.loads((HERE/'CONFIG_v2.json').read_text())
import hashlib
def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
assert pin(cfg['helper']['path']) == cfg['helper']
spec = importlib.util.spec_from_file_location('independent_conic_arithmetic', cfg['helper']['path'])
h = importlib.util.module_from_spec(spec); spec.loader.exec_module(h)
req, read, near = h.require, h.read, h.near
REMOVED = {'block/'+k for k in h.OLD}


def matrix(ref):
    req(pin(ref['path']) == ref, 'Exact sparse file')
    with np.load(ref['path'], allow_pickle=False) as a:
        return (sparse.csc_matrix((a['A_data'], a['A_indices'], a['A_indptr']), shape=tuple(a['A_shape'])),
                a['b'].copy(), a['c'].copy())


def certificate(p, slot, d, q, success):
    r = h.scalar_report(p, d, q, success); c = slot['certificate']
    for key, value in [('residual_constraints', r['margins']), ('blocks', r['blocks']), ('original25', r['original25'])]:
        req(c[key] == value, 'Exact independent scalar '+key)
    req(len(r['margins']) == 851 and len(r['original25']) == 25, '851 scalar constraints/25 gates')
    req(slot['parameters'] == r['parameters'] and slot['key'] == h.digest(r['parameters']), 'Exact physical point/key')
    req(slot['eligible'] is r['eligible'], 'Original eligibility preserved')
    expected = {'method': 'independent_scalar_fsum', 'solver_success': success,
        'failed_solver_status_remains_failed': not success, 'q_fraction': q, 'q_current': p['q_current'],
        'strict_q_improvement': q < p['q_current'], 'in_box_and_trust': r['valid_box'],
        'minimum_box_margin': r['minimum_box_margin'], 'minimum_residual_constraint_margin': r['minimum_margin'],
        'minimum_epigraph_margin': min(r['margins']['block/'+k] for k in h.COMPONENTS),
        'all_residual_margins_strictly_positive': r['minimum_margin'] > 0,
        'all25_predicted': r['all25'], 'fixed_old_ceiling_pass': r['ceilings']}
    for k, v in expected.items(): req(c[k] == v, 'Exact certificate summary '+k)
    # Bare CORE slots contain only fraction/parameters/key/eligible/certificate.
    # Tradeoff comparison epsilon and native_evaluated are checked on each outer witness.
    return r


def decimal_margins(p, d, q):
    """Higher precision of exact saved binary64 inputs; no feasibility tolerance."""
    with localcontext() as ctx:
        ctx.prec = 70
        D = lambda x: Decimal.from_float(float(x))
        ds = list(map(D,d)); margins = {}
        def affine(block):
            return [D(r)+sum((D(a)*v for a,v in zip(row,ds)),Decimal(0)) for r,row in zip(block['r'],block['J'])]
        for name,block in p['blocks'].items():
            values=affine(block)
            margins['block/'+name]=(D(q) if name in h.COMPONENTS else Decimal(1))-sum(v*v for v in values)/D(block['denominator'])
        for name,v in zip(p['energy']['names'],affine(p['energy'])):
            margins['energy/upper/'+name]=1-v; margins['energy/lower/'+name]=1+v
        for name,block in p['forces'].items():
            values=affine(block)
            margins['force/RMS/'+name]=1-sum(D(w)*v*v for w,v in zip(block['weights'],values))/4
            for i,v in enumerate(values):
                margins['force/upper/'+name+'/'+str(i)]=1-v/5
                margins['force/lower/'+name+'/'+str(i)]=1+v/5
        retained={k:v for k,v in margins.items() if k not in REMOVED}
        bounds=[v-D(lo) for v,lo in zip(ds,p['bounds_lower'][:12])]+[D(hi)-v for v,hi in zip(ds,p['bounds_upper'][:12])]
        physical=[D(c)+D(w)*v for c,w,v in zip(p['center'],p['width'],ds)]
        bounds += [v-D(lo) for v,lo in zip(physical,p['lower'])]+[D(hi)-v for v,hi in zip(physical,p['upper'])]
        bounds += [D(q),D(p['q_current'])-D(q)]
        req(len(retained)==847 and min(retained.values())>0 and min(bounds)>=0, 'Decimal retained signs and domain')
        return {'minimum_retained_margin':str(min(retained.values())), 'limiting_constraint':min(retained,key=retained.get),
                'minimum_domain_margin':str(min(bounds)), 'all847_strictly_positive':True, 'precision_decimal_digits':70,
                'scope':'Diagnostic high-precision arithmetic on saved binary64 coefficients; not an interval or nonlinear certificate'}


def run():
    tick,cpu=time.monotonic(),time.process_time()
    out=Path(cfg['output']);out.mkdir(parents=True,exist_ok=False)
    report={'status':'started','implementation':pin(__file__),'configuration':pin(HERE/'CONFIG_v2.json'),
            'result':cfg['result'],'helper':cfg['helper'],'new_solver_calls':0,'new_native_calls':0,'error':None,
            'authorship_disclosure':'Reviewer authored original conic producer; this new tradeoff was authored by Arendt. Reused independent reviewer scalar helpers are pinned and do not import producer/CORE/optimizer.'}
    try:
        req(all(os.environ.get(k)=='1' for k in h.THREADS),'One numerical thread')
        r=read(cfg['result']);c=read(r['configuration']);inputs=read(c['inputs'])
        req(r['status']=='completed_saved_ceiling_tradeoff' and r['error'] is None,'Complete actual result')
        req(pin(r['implementation']['path'])==r['implementation']==c['implementation'],'Producer pin')
        source=read(c['source_review']);release=read(c['root_release'])
        fp=h.digest({k:v for k,v in c.items() if k not in ('source_review','root_release')})
        req(fp==r['configuration_fingerprint']==source['configuration_fingerprint']==release['configuration_fingerprint'],'Fingerprint joins')
        req(source['status']=='passed_ceiling_tradeoff_source_review' and source['implementation']==r['implementation'] and source['inputs']==r['inputs']==c['inputs'],'Source admission')
        req(release['status']=='released_one_saved_ceiling_tradeoff' and release['source_review']==c['source_review'] and release['baseline_certificate']==r['baseline_certificate']==c['baseline_certificate'],'Exact release')
        req(release['maximum_solver_calls']==r['solver_calls']==1 and release['new_native_calls']==r['new_native_calls']==r['new_candidate_evaluations']==0 and r['selected_candidate'] is None,'Counts/no promotion')
        for ref in inputs.values():
            if isinstance(ref,dict) and set(ref)=={'path','sha256'}: req(pin(ref['path'])==ref,'Input pin')
        req(pin(c['witness_contract']['path'])==c['witness_contract'],'Auxiliary witness contract')
        full=read(inputs['full_result']);admission=read(inputs['full_actual_review']);baseline=read(c['baseline_certificate'])
        req(admission['status']=='passed_saved_conic_actual_review' and admission['certificate_checks_passed'] is True and admission['result']==inputs['full_result']==r['full_comparison'] and admission['implementation']==cfg['helper'],'Original arithmetic authority')
        req(full['solver_status']=='AlmostSolved' and full['solver_success'] is False,'Original reduced status preserved')
        req(all(r[k]==full[k] for k in ('origin_key','derivative_keys','linearized_problem','effective_settings')),'Exact origin/probes/settings/affine model')
        fcfg=read(full['configuration']);req(r['effective_settings']==read(fcfg['settings']),'Pinned effective settings')
        marker=json.loads(Path(c['one_attempt_marker']).read_text())
        req(marker['configuration']==r['configuration'] and marker['baseline_certificate']==c['baseline_certificate'],'One-attempt marker')
        p=read(r['linearized_problem']);h.validate_problem(p)
        A0,b0,c0=matrix(full['conic_data']);A,b,obj=matrix(r['conic_data'])
        layout0=read(full['row_layout']);layout=read(r['row_layout'])
        removed=layout0[:4]
        req([(v['name'],v['start'],v['stop']) for v in removed]==[('block/old_EF',0,414),('block/pair_mu',414,517),('block/C_mu',517,530),('block/chi',530,535)],'Exact four SOC removal')
        req(r['removed_cones']==removed and all(v['kind']=='SOC' for v in removed),'Removed cone metadata')
        req(r['retained_original_row_indices']==list(range(535,1817)),'All1282 original remaining rows retained')
        expected=A0[535:,:].tocsc()
        req(A.shape==(1282,13) and all(np.array_equal(getattr(A,k),getattr(expected,k)) for k in ('data','indices','indptr')) and np.array_equal(b,b0[535:]) and np.array_equal(obj,c0),'Exact sparse slice, no other matrix change')
        req(layout==[dict(v,start=v['start']-535,stop=v['stop']-535) for v in layout0[4:]],'Exact retained row layout')
        req(sum(v['kind']=='SOC' for v in layout)==17 and sum(v['kind']=='nonnegative' for v in layout)==856,'17 SOCs and856 linear sides')
        x,s,z=(np.asarray(r['solution'][k],dtype=float) for k in ('x','s','z'))
        req(x.shape==(13,) and s.shape==z.shape==(1282,) and all(np.isfinite(v).all() for v in (x,s,z)) and x[-1]>=0,'Finite raw solution')
        req(r['solver_status']==r['solution']['status']==r['info']['status']=='AlmostSolved' and r['solver_success'] is False,'AlmostSolved remains false success')
        req(r['returned_q']==float(x[-1]**2),'Raw q is t squared')
        primal=A@x+s-b;dual=A.T@z+obj;pc=float(obj@x);dc=float(-b@z)
        residuals={'primal_equation_residual':primal,'dual_stationarity_residual':dual,'primal_residual_inf':max(abs(primal)),
                   'dual_residual_inf':max(abs(dual)),'primal_objective_t':pc,'dual_objective_t':dc,'signed_gap_t':pc-dc,'complementarity_dot':float(s@z)}
        for k,v in residuals.items():near(v,r['independent_linear_residuals'][k],k)
        for name,vec in [('primal_affine_cone_margins',b-A@x),('raw_slack_cone_margins',s),('dual_cone_margins',z)]:
            observed=h.margins(layout,vec);stored=r['independent_linear_residuals'][name]
            req(set(observed)==set(stored),'Cone margin names')
            near(list(observed.values()),[stored[k] for k in observed],name)
        affine=b-A@x; squared={}
        for row in layout:
            values=affine[row['start']:row['stop']]
            if row['kind']=='SOC': squared[row['name']]=float(values[0]**2-values[1:]@values[1:])
            elif not row['name'].startswith('bound/'): squared[row['name']]=float(values[0])
        raw_independent=h.scalar_report(p,list(map(float,x[:12])),r['returned_q'],False)
        retained_raw={k:v for k,v in raw_independent['margins'].items() if k not in REMOVED}
        req(len(squared)==len(retained_raw)==847 and set(squared)==set(retained_raw),'All847 retained cone/scalar names')
        near(list(squared.values()),[retained_raw[k] for k in squared],'All847 reduced cone/scalar correspondence')
        corr=r['retained_correspondence']
        req(corr['retained_source_inequalities']==847 and corr['nonnegative_t'] is True,'Correspondence inventory')
        near(corr['minimum_retained_source_margin'],min(retained_raw.values()),'Reported retained source minimum')
        req(r['independent_linear_residuals']['rigorously_rounded_dual_certificate'] is False,'Numerical certificate limit')
        base=baseline['baseline'];xb=full['solution']['x'];db=xb[:12]
        req(baseline['status']=='passed_saved_auxiliary_epigraph_baseline' and baseline['error'] is None and baseline['solver_calls']==baseline['native_calls']==0 and baseline['inputs']==c['inputs'],'Saved zero-solve baseline')
        br=certificate(p,base['raw_certificate'],db,float(xb[-1]**2),False)
        rho0=max(br['blocks'][k]['ratio'] for k in h.COMPONENTS);q0=rho0+1e-8*max(1.,abs(rho0))
        bc=certificate(p,base['certificate'],db,q0,False)
        req(bc['minimum_margin']>0 and bc['valid_box'] and bc['all25'] and bc['ceilings'],'Full851 strict baseline')
        req(rho0==base['actual_component_rho']==r['baseline_actual_component_rho'] and q0==base['auxiliary_q']==r['baseline_auxiliary_q'],'Baseline rho/q')
        req([w['fraction'] for w in r['witnesses']]==[1.,.5,.25,.125],'Four declared fractions')
        summaries=[]
        for w in r['witnesses']:
            f=w['fraction'];d=[f*float(v) for v in x[:12]];rawq=(1-f)*p['q_current']+f*r['returned_q']
            raw=certificate(p,w['raw_source_certificate'],d,rawq,False)
            rho=max(raw['blocks'][k]['ratio'] for k in h.COMPONENTS);cushion=1e-8*max(1.,abs(rho));q=rho+cushion
            actual=certificate(p,w['constructed_source_certificate'],d,q,False)
            req(w['raw_source_certificate']['fraction']==f and w['constructed_source_certificate']['fraction']==1.,'Auxiliary construction fraction metadata')
            req(actual['parameters']==raw['parameters'] and w['physical_fraction_changed'] is False,'Only q changed')
            retained={k:v for k,v in actual['margins'].items() if k not in REMOVED}
            feasible=actual['valid_box'] and actual['all25'] and all(v>0 for v in retained.values())
            eps=1e-10*max(1.,abs(rho0));gain=rho0-rho
            expected_fields={'raw_auxiliary_q':rawq,'raw_minimum_epigraph_margin':min(raw['margins']['block/'+k] for k in h.COMPONENTS),
              'auxiliary_cushion':cushion,'retained_scalar_margin_count':847,'minimum_retained_margin':min(retained.values()),
              'strict_retained_feasible':feasible,'actual_component_rho':rho,'auxiliary_q':q,'baseline_actual_component_rho':rho0,
              'gain_vs_baseline_actual_rho':gain,'gain_comparison_epsilon':eps,'feasible_gain_exceeds_numerical_comparison':feasible and rho<rho0-eps,
              'removed_ceiling_ratios':{k:actual['blocks'][k]['ratio'] for k in h.OLD},'native_evaluated':False,'candidate_admitted':False,'nonlinear_or_global_optimality_claim':False}
            for k,v in expected_fields.items():req(w[k]==v,'Witness field '+k)
            req(len(retained)==847 and feasible,'Strict original847 and25/domain')
            dec=decimal_margins(p,d,q)
            summaries.append({'fraction':f,'rho':rho,'auxiliary_q':q,'gain':gain,'strict_retained_feasible':feasible,
              'minimum_retained_margin':min(retained.values()),'minimum_domain_margin':actual['minimum_box_margin'],
              'all25_passed':actual['all25'],'old_score_ratios':expected_fields['removed_ceiling_ratios'],
              'six_squared_norms':{k:v['predicted_squared_norm'] for k,v in actual['blocks'].items()},
              'physical_key':w['constructed_source_certificate']['key'],'decimal_confirmation':dec})
        policy=read(cfg['nonlinear_policy'])
        req(policy['created_utc']<release['created_utc'],'Followup policy prospectively predates solve release')
        rule=policy['selection_rule'];req(rule['order']==[1.,.5,.25,.125] and rule['maximum_selected_vectors']==1 and rule['minimum_absolute_gain_in_component_q']==.05,'Exact prospective trigger')
        eligible=[v for v in summaries if v['strict_retained_feasible'] and v['gain']>=.05]
        req(eligible and eligible[0]['fraction']==1.,'Full fraction first satisfies prospective trigger')
        costs=r['cost'];req(costs['numerical_threads']==1 and costs['scheduler_job_id'] is None and all(math.isfinite(costs[k]) and costs[k]>0 for k in ('wall_seconds','CPU_seconds')),'Actual truthful cost envelope')
        req(r['info']['iterations']==r['solution']['iterations']==17 and r['info']['linsolver']['threads']==1,'Actual solver iteration/thread metadata')
        req('clarabel' not in sys.modules,'No optimizer imported')
        report.update(status='passed_ceiling_tradeoff_actual_arithmetic_review',source_review=c['source_review'],
          certificate_checks_passed=True,solver_status=r['solver_status'],solver_success=False,
          removed_rows=535,retained_rows=1282,retained_constraints_per_witness=847,
          baseline_rho=rho0,returned_q=r['returned_q'],witnesses=summaries,
          prospective_policy=cfg['nonlinear_policy'],prospective_trigger_fraction=eligible[0]['fraction'],
          prospective_trigger_witness={'fraction':eligible[0]['fraction'],'physical_key':eligible[0]['physical_key'],
            'parameters':r['witnesses'][0]['constructed_source_certificate']['parameters']},
          raw_primal_residual_inf=float(max(abs(primal))),raw_dual_residual_inf=float(max(abs(dual))),signed_gap_t=pc-dc,
          actual_solver_cost=costs,raw_solver_seconds=r['solution']['solve_time'],
          interpretation='Exact saved local tradeoff and strict numerical witnesses verified; AlmostSolved remains AlmostSolved. Full fraction meets the predeclared followup trigger. No candidate/native release, nonlinear accuracy, optimality or physical admission.')
    except Exception as exc:
        report.update(status='failed_ceiling_tradeoff_actual_arithmetic_review',error=repr(exc),traceback=traceback.format_exc())
    report['checks']=h.CHECKS;report['review_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}
    dest=out/'REVIEW.json';dest.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],'review':pin(dest),'error':report['error'],'checks':h.CHECKS}))
    return report

if __name__=='__main__':raise SystemExit(bool(run()['error']))
