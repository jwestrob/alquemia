"""Gated two-coordinate mechanics assessment; never enables an unvalidated score."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new
from affordable_response import source_key
from ggr_sensitivity import executed
from mace_hybrid import EV_TO_KCAL
from mace_curvature import actual_collection,dft_source
from mace_mechanics import H,UNAVAILABLE
from mace_mechanics_run import preparation,validate


def matrix(values,half=False):
    prefix='half_' if half else '';step=.5 if half else 1.
    def e(label):return float(values[prefix+label])-float(values['center'])
    a=(e('s_plus')+e('s_minus'))/step**2
    b=(e('theta_plus')+e('theta_minus'))/step**2
    c=(e('same_plus')+e('same_minus')-e('opposite_a')-e('opposite_b'))/(4*step**2)
    result=np.array([[a,c],[c,b]])
    if not np.isfinite(result).all():raise InvalidArtifact('nonfinite mechanical matrix')
    return result


def square_extrema(a):
    """Exact range of a homogeneous two-dimensional quadratic on [-1,1]^2."""
    a=np.asarray(a,dtype=float)
    if a.shape!=(2,2) or not np.isfinite(a).all() or not np.array_equal(a,a.T):
        raise InvalidArtifact('finite symmetric2x2 matrix required')
    points=[np.array([x,y]) for x in (-1.,1.) for y in (-1.,1.)]+[np.zeros(2)]
    for fixed in (-1.,1.):
        if a[1,1]!=0:
            y=-a[0,1]*fixed/a[1,1]
            if abs(y)<=1:points.append(np.array([fixed,y]))
        if a[0,0]!=0:
            x=-a[0,1]*fixed/a[0,0]
            if abs(x)<=1:points.append(np.array([x,fixed]))
    values=[float(.5*z@a@z) for z in points]
    return {'minimum_kcal_mol':min(values),'maximum_kcal_mol':max(values),
            'max_absolute_kcal_mol':max(map(abs,values)),
            'checked_points':[z.tolist() for z in points]}


def stationary_point(a,b,coarse):
    a=np.asarray(a);b=np.asarray(b);coarse=np.asarray(coarse)
    if a.shape!=(2,2) or b.shape!=(2,) or not np.isfinite(a).all() or not np.isfinite(b).all() or not np.isfinite(coarse).all():
        raise InvalidArtifact('nonfinite/invalid bounded mechanical model')
    fine_eigen=np.linalg.eigvalsh(a);coarse_eigen=np.linalg.eigvalsh(coarse)
    condition=float(np.linalg.cond(a))
    result={'fine_eigenvalues_kcal_mol':fine_eigen.tolist(),'coarse_eigenvalues_kcal_mol':coarse_eigen.tolist(),
            'coordinate_measure':'z=(s/0.02Angstrom,theta/(pi/180radian))',
            'condition_number':condition if math.isfinite(condition) else None,'normalized_optimum':None,
            'predicted_relaxation_kcal_mol':None,'status':'unstable_or_unresolved_curvature'}
    if min(fine_eigen)<=0 or min(coarse_eigen)<=0:return result
    z=-np.linalg.solve(a,b)
    residual=float(np.linalg.norm(a@z+b)/max(float(np.linalg.norm(b)),np.finfo(float).tiny))
    result.update(normalized_optimum=z.tolist(),physical_optimum=(H*z).tolist(),solve_relative_residual=residual)
    if not np.isfinite(z).all() or residual>1e-10:
        result['status']='unreliable_linear_solve';return result
    if np.max(np.abs(z))>1.+1e-12:
        result['status']='trust_region_exceeded';return result
    result.update(status='eligible_for_independent_DFT_minimum_validation',
                  predicted_relaxation_kcal_mol=float(.5*b@z))
    return result


def verified_new_dft(path):
    c=read_json(path);mp=verify(c['manifest']);m,rows=executed(mp)
    if c['status']!='complete' or any(r['status']!='complete' for r in rows.values()):
        raise InvalidArtifact('all declared new DFT endpoints must complete before assessment')
    for key,r in rows.items():
        if any(c['rows'][key].get(k)!=v for k,v in r.items()):
            raise InvalidArtifact('DFT collected result differs from verified execution')
    for g in c['gradients'].values():
        for ref in g['artifacts'].values():verify(ref)
    return c,m


def all_core_rows(c,p,source_reference):
    original,_=actual_collection(source_reference)
    rows=dict(c['rows'])
    for key,r in p['reused_GGR'].items():
        stored=c['reused_rows'][key]
        actual=original['rows'][r['source_task_id']]
        if stored['result']!=actual or stored['source_collection']!=record(source_reference):
            raise InvalidArtifact('reused component does not match actual archive')
        rows[key]=actual
    if len(rows)!=136:raise InvalidArtifact('incomplete core point inventory')
    return rows


def assess(prepared,core,gb,short,dft,normal_manifest,convention,output):
    p=preparation(prepared);mc,mm=actual_collection(core);bc,bm=actual_collection(gb);sc,sm=actual_collection(short)
    for c,m in ((mc,mm),(bc,bm),(sc,sm)):
        validate(verify(c['manifest']))
        if m['preparation']!=record(prepared):raise InvalidArtifact('calculation preparation mismatch')
        submission=read_json(Path(verify(c['manifest'])).parent/'submission.json')
        if (submission['manifest']!=c['manifest'] or submission['returncode']!=0 or
                submission['numerical_convention']!=record(convention)):
            raise InvalidArtifact('numerical convention differs from the pre-execution record')
    if bm['source_mace_collection']!=record(core):raise InvalidArtifact('GB does not use these MACE densities')
    dc,dm=verified_new_dft(dft)
    if record(verify(dc['manifest']))!=p['DFT_tasks']:raise InvalidArtifact('wrong DFT validation panel')
    old_dft,old_dft_m=dft_source(verify(p['DFT_reference']))
    normal_m,normal_rows=executed(normal_manifest)
    short_gate,_=actual_collection(verify(sm['short_gate']))
    if not short_gate['numerical_checks_pass']:raise InvalidArtifact('short analytic derivative gate failed')
    mr=all_core_rows(mc,p,verify(mm['reference_mace']));br=all_core_rows(bc,p,verify(mm['reference_gb']))
    rows={};actual_dft={};gradients={};normal_bridge={}
    for name,ref in p['cases'].items():
        c=read_json(verify(ref));cj=np.load(verify(c['core_jacobians']));fj=np.load(verify(c['full_jacobians']))
        physical=read_json(verify(c['physical_atoms']));values={};b_vectors={};g_parts={};case_bridge={}
        for metal in ('La','Ca'):
            key=f'{name}_center_{metal}';actual_dft[key]=None
            if key in p['reused_GGR']:
                original=p['reused_GGR'][key]['source_task_id'];d0=old_dft['energies'][original]['energy_hartree'];gradient=old_dft['gradients'][original]
                short_core=short_gate['rows']['core__'+original]
            else:
                d0=dc['rows'][key]['energy_hartree'];gradient=dc['gradients'][key]
                short_core=sc['rows']['core__'+key]
            center=c['grids']['center']['endpoints'][metal]
            normal_matches=[t for t in normal_m['tasks'] if t['xyz']['sha256']==center['xyz']['sha256']]
            if len(normal_matches)!=1:raise InvalidArtifact('exact normal-SCF center reference unavailable/ambiguous')
            nr=normal_rows[normal_matches[0]['task_id']]
            if nr['status']!='complete':raise InvalidArtifact('normal center reference not complete')
            bridge=(d0-nr['energy_hartree'])*HA_TO_KCAL
            case_bridge[metal]={'difference_kcal_mol':bridge,'pass':abs(bridge)<=.05,'source':nr}
            g=np.array(gradient['gradient_kcal_mol_per_A']);gd=np.einsum('nij,ij->n',cj,g)
            source_g=gradient['mapped_source_gradient_kcal_mol_per_A']
            moving={key for spec in c['specs'] for key in spec['moving_source_keys']}
            if not moving.issubset(source_g):raise InvalidArtifact('mapped gradient is missing a moving physical atom')
            # Atoms outside the quantum energy have exactly zero contribution
            # to its derivative; the separate J gradient supplies the scaffold.
            mapped=np.array([source_g.get(a['source_key'],[0.,0.,0.]) for a in physical])
            if not np.allclose(gd,np.einsum('nij,ij->n',fj,mapped),atol=1e-7,rtol=0):
                raise InvalidArtifact('physical/QM DFT derivative projection differs')
            full_force=np.load(verify(sc['rows'][f'full__{c["global_id"]}_center_{metal}']['forces']))
            core_force=np.load(verify(short_core['forces']))
            gj=(-np.einsum('nij,ij->n',fj,full_force)+np.einsum('nij,ij->n',cj,core_force))*EV_TO_KCAL
            b_vectors[metal]=H*(gd+gj);g_parts[metal]={'DFT_core':gd.tolist(),'J_scaffold':gj.tolist(),'combined':(gd+gj).tolist()}
            values[metal]={'C_core':{},'J':{},'DFT_changes':{},'DFT_center_hartree':d0}
            g_parts[metal]['DFT_scaled_b']=(H*gd).tolist();gradients[key]=gradient['artifacts']
            for point in c['grids']:
                logical=f'{name}_{point}_{metal}';full_id=f'full__{c["global_id"]}_{point}_{metal}'
                center_full=f'full__{c["global_id"]}_center_{metal}'
                values[metal]['C_core'][point]=(mr[logical]['energy_eV']-mr[key]['energy_eV'])*EV_TO_KCAL+(br[logical]['GB_reaction_kcal_mol']-br[key]['GB_reaction_kcal_mol'])
                values[metal]['J'][point]=((sc['rows'][full_id]['energy_eV']-sc['rows'][center_full]['energy_eV'])-
                    (mr[logical]['energy_components_eV']['interaction_energy']-mr[key]['energy_components_eV']['interaction_energy']))*EV_TO_KCAL
                if logical in p['reused_GGR']:
                    actual=old_dft['energies'][p['reused_GGR'][logical]['source_task_id']]['energy_hartree']
                elif logical in dc['rows']:actual=dc['rows'][logical]['energy_hartree']
                else:continue
                values[metal]['DFT_changes'][point]=(actual-d0)*HA_TO_KCAL
            actual_dft[key]=d0
        case_bridge['R']={'difference_kcal_mol':case_bridge['Ca']['difference_kcal_mol']-case_bridge['La']['difference_kcal_mol']}
        case_bridge['R']['pass']=abs(case_bridge['R']['difference_kcal_mol'])<=.05
        normal_bridge[name]=case_bridge
        matrices={};convergence={};local_checks={};short_checks={};stationary={}
        for metal in ('La','Ca','R'):
            if metal=='R':
                values[metal]={term:{pt:values['Ca'][term][pt]-v for pt,v in values['La'][term].items()}
                               for term in ('C_core','J','DFT_changes')}
                g_parts[metal]={k:(np.array(g_parts['Ca'][k])-g_parts['La'][k]).tolist() for k in g_parts['La']}
            matrices[metal]={};convergence[metal]={}
            for component in ('C_core','J'):
                matrices[metal][component]={'coarse':matrix(values[metal][component]),'fine':matrix(values[metal][component],True)}
            matrices[metal]['combined']={scale:matrices[metal]['C_core'][scale]+matrices[metal]['J'][scale] for scale in ('coarse','fine')}
            for component,pair in matrices[metal].items():
                check=square_extrema(pair['fine']-pair['coarse']);check['pass']=check['max_absolute_kcal_mol']<=.01
                convergence[metal][component]=check
            local_checks[metal]=[];short_checks[metal]=[]
            for plus,minus,z in [('s_plus','s_minus',np.array([1.,0.])),('theta_plus','theta_minus',np.array([0.,1.])),('same_plus','same_minus',np.ones(2))]:
                de=values[metal]['DFT_changes'];observed_even=(de[plus]+de[minus])/2
                predicted_even=float(.5*z@matrices[metal]['C_core']['fine']@z)
                err=predicted_even-observed_even;even_tol=max(.005,.25*abs(observed_even));energy_tol=max(.02,.25*abs(observed_even))
                b=np.array(g_parts[metal]['DFT_scaled_b'])
                errors={'plus':float(b@z)+predicted_even-de[plus],'minus':-float(b@z)+predicted_even-de[minus]}
                local_checks[metal].append({'direction':plus.removesuffix('_plus'),'DFT_even_kcal_mol':observed_even,
                    'predicted_even_kcal_mol':predicted_even,'even_error_kcal_mol':err,'even_tolerance_kcal_mol':even_tol,
                    'anchored_errors_kcal_mol':errors,'anchored_tolerance_kcal_mol':energy_tol,
                    'pass':abs(err)<=even_tol and max(map(abs,errors.values()))<=energy_tol})
                odd=(values[metal]['J'][plus]-values[metal]['J'][minus])/2
                predicted=float((H*np.array(g_parts[metal]['J_scaffold']))@z)
                tolerance=max(.02,.05*abs(predicted))
                short_checks[metal].append({'direction':plus.removesuffix('_plus'),'odd_error_kcal_mol':odd-predicted,
                                            'tolerance_kcal_mol':tolerance,'pass':abs(odd-predicted)<=tolerance})
        gates={'normal_SCF_bridge':all(r['pass'] for r in case_bridge.values()),
               'local_DFT_curvature':all(r['pass'] for rows0 in local_checks.values() for r in rows0),
               'short_analytic_gradient':all(r['pass'] for rows0 in short_checks.values() for r in rows0),
               'grid_refinement':all(r['pass'] for components in convergence.values() for r in components.values())}
        for metal in ('La','Ca'):
            combo=matrices[metal]['combined'];stationary[metal]=stationary_point(combo['fine'],b_vectors[metal],combo['coarse'])
            if not all(gates.values()):stationary[metal]['status']='prerequisite_gate_failed'
        rows[name]={'preparation':ref,'evidence':c['evidence'],'DFT_R_kcal_mol':(actual_dft[f'{name}_center_Ca']-actual_dft[f'{name}_center_La'])*HA_TO_KCAL,
                    'gradient_components':g_parts,'scaled_gradient_b':{k:v.tolist() for k,v in b_vectors.items()},
                    'matrices_kcal_mol':{metal:{component:{scale:a.tolist() for scale,a in pair.items()} for component,pair in comps.items()} for metal,comps in matrices.items()},
                    'grid_energies_relative_kcal_mol':values,'convergence':convergence,'local_DFT_checks':local_checks,
                    'short_gradient_checks':short_checks,'gates':gates,'stationary_points':stationary,**UNAVAILABLE}
    sources={k:record(v) for k,v in [('prepared',prepared),('core',core),('gb',gb),('short',short),('DFT',dft),('normal_manifest',normal_manifest),('numerical_convention',convention)]}
    result={'status':'complete','sources':sources,'implementation':record(__file__),'rows':rows,
            'normal_SCF_bridges':normal_bridge,'gradient_sources':gradients,
            'minimum_validation_status':'not_run','predicted_minima_are_not_validated_scores':True,**UNAVAILABLE}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    lines=['# Coupled mechanical-response assessment','',
           'Scores remain unavailable until eligible predicted minima pass independent DFT validation.',
           'Both alpha structures remain one biological group; all cases are consumed development.','',
           '| Case | All prerequisite gates | La stationary status | Ca stationary status |',
           '|---|---|---|---|']
    for case,r in rows.items():lines.append(f'| {case} | {all(r["gates"].values())} | {r["stationary_points"]["La"]["status"]} | {r["stationary_points"]["Ca"]["status"]} |')
    lines.extend(['','All matrices, negative eigenvalues, gradients, residuals and sources are in result.json.',
                  'No clipping, inherited threshold, zero substitution or baseline promotion.'])
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('prepared','core','gb','short','dft','normal-manifest','convention','output'):p.add_argument('--'+key,required=True)
    a=p.parse_args();result=assess(**vars(a));print(json.dumps({k:v['gates'] for k,v in result['rows'].items()},indent=2))
