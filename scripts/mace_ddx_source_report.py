"""Audit actual ddX coefficients and predeclared component-convergence checks."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
import shutil
import sys
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,read_json,record,verify,write_new
from mace_ddx_source_self import CASES,GRIDS,TOL,validate,potential_error


def report(collection,output):
    start=time.monotonic();r=read_json(collection);mp=verify(r['manifest'])
    recovery=r['protocol'].endswith('_independent_iteration_recovery_v3')
    conductor=r['protocol'] in ('fixed_source_full_protein_ddCPCM_component_v1','fixed_source_full_protein_ddCPCM_component_refinement_v2')
    if conductor:
        from mace_ddx_cpcm import preflight
        from mace_ddx_recovery import retained_state
        m=preflight(mp)
    elif recovery:
        from mace_ddx_recovery import preflight,retained_state
        m=preflight(mp)
    else:m=validate(mp)
    checks=[];rows={};lookup={};arrays={};factor=m.get('energy_prefactor',1.)
    def check(name,error,tolerance,**meta):
        checks.append(dict(name=name,error=error,tolerance=tolerance,pass_=error is not None and abs(error)<=tolerance,**meta))
    for g in m['groups']:
        gid=g['group_id'];result=r['groups'][gid];d=read_json(verify(g['source']));lookup[(g['case_id'],g['variant'],g['grid'])]=gid
        if result['cache_key']!=g['cache_key'] or result['source']!=g['source']:raise InvalidArtifact('native group provenance changed')
        expected_parameters={**m['model'],**m['grids'][g['grid']]}
        for name,expected in expected_parameters.items():
            if name in ('model','enable_force'):continue
            if recovery or conductor:continue # audited per state below, including original reused iteration settings
            observed=result.get('actual_parameters',{}).get(name)
            check('runtime_parameter_'+name,0. if observed==expected else None,0.,group=gid)
        with np.load(verify(result['arrays']),allow_pickle=False) as data:arrays[gid]={k:data[k] for k in data.files}
        if set(result['rows'])!=set(g['states']):
            check('all_declared_states_available',None,0.,group=gid)
        for label in g['states']:
            key=gid+'/'+label;row=result['rows'].get(label,dict(status='unavailable',failure_reason=result.get('failure_reason','state absent')))
            rows[key]=dict(row)
            if recovery or conductor:
                expected={**expected_parameters,'maxiter':g['reuse'][label]['maxiter'] if label in g['reuse'] else 1200}
                for name,value in expected.items():
                    if name in ('model','enable_force'):continue
                    check('runtime_parameter_'+name,0. if row.get('actual_parameters',{}).get(name)==value else None,0.,group=gid,state=label)
                if label in g['reuse']:
                    original,old_arrays=retained_state(g['reuse'][label])
                    check('reuse_receipt',0. if row.get('reused_from')==g['reuse'][label] and row.get('execution')=='reused' and not row.get('new_forward_solve_started') else None,0.,group=gid,state=label)
                    for part in ('phi','psi','x'):
                        check('reuse_'+part,0. if np.array_equal(arrays[gid].get(label+'_'+part),old_arrays[part]) else None,0.,group=gid,state=label)
                    check('reuse_energy_kcal',row.get('energy_kcal_mol')-original['energy_kcal_mol'] if row.get('energy_kcal_mol') is not None else None,0.,group=gid,state=label)
            if row['status']!='computed':
                reason=row.get('failure_reason','')
                rows[key]['execution_failure_class']=('forward_solver_failed' if row.get('forward_solve_started') else
                    'not_started_with_native_error_flag' if 'Jacobi solver' in reason else
                    'source_preflight_failed' if 'source potential/integral gate' in reason else 'unavailable_or_setup_failed')
            if row['status']!='computed':check('state_computed',None,0.,group=gid,state=label);continue
            a=arrays[gid];psi=a[label+'_psi'];x=a[label+'_x'];phi=a[label+'_phi'];q=np.zeros(len(d['physical_ids'])) if label=='zero' else np.array(d['endpoints'][label.split('_')[0]]['charge_e'])
            energy=factor*.5*math.fsum((psi*x).ravel());expected=np.zeros_like(psi);expected[0]=np.sqrt(4*np.pi)*q
            if conductor:
                check('conductor_energy_scaling',row.get('energy_prefactor')-factor if row.get('energy_prefactor') is not None else None,0.,group=gid,state=label)
                raw=row.get('raw_native_energy_hartree')
                check('retained_raw_conductor_energy_kcal',(raw*factor-row['energy_hartree'])*HA_TO_KCAL if raw is not None else None,TOL['energy_kcal'],group=gid,state=label)
            check('energy_contraction_kcal',(energy-row['energy_hartree'])*HA_TO_KCAL,TOL['energy_kcal'],group=gid,state=label)
            check('single_unit_conversion',row['energy_kcal_mol']-row['energy_hartree']*HA_TO_KCAL,1e-12,group=gid,state=label)
            check('source_psi',float(np.max(abs(psi-expected))),TOL['psi'],group=gid,state=label)
            perr=potential_error(a['cavity_bohr'],np.array(d['coordinates_A']).T/BOHR_TO_A,q,phi)
            check('source_potential_au',perr,TOL['potential_au'],group=gid,state=label)
            check('passive_energy_kcal',max(0.,row['energy_kcal_mol']),TOL['energy_kcal'],group=gid,state=label)
            if label=='zero':check('zero_source_kcal',row['energy_kcal_mol'],TOL['identity_kcal'],group=gid)
            if label.endswith('_repeat'):
                original=result['rows'].get(label.split('_')[0],{})
                delta=row['energy_kcal_mol']-original['energy_kcal_mol'] if original.get('status')=='computed' else None
                check('repeat_kcal',delta,TOL['repeat_kcal'],group=gid,state=label)
        if all(result['rows'].get(x,{}).get('status')=='computed' for x in ('Ca','La')):
            a=arrays[gid];err=factor*.5*abs(math.fsum((a['Ca_psi']*a['La_x']).ravel())-math.fsum((a['La_psi']*a['Ca_x']).ravel()))*HA_TO_KCAL
            check('reciprocity_kcal',err,TOL['reciprocity_kcal'],group=gid)
        else:check('reciprocity_kcal',None,TOL['reciprocity_kcal'],group=gid)
        arrays[gid].clear()
    def values(case,variant,grid):
        group=r['groups'][lookup[(case,variant,grid)]];values={}
        for metal in ('Ca','La'):
            row=group['rows'].get(metal,{});values[metal]=row['energy_kcal_mol'] if row.get('status')=='computed' else None
        values['R']=values['Ca']-values['La'] if None not in (values['Ca'],values['La']) else None
        return values
    results={};differences={};coarse_changes={};native_differences={}
    for case in CASES:
        results[case]={grid:values(case,'primary',grid) for grid in GRIDS};primary=results[case]['primary'];refined=results[case]['refined'];rigid=values(case,'rigid','primary');coarse=results[case]['coarse']
        coarse_changes[case]={}
        for term in ('Ca','La','R'):
            check('refinement_kcal',refined[term]-primary[term] if None not in (refined[term],primary[term]) else None,TOL['refinement_kcal'],case=case,term=term)
            check('rigid_kcal',rigid[term]-primary[term] if None not in (rigid[term],primary[term]) else None,TOL['rigid_kcal'],case=case,term=term)
            coarse_changes[case][term]=primary[term]-coarse[term] if None not in (primary[term],coarse[term]) else None
        source=read_json(verify(next(g['source'] for g in m['groups'] if g['case_id']==case and g['variant']=='primary')))
        native={metal:e['source_self_GK_primary_kcal'] for metal,e in source['endpoints'].items()};native['R']=native['Ca']-native['La']
        native_differences[case]={'native_GK':native,('CPCM_primary_minus_GK' if conductor else 'PCM_primary_minus_GK'):{k:primary[k]-native[k] if primary[k] is not None else None for k in native}}
    left,right=CASES
    for grid in GRIDS:
        differences[grid]={term:results[left][grid][term]-results[right][grid][term] if None not in (results[left][grid][term],results[right][grid][term]) else None for term in ('Ca','La','R')}
    for term in ('Ca','La','R'):
        primary=differences['primary'][term];refined=differences['refined'][term]
        check('between_structure_refinement_kcal',refined-primary if None not in (primary,refined) else None,TOL['refinement_kcal'],term=term)
        a,b=values(left,'rigid','primary')[term],values(right,'rigid','primary')[term]
        check('between_structure_rigid_kcal',(a-b)-primary if None not in (a,b,primary) else None,TOL['rigid_kcal'],term=term)
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,root/Path(__file__).name)
    out=dict(protocol=r['protocol']+'_numerical_report',collection=record(collection),manifest=r['manifest'],implementation=record(root/Path(__file__).name),
             dependencies={n:record(sys.modules[n].__file__) for n in ('affordable_common','mace_ddx_source_self')},
             checks=checks,checks_pass=all(c['pass_'] for c in checks),complete=all(row['status']=='computed' for row in rows.values()),
             source_self_kcal=results,GGR_2FW0_minus_2FVY_kcal=differences,coarse_to_primary_changes_kcal=coarse_changes,native_GK_comparison=native_differences,
             failed_states={k:v for k,v in rows.items() if v['status']!='computed'},new_solver_calls=0,new_score=None,biological_classification=None,baseline_changed=False,
             scope='source-only fixed-charge reaction energies on a full physical protein cavity; no full hybrid or exact-density calculation',wall_seconds=time.monotonic()-start)
    if recovery:out['dependencies']['mace_ddx_recovery']=record(sys.modules['mace_ddx_recovery'].__file__)
    if conductor:
        out['dependencies']['mace_ddx_cpcm']=record(sys.modules['mace_ddx_cpcm'].__file__)
        out['energy_prefactor']=factor
    write_new(root/'result.json',out);return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--collection',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=report(a.collection,a.output)
    print({k:r[k] for k in ('complete','checks_pass','source_self_kcal','GGR_2FW0_minus_2FVY_kcal','wall_seconds')})
