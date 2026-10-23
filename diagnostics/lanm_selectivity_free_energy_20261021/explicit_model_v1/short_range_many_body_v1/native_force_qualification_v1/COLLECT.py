#!/usr/bin/env python3
"""Saved-only finite directional force comparison; never dispatches native work."""
import argparse,json,math,platform,time
from pathlib import Path
from PREPARE import pin,read,verify,save

def slope(plus,minus,h,gradient,direction):
    analytic=math.fsum(g*d for ga,da in zip(gradient,direction) for g,d in zip(ga,da));fd=(plus-minus)/(2*h);error=fd-analytic;limit=5e-5+1e-6*abs(analytic)
    assert all(math.isfinite(x) for x in (analytic,fd,error,limit))
    return {'analytic_kcal_mol_A':analytic,'finite_difference_kcal_mol_A':fd,'signed_error_kcal_mol_A':error,'absolute_error_kcal_mol_A':abs(error),'limit_kcal_mol_A':limit,'passed':abs(error)<=limit}

def collect(scout_pin,remaining_pin,out):
    t=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False);scout=read(scout_pin);panel=read(scout['panel']);full=remaining_pin is not None
    result={'status':'started','implementation':pin(__file__),'scout_result':scout_pin,'remaining_result':remaining_pin,'panel':scout['panel'],'mode':'full' if full else 'scout','checks':{},'error':None,'new_native_calls':0}
    try:
        assert scout['status']=='native_force_scout_complete' and scout['mode']=='scout' and scout['native_program_invocations']==scout['native_program_attempts']==6 and set(scout['rows'])==set(panel['scout_ids'])
        maps=dict(scout['rows']);runs=[scout]
        if full:
            remaining=read(remaining_pin);assert remaining['status']=='native_force_remaining_complete' and remaining['mode']=='remaining' and remaining['native_program_invocations']==remaining['native_program_attempts']==96
            for key in ('implementation','panel','native_binary','kernel','combiner','build'):assert remaining[key]==scout[key]
            assert remaining['reused_scout_rows']==scout['rows'] and set(remaining['rows'])==set(panel['remaining_ids']) and set(maps).isdisjoint(remaining['rows']);maps.update(remaining['rows']);runs.append(remaining)
        expected=set(panel['states']) if full else set(panel['scout_ids']);assert set(maps)==expected
        rows={};regressions={};calls={'direct_initializations':0,'mutual_operator_applications':0,'baseline_gradient_calls':0,'induced_solves':0}
        for id,p in maps.items():
            r=read(p);assert r['status']=='native_force_state_complete' and r['error'] is None and r['returncode']==0 and r['native_program_attempts']==r['native_program_invocations']==1 and r['state']==panel['states'][id]
            assert r['kernel']==scout['kernel'] and r['native_binary']==scout['native_binary'] and r['combination_source']==scout['combiner']
            for key in ('baseline_state','supplied_uv','force_bridge_result'):verify(r[key])
            for pinput in r['state']['inputs'].values():verify(pinput)
            assert math.isfinite(r['energy_kcal_mol']) and r['minimum_corrected_eigenvalue']>0 and math.isfinite(r['solve_residual'])
            assert all(math.isfinite(x) for a in r['gradient_kcal_mol_A'] for x in a)
            assert r['raw_native']['supplied_uv_unchanged'] and r['counts']=={'direct_initializations':1,'mutual_operator_applications':3*len(r['state']['geometry']['xyz_A']),'baseline_gradient_calls':1,'induced_solves':0}
            for k in calls:calls[k]+=r['counts'][k]
            if r['zero_regression'] is not None:regressions[id]=r['zero_regression'];assert r['zero_regression']['passed']
            rows[id]=r
        for id,r in rows.items():
            s=r['state']
            if s['sign']!=1:continue
            minus_id=id[:-4]+'minus';center_id=s['case_id']+'__'+s['coefficient_mode']+'__center';assert minus_id in rows and center_id in rows
            check=slope(r['energy_kcal_mol'],rows[minus_id]['energy_kcal_mol'],s['step_A'],rows[center_id]['gradient_kcal_mol_A'],panel['directions'][s['case_id']][s['direction']])
            check.update(case=s['case_id'],coefficient_mode=s['coefficient_mode'],direction=s['direction'],step_tag=s['step_tag'],step_A=s['step_A'],plus=maps[id],minus=maps[minus_id],center=maps[center_id]);result['checks'][id[:-6]]=check
        assert len(result['checks'])==(48 if full else 2)
        result.update(status=('passed_native_force_full_panel' if full else 'passed_native_force_scout') if all(r['passed'] for r in result['checks'].values()) else 'failed_native_force_directional_gate',native_states=len(rows),native_programs_observed=sum(r['native_program_invocations'] for r in runs),native_calls_by_kind=calls,zero_regressions=regressions,source_runs=[{'implementation':r['implementation'],'configuration':r['configuration'],'root_release':r['root_release'],'cost':r['cost']} for r in runs],source_cost_totals={'wall_seconds_serial_runs':sum(r['cost']['wall_seconds'] for r in runs),'controller_CPU_seconds':sum(r['cost']['controller_CPU_seconds'] for r in runs),'native_child_CPU_seconds':sum(r['cost']['native_child_CPU_seconds'] for r in runs)},limitations=['Finite vacuum derivative implementation only; no physical accuracy, periodic, stability, fit or exchange-free-energy validation.','Missing archived zero-center water full gradients remain explicitly unavailable.','Raw native virial lacks the added local correction and is not a qualified total virial.'])
    except Exception as exc:result['status']='failed_saved_collection';result['error']=repr(exc)
    result['collection_cost']={'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1,'host':platform.node()};save(out/'RESULT.json',result);return result

def main():
    a=argparse.ArgumentParser();a.add_argument('--scout',type=Path,required=True);a.add_argument('--remaining',type=Path);a.add_argument('--out',type=Path,required=True);c=a.parse_args();r=collect(pin(c.scout),pin(c.remaining) if c.remaining else None,c.out);print(json.dumps({'status':r['status'],'checks':len(r['checks']),'error':r['error']}));return 0 if r['status'].startswith('passed_') else 1
if __name__=='__main__':raise SystemExit(main())
