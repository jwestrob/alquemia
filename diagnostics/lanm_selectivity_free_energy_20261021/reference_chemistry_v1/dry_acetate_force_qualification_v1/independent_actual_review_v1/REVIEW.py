#!/usr/bin/env python3
"""Independent saved-token/gradient arithmetic and exact native lineage review."""
import hashlib, importlib.util, json, math, re, time
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
WS=Path('workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_force_qualification_v1').resolve()
SOURCE=HERE.parent/'REFERENCE.py'
spec=importlib.util.spec_from_file_location('dry_native_admission',SOURCE)
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
F=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
HA=627.509474;BOHR=.529177210903

def read(p):return json.loads(Path(p).read_text())
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verified(p):
    assert pin(p['path'])==p
    return Path(p['path'])
def close(a,b,tol=1e-10):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    assert a.shape==b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b))
    err=float(np.max(np.abs(a-b))) if a.size else 0.
    assert err<=tol,(err,tol)
    return err
def energy_tokens(raw):
    native=raw.split('****END OF INPUT****')[-1]
    patterns=dict(electronic=r'^Total Energy\s*:\s*('+F+r')\s+Eh',D4=r'^Dispersion correction\s+('+F+r')',PBE0_D4=r'^FINAL SINGLE POINT ENERGY\s+('+F+r')')
    result={}
    for key,pattern in patterns.items():
        matches=list(re.finditer(pattern,native,re.M));assert len(matches)==1,(key,len(matches))
        token=matches[0].group(1);value=float(token.replace('D','E'));assert math.isfinite(value)
        result[key]=dict(token=token,value_hartree=value)
    assert abs(result['electronic']['value_hartree']+result['D4']['value_hartree']-result['PBE0_D4']['value_hartree'])<1e-8
    return result
def gradient_file(path):
    lines=[s.strip() for s in Path(path).read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    n=int(lines[0]);assert n==8 and len(lines)==2+3*n+n
    gradient=np.array([float(v) for v in lines[2:2+3*n]]).reshape(n,3)
    xyz=np.array([[float(v) for v in line.split()[1:4]] for line in lines[2+3*n:]])*BOHR
    numbers=[int(line.split()[0]) for line in lines[2+3*n:]]
    return float(lines[1]),gradient,xyz,numbers

def main():
    tick,cpu=time.monotonic(),time.process_time();out=WS/'independent_actual_review_v1';out.mkdir(exist_ok=False)
    result_path=WS/'remaining_submission_v1/collection_v1/RESULT.json';result=read(result_path)
    assert result['status']=='dry_acetate_force_complete' and result['adapter']==pin(SOURCE) and result['admitted_cells']==16
    panel_path=verified(result['panel']);p=h.load(panel_path);p['_path']=str(panel_path)
    assert h.b.HA==HA and h.b.BOHR==BOHR
    inventory=read(verified(result['inventory']));expected={f'pbe0__{m}__{d}__{s}__complex' for m in ['La','Dy'] for d in ['radial','rotation'] for s in ['minus_010','minus_005','plus_005','plus_010']}
    assert set(inventory)==set(result['slots'])==expected and len(set(inventory.values()))==16
    assert not result['failures'] and result['new_native_calls']==0
    centers={};frames={};analytic={};central_geometries={}
    for metal,center in p['central'].items():
        saved=read(verified(center['saved_result']));st=read(verified(center['saved_start']));g=center['original_cell']['geometry']
        assert saved['status']=='complete' and all(saved['checks'].values())
        assert st['cell']==center['original_cell'] and st['panel']==center['original_panel'] and st['input']==center['executed_input']
        parent=h.b.load_panel(verified(center['original_panel']))
        assert verified(center['executed_input']).read_text()==h.b.deck(parent,center['original_cell'],st['ranks'],st['maxcore_MiB'])
        raw=verified(center['output']).read_text();tokens=energy_tokens(raw)
        e,grad,gxyz,z=gradient_file(verified(center['engrad']));x=np.array(g['positions_A']);force=-grad*HA/BOHR
        close(e,tokens['PBE0_D4']['value_hartree'],1e-10);close(force,saved['forces_kcal_mol_A'],1e-12);close(gxyz,x,1e-7)
        assert z==[57 if metal=='La' else 66,6,6,8,8,1,1,1]
        u=x[3]-x[0];u=u/np.linalg.norm(u);centroid=x[1:].mean(axis=0);axis=np.cross(u,x[4]-x[3]);axis=axis/np.linalg.norm(axis)
        net=force[1:].sum(axis=0);torque=np.cross(x[1:]-centroid,force[1:]).sum(axis=0)
        projection=dict(radial=float(-(net@u)),rotation=float(-(torque@axis)))
        for d,v in projection.items():close(v,center['analytic'][d])
        strict=h.w.parse(raw,h.h.metadata(g),saved['returncode']);strict['dft_xc_grid_metadata']=h.xc_metadata(raw)
        assert strict==center['strict_native_replay'] and strict['status']=='complete'
        frames[metal]=(u,centroid,axis);analytic[metal]=projection;central_geometries[metal]=x
        centers[metal]=dict(result=center['saved_result'],output=center['output'],engrad=center['engrad'],input=center['executed_input'],
            independent_energy_tokens=tokens,gradient_forces_max_replay_error=close(force,saved['forces_kcal_mol_A']),analytic=projection,
            net_force_kcal_mol_A=force.sum(axis=0).tolist(),strict_native=strict,center_rerun=False)
    assert np.array_equal(central_geometries['La'],central_geometries['Dy'])
    native={};energies={};jobs=[];geometry_errors=[];native_cost=[]
    for c in p['cells']:
        cid=c['cell_id'];path=Path(inventory[cid]);assert result['slots'][cid]['result']==pin(path)
        r=h.admit(p,result['panel'],path);raw=verified(r['output']).read_text();tokens=energy_tokens(raw)
        for label,field in [('electronic','electronic_energy_hartree'),('D4','D4_energy_hartree'),('PBE0_D4','energy_hartree')]:
            assert tokens[label]['value_hartree']==r[field]
        x=central_geometries[c['metal']];u,centroid,axis=frames[c['metal']];q=c['q_A'];assert q in [-.01,-.005,.005,.01]
        expected_xyz=x.copy()
        if c['direction']=='radial':expected_xyz[1:]+=q*u
        else:
            v=x[1:]-centroid
            expected_xyz[1:]=centroid+v*np.cos(q)+np.cross(axis,v)*np.sin(q)+np.outer(v@axis,axis)*(1-np.cos(q))
        actual=np.array(c['positions_A']);geomerr=close(actual,expected_xyz,1e-12);geometry_errors.append(geomerr)
        assert np.array_equal(actual[0],x[0]) and c['charge']==2 and c['multiplicity']==1 and c['expected_electrons']==40
        assert c['physical_multiplicity_of_parent']==(1 if c['metal']=='La' else 6) and c['expected_ecp_core']==(46 if c['metal']=='La' else 55)
        assert c['expected_basis_functions']==189 and not any(c['ghost_mask'])
        close(np.linalg.norm(actual[1:,None,:]-actual[None,1:,:],axis=2),np.linalg.norm(x[1:,None,:]-x[None,1:,:],axis=2),1e-9)
        original=p['central'][c['metal']];resources=original['original_resource_values']
        assert h.deck(p,c,resources['ranks'],resources['maxcore_MiB'])==verified(original['executed_input']).read_text().replace(' EnGrad','')
        start=read(path.parent/'START.json');assert start['job_id']==r['job_id'] and start['command']==[p['orca']['path'],'endpoint.inp','--bind-to none']
        assert r['new_native_calls']==1 and r['returncode']==0 and not r['native_energy_stages']['has_analytic_gradient_stage']
        assert Path(path.parent/'endpoint.engrad').exists() is False
        verified(r['stderr']);jobs.append(r['job_id']);energies[cid]={k:v['value_hartree'] for k,v in tokens.items()}
        native_cost.append(dict(cell_id=cid,job_id=r['job_id'],elapsed_seconds=r['elapsed_seconds'],allocated_cpu_seconds=r['allocated_cpu_seconds']))
        prior=next(z for z in result['native_reviews'] if z['cell_id']==cid)
        assert prior['printed_scf_criteria']==r['printed_scf_criteria'] and prior['energy_stages']==r['native_energy_stages'] and prior['dft_xc_grid_metadata']==r['dft_xc_grid_metadata']
        assert prior['center_stage_comparison']==h.w.compare_stages(original['strict_native_replay']['native_energy_stages'],r['native_energy_stages'])
        native[cid]=dict(result=pin(path),output=r['output'],input=r['input'],geometry=r['geometry'],start=pin(path.parent/'START.json'),release=r['release'],
            independent_energy_tokens=tokens,transform_max_error_A=geomerr,all_native_checks_passed=True,native_checks=len(r['checks']),
            native_caveats=prior)
    assert len(set(jobs))==16
    scout=read(WS/'scout_submission_v1/collection_v1/RESULT.json')
    assert all(result['slots'][cid]==scout['slots'][cid] for cid in h.SCOUT)
    # Independent finite differences and torque projection; do not call the producer's score().
    checks={};table={}
    for d in ['radial','rotation']:
        channels={}
        for metal in ['La','Dy']:
            fd={}
            for step,suffix in [(.01,'010'),(.005,'005')]:
                plus=energies[f'pbe0__{metal}__{d}__plus_{suffix}__complex'];minus=energies[f'pbe0__{metal}__{d}__minus_{suffix}__complex']
                fd[str(step)]={k:(plus[k]-minus[k])*HA/(2*step) for k in plus}
            channels[metal]=dict(analytic=analytic[metal][d],derivatives=fd)
        channels['paired']=dict(analytic=channels['Dy']['analytic']-channels['La']['analytic'],
            derivatives={step:{k:channels['Dy']['derivatives'][step][k]-channels['La']['derivatives'][step][k] for k in ['electronic','D4','PBE0_D4']} for step in ['0.01','0.005']})
        for metal,values in channels.items():
            saved=result['scoring']['directions'][d][metal];close(values['analytic'],saved['analytic_kcal_mol_A'])
            errors={s:values['derivatives'][s]['PBE0_D4']-values['analytic'] for s in ['0.01','0.005']}
            delta={k:values['derivatives']['0.005'][k]-values['derivatives']['0.01'][k] for k in ['electronic','D4','PBE0_D4']}
            for s in errors:
                for k in delta:close(values['derivatives'][s][k],saved['derivatives_kcal_mol_A'][s][k])
                close(errors[s],saved['total_FD_minus_analytic'][s])
            for k in delta:close(delta[k],saved['smaller_minus_larger_step'][k])
            for name,error in [('two_step',delta['PBE0_D4']),('analytic_0.01',errors['0.01']),('analytic_0.005',errors['0.005'])]:
                key=f'{d}/{metal}/{name}';passed=abs(error)<=.1
                close(error,result['scoring']['primary_checks'][key]['error_kcal_mol_A']);assert result['scoring']['primary_checks'][key]['passed']==passed
                checks[key]=dict(error_kcal_mol_A=error,passed=passed,resolution_kcal_mol_A=.1)
            values.update(errors=errors,two_step_difference=delta)
        table[d]=channels
    assert len(checks)==18 and all(v['passed'] for v in checks.values()) and result['stage_numerical_checks_passed']
    receipts={};scheduler={};support={}
    for stage in ['scout_submission_v1','remaining_submission_v1']:
        for suffix,is_native in [('completion_event.json',True),('collector_submission/completion_event.json',False)]:
            path=WS/stage/suffix;receipt=read(path);receipts[str(path)]=pin(path)
            for job,rec in receipt['records'].items():
                assert rec['state']=='COMPLETED' and rec['exit_code']=='0:0'
                (scheduler if is_native else support)[job]=rec
    assert set(scheduler)==set(jobs) and all(int(x['allocated_cpus'])==4 for x in scheduler.values())
    native_cpu=sum(v['allocated_cpu_seconds'] for v in native_cost);close(native_cpu,result['native_allocated_cpu_seconds_sum'])
    maxerr=max(abs(v['error_kcal_mol_A']) for k,v in checks.items() if '/analytic_' in k)
    maxstep=max(abs(v['error_kcal_mol_A']) for k,v in checks.items() if k.endswith('/two_step'))
    q=dict(status='passed_dry_acetate_force_actual_review',implementation=pin(__file__),reviewed_implementation=pin(SOURCE),result=pin(result_path),panel=result['panel'],
        complete=16,primary_checks=checks,all_18_numerical_checks_passed=True,independent_table=table,
        maximum_derivative_discrepancy_kcal_mol_A=maxerr,maximum_two_step_difference_kcal_mol_A=maxstep,
        central_raw_replays=centers,displaced_raw_replays=native,unique_native_jobs=jobs,exact_scout_reuse=4,new_center_calls=0,
        independent_geometry_max_error_A=max(geometry_errors),raw_energy_gradient_coordinates_and_state_verified=True,
        historical_centers_reused=True,physical_state_qualified=False,method_accuracy_qualified=False,
        caveats=['Only two prescribed directions on two consumed exact dry mono-acetate centers are qualified at 0.1 kcal/mol/Angstrom.','SCF Energy Check stopping and failed printed density/DIIS criteria remain explicit.','COSX stage signatures and DFT/XC settings match, while coordinate-dependent grid counts and EnGrad-to-SP distinction do not establish exact numerical-surface identity.','No proof of other donor forces, unique electronic state, basis/method accuracy or physical preference follows.'],
        cost=dict(native_records=native_cost,native_wall_seconds_sum=sum(x['elapsed_seconds'] for x in native_cost),native_allocated_CPU_seconds=native_cpu,
            scheduler_records=scheduler,scheduler_native_CPU_seconds=sum(int(x['elapsed_seconds'])*int(x['allocated_cpus']) for x in scheduler.values()),
            collector_scheduler_records=support,collector_scheduler_CPU_seconds=sum(int(x['elapsed_seconds'])*int(x['allocated_cpus']) for x in support.values()),
            saved_scheduler_receipts=receipts,original_collection_measured_cost=result['cost'],review_wall_seconds=time.monotonic()-tick,review_CPU_seconds=time.process_time()-cpu,threads=1,new_native_calls=0))
    dest=out/'REVIEW.json'
    with dest.open('x') as f:json.dump(q,f,indent=2,allow_nan=False);f.write('\n')
    lines=['# Independent actual dry mono-acetate derivative review','',
        '**Passed: all 16 displaced outputs and 18 prescribed numerical checks.** Both original central gradients replay from `.engrad`; all four scout outputs are reused exactly. No center rerun or omitted cell.','',
        '| Direction | Comparison | Analytic | FD 0.010 | FD 0.005 |','|---|---|---:|---:|---:|']
    for direction,channels in table.items():
        for metal,v in channels.items():lines.append(f"| {direction} | {metal} | {v['analytic']:.8f} | {v['derivatives']['0.01']['PBE0_D4']:.8f} | {v['derivatives']['0.005']['PBE0_D4']:.8f} |")
    lines+=['',f'Units are kcal/mol/Angstrom, including the fixed 1 Angstrom rotational coordinate. Maximum FD/analytic discrepancy is **{maxerr:.8f}**; maximum two-step difference is **{maxstep:.8f}**, both below the unchanged 0.1 resolution. Independent raw final/electronic/D4 tokens, gradient signs, torque projections, geometry transformations, input/state lineage and all component derivatives agree with the saved result.','',
        'Every displaced calculation uses Energy Check stopping. Printed density/DIIS failures remain; COSX stage signatures and DFT/XC settings agree with centers, but coordinate-dependent grid counts and central EnGrad versus displaced SP remain explicit. This qualifies only the sampled directions and resolution, not physical reference accuracy, other donor forces or a unique electronic state.','',
        f"Sixteen native calls: {q['cost']['native_wall_seconds_sum']:.6f} summed wall seconds and **{native_cpu:.6f} native allocated CPU-s**. Saved scheduler totals: **{q['cost']['scheduler_native_CPU_seconds']} CPU-s**, plus **{q['cost']['collector_scheduler_CPU_seconds']} CPU-s** for the two collectors. Historical center costs are not charged again. This review made zero native calls and used one local thread."]
    for folder in (out,HERE):
        with (folder/'REPORT.md').open('x') as f:f.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(review=pin(dest),report=pin(HERE/'REPORT.md'),max_derivative_error=maxerr,max_step_error=maxstep,cost={k:v for k,v in q['cost'].items() if k not in ['native_records','scheduler_records','collector_scheduler_records','saved_scheduler_receipts']})))

if __name__=='__main__':main()
