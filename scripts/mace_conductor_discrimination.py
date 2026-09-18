"""Complete frozen-response ddCPCM versus GK discrimination on archived states."""
from pathlib import Path
import argparse
import fcntl
import gc
import json
import math
import os
import resource
import shutil
import subprocess
import time
import numpy as np
from affordable_common import BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new
from mace_frozen_response import source_arrays, parse as parse_replay
from mace_ddx_multipole_bridge import distributions, coefficients, cartesian_potential

PROTOCOL = 'trial_density_MACE_AMOEBA_frozen_response_ddCPCM_v1'
CASES = ('GGR_extended', 'GGR_2FW0', 'GGR_2FVY', 'ALPHA_1F6S', 'ALPHA_6IP9', 'PARV_4CPV_CD', 'PARV_4CPV_EF')
GRIDS = dict(primary=dict(lmax=18,n_lebedev=974), refined=dict(lmax=24,n_lebedev=2030))
TOL = dict(source_au=1e-10, psi=1e-12, stock_array=1e-12, solver=1e-10,
           contraction_kcal=1e-8, refinement_kcal=.1, direction_kcal=.02)


def corrected_contrast(old, cpcm, gk):
    return old + (cpcm['Ca']-gk['Ca']) - (cpcm['La']-gk['La'])


def prepare(config, output):
    c=read_json(config); parent=read_json(verify(c['conductor'])); build=read_json(verify(c['software']))
    sw=read_json(verify(c['native_software'])); functional=read_json(verify(c['functional']))
    if build['status']!='built_and_imported' or sw['returncode'] or not functional['checks_pass']:
        raise InvalidArtifact('completed software and native functional required')
    env=Path(c['software']['path']).parent/'venv'; modules=list((env/'lib/python3.11/site-packages').glob('pyddx*.so'))
    if len(modules)!=1: raise InvalidArtifact('isolated module unavailable')
    root=Path(output).resolve(); root.mkdir(parents=True,exist_ok=False); impl=root/'implementation'; impl.mkdir()
    # Reuse the existing dependency closure, with only this driver and its helpers updated.
    bm=read_json(verify(c['bridge_manifest']))
    for name,pin in bm['implementation'].items(): shutil.copyfile(verify(pin),impl/name)
    for name in ('mace_conductor_discrimination.py','mace_ddx_multipole_bridge.py','mace_frozen_response.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    bridge=read_json(verify(c['bridge_collection']))
    reference=bridge['groups']['GGR_2FW0_primary']
    if not all(x['checks_pass'] for x in reference['states'].values()): raise InvalidArtifact('real stock source checks failed')
    py=env/'bin/python'
    m=dict(protocol=PROTOCOL,config=record(config),plan=c['plan'],software=c['software'],native_software=c['native_software'],
           functional=c['functional'],conductor=c['conductor'],python={**record(py),'path':str(py)},module=record(modules[0]),
           model=parent['model'],energy_prefactor=parent['energy_prefactor'],grids=GRIDS,tolerances=TOL,
           stock_fixture=dict(arrays=reference['arrays'],physical=next(g['physical'] for g in bm['groups'] if g['group_id']=='GGR_2FW0_primary'),grid=bm['grid']),
           implementation={p.name:record(p) for p in impl.glob('*.py')},groups=[],requested_forward_solves=56,
           new_DFT_calls=0,new_MACE_calls=0,baseline_changed=False,reference=None,calibrated_class=None)
    fm=read_json(verify(functional['manifest']))
    for case in CASES:
        co=read_json(verify(c['collections'][case])); pm=read_json(verify(co['manifest']))
        if not co['complete'] or not co['numerical_pass']: raise InvalidArtifact('incomplete native archive')
        if read_json(verify(pm['software']))['library']!=sw['library']: raise InvalidArtifact('native library differs')
        old=co['variants']['primary']['cases'][case]; a,b=(old['endpoints'][x] for x in ('Ca','La'))
        if any(a[k]!=b[k] for k in ('induction_environment_only_kcal','native_environment_static_components')):
            raise InvalidArtifact('environment-only reference does not cancel')
        g=dict(case_id=case,collection=c['collections'][case],old=old,endpoints={})
        for metal in ('Ca','La'):
            st=next(t for t in pm['static_tasks'] if t['case_id']==case and t['state']==metal and t['variant']=='primary')
            rt=next(t for t in pm['response_tasks'] if t['case_id']==case and t['state']==metal and t['variant']=='primary' and not t['task_id'].endswith('_standard'))
            sp=co['tasks'][st['task_id']]['result']; rp=co['tasks'][rt['task_id']]['result']
            s=read_json(verify(sp)); r=read_json(verify(rp)); arrays=source_arrays(s,r)
            frozen=st['frozen_indices']; ix=np.array(frozen)-1
            if np.any(arrays[1][ix]) or np.any(arrays[2][ix]): raise InvalidArtifact('source dipoles not frozen')
            physical=dict(coordinates_A=[x['xyz_A'] for x in s['parameters']['atoms']],radii_A=[x['radius_A'] for x in s['parameters']['atoms']],
                          physical_ids=read_json(verify(st['boundary']))['physical_ids'],frozen_indices=frozen)
            if 'physical' in g and physical!=g['physical']: raise InvalidArtifact('paired physical cavity/mapping differs')
            g['physical']=physical
            permanent=np.array([x['global_'] for x in s['moments']]); permanent[ix]=0
            if metal=='Ca': environment=permanent
            elif not np.array_equal(environment,permanent): raise InvalidArtifact('paired permanent environment differs')
            d=root/'inputs'/(case+'_'+metal);d.mkdir(parents=True)
            t=dict(task_id=d.name,static=sp,response=rp,frozen_indices=frozen)
            for name in ('xyz','key','mask','overrides'):
                target=d/Path(st[name]['path']).name;shutil.copyfile(verify(st[name]),target);t[name]=record(target)
            rows=np.column_stack((np.arange(1,len(arrays[0])+1),*arrays)); target=d/'frozen_response.dat'
            target.write_text(str(len(rows))+'\n'+''.join(str(int(row[0]))+' '+' '.join(format(v,'.17g') for v in row[1:])+'\n' for row in rows))
            t['input']=record(target);t['reuse']=None
            for ft in fm['tasks']:
                if ft['static']==sp and ft['response']==rp and ft['frozen_indices']==frozen:
                    rr=functional['rows'][ft['task_id']]
                    if rr['status']=='computed': t['reuse']=rr['result']
            g['endpoints'][metal]=t
        m['groups'].append(g)
    m['cache_key']=cache_key(m);write_new(root/'manifest.json',m);return preflight(root/'manifest.json')


def preflight(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['grids']!=GRIDS or m['tolerances']!=TOL or tuple(g['case_id'] for g in m['groups'])!=CASES:
        raise InvalidArtifact('comparison scope changed')
    if m['cache_key']!=cache_key({k:v for k,v in m.items() if k!='cache_key'}): raise InvalidArtifact('manifest changed')
    parent=read_json(verify(m['conductor']))
    if m['model']!=parent['model'] or m['energy_prefactor']!=parent['energy_prefactor']: raise InvalidArtifact('physical model changed')
    for p in (m['config'],m['plan'],m['software'],m['native_software'],m['functional'],m['python'],m['module'],m['stock_fixture']['arrays'],*m['implementation'].values()):verify(p)
    for g in m['groups']:
        verify(g['collection'])
        for t in g['endpoints'].values():
            for k in ('static','response','xyz','key','mask','overrides','input'):verify(t[k])
            if t['reuse']:verify(t['reuse'])
    return m


def source_model(pyddx,physical,settings,grid,dense):
    return pyddx.Model(sphere_centres=np.asfortranarray(np.array(physical['coordinates_A']).T/BOHR_TO_A),
                       sphere_radii=np.array(physical['radii_A'])/BOHR_TO_A,**{**settings,'enable_fmm':not dense},**grid)


def stock_check(pyddx,m):
    tic=time.monotonic(); f=m['stock_fixture']; model=source_model(pyddx,f['physical'],m['model'],f['grid'],True)
    with np.load(verify(f['arrays']),allow_pickle=False) as a:
        if not np.array_equal(model.cavity,a['cavity_bohr']):raise InvalidArtifact('stock source cavity differs')
        coeff=np.asfortranarray(a['Ca_average_coefficients']);phi=model.multipole_electrostatics(coeff,derivative_order=0)['phi'];psi=model.multipole_psi(coeff)
        errors=dict(phi=float(np.max(abs(phi-a['Ca_average_phi']))),psi=float(np.max(abs(psi-a['Ca_average_psi']))))
    if max(errors.values())>TOL['stock_array']:raise InvalidArtifact('parallel source full-array parity failed')
    return dict(errors=errors,tolerance=TOL['stock_array'],pass_=True,wall_seconds=time.monotonic()-tic)


def replay(t,sw,d):
    if t['reuse']:
        r=read_json(verify(t['reuse']));return dict(**r,reused=True,reused_from=t['reuse'])
    tic=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    cmd=[str(verify(sw['executable'])),str(verify(t['xyz'])),str(verify(t['mask'])),str(verify(t['overrides'])),'replay',str(verify(t['input']))]
    log=d/'native.log'
    with log.open('x') as out:p=subprocess.run(cmd,cwd=Path(t['xyz']['path']).parent,stdout=out,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    r=dict(status='failed',reused=False,receipt=dict(command=cmd,returncode=p.returncode,log=record(log),wall_seconds=time.monotonic()-tic,
           CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
    try:
        if p.returncode:raise InvalidArtifact('native replay failed')
        r.update(parse_replay(log.read_text(),t))
    except Exception as exc:r['failure_reason']=str(exc)
    return r


def solve_grid(pyddx,m,g,gridname,d):
    grid=m['grids'][gridname];physical=g['physical'];rows={};tic=time.monotonic();cpu=time.process_time()
    centers=np.array(physical['coordinates_A'])/BOHR_TO_A;radii=np.array(physical['radii_A'])/BOHR_TO_A
    dense=source_model(pyddx,physical,m['model'],grid,True);cavity=np.array(dense.cavity)
    probes=np.unique(np.linspace(0,cavity.shape[1]-1,min(256,cavity.shape[1]),dtype=int))
    # Save sources before destroying dense model; solve fresh FMM models to avoid persistent error flags.
    arrays=dict(cavity_bohr=cavity,probe_indices=probes)
    for metal,t in g['endpoints'].items():
        s=read_json(verify(t['static']));r=read_json(verify(t['response']))
        for kind,poles in distributions(s,r).items():
            label=metal+'_'+kind;start=time.monotonic();row=dict(status='failed',forward_solve_started=False)
            try:
                coeff=coefficients(poles);phi=np.array(dense.multipole_electrostatics(coeff,derivative_order=0)['phi']);psi=np.array(dense.multipole_psi(coeff))
                expected=np.zeros_like(psi)
                for l in range(3):expected[l*l:(l+1)**2]=4*np.pi*coeff[l*l:(l+1)**2]/((2*l+1)*radii**l)
                pe=float(np.max(abs(phi[probes]-cartesian_potential(cavity[:,probes].T,centers,poles))));se=float(np.max(abs(psi-expected)))
                if not np.isfinite(phi).all() or not np.isfinite(psi).all() or pe>TOL['source_au'] or se>TOL['psi']:raise InvalidArtifact('multipole source check failed')
                row.update(status='source_ready',potential_error_au=pe,psi_error=se,charge_sum_e=math.fsum(poles[:,0]))
                arrays.update({label+'_phi':phi,label+'_psi':psi})
            except Exception as exc:row['failure_reason']=str(exc)
            row['source_wall_seconds']=time.monotonic()-start;rows[label]=row
    dense=None;gc.collect();np.savez_compressed(d/'sources.npz',**arrays)
    for label,row in rows.items():
        if row['status']=='source_ready':
            model=state=None;start=time.monotonic()
            try:
                model=source_model(pyddx,physical,m['model'],grid,False)
                if not np.array_equal(model.cavity,cavity):raise InvalidArtifact('dense-source/FMM cavity differs')
                psi=np.asfortranarray(arrays[label+'_psi']);phi=arrays[label+'_phi']
                state=pyddx.State(model,psi,phi);row['forward_solve_started']=True;solve_start=time.monotonic();state.solve(tol=TOL['solver'])
                e=float(state.energy());x=np.array(state.x);factor=m['energy_prefactor'];err=abs(e-.5*math.fsum((psi*x).ravel()))*HA_TO_KCAL*factor
                if not state.is_solved or not math.isfinite(e) or not np.isfinite(x).all():raise InvalidArtifact('no finite solved state')
                row.update(status='computed',conductor_energy_hartree=e,energy_prefactor=factor,energy_kcal_mol=e*factor*HA_TO_KCAL,
                           contraction_error_kcal=err,contraction_pass=err<=TOL['contraction_kcal'],passive_energy_pass=e*factor*HA_TO_KCAL<=TOL['contraction_kcal'],
                           iterations=state.x_n_iter,solve_wall_seconds=time.monotonic()-solve_start)
            except Exception as exc:row.update(status='failed',failure_reason=str(exc))
            row['forward_wall_seconds']=time.monotonic()-start;state=model=None;gc.collect()
        write_new(d/(label+'.json'),row);print(g['case_id'],gridname,label,row['status'],row.get('energy_kcal_mol',row.get('failure_reason')),flush=True)
    out=dict(rows=rows,sources=record(d/'sources.npz'),wall_seconds=time.monotonic()-tic,CPU_seconds=time.process_time()-cpu,
             complete=all(r['status']=='computed' for r in rows.values()))
    if out['complete']:out['C_CPCM_kcal']={s:rows[s+'_average']['energy_kcal_mol']-rows[s+'_difference']['energy_kcal_mol'] for s in ('Ca','La')}
    write_new(d/'result.json',out);return out


def execute(path,case):
    import pyddx
    m=preflight(path);root=Path(path).resolve().parent;g=next(g for g in m['groups'] if g['case_id']==case)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_PER_TASK','0'))!=64 or record(pyddx.__file__)!=m['module']:
        raise InvalidArtifact('isolated module and 64-CPU allocation required')
    base=root/'cases'/case;base.mkdir(parents=True,exist_ok=True)
    with (base/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        d=base/'attempt_0001';d.mkdir();tic=time.monotonic();cpu=time.process_time()
        row=dict(case_id=case,manifest=record(path),status='failed',slurm_job_id=os.environ['SLURM_JOB_ID'],grids={},native={})
        try:
            row['stock_check']=stock_check(pyddx,m);sw=read_json(verify(m['native_software']))
            for metal,t in g['endpoints'].items():
                nd=d/('GK_'+metal);nd.mkdir();r=replay(t,sw,nd);write_new(nd/'result.json',r);row['native'][metal]=r
            if not all(x['status']=='computed' for x in row['native'].values()):raise InvalidArtifact('native GK functional unavailable')
            gk={s:r['GK_cross_kcal'] for s,r in row['native'].items()}
            for gridname in m['grids']:
                gd=d/gridname;gd.mkdir();result=solve_grid(pyddx,m,g,gridname,gd)
                if result['complete']:
                    result['R_old_kcal']=g['old']['R_kcal'];result['C_GK_kcal']=gk
                    result['R_new_kcal']=corrected_contrast(g['old']['R_kcal'],result['C_CPCM_kcal'],gk)
                    result['correction_R_kcal']=result['R_new_kcal']-result['R_old_kcal']
                row['grids'][gridname]=result;write_new(gd/'score.json',result)
            row['status']='computed' if all(x['complete'] for x in row['grids'].values()) else 'partial_failure'
        except Exception as exc:row['failure_reason']=str(exc)
        row.update(wall_seconds=time.monotonic()-tic,CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_new(d/'result.json',row);return row


def collect(path,output):
    m=preflight(path);root=Path(path).resolve().parent;cases={};margins={};checks=[]
    for g in m['groups']:
        case=g['case_id'];p=root/'cases'/case/'attempt_0001/result.json'
        r=read_json(p) if p.exists() else dict(status='pending',grids={},native={})
        if not p.exists():
            # Workers publish each completed grid before the next refinement.
            # Retain its receipt and report useful partial scores immediately.
            for grid in GRIDS:
                partial=p.parent/grid/'score.json'
                if partial.exists():
                    r['grids'][grid]={**read_json(partial),'partial_receipt':record(partial)}
                    r['status']='running_with_partial_results'
            for metal in ('Ca','La'):
                native=p.parent/('GK_'+metal)/'result.json'
                if native.exists():r['native'][metal]=read_json(native)
        if p.exists() and r['manifest']!=record(path):raise InvalidArtifact('result belongs to another manifest')
        cases[case]=dict(status=r['status'],old_R_kcal=g['old']['R_kcal'],evidence=g['old']['evidence'],old_components=g['old']['components_R_kcal'],grids=r['grids'],result=record(p) if p.exists() else None)
        for metal,x in r.get('native',{}).items():checks.extend(dict(case=case,metal=metal,**c) for c in x.get('checks',[]))
        if all(r.get('grids',{}).get(k,{}).get('complete') for k in GRIDS):
            a=r['grids']['primary'];b=r['grids']['refined'];delta=b['R_new_kcal']-a['R_new_kcal']
            cases[case]['refinement_R_kcal']=delta;checks.append(dict(case=case,name='paired_refinement',error=delta,tolerance=TOL['refinement_kcal'],pass_=abs(delta)<=TOL['refinement_kcal']))
    for grid in GRIDS:
        rows=[]
        for ln in CASES[3:]:
            for ca in CASES[:3]:
                a=cases[ln];b=cases[ca];old=a['old_R_kcal']-b['old_R_kcal'];aa=a['grids'].get(grid,{});bb=b['grids'].get(grid,{})
                new=aa.get('R_new_kcal')-bb.get('R_new_kcal') if 'R_new_kcal' in aa and 'R_new_kcal' in bb else None
                change=a.get('refinement_R_kcal',0)-b.get('refinement_R_kcal',0) if 'refinement_R_kcal' in a and 'refinement_R_kcal' in b else None
                rows.append(dict(La_case=ln,Ca_case=ca,old_margin_kcal=old,new_margin_kcal=new,old_direction_pass=old>TOL['direction_kcal'],new_direction_pass=new>TOL['direction_kcal'] if new is not None else None,margin_refinement_kcal=change))
        margins[grid]=rows
    result=dict(protocol=PROTOCOL,manifest=record(path),reporter=record(__file__),cases=cases,margins=margins,checks=checks,complete=all(x['status']=='computed' for x in cases.values()),
                interpretation='Exploratory frozen-response contrasts; no absolute calibration; replicates/sites grouped, not independent biological observations.',
                inherited_qualification='Conductor endpoint/rotation qualification previously failed; not waived.',baseline_changed=False,new_DFT_calls=0,new_MACE_calls=0)
    write_new(output,result);return result


def report(collection,output):
    r=read_json(collection);m=read_json(verify(r['manifest']))
    lines=['# Complete conductor transfer: discrimination comparison','',
           'Protocol: `'+PROTOCOL+'`. Production baseline unchanged. This compares the prior',
           'GK hybrid with its frozen-response conductor correction, not two calibrated classifiers.',
           'All structures are consumed development cases; crystal replicates and sites are grouped.',
           '', '| Comparison | Prior GK | Conductor primary | Conductor refined | Refinement change |',
           '|---|---:|---:|---:|---:|']
    fmt=lambda v:'unavailable' if v is None else f'{v:+.6f}'
    for a,b in zip(r['margins']['primary'],r['margins']['refined']):
        lines.append('| '+a['La_case']+' − '+a['Ca_case']+' | '+' | '.join(fmt(v) for v in
          (a['old_margin_kcal'],a['new_margin_kcal'],b['new_margin_kcal'],b['margin_refinement_kcal']))+' |')
    lines+=['','Values are kcal/mol differences in R = E_Ca − E_La. The predeclared directional',
            'screen is >0.02. These twelve comparisons are not twelve independent biological tests.',
            'GGR has direct same-assay Ca-favoring evidence; alpha is a qualified strong-site',
            'comparison and parvalbumin is supporting cross-study evidence. No absolute bands.', '', '## Directional summary','']
    for family,subset in [('alpha',CASES[3:5]),('parvalbumin',CASES[5:])]:
        rows=[x for x in r['margins']['primary'] if x['La_case'] in subset]
        fine=[x for x in r['margins']['refined'] if x['La_case'] in subset]
        lines.append(f'- {family}: prior {sum(x["old_direction_pass"] for x in rows)}/{len(rows)}; '
                     f'primary {sum(x["new_direction_pass"] is True for x in rows)}/{len(rows)}; '
                     f'refined {sum(x["new_direction_pass"] is True for x in fine)}/{len(fine)} '
                     f'({sum(x["new_margin_kcal"] is not None for x in fine)} available).')
    lines+=['','## Numerical qualification and cost','',
            r['inherited_qualification'],
            'New native-functional and paired-refinement checks: '
            f'{sum(c["pass_"] for c in r["checks"])}/{len(r["checks"])} passing.',
            'Source representation, solved-state and contraction checks remain in each endpoint receipt.',
            'A small paired grid change does not erase earlier endpoint or rigid-transform failures.',
            'No new DFT or MACE inference. Native GK evaluations reuse four exact prior endpoints.',
            'Scheduler allocation cost is recorded separately from worker timings; cache preparation',
            'and prior development costs are not implied to be free.', '', '## Reproduction','',
            'From the repository root (choose a fresh output filename):','', '```bash',
            '/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_conductor_discrimination.py collect --manifest '+r['manifest']['path']+' --output '+str(Path(collection).with_name('recollected.json')),
            '```','', 'Collection: `'+str(Path(collection).resolve())+'`.',
            'Manifest SHA256: `'+r['manifest']['sha256']+'`. Unrounded values and components are retained.','']
    with Path(output).open('x') as f:f.write('\n'.join(lines))
    return dict(status='reported',collection=record(collection),report=record(output))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('report');q.add_argument('--collection',required=True);q.add_argument('--output',required=True)
    for op in ('dry-run','execute','collect'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op=='execute':q.add_argument('--case',choices=CASES,required=True)
        if op=='collect':q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output)
    elif a.command=='report':r=report(a.collection,a.output)
    elif a.command=='dry-run':r=preflight(a.manifest)
    elif a.command=='execute':r=execute(a.manifest,a.case)
    else:r=collect(a.manifest,a.output)
    print(json.dumps({k:r[k] for k in ('protocol','status','complete','failure_reason') if k in r}))
    if a.command=='execute' and r['status']!='computed':raise SystemExit(1)
