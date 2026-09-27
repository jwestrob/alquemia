"""Finite real-source native derivative checks; no protein relaxation or classifier.

Twenty endpoints, two metals, two physical source-coordinate derivative modes,
two step sizes, exact repeats and joint rigid transforms. Root owns submission.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import numpy as np
from affordable_common import (BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, cache_key,
                               read_json, record, verify, write_new)
from affordable_workflow import dry_run, execute
from mace_omol_vacuum import METHOD, embedded_input, parse_endpoint
from metal_environment_reference import read_pcgrad

PROTOCOL = 'nikasha_native_embedded_selected_derivatives_v1'
STEPS = (.001, .0005)
TOLERANCES = dict(repeat_energy_Eh=1e-7, repeat_gradient_Eh_bohr=1e-6,
                  rigid_energy_Eh=1e-5, rigid_gradient_Eh_bohr=1e-4,
                  fd_absolute_kcal_mol_A=.05, fd_relative=.005)


def check_pins(value):
    if isinstance(value, dict):
        if 'path' in value and 'sha256' in value:
            verify(value)
        else:
            for child in value.values(): check_pins(child)
    elif isinstance(value, list):
        for child in value: check_pins(child)


def xyz_data(path):
    lines = Path(path).read_text().splitlines(keepends=True)
    n = int(lines[0]); rows = [line.split() for line in lines[2:]]
    if len(rows) != n or any(len(row) != 4 for row in rows):
        raise InvalidArtifact('invalid XYZ inventory')
    return lines, [r[0] for r in rows], np.array([[float(v) for v in r[1:]] for r in rows])


def pc_data(path):
    lines = Path(path).read_text().splitlines(keepends=True)
    rows = np.array([[float(v) for v in line.split()] for line in lines[1:]])
    if rows.shape != (int(lines[0]),4) or not np.isfinite(rows).all():
        raise InvalidArtifact('invalid point-charge inventory')
    return lines, rows[:,0], rows[:,1:]


def rotation():
    c,s = np.cos(.37), np.sin(.37)
    return np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])


def construct_modes(config):
    """Derive modes from pinned source maps and actual serialized A coordinates."""
    mapping = read_json(verify(config['core_mapping']))
    atoms = read_json(verify(config['environments']['A']['atoms']))
    _, q, y = pc_data(verify(config['environments']['A']['pointcharges']))
    if len(atoms)!=len(y) or not np.allclose(y,[a['xyz_A'] for a in atoms],atol=1e-12,rtol=0):
        raise InvalidArtifact('field coordinate mapping differs')
    if not np.allclose(q,[a['charge_e'] for a in atoms],atol=1e-14,rtol=0):
        raise InvalidArtifact('field charge mapping differs')
    def index(identifier):
        selected=[i for i,a in enumerate(atoms) if a['id']==identifier]
        if len(selected)!=1: raise InvalidArtifact('missing/duplicate mapped source: '+identifier)
        return selected[0]
    hi,oi,ci = [index('A/159/ /'+name) for name in ('HG1','OG1','CB')]
    axis = y[oi]-y[ci]; axis/=np.linalg.norm(axis)
    tangent = np.cross(axis,y[hi]-y[oi]); radius=np.linalg.norm(tangent)
    if radius < 1e-6: raise InvalidArtifact('degenerate hydroxyl rotation')
    cap = [a for a in mapping if a.get('id')=='cap/A:GLU177']
    if len(cap)!=1: raise InvalidArtifact('missing/duplicate Glu177 cap')
    cap=cap[0]; retained=np.asarray(cap['retained_xyz_A']); omitted=np.asarray(cap['omitted_xyz_A'])
    d=omitted-retained; length=np.linalg.norm(d); u=d/length
    direction=np.eye(3)[np.argmin(abs(u))]; direction-=u*np.dot(u,direction); direction/=np.linalg.norm(direction)
    jac=cap['length_A']/length*(np.eye(3)-np.outer(u,u))
    if not np.allclose(jac,cap['jacobian_omitted'],atol=1e-12,rtol=0):
        raise InvalidArtifact('source cap Jacobian differs')
    mm=[i for i,a in enumerate(atoms) if a['id']==cap['omitted_source_id']]
    if len(mm)>1: raise InvalidArtifact('duplicate omitted CA field mapping')
    _,_,x=xyz_data(verify(config['endpoints']['Ca']['xyz']))
    if not np.allclose(x[cap['qm_index']],cap['xyz_A'],atol=1e-12,rtol=0):
        raise InvalidArtifact('serialized cap mapping differs')
    offset=x[cap['qm_index']]-(retained+cap['length_A']*u)
    return dict(mm=dict(h_index=hi,o_index=oi,c_index=ci,axis=axis.tolist(),
                        radius_A=float(radius), tangent=(tangent/radius).tolist()),
                boundary=dict(cap_index=cap['qm_index'],retained_A=retained.tolist(),
                              omitted_A=omitted.tolist(),length_A=cap['length_A'],
                              direction=direction.tolist(),jacobian_omitted=jac.tolist(),
                              cap_tangent=(jac@direction).tolist(),offset_A=offset.tolist(),
                              mm_index=mm[0] if mm else None,
                              source_id=cap['omitted_source_id']))


def displaced(x, y, modes, mode, h=0.):
    """Return coordinates with only the declared source-induced changes."""
    xx,yy=x.copy(),y.copy()
    if mode=='repeat': return xx,yy
    if mode=='rigid': return x@rotation().T+(.173,.117,.231), y@rotation().T+(.173,.117,.231)
    if mode=='mm':
        m=modes['mm'];axis=np.asarray(m['axis']); v=y[m['h_index']]-y[m['o_index']]
        angle=h/m['radius_A']
        yy[m['h_index']]=y[m['o_index']]+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+axis*np.dot(axis,v)*(1-np.cos(angle))
    elif mode=='boundary':
        m=modes['boundary'];r=np.asarray(m['retained_A']); direction=np.asarray(m['direction'])
        d=np.asarray(m['omitted_A'])+h*direction-r
        xx[m['cap_index']]=r+m['length_A']*d/np.linalg.norm(d)+m['offset_A']
        if m['mm_index'] is not None: yy[m['mm_index']]+=h*direction
    else: raise InvalidArtifact('unknown declared mode')
    return xx,yy


def rendered_geometry(config, metal, modes, mode, h):
    xl,symbols,x=xyz_data(verify(config['endpoints'][metal]['xyz']))
    pl,q,y=pc_data(verify(config['environments']['A']['pointcharges']))
    xx,yy=displaced(x,y,modes,mode,h)
    for i in np.flatnonzero(np.any(xx!=x,axis=1)):
        xl[i+2]=symbols[i]+' '+' '.join(f'{v:.16g}' for v in xx[i])+'\n'
    for i in np.flatnonzero(np.any(yy!=y,axis=1)):
        pl[i+1]=' '.join(f'{v:.17g}' for v in (q[i],*yy[i]))+'\n'
    return ''.join(xl),''.join(pl)


def inventory():
    return [(metal,mode,h) for metal in ('Ca','La') for mode in ('mm','boundary')
            for step in STEPS for h in (-step,step)]+[(metal,mode,0.) for metal in ('Ca','La') for mode in ('repeat','rigid')]


def task_id(metal,mode,h):
    return f'{metal}_{mode}_'+('origin' if not h else ('plus' if h>0 else 'minus')+str(abs(h)).replace('.','p'))


def centers(collection, inputs):
    from ggr_sensitivity import executed
    c=read_json(collection); cm,actual=executed(verify(c['manifest']))
    if cm['inputs']['sha256']!=record(inputs)['sha256'] or cm['method']!=METHOD:
        raise InvalidArtifact('center source/method differs')
    rows={}
    for metal in ('Ca','La'):
        row=c['rows'][metal+'_A']; task=next(t for t in cm['tasks'] if t['task_id']==metal+'_A')
        if row['status']!='complete': raise InvalidArtifact('center endpoint unavailable')
        observed=actual[metal+'_A']
        if observed['status']!='complete' or any(observed[k]!=row[k] for k in ('receipt','output','energy_hartree')):
            raise InvalidArtifact('center receipt failed original-manifest verification')
        for key in ('output','engrad','receipt','pointcharge_gradient'): verify(row[key])
        parsed=parse_endpoint(task,verify(row['output']),verify(row['engrad']),permanent_field=True)
        if parsed['energy_hartree']!=row['energy_hartree'] or parsed['gradient_kcal_mol_per_A']!=row['gradient_kcal_mol_per_A']:
            raise InvalidArtifact('center parsed energy/gradient differs')
        rows[metal]=row
    return rows,cm


def prepare(inputs, collection, plan, output, *, workers, mpi_ranks):
    if workers<1 or workers>20 or mpi_ranks<1: raise InvalidArtifact('one to twenty workers and positive ranks required')
    config=read_json(inputs);check_pins(config); center,cm=centers(collection,inputs); modes=construct_modes(config)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir()
    pins={}
    for source in sorted(Path(__file__).parent.glob('*.py')):
        shutil.copyfile(source,impl/source.name); pins[source.name]=record(impl/source.name)
    shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py')
    shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
    for name in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[name]=record(impl/name)
    tasks=[]
    for metal,mode,h in inventory():
        td=out/task_id(metal,mode,h);td.mkdir();xtext,ptext=rendered_geometry(config,metal,modes,mode,h)
        (td/'core.xyz').write_text(xtext);(td/'environment.pc').write_text(ptext)
        charge=config['endpoints'][metal]['charge'];(td/'endpoint.inp').write_text(embedded_input(charge))
        tasks.append(dict(task_id=td.name,metal=metal,mode=mode,displacement_A=h,charge=charge,multiplicity=1,
                          input=record(td/'endpoint.inp'),xyz=record(td/'core.xyz'),pointcharges=record(td/'environment.pc'),
                          output_path=str(td/'endpoint.out'),engrad_path=str(td/'endpoint.engrad'),task_type='analytic_gradient'))
    m=dict(protocol_id=PROTOCOL,method=METHOD,inputs=record(inputs),center_collection=record(collection),agreement=record(plan),
           tasks=tasks,modes=modes,tolerances=TOLERANCES,orca=cm['orca'],implementation=pins,
           execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
           execution_resources={'mpi_ranks':mpi_ranks,'concurrent_tasks':workers},
           energy_scope='embedded_electronic_component_only',compute_budget=None,wall_time_limit=None)
    for t in tasks:t['cache_key']=cache_key(dict(task=t,protocol=PROTOCOL,method=METHOD,orca=m['orca'],tolerances=TOLERANCES,inputs=m['inputs'],center=m['center_collection'],implementation=pins))
    write_new(out/'manifest.json',m)
    return dict(status='prepared',manifest=record(out/'manifest.json'),task_count=20)


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['method']!=METHOD or m['tolerances']!=TOLERANCES:raise InvalidArtifact('protocol/method/tolerances changed')
    check_pins(m['implementation']);verify(m['agreement']);verify(m['orca'])
    config=read_json(verify(m['inputs']));check_pins(config);centers(verify(m['center_collection']),verify(m['inputs']))
    modes=construct_modes(config)
    if modes!=m['modes']:raise InvalidArtifact('physical mappings changed')
    if len(m['tasks'])!=20 or {(t['metal'],t['mode'],t['displacement_A']) for t in m['tasks']}!=set(inventory()):
        raise InvalidArtifact('twenty-endpoint inventory changed')
    for t in m['tasks']:
        x,p=rendered_geometry(config,t['metal'],modes,t['mode'],t['displacement_A'])
        if verify(t['xyz']).read_text()!=x or verify(t['pointcharges']).read_text()!=p:raise InvalidArtifact('geometry/field displacement changed')
        if t['charge']!=config['endpoints'][t['metal']]['charge'] or t['multiplicity']!=1 or verify(t['input']).read_text()!=embedded_input(t['charge']):raise InvalidArtifact('electronic state/Hamiltonian changed')
        bare={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key(dict(task=bare,protocol=PROTOCOL,method=METHOD,orca=m['orca'],tolerances=TOLERANCES,inputs=m['inputs'],center=m['center_collection'],implementation=m['implementation'])):raise InvalidArtifact('cache identity changed')
    return dry_run(manifest)


def gradients(row,count):
    return np.array(row['gradient_kcal_mol_per_A']),read_pcgrad(verify(row['pointcharge_gradient']),count)*HA_TO_KCAL/BOHR_TO_A


def collect(manifest):
    from ggr_sensitivity import executed
    validate(manifest);m,rows=executed(manifest);config=read_json(verify(m['inputs']));orig,_=centers(verify(m['center_collection']),verify(m['inputs']))
    count=config['environments']['A']['count'];modes=m['modes'];checks=[]
    for t in m['tasks']:
        row=rows[t['task_id']]
        if row['status']!='complete':continue
        try:
            row.update(parse_endpoint(t,verify(row['output']),t['engrad_path'],permanent_field=True))
            pc=Path(t['output_path']).with_name('endpoint.runtime.pcgrad');read_pcgrad(pc,count)
            row['pointcharge_gradient']=record(pc)
        except (ValueError,OSError) as exc:row.update(status='invalid',reason=str(exc),energy_hartree=None)
    for metal in ('Ca','La'):
        gx,gy=gradients(orig[metal],count)
        for mode in ('mm','boundary'):
            if mode=='mm': projection=float(gy[modes['mm']['h_index']]@modes['mm']['tangent'])
            else:
                b=modes['boundary'];projection=float(gx[b['cap_index']]@b['cap_tangent'])
                if b['mm_index'] is not None:projection+=float(gy[b['mm_index']]@b['direction'])
            tolerance=.05+.005*abs(projection);fd=[]
            for h in STEPS:
                minus,plus=(rows[task_id(metal,mode,s*h)] for s in (-1,1))
                value=(plus['energy_hartree']-minus['energy_hartree'])*HA_TO_KCAL/(2*h) if all(r['status']=='complete' for r in (minus,plus)) else None
                fd.append(value)
                checks.append(dict(metal=metal,mode=mode,step_A=h,analytic_kcal_mol_A=projection,fd_kcal_mol_A=value,
                                   residual_kcal_mol_A=None if value is None else value-projection,tolerance=tolerance,
                                   status='unavailable' if value is None else ('passed' if abs(value-projection)<=tolerance else 'failed')))
            residual=None if None in fd else fd[0]-fd[1]
            checks.append(dict(metal=metal,mode=mode+'_step_convergence',residual_kcal_mol_A=residual,tolerance=tolerance,
                               status='unavailable' if residual is None else ('passed' if abs(residual)<=tolerance else 'failed')))
        for mode in ('repeat','rigid'):
            row=rows[task_id(metal,mode,0.)]
            if row['status']!='complete':checks.append(dict(metal=metal,mode=mode,status='unavailable'));continue
            xx,yy=gradients(row,count);rotation_matrix=rotation() if mode=='rigid' else np.eye(3)
            er=abs(row['energy_hartree']-orig[metal]['energy_hartree'])
            qr=float(abs(xx-gx@rotation_matrix.T).max()*BOHR_TO_A/HA_TO_KCAL)
            mr=float(abs(yy-gy@rotation_matrix.T).max()*BOHR_TO_A/HA_TO_KCAL)
            ok=er<=TOLERANCES[mode+'_energy_Eh'] and max(qr,mr)<=TOLERANCES[mode+'_gradient_Eh_bohr']
            checks.append(dict(metal=metal,mode=mode,energy_residual_Eh=er,core_gradient_residual_Eh_bohr=qr,mm_gradient_residual_Eh_bohr=mr,status='passed' if ok else 'failed'))
    return dict(protocol_id=PROTOCOL,manifest=record(manifest),rows=rows,checks=checks,
                status='qualified_selected_directions_only' if all(c['status']=='passed' for c in checks) else 'not_qualified',
                full_hybrid_status='unsupported',classification=None,affinity=None,
                measured_execution_events=str(Path(manifest).parent/'budget_events.jsonl'))


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    a.add_argument('--inputs')
    a.add_argument('--source-manifest')
    a.add_argument('--collection','--source-collection',dest='collection',required=True)
    a.add_argument('--plan','--agreement',dest='plan',required=True)
    a.add_argument('--output',required=True)
    for k in ('workers','mpi-ranks'):a.add_argument('--'+k,type=int,required=True)
    for op in ('dry-run','collect','execute'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=p.parse_args()
    if a.op=='prepare':
        if bool(a.inputs)==bool(a.source_manifest):p.error('exactly one of --inputs or --source-manifest required')
        if a.source_manifest:
            source=record(a.source_manifest)
            if source!=read_json(a.collection)['manifest']:raise InvalidArtifact('source manifest differs from center collection')
            inputs=verify(read_json(a.source_manifest)['inputs'])
        else:inputs=a.inputs
        result=prepare(inputs,a.collection,a.plan,a.output,workers=a.workers,mpi_ranks=a.mpi_ranks)
    elif a.op=='dry-run':result=validate(a.manifest)
    elif a.op=='collect':result=collect(a.manifest)
    else:
        validate(a.manifest);resources=read_json(a.manifest)['execution_resources']
        if int(os.environ.get('SLURM_NTASKS','0'))<resources['concurrent_tasks']*resources['mpi_ranks']:raise InvalidArtifact('requested MPI layout exceeds actual slots')
        os.environ['METAL_ENV_WORKERS']=str(resources['concurrent_tasks']);result=execute(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
