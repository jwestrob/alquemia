"""Native compact-context coordination proposals; original-core DFT adjudication."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import resource
import shutil
import time
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from hydration_square import endpoint
from mace_hybrid import EV_TO_KCAL, check_atoms, write_xyz
from mace_site_kinematics import Kinematics

CASES=('1H4I','4MAE','1GLG','1F6S','6IP9')

PROTOCOL = 'native_OMOL_context_physical_coordination_proposal_v1'
SETTINGS = {'optimizer':'SLSQP', 'maxiter':200, 'ftol_eV':1e-9,
            'physical_heavy_displacement_A':.20, 'angular_bound_rad':.20,
            'selection':'final_feasible_nonincreasing_context_energy_iterate',
            'fixed':'PQQ_ARG_second_shell_waters_scaffold',
            'native_score':'original_core_r2SCAN3c_CPCM_Water_DefGrid3_SP'}


def donor_units(state):
    p=state['parent'];g=state['graph'];result=[]
    if 'fixed_core' in p:
        for f in p['qm_fragments']:
            if f['kind']!='fixed_core_protein_sidechain':continue
            ch,rn,num,ic=re.fullmatch(r'([^:]+):([A-Z]+)(-?\d+)([A-Za-z]?)',f['id']).groups()
            r=g.locate(dict(chain=ch,resname=rn,resnum=int(num),insertion_code=ic))
            result.append(('sidechain',r.key))
    else:
        for f in p['charge_ledger']:
            if f['kind']!='water':result.append((f['kind'],g.locate(f['source']).key))
    return sorted(result,key=lambda z:(z[1],z[0]))


def geometry(state,prep,context_rows,core_rows):
    """Reuse established exact rotations, with both compact and core maps."""
    from mace_site_coordinates import CHI
    from second_shell_context import atom_key
    g=state['graph'];adj={k:set(v for v in g.required[k] if v in g.atoms) for k in g.atoms}
    for h,p in g.hparents.items():adj[h].add(p);adj[p].add(h)
    modes=[];keys=set()
    for kind,rkey in donor_units(state):
        r=g.residues[rkey]
        if kind=='sidechain':
            if r.canonical_resname not in CHI:raise InvalidArtifact('unsupported original donor torsions')
            for j,(an,bn) in enumerate(CHI[r.canonical_resname],1):
                a,b=g.key(rkey,an),g.key(rkey,bn);seen=set();stack=[b]
                while stack:
                    k=stack.pop()
                    if k in seen:continue
                    seen.add(k);stack.extend(v for v in adj[k] if {v,k}!={a,b})
                if a in seen or any(k[:2]!=rkey for k in seen):raise InvalidArtifact('ring or nonlocal donor torsion')
                modes.append(dict(id=f'{r.chain}/{r.resnum}/chi{j}',kind='sidechain_torsion',unit='radian',
                                  axis_keys=[a,b],moving_keys=sorted(seen)))
        elif kind=='backbone_carbonyl':
            n=g.amide_next.get(rkey)
            if n is None or g.residues[n[:2]].canonical_resname=='PRO':raise InvalidArtifact('unsupported peptide crankshaft')
            hs=[k for k in adj[n] if g.meta[k]['element']=='H']
            if len(hs)!=1:raise InvalidArtifact('peptide crankshaft requires one physical N-H')
            modes.append(dict(id=f'{r.chain}/{r.resnum}/peptide_crankshaft',kind='peptide_crankshaft',unit='radian',
                axis_keys=[g.key(rkey,'CA'),g.key(n[:2],'CA')],moving_keys=[g.key(rkey,'C'),g.key(rkey,'O'),n,*hs]))
        else:raise InvalidArtifact('unsupported donor ledger kind')
    for m in modes:keys.update(m['axis_keys']);keys.update(m['moving_keys'])
    moving=set(k for m in modes for k in m['moving_keys'])
    for k in moving:keys.update(adj[k])
    for mapping in (state['maps']+state['opaque'],prep['mapping']['source_to_qm']):
        for a in mapping:
            for name in ('source','retained','omitted'):
                if name in a:keys.add(atom_key(a[name]))
    ordered=sorted(keys);lookup={k:i+1 for i,k in enumerate(ordered)}
    positions=[list(core_rows[0][1:])]+[[g.atoms[k].pos.x,g.atoms[k].pos.y,g.atoms[k].pos.z] for k in ordered]
    # Endpoint-specific archived H preparation replaces source positions exactly.
    seen={}
    for mapping,rows in ((prep['mapping']['source_to_qm'],context_rows),(state['maps']+state['opaque'],core_rows)):
        for a in mapping:
            if 'source' not in a:continue
            key=atom_key(a['source']);v=list(rows[a['qm_index']][1:])
            if key in seen and not np.allclose(v,seen[key],rtol=0,atol=1e-9):raise InvalidArtifact('core/context source mismatch')
            positions[lookup[key]]=v;seen[key]=v
    physical=np.array(positions)
    for m in modes:
        m['axis_indices']=[lookup[k] for k in m.pop('axis_keys')]
        m['moving_indices']=[lookup[k] for k in m.pop('moving_keys')]
    modes=[dict(id='metal_'+s,kind='metal_translation',unit='angstrom',axis=v,moving_indices=[0])
           for s,v in zip('xyz',np.eye(3).tolist())]+modes
    bonds=sorted({tuple(sorted((lookup[k],lookup[v]))) for k in ordered for v in adj[k] if v in lookup})
    data={'positions_A':positions,'modes':modes,'metal_index':0,
          'heavy_indices':[0]+[lookup[k] for k in ordered if g.meta[k]['element'] not in ('H','D')],
          'physical_ids':['metal']+[list(k) for k in ordered],'bonds':bonds,
          'source_atom_metadata':[None]+[g.meta[k] for k in ordered]}
    def mapped(mapping,rows):
        links=[]
        for a in mapping:
            i=a['qm_index']
            if 'source' in a:links.append([i,'source',lookup[atom_key(a['source'])]])
            elif a['kind']=='sigma_link_H':
                x,y=[lookup[atom_key(a[s])] for s in ('retained','omitted')]
                length=float(a.get('length_A',np.linalg.norm(np.array(rows[i][1:])-physical[x])))
                expected=physical[x]+length*(physical[y]-physical[x])/np.linalg.norm(physical[y]-physical[x])
                links.append([i,'cap',x,y,length,expected.tolist()])
            elif a['kind']!='opaque_cofactor':raise InvalidArtifact('unknown source mapping')
        return data|{'core_positions_A':[list(a[1:]) for a in rows],'core_links':links}
    original=mapped(state['maps']+state['opaque'],core_rows)
    context=mapped(prep['mapping']['source_to_qm'],context_rows)
    q=np.zeros(len(modes));probe=np.full(len(modes),.001)
    for d,rows in ((original,core_rows),(context,context_rows)):
        kin=Kinematics(d)
        if not np.allclose(kin.evaluate(q)[1],np.array([a[1:] for a in rows]),atol=1e-10,rtol=0):raise InvalidArtifact('nonidentity map at origin')
        checks=kin.check(probe)
        if not checks['pass']:raise InvalidArtifact('physical coordinate map failed: '+str(checks))
    return {'core':original,'context':context,'checks':checks,
            'fixed_core_indices':[a['qm_index'] for a in state['maps'] if 'source' in a and atom_key(a['source']) in state['water_keys']]}


def prepare(static_manifest,pqq_collection,pqq_dft,agreement,output):
    from second_shell_context import parent_state
    src=read_json(static_manifest);hcol=read_json(pqq_collection);dft=read_json(pqq_dft)
    hp={r['task_id']:r for pin in hcol['results'] for r in [read_json(verify(pin))]}
    dm=read_json(verify(dft['manifest']));cfg=read_json(verify(src['configuration']))
    old=read_json(Path(static_manifest).parent/'mace_manifest.json')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir()
    pins={}
    for p in Path(__file__).parent.glob('*.py'):
        z=impl/p.name;shutil.copyfile(p,z);pins[p.name]=record(z)
    tasks=[]
    for item in cfg['cases']:
        case=item['case'];state=parent_state(item,cfg['topology'])
        p=next(s['preparation'] for s in src['states'] if s['case']==case);prep=read_json(verify(p))
        for metal in ('Ca','La'):
            tid=case+'__'+metal;t=next(t for t in src['tasks'] if t['case']==case and t['metal']==metal)
            core=item['endpoints'][metal];context=t['xyz']
            if case in CASES[:2]:
                h=hp[tid];dt=next(t for t in dm['tasks'] if t['task_id']==tid)
                core={'xyz':dt['xyz'],'input':dt['input'],'output':record(dt['output_path']),
                      'receipt':record(dt['output_path']+'.execution.json'),'charge':dt['charge']}
                core['energy_hartree']=endpoint(core['output'],core['receipt'],core['xyz'],core['input'])['energy_hartree']
                context=h['context_xyz']
            rows=xyz(verify(context));cr=xyz(verify(core['xyz']));geo=geometry(state,prep,rows,cr)
            gp=out/(tid+'_geometry.json');write_new(gp,geo)
            task={'task_id':tid,'case':case,'metal':metal,'context_xyz':context,'charge':t['charge'],
                  'core':core,'geometry':record(gp),'preparation':p,'state':check_atoms(rows,t['charge']),
                  'group':item['group'],'evidence_stratum':item['evidence_stratum']}
            task['cache_key']=cache_key(task|{'protocol':PROTOCOL,'settings':SETTINGS,'model':old['model'],'implementation':pins})
            tasks.append(task)
    m={'protocol_id':PROTOCOL,'settings':SETTINGS,'tasks':tasks,'source':record(static_manifest),
       'pqq_preparation':record(pqq_collection),'pqq_dft':record(pqq_dft),'agreement':record(agreement),
       'model':old['model'],'software':old['software'],'implementation':pins,
       'orca':src['orca'],'execution_policy':src['execution_policy'],'compute_budget':None}
    write_new(out/'manifest.json',m)
    return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or len(m['tasks'])!=10:raise InvalidArtifact('scope differs')
    for pin in m['implementation'].values():verify(pin)
    for name in ('source','pqq_preparation','pqq_dft','agreement','software'):verify(m[name])
    modes={}
    for t in m['tasks']:
        rows=xyz(verify(t['context_xyz']));verify(t['core']['xyz']);geo=read_json(verify(t['geometry']))
        if check_atoms(rows,t['charge'])!=t['state']:raise InvalidArtifact('charge/state differs')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key(payload|{'protocol':PROTOCOL,'settings':SETTINGS,'model':m['model'],'implementation':m['implementation']}):
            raise InvalidArtifact('cache differs')
        modes[t['task_id']]=len(geo['core']['modes'])
    return {'status':'pass','tasks':10,'modes':modes,'new_energy_calls':0}


def constraints(kin,q):
    p,_,j,_=kin.evaluate(q);ids=kin.heavy;d=p[ids]-kin.positions[ids]
    return SETTINGS['physical_heavy_displacement_A']**2-np.sum(d*d,axis=1),-2*np.einsum('ij,mij->im',d,j[:,ids])


def proposal_checks(geo,q,initial_context,initial_core):
    ck=Kinematics(geo['context']);ok=Kinematics(geo['core']);p,c,j,cj=ck.evaluate(q);_,core,_,_=ok.evaluate(q)
    extent=np.linalg.norm(p[ck.heavy]-ck.positions[ck.heavy],axis=1)
    if max(extent)>.20+1e-7 or max(abs(np.array(q)[3:]))>.20+1e-7:raise InvalidArtifact('proposal outside physical domain')
    checks=ck.check(q)
    if not checks['pass'] or not ok.check(q)['pass']:raise InvalidArtifact('proposal violates physical map')
    water=geo['fixed_core_indices']
    if not np.array_equal(core[water],np.array([a[1:] for a in initial_core])[water]):raise InvalidArtifact('water coordinates changed')
    distance=np.linalg.norm(c[:,None]-c[None,:],axis=2);np.fill_diagonal(distance,np.inf)
    if distance.min()<.45:raise InvalidArtifact('proposal atom/cap overlap')
    return checks|{'maximum_physical_heavy_displacement_A':float(max(extent)),
                  'boundary_heavy_atoms':int(sum(extent>=.20-1e-5)),
                  'water_coordinates_unchanged':True}


def execute(manifest):
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import input_batch
    from scipy.optimize import minimize
    if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():raise InvalidArtifact('allocated GPU required')
    validate(manifest);m=read_json(manifest);root=Path(manifest).parent;began=time.monotonic()
    torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']));torch.set_default_dtype(torch.float64)
    torch.cuda.reset_peak_memory_stats();calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cuda',default_dtype='float64')
    model=calc.models[0]
    if type(model).__name__!='ScaleShiftMACE' or dict(model.embedding_specs)!=m['model']['embedding_specs']:raise InvalidArtifact('native model differs')
    for p in model.parameters():p.requires_grad_(False)
    versions={k:p._version for k,p in model.named_parameters()};results=[]
    for t in m['tasks']:
        dest=root/'tasks'/t['task_id'];rp=dest/'result.json'
        if rp.exists():
            previous=read_json(rp)
            if previous['cache_key']!=t['cache_key']:raise InvalidArtifact('partial result cache mismatch')
            results.append(record(rp));continue
        dest.mkdir(parents=True,exist_ok=False);ts=time.monotonic()
        rows=xyz(verify(t['context_xyz']));cr=xyz(verify(t['core']['xyz']));geo=read_json(verify(t['geometry']))
        kin=Kinematics(geo['context']);original=Kinematics(geo['core']);n=len(kin.modes);zero=np.zeros(n)
        atoms=Atoms([a[0] for a in rows],positions=[a[1:] for a in rows],pbc=False)
        atoms.info.update(charge=t['charge'],spin=1);atoms.calc=calc;batch=input_batch(calc,atoms,t['charge'],1)
        trace=[];last={}
        def objective(q):
            p,c,j,cj=kin.evaluate(q);atoms.set_positions(c)
            energy=float(atoms.get_potential_energy());forces=np.asarray(atoms.get_forces(),dtype=float)
            gradient=-np.einsum('mij,ij->m',cj,forces)
            trace.append({'energy_eV':energy,'q':q.tolist(),'gradient_max_eV_per_unit':float(max(abs(gradient)))})
            last.update(q=q.copy(),energy=energy,forces=forces,gradient=gradient,context=c)
            return energy,gradient
        objective(zero);initial_energy=last['energy']
        # Subtract one fixed energy constant to avoid premature precision loss
        # from atomic baseline energies in SLSQP's convergence criterion.
        def shifted(q):
            e,g=objective(q);return e-initial_energy,g
        result=minimize(shifted,zero,jac=True,method='SLSQP',bounds=[(-.20,.20)]*n,
            constraints=[{'type':'ineq','fun':lambda q:constraints(kin,q)[0],
                          'jac':lambda q:constraints(kin,q)[1]}],
            options={'maxiter':SETTINGS['maxiter'],'ftol':SETTINGS['ftol_eV']})
        if not np.array_equal(result.x,last['q']):objective(result.x)
        r={'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'cache_key':t['cache_key'],
           'manifest':record(manifest),'native_batch':batch,'q':result.x.tolist(),
           'iterations':int(result.nit),'objective_evaluations':len(trace),'scipy_success':bool(result.success),
           'termination':str(result.message),'source_context_energy_eV':initial_energy,
           'context_energy_eV':last['energy'],'trace':trace,'slurm_job_id':os.environ['SLURM_JOB_ID']}
        try:
            r['geometry_checks']=proposal_checks(geo,result.x,rows,cr)
            if last['energy']>initial_energy+1e-7:raise InvalidArtifact('context energy increased')
            # Tangent-cone projected gradient diagnostic in declared physical
            # coordinate units. This is a stationarity diagnostic, not curvature.
            values,cjac=constraints(kin,result.x);active=values<=4e-6
            A=list(cjac[active])
            for i,q in enumerate(result.x):
                if abs(q)>=.20-1e-5:
                    v=np.zeros(n);v[i]=-np.sign(q);A.append(v)
            A=np.array(A).reshape((-1,n));grad=last['gradient']
            proj=minimize(lambda v:(.5*np.sum((v+grad)**2),v+grad),np.zeros(n),jac=True,method='SLSQP',
                constraints=[{'type':'ineq','fun':lambda v:A@v,'jac':lambda v:A}],options={'ftol':1e-12,'maxiter':200})
            r['tangent_projected_gradient_max_eV_per_unit']=float(max(abs(proj.x)))
            r['stationarity_projection_success']=bool(proj.success)
            cp=dest/'context_proposed.xyz';write_xyz(cp,[(a[0],*v) for a,v in zip(rows,last['context'])])
            core=original.evaluate(result.x)[1];xp=dest/'core_proposed.xyz';core_rows=[(a[0],*v) for a,v in zip(cr,core)];write_xyz(xp,core_rows)
            ca=Atoms([a[0] for a in cr],positions=core,pbc=False);ca.info.update(charge=t['core']['charge'],spin=1);ca.calc=calc
            r.update(context_xyz=record(cp),core_xyz=record(xp),core_native_batch=input_batch(calc,ca,t['core']['charge'],1),
                     core_energy_eV=float(ca.get_potential_energy()),status='computed')
        except (ValueError,KeyError) as exc:r.update(status='unsupported_proposal',reason=str(exc))
        r['parameter_versions_unchanged']=versions=={k:p._version for k,p in model.named_parameters()}
        if not r['parameter_versions_unchanged']:raise InvalidArtifact('model parameters changed')
        r['wall_seconds']=time.monotonic()-ts;write_new(rp,r);results.append(record(rp))
        print(json.dumps({k:r[k] for k in ('task_id','status','iterations','objective_evaluations','wall_seconds')}),flush=True)
    receipt={'manifest':record(manifest),'status':'complete','results':results,'wall_seconds':time.monotonic()-began,
             'slurm_job_id':os.environ['SLURM_JOB_ID'],'allocated_cpus':int(os.environ['SLURM_CPUS_PER_TASK']),
             'peak_host_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'peak_cuda_allocated_bytes':torch.cuda.max_memory_allocated()}
    write_new(root/('collection_'+os.environ['SLURM_JOB_ID']+'.json'),receipt);return receipt


def prepare_dft(collection,output):
    col=read_json(collection);m=read_json(verify(col['manifest']));out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';shutil.copytree(Path(m['implementation'][Path(__file__).name]['path']).parent,impl);tasks=[];failures=[]
    for pin in col['results']:
        r=read_json(verify(pin));t=next(t for t in m['tasks'] if t['task_id']==r['task_id'])
        if r['status']!='computed':failures.append(pin);continue
        d=out/'tasks'/t['task_id'];d.mkdir(parents=True);xp=d/'core.xyz';xp.write_bytes(verify(r['core_xyz']).read_bytes())
        ip=d/'endpoint.inp';ip.write_text('! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3\n'+f"* xyzfile {t['core']['charge']} 1 core.xyz\n")
        tasks.append({'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'charge':t['core']['charge'],'multiplicity':1,
                      'xyz':record(xp),'input':record(ip),'output_path':str(d/'endpoint.out'),'proposal':pin,'original_core':t['core']})
    mm={'protocol_id':'native_r2scan3c_context_proposed_coordination_v1','tasks':tasks,'unscorable_proposals':failures,
        'agreement':m['agreement'],'orca':m['orca'],'execution_policy':m['execution_policy'],
        'execution_resources':{'mpi_ranks':16,'concurrent_tasks':4},'source_collection':record(collection),
        'implementation':{p.name:record(p) for p in impl.glob('*.py')},'compute_budget':None,'reference':None,'calibrated_decision':None}
    write_new(out/'manifest.json',mm)
    from affordable_workflow import dry_run
    return dry_run(out/'manifest.json')


def collect(manifest,output):
    m=read_json(manifest);rows=[];failed=[]
    for t in m['tasks']:
        try:
            e=endpoint(record(t['output_path']),record(t['output_path']+'.execution.json'),t['xyz'],t['input'])
            old=t['original_core'];o=endpoint(old['output'],old['receipt'],old['xyz'],old['input'])
            rows.append({'case':t['case'],'metal':t['metal'],'energy_hartree':e['energy_hartree'],
                         'original_energy_hartree':o['energy_hartree'],
                         'change_kcal_mol':(e['energy_hartree']-o['energy_hartree'])*HA_TO_KCAL,
                         'endpoint':e,'proposal':t['proposal']})
        except (ValueError,KeyError,OSError) as exc:failed.append({'task_id':t['task_id'],'reason':str(exc)})
    contrasts=[]
    for case in CASES:
        by={r['metal']:r for r in rows if r['case']==case}
        if set(by)!=set(('Ca','La')):continue
        contrasts.append({'case':case,'before_R_kcal_mol':(by['Ca']['original_energy_hartree']-by['La']['original_energy_hartree'])*HA_TO_KCAL,
                          'after_R_kcal_mol':(by['Ca']['energy_hartree']-by['La']['energy_hartree'])*HA_TO_KCAL})
    by={r['case']:r for r in contrasts};comparisons=[]
    for la,ca in [('4MAE','1H4I'),('1F6S','1GLG'),('6IP9','1GLG')]:
        if la in by and ca in by:
            comparisons.append({'La_like':la,'Ca_like':ca,**{s+'_gap_kcal_mol':by[la][s+'_R_kcal_mol']-by[ca][s+'_R_kcal_mol'] for s in ('before','after')}})
    result={'manifest':record(manifest),'rows':rows,'failed':failed,'contrasts':contrasts,'comparisons':comparisons,
            'calibrated_decision':None,'reference':None,'all_cases_consumed_development':True,'baseline_changed':False}
    write_new(output,result);return result


def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare')
    for name in ('static-manifest','pqq-collection','pqq-dft','agreement','output'):p.add_argument('--'+name,required=True)
    for command in ('validate','execute'):
        p=sub.add_parser(command);p.add_argument('--manifest',required=True)
    p=sub.add_parser('prepare-dft');p.add_argument('--collection',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('collect');p.add_argument('--manifest',required=True);p.add_argument('--output',required=True)
    a=vars(ap.parse_args());command=a.pop('command').replace('-','_');result=globals()[command](**a);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
