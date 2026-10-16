#!/usr/bin/env python3
"""Finite Mex interface check. --stage run|mm --config CONFIG --output FRESH.
Main config pins implementation,plan,static_result,classical_helper; workers=4.
MM stage is an internal pinned-CPU-interpreter subprocess, no origin evaluations.
"""
from concurrent.futures import ThreadPoolExecutor
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import socket
import sys
import time
import traceback

ROOT=Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
STEPS=(-.0001,-.00005,.00005,.0001)
FIXTURES=('external_Leu97_O','boundary_Leu83_CA')
COMPONENTS=('MACE_vacuum','GFN2_field_difference','MM_ledger','fixed_FF_charge_ALPB')


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path']);return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def write(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def inputs(cfg,out):
    verify(cfg['plan']);require(cfg['workers']==4,'Four-worker scope differs')
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Implementation differs')
    previous=read(cfg['static_result'])
    require(previous['status']=='complete' and previous['slurm_job_id']=='12685'
        and previous['new_MACE_calls_attempted']==previous['new_ALPB_calls_attempted']==2
        and previous['new_SCC_calls']==previous['new_MM_calls']==0,'Admitted static parent differs')
    old=module('mex_force_static',verify(previous['implementation']));origin_replay=out/'origin_replay';origin_replay.mkdir()
    sc,cells,solcfg,mcfg,native_helper,mapper,np=old.inputs(verify(previous['config']),origin_replay)
    field=read(sc['field_result']);local=module('mex_force_local',verify(field['implementation']))
    world=local.inputs(verify(field['config']))
    source=read(cells['La']['mapping']['inputs']);export=read(source['export'])
    helper=verify(cfg['classical_helper']);origins={}
    require(len(previous['rows'])==2 and {r['metal'] for r in previous['rows']}=={'La','Dy'},'Static denominator')
    for row in previous['rows']:
        require(row['status']=='complete','Static cell unavailable');a=read(row['assembled']);metal=row['metal'];c=cells[metal]
        require(a['source_ids']==c['ids'] and a['source_xyz']==c['xyz'] and a['physical_mapping']==c['mapping_pin']
             and a['MACE']==row['MACE'] and a['ALPB']==row['ALPB'] and a['classical_reused']==sc['classical'][metal]
             and a['GFN2_reused']==c['field_rows'],'Static source/component binding differs')
        require(np.array_equal([r['xyz_A'] for r in a['physical_atoms']],c['x']),'Origin coordinates differ')
        for p in (a['MACE'],a['ALPB'],a['classical_reused'],*a['GFN2_reused'].values()):verify(p)
        for g in [a['gradient_kcal_mol_A'],*a['component_gradients_kcal_mol_A'].values()]:
            g=np.asarray(g);require(g.shape==(2049,3) and np.isfinite(g).all(),'Origin gradient invalid')
        require(set(a['component_energies_kcal_mol'])==set(a['component_gradients_kcal_mol_A'])==set(COMPONENTS),'Origin components differ')
        ledger=read(c['classical']['ledger']);require(verify(ledger['reused_builder'])==helper,'Classical helper differs')
        origins[metal]={'assembled':a,'pin':row['assembled'],'MACE':read(a['MACE']),
                        'GFN2':{k:read(p) for k,p in a['GFN2_reused'].items()},'ledger':ledger}
    verify(previous['build']['binary'])
    require(np.array_equal(cells['La']['x'],cells['Dy']['x']) and np.array_equal(cells['La']['q'],cells['Dy']['q']),
            'Paired physical origins/charges differ')
    return old,cells,solcfg,mcfg,native_helper,mapper,np,local,world,source,export,origins,previous


def prepare(cfg,out,world):
    old,cells,solcfg,mcfg,native_helper,mapper,np,local,field_world,source,export,origins,previous=world
    c=cells['La'];ids=c['ids'];ix={s:i for i,s in enumerate(ids)};x=c['x']
    require(ix['A/97//O']==1019 and ix['A/83//CA']==801 and ix['A/83//C']==803 and ix['A/83//N']==799,'Fixture identity differs')
    bonds={frozenset(v) for v in read(export['bonds'])}
    for pair in [('A/97//C','A/97//O'),('A/83//C','A/83//CA'),('A/83//CA','A/83//N')]:
        require(frozenset(pair) in bonds,'Fixture frame bond missing')
    field=c['mapping']['field'];require(field[832]['id']=='A/97//O' and field[832]['charge_e']==-.5679
        and all(a['id']!='A/83//CA' for a in field),'Fixture field membership differs')
    cap=next(a for a in c['mapping']['core'] if a['qm_index']==205)
    require(cap['kind']=='cap' and cap['retained_source_id']=='A/83//C' and cap['omitted_source_id']=='A/83//CA'
        and cap['length_A']==1.09,'Fixture cap relation differs')
    def unit(v):
        n=float(np.linalg.norm(v));require(n>1e-8,'Degenerate source-selected direction');return v/n
    u=unit(x[801]-x[803]);w=x[799]-x[801]
    fixtures={FIXTURES[0]:{'index':1019,'source_id':'A/97//O','direction':unit(x[1019]-x[ix['A/97//C']]).tolist()},
              FIXTURES[1]:{'index':801,'source_id':'A/83//CA','direction':unit(w-u*np.dot(u,w)).tolist()}}
    native=read(export['atoms']);boundary=read(source['boundary_mapping']);trials=[];tasks=[]
    directory=out/'fixtures';directory.mkdir()
    for name in FIXTURES:
        f=fixtures[name]
        for h in STEPS:
            y=x.copy();y[f['index']]+=h*np.asarray(f['direction']);changed=np.flatnonzero(np.any(y!=x,axis=1))
            require(changed.tolist()==[f['index']] and abs(float(np.linalg.norm(y-x,axis=1).max())-abs(h))<1e-12,'Single physical displacement differs')
            for metal in ('La','Dy'):
                c=cells[metal];tid=name+'_'+metal+'_'+('%+.5f'%h).replace('+','p').replace('-','m').replace('.','_')
                d=directory/tid;d.mkdir();physical=copy.deepcopy(c['physical'])
                for a in physical:a['xyz_A']=y[ix[a['id']]].tolist()
                write(d/'physical.json',physical)
                atoms=copy.deepcopy(native)
                for a in atoms:a['xyz_A']=y[ix[a['id']]].tolist()
                write(d/'atoms_charges.json',atoms);ex=copy.deepcopy(export);ex.update(atoms=pin(d/'atoms_charges.json'),origin_export=source['export'])
                write(d/'EXPORT.json',ex)
                env=copy.deepcopy(c['mapping']['field'])
                for a in env:
                    if a['id']==f['source_id']:
                        a['xyz_A']=y[ix[a['id']]].tolist()
                        if 'source_xyz_A' in a:a['source_xyz_A']=a['xyz_A']
                        if 'source' in a:a['source']['xyz_A']=a['xyz_A']
                write(d/'environment_atoms.json',env)
                pc=verify(c['mapping']['pointcharges']).read_text().splitlines(keepends=True)
                if name==FIXTURES[0]:
                    tokens=pc[833].split();require(len(tokens)==4,'Original four-column field layout')
                    pc[833]=tokens[0]+' '+ ' '.join('%.17g'%v for v in y[1019])+'\n'
                (d/'environment.pc').write_text(''.join(pc))
                if name==FIXTURES[1]:require(pin(d/'environment.pc')['sha256']==c['mapping']['pointcharges']['sha256'],'Boundary field changed')
                bd=copy.deepcopy(boundary);bd['original_boundary_mapping']=source['boundary_mapping']
                bd['dipole_diagnostics_scope']='historical origin only; fixed charge-shift membership at perturbed geometry'
                for sid in bd['retained_physical_CA_coordinates']:bd['retained_physical_CA_coordinates'][sid]=y[ix[sid]].tolist()
                write(d/'boundary_mapping.json',bd)
                core=copy.deepcopy(c['mapping']['core']);qm=c['qm_positions'].copy()
                if name==FIXTURES[1]:
                    a=next(a for a in core if a['qm_index']==205);r=y[ix[a['retained_source_id']]];o=y[ix[a['omitted_source_id']]]
                    v=o-r;dist=float(np.linalg.norm(v));require(dist>0,'Coincident cap anchors');u=v/dist
                    jo=a['length_A']/dist*(np.eye(3)-np.outer(u,u));jr=np.eye(3)-jo;qm[205]=r+a['length_A']*u
                    a.update(xyz_A=qm[205].tolist(),retained_xyz_A=r.tolist(),omitted_xyz_A=o.tolist(),jacobian_retained=jr.tolist(),jacobian_omitted=jo.tolist())
                    require(np.flatnonzero(np.any(qm!=c['qm_positions'],axis=1)).tolist()==[205],'Boundary changed unexpected local atom')
                write(d/'core_mapping.json',core)
                lines=verify(c['xyz']).read_text().splitlines(keepends=True)
                if name==FIXTURES[1]:lines[207]='H '+' '.join('%.17g'%v for v in qm[205])+'\n'
                (d/'core.xyz').write_text(''.join(lines))
                if name==FIXTURES[0]:require(pin(d/'core.xyz')['sha256']==c['xyz']['sha256'],'External fixture changed local XYZ')
                symbols,actual_qm=local.xyz(pin(d/'core.xyz'));require(symbols==c['symbols'] and np.array_equal(actual_qm,qm),'Local serialization differs')
                endpoint=copy.deepcopy(source['configurations']['A']['endpoints'][metal]);endpoint['xyz']=pin(d/'core.xyz')
                conf={'endpoints':{metal:endpoint},'core_mapping':pin(d/'core_mapping.json'),'pointcharges':pin(d/'environment.pc')}
                inp=copy.deepcopy(source);inp.update(origin_inputs=c['mapping']['inputs'],export=pin(d/'EXPORT.json'),
                    environment_atoms=pin(d/'environment_atoms.json'),boundary_mapping=pin(d/'boundary_mapping.json'),configurations={'A':conf},
                    physical_source=pin(d/'atoms_charges.json'),physical_sources={'A':pin(d/'atoms_charges.json')},core_mapping={'A':pin(d/'core_mapping.json')})
                write(d/'INPUTS.json',inp);mp=mapper.prepare_mapping(d/'INPUTS.json','A',metal)
                require([a['id'] for a in mp['physical_atoms']]==ids and np.array_equal([a['xyz_A'] for a in mp['physical_atoms']],y),
                        'Derived mapping changed physical order/coordinates')
                write(d/'MAPPING.json',mp)
                fdesc=copy.deepcopy(read(origins[metal]['GFN2']['embedded']['task']['field']))
                fdesc.update(origin_field=origins[metal]['GFN2']['embedded']['task']['field'],source=pin(d/'environment.pc'),
                    source_mapping=pin(d/'environment_atoms.json'),boundary_mapping=pin(d/'boundary_mapping.json'))
                for row in fdesc['rows']:row['xyz_A']=y[ix[row['source_id']]].tolist()
                write(d/'FIELD.json',fdesc)
                trial={'task_id':tid,'fixture':name,'metal':metal,'h_A':h,'physical':pin(d/'physical.json'),'mapping':pin(d/'MAPPING.json'),
                       'xyz':pin(d/'core.xyz'),'field':pin(d/'FIELD.json'),'source_inputs':pin(d/'INPUTS.json'),'GFN2':{}}
                for kind in ('embedded','isolated'):
                    origin=origins[metal]['GFN2'][kind]
                    if name==FIXTURES[0] and kind=='isolated':
                        trial['GFN2'][kind]={'reused':origin['result'] if 'result' in origin else c['field_rows'][kind]};continue
                    gd=d/('GFN2_'+kind);gd.mkdir();task=copy.deepcopy(origin['task']);staged={}
                    for fname,p in task['staged_files'].items():
                        content=verify(p).read_bytes()
                        if fname=='core.xyz':content=(d/'core.xyz').read_bytes()
                        if fname=='external.pc' and name==FIXTURES[0]:
                            raw=content.decode().splitlines(keepends=True);tokens=raw[833].split()
                            require(len(tokens)==5,'Explicit GAM field layout differs')
                            raw[833]=' '.join([tokens[0],*('%.17g'%v for v in y[1019]),tokens[4]])+'\n';content=''.join(raw).encode()
                        (gd/fname).write_bytes(content);staged[fname]=pin(gd/fname)
                    if kind=='embedded' and name==FIXTURES[1]:
                        require(staged['external.pc']['sha256']==task['staged_files']['external.pc']['sha256'],'Boundary changed embedded field')
                    task.update(task_id=tid+'_'+kind,xyz=staged['core.xyz'],source_xyz=trial['xyz'],mapping=trial['mapping'],field=trial['field'],
                        source_pointcharges=pin(d/'environment.pc'),staged_files=staged,origin_task=origin['result'] if 'result' in origin else c['field_rows'][kind])
                    task.pop('preparation_task',None);task.pop('task_record',None)
                    write(gd/'TASK.json',task);trial['GFN2'][kind]={'task':pin(gd/'TASK.json')};tasks.append((task,gd))
                write(d/'TRIAL.json',trial);trial['record']=pin(d/'TRIAL.json');trials.append(trial)
    require(len(trials)==16 and len(tasks)==24,'Finite matrix count differs')
    for name in FIXTURES:
        for h in STEPS:
            pair=[t for t in trials if t['fixture']==name and t['h_A']==h]
            require(np.array_equal([a['xyz_A'] for a in read(pair[0]['physical'])],[a['xyz_A'] for a in read(pair[1]['physical'])]),'Paired trial coordinates differ')
    write(out/'PRECHECK.json',{'status':'passed','fixtures':fixtures,'trials':[t['record'] for t in trials],
          'declared_combined_cells':16,'new_MACE_calls_max':8,'new_GFN2_calls_max':24,'new_MM_configurations_max':16,
          'new_MM_group_queries_max':64,'new_ALPB_calls_max':16,'new_molecular_calls':0,'static_result':cfg['static_result']})
    return fixtures,trials,tasks


def mm_stage(cfg,out,result):
    require(pin(sys.executable)['sha256']==cfg['python']['sha256'],'MM interpreter binary differs')
    verify(cfg['common_helper']);sys.path.insert(0,str(ROOT/'scripts'))
    builder=module('mex_fd_classical',verify(cfg['classical_helper']))
    import numpy as np
    import openmm as mm
    from openmm import unit
    require(mm.__version__==cfg['openmm_version'],'Original OpenMM version differs')
    platform=mm.Platform.getPlatformByName('Reference');prepared={};trial_data={}
    require(len(cfg['trials'])==16,'MM denominator differs')
    # All imports, systems and coordinates are admitted before the first energy.
    for metal in ('La','Dy'):
        ledger=read(cfg['ledgers'][metal]);require(ledger['reused_builder']==cfg['classical_helper'],'MM builder lineage')
        artifacts={k:read(p) for k,p in ledger['artifacts'].items()}
        system,counts=builder.build_system(ledger,artifacts,metal)
        require(system.getNumParticles()==2049 and system.getNumConstraints()==0,'MM physical inventory/constraints')
        ids=[cfg['target_source_id'] if a['id']=='target_EF3' else a['id'] for a in artifacts['particles.json']]
        require(len(ids)==len(set(ids))==2049,'MM particle identity differs')
        prepared[metal]=(system,ids,counts)
    for p in cfg['trials']:
        t=read(p);physical=read(t['physical']);by={a['id']:a for a in physical};system,ids,counts=prepared[t['metal']]
        require(len(by)==2049 and set(by)==set(ids),'MM trial/source ID join differs')
        x=np.asarray([by[s]['xyz_A'] for s in ids]);require(x.shape==(2049,3) and np.isfinite(x).all(),'MM trial coordinate shape')
        trial_data[t['task_id']]=(t,x)
    result.update(runtime={'python':sys.executable,'numpy':np.__version__,'openmm':mm.__version__,
                          'openmm_full_version':mm.version.full_version,'platform':platform.getName()},rows=[],
                  configurations_attempted=0,group_queries_attempted=0)
    contexts={}
    for metal,(system,ids,counts) in prepared.items():
        integrator=mm.VerletIntegrator(.001);contexts[metal]=(mm.Context(system,integrator,platform),integrator)
    for tid,(t,x) in trial_data.items():
        row={'task_id':tid,'metal':t['metal'],'status':'unavailable','trial':next(p for p in cfg['trials'] if read(p)['task_id']==tid)}
        try:
            context=contexts[t['metal']][0];context.setPositions(x*.1*unit.nanometer);energies={}
            result['configurations_attempted']+=1
            for i,name in enumerate(builder.GROUPS):
                result['group_queries_attempted']+=1
                energies[name]=context.getState(getEnergy=True,groups=1<<i).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
            require(all(math.isfinite(v) for v in energies.values()),'Nonfinite MM energy')
            row.update(status='complete',components_kcal_mol=energies,total_kcal_mol=math.fsum(energies.values()))
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
        result['rows'].append(row)
    result['status']='complete' if all(r['status']=='complete' for r in result['rows']) else 'unavailable'


def compare(fixtures,trials,rows,world):
    old,cells,solcfg,mcfg,helper,mapper,np,local,fw,source,export,origins,previous=world
    by={r['task_id']:r for r in rows};summary=[];derived=[]
    for name,f in fixtures.items():
        i=f['index'];v=np.asarray(f['direction']);states={}
        for metal in ('La','Dy'):
            a=origins[metal]['assembled'];c=cells[metal]
            analytic={k:float(np.dot(np.asarray(a['component_gradients_kcal_mol_A'][k])[i],v)) for k in COMPONENTS}
            analytic['total']=float(np.dot(np.asarray(a['gradient_kcal_mol_A'])[i],v))
            analytic.update({'GFN2_'+kind:float(np.dot(np.asarray(origins[metal]['GFN2'][kind]['mapped_physical_gradient']['gradient_kcal_mol_A'])[i],v)) for kind in ('embedded','isolated')})
            ci=c['classical_order'][i]
            group_a={k:float(np.dot(np.asarray(g)[ci],v)) for k,g in c['classical']['component_gradients_kcal_mol_A'].items()}
            row={'fixture':name,'metal':metal,'status':'unavailable','analytic_kcal_mol_A':analytic,'MM_group_analytic_kcal_mol_A':group_a,'steps':[]}
            selected={t['h_A']:by[t['task_id']] for t in trials if t['fixture']==name and t['metal']==metal}
            if all(r['status']=='complete' for r in selected.values()):
                for h in (.0001,.00005):
                    p,m=selected[h]['native_energies'],selected[-h]['native_energies'];den=2*h
                    fd={'MACE_vacuum':(p['MACE_eV']-m['MACE_eV'])*old.EV/den,
                        'GFN2_embedded':(p['GFN2_embedded_Eh']-m['GFN2_embedded_Eh'])*old.HA/den,
                        'GFN2_isolated':(p['GFN2_isolated_Eh']-m['GFN2_isolated_Eh'])*old.HA/den,
                        'fixed_FF_charge_ALPB':(p['ALPB_Eh']-m['ALPB_Eh'])*old.HA/den}
                    groups={k:(p['MM_groups_kcal_mol'][k]-m['MM_groups_kcal_mol'][k])/den for k in group_a}
                    fd['MM_ledger']=math.fsum(groups.values());fd['GFN2_field_difference']=math.fsum([fd['GFN2_embedded'],-fd['GFN2_isolated']])
                    fd['total']=math.fsum(fd[k] for k in COMPONENTS)
                    checks={k:{'analytic':analytic[k],'finite_difference':fd[k],'residual':fd[k]-analytic[k],
                        'tolerance':.2+.05*abs(analytic[k]),'passed':abs(fd[k]-analytic[k])<=.2+.05*abs(analytic[k])} for k in (*COMPONENTS,'total')}
                    group_checks={k:{'analytic':group_a[k],'finite_difference':groups[k],'residual':groups[k]-group_a[k],
                        'tolerance':.001+.0001*abs(group_a[k]),'passed':abs(groups[k]-group_a[k])<=.001+.0001*abs(group_a[k])} for k in groups}
                    row['steps'].append({'h_A':h,'finite_difference_kcal_mol_A':fd,'checks':checks,'MM_group_checks':group_checks,
                        'raw_GFN2_residuals':{k:fd[k]-analytic[k] for k in ('GFN2_embedded','GFN2_isolated')}})
                step_checks={k:{'difference':row['steps'][0]['finite_difference_kcal_mol_A'][k]-row['steps'][1]['finite_difference_kcal_mol_A'][k],
                     'tolerance':.2+.05*abs(analytic[k])} for k in (*COMPONENTS,'total')}
                for z in step_checks.values():z['passed']=abs(z['difference'])<=z['tolerance']
                row['step_agreement']=step_checks
                row['unchanged_local_exact_zero']=name!=FIXTURES[0] or all(
                    step['finite_difference_kcal_mol_A'][k]==0 and abs(analytic[k])<=1e-10
                    for step in row['steps'] for k in ('MACE_vacuum','GFN2_isolated'))
                row['status']='passed' if (all(z['passed'] for s in row['steps'] for group in ('checks','MM_group_checks') for z in s[group].values())
                    and all(z['passed'] for z in step_checks.values()) and row['unchanged_local_exact_zero']) else 'failed'
            row['native_PC_rounding_projection_bound_kcal_mol_A']=float(5e-9*old.HA/old.BOHR*np.abs(v).sum()) if name==FIXTURES[0] else 0.
            summary.append(row);states[metal]=row
        for state,weights in (('common',(.5,.5)),('Dy_minus_La',(-1.,1.))):
            a=math.fsum(weights[j]*states[m]['analytic_kcal_mol_A']['total'] for j,m in enumerate(('La','Dy')))
            r={'fixture':name,'state':state,'analytic_total_kcal_mol_A':a,'status':'unavailable'}
            if all(states[m]['status']!='unavailable' for m in states):
                vals=[math.fsum(weights[j]*states[m]['steps'][i]['finite_difference_kcal_mol_A']['total'] for j,m in enumerate(('La','Dy'))) for i in range(2)]
                tol=.2+.05*abs(a);r.update(finite_difference_totals_kcal_mol_A=vals,residuals=[v-a for v in vals],
                    step_difference=vals[0]-vals[1],tolerance=tol,status='passed' if all(abs(v-a)<=tol for v in vals) and abs(vals[0]-vals[1])<=tol else 'failed')
            derived.append(r)
    return {'directions':summary,'derived_total_states':derived,'all_passed':all(r['status']=='passed' for r in summary+derived),
            'interpretation':'Two local derivative interfaces only; no equilibrium, affinity, solvent-accuracy or classifier qualification'}


def run(cfg,out,result):
    require(int(os.environ['SLURM_CPUS_PER_TASK'])>=4,'Four one-thread workers need allocated CPUs')
    world=inputs(cfg,out);old,cells,solcfg,mcfg,helper,mapper,np,local,fw,source,export,origins,previous=world
    fixtures,trials,tasks=prepare(cfg,out,world);result['precheck']=pin(out/'PRECHECK.json');result['origin_replay']=pin(out/'origin_replay/PRECHECK.json')
    backend,env=fw[6],fw[9];mmcfg={'implementation':cfg['implementation'],'plan':cfg['plan'],'static_result':cfg['static_result'],
        'python':fw[1]['python'],'common_helper':fw[1]['common_helper'],'classical_helper':cfg['classical_helper'],
        'openmm_version':export['openmm_version'],'ledgers':{m:cells[m]['classical']['ledger'] for m in cells},
        'target_source_id':cells['La']['target_id'],'trials':[t['record'] for t in trials]}
    write(out/'MM_CONFIG.json',mmcfg);mmdir=out/'mm_driver';mmdir.mkdir()
    command=[str(verify(mmcfg['python'])),str(verify(cfg['implementation'])),'--stage','mm','--config',str(out/'MM_CONFIG.json'),'--output',str(out/'MM')]
    execution=backend.call(command,mmdir,env);result['MM_execution']=pin(mmdir/'EXECUTION.json')
    mmresult=read(pin(out/'MM/RESULT.json'));result['MM_result']=pin(out/'MM/RESULT.json')
    result['new_MM_configurations_attempted']=mmresult.get('configurations_attempted',0)
    result['new_MM_group_queries_attempted']=mmresult.get('group_queries_attempted',0)
    require(execution['returncode']==0 and execution['owned_process_group_empty'] and mmresult['status']=='complete'
        and mmresult['config']==pin(out/'MM_CONFIG.json'),'MM stage unavailable')
    mmby={r['task_id']:r for r in mmresult['rows']};require(len(mmby)==16,'MM result denominator')
    native=helper.NativeMACE(mcfg);ml={}
    for t in trials:
        if t['fixture']==FIXTURES[0]:ml[t['task_id']]={'status':'reused','value':origins[t['metal']]['MACE'],'receipt':origins[t['metal']]['assembled']['MACE']};continue
        row={'status':'unavailable'};ml[t['task_id']]=row;d=Path(t['record']['path']).parent
        try:
            c=dict(cells[t['metal']]);c.update(qm_positions=np.asarray(local.xyz(t['xyz'])[1]),xyz=t['xyz'])
            result['new_MACE_calls_attempted']+=1;value=old.mace_cell(native,c,t['metal'])
            write(d/'MACE.json',value);row.update(status='complete',value=value,receipt=pin(d/'MACE.json'))
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
    result.update(native_model_calls=native.calls,native_device=native.device,model=mcfg['model'],versions=native.versions)
    def gfn_cell(item):
        task,d=item;row={'task':task,'status':'unavailable','attempted':False}
        try:
            for p in [fw[4]['executable'],fw[4]['parameters'],*fw[4]['runtime_libraries']]:verify(p)
            for p in task['staged_files'].values():verify(p)
            row['attempted']=True;execution=backend.call(task['command'],d,env);row['execution']=pin(d/'EXECUTION.json')
            for p in task['staged_files'].values():verify(p)
            row.update(local.admit(d,task,execution,fw),status='accepted')
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
        row['artifacts']={p.name:pin(p) for p in d.iterdir() if p.is_file()};write(d/'RESULT.json',row);row['result']=pin(d/'RESULT.json');return row
    with ThreadPoolExecutor(max_workers=4) as pool:gfns=list(pool.map(gfn_cell,tasks))
    gby={r['task']['task_id']:r for r in gfns};result['GFN2_results']=[r['result'] for r in gfns]
    result['new_GFN2_calls_attempted']=sum(r['attempted'] for r in gfns)
    executions=[read(r['execution']) for r in gfns if 'execution' in r]
    result['observed_GFN2_child_CPU_seconds']=math.fsum(r['cpu_user_seconds']+r['cpu_system_seconds'] for r in executions)
    alpb_env=dict(env,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=solcfg['library_directory'])
    binary=verify(previous['build']['binary']);rows=[];result['rows']=rows
    for t in trials:
        d=Path(t['record']['path']).parent;row={'task_id':t['task_id'],'fixture':t['fixture'],'metal':t['metal'],'h_A':t['h_A'],'status':'unavailable','trial':t['record']};rows.append(row)
        try:
            c=dict(cells[t['metal']]);physical=read(t['physical']);c.update(physical=physical,x=np.asarray([a['xyz_A'] for a in physical]),mapping_pin=t['mapping'])
            result['new_ALPB_calls_attempted']+=1;sol=old.solvent_cell(binary,c,d/'solvent',alpb_env,np);row['ALPB']=pin(d/'solvent/RESULT.json')
            m=ml[t['task_id']];require(m['status'] in ('complete','reused'),'MACE unavailable');row['MACE']=m['receipt'];row['MACE_reused']=m['status']=='reused'
            g={}
            for kind in ('embedded','isolated'):
                entry=t['GFN2'][kind]
                if 'reused' in entry:g[kind]=read(entry['reused']);row['GFN2_'+kind]={'result':entry['reused'],'reused':True}
                else:g[kind]=gby[read(entry['task'])['task_id']];row['GFN2_'+kind]={'result':g[kind]['result'],'reused':False}
                require(g[kind]['status']=='accepted','GFN2 component unavailable')
            mmrow=mmby[t['task_id']];require(mmrow['status']=='complete','MM component unavailable')
            values={'MACE_eV':m['value']['energy_eV'],'GFN2_embedded_Eh':g['embedded']['energy_hartree'],
                'GFN2_isolated_Eh':g['isolated']['energy_hartree'],'MM_groups_kcal_mol':mmrow['components_kcal_mol'],'ALPB_Eh':sol['total_hartree']}
            row.update(status='complete',native_energies=values,MM_result=result['MM_result'])
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
        write(d/'CELL.json',row);row['cell']=pin(d/'CELL.json')
    result['comparison']=compare(fixtures,trials,rows,world)
    require(result['new_MACE_calls_attempted']<=8 and result['new_GFN2_calls_attempted']<=24
        and result['new_MM_configurations_attempted']<=16 and result['new_MM_group_queries_attempted']<=64
        and result['new_ALPB_calls_attempted']<=16,'Declared call ceiling exceeded')
    result['status']='passed' if len(rows)==16 and all(r['status']=='complete' for r in rows) and result['comparison']['all_passed'] else 'failed'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=('run','mm'),default='run')
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);cfg=json.loads(a.config.read_text());start=time.monotonic();cpu=time.process_time()
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Implementation differs');verify(cfg['plan'])
    result={'status':'running','stage':a.stage,'implementation':pin(__file__),'config':pin(a.config),'plan':cfg['plan'],
        'static_result':cfg['static_result'],'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),
        'new_MACE_calls_attempted':0,'new_GFN2_calls_attempted':0,'new_ALPB_calls_attempted':0,
        'new_MM_configurations_attempted':0,'new_MM_group_queries_attempted':0,'optimizer_steps':0,'new_DFT_calls':0,'new_wholeprotein_SCC_calls':0}
    try:
        if a.stage=='mm':mm_stage(cfg,out,result)
        else:run(cfg,out,result)
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-start;result['process_CPU_seconds']=time.process_time()-cpu
    if a.stage=='run':
        result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ['SLURM_CPUS_PER_TASK'])
        result['application_allocated_GPU_seconds']=result['elapsed_seconds']
    write(out/'RESULT.json',result)
    return 0 if result['status'] in ('passed','complete') else 1


if __name__=='__main__':raise SystemExit(main())
