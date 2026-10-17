#!/usr/bin/env python3
"""Four fixed LBT starts: native MM class and complete Phe1-CA term attribution."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import traceback
import xml.etree.ElementTree as ET

PROTOCOL = 'lbt_four_start_MM_Phe1CA_term_audit_v1'
AUDIT_SHA = '46d01517a6bb3c5dd9faf37ebe91724d433c314c50341d0aeacb60a141541549'
MM_SHA = '686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661'
BUILDER_SHA = '2a645bfcf766916665e2b9090a9a0c71f631f518d6a1fe460e3a93754473a5d1'
GROUPS = ('retained_bonded', 'MM_LJ', 'MM_Coulomb', 'QM_MM_LJ')
CLASSES = ('HarmonicBondForce', 'HarmonicAngleForce', 'PeriodicTorsionForce', 'CMAPTorsionForce')
TERM_API = {'HarmonicBondForce': ('Bonds','getNumBonds','getBondParameters','addBond',2),
            'HarmonicAngleForce': ('Angles','getNumAngles','getAngleParameters','addAngle',3),
            'PeriodicTorsionForce': ('Torsions','getNumTorsions','getTorsionParameters','addTorsion',4),
            'CMAPTorsionForce': ('Torsions','getNumTorsions','getTorsionParameters','addTorsion',8)}
SETTINGS = {'states': [['7CCO','La'],['7CCO','Dy'],['7CCN','La'],['7CCN','Dy']],
    'target_id':'A/1//CA', 'target_index':3, 'physical_atoms':256, 'region':'baseline', 'platform':'Reference',
    'energy_absolute_kcal_mol':1e-8, 'force_component_absolute_kcal_mol_A':1e-5,
    'local_term_counts':dict(zip(CLASSES,(4,14,34,0))), 'term_batch_size':31,
    'maximum_contexts_per_state':4, 'maximum_queries_per_state':80}


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])==record,'Changed pinned artifact: '+str(record)); return Path(record['path'])


def read(record): return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value,handle,indent=2,allow_nan=False); handle.write('\n')
    return pin(path)


def force_xml(mm, force): return ET.fromstring(mm.XmlSerializer.serialize(force))


def term_attributes(mm, force):
    return [dict(x.attrib) for x in force_xml(mm,force).find(TERM_API[type(force).__name__][0])]


def term_indices(force, index):
    kind=type(force).__name__; api=TERM_API[kind]; p=getattr(force,api[2])(index)
    return list(map(int,p[1:] if kind=='CMAPTorsionForce' else p[:api[4]]))


def input_state(cfg, audit, state, row, np, mm):
    source,metal=state; saved=next(r for r in audit['rows'] if [r['source_id'],r['metal']]==state)
    row['input_pins']={'force_audit_result':cfg['force_audit_result'],'POINT':saved['hybrid_point']['point'],
                       'force_table':saved['full_tables']}
    require(saved['status']=='available','Saved force state unavailable')
    point=read(saved['hybrid_point']['point']); receipts=point['component_receipts']; archived=read(receipts['MM'])
    execution=read(receipts['MM_execution']); mc=read(archived['config']); mapping=read(point['mapping'])
    table=read(saved['full_tables']); task=read(point['task']); ledger_pin=point['source_task']['ledger']; ledger=read(ledger_pin)
    require(point['status']=='complete' and point['state']['source']==source and point['state']['metal']==metal
            and point['state']['region']=='baseline' and point['source_task']['source_id']==source
            and point['source_task']['metal']==metal and point['source_task']['region']=='baseline', 'POINT state differs')
    require(task['ledger']==ledger_pin and task['partition']==point['source_task']['partition']==ledger['partition'], 'Task/ledger differs')
    partition=read(ledger['partition']); require(partition['source_id']==source and partition['region_name']=='baseline','Partition differs')
    require(archived['status']=='complete' and archived['stage']=='mm' and archived['group_order']==list(GROUPS)
            and archived['platform']=='Reference' and archived['openmm_version']==cfg['openmm_version']==mc['version']
            and archived['system']==mc['system'] and archived['positions']==mc['positions']==point['positions'], 'Archived MM differs')
    require(archived['implementation']==mc['implementation'] and mc['implementation']['sha256']==MM_SHA
            and mc['builder']['sha256']==BUILDER_SHA and mc['python']==cfg['python'], 'MM implementation/interpreter differs')
    for name in ('implementation','builder','python'): verify(mc[name])
    expected=[str(verify(mc['python'])),str(verify(mc['implementation'])),'--stage','mm','--config',str(verify(archived['config'])),
              '--output',str(verify(receipts['MM']).parent)]
    require(execution['returncode']==0 and execution['owned_process_group_empty'] and execution['command']==expected,
            'Archived MM execution differs')
    for name in ('stdout','stderr'): verify(execution[name])
    physical=mapping['physical_atoms']; ids=[a['id'] for a in physical]; target=SETTINGS['target_index']
    x=np.asarray(read(point['positions'])['positions_A'],dtype=float)
    require(x.shape==(256,3) and np.isfinite(x).all() and np.array_equal(x,[a['xyz_A'] for a in physical])
            and len(set(ids))==256 and ids[target]==SETTINGS['target_id'], 'Physical order/coordinates differ')
    require([a['id'] for a in table['atom_rows']]==ids and all(a['index']==i for i,a in enumerate(table['atom_rows']))
            and table['atom_rows'][target]['category']=='MM_other' and table['atom_rows'][target]['boundary_role'] is None,
            'Force table/target identity differs')
    artifacts={k:read(v) for k,v in ledger['artifacts'].items()}; particles=artifacts['particles.json']
    require([a['id'] for a in particles]==ids and particles[target]['region']=='MM','Ledger physical/target identity differs')
    verify(ledger['forcefield']); parent=mm.XmlSerializer.deserialize(verify(ledger['native_system']).read_text())
    system=mm.XmlSerializer.deserialize(verify(mc['system']).read_text())
    require(system.getNumParticles()==256 and system.getNumConstraints()==0,'MM inventory/constraints differ')
    forces=list(system.getForces()); require(len(forces)==8,'Expected four bonded and four nonbonded Force objects')
    by_class={}; native={}
    for kind in CLASSES:
        a=[f for f in forces if type(f).__name__==kind]; b=[f for f in parent.getForces() if type(f).__name__==kind]
        require(len(a)==len(b)==1 and a[0].getForceGroup()==0 and a[0].getName()=='retained_'+kind,'Bonded Force identity differs')
        by_class[kind]=a[0]; native[kind]=b[0]
    expected_nonbonded=[('MM_LJ','NonbondedForce',1),('MM_Coulomb','NonbondedForce',2),
        ('QM_MM_LJ_generic','CustomNonbondedForce',3),('QM_MM_LJ_native_exceptions','CustomBondForce',3)]
    require([(f.getName(),type(f).__name__,f.getForceGroup()) for f in forces if type(f).__name__ not in CLASSES]==expected_nonbonded,
            'Nonbonded ledger/group identity differs')
    decisions=read(ledger['bonded_decisions']); selected=artifacts['bonded_terms.json']; local=[]
    for kind in CLASSES:
        indices=sorted(t['index'] for t in selected if t['force']==kind)
        require(len(indices)==len(set(indices)),'Duplicate selected native term')
        raw=term_attributes(mm,native[kind]); actual=term_attributes(mm,by_class[kind])
        require(actual==[raw[i] for i in indices],'Native-to-retained term parameters differ: '+kind)
        if kind=='CMAPTorsionForce':
            require(ET.tostring(force_xml(mm,native[kind]).find('Maps'))==ET.tostring(force_xml(mm,by_class[kind]).find('Maps')),
                    'CMAP tables differ')
        count=0
        for local_index,native_index in enumerate(indices):
            atoms=term_indices(by_class[kind],local_index)
            if target not in atoms: continue
            d=next(d for d in decisions if d['force']==kind and d['index']==native_index)
            require(d['keep'] and d['reason']=='wholly_MM' and d['atoms']==atoms
                    and d['source_ids']==[ids[i] for i in atoms], 'Target local ledger decision differs')
            local.append({'force_class':kind,'retained_index':local_index,'native_index':native_index,
                'atoms':atoms,'atom_ids':[ids[i] for i in atoms],'reason':d['reason'],
                'native_parameter_XML':actual[local_index]});count+=1
        require(count==SETTINGS['local_term_counts'][kind],'Declared target-term denominator differs: '+kind)
    require(len(local)==52,'All52 target terms required')
    archived_force=-np.asarray(archived['gradient_kcal_mol_A'],float)
    require(archived_force.shape==(256,3) and np.isfinite(archived_force).all(),'Archived full force unavailable')
    table_force=np.asarray([a['components']['hybrid']['MM_ledger'] for a in table['atom_rows']],float)
    require(np.array_equal(archived_force,table_force),'Force-table MM vector differs')
    row['provenance']={'audit_row':saved,'POINT':saved['hybrid_point']['point'],'MM_result':receipts['MM'],
        'MM_execution':receipts['MM_execution'],'MM_config':archived['config'],'system':mc['system'],'positions':point['positions'],
        'mapping':point['mapping'],'task':point['task'],'ledger':ledger_pin,'native_system':ledger['native_system'],
        'forcefield':ledger['forcefield'],'bonded_decisions':ledger['bonded_decisions'],'artifacts':ledger['artifacts']}
    row['input_status']='passed';row['local_term_denominator']=len(local)
    return system,x,archived,by_class,local,ids


def query_system(system,x,queries,label,out,row,np,mm,unit):
    require(row['contexts_attempted']<4,'Context ceiling exceeded');row['contexts_attempted']+=1
    directory=out/label;directory.mkdir(); xml=directory/'SYSTEM.xml';xml.write_text(mm.XmlSerializer.serialize(system))
    row['diagnostic_systems'][label]=pin(xml)
    integrator=mm.VerletIntegrator(.001); context=None; results={}
    try:
        context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'))
        context.setPositions(x*.1*unit.nanometer)
        for name,mask in queries:
            require(row['queries_attempted']<80,'Query ceiling exceeded');row['queries_attempted']+=1
            state=context.getState(getEnergy=True,getForces=True,groups=mask)
            energy=state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
            force=np.asarray(state.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom))
            require(math.isfinite(energy) and force.shape==(256,3) and np.isfinite(force).all(),'Nonfinite diagnostic query')
            value={'energy_kcal_mol':energy,'force_kcal_mol_A':force.tolist()}
            record=write(directory/(name+'.json'),dict(value,system=pin(xml),group_mask=mask,physical_coordinates_unchanged=True))
            row['queries'].append({'label':label,'name':name,'groups':mask,'result':record});results[name]=value
        return results
    finally:
        if context is not None:del context
        del integrator


def close(label,observed,expected,row,np,energy=True):
    residual=float(np.max(np.abs(np.asarray(observed['force_kcal_mol_A'])-np.asarray(expected['force_kcal_mol_A']))))
    check={'name':label,'maximum_force_component_residual_kcal_mol_A':residual,
           'force_tolerance_kcal_mol_A':SETTINGS['force_component_absolute_kcal_mol_A']}
    passed=residual<=SETTINGS['force_component_absolute_kcal_mol_A']
    if energy:
        delta=abs(observed['energy_kcal_mol']-expected['energy_kcal_mol'])
        check.update(absolute_energy_residual_kcal_mol=delta,energy_tolerance_kcal_mol=SETTINGS['energy_absolute_kcal_mol'])
        passed=passed and delta<=SETTINGS['energy_absolute_kcal_mol']
    check['passed']=bool(passed);row['checks'].append(check);require(passed,'Numerical closure failed: '+label)


def subset_force(force,indices,name,group,mm):
    kind=type(force).__name__; require(kind in CLASSES[:3],'Only Bond/Angle/Periodic containers are split')
    api=TERM_API[kind];new=type(force)();new.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
    new.setName(name);new.setForceGroup(group)
    for i in indices:getattr(new,api[3])(*getattr(force,api[2])(i))
    original=term_attributes(mm,force)
    require(term_attributes(mm,new)==[original[i] for i in indices],'Split term parameters differ')
    return new


def split_system(original,batch,mm):
    system=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(original))
    for i in reversed(range(system.getNumForces())):system.removeForce(i)
    for force in original.getForces():
        kind=type(force).__name__; chosen=[(g,t) for g,t in enumerate(batch) if t['force_class']==kind]
        if not chosen:
            cloned=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(force));cloned.setForceGroup(31);system.addForce(cloned);continue
        n=getattr(force,TERM_API[kind][1])(); indices=[t['retained_index'] for _,t in chosen]
        require(len(indices)==len(set(indices)),'Duplicate selected batch term')
        remainder=[i for i in range(n) if i not in indices]
        system.addForce(subset_force(force,remainder,force.getName()+'_remainder',31,mm))
        for group,term in chosen:
            system.addForce(subset_force(force,[term['retained_index']],force.getName()+'_term_'+str(term['retained_index']),group,mm))
        require(sorted(indices+remainder)==list(range(n)),'Terms not represented exactly once')
    require(system.getNumParticles()==original.getNumParticles() and system.getNumConstraints()==0,'Split physical system differs')
    return system


def term_description(term,force,x,np,unit):
    kind=term['force_class'];p=getattr(force,TERM_API[kind][2])(term['retained_index']);atoms=term['atoms']
    if kind=='HarmonicBondForce':
        return {'parameter_units':{'d':'nm','k':'kJ/mol/nm^2','p1':'physical_index','p2':'physical_index'},
                'actual_bond_length_A':float(np.linalg.norm(x[atoms[1]]-x[atoms[0]])),
                'reference_bond_length_A':p[2].value_in_unit(unit.angstrom),
                'stiffness_kcal_mol_A2':p[3].value_in_unit(unit.kilocalorie_per_mole/unit.angstrom**2)}
    if kind=='HarmonicAngleForce':
        u=x[atoms[0]]-x[atoms[1]];v=x[atoms[2]]-x[atoms[1]];den=float(np.linalg.norm(u)*np.linalg.norm(v))
        require(den>0,'Singular angle');angle=math.acos(max(-1.,min(1.,float(np.dot(u,v))/den)))
        return {'parameter_units':{'a':'radian','k':'kJ/mol/radian^2','p1':'physical_index','p2':'physical_index','p3':'physical_index'},
                'actual_angle_radians':angle,'actual_angle_degrees':math.degrees(angle),
                'reference_angle_radians':p[3].value_in_unit(unit.radian),
                'stiffness_kcal_mol_radian2':p[4].value_in_unit(unit.kilocalorie_per_mole/unit.radian**2)}
    return {'parameter_units':{'periodicity':'integer','phase':'radian','k':'kJ/mol','p1':'physical_index',
                'p2':'physical_index','p3':'physical_index','p4':'physical_index'},
            'periodicity':int(p[4]),'phase_radians':p[5].value_in_unit(unit.radian),
            'amplitude_kcal_mol':p[6].value_in_unit(unit.kilocalorie_per_mole)}


def run_state(cfg,audit,state,out,row,np,mm,unit):
    row['active_stage']='input';system,x,archive,by_class,terms,ids=input_state(cfg,audit,state,row,np,mm)
    expected={name:{'energy_kcal_mol':archive['component_energies_kcal_mol'][name],
        'force_kcal_mol_A':(-np.asarray(archive['component_gradients_kcal_mol_A'][name])).tolist()} for name in GROUPS}
    expected_total={'energy_kcal_mol':archive['energy_kcal_mol'],'force_kcal_mol_A':(-np.asarray(archive['gradient_kcal_mol_A'])).tolist()}
    row['active_stage']='replay';baseline=query_system(system,x,[(n,1<<i) for i,n in enumerate(GROUPS)]+[('total',-1)],'replay',out,row,np,mm,unit)
    for name in GROUPS:close('archived_'+name,baseline[name],expected[name],row,np)
    close('archived_total',baseline['total'],expected_total,row,np);row['replay_status']='passed'
    row['active_stage']='class';classified=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))
    for f in classified.getForces():
        if type(f).__name__ in CLASSES:f.setForceGroup(4+CLASSES.index(type(f).__name__))
    classes=query_system(classified,x,[(name,1<<(4+i)) for i,name in enumerate(CLASSES)]+[('total',-1)],'classes',out,row,np,mm,unit)
    summed={'energy_kcal_mol':math.fsum(classes[n]['energy_kcal_mol'] for n in CLASSES),
        'force_kcal_mol_A':sum(np.asarray(classes[n]['force_kcal_mol_A']) for n in CLASSES).tolist()}
    close('bonded_class_sum',summed,baseline['retained_bonded'],row,np)
    close('class_system_total',classes['total'],baseline['total'],row,np);row['class_status']='passed'
    row['full_class_force_table']=write(out/'CLASS_FORCES.json',{'physical_ids':ids,'positions':row['provenance']['positions'],
        'force_units':'kcal/mol/angstrom','classes':{n:classes[n] for n in CLASSES},'all_classes_retained':True})
    row['active_stage']='local';target=SETTINGS['target_index'];local_values=[]
    for batch_index,start in enumerate(range(0,len(terms),31)):
        batch=terms[start:start+31];split=split_system(system,batch,mm)
        values=query_system(split,x,[('term_'+str(start+i),1<<i) for i in range(len(batch))]+[('total',-1)],
            'batch_'+str(batch_index),out,row,np,mm,unit)
        close('batch_'+str(batch_index)+'_total',values['total'],baseline['total'],row,np)
        for i,t in enumerate(batch):
            q=values['term_'+str(start+i)];local_values.append(dict(t,**term_description(t,by_class[t['force_class']],x,np,unit),
                energy_kcal_mol=q['energy_kcal_mol'],target_force_kcal_mol_A=q['force_kcal_mol_A'][target]))
        row['local_terms']=list(local_values)
    require(len(local_values)==52,'Incomplete local term denominator')
    for kind in CLASSES:
        rows=[v for v in local_values if v['force_class']==kind]
        summed=sum((np.asarray(v['target_force_kcal_mol_A']) for v in rows),np.zeros(3))
        close('target_local_sum_'+kind,{'force_kcal_mol_A':summed.tolist()},
              {'force_kcal_mol_A':classes[kind]['force_kcal_mol_A'][target]},row,np,energy=False)
    row['local_status']='passed';row['status']='available';row['active_stage']=None
    row['target_summary']={'id':ids[target],'index':target,'archived_force_kcal_mol_A':expected_total['force_kcal_mol_A'][target],
        'class_forces_kcal_mol_A':{k:classes[k]['force_kcal_mol_A'][target] for k in CLASSES},
        'nonbonded_group_forces_kcal_mol_A':{k:baseline[k]['force_kcal_mol_A'][target] for k in GROUPS[1:]},
        'local_terms':len(local_values),'not_a_residue_energy_decomposition':True,'physical_accuracy_established':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);started=time.monotonic();cpu=time.process_time()
    result={'protocol':PROTOCOL,'status':'execution_failure','implementation':pin(__file__),'config':pin(args.config),
        'slurm_job_id':os.environ['SLURM_JOB_ID'],'declared_states':SETTINGS['states'],'declared_local_terms':208,
        'new_DFT_calls':0,'new_MACE_calls':0,'new_GFN_calls':0,'new_ALPB_calls':0,'new_geometries':0,'optimizer_steps':0,
        'rows':[{'source_id':s,'metal':m,'status':'unavailable','input_status':'unavailable','replay_status':'unavailable',
                 'class_status':'unavailable','local_status':'unavailable','contexts_attempted':0,'queries_attempted':0,
                 'queries':[],'checks':[],'diagnostic_systems':{},'local_terms':[]} for s,m in SETTINGS['states']]}
    try:
        cfg=json.loads(args.config.read_text());require(cfg['protocol']==PROTOCOL and cfg['settings']==SETTINGS,'Frozen diagnostic settings differ')
        require(pin(__file__)['sha256']==cfg['implementation']['sha256'] and pin(sys.executable)==cfg['python'],'Source/Python differs')
        verify(cfg['implementation']);verify(cfg['specification'])
        require(cfg['force_audit_result']['sha256']==AUDIT_SHA,'Exact force audit required');audit=read(cfg['force_audit_result'])
        require(audit['status']=='complete' and audit['slurm_job_id']=='13032' and len(audit['rows'])==4
                and [[r['source_id'],r['metal']] for r in audit['rows']]==SETTINGS['states'],'Four original states required')
        for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
            require(os.environ.get(name)=='1','Single-thread native environment required: '+name)
        import numpy as np
        import openmm as mm
        from openmm import unit
        require(mm.__version__==cfg['openmm_version']=='8.5.1','OpenMM version differs')
        result.update(settings=SETTINGS,force_audit_result=cfg['force_audit_result'],
            runtime={'python':pin(sys.executable),'numpy_version':np.__version__,'openmm_version':mm.__version__,'platform':'Reference'})
        for state,row in zip(SETTINGS['states'],result['rows']):
            directory=out/(state[0]+'_'+state[1]);directory.mkdir();tick=time.monotonic()
            try:run_state(cfg,audit,state,directory,row,np,mm,unit)
            except Exception as exc:
                row.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
                if row.get('active_stage'):row[row['active_stage']+'_status']='failed'
            row['elapsed_seconds']=time.monotonic()-tick;row['result']=write(directory/'RESULT.json',row)
        verify(result['config']);result['available_states']=sum(r['status']=='available' for r in result['rows'])
        result['status']='complete' if result['available_states']==4 else 'execution_failure'
    except Exception as exc:result.update(failure=repr(exc),traceback=traceback.format_exc())
    result.update(contexts_attempted=sum(r['contexts_attempted'] for r in result['rows']),
        MM_diagnostic_queries_attempted=sum(r['queries_attempted'] for r in result['rows']),
        elapsed_seconds=time.monotonic()-started,process_CPU_seconds=time.process_time()-cpu,
        peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    record=write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':record}))
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
