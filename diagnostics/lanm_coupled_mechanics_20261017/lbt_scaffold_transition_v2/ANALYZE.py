#!/usr/bin/env python3
"""Fixed four-cell scaffold transition: saved references and finite native MM attribution."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import sys
import time
import traceback
import xml.etree.ElementTree as ET

PROTOCOL='lbt_fixed_scaffold_transition_MM_attribution_v2'
IDS=('before_La','after_La','before_Dy','after_Dy')
GROUPS=('retained_bonded','MM_LJ','MM_Coulomb','QM_MM_LJ')
CLASSES=('HarmonicBondForce','HarmonicAngleForce','PeriodicTorsionForce','CMAPTorsionForce')
COMPONENTS=('MACE_vacuum','GFN2_field_difference','MM_ledger','fixed_FF_charge_ALPB')
TERMINAL=list(range(22))+list(range(244,255))
SETTINGS={'physical_atoms':256,'source':'7CCO','region':'baseline','states':list(IDS),
 'terminal_indices':TERMINAL,'named_indices':{'Phe1_CA':3,'Ala17_C':246},'platform':'Reference',
 'retained_term_counts':dict(zip(CLASSES,(87,159,247,3))),
 'terminal_term_counts':dict(zip(CLASSES,(34,61,90,2))), 'CMAP_maps':16,
 'batch_terms':31,'contexts_per_state':9,'queries_per_state':204,
 'maximum_contexts':36,'maximum_queries':816,
 'energy_absolute_kcal_mol':1e-8,'force_component_absolute_kcal_mol_A':1e-5}
EV,HA,BOHR=23.06054783061903,627.509474,.529177210903
FROZEN={'reference_collector':'06e40d37a44c16446a062bbc479b6a49db70c72e6d16057a2fd7f7278916ee0f',
 'unit_constant_source':'4b683f86d4c4f81c6246d5cecd0acfbc1d2565db157321e3d8fbc8073502fde5',
 'whole_model_implementation':'f35e2225473b1df10288c532a930bdd49186eaaf2e3bab40c65026c44024d8bc',
 'reference_result':'9229f7bebb533dd44a024bb80eca10eded2b3a64cc3972b245d65d82b5ecda97',
 'root_scope':'e7493639861f02a5309dda574f7780147d81f2fb8ff209b780eeb06e79ba6a40',
 'model_result':'64f27e852d2403600ef7a73616c7e52cf50cfbace75689a5b932fde926a06ef0',
 'audit_helper':'a7154cc7b62a506eb92a01e29e063927c952c7cbdeecd158cb1d15a6a0098b08',
 'partition_helper':'0d02c79ff22086aa40f366d11951f3f98fcd58b786f4773c6fb7f48d3e340da7'}


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(pin(record['path'])==record,'Changed artifact: '+str(record));return Path(record['path'])


def read(record):return json.loads(verify(record).read_text())


def write(path,obj):
    with Path(path).open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
    return pin(path)


def module(name,record):
    spec=importlib.util.spec_from_file_location(name,verify(record));obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj);return obj


def vectors(a,np):
    x=np.asarray(a,dtype=float);require(x.shape==(256,3) and np.isfinite(x).all(),'Invalid full256 vectors');return x


def summed(values,np):
    return {'energy_kcal_mol':math.fsum(v['energy_kcal_mol'] for v in values),
            'force_kcal_mol_A':sum((np.asarray(v['force_kcal_mol_A']) for v in values),np.zeros((256,3))).tolist()}


def describe(f,ids,directions,np):
    f=vectors(f,np);norm=np.linalg.norm(f,axis=1);i=int(np.argmax(norm))
    return {'maximum_atom_force_kcal_mol_A':float(norm[i]),'maximum_atom_index':i,'maximum_atom_id':ids[i],
      'atom_vector_RMS_kcal_mol_A':float(np.sqrt(np.mean(norm**2))),
      'projections_kcal_mol_A':{k:float(np.sum(f*d)) for k,d in directions.items()},
      'named_atom_forces':{k:f[i].tolist() for k,i in SETTINGS['named_indices'].items()},
      'named_atom_projection_contributions':{k:{n:float(np.sum(f[i]*d[i])) for n,d in directions.items()}
                                             for k,i in SETTINGS['named_indices'].items()}}


def saved_inputs(cfg,np):
    ref,model=read(cfg['reference_result']),read(cfg['model_result'])
    require(ref['implementation']==cfg['reference_collector'] and model['implementation']==cfg['whole_model_implementation']
            and cfg['unit_constants']=={'eV_to_kcal_mol':EV,'Hartree_to_kcal_mol':HA,'Bohr_to_A':BOHR},
            'Original collector/conversion constants differ')
    require(ref['status']==model['status']=='complete' and ref['full_four_available']
            and len(model['rows'])==4 and set(ref['rows'])==set(IDS),'Complete four-cell parents required')
    require(ref['model_completion_result']['sha256']==cfg['model_result']['sha256'],'Reference/model parent differs')
    require(ref['directions']==model['directions']==cfg['directions'],'Frozen directions differ')
    d=read(cfg['directions']);directions={k:vectors(v,np) for k,v in d['directions'].items()}
    require(d['frozen_before_new_reference_calls'] and set(directions)=={'full_deformation_chord','initial_hybrid_Phe1_CA_force'}
            and all(abs(float(np.linalg.norm(v))-1)<=1e-12 for v in directions.values()),'Direction identity/normalization differs')
    cells={}
    for tid in IDS:
        r=next(r for r in model['rows'] if r['cell_id']==tid);s=next(s for s in ref['original_signatures'] if s['task_id']==tid)
        require(r['status']=='available' and r['source_id']=='7CCO' and r['metal']==tid.split('_')[1]
                and r['charge']==s['charge']==-2 and r['physical_multiplicity']==s['physical_multiplicity']==(1 if r['metal']=='La' else 6)
                and s['target_index']==r['target_index']==255 and s['physical_parent']==r['physical_source']
                and s['positions']==r['positions'] and s['physical_ids']==r['physical_ids']==d['physical_ids'], 'Frozen cell identity differs')
        require(cfg['cells'][tid]=={'physical_source':r['physical_source'],'positions':r['positions'],
                'hybrid_point':r['hybrid_point'],'reference_forces':ref['force_reports'][tid],'ALPB':r['ALPB']},'Cell pins differ')
        atoms=read(r['physical_source']);ids=[a['id'] for a in atoms];positions=read(r['positions']);x=vectors(positions['positions_A'],np)
        require(ids==r['physical_ids'] and len(set(ids))==256 and np.array_equal(x,[a['xyz_A'] for a in atoms])
                and ids[3]=='A/1//CA' and ids[246]=='A/17//C' and ids[255]=='metal'
                and [i for i,a in enumerate(atoms) if a['id'].split('/')[:2] in (['A','1'],['A','17'])]==TERMINAL,
                'Physical order/terminal inventory differs')
        report=read(ref['force_reports'][tid]);table=report['force_table']
        require(report['quantity']=='force=-gradient' and report['units']=='kcal/mol/Angstrom'
                and report['fixed_physical_input']==r['positions']
                and [a['physical_id'] for a in table]==ids and [a['index'] for a in table]==list(range(256)), 'Reference force table order differs')
        models={k:vectors([a['force_kcal_mol_A'][k] for a in table],np) for k in table[0]['force_kcal_mol_A']}
        expected={'reference_PBE0_vacuum','reference_PBE0_plus_fixed_ALPB'}
        if tid!='before_Dy':expected.add('whole_MACE_composite')
        if tid!='after_Dy':expected.add('hybrid_composite')
        require(set(models)==expected,'Missing comparator policy differs')
        sr=ref['solvent_admission'][tid];sm=model['solvent_admission'][tid]
        raw_solvent=read(r['ALPB'])
        require(sr['status']=='available' and sr['receipt']==r['ALPB']
                and sr['energy_hartree']==raw_solvent['total_hartree']==r['native_energies']['fixed_ALPB_hartree']
                and sr['components_hartree']==raw_solvent['components_hartree'],'Fixed ALPB raw receipt differs')
        if tid!='after_Dy':
            require(sr['receipt']==sm['receipt'] and sr['energy_hartree']==sm['energy_hartree']
                    and sr['gradient_kcal_mol_A']==sm['gradient_kcal_mol_A']
                    and sr['gradient_hartree_per_bohr']==sm['gradient_hartree_per_bohr'],
                    'Fixed ALPB differs between saved models')
        # Model collection's preparation admission predates its new after-Dy call.
        # The actual completed row/raw receipt and subsequent full-reference admission are authoritative there.
        require(ref['rows'][tid]['status']=='available' and ref['rows'][tid]['electronic_state']==s['electronic_state'],
                'Admitted native reference state differs')
        g=vectors(ref['rows'][tid]['gradient_kcal_mol_per_A'],np);gs=vectors(sr['gradient_kcal_mol_A'],np)
        require(np.array_equal(models['reference_PBE0_vacuum'],-g)
                and np.array_equal(models['reference_PBE0_plus_fixed_ALPB'],-(g+gs)), 'Admitted reference force replay differs')
        task=ref['rows'][tid]['task']
        require(all(task[k]==s[k] for k in ('task_id','source_id','geometry','metal','charge','multiplicity',
              'physical_multiplicity','electronic_state','physical_parent','positions','physical_ids','atom_count','target_index','ALPB'))
              and task['xyz']['sha256']==s['xyz_sha256'] and task['input']['sha256']==s['input_sha256']
              and task['task_type']=='analytic_gradient', 'Reference inline task/signature differs')
        for name in ('xyz','input','physical_parent','positions','ALPB'):verify(task[name])
        for name in ('receipt','engrad','output'):verify(ref['rows'][tid][name])
        hybrid=None
        if r['hybrid_point'] is not None:
            hybrid=read(r['hybrid_point']);p=hybrid['point'];mapping=read(hybrid['mapping'])
            require(hybrid['status']=='complete' and hybrid['state']['source']=='7CCO'
                    and hybrid['state']['region']=='baseline' and hybrid['state']['metal']==r['metal']
                    and hybrid['state']['whole_charge']==-2 and hybrid['state']['MACE_multiplicity']==r['physical_multiplicity']
                    and hybrid['positions']==r['positions'] and mapping['physical_atoms']==atoms
                    and np.array_equal(p['positions_A'],x),'Hybrid physical/state replay differs')
            components={k:vectors(p['component_gradients_kcal_mol_A'][k],np) for k in COMPONENTS}
            require(set(p['component_gradients_kcal_mol_A'])==set(COMPONENTS)
                    and np.array_equal(sum(components[k] for k in COMPONENTS),p['gradient_kcal_mol_A'])
                    and np.array_equal(-vectors(p['gradient_kcal_mol_A'],np),models['hybrid_composite'])
                    and np.array_equal(components['fixed_FF_charge_ALPB'],gs),'Archived hybrid component/ALPB closure differs')
            require(math.fsum(p['component_energies_kcal_mol'].values())==p['total_energy_kcal_mol'], 'Hybrid energy closure differs')
        else:require(tid=='after_Dy','Unexpected absent hybrid')
        if r['force_arrays'] is not None:
            arrays=read(r['force_arrays'])
            require(arrays['physical_ids']==ids and arrays['units']=='kcal/mol/Angstrom'
                    and arrays['quantity']=='force = negative derivative','Whole force provenance differs')
            # Exact original CONTINUATION.force_report / MODEL_COMPLETION.package arithmetic.
            # REFERENCE.admit_alpb stores gs=(g*HA)/BOHR; whole force stores sf=-g*(HA/BOHR).
            # Each representation retains its own exact check against the same admitted raw vector.
            raw=read(r['MACE']);mf=vectors(raw['forces_eV_A'],np)*EV
            gh=vectors(sr['gradient_hartree_per_bohr'],np);sf=-gh*(HA/BOHR)
            require(np.array_equal(gs,gh*HA/BOHR),'Reference ALPB conversion replay differs')
            require(np.array_equal(arrays['components']['MACE_vacuum'],mf)
                    and np.array_equal(arrays['components']['fixed_FF_ALPB'],sf)
                    and np.array_equal(arrays['components']['composite'],mf+sf)
                    and np.array_equal(arrays['components']['composite'],models['whole_MACE_composite']),
                    'Whole raw force/collector replay differs')
        else:require(tid=='before_Dy' and r['force_status']=='unavailable_energy_only_archived_request','Unexpected missing whole force')
        cells[tid]={'model_row':r,'signature':s,'atoms':atoms,'ids':ids,'x':x,'models':models,'hybrid':hybrid,
                    'reference_report':ref['force_reports'][tid],'solvent':sr}
    for prefix in ('before','after'):
        a,b=(cells[prefix+'_'+m] for m in ('La','Dy'))
        require(a['atoms']==b['atoms'] and np.array_equal(a['x'],b['x']), 'Shared metal-pair geometry differs')
    return ref,model,cells,directions


def mm_inputs(cfg,cell,row,h,np,mm):
    metal=cell['model_row']['metal'];p=cell['hybrid'];authority=p or read(cfg['cells']['before_Dy']['hybrid_point'])
    ledger_pin=authority['source_task']['ledger'];ledger=read(ledger_pin);part=read(ledger['partition'])
    require(part['source_id']=='7CCO' and part['region_name']=='baseline' and ledger_pin==cfg['ledger'],'Ledger identity differs')
    artifacts={k:read(v) for k,v in ledger['artifacts'].items()};particles=artifacts['particles.json'];ids=cell['ids']
    require([a['id'] for a in particles]==ids and all(particles[i]['region']=='MM' for i in TERMINAL),'Terminal MM identity differs')
    verify(ledger['forcefield']);verify(ledger['ion_parameters'])
    system_pin=cfg['systems'][metal];system=mm.XmlSerializer.deserialize(verify(system_pin).read_text())
    parent=mm.XmlSerializer.deserialize(verify(ledger['native_system']).read_text());archive=None
    if p is not None:
        receipts=p['component_receipts'];archive=read(receipts['MM']);mc=read(archive['config']);execution=read(receipts['MM_execution'])
        require(archive['status']=='complete' and archive['stage']=='mm' and archive['group_order']==list(GROUPS)
          and archive['system']==mc['system']==system_pin and archive['positions']==mc['positions']==cell['model_row']['positions']
          and archive['openmm_version']==mc['version']==cfg['openmm_version'] and archive['platform']=='Reference'
          and archive['implementation']==mc['implementation'] and mc['implementation']['sha256']==h.MM_SHA
          and mc['builder']['sha256']==h.BUILDER_SHA and mc['python']==cfg['python'],'Archived MM recipe differs')
        expected=[str(verify(mc['python'])),str(verify(mc['implementation'])),'--stage','mm','--config',str(verify(archive['config'])),
                  '--output',str(verify(receipts['MM']).parent)]
        verify(mc['builder']);require(execution['returncode']==0 and execution['owned_process_group_empty']
                 and execution['command']==expected,'Archived MM execution differs')
        for k in ('stdout','stderr'):verify(execution[k])
        require(np.array_equal(archive['gradient_kcal_mol_A'],p['point']['component_gradients_kcal_mol_A']['MM_ledger'])
                and archive['energy_kcal_mol']==p['point']['component_energies_kcal_mol']['MM_ledger'],'MM POINT component differs')
        task=read(p['task']);require(task['ledger']==ledger_pin and task['partition']==ledger['partition'],'Hybrid task/ledger differs')
    forces=list(system.getForces());require(system.getNumParticles()==256 and system.getNumConstraints()==0 and len(forces)==8,'MM system inventory differs')
    expect=[('MM_LJ','NonbondedForce',1),('MM_Coulomb','NonbondedForce',2),
            ('QM_MM_LJ_generic','CustomNonbondedForce',3),('QM_MM_LJ_native_exceptions','CustomBondForce',3)]
    require([(f.getName(),type(f).__name__,f.getForceGroup()) for f in forces if type(f).__name__ not in CLASSES]==expect,'Nonbonded ledger differs')
    decisions=read(ledger['bonded_decisions']);selected=artifacts['bonded_terms.json'];byclass={};terms=[]
    for kind in CLASSES:
        fs=[f for f in forces if type(f).__name__==kind];ns=[f for f in parent.getForces() if type(f).__name__==kind]
        require(len(fs)==len(ns)==1 and fs[0].getForceGroup()==0 and fs[0].getName()=='retained_'+kind,'Bonded class differs')
        f=fs[0];byclass[kind]=f;indices=sorted(t['index'] for t in selected if t['force']==kind)
        raw=h.term_attributes(mm,ns[0]);actual=h.term_attributes(mm,f)
        require(len(indices)==len(set(indices))==SETTINGS['retained_term_counts'][kind]
                and actual==[raw[i] for i in indices],'Native/retained parameter identity differs')
        if kind=='CMAPTorsionForce':
            require(f.getNumMaps()==SETTINGS['CMAP_maps'] and ET.tostring(h.force_xml(mm,f).find('Maps'))
                    ==ET.tostring(h.force_xml(mm,ns[0]).find('Maps')),'CMAP maps differ')
        local=[]
        for i,native_i in enumerate(indices):
            atoms=h.term_indices(f,i)
            if not set(atoms).intersection(TERMINAL):continue
            d=next(d for d in decisions if d['force']==kind and d['index']==native_i)
            require(d['keep'] and d['atoms']==atoms and d['source_ids']==[ids[j] for j in atoms],'Term ledger differs')
            local.append({'force_class':kind,'retained_index':i,'native_index':native_i,'atoms':atoms,
              'atom_ids':[ids[j] for j in atoms],'reason':d['reason'],'native_parameter_XML':actual[i],
              'terminal_atoms':sorted(set(atoms).intersection(TERMINAL))})
        require(len(local)==SETTINGS['terminal_term_counts'][kind],'Terminal term denominator differs');terms.extend(local)
    require(len(terms)==187,'187 unique terminal-touching interactions required')
    row['provenance']={'cell':cfg['cells'][row['cell_id']],'system':system_pin,'ledger':ledger_pin,
      'native_system':ledger['native_system'],'forcefield':ledger['forcefield'],'ion_parameters':ledger['ion_parameters'],
      'artifacts':ledger['artifacts'],'bonded_decisions':ledger['bonded_decisions'],
      'MM_result':None if p is None else p['component_receipts']['MM'],
      'MM_execution':None if p is None else p['component_receipts']['MM_execution'],
      'new_after_Dy_component_only':p is None,'native_to_retained_parameters_exact':True}
    row['physical_classification']=[{'index':i,'physical_id':a['id'],'region':a['region'],
       'terminal_residue':i in TERMINAL,'boundary_roles':[
          role for cap in part['caps'] for role,key in (('retained','retained_source_id'),('omitted','omitted_source_id'))
          if cap[key]==a['id']]} for i,a in enumerate(particles)]
    row['input_status']='passed';return system,archive,byclass,terms,particles


def subset(force,indices,name,group,h,mm):
    if type(force).__name__!='CMAPTorsionForce':return h.subset_force(force,indices,name,group,mm)
    new=mm.CMAPTorsionForce();new.setName(name);new.setForceGroup(group)
    new.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
    for i in range(force.getNumMaps()):new.addMap(*force.getMapParameters(i))
    for i in indices:new.addTorsion(*force.getTorsionParameters(i))
    require(h.term_attributes(mm,new)==[h.term_attributes(mm,force)[i] for i in indices]
      and ET.tostring(h.force_xml(mm,new).find('Maps'))==ET.tostring(h.force_xml(mm,force).find('Maps')),'CMAP split differs')
    return new


def split_system(original,batch,h,mm):
    system=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(original))
    for i in reversed(range(system.getNumForces())):system.removeForce(i)
    for f in original.getForces():
        kind=type(f).__name__;chosen=[(g,t) for g,t in enumerate(batch) if t['force_class']==kind]
        if not chosen:
            copy=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(f));copy.setForceGroup(31);system.addForce(copy);continue
        n=getattr(f,h.TERM_API[kind][1])();indices=[t['retained_index'] for _,t in chosen];rest=[i for i in range(n) if i not in indices]
        require(len(indices)==len(set(indices)) and sorted(indices+rest)==list(range(n)),'Split accounting differs')
        system.addForce(subset(f,rest,f.getName()+'_remainder',31,h,mm))
        for group,t in chosen:system.addForce(subset(f,[t['retained_index']],f.getName()+'_term_'+str(t['retained_index']),group,h,mm))
    require(system.getNumParticles()==256 and system.getNumConstraints()==0,'Split physical system differs');return system


def query(system,x,queries,label,out,row,np,mm,unit):
    require(row['contexts_attempted']<9,'Context ceiling');row['contexts_attempted']+=1
    directory=out/label;directory.mkdir();xml=directory/'SYSTEM.xml';xml.write_text(mm.XmlSerializer.serialize(system))
    row['diagnostic_systems'][label]=pin(xml);integrator=mm.VerletIntegrator(.001);context=None;values={}
    try:
        context=mm.Context(system,integrator,mm.Platform.getPlatformByName('Reference'));context.setPositions(x*.1*unit.nanometer)
        for name,mask in queries:
            require(row['queries_attempted']<204,'Query ceiling');row['queries_attempted']+=1
            st=context.getState(getEnergy=True,getForces=True,groups=mask)
            energy=float(st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
            f=vectors(st.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole/unit.angstrom),np)
            require(math.isfinite(energy),'Nonfinite MM energy');v={'energy_kcal_mol':energy,'force_kcal_mol_A':f.tolist()}
            rec=write(directory/(name+'.json'),dict(v,system=pin(xml),groups=mask,positions=row['provenance']['cell']['positions']))
            row['queries'].append({'label':label,'name':name,'groups':mask,'result':rec});values[name]=v
        return values
    finally:
        if context is not None:del context
        del integrator


def run_cell(cell,prepared,row,out,h,np,mm,unit,directions,common_after):
    system,archive,byclass,terms,_=prepared;x=cell['x'];row['active_stage']='replay'
    require(archive is not None or common_after is not None,'After-Dy requires admitted after-La control before queries')
    base=query(system,x,[(n,1<<i) for i,n in enumerate(GROUPS)]+[('total',-1)],'original',out,row,np,mm,unit)
    h.close('original_group_sum',summed([base[n] for n in GROUPS],np),base['total'],row,np)
    if archive is not None:
        for n in GROUPS:
            h.close('archived_'+n,base[n],{'energy_kcal_mol':archive['component_energies_kcal_mol'][n],
               'force_kcal_mol_A':(-np.asarray(archive['component_gradients_kcal_mol_A'][n])).tolist()},row,np)
        h.close('archived_total',base['total'],{'energy_kcal_mol':archive['energy_kcal_mol'],
                 'force_kcal_mol_A':(-np.asarray(archive['gradient_kcal_mol_A'])).tolist()},row,np)
        row['replay_status']='passed'
    else:
        require(common_after is not None,'After-Dy requires actual admitted after-La group control')
        for n in GROUPS[:3]:h.close('same_geometry_common_'+n,base[n],common_after[n],row,np)
        row['replay_status']='new_component_common_groups_passed_no_archived_after_Dy';row['full_hybrid_available']=False
    row['original_groups']=base;row['active_stage']='class'
    classified=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))
    for f in classified.getForces():
        if type(f).__name__ in CLASSES:f.setForceGroup(4+CLASSES.index(type(f).__name__))
    classes=query(classified,x,[(n,1<<(4+i)) for i,n in enumerate(CLASSES)]+[('total',-1)],'classes',out,row,np,mm,unit)
    h.close('bonded_class_sum',summed([classes[n] for n in CLASSES],np),base['retained_bonded'],row,np)
    h.close('class_system_total',classes['total'],base['total'],row,np);row['class_status']='passed';row['classes']=classes
    row['active_stage']='terms';terminal_sums={k:np.zeros((len(TERMINAL),3)) for k in CLASSES}
    for batch_i,start in enumerate(range(0,len(terms),31)):
        batch=terms[start:start+31];split=split_system(system,batch,h,mm)
        values=query(split,x,[('term_'+str(start+i),1<<i) for i in range(len(batch))]+[('total',-1)],
                     'batch_'+str(batch_i),out,row,np,mm,unit)
        h.close('batch_'+str(batch_i)+'_full_system',values['total'],base['total'],row,np)
        for i,t in enumerate(batch):
            value=values['term_'+str(start+i)];f=np.asarray(value['force_kcal_mol_A']);terminal_sums[t['force_class']]+=f[TERMINAL]
            desc=h.term_description(t,byclass[t['force_class']],x,np,unit) if t['force_class']!='CMAPTorsionForce' else {
              'map_index':int(byclass[t['force_class']].getTorsionParameters(t['retained_index'])[0]),
              'parameter_units':{'map_energies':'kJ/mol','torsion_atoms':'physical_indices'},'all_16_maps_in_query_system':True}
            row['terms'].append(dict(t,**desc,energy_kcal_mol=value['energy_kcal_mol'],
                 terminal_forces_kcal_mol_A=f[TERMINAL].tolist(),force_summary=describe(f,cell['ids'],directions,np),
                 query=next(q['result'] for q in row['queries'] if q['label']=='batch_'+str(batch_i) and q['name']=='term_'+str(start+i))))
    require(len(row['terms'])==187,'Terminal term denominator incomplete')
    for n in CLASSES:h.close('all_terminal_atoms_term_sum_'+n,{'force_kcal_mol_A':terminal_sums[n].tolist()},
                    {'force_kcal_mol_A':np.asarray(classes[n]['force_kcal_mol_A'])[TERMINAL].tolist()},row,np,energy=False)
    row['terms_status']='passed';row['active_stage']=None;row['status']='available'
    row['force_summaries']={n:describe(v['force_kcal_mol_A'],cell['ids'],directions,np)
                            for n,v in dict(base,**{k:classes[k] for k in CLASSES}).items()}


def comparison(ref,model,cells,rows,directions,out,np):
    ids=cells['before_La']['ids'];result={'work_authority':{'reference':ref['comparison'],'whole_model':model['comparison']},
      'per_metal':{},'force_error_definition':'model minus PBE0 plus identical fixed-FF ALPB',
      'force_change_definition':'after minus before at the two fixed geometries',
      'not_a_Hessian_or_unique_stiffness_fit':True,'no_quantum_interaction_energy_attribution':True}
    allvectors={}
    for metal in ('La','Dy'):
        a,b=(cells[p+'_'+metal] for p in ('before','after'));ra,rb=(next(r for r in rows if r['cell_id']==p+'_'+metal) for p in ('before','after'))
        entry={'force_models':{},'hybrid_components':{},'MM_groups':{},'MM_classes':{},'terminal_interaction_works':[]}
        for name in ('reference_PBE0_plus_fixed_ALPB','whole_MACE_composite','hybrid_composite'):
            available=name in a['models'] and name in b['models'];v={'status':'available' if available else 'unavailable',
              'before_available':name in a['models'],'after_available':name in b['models']}
            fields={}
            for label,c in (('before',a),('after',b)):
                if name not in c['models']:continue
                f=c['models'][name];fields['force_'+label]=f
                if name!='reference_PBE0_plus_fixed_ALPB':fields['error_'+label]=f-c['models']['reference_PBE0_plus_fixed_ALPB']
            if available:
                fields['force_change']=fields['force_after']-fields['force_before']
                if name!='reference_PBE0_plus_fixed_ALPB':fields['error_change']=fields['error_after']-fields['error_before']
            v['summaries']={k:describe(f,ids,directions,np) for k,f in fields.items()};entry['force_models'][name]=v
            allvectors[metal+'_'+name]={k:f.tolist() for k,f in fields.items()}
        for name in COMPONENTS:
            fields={}
            for label,c in (('before',a),('after',b)):
                if c['hybrid'] is not None:
                    fields[label]=-vectors(c['hybrid']['point']['component_gradients_kcal_mol_A'][name],np)
            if len(fields)==2:fields['change']=fields['after']-fields['before']
            entry['hybrid_components'][name]={'before_available':'before' in fields,'after_available':'after' in fields,
                'summaries':{k:describe(v,ids,directions,np) for k,v in fields.items()},
                'component_only_not_quantum_bond_attribution':True}
            allvectors[metal+'_hybrid_component_'+name]={k:v.tolist() for k,v in fields.items()}
        for key,source,names in (('MM_groups','original_groups',GROUPS+('total',)),('MM_classes','classes',CLASSES)):
            for n in names:
                va,vb=ra[source][n],rb[source][n];df=np.asarray(vb['force_kcal_mol_A'])-np.asarray(va['force_kcal_mol_A'])
                entry[key][n]={'work_kcal_mol':vb['energy_kcal_mol']-va['energy_kcal_mol'],
                              'force_change':describe(df,ids,directions,np)}
                allvectors[metal+'_'+key+'_'+n]={'before':va['force_kcal_mol_A'],'after':vb['force_kcal_mol_A'],'change':df.tolist()}
        for ta,tb in zip(ra['terms'],rb['terms']):
            identity=('force_class','retained_index','native_index','atoms','atom_ids','reason','native_parameter_XML','terminal_atoms')
            require(all(ta[k]==tb[k] for k in identity),'Before/after interaction identity differs')
            entry['terminal_interaction_works'].append({**{k:ta[k] for k in identity},
                'before_energy_kcal_mol':ta['energy_kcal_mol'],'after_energy_kcal_mol':tb['energy_kcal_mol'],
                'work_kcal_mol':tb['energy_kcal_mol']-ta['energy_kcal_mol'],
                'whole_interaction_energy_not_residue_energy':True,'before_query':ta['query'],'after_query':tb['query']})
        entry['terminal_touching_work_sums']={}
        for label,indices in (('union',set(TERMINAL)),('Phe1',set(range(22))),('Ala17',set(range(244,255))),
                              ('both_termini_overlap',None)):
            if indices is None:
                selected=[t for t in entry['terminal_interaction_works'] if
                  set(t['atoms']) & set(range(22)) and set(t['atoms']) & set(range(244,255))]
            else:selected=[t for t in entry['terminal_interaction_works'] if set(t['atoms']) & indices]
            entry['terminal_touching_work_sums'][label]={'count':len(selected),'whole_interaction_not_residue_energy':True,
              'by_class':{k:{'count':sum(t['force_class']==k for t in selected),
                'work_kcal_mol':math.fsum(t['work_kcal_mol'] for t in selected if t['force_class']==k)} for k in CLASSES},
              'work_kcal_mol':math.fsum(t['work_kcal_mol'] for t in selected)}
        result['per_metal'][metal]=entry
    result['full256_force_offset_and_change']=write(out/'FORCE_TRANSITIONS.json',{'physical_ids':ids,'units':'kcal/mol/angstrom',
        'frozen_directions':ref['directions'],'arrays':allvectors})
    return result


def report(result,out):
    lines=['# Fixed scaffold transition audit','',
      'Four fixed reference cells; both metals use the same La-before and La-proposed after coordinates.',
      'All energies below are after minus before (kcal/mol); forces are kcal/mol/Angstrom.',
      'MM interaction energies are whole-interaction contributions, not residue energies or quantum bond energies.','']
    if result['status']!='complete':
        lines+=['**Incomplete: no complete transition interpretation is admitted.**','']
        for r in result['rows']:lines.append('- '+r['cell_id']+': '+r['status']+'; '+r.get('failure','not reached'))
    else:
        c=result['comparison'];lines+=['## MM work','', '| Group/class | La | Dy |','|---|---:|---:|']
        for key,names in (('MM_groups',GROUPS+('total',)),('MM_classes',CLASSES)):
            for n in names:lines.append('| '+n+' | '+' | '.join(f"{c['per_metal'][m][key][n]['work_kcal_mol']:.6f}" for m in ('La','Dy'))+' |')
        lines+=['','## Terminal-touching interaction work','',
          'The union counts every interaction once. Per-terminus sets can overlap; overlap is listed explicitly.','',
          '| Set / class | Count | La | Dy |','|---|---:|---:|---:|']
        for label in ('union','Phe1','Ala17','both_termini_overlap'):
            values=[c['per_metal'][m]['terminal_touching_work_sums'][label] for m in ('La','Dy')]
            lines.append(f"| {label} / total | {values[0]['count']} | {values[0]['work_kcal_mol']:.6f} | {values[1]['work_kcal_mol']:.6f} |")
            for n in CLASSES:lines.append(f"| {label} / {n} | {values[0]['by_class'][n]['count']} | {values[0]['by_class'][n]['work_kcal_mol']:.6f} | {values[1]['by_class'][n]['work_kcal_mol']:.6f} |")
        paired=list(zip(c['per_metal']['La']['terminal_interaction_works'],c['per_metal']['Dy']['terminal_interaction_works']))
        paired.sort(key=lambda pair:(-max(abs(t['work_kcal_mol']) for t in pair),pair[0]['force_class'],pair[0]['native_index']))
        lines+=['','## Fifteen largest absolute interaction works','',
          'Descriptive ranking by the larger absolute La/Dy work; all 187 interactions remain in RESULT.json.','',
          '| Class / native index | Atoms | La | Dy |','|---|---|---:|---:|']
        for a,b in paired[:15]:lines.append(f"| {a['force_class']} / {a['native_index']} | {', '.join(a['atom_ids'])} | {a['work_kcal_mol']:.6f} | {b['work_kcal_mol']:.6f} |")
        lines+=['','## Force offsets and changes','',
          'Errors are model minus PBE0 plus the identical fixed-charge ALPB. Change is after minus before.','',
          '| Metal / model / quantity | Full256 vector RMS | Phe1-direction projection | Deformation-chord projection |',
          '|---|---:|---:|---:|']
        for metal in ('La','Dy'):
            for name,v in c['per_metal'][metal]['force_models'].items():
                for key,s in v['summaries'].items():
                    lines.append(f"| {metal} / {name} / {key} | {s['atom_vector_RMS_kcal_mol_A']:.6f} | {s['projections_kcal_mol_A']['initial_hybrid_Phe1_CA_force']:.6f} | {s['projections_kcal_mol_A']['full_deformation_chord']:.6f} |")
                if v['status']!='available':lines.append(f"| {metal} / {name} / transition unavailable | — | — | — |")
        lines+=['','Before-Dy whole-MACE forces and after-Dy hybrid forces/work remain unavailable. After-Dy MM is a component-only diagnostic.',
          'Two endpoints do not measure a Hessian or separate equilibrium offsets from stiffness uniquely. No parameter fit, accuracy, affinity or preference claim follows.',
          'Inherited limitations: conditional terminal chemistry, vacuum dianion, frozen-f representation, SCF evidence, nonstationary endpoints, and no independent numerical qualification of these whole-reference forces.']
    lines+=['',f"Actual MM diagnostics: {result['contexts_attempted']} contexts; {result['MM_diagnostic_queries_attempted']} energy/full-force queries.",'']
    path=out/'REPORT.md';path.write_text('\n'.join(lines));return pin(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm allocation required')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    result={'protocol':PROTOCOL,'status':'execution_failure','config':pin(args.config),'implementation':pin(__file__),
      'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'declared_states':list(IDS),'settings':SETTINGS,
      'new_DFT_calls':0,'new_MACE_calls':0,'new_GFN_calls':0,'new_ALPB_calls':0,'optimizer_steps':0,'new_geometries':0,
      'rows':[{'cell_id':tid,'status':'unavailable','input_status':'unavailable','replay_status':'unavailable',
               'class_status':'unavailable','terms_status':'unavailable','contexts_attempted':0,'queries_attempted':0,
               'checks':[],'queries':[],'diagnostic_systems':{},'terms':[]} for tid in IDS]}
    try:
        cfg=read(result['config'])
        previous=read(cfg['preserved_v1_failure']);verify(cfg['previous_implementation']);verify(cfg['previous_specification'])
        require(previous['status']=='execution_failure' and previous['slurm_job_id']=='14785'
                and previous['implementation']==cfg['previous_implementation']
                and previous['contexts_attempted']==previous['MM_diagnostic_queries_attempted']==0
                and previous['failure']=="ValueError('Whole comparator/ALPB force replay differs')",
                'Preserved zero-query V1 failure differs')
        result['technical_recovery']={'preserved_failure':cfg['preserved_v1_failure'],
            'original_collector':cfg['reference_collector'],'unit_constant_source':cfg['unit_constant_source'],
            'same_raw_ALPB_two_original_arithmetic_representations':True,'tolerances_changed':False}
        require(cfg['protocol']==PROTOCOL and cfg['settings']==SETTINGS,'Frozen scope differs')
        require(pin(__file__)['sha256']==cfg['implementation']['sha256'] and pin(sys.executable)==cfg['python'],'Source/interpreter differs')
        for key in ('implementation','specification','root_scope','plan','directions','ledger'):verify(cfg[key])
        for key,sha in FROZEN.items():require(cfg[key]['sha256']==sha,'Frozen parent/helper differs');verify(cfg[key])
        require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'One CPU required')
        for n in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
            require(os.environ.get(n)=='1','Single-thread environment required: '+n)
        import numpy as np
        import openmm as mm
        from openmm import unit
        require(mm.__version__==cfg['openmm_version']=='8.5.1','OpenMM version differs')
        h=module('scaffold_old_term_audit',cfg['audit_helper']);prep=module('scaffold_partition_definition',cfg['partition_helper']);prep.mm=mm
        result['runtime']={'python':pin(sys.executable),'numpy':np.__version__,'openmm':mm.__version__,'platform':'Reference'}
        ref,model,cells,directions=saved_inputs(cfg,np);result['reference_caveats']=ref['caveats']
        result['reference_native_evidence']={tid:{k:ref['rows'][tid][k] for k in
          ('electronic_state','SCF_evidence','receipt','output','engrad')} for tid in IDS}
        prepared={}
        for row in result['rows']:
            row['active_stage']='input'
            try:prepared[row['cell_id']]=mm_inputs(cfg,cells[row['cell_id']],row,h,np,mm)
            except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc(),input_status='failed')
        require(len(prepared)==4,'All four input bindings required before any MM query')
        common,proof=prep.common_mm_definitions({m:{'system':cfg['systems'][m]} for m in ('La','Dy')},prepared['before_La'][4])
        require(all(common.values()),'Common MM definitions differ');result['common_MM_definition_proof']=proof
        for row in result['rows']:
            tid=row['cell_id'];directory=out/tid;directory.mkdir();tick=time.monotonic()
            try:
                after=next(r for r in result['rows'] if r['cell_id']=='after_La')
                control=after.get('original_groups') if after['status']=='available' else None
                run_cell(cells[tid],prepared[tid],row,directory,h,np,mm,unit,directions,control)
            except Exception as exc:
                row.update(status='unavailable',failure=repr(exc),traceback=traceback.format_exc())
                if row.get('active_stage'):row[row['active_stage']+'_status']='failed'
            row['elapsed_seconds']=time.monotonic()-tick;row['result']=write(directory/'RESULT.json',row)
        require(all(r['status']=='available' for r in result['rows']),'Incomplete four-state audit')
        result['comparison']=comparison(ref,model,cells,result['rows'],directions,out,np)
        verify(result['config']);result['status']='complete'
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result.update(contexts_attempted=sum(r['contexts_attempted'] for r in result['rows']),
       MM_diagnostic_queries_attempted=sum(r['queries_attempted'] for r in result['rows']),
       elapsed_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu,
       peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    try:result['report']=report(result,out)
    except Exception as exc:result.update(status='execution_failure',report_failure=repr(exc),report_traceback=traceback.format_exc())
    record=write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':record}))
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
