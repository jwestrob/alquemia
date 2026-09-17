"""Physical coupled-coordinate preparation for the gated MACE response pilot."""
from __future__ import annotations
import argparse
import copy
import itertools
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz,paired
from affordable_peptide import SourceGraph
from affordable_response import source_key,extract
from ggr_sensitivity import (atom_key,physical_velocity,moved_source,materialize_motion,membership,
                             METHOD,ORCA,executed,endpoint_output_validation)
from mace_hybrid import check_atoms,write_xyz
from mace_global_benchmark import physical_preparations,snapshot
from mace_curvature import dft_source

ROOT=Path(__file__).resolve().parents[1]
POLICY='source_original_H_coupled_metal_peptide_scaffold_response_v1'
CASES={'GGR_extended':('GGR_1GLG','ggr_1glg_nma'),
       'GGR_connected':('GGR_1GLG','ggr_1glg_connected'),
       'ALPHA_1F6S':('ALPHA_1F6S','alacta_1f6s_nma'),
       'ALPHA_6IP9':('ALPHA_6IP9','alacta_6ip9_nma')}
H=np.array([.02,math.pi/180])
POINTS={'center':(0.,0.),'s_minus':(-1.,0.),'s_plus':(1.,0.),
        'theta_minus':(0.,-1.),'theta_plus':(0.,1.),'same_minus':(-1.,-1.),
        'same_plus':(1.,1.),'opposite_a':(-1.,1.),'opposite_b':(1.,-1.)}
GRID={**POINTS,**{'half_'+k:tuple(.5*x for x in v) for k,v in POINTS.items() if k!='center'}}
OLD_NAMES={'center':'center','s_minus':'metal_minus','s_plus':'metal_plus',
           'theta_minus':'peptide_minus','theta_plus':'peptide_plus'}
UNAVAILABLE={'response_status':'response_model_not_validated','relaxation_correction_kcal_mol':None,
             'entropy_correction_kcal_mol':None,'calibrated_class':None}


def physical_model(repair):
    graph=SourceGraph(verify(repair['source_structure']),verify(repair['topology_definition']))
    metal=np.array(repair['selected_site']['xyz_A']);options=[]
    for entry in repair['charge_ledger']:
        if entry['kind']!='backbone_carbonyl':continue
        res=graph.locate(entry['source']);oxygen=graph.key(res.key,'O')
        options.append((float(np.linalg.norm(np.array(tuple(graph.atoms[oxygen].pos))-metal)),oxygen,res))
    if not options:raise InvalidArtifact('no declared backbone donor for physical motion')
    distance,oxygen,res=min(options,key=lambda v:(v[0],v[1]))
    carbon=graph.key(res.key,'C');alpha=graph.key(res.key,'CA');nitrogen=graph.amide_next.get(res.key)
    if nitrogen is None:raise InvalidArtifact('actual amide bond missing')
    nextalpha=graph.key(nitrogen[:2],'CA');hydrogens=[h for h,p in graph.hparents.items() if p==nitrogen]
    if len(hydrogens)!=1:raise InvalidArtifact('supported amide needs one source N hydrogen')
    moving=[carbon,oxygen,nitrogen,hydrogens[0]]
    selected={atom_key(a['source']) for a in repair['atom_graph']['source_to_qm'] if a['kind']=='source'}
    if not set(moving+[alpha,nextalpha]).issubset(selected):
        raise InvalidArtifact('complete physical peptide and both alpha anchors required')
    coords={source_key(graph.meta[k]):np.array(tuple(a.pos)) for k,a in graph.atoms.items()}
    origin=coords[source_key(graph.meta[alpha])]
    axis=coords[source_key(graph.meta[nextalpha])]-origin;axis/=np.linalg.norm(axis)
    direction=coords[source_key(graph.meta[oxygen])]-metal;direction/=np.linalg.norm(direction)
    specs=[{'units':'angstrom','axis':direction.tolist(),'moving_source_keys':['metal']},
           {'units':'radian','axis':axis.tolist(),'origin_A':origin.tolist(),
            'moving_source_keys':[source_key(graph.meta[k]) for k in moving],
            'fixed_anchor_keys':[source_key(graph.meta[k]) for k in (alpha,nextalpha)]}]
    mk=graph.key(graph.locate(repair['selected_site']).key,repair['selected_site']['atom'])
    selection={'rule':'nearest_declared_backbone_carbonyl_O_then_source_identity',
               'oxygen':graph.meta[oxygen],'amide_nitrogen':graph.meta[nitrogen],'distance_A':distance}
    return graph,coords,metal,mk,specs,selection


def physical_id(meta):
    return f'{meta["chain"]}/{meta["resnum"]}/{meta["insertion_code"]}/{meta["atom"]}'


def original_full(prep,graph,coords,metal,mk):
    if prep['preparation_details']['terminal_additions'] or prep['cofactor_charge_e']!=0:
        raise InvalidArtifact('this pilot excludes PQQ/terminal additions')
    physical=copy.deepcopy(prep['physical_atoms'])
    original={a['id']:a['xyz_A'] for a in prep['preparation_details']['original_protein_atoms']}
    original.update({a['id']:a['before_A'] for a in prep['water_H_moves']})
    lookup={physical_id(meta):source_key(meta) for meta in graph.meta.values()}
    for a in physical:
        if a['id'] in original:a['xyz_A']=original[a['id']]
        if a['id']=='metal':
            if not np.allclose(a['xyz_A'],metal,atol=1e-9,rtol=0):raise InvalidArtifact('full/source metal mismatch')
            a['source_key']='metal'
        else:
            key=lookup.get(a['id'])
            if key is None or not np.allclose(a['xyz_A'],coords[key],atol=1e-9,rtol=0):
                raise InvalidArtifact('full atom lacks exact original source match: '+a['id'])
            a['source_key']=key
    if len({a['source_key'] for a in physical})!=len(physical):raise InvalidArtifact('duplicate full source atoms')
    return physical


def displaced(repair,graph,coords,metal,mk,specs,z,physical,original_membership):
    moved={k:v.copy() for k,v in coords.items()};newmetal=metal.copy()
    for spec,amount in zip(specs,H*np.array(z)):
        moved,newmetal=moved_source(moved,newmetal,spec,float(amount))
    maximum=max(float(np.linalg.norm(newmetal-metal)),max(float(np.linalg.norm(moved[source_key(meta)]-coords[source_key(meta)]))
                  for meta in graph.meta.values() if meta['element']!='H'))
    if maximum>.05+1e-12:raise InvalidArtifact('trust displacement exceeds0.05Angstrom')
    if membership(graph,moved,newmetal,mk)!=original_membership:raise InvalidArtifact('physical donor membership changed')
    core,_=materialize_motion(repair,moved,newmetal,{})
    full=[('La' if a['id']=='metal' else a['element'],*(newmetal if a['id']=='metal' else moved[a['source_key']])) for a in physical]
    return core,full,{'maximum_heavy_displacement_A':maximum,'typed_donors':'unchanged',
                      'protonation':'unchanged','water_inventory':'unchanged','cap_anchors':'fixed'}


def prepare(source_root,global_preparation,ggr_dft,agreement,output):
    wall_start=time.monotonic();cpu_start=time.process_time()
    global_m,_=physical_preparations(global_preparation)
    gc,gm=dft_source(ggr_dft,True)
    globals_by_id={r['case_id']:r['preparation'] for r in global_m['rows']}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    names=('mace_mechanics.py','ggr_sensitivity.py','ggr_preparation.py','affordable_peptide.py',
           'affordable_response.py','carve_generic.py','coordination_policy.py','run_orca_task_manifest.py',
           'render_orca_runtime_input.py','result_protocol.py')
    pins=snapshot(out,names);cases={};physical_cases={};reused={};dft_tasks=[]
    qroot=out/'quantum';qroot.mkdir();physical_root=out/'physical_systems';physical_root.mkdir()
    for label,(global_id,folder) in CASES.items():
        source=Path(source_root)/folder/'preparation_manifest.json';repair=read_json(source)
        if repair['protocol_id'] not in ('generic_peptide_alpha_caps_native_r2scan3c_dev_v1','ggr_connected_segment_native_r2scan3c_dev_v1'):
            raise InvalidArtifact('unexpected core preparation protocol')
        for metal in ('La','Ca'):
            verify(repair['outputs'][metal]['xyz']);verify(repair['outputs'][metal]['input'])
        paired(verify(repair['outputs']['La']['xyz']),verify(repair['outputs']['Ca']['xyz']),
               repair['outputs']['La']['charge'],repair['outputs']['Ca']['charge'])
        graph,coords,metalpos,mk,specs,selection=physical_model(repair)
        parent=read_json(verify(globals_by_id[global_id]));physical=original_full(parent,graph,coords,metalpos,mk)
        first_membership=membership(graph,coords,metalpos,mk)
        directory=out/label;directory.mkdir();core_jac=[];full_jac=[]
        for spec in specs:
            vel=physical_velocity(coords,spec);_,jac=materialize_motion(repair,coords,metalpos,vel)
            core_jac.append(jac);full_jac.append([vel.get(a['source_key'],[0.,0.,0.]) for a in physical])
        np.save(directory/'core_jacobians.npy',np.array(core_jac));np.save(directory/'full_jacobians.npy',np.array(full_jac))
        physical_ref=None
        if global_id not in physical_cases:
            pd=physical_root/global_id;pd.mkdir(exist_ok=True)
            write_new(pd/'physical_atoms.json',physical);physical_ref=record(pd/'physical_atoms.json')
            physical_cases[global_id]={'physical_atoms':physical_ref,'original_global_preparation':globals_by_id[global_id],
                                      'specs':specs,'grids':{},'endpoints':parent['endpoints'],
                                      'evidence':parent['evidence'],'assembly':parent['assembly'],'microstate':parent['microstate']}
        else:
            prev=physical_cases[global_id]
            if read_json(verify(prev['physical_atoms']))!=physical or prev['specs']!=specs:
                raise InvalidArtifact('GGR representations change physical protein or coordinate measure')
            physical_ref=prev['physical_atoms']
        grids={}
        for point,z in GRID.items():
            core,full,valid=displaced(repair,graph,coords,metalpos,mk,specs,z,physical,first_membership)
            grid={'normalized_coordinates':list(z),'physical_coordinates':(H*np.array(z)).tolist(),
                  'validation':valid,'endpoints':{}}
            for metal in ('La','Ca'):
                core_rows=[(metal if i==0 else a[0],*a[1:]) for i,a in enumerate(core)]
                cp=directory/f'{point}_{metal}.xyz';old_task=None
                if label.startswith('GGR_') and point in OLD_NAMES:
                    key=f'{label.removeprefix("GGR_")}_{OLD_NAMES[point]}_{metal}'
                    old_task=next(t for t in gm['tasks'] if t['task_id']==key)
                    old_rows=xyz(verify(old_task['xyz']))
                    if [a[0] for a in old_rows]!=[a[0] for a in core_rows] or not np.allclose([a[1:] for a in old_rows],[a[1:] for a in core_rows],atol=1e-9,rtol=0):
                        raise InvalidArtifact('regenerated GGR motion differs from exact archive')
                    shutil.copyfile(verify(old_task['xyz']),cp)
                    reused[f'{label}_{point}_{metal}']={'source_task_id':key,'DFT_collection':record(ggr_dft),'xyz':record(cp)}
                elif point=='center':
                    old_rows=xyz(verify(repair['outputs'][metal]['xyz']))
                    if not np.allclose([a[1:] for a in old_rows],[a[1:] for a in core_rows],atol=1e-9,rtol=0):
                        raise InvalidArtifact('regenerated center differs from archive')
                    shutil.copyfile(verify(repair['outputs'][metal]['xyz']),cp)
                else:write_xyz(cp,core_rows)
                charge=repair['outputs'][metal]['charge']
                grid['endpoints'][metal]={'xyz':record(cp),'charge':charge,'spin_multiplicity':1,
                                          'state':check_atoms(xyz(cp),charge)}
                wanted_dft=point in ('same_minus','same_plus') or (not label.startswith('GGR_') and point in OLD_NAMES)
                if wanted_dft:
                    tid=f'{label}_{point}_{metal}';td=qroot/tid;td.mkdir();xp=td/'core.xyz';shutil.copyfile(cp,xp)
                    ip=td/'endpoint.inp';is_grad=point=='center'
                    ip.write_text(f'! {METHOD}{" EnGrad" if is_grad else ""}\n* xyzfile {charge} 1 core.xyz\n')
                    settings={'protocol_id':POLICY,'source_preparation':record(source),'xyz':record(xp),
                              'input':record(ip),'point':point,'physical_coordinates':grid['physical_coordinates'],
                              'charge':charge,'spin_multiplicity':1,'metal':metal,'case_id':label}
                    dft_tasks.append({'task_id':tid,'case_id':label,'point':point,'metal':metal,'charge':charge,'multiplicity':1,
                                      'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),
                                      'engrad_path':str(td/'endpoint.engrad') if is_grad else None,
                                      'task_type':'analytic_gradient' if is_grad else 'single_point','scientific_settings':settings,
                                      'cache_key':cache_key(settings)})
                if point not in physical_cases[global_id]['grids']:
                    physical_cases[global_id]['grids'][point]={'normalized_coordinates':list(z),'endpoints':{}}
                full_end=physical_cases[global_id]['grids'][point]['endpoints']
                full_rows=[(metal if a[0]=='La' else a[0],*a[1:]) for a in full]
                if metal not in full_end:
                    fp=physical_root/global_id/f'{point}_{metal}.xyz';write_xyz(fp,full_rows)
                    full_end[metal]={'xyz':record(fp),'charge':parent['endpoints'][metal]['charge'],'spin_multiplicity':1,
                                    'state':check_atoms(full_rows,parent['endpoints'][metal]['charge'])}
                elif xyz(verify(full_end[metal]['xyz']))!=full_rows:
                    raise InvalidArtifact('same physical protein was moved differently across representations')
            grids[point]=grid
        case={'case_id':label,'global_id':global_id,'source_preparation':record(source),'selected_coordinate':selection,
              'specs':specs,'physical_atoms':physical_ref,'core_jacobians':record(directory/'core_jacobians.npy'),
              'full_jacobians':record(directory/'full_jacobians.npy'),'grids':grids,
              'explicit_waters':repair['explicit_water_inventory'],'source_graph':repair['atom_graph'],
              'evidence':parent['evidence'],'policy_id':POLICY}
        write_new(directory/'mapping.json',case);cases[label]=record(directory/'mapping.json')
    qmanifest={'schema_version':'alquemia.mace_mechanics_DFT.v1','protocol_id':POLICY,'agreement':record(agreement),
               'orca':record(ORCA),'tasks':dft_tasks,'cases':cases,'implementation':pins,
               'execution_policy':{'task_runner':record(ROOT/'scripts/run_orca_task_manifest.py'),
                                   'runtime_renderer':record(ROOT/'scripts/render_orca_runtime_input.py')},**UNAVAILABLE}
    if len(dft_tasks)!=36 or len(reused)!=20:raise InvalidArtifact('unexpected DFT/reuse inventory')
    write_new(qroot/'manifest.json',qmanifest)
    manifest={'schema_version':'alquemia.mace_mechanics_preparation.v1','policy_id':POLICY,
              'agreement':record(agreement),'implementation':pins,'cases':cases,'physical_cases':physical_cases,
              'DFT_tasks':record(qroot/'manifest.json'),'reused_GGR':reused,'DFT_reference':record(ggr_dft),
              'run_inventory':{'new_DFT_initial':36,'conditional_DFT_optimum_validation_maximum':8,
                               'new_MACE_core':116,'new_GB_core':116,'new_short_component':106},
              'preparation_wall_seconds':time.monotonic()-wall_start,'preparation_CPU_seconds':time.process_time()-cpu_start,**UNAVAILABLE}
    write_new(out/'manifest.json',manifest)
    return {'status':'prepared','manifest':record(out/'manifest.json'),'new_DFT_tasks':36,'reused_DFT':20}


def collect_dft(manifest,output):
    m,rows=executed(manifest);gradients={}
    for t in m['tasks']:
        r=rows[t['task_id']]
        if r['status']!='complete':continue
        try:
            is_gradient=t['task_type']=='analytic_gradient'
            r['native_validation']=endpoint_output_validation(t,verify(r['output']),require_gradient=is_gradient)
            if is_gradient:
                case=read_json(verify(m['cases'][t['case_id']]))
                gradients[t['task_id']]=extract(t['engrad_path'],verify(r['output']),verify(t['input']),verify(t['xyz']),verify(case['source_preparation']))
        except (InvalidArtifact,ValueError,OSError) as exc:
            rows[t['task_id']]={**r,'status':'invalid','reason':str(exc)}
    result={'status':'complete' if all(r['status']=='complete' for r in rows.values()) else 'incomplete',
            'manifest':record(manifest),'rows':rows,'gradients':gradients,'implementation':record(__file__),**UNAVAILABLE}
    write_new(output,result);return {'status':result['status'],'completed':sum(r['status']=='complete' for r in rows.values()),'gradients':len(gradients)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for k in ('source-root','global-preparation','ggr-dft','agreement','output'):a.add_argument('--'+k,required=True)
    a=sub.add_parser('collect-dft')
    for k in ('manifest','output'):a.add_argument('--'+k,required=True)
    values=vars(p.parse_args());command=values.pop('command')
    print(json.dumps({'prepare':prepare,'collect-dft':collect_dft}[command](**values),indent=2))
