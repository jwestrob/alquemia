#!/usr/bin/env python3
"""Prepare the fixed six-endpoint response-transfer continuation in existing runners."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
import shutil

from affordable_common import InvalidArtifact,record,verify,read_json,write_new,xyz,cache_key
from mace_hybrid import check_atoms
from response_probe_archive import paired_coordinates,canonical_mapping,generic_mapping,selected_donors


def prepare(root, agreement, output):
    root=Path(root).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    paths={f'GGR_{n.upper()}':root/f'workspaces/ggr_mechanism_20260915/stage_b_prepared_v1/{n}/alpha_caps/preparation_manifest.json'
           for n in ('2fw0','2fvy')}
    paths['8GY2']=root/'workspaces/benchmark_augmentation_20260918/prepared_8gy2_v1/01_8GY2/pmdh_fc_holdout_8gy2_carve_manifest.json'
    cases={}
    for name,path in paths.items():
        prep=read_json(path);states={};is_pqq=name=='8GY2'
        for metal in ('Ca','La'):
            if is_pqq:
                ep=prep['outputs'];xp=ep[metal+'_xyz'];ip=ep[metal+'_input']
                q=-2 if metal=='Ca' else -1
            else:
                ep=prep['outputs'][metal];xp=ep['xyz'];ip=ep['input'];q=ep['charge']
            verify(xp);verify(ip)
            states[metal]={'xyz':xp,'charge':q,'spin_multiplicity':1,'state':check_atoms(xyz(verify(xp)),q),
                           'baseline_input':ip}
        _,coords=paired_coordinates(states)
        mapper=canonical_mapping if is_pqq else generic_mapping
        mapping,physical,cofactor,source=mapper(prep,coords)
        donors=selected_donors(mapping,physical,coords,prep)
        cases[name]={'preparation':record(path),'states':states,'source_structure':source,'mapping':mapping,
            'donors':donors,'cofactor_indices':cofactor,'source_protocol':prep['protocol_id'],
            'evidence_stratum':'PQQ_enzyme_structural_association' if is_pqq else 'direct_direction_condition_qualified',
            'expected_label':'Ca','biological_group':'O05542' if is_pqq else 'GGR_MglB',
            'DFT_policy':'unchanged_canonical_input' if is_pqq else 'native_TightSCF_same_source_coordinates',
            'whole_protein_supported':False if is_pqq else None}
    pp=out/'preparation.json'
    write_new(pp,{'cases':cases,'agreement':record(agreement),'new_MACE_endpoints':6,'new_DFT_endpoints':6,
                  'baseline_changed':False,'implementation':record(__file__)})
    parent_path=root/'workspaces/mace_canonical_20260916/mace_v2/medium/manifest.json'
    parent=read_json(parent_path);md=out/'mace';md.mkdir();impl=md/'implementation';impl.mkdir()
    for name,pin in parent['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    # Preserve the qualified worker/execute implementation unchanged; adapter
    # adds only the bounded validator and ordinary result collection.
    shutil.copyfile(impl/'mace_hybrid.py',impl/'response_probe_backend.py')
    shutil.copyfile(Path(__file__).with_name('response_probe_run.py'),impl/'mace_hybrid.py')
    pins={p.name:record(p) for p in impl.glob('*.py')}
    model=copy.deepcopy(parent['model']);model['preparation_policy']='explicit_frozen_source_core_response_v1'
    mt=[]
    for name,c in cases.items():
        for metal,s in c['states'].items():
            t={'task_id':name+'_'+metal,'case_id':name,'metal':metal,'kind':'core','variant':'primary',**s}
            t['cache_key']=cache_key({'task':t,'model':model,'software':parent['software'],'implementation':pins});mt.append(t)
    mm={'schema_version':'alquemia.response_probe_core.v1','protocol_id':'mace_polar_archived_donor_response_probe_v1',
        'preparation':record(pp),'agreement':record(agreement),'parent_manifest':record(parent_path),'model':model,
        'software':parent['software'],'implementation':pins,'tasks':mt,'tolerances':parent['tolerances']}
    write_new(md/'manifest.json',mm)
    old=read_json(root/'workspaces/ggr_mechanism_20260915/stage_b_tasks_v1/manifest.json')
    dd=out/'dft';dd.mkdir();dt=[]
    for name,c in cases.items():
        for metal,s in c['states'].items():
            td=dd/(name+'_'+metal);td.mkdir();xp=td/'core.xyz';shutil.copyfile(verify(s['xyz']),xp)
            if name=='8GY2':
                text=verify(s['baseline_input']).read_text()
                text=re.sub(r'(?m)^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+\S+\s*$',
                            lambda m:f'* xyzfile {m[1]} {m[2]} core.xyz',text)
            else:text=f"! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF\n* xyzfile {s['charge']} 1 core.xyz\n"
            ip=td/'endpoint.inp';ip.write_text(text.rstrip()+'\n')
            dt.append({'task_id':name+'_'+metal,'case':name,'metal':metal,'charge':s['charge'],'multiplicity':1,
                'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),
                'source_input':s['baseline_input'],'source_xyz':s['xyz'],'source_manifest':c['preparation'],
                'task_type':'single_point','protocol_id':c['source_protocol']})
    dm={'schema_version':'alquemia.response_probe_dft.v1','protocol_id':'response_probe_matched_native_endpoints_v1',
        'agreement':record(agreement),'preparation':record(pp),'orca':old['orca'],
        'execution_policy':{name:record(root/'scripts'/filename) for name,filename in (
            ('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py'))},
        'execution_resources':{'mpi_ranks':16,'concurrent_tasks':4},'tasks':dt,'baseline_changed':False}
    write_new(dd/'manifest.json',dm)
    from affordable_workflow import dry_run
    return {'preparation':record(pp),'MACE_manifest':record(md/'manifest.json'),'DFT':dry_run(dd/'manifest.json'),
            'tasks':{'MACE':len(mt),'DFT':len(dt)}}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','agreement','output'):p.add_argument('--'+name,required=True)
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
