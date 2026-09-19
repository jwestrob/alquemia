#!/usr/bin/env python3
"""Opt-in scoring of explicit prepared cores with frozen PQQ response weights."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import time

from affordable_common import (InvalidArtifact,HA_TO_KCAL,record,verify,read_json,write_new,
                               xyz,digest,energy,classify_raw)
from response_probe_archive import (PROTOCOL,PRIMARY,paired_coordinates,native_gradient,
    canonical_mapping,generic_mapping,selected_donors,describe_pair)
from site_classifier import predict

CANONICAL='pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3'


def method_text(path):
    return re.sub(r'(?m)^\s*\*\s+xyzfile.*$','',Path(path).read_text()).strip()


def prepared_states(preparation):
    p=read_json(preparation);states={}
    for metal in ('Ca','La'):
        if p['protocol_id']==CANONICAL:
            output=p['outputs'];xp=output[metal+'_xyz'];ip=output[metal+'_input']
            ip={**ip,'path':str((Path(preparation).resolve().parent/ip['path']).resolve())}
            xp={**xp,'path':str((Path(preparation).resolve().parent/xp['path']).resolve())}
            match=re.search(r'^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+\S+',verify(ip).read_text(),re.M)
            if match is None:raise InvalidArtifact('missing canonical charge/multiplicity declaration')
            charge,mult=int(match[1]),int(match[2])
        else:
            output=p['outputs'][metal];xp=output['xyz'];ip=output['input']
            charge=output['charge'];mult=output['multiplicity']
        states[metal]={'xyz':xp,'charge':charge,'spin_multiplicity':mult,'baseline_input':ip}
    _,coordinates=paired_coordinates(states)
    if p['protocol_id']==CANONICAL:
        if p['pqq']['microstate_id']!='pqq_ox_3minus_v1' or p['fixed_core']['water_policy']!='dry_exclude_all_source_and_synthetic_waters':
            raise InvalidArtifact('unsupported canonical microstate/water policy')
        mapping,physical,cofactor,source=canonical_mapping(p,coordinates)
    else:mapping,physical,cofactor,source=generic_mapping(p,coordinates)
    donors=selected_donors(mapping,physical,coordinates,p)
    return p,states,coordinates,mapping,physical,cofactor,source,donors


def dft_pair(manifest, case_id, states, protocol):
    from run_orca_task_manifest import load_manifest_tasks,_completed_attempt_is_valid
    mp=Path(manifest).resolve();m,normalized=load_manifest_tasks(mp)
    normal={t['task_id']:t for t in normalized};result={}
    selected=[t for t in m['tasks'] if t.get('case',t.get('case_id'))==case_id]
    if len(selected)!=2 or {t['metal'] for t in selected}!={'Ca','La'}:
        raise InvalidArtifact('DFT manifest must supply exact paired case')
    for t in selected:
        metal=t['metal'];s=states[metal]
        if xyz(verify(t['xyz']))!=xyz(verify(s['xyz'])) or t['charge']!=s['charge'] or t['multiplicity']!=1:
            raise InvalidArtifact('DFT/MACE prepared-state mismatch')
        if protocol==CANONICAL:
            if t['protocol_id']!=CANONICAL or method_text(verify(t['input']))!=method_text(verify(s['baseline_input'])):
                raise InvalidArtifact('DFT not the canonical prepared Hamiltonian/input policy')
        elif method_text(verify(t['input']))!='! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF':
            raise InvalidArtifact('direct transfer requires the declared native TightSCF policy')
        op=Path(t['output_path']);rp=Path(str(op)+'.execution.json')
        if not rp.exists() or not op.exists():raise InvalidArtifact('DFT endpoint not completed: '+t['task_id'])
        if not _completed_attempt_is_valid(rp,op,manifest_sha256=digest(mp),task=normal[t['task_id']],
                    runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
            raise InvalidArtifact('DFT execution receipt failed compatibility validation')
        receipt=read_json(rp)
        if receipt['orca_version']!='6.1.1':raise InvalidArtifact('unqualified ORCA version')
        result[metal]={'energy_hartree':energy(op),'output':record(op),'execution_receipt':record(rp),
                       'input':t['input'],'xyz':t['xyz'],'method':method_text(verify(t['input']))}
    return result


def score(preparation,mace_collection,dft_manifest,case_id,evaluation,baseline_release,output):
    start,cpu=time.monotonic(),time.process_time()
    result={'protocol':PROTOCOL,'case_id':case_id,'status':'unavailable','baseline':None,'response_heads':None,
            'affinity_probability':None,'baseline_changed':False,'new_fits':0,'new_molecular_calls':0,
            'sources':{k:record(p) for k,p in [('preparation',preparation),('MACE_collection',mace_collection),
                ('DFT_manifest',dft_manifest),('evaluation',evaluation),('baseline_release',baseline_release)]}}
    try:
        p,states,coords,mapping,physical,cofactor,source,donors=prepared_states(preparation)
        fitted=read_json(evaluation);f=read_json(verify(fitted['features']))
        if fitted['protocol']!=PROTOCOL or f['protocol']!=PROTOCOL:raise InvalidArtifact('unqualified frozen model')
        dft=dft_pair(dft_manifest,case_id,states,p['protocol_id'])
        contrast=(dft['Ca']['energy_hartree']-dft['La']['energy_hartree'])*HA_TO_KCAL
        result['DFT_endpoints']=dft;result['features']={'DFT_R_kcal_mol':contrast}
        release=read_json(baseline_release)
        result['baseline']={'R_kcal_mol':contrast,'S_kcal_mol':None,'decision':'uncalibrated_protocol',
                            'source_protocol':p['protocol_id'],'reference_status':'unavailable_for_this_protocol'}
        canonical=p['protocol_id']==CANONICAL
        if canonical:
            if release['protocol_id']!=CANONICAL:raise InvalidArtifact('baseline reference protocol mismatch')
            gauge=release['aquo_reporting_gauge'];verify({'path':gauge['path'],'sha256':gauge['sha256']})
            result['baseline'].update(S_kcal_mol=contrast-gauge['A_kcal_mol'],
                decision=classify_raw(contrast,release,CANONICAL),reference_status='released_canonical_reporting_gauge',
                released_bands=release['calibration']['released_supported_bands'])
        collection=read_json(mace_collection);mm=read_json(verify(collection['manifest']))
        tasks={t['metal']:t for t in mm['tasks'] if t.get('case_id')==case_id}
        if set(tasks)!={'Ca','La'}:raise InvalidArtifact('MACE manifest must supply exact paired case')
        gradients={};forces={}
        for metal,t in tasks.items():
            native=collection['rows'][t['task_id']]
            gradients[metal]=native_gradient(native,states[metal],f['model'],f['software'])
            forces[metal]=native['forces']
        description=describe_pair(gradients,mapping,physical,coords,donors,cofactor)
        result['features'].update(description['features']);result['response_measurement']=description
        result.update(donors=donors,physical_mapping=mapping,force_artifacts=forces,source_structure=source)
        if canonical:
            row={'features':result['features'],'target':'PQQ_functional_class','DFT_protocol':CANONICAL,
                 'MACE_protocol':PROTOCOL}
            result['response_heads']={arm:{**predict(fitted['models'][arm],row),'logit_threshold':0.,
                'probability':None,'qualification':'research_PQQ_functional_or_association_transfer'} for arm in ('DFT','DFT_radial')}
        else:
            result['response_heads']={'status':'absolute_cross_target_classification_unavailable',
                                     'class':None,'logit':None,'probability':None}
        result.update(status='complete',source_protocol=p['protocol_id'],
                      target='PQQ_functional_or_association_transfer' if canonical else 'direct_cross_target_ordering_only',
                      relaxation_correction_kcal_mol=None,entropy_correction_kcal_mol=None)
    except (InvalidArtifact,KeyError,ValueError,StopIteration) as error:
        result.update(status='unavailable',reason=str(error))
    result.update(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,implementation=record(__file__))
    write_new(output,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('preparation','mace-collection','dft-manifest','case-id','evaluation','baseline-release','output'):
        p.add_argument('--'+name,required=True)
    result=score(**vars(p.parse_args()))
    print(json.dumps({k:result.get(k) for k in ('status','reason','case_id','features','baseline','response_heads')},indent=2))
    if result['status']!='complete':raise SystemExit(1)
