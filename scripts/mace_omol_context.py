"""Cached partition prerequisite for a static DFT/masked-MACE context descriptor."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new,xyz
from affordable_response import source_key
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL
from mace_omol_response import validate,collect,CASES

PROTOCOL='masked_omol_subtractive_context_DFT_CPCM_v1'
PARTITION_TOLERANCE=2.


@cached_file_checks
def partition(response_report,agreement,output):
    saved=read_json(response_report);mp=verify(saved['manifest']);validate(mp);actual=collect(mp)
    if (saved['status']!='complete' or not saved['gates']['analytic_derivative']
            or any(saved[k]!=v for k,v in actual.items())):
        raise InvalidArtifact('actual numerically qualified response receipts required')
    m=read_json(mp);p=read_json(verify(m['prepared']));tasks={t['task_id']:t for t in m['tasks']}
    rows={};mapping={}
    for name in CASES:
        c=read_json(verify(p['cases'][name]));physical=read_json(verify(c['physical_atoms']))
        full=p['physical_cases'][c['global_id']];atoms_by_key={a['source_key']:i for i,a in enumerate(physical)}
        if len(atoms_by_key)!=len(physical):raise InvalidArtifact('nonbijective physical inventory')
        endpoints={};errors=[];caps=0;source_atoms=0
        for metal in ('Ca','La'):
            key=f'{name}_center_{metal}';t=tasks[key];q=xyz(verify(t['xyz']))
            f=xyz(verify(full['grids']['center']['endpoints'][metal]['xyz']))
            if len(f)!=len(physical):raise InvalidArtifact('full atom inventory changed')
            for a,v in zip(physical,f):
                if ((metal if a['element']=='M' else a['element'])!=v[0]
                        or np.max(np.abs(np.array(a['xyz_A'])-v[1:]))>1e-6):
                    raise InvalidArtifact('full coordinates differ from source physical inventory')
            errors.append(float(np.max(np.abs(np.array(q[0][1:])-f[atoms_by_key['metal']][1:]))))
            for a in c['source_graph']['source_to_qm']:
                i=a['qm_index']
                if a['kind']=='source':
                    j=atoms_by_key[source_key(a['source'])];wanted=np.array(f[j][1:])
                    if q[i][0]!=f[j][0]:raise InvalidArtifact('shared source element differs')
                    source_atoms+=1
                elif a['kind']=='sigma_link_H':
                    x,y=(np.array(f[atoms_by_key[source_key(a[k])]][1:]) for k in ('retained','omitted'))
                    wanted=x+a['length_A']*(y-x)/np.linalg.norm(y-x);caps+=1
                    if q[i][0]!='H':raise InvalidArtifact('boundary cap is not hydrogen')
                else:raise InvalidArtifact('unsupported core boundary mapping')
                errors.append(float(np.max(np.abs(np.array(q[i][1:])-wanted))))
            if max(errors)>1e-6:raise InvalidArtifact('core/full source or cap coordinate mismatch')
            endpoints[metal]={'DFT_energy_hartree':t['DFT_source']['energy_hartree'],
                'DFT_source':t['DFT_source'],'masked_energy_eV':actual['rows'][key]['energy_eV'],
                'masked_result':actual['rows'][key], 'core_xyz':t['xyz'],
                'full_xyz':full['grids']['center']['endpoints'][metal]['xyz']}
        dr=(endpoints['Ca']['DFT_energy_hartree']-endpoints['La']['DFT_energy_hartree'])*HA_TO_KCAL
        tr=(endpoints['Ca']['masked_energy_eV']-endpoints['La']['masked_energy_eV'])*EV_TO_KCAL
        rows[name]={'endpoints':endpoints,'DFT_core_R_kcal_mol':dr,'masked_core_R_model_kcal':tr,
                    'core_anchor_R_kcal_scale':dr-tr,'evidence':c['evidence'],'global_id':c['global_id']}
        mapping[name]={'physical_atoms':c['physical_atoms'],'source_mapping':p['cases'][name],
                       'source_atoms_checked_across_both_endpoints':source_atoms,'caps_checked_across_both_endpoints':caps,
                       'maximum_coordinate_error_A':max(errors),'original_source_H_preserved':True}
    for metal in ('Ca','La'):
        if rows['GGR_extended']['endpoints'][metal]['full_xyz']!=rows['GGR_connected']['endpoints'][metal]['full_xyz']:
            raise InvalidArtifact('GGR partitions do not share exactly the same physical protein')
    d=rows['GGR_connected']['DFT_core_R_kcal_mol']-rows['GGR_extended']['DFT_core_R_kcal_mol']
    t=rows['GGR_connected']['masked_core_R_model_kcal']-rows['GGR_extended']['masked_core_R_model_kcal']
    result={'protocol_id':PROTOCOL,'status':'complete','response_report':record(response_report),
        'agreement':record(agreement),'implementation':record(__file__),'rows':rows,'mapping_checks':mapping,
        'partition_convention':'connected_minus_extended','DFT_partition_shift_kcal_mol':d,
        'masked_core_partition_shift_model_kcal':t,'hybrid_partition_shift_kcal_scale':d-t,
        'partition_tolerance_kcal_scale':PARTITION_TOLERANCE,'partition_gate_pass':abs(d-t)<=PARTITION_TOLERANCE,
        'whole_inference_eligible':abs(d-t)<=PARTITION_TOLERANCE,'new_model_calls':0,'new_DFT_calls':0,
        'whole_context_scores':None,'calibrated_class':None,'predictive_improvement_claimed':False,'baseline_changed':False}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('response-report','agreement','output'):p.add_argument('--'+key,required=True)
    r=partition(**vars(p.parse_args()))
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','mapping_checks')},indent=2))
