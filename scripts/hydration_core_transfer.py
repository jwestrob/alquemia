"""Score improved water hydrogens in the unchanged original amide-v3 cores."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
import hydration_square as square
from hydration_network import source_id, write_xyz

PROTOCOL='amide_v3_native_r2scan3c_context_prepared_water_H_v1'


def transfer_water_hydrogens(parent,original,context,proposed):
    """Source-atom identity maps water H; all other coordinates stay unchanged."""
    result=list(original)
    lookup={source_id(a['source']):a['qm_index'] for a in context['atom_graph']['source_to_qm'] if a['kind']=='source'}
    changed=[]
    for a in parent['atom_graph']['source_to_qm']:
        if a['kind']!='source' or a['source']['canonical_resname']!='HOH' or a['source']['element']!='H':continue
        key=source_id(a['source'])
        if key not in lookup:raise InvalidArtifact('water H missing from physical context')
        i=a['qm_index'];j=lookup[key]
        if original[i][0]!='H' or proposed[j][0]!='H':raise InvalidArtifact('water mapping changed element')
        result[i]=proposed[j];changed.append(i)
    for w in square.water_groups(parent,original):
        key=source_id(next(a['source'] for a in parent['atom_graph']['source_to_qm'] if a['qm_index']==w['oxygen_index']))
        if original[w['oxygen_index']]!=proposed[lookup[key]]:raise InvalidArtifact('physical context moved water oxygen')
    if any(result[i]!=original[i] for i in range(len(original)) if i not in changed):
        raise InvalidArtifact('non-water-H coordinates changed')
    return result,changed


def prepare(config,adjudication,ggr_release,agreement,output):
    cfg=read_json(config);source=read_json(verify(cfg['source_config']))
    a=read_json(adjudication)
    if a['status']!='complete' or not a['all_proposals_lower_both_DFT_starts']:
        raise InvalidArtifact('actual DFT adjudication must support improved starting geometry')
    am=read_json(verify(a['manifest']));mc=read_json(verify(am['mace_collection']));mm=read_json(verify(mc['manifest']))
    contexts={}
    for pin in mm['source_manifests']:
        m=read_json(verify(pin))
        for p in m['preparations']:
            c=read_json(verify(p));contexts[c['case']]=c
    tasks=[];out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for item in source['parents']:
        parent,_,waters=square.checked_parent(item);case=item['case']
        for metal in ('Ca','La'):
            target=next(t for t in am['tasks'] if t['case']==case and t['metal']==metal)
            original=xyz(verify(parent['outputs'][metal]['xyz']));proposed=xyz(verify(target['xyz']))
            rows,mobile=transfer_water_hydrogens(parent,original,contexts[case],proposed)
            if len(mobile)!=2*len(waters):raise InvalidArtifact('unexpected water H inventory')
            d=out/'tasks'/(case+'__'+metal);d.mkdir(parents=True);xp=d/'core.xyz';write_xyz(xp,rows,PROTOCOL)
            q=parent['outputs'][metal]['charge'];ip=d/'endpoint.inp';ip.write_text(square.METHOD+f'\n* xyzfile {q} 1 core.xyz\n')
            tasks.append({'task_id':case+'__'+metal,'case':case,'metal':metal,'charge':q,'multiplicity':1,
                'input':record(ip),'xyz':record(xp),'output_path':str(d/'endpoint.out'),
                'original_xyz':parent['outputs'][metal]['xyz'],'parent':item['preparation'],
                'transferred_context_xyz':target['xyz'],'water_H_indices':mobile,
                'atom_count':len(rows),'original_endpoints':item['endpoints'],
                'cache_key':cache_key({'protocol':PROTOCOL,'input':record(ip),'xyz':record(xp),'context':target['xyz']})})
    release=read_json(ggr_release);ggr=next(r for r in release['scores'] if r['case']=='ggr_1glg_GGR' and r['lane']=='repaired')
    gp=read_json(verify(ggr['source_manifest']))
    if gp['explicit_water_inventory']:raise InvalidArtifact('declared GGR dry identity does not hold')
    for metal in ('Ca','La'):
        t=ggr['endpoints'][metal];square.endpoint(t['output'],t['receipt'],gp['outputs'][metal]['xyz'],gp['outputs'][metal]['input'])
    impl=out/'implementation';impl.mkdir();pins={}
    for name in ('hydration_core_transfer.py','hydration_network.py','hydration_square.py','affordable_peptide.py'):
        p=impl/name;shutil.copyfile(Path(__file__).with_name(name),p);pins[name]=record(p)
    result={'protocol_id':PROTOCOL,'tasks':tasks,'agreement':record(agreement),'configuration':record(config),
        'orca':source['orca'],'execution_policy':source['execution_policy'],'adjudication':record(adjudication),
        'ggr_release':record(ggr_release),'ggr_reused':ggr,'implementation':pins,
        'baseline_changed':False,'PQQ_inputs_changed':False,'new_quantum_endpoints':4}
    write_new(out/'manifest.json',result)
    from affordable_workflow import dry_run
    return dry_run(out/'manifest.json')


def compare_ggr_replicates(collection,ggr_study,output):
    result=read_json(collection);study=read_json(ggr_study);comparisons=[];controls=[]
    if result['status']!='complete':raise InvalidArtifact('complete prepared-core result required')
    for entry in study['source_by_representation']:
        if entry['representation']!='formamide':continue
        archive=read_json(verify(entry['collection']))
        row=next(r for r in archive['rows'] if r['case']==entry['case'] and r['protocol_id']=='generic_peptide_amide_vertical_native_r2scan3c_v3')
        parent=read_json(verify(row['source_manifest']))
        if parent['explicit_water_inventory']:raise InvalidArtifact('GGR dry preparation identity failed')
        energies={metal:square.endpoint(e['output'],e['receipt'],parent['outputs'][metal]['xyz'],parent['outputs'][metal]['input'])['energy_hartree'] for metal,e in row['endpoints'].items()}
        r=energies['Ca']-energies['La']
        if abs(r-entry['R_hartree'])>1e-10:raise InvalidArtifact('archived GGR contrast differs')
        controls.append({'source':entry['source'],'R_hartree':r,'collection':entry['collection'],'source_manifest':row['source_manifest']})
        for a in result['cases']:
            comparisons.append({'alpha_structure':a['case'],'GGR_structure':entry['source'],
                'before_kcal_mol':(a['original_R_hartree']-r)*HA_TO_KCAL,
                'after_kcal_mol':(a['prepared_R_hartree']-r)*HA_TO_KCAL,
                'direction_reproduced':a['prepared_R_hartree']>r})
    if {c['source'] for c in controls}!={'1GLG','2FW0','2FVY'}:raise InvalidArtifact('all three archived GGR sources required')
    combined={'collection':record(collection),'ggr_study':record(ggr_study),'controls':controls,'comparisons':comparisons,
        'correct_before':sum(c['before_kcal_mol']>0 for c in comparisons),'correct_after':sum(c['direction_reproduced'] for c in comparisons),
        'comparisons_count':len(comparisons),'independent_biological_comparisons':1,
        'minimum_after_margin_kcal_mol':min(c['after_kcal_mol'] for c in comparisons),
        'evidence_role':'additional previously inspected structural robustness controls; retrospective development',
        'new_quantum_calls':0,'threshold_fit':False}
    write_new(output,combined);return combined


def collect(manifest,output):
    m=read_json(manifest);rows=[];g=m['ggr_reused']
    gr={metal:square.endpoint(e['output'],e['receipt'])['energy_hartree'] for metal,e in g['endpoints'].items()}
    ggr_r=gr['Ca']-gr['La']
    for t in m['tasks']:
        try:
            op=Path(t['output_path']);r=square.endpoint(record(op),record(str(op)+'.execution.json'),t['xyz'],t['input'])
            rows.append({'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'status':'complete','result':r})
        except (OSError,ValueError,KeyError) as exc:rows.append({'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'status':'unavailable','failure':str(exc)})
    cases=[]
    for case in ('1F6S','6IP9'):
        pair={r['metal']:r for r in rows if r['case']==case}
        if len(pair)!=2 or any(r['status']!='complete' for r in pair.values()):continue
        r=pair['Ca']['result']['energy_hartree']-pair['La']['result']['energy_hartree']
        original=next(t['original_endpoints'] for t in m['tasks'] if t['case']==case)
        old={metal:square.endpoint(e['output'],e['receipt'])['energy_hartree'] for metal,e in original.items()}
        old_r=old['Ca']-old['La']
        cases.append({'case':case,'original_R_hartree':old_r,'prepared_R_hartree':r,
            'water_preparation_delta_R_kcal_mol':(r-old_r)*HA_TO_KCAL,
            'original_alpha_minus_GGR_kcal_mol':(old_r-ggr_r)*HA_TO_KCAL,
            'prepared_alpha_minus_GGR_kcal_mol':(r-ggr_r)*HA_TO_KCAL,
            'direction_reproduced':r>ggr_r})
    result={'manifest':record(manifest),'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
            'rows':rows,'cases':cases,'ggr_R_hartree':ggr_r,'both_development_directions_reproduced':len(cases)==2 and all(c['direction_reproduced'] for c in cases),
            'biological_groups':'one consumed alpha-lactalbumin group versus one GGR group; not two independent successes',
            'baseline_changed':False,'calibrated_decision':None,'occupancy_probabilities':None}
    write_new(output,result);return {'status':result['status'],'cases':cases,'both_development_directions_reproduced':result['both_development_directions_reproduced']}


if __name__=='__main__':
    import json
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    q=s.add_parser('prepare')
    for name in ('config','adjudication','ggr-release','agreement','output'):q.add_argument('--'+name,required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('compare-ggr-replicates')
    for name in ('collection','ggr-study','output'):q.add_argument('--'+name,required=True)
    a=vars(p.parse_args());op=a.pop('op').replace('-','_');print(json.dumps(globals()[op](**a),indent=2))
