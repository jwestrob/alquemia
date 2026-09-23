"""Existing-output A8 audit: all four triples, no scientific execution."""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import re
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from accommodation_nonlinear import contrast_components,relative_components
from structure_informed_starts import scf_details
from consistent_context import fragment_id


@lru_cache(None)
def loaded(path):return read_json(path)


def load(pin):return loaded(str(verify(pin)))


def logical_pool(collection,cid):
    d=load(collection)
    return next(r for r in d['cases'] if r['source'].get('actual_union_case_id',r['source'].get('source_case_id',r['case_id']))==cid)


def scalar(cell,metal,candidate,medium):
    low=cell['low'][medium];mp=load(low['manifest']);t=next(t for t in mp.get('all_tasks',mp['tasks']) if t['task_id']==low['task_id'])
    text=verify(low['output']).read_text();details=scf_details(text);receipt=load(low['receipt']);root=verify(low['output']).parent
    trace=re.findall(r'^\s*\d+\s+[-+]?\d+\.\d{10,}\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*$',text,re.M)
    seeds={k:record(root/('endpoint.runtime.'+k)) if (root/('endpoint.runtime.'+k)).is_file() else None for k in ('gbw','xtbw')}
    return {'metal':metal,'candidate':candidate,'medium':medium,'actual':low,'xyz':t['xyz'],'charge':t['charge'],
        'multiplicity':t['multiplicity'],'atoms':len(xyz(verify(t['xyz']))),'input':t['input'],
        'raw_input':verify(t['input']).read_text(),'normal_termination':'ORCA TERMINATED NORMALLY' in text,
        'energy_check_stop':'Energy Check signals convergence' in text,'SCF':details,'last_five_printed_iterations':trace[-5:],
        'retained_restart_seeds':seeds,'runtime_input':receipt['artifacts']['runtime_input']}


def summarize(c,collection):
    m=load(load(collection)['manifest']);sm=load(m['source_manifest'])
    cid=c['case_id'];tasks=[t for t in sm['tasks'] if t['case_id']==cid]
    proposals={}
    for candidate in c['candidates']:
        if candidate['id']=='origin':continue
        rp=verify(candidate['xyz']).parents[2]/'result.json';r=loaded(str(rp))
        if r['proposal']['coordinate']!=candidate['xyz']:raise InvalidArtifact('proposal coordinate mismatch')
        metal=candidate['id'].split('_')[-1];task=next(t for t in tasks if t['metal']==metal)
        proposals[candidate['id']]={'receipt':record(rp),'full_q':r['proposal']['full_q'],
            'active_mode_ids':task['active_mode_ids'],'active_q_radian':r['proposal']['active_q_radian'],
            'mapping':task['mapping'],'boundary_flag':r['boundary_flag'],'final_geometry':r['final_geometry']}
    return {'collection':collection,'case_id':cid,'source_manifest':m['source_manifest'],'proposals':proposals,
        'selected_modes':{t['metal']:t['active_mode_ids'] for t in tasks},'pool':c['pool'],'matrix':c['matrix'],
        'origin_components':contrast_components(c['matrix']['Ca']['origin']['components'],c['matrix']['La']['origin']['components']),
        'candidate_contrasts':{name:contrast_components(c['matrix']['Ca'][name]['components'],c['matrix']['La'][name]['components']) for name in c['matrix']['Ca']},
        'selected_work':{z:relative_components(c['matrix'][z][c['pool']['rows'][z]['operational_candidate']]['components'],c['matrix'][z]['origin']['components']) for z in ('Ca','La')},
        'scalar_cells':[scalar(cell,z,name,medium) for z,row in c['matrix'].items() for name,cell in row.items() for medium in ('vacuum','alpb')]}


def audit(comparison,continuing_inventory,output):
    d=read_json(comparison);a=load(d['audit']);prior=load(d['prior']);continuing=read_json(continuing_inventory)
    rows=[r for r in d['rows'] if r['case_id'].startswith('a8r3s4-')];triples=[t for t in d['triples'] if t['protein_id']=='a8r3s4-pqq-la_model']
    if len(rows)!=7 or len(triples)!=4:raise InvalidArtifact('all seven pairs/four A8 triples required')
    prepared={(r['selection_id'],r['case_id']):r['source'] for r in a['rows']};oldrows={r['case_id']:r for r in prior['rows']};saved={};results=[]
    for r in rows:
        cid=r['case_id'];source=prepared[r['selection_id'],cid];old=oldrows[cid]
        if cid not in saved:
            c=logical_pool(old['precision_collection'],cid);saved[cid]=summarize(c,old['precision_collection'])
        newc=next(c for c in load(r['source_collection'])['cases'] if c['case_id']==r['pair_id']);new=summarize(newc,r['source_collection'])
        oldc=logical_pool(old['precision_collection'],cid)
        before={fragment_id(x) for x in load(oldc['source']['union']['union'])['fragments']}
        after={fragment_id(x) for x in load(source['union'])['fragments']}
        result={'pair_id':r['pair_id'],'case_id':cid,'selection_id':r['selection_id'],'expected_class':r['expected_class'],
            'before_atoms':r['tenfold_atom_count'],'after_atoms':r['atom_count'],'before_La_charge':source['tenfold_comparison']['La_charge'],
            'after_La_charge':r['La_charge'],'removed_fragments':sorted(before-after),'added_fragments':sorted(after-before),
            'source_preparation':source,'before':saved[cid],'after':new,'methods':r['methods'],
            'R_component_changes':{k:new['pool']['operational'][k]-saved[cid]['pool']['operational'][k] for k in ('native_R_model_kcal_mol','GFN2_vacuum_R_kcal_mol','GFN2_ALPB_R_kcal_mol','solvation_delta_R_kcal_mol','composite_R_model_kcal_mol')}}
        results.append(result)
    ref=load(d['reference']);canon=next(r for r in ref['rows'] if r['case_id'].startswith('a8r3'))
    cr=next(c for c in load(canon['source_collection'])['cases'] if c['case_id']==canon['case_id']);canonical=summarize(cr,canon['source_collection'])
    ongoing=[r for r in continuing['rows'] if r['case_id']=='a8r3s4-pqq-la_model']
    result={'comparison':record(comparison),'rows':results,'all_four_triples':triples,'changed_threefold_canonical':canonical,
        'ongoing_inventory':record(continuing_inventory),'ongoing_A8_cells':ongoing,
        'ongoing_A8_atom_counts':sorted({len(xyz(verify(r['xyz']))) for r in ongoing}),
        'changed_canonical_atom_count':canon['source_preparation']['atom_count'],
        'ongoing_case_is_same_as_any_changed_pool':any(r['xyz']['sha256']==x['xyz']['sha256'] for r in ongoing for c in results for x in c['after']['scalar_cells']),
        'new_molecular_calls':0,'changed_parameters':False,'production_changed':False,'implementation':record(__file__)}
    write_new(output,result);return {'output':record(output),'all_four_triples':len(triples),'source_membership_pairs':len(results),
        'ongoing_atom_counts':result['ongoing_A8_atom_counts'],'changed_canonical_atoms':result['changed_canonical_atom_count'],'new_molecular_calls':0}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('comparison','continuing_inventory','output'):p.add_argument('--'+k.replace('_','-'),required=True)
    print(json.dumps(audit(**vars(p.parse_args())),indent=2))
