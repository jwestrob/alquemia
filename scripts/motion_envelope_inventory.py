"""Read-only envelope archive identity audit and finite native-MACE origin inputs."""
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path
import time
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from adaptive_origin_recovery import snapshot
import consistent_context as context
import motion_envelope as envelope
import union_adaptive


def inventory(preparation,archive_spec,output):
    start=time.monotonic();p=read_json(preparation);m=envelope.validate(verify(p['selection']));spec=read_json(archive_spec)
    source=read_json(verify(m['source_inputs']));old_selection=read_json(verify(m['prior_selection']))
    oldten=read_json(verify(old_selection['tenfold_preparation']));ten={r['case_id']:r for r in oldten['cases']}
    canonical=read_json(verify(m['canonical_selection']));cr=read_json(verify(canonical['crystals']));ten.update({r['case_id']:r for r in cr['cases']})
    root_actual={x['protein_id']:x['case_id'] for x in canonical['canonical_sources']}
    origins=defaultdict(list);pools=defaultdict(list);archive_counts={}
    for pin in spec['pool_collections']:
        d=read_json(verify(pin));n=0
        for c in d.get('cases',[]):
            old=c.get('source',{}).get('union')
            if not old or c.get('pool',{}).get('status')!='available':continue
            cid=old['case_id'];pools[cid].append((old,c,pin));n+=1
            ep=c['source'].get('origin_row',{}).get('native_endpoints',{})
            if len(ep)==2:origins[cid].append((ep,pin,old['state_key']))
        archive_counts[pin['path']]=n
    for pin in spec['origin_collections']:
        d=read_json(verify(pin))
        for r in d['rows']:
            if r.get('status')!='complete' or len(r.get('native_endpoints',{}))!=2:continue
            cid=r.get('source_case_id',r['case_id']);cid=root_actual.get(cid,cid)
            state=r.get('source',{}).get('state_key');origins[cid].append((r['native_endpoints'],pin,state))
    for pin in spec.get('local_inventories',[]):
        d=read_json(verify(pin))
        for c in d['cases']:
            ep=c.get('representations',{}).get('context',{}).get('endpoints',{})
            if len(ep)==2:origins[root_actual.get(c['case_id'],c['case_id'])].append((ep,pin,None))
    rows=[]
    for r in p['cases']:
        cid=r['case_id'];entry={'case_id':cid,'selection_id':r['selection_id'],'roles':r['roles'],'status':r['status'],
            'origin_reuse':None,'pool_reuse':None,'candidate_coordinate_matches_requiring_further_state_audit':[]}
        if r['status']=='prepared':
            eps=r['representations']['context']['endpoints'];old=ten.get(cid)
            now={context.fragment_id(f) for f in read_json(verify(r['union']))['fragments']}
            previous={context.fragment_id(f) for f in read_json(verify(old['union']))['fragments']} if old else set()
            entry.update(atom_count=r['atom_count'],La_charge=r['new_La_charge'],cap_count=r['cap_count'],
                old_local_atom_count=r['original_atom_count'],old_local_La_charge=r['original_La_charge'],
                tenfold_atom_count=old.get('atom_count') if old else None,tenfold_La_charge=old.get('new_La_charge') if old else None,
                fragments_added_to_tenfold=sorted(now-previous),fragments_removed_from_tenfold=sorted(previous-now),
                tenfold_same_state=old is not None and old.get('state_key')==r['state_key'],maps=r['maps'])
            for ep,pin,state in origins[cid]:
                if all(context.reusable_state(eps[z],ep[z]) for z in ('Ca','La')):
                    if state!=r['state_key']:
                        entry['candidate_coordinate_matches_requiring_further_state_audit'].append(pin);continue
                    actual={}
                    for z in ('Ca','La'):
                        receipt,forces=union_adaptive.origin(ep[z],m['config']['model']);actual[z]=ep[z]
                    entry['origin_reuse']={'source':pin,'endpoints':actual,'exact_graph_state':True};break
            for old,c,pin in pools[cid]:
                if old['state_key']==r['state_key'] and old['original_core']==r['original_core'] and all(context.reusable_state(eps[z],old['representations']['context']['endpoints'][z])for z in ('Ca','La')):
                    entry['pool_reuse']={'status':'physical_match_requires_selector_and_actual_candidate_audit','source':pin,'source_case_id':c['case_id']};break
        rows.append(entry)
    ri={(r['selection_id'],r['case_id']):r for r in rows};pilot=[]
    for item in p['pilot']:
        r=ri.get((item['selection_id'],item['case_id']));pilot.append({**item,'audit':r})
    supported=[x for x in pilot if x['status']=='prepared'];unresolved=sum(bool(x['audit']['candidate_coordinate_matches_requiring_further_state_audit'])or bool(x['audit']['pool_reuse']) for x in supported)
    if unresolved:raise InvalidArtifact('exact coordinate/pool matches need further audit before finite new-call counts')
    counts={'pilot_denominator':34,'pilot_supported':len(supported),'whole_pool_reuses':0,
        'origin_pairs_reused':sum(x['audit']['origin_reuse'] is not None for x in supported),
        'new_origin_MACE_calls':sum(2 for x in supported if not x['audit']['origin_reuse']),
        'new_optimizer_starts':2*len(supported),'maximum_cross_MACE_calls':2*len(supported),
        'prospective_origin_GFN2_cells_without_scalar_reuse':4*len(supported),'prospective_candidate_GFN2_cells_without_scalar_reuse':8*len(supported)}
    result={'protocol_id':envelope.PROTOCOL,'preparation':record(preparation),'selection':p['selection'],'archive_spec':record(archive_spec),
        'archive_pool_counts':archive_counts,'rows':rows,'pilot':pilot,'counts':counts,
        'GFN2_reuse_status':'not_authorized_under_unselected_new_scalar_numerical_profile','new_molecular_calls':0,
        'new_reference':None,'production_changed':False,'wall_seconds':time.monotonic()-start,'implementation':record(__file__)}
    write_new(output,result);return {'inventory':record(output),'counts':counts}


def mace_inputs(inventory,template,output):
    audit=read_json(inventory);prep=read_json(verify(audit['preparation']));sel=envelope.validate(verify(audit['selection']));base=read_json(template)
    if sel['config']['model']!=base['model'] or sel['config']['software']!=base['software']:raise InvalidArtifact('native checkpoint/runtime differs')
    if audit['counts']['pilot_denominator']!=34:raise InvalidArtifact('fixed pilot scope differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=snapshot(Path(__file__).parent,out/'implementation')
    selected={(x['selection_id'],x['case_id']):x for x in prep['cases']};tasks=[];sources=[]
    for p in audit['pilot']:
        if p['status']!='prepared':continue
        r=selected[p['selection_id'],p['case_id']];sources.append(r)
        if p['audit']['origin_reuse']:continue
        cid=r['case_id']+'__envelope_'+r['selection_id'].rsplit('_',1)[-1]
        for z,ep in r['representations']['context']['endpoints'].items():
            task={'task_id':cid+'__context__'+z,'case_id':cid,'source_case_id':r['case_id'],'selection_id':r['selection_id'],
                'representation':'context','metal':z,'metal_index':0,'kind':'core','xyz':ep['xyz'],'charge':ep['charge'],
                'spin_multiplicity':ep['multiplicity'],'energy_component':'MACE_OMOL_total_vacuum_energy'}
            task['cache_key']=cache_key({'task':task,'model':base['model'],'software':base['software'],'implementation':impl});tasks.append(task)
    if len(tasks)!=audit['counts']['new_origin_MACE_calls'] or len({t['task_id'] for t in tasks})!=len(tasks):raise InvalidArtifact('finite native origin count differs')
    inputs={'protocol_id':envelope.PROTOCOL,'inventory':record(inventory),'preparation':audit['preparation'],'agreement':sel['agreement'],
        'cases':sources,'model':base['model'],'software':base['software'],'implementation':impl,'new_GFN2_calls':0,
        'composite_status':'unavailable_until_real_scalar_endpoints','reference':None}
    write_new(out/'INPUTS.json',inputs)
    manifest={'protocol_id':'mace_omol_0_100m_vacuum_descriptor_v1','stage':envelope.PROTOCOL,'preparation':record(out/'INPUTS.json'),
        'agreement':sel['agreement'],'software':base['software'],'model':base['model'],'implementation':impl,
        'tasks':tasks,'reused':{},'verification_references':{}}
    write_new(out/'mace'/'manifest.json',manifest)
    from mace_omol import check_atoms
    for t in tasks:check_atoms(xyz(verify(t['xyz'])),t['charge'])
    ready={'inputs':record(out/'INPUTS.json'),'MACE_manifest':record(out/'mace'/'manifest.json'),'new_origin_MACE_calls':len(tasks),
        'new_GFN2_calls':0,'new_optimizer_starts':0,'gpu_python':base['gpu_python'],'cpu_python':base['cpu_python'],
        'actual_energy_execution':False,'allocation':{'CPU':32,'memory_MiB':200000,'GPUs':1},'production_changed':False}
    write_new(out/'READY.json',ready);return ready

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    for op,keys in {'inventory':('preparation','archive_spec','output'),'mace_inputs':('inventory','template','output')}.items():
        q=s.add_parser(op)
        for k in keys:q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());cmd=a.pop('command');print(json.dumps(globals()[cmd](**a),indent=2))
