"""Opt-in fixed4.3Å source-fragment envelope; preparation and archive inventory only."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import copy
import json
from pathlib import Path
import time
import types

from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
import consistent_context as context
import second_shell_context as shell
from union_triple_transfer import mapped,mapping_equivalence

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL='Nikasha_three_La_fixed_polar_motion_envelope4p3_preparation_v1'
POLICY={**shell.POLICY,'polar_neighbor_cutoff_A':4.3}
PROBES=[('a0a3f2yly8','Ca',1),('a0a3f2yly8','Ca',3),('a0acd6b9f2','Ca',4),('a0acd6b9f2','La',4),('a8r3s4','La',1),('a8r3s4','La',3)]


def isolated_discover(state):
    """Exact existing kernel, private globals; never mutate shared module POLICY."""
    globals_=dict(shell.discover.__globals__);globals_['POLICY']=copy.deepcopy(POLICY)
    return types.FunctionType(shell.discover.__code__,globals_,'motion_envelope_discover')(state)


def discover_source(row,config):
    result={'case_id':row['case_id'],'status':'unsupported','source':row['source']}
    try:
        state=shell.parent_state(row['core'],config['topology'],require_endpoint_receipts=False)
        old=read_json(verify(row['representations']['context']['preparation']))
        anchors,fragments,contacts,excluded=isolated_discover(state)
        actual={shell.atom_key(x) for x in old['anchors']}
        if anchors!=actual:raise InvalidArtifact('original donor-anchor identity does not replay exactly')
        old_keys={context.fragment_id(x) for x in old['added_fragments']}
        records=[{'kind':kind,'source':state['graph'].residues[rkey].source_dict()} for kind,rkey in sorted(fragments)]
        new_keys={context.fragment_id(x) for x in records}
        if not old_keys<=new_keys:raise InvalidArtifact('larger source envelope lost old selected fragment')
        result.update(status='discovered',anchors=[state['graph'].meta[k] for k in sorted(anchors)],
            original_preparation=row['representations']['context']['preparation'],fragments=records,
            added_fragment_ids=sorted(new_keys-old_keys),contacts=contacts,excluded_contacts=excluded,
            original_radius_A=3.5,envelope_radius_A=4.3,anchor_identity_replayed=True)
    except Exception as exc:result.update(reason=str(exc),exception_type=type(exc).__name__)
    return result


def select(prior_selection,canonical_selection,agreement,output,workers=4):
    start=time.monotonic();prior=read_json(prior_selection);canonical=read_json(canonical_selection)
    original=read_json(verify(prior['source_preparation']));config=original['config']
    if len(prior['triples'])!=100 or canonical['denominator']!=28 or canonical['primary_selection']!=record(prior_selection):
        raise InvalidArtifact('original100/canonical28 selection linkage differs')
    if sum(bool(t['missing_members']) for t in prior['triples'])!=6:raise InvalidArtifact('six original missing triples required')
    src={r['case_id']:r for r in original['cases']}
    oldcr=read_json(verify(canonical['crystals']));rawcr=read_json(verify(oldcr['source_preparation']))
    crystals={r['case_id']:copy.deepcopy(r) for r in rawcr['cases'] if r['case_id'] in ('1H4I','4MAE','1KB0')}
    if len(crystals)!=3:raise InvalidArtifact('original three crystal sources missing')
    oldcri={r['case_id']:r for r in oldcr['cases']}
    for cid,r in crystals.items():r['source']=copy.deepcopy(oldcri[cid]['source'])
    src.update(crystals)
    ids=sorted({c for t in prior['triples'] if not t['missing_members'] for c in t['members']}|set(crystals))
    out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('outputs must be under workspaces')
    out.mkdir(parents=True,exist_ok=False)
    discovery={}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fs={ex.submit(discover_source,src[cid],config):cid for cid in ids}
        for f in as_completed(fs):
            r=f.result();discovery[r['case_id']]=r;print(json.dumps({'discovery':r['case_id'],'status':r['status'],'reason':r.get('reason')}),flush=True)
    write_new(out/'DISCOVERY.json',{'protocol_id':PROTOCOL,'agreement':record(agreement),'cases':[discovery[c] for c in ids],
        'original_policy':shell.POLICY,'envelope_policy':POLICY,'new_molecular_calls':0})
    groups={};triples=[]
    def group_for(protein,members,tid):
        keys=set().union(*({context.fragment_id(x) for x in discovery[c]['fragments']} for c in members))
        sid=protein+'__envelope_'+cache_key(sorted(keys))[:16]
        if sid not in groups:groups[sid]={'selection_id':sid,'protein_id':protein,'root_case_id':protein,
            'fragments':[context.fragment_record(k) for k in sorted(keys)],'members':[],'triple_ids':[],'requested_case_ids':[]}
        g=groups[sid];g['members']=sorted(set(g['members'])|set(members))
        if tid is not None:g['triple_ids'].append(tid)
        return sid
    for t in prior['triples']:
        r={k:copy.deepcopy(t[k]) for k in ('triple_id','protein_id','members','missing_members','source_statuses')}
        newmissing=[c for c in t['members'] if c in discovery and discovery[c]['status']!='discovered']
        r.update(envelope_discovery_missing_members=newmissing,selection_id=None,
            status='prior_preparation_unavailable' if t['missing_members'] else 'envelope_discovery_unavailable' if newmissing else 'selected')
        if r['status']=='selected':r['selection_id']=group_for(t['protein_id'],t['members'],t['triple_id'])
        triples.append(r)
    chosen={x['protein_id']:x for x in canonical['canonical_sources']};bytriple={t['triple_id']:t for t in triples};requests={};pilot=[]
    def request(cid,sid,role,tid=None):
        if sid is None:return None
        key=(sid,cid)
        if key not in requests:requests[key]={'case_id':cid,'selection_id':sid,'roles':[],'triple_ids':[]}
        r=requests[key]
        if role not in r['roles']:r['roles'].append(role)
        if tid is not None and tid not in r['triple_ids']:r['triple_ids'].append(tid)
        if cid not in groups[sid]['requested_case_ids']:groups[sid]['requested_case_ids'].append(cid)
        return key
    for t in triples:
        if t['status']=='selected':
            for cid in t['members']:request(cid,t['selection_id'],'primary_triple_member',t['triple_id'])
    for x in canonical['canonical_sources']:
        t=bytriple[x['selected_triple']['triple_id']]
        if t['members']!=x['selected_triple']['members']:raise InvalidArtifact('predeclared canonical triple changed')
        request(x['case_id'],t['selection_id'],'canonical_calibration',t['triple_id'])
        pilot.append({'case_id':x['case_id'],'selection_id':t['selection_id'],'role':'canonical_calibration','expected_class':x['expected_class'],
            'declared_triple_id':t['triple_id'],'status':t['status']})
    for cid in ('1H4I','4MAE','1KB0'):
        sid=group_for(cid,[cid],None) if discovery[cid]['status']=='discovered' else None
        request(cid,sid,'crystal_transfer');pilot.append({'case_id':cid,'selection_id':sid,'role':'crystal_transfer',
            'expected_class':src[cid]['expected_class'],'status':'selected' if sid else 'envelope_discovery_unavailable'})
    for accession,metal,sample in PROBES:
        gid=accession+'-pqq-la_model';cid=gid+f'__conditioned_{metal}__seed-1_sample-{sample}'
        x=chosen[gid];t=bytriple[x['selected_triple']['triple_id']]
        if src[cid]['status']!='prepared':raise InvalidArtifact('predeclared probe source not previously prepared')
        request(cid,t['selection_id'],'separate_probe',t['triple_id']);pilot.append({'case_id':cid,'selection_id':t['selection_id'],
            'role':'separate_probe','expected_class':src[cid]['expected_class'],'declared_triple_id':t['triple_id'],'status':t['status']})
    if len(pilot)!=34:raise InvalidArtifact('canonical25+crystal3+probe6 denominator changed')
    for g in groups.values():
        g['requested_case_ids'].sort();p=out/'selections'/g['selection_id']/'UNION.json'
        write_new(p,{**g,'protocol_id':PROTOCOL,'agreement':record(agreement),'rule':'fixed_complete_fragments_4p3_A_over_exact_declared_sources',
            'discovery':record(out/'DISCOVERY.json')});g['union']=record(p)
    selected_sources={cid:src[cid] for g in groups.values() for cid in g['requested_case_ids']}
    write_new(out/'SOURCE_INPUTS.json',{'source_preparation':prior['source_preparation'],'crystal_preparation':oldcr['source_preparation'],
        'sources':selected_sources,'config':config})
    result={'protocol_id':PROTOCOL,'prior_selection':record(prior_selection),'canonical_selection':record(canonical_selection),
        'agreement':record(agreement),'discovery':record(out/'DISCOVERY.json'),'source_inputs':record(out/'SOURCE_INPUTS.json'),
        'config':config,'original_policy':shell.POLICY,'envelope_policy':POLICY,'groups':list(groups.values()),'triples':triples,
        'requests':list(requests.values()),'pilot':pilot,'counts':{'triples':100,'prior_unavailable_triples':6,
            'selected_triples':sum(t['status']=='selected' for t in triples),'unique_groups':len(groups),'unique_source_selection_requests':len(requests),
            'primary_requests':sum('primary_triple_member' in r['roles'] for r in requests.values()),'pilot_denominator':34},
        'new_molecular_calls':0,'new_reference':None,'production_changed':False,'wall_seconds':time.monotonic()-start,
        'implementation':{n:record(ROOT/'scripts'/n) for n in ('motion_envelope.py','second_shell_context.py','consistent_context.py',
            'environment_context_chemistry.py','union_triple_transfer.py','coordination_preparation_context.py','mace_site_kinematics.py')}}
    write_new(out/'SELECTION.json',result);return {'selection':record(out/'SELECTION.json'),**result['counts']}


def validate(path):
    m=read_json(path)
    if m['protocol_id']!=PROTOCOL or m['envelope_policy']!=POLICY or m['original_policy']!=shell.POLICY:
        raise InvalidArtifact('frozen envelope protocol/policy differs')
    for pin in m['implementation'].values():verify(pin)
    for k in ('prior_selection','canonical_selection','agreement','discovery','source_inputs'):verify(m[k])
    for g in m['groups']:
        p=read_json(verify(g['union']))
        if any(p[k]!=g[k] for k in ('fragments','members','requested_case_ids','triple_ids')):raise InvalidArtifact('frozen union membership changed')
    return m


def prepare_task(row,group,config,out):
    r=context.prepare_one(row,group,config,out);r['preparation_profile']=PROTOCOL
    if r['status']=='prepared_awaiting_group_check':
        try:
            maps=mapped(r,config);pins={}
            for z,geometry in maps.items():
                p=Path(out)/'cases'/r['case_id']/(z+'_mapping.json');write_new(p,geometry);pins[z]=record(p)
            r['maps']=pins;r['paired_maps_equal']=maps['Ca']==maps['La'];r['mapping_eligible']=True
        except Exception as exc:r.update(status='mapping_unavailable',mapping_eligible=False,reason=str(exc))
    return r


def finalize(m,results):
    """Strict state anchors and complete declared triple membership, no substitution."""
    groups={g['selection_id']:g for g in m['groups']}
    for (sid,cid),r in results.items():
        anchor=results[sid,groups[sid]['members'][0]];r['state_anchor_source']=groups[sid]['members'][0]
        if r['status']=='prepared_awaiting_group_check':
            if not anchor.get('mapping_eligible') or anchor['status'] not in ('prepared_awaiting_group_check','prepared'):
                r.update(status='state_anchor_unavailable',reason='declared state anchor failed preparation or physical mapping')
            elif anchor.get('state_key')!=r['state_key'] or anchor.get('protein_key')!=r['protein_key']:
                r.update(status='source_state_mismatch',reason='selected identities/H parents/caps/bonds/charges/protein differ across source group')
            else:r['status']='prepared'
    triples=[]
    for t in m['triples']:
        r=copy.deepcopy(t)
        if r['status']=='selected':
            r['new_preparation_missing_members']=[c for c in r['members'] if results[r['selection_id'],c]['status']!='prepared']
            r['status']='prepared' if not r['new_preparation_missing_members'] else 'envelope_preparation_unavailable'
        triples.append(r)
    pilot=[]
    final_triples={t['triple_id']:t for t in triples}
    for p in m['pilot']:
        r=copy.deepcopy(p);prepared=results.get((p['selection_id'],p['case_id']))
        r['target_preparation_status']=prepared['status'] if prepared else p['status']
        r['status']=r['target_preparation_status'];r['reason']=prepared.get('reason') if prepared else 'declared union discovery unavailable'
        if p.get('declared_triple_id'):
            triple=final_triples[p['declared_triple_id']];r['declared_triple_preparation_status']=triple['status']
            if triple['status']!='prepared':
                r.update(status='declared_triple_unavailable',reason='one or more exact defining source members is unavailable')
        pilot.append(r)
    ordered=[results[r['selection_id'],r['case_id']] for r in m['requests']]
    return ordered,triples,pilot


def prepare(manifest,workers=4):
    start=time.monotonic();m=validate(manifest);out=Path(manifest).parent;src=read_json(verify(m['source_inputs']))['sources']
    groups={g['selection_id']:g for g in m['groups']};results={}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fs={ex.submit(prepare_task,src[r['case_id']],groups[r['selection_id']],m['config'],out/'prepared'/r['selection_id']):r for r in m['requests']}
        for f in as_completed(fs):
            request=fs[f];r=f.result();r.update(selection_id=request['selection_id'],roles=request['roles'],triple_ids=request['triple_ids'])
            results[r['selection_id'],r['case_id']]=r;print(json.dumps({'prepared':r['case_id'],'selection':r['selection_id'],'status':r['status'],'reason':r.get('reason')}),flush=True)
    ordered,triples,pilot=finalize(m,results)
    result={'protocol_id':PROTOCOL,'selection':record(manifest),'cases':ordered,'triples':triples,'pilot':pilot,
        'counts':{'triple_denominator':100,'prepared_triples':sum(t['status']=='prepared' for t in triples),
            'source_selection_pairs':len(ordered),'prepared_pairs':sum(r['status']=='prepared' for r in ordered),'pilot_denominator':34,
            'prepared_pilot':sum(r['status']=='prepared' for r in pilot),'prior_unavailable_triples':6},
        'new_molecular_calls':0,'new_reference':None,'new_protonation':False,'new_waters':False,'production_changed':False,
        'wall_seconds':time.monotonic()-start,'workers':workers}
    write_new(out/'PREPARATION.json',result);return {k:v for k,v in result.items() if k not in ('cases','triples','pilot')}


def reconcile(preparation,prior_implementation,output):
    """Apply strict failure accounting to completed real geometry without rerunning it."""
    p=read_json(preparation);m=read_json(verify(p['selection']))
    if record(prior_implementation)['sha256']!=m['implementation']['motion_envelope.py']['sha256']:
        raise InvalidArtifact('original implementation snapshot differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    m['implementation']['motion_envelope.py']=record(__file__)
    m['prior_geometry_implementation']=record(prior_implementation);m['prior_preparation']=record(preparation)
    write_new(out/'SELECTION.json',m);validate(out/'SELECTION.json')
    rows=copy.deepcopy(p['cases'])
    for r in rows:
        if r['status']=='prepared':r['status']='prepared_awaiting_group_check'
        for pin in r.get('maps',{}).values():verify(pin)
        for k in ('state_signature','protein_identity'):
            if k in r:verify(r[k])
    ordered,triples,pilot=finalize(m,{(r['selection_id'],r['case_id']):r for r in rows})
    p.update(selection=record(out/'SELECTION.json'),cases=ordered,triples=triples,pilot=pilot,
        prior_preparation=record(preparation),geometry_rerun=False,reconciliation='strict_target_and_defining_triple_availability')
    p['counts'].update(prepared_triples=sum(t['status']=='prepared' for t in triples),
        prepared_pairs=sum(r['status']=='prepared' for r in ordered),prepared_pilot=sum(r['status']=='prepared' for r in pilot))
    write_new(out/'PREPARATION.json',p);return {'preparation':record(out/'PREPARATION.json'),**p['counts']}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    q=s.add_parser('select')
    for k in ('prior_selection','canonical_selection','agreement','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    q.add_argument('--workers',type=int,default=4)
    q=s.add_parser('prepare');q.add_argument('--manifest',required=True);q.add_argument('--workers',type=int,default=4)
    q=s.add_parser('reconcile')
    for k in ('preparation','prior_implementation','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());cmd=a.pop('command');print(json.dumps({'select':select,'prepare':prepare,'reconcile':reconcile}[cmd](**a),indent=2))
if __name__=='__main__':main()
