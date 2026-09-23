"""All declared La triples: source-graph preparation and exact archive mapping only."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools
import json
from pathlib import Path
import time

from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
import consistent_context as context
from accommodation_folds_compare import checked_scoped_transfer, source_groups
from second_shell_context import parent_state
import union_adaptive as union

PROTOCOL='Nikasha_three_La_source_fragment_union_preparation_v1'
ROOT=Path(__file__).resolve().parents[1]
DEFAULT_PLAN=ROOT/'diagnostics/union_triple_preparation_20260923/PLAN.md'
STRESS='a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1'


def fragment_set(row):
    return {context.fragment_id(x) for x in read_json(verify(row['representations']['context']['preparation']))['added_fragments']}


def selection(inventory, output, agreement=DEFAULT_PLAN):
    start=time.monotonic();d=read_json(inventory);prep=read_json(verify(d['source_preparation']))
    ten=read_json(verify(d['union_preparation']));source=read_json(verify(prep['source_manifest']))
    source_groups(source);original={r['case_id']:r for r in prep['cases']}
    ten_groups={g['root_case_id']:g for g in ten['groups']};by_group=defaultdict(list)
    for c in source['cases']:
        if c['source_conditioning_metal']=='La' and not c['canonical_coordinate_match']:
            by_group[c['root_case_id']].append(c['case_id'])
    expected={(gid,tuple(sorted(triple))) for gid,ids in by_group.items() for triple in itertools.combinations(sorted(ids),3)}
    actual={(t['group'],tuple(sorted(t['sources']))) for t in d['triples']}
    if len(d['triples'])!=100 or len(expected)!=100 or actual!=expected:
        raise InvalidArtifact('all100 declared La triples must remain, without substitutes')
    out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces/')
    out.mkdir(parents=True,exist_ok=False);groups={};triples=[]
    for n,t in enumerate(d['triples']):
        ids=t['sources'];missing=[c for c in ids if original[c]['status']!='prepared']
        if missing!=t['missing_sources']:raise InvalidArtifact('old preparation exclusions changed')
        row={'triple_id':f'triple_{n:03d}','protein_id':t['group'],'members':ids,'missing_members':missing,
            'status':'prior_preparation_unavailable' if missing else 'selected','selection_id':None,
            'source_statuses':{c:original[c]['status'] for c in ids}}
        if not missing:
            fragments=set().union(*(fragment_set(original[c]) for c in ids))
            old={context.fragment_id(f) for f in ten_groups[t['group']]['fragments']}
            if not fragments<=old:raise InvalidArtifact('triple selection is not a subset of original union')
            status='identical_fragment_selection' if fragments==old else 'smaller_fragment_selection'
            if status!=t['status']:raise InvalidArtifact('read-only selection inventory no longer matches')
            sid=t['group']+'__'+cache_key(sorted(fragments))[:16];row['selection_id']=sid
            if sid not in groups:
                groups[sid]={'selection_id':sid,'protein_id':t['group'],'root_case_id':t['group'],
                    'fragments':[context.fragment_record(k) for k in sorted(fragments)],
                    'members':[],'selection_status':status,'triple_ids':[],'stress_members':[]}
            g=groups[sid];g['members']=sorted(set(g['members'])|set(ids));g['triple_ids'].append(row['triple_id'])
        triples.append(row)
    stress_group=[g for g in groups.values() if g['protein_id']=='a0a3f2yly8-pqq-la_model']
    if len(stress_group)!=1 or original[STRESS]['source']['source_conditioning_metal']!='Ca':
        raise InvalidArtifact('single predeclared A0A3 La-selection stress probe unavailable')
    stress_group[0]['stress_members']=[STRESS]
    for g in groups.values():
        path=out/'selections'/g['selection_id']/'UNION.json'
        write_new(path,{**g,'membership_contributors_are_La_only':True,'agreement':record(agreement),
            'rule':'union_of_each_declared_three_La_sources_before_any_energy'})
        g['union']=record(path)
    manifest={'protocol_id':PROTOCOL,'inventory':record(inventory),'source_preparation':d['source_preparation'],
        'tenfold_preparation':d['union_preparation'],'source_manifest':prep['source_manifest'],
        'config':prep['config'],'agreement':record(agreement),'groups':list(groups.values()),'triples':triples,
        'triple_denominator':100,'complete_triples':sum(not t['missing_members'] for t in triples),
        'unique_selections':len(groups),'unique_source_selection_pairs':sum(len(g['members']) for g in groups.values()),
        'stress_probe':{'source':STRESS,'selection_id':stress_group[0]['selection_id'],'separate_from100':True},
        'new_molecular_calls':0,'production_changed':False,'selection_wall_seconds':time.monotonic()-start,
        'implementation':{n:record(ROOT/'scripts'/n) for n in ('union_triple_preparation.py','consistent_context.py','second_shell_context.py','environment_context_chemistry.py','union_adaptive.py','compact_solvation.py','accommodation_folds_compare.py')}}
    write_new(out/'SELECTION.json',manifest);return {'selection':record(out/'SELECTION.json'),**{k:manifest[k] for k in ('triple_denominator','complete_triples','unique_selections','unique_source_selection_pairs','new_molecular_calls')}}


def validate_selection(path):
    m=read_json(path)
    if m['protocol_id']!=PROTOCOL or m['triple_denominator']!=100 or len(m['triples'])!=100 or m['new_molecular_calls']:
        raise InvalidArtifact('preparation-only scope changed')
    for p in m['implementation'].values():verify(p)
    for k in ('inventory','source_preparation','tenfold_preparation','source_manifest','agreement'):verify(m[k])
    for g in m['groups']:
        pinned=read_json(verify(g['union']))
        if any(pinned[k]!=g[k] for k in ('fragments','members','stress_members','triple_ids')):
            raise InvalidArtifact('frozen triple selection changed')
    return m


def prepare(manifest,workers=4):
    start=time.monotonic();m=validate_selection(manifest);out=Path(manifest).parent
    source={r['case_id']:r for r in read_json(verify(m['source_preparation']))['cases']}
    ten={r['case_id']:r for r in read_json(verify(m['tenfold_preparation']))['cases']}
    results={};pending=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for g in m['groups']:
            for cid in g['members']+g['stress_members']:
                pending.append((pool.submit(context.prepare_one,source[cid],g,m['config'],out/'prepared'/g['selection_id']),g,cid))
        lookup={f:(g,cid) for f,g,cid in pending}
        for f in as_completed(lookup):
            g,cid=lookup[f];r=f.result();r['selection_id']=g['selection_id'];r['protein_id']=g['protein_id']
            r['is_separate_stress_probe']=cid in g['stress_members'];results[g['selection_id'],cid]=r
            print(json.dumps({'source':cid,'selection':g['selection_id'],'status':r['status']}),flush=True)
    for g in m['groups']:
        anchor=results[g['selection_id'],g['members'][0]]
        for cid in g['members']+g['stress_members']:
            r=results[g['selection_id'],cid]
            if r['status']!='prepared_awaiting_group_check':continue
            r['state_anchor_source']=g['members'][0]
            if 'state_key' not in anchor or (r['state_key'],r['protein_key'])!=(anchor['state_key'],anchor['protein_key']):
                r.update(status='source_state_mismatch',reason='actual selected source atom/H-parent/cap/charge/protein identity differs')
                continue
            r['status']='prepared';old=ten[cid]
            same_state=r['state_key']==old.get('state_key')
            same_geometry=all(context.reusable_state(r['representations']['context']['endpoints'][z],old['representations']['context']['endpoints'][z]) for z in ('Ca','La')) if old['status']=='prepared' else False
            r['tenfold_comparison']={'same_state':same_state,'same_geometry_and_charge':same_geometry,
                'exact_reuse_eligible':same_state and same_geometry,'tenfold_state_signature':old.get('state_signature'),
                'atom_count':old.get('atom_count'),'cap_count':old.get('cap_count'),'La_charge':old.get('new_La_charge'),
                'atom_count_difference':r['atom_count']-old['atom_count'],'cap_count_difference':r['cap_count']-old['cap_count']}
    triples=[]
    for t in m['triples']:
        row=dict(t)
        if not row['missing_members']:
            available=[results[t['selection_id'],c] for c in t['members']]
            row['preparation_missing_members']=[r['case_id'] for r in available if r['status']!='prepared']
            row['status']='prepared' if not row['preparation_missing_members'] else 'preparation_unavailable'
            row['all_members_exact_tenfold']=all(r.get('tenfold_comparison',{}).get('exact_reuse_eligible',False) for r in available)
        triples.append(row)
    ordered=[results[g['selection_id'],c] for g in m['groups'] for c in g['members']+g['stress_members']]
    result={'protocol_id':PROTOCOL,'selection':record(manifest),'triples':triples,'cases':ordered,
        'triple_denominator':100,'prepared_triples':sum(t['status']=='prepared' for t in triples),
        'unique_source_selection_pairs':m['unique_source_selection_pairs'],'separate_stress_probe_count':1,
        'supported_pairs_including_stress':sum(r['status']=='prepared' for r in ordered),
        'new_molecular_calls':0,'new_protonation':False,'wall_seconds':time.monotonic()-start,'workers':workers}
    write_new(out/'PREPARATION.json',result);return {k:v for k,v in result.items() if k not in ('cases','triples')}


def reuse(preparation,tenfold_collection,local_inventory,local_comparison,output):
    """Audit actual origins. No adaptive map/geometry or score is substituted."""
    start=time.monotonic();p=read_json(preparation);m=validate_selection(verify(p['selection']))
    ten=read_json(tenfold_collection);union_rows={r['case_id']:r for r in ten['rows']}
    local=read_json(local_inventory);local_rows={r['case_id']:r for r in local['cases']}
    source={r['case_id']:r for r in read_json(verify(m['source_preparation']))['cases']}
    comparison=read_json(local_comparison)
    if local['source_preparation']!=m['source_preparation'] or comparison['inventory']!=record(local_inventory):
        raise InvalidArtifact('local source archive differs')
    low={}
    for pin in comparison['collections']:
        c=read_json(verify(pin))
        for row in c['rows']:
            if row['representation']=='context' and row['solver']=='native':low[row['case_id'],row['metal']]=(row,c['manifest'],pin)
    native_cache={};rows=[]
    def audit_native(ep):
        key=ep['native_MACE_receipt']['sha256']
        if key not in native_cache:
            r,_=union.origin(ep,m['config']['model']);native_cache[key]={'energy_eV':r['energy_eV'],
                'receipt':ep['native_MACE_receipt'],'forces':r['forces'],'force_units':'eV/angstrom','force_definition':r['force_definition'],
                'model':m['config']['model'],'source_manifest':ep['source_manifest'],'source_task_id':ep['source_task_id']}
        return native_cache[key]
    for r in p['cases']:
        cid=r['case_id'];row={'case_id':cid,'selection_id':r['selection_id'],'is_separate_stress_probe':r['is_separate_stress_probe'],
            'preparation_status':r['status'],'native_endpoints':{},'solvent_endpoints':{},'reuse_source':None,
            'candidate_pool_reuse':'not_audited_requires_exact_physical_maps_and_settings','classification':None}
        if r['status']=='prepared':
            prepared=r['representations']['context']['endpoints'];old=union_rows.get(cid)
            if r['tenfold_comparison']['exact_reuse_eligible'] and old and old['status']=='complete':
                if old['group_state_signature']['sha256']!=r['tenfold_comparison']['tenfold_state_signature']['sha256']:
                    raise InvalidArtifact('tenfold collection graph state differs')
                for z in ('Ca','La'):
                    ep=old['native_endpoints'][z]
                    if not context.reusable_state(prepared[z],ep):raise InvalidArtifact('tenfold origin coordinate/state differs')
                    row['native_endpoints'][z]=audit_native(ep)
                    for medium in ('vacuum','alpb'):
                        cell=old['solvent_endpoints'][z][medium]
                        from compact_solvation import completed,input_text,METHOD,PROTOCOL as SOLVENT_PROTOCOL
                        actual=completed(verify(cell['manifest']),cell['task_id'])
                        if actual is None or any(actual[k]!=cell[k] for k in actual):raise InvalidArtifact('tenfold actual solvent receipt differs')
                        mm=read_json(verify(cell['manifest']));task=next(t for t in mm.get('all_tasks',mm['tasks']) if t['task_id']==cell['task_id'])
                        if not context.reusable_state(prepared[z],task) or task['medium']!=medium or task['solver']!='native':
                            raise InvalidArtifact('tenfold solvent physical state differs')
                        if (mm['method_id']!=METHOD or mm['protocol_id']!=SOLVENT_PROTOCOL or
                                verify(task['input']).read_text().replace(' MaxIter 500\n','')!=input_text(task['charge'],task['multiplicity'],medium,'native')):
                            raise InvalidArtifact('tenfold solvent Hamiltonian or numerical recipe differs')
                        row['solvent_endpoints'].setdefault(z,{})[medium]=actual
                row['reuse_source']='tenfold_exact_origin'
            else:
                oldlocal=source[cid];state=parent_state(oldlocal['core'],m['config']['topology'],require_endpoint_receipts=False)
                audit=read_json(verify(oldlocal['representations']['context']['preparation']))
                # Historical local audits predate the explicit formal-charge ledger.
                # Replay their exact source-backed fragment list through the same
                # graph implementation to derive the complete physical signature.
                reconstructed,audit=context.forced_expansion(state,audit['added_fragments'])
                if not all(context.same_geometry(reconstructed[z],xyz(verify(oldlocal['representations']['context']['endpoints'][z]['xyz']))) for z in ('Ca','La')):
                    raise InvalidArtifact('local context graph reconstruction changed original coordinates')
                signature=context.state_signature(state,audit,{z:oldlocal['representations']['context']['endpoints'][z]['charge'] for z in ('Ca','La')})
                endpoints=local_rows.get(cid,{}).get('representations',{}).get('context',{}).get('endpoints',{})
                same=cache_key(signature)==r['state_key'] and len(endpoints)==2 and all(context.reusable_state(prepared[z],endpoints[z]) for z in ('Ca','La'))
                if same:
                    for z in ('Ca','La'):
                        row['native_endpoints'][z]=audit_native(endpoints[z])
                        oldlow=low.get((cid,z))
                        if oldlow:
                            lr,lm,lp=oldlow;checked_scoped_transfer(lr,endpoints[z],record(local_inventory),lm)
                            row['solvent_endpoints'][z]={k:lr['endpoints'][k] for k in ('vacuum','alpb') if lr['endpoints'][k]['status']=='complete'}
                    row['reuse_source']='original_local_context_exact_origin'
                else:row['reuse_source']='no_exact_archived_origin'
        row['missing_native_endpoints']=[z for z in ('Ca','La') if z not in row['native_endpoints']]
        row['missing_solvent_cells']=[z+'_'+s for z in ('Ca','La') for s in ('vacuum','alpb') if s not in row['solvent_endpoints'].get(z,{})]
        rows.append(row)
    primary=[r for r in rows if not r['is_separate_stress_probe']]
    summary={'source_selection_pairs':len(primary),'reuse_sources':dict(Counter(r['reuse_source'] for r in primary)),
        'missing_native_endpoints':sum(len(r['missing_native_endpoints']) for r in primary),
        'missing_solvent_cells':sum(len(r['missing_solvent_cells']) for r in primary)}
    result={'protocol_id':PROTOCOL,'preparation':record(preparation),'tenfold_collection':record(tenfold_collection),
        'local_inventory':record(local_inventory),'local_comparison':record(local_comparison),'rows':rows,'summary':summary,
        'new_molecular_calls':0,'new_reference':None,'production_changed':False,'wall_seconds':time.monotonic()-start,
        'limits':'Only original scalar/force receipts mapped; changed representation has no new calibrated decision or inherited adaptive pool.'}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    q=s.add_parser('select');q.add_argument('--inventory',required=True);q.add_argument('--output',required=True);q.add_argument('--agreement',default=str(DEFAULT_PLAN))
    q=s.add_parser('prepare');q.add_argument('--manifest',required=True);q.add_argument('--workers',type=int,default=4)
    q=s.add_parser('reuse');q.add_argument('--preparation',required=True);q.add_argument('--tenfold-collection',required=True);q.add_argument('--local-inventory',required=True);q.add_argument('--local-comparison',required=True);q.add_argument('--output',required=True)
    a=vars(p.parse_args());cmd=a.pop('command');print(json.dumps({'select':selection,'prepare':prepare,'reuse':reuse}[cmd](**a),indent=2))


if __name__=='__main__':main()
