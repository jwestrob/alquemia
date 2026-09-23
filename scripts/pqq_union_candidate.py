"""Opt-in declared ten-fold PQQ union preparation and exact scientific replay.

This adapter performs no molecular energy/force calls and changes no defaults.
"""
from __future__ import annotations
import argparse
import copy
import itertools
import json
from pathlib import Path
import re
import statistics
import time

import gemmi
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
import accommodation_folds as folds
from accommodation_fold_proposals import strict_summary
from accommodation_folds_compare import decision
import consistent_context as context
import pqq_adaptive_candidate as candidate
import pqq_fast_prepare as source_prepare
from pqq_ensemble import resolve_residue
import pqq_standard as standard
import union_adaptive as union
import nikasha_pool as pool

PROTOCOL = 'Nikasha_opt_in_declared_ten_fold_union_adaptive_replay_v1'
PROFILE = 'frozen_SLSQP200_native_GFN2_8rank_MaxIter500_v1'
CANDIDATES = ['origin', 'adaptive_Ca', 'adaptive_La']
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PLAN = ROOT/'diagnostics/pqq_union_candidate_20260923/PLAN.md'


def archive(preparation, crystals, comparison, reference, output,
            release=standard.DEFAULT_RELEASE):
    """Make explicit, portable pins to existing compatible experiment records."""
    r=standard.release(release); p=read_json(preparation); cr=read_json(crystals)
    comp=read_json(comparison); ref=read_json(reference)
    if (p['protocol_id']!=context.PROTOCOL or cr['protocol_id']!=context.PROTOCOL or
            comp['protocol_id']!=union.POOL_PROTOCOL or comp['reference']!=record(reference) or
            ref['optimizer_settings']!=union.SETTINGS or ref['settings']!=union.POOL_SETTINGS or
            not ref['complete']):raise InvalidArtifact('incompatible completed union archive')
    selection=read_json(verify(comp['selection']))
    pins=comp['collections']+ref['collections']+[selection['pilot_collection']]
    unique={x['sha256']:x for x in pins}
    result={'protocol_id':PROTOCOL,'profile':PROFILE,'union_preparations':[record(preparation),record(crystals)],
        'source_preparations':[p['source_preparation'],cr['source_preparation']],
        'comparison':record(comparison),'reference':record(reference),'release':record(release),
        'original_calibration':r['artifacts']['calibration'],'pool_collections':list(unique.values()),
        'config':p['config'],'new_molecular_calls':0}
    if p['config']!=cr['config'] or p['config']!=r['source_configuration']:
        raise InvalidArtifact('source configuration differs')
    candidate.put(output,result);return {'archive':record(output),'new_molecular_calls':0}


def source_identity(case, config):
    st=gemmi.read_structure(str(verify(case['source_structure'])))
    if case['normalization']=='protenix_generic_PQQ':
        if len(st)!=1 or case['source_conditioning_metal'] not in ('Ca','La'):
            raise InvalidArtifact('one-model declared Ca/La fold required')
        raw=case['raw_source_metal'];res=resolve_residue(st[0],raw)
        if (raw['element']!=case['source_conditioning_metal'] or len(res)!=1 or
                res[0].name!=raw['atom'] or res[0].element.name!=raw['element']):
            raise InvalidArtifact('actual conditioning metal differs')
        if any(raw.get(k,'')!=case['metal'].get(k,'') for k in ('chain','resnum','icode','atom')):
            raise InvalidArtifact('native/normalized atom selectors disagree')
        model=st[0];chains=list(model)
    elif case['normalization']=='selected_crystal_chain':
        _,_,crystal,_,_=source_prepare.helpers(config);a=case['assembly']
        model=crystal.exact_model(st,a['model']);chains=[crystal.exact_chain(model,a['chain'])]
        crystal.parse_atom_selector(a['raw_metal']);crystal.parse_residue_selector(a['raw_pqq'])
    else:raise InvalidArtifact('unsupported source normalization')
    sequence=[]
    for chain in chains:
        aa=[(r.seqid.num,r.seqid.icode.strip(),r.name) for r in chain
            if gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        if aa:sequence.append((chain.name,aa))
    if not sequence or len({c for c,_ in sequence})!=len(sequence):
        raise InvalidArtifact('missing/ambiguous actual protein chains')
    for role in case['roles'].values():resolve_residue(model,role)
    return {'sequence_and_numbering':sequence,'site_roles':case['roles'],
            'pqq_selector':case['pqq'],'assembly':case['assembly']}


def check_request(req):
    if req.get('protocol_id')!=PROTOCOL or req.get('profile')!=PROFILE:
        raise InvalidArtifact('unsupported profile; cheaper numerical profiles are not enabled')
    release=standard.release(verify(req['release']))
    if req['config']!=release['source_configuration']:raise InvalidArtifact('source method/configuration differs')
    ref=read_json(verify(req['reference']))
    arc=read_json(verify(req['archive']))
    if (arc['protocol_id']!=PROTOCOL or arc['profile']!=PROFILE or arc['reference']!=req['reference']
            or arc['release']!=req['release'] or arc['config']!=req['config'] or ref['model']!=req['config']['model']):
        raise InvalidArtifact('source/archive/reference method identity differs')
    if ref['optimizer_settings']!=union.SETTINGS or ref['settings']!=union.POOL_SETTINGS or not ref['complete']:
        raise InvalidArtifact('frozen numerical/reference compatibility failed')
    historical_sources={x['case_id']:x['source'] for pin in arc['source_preparations']
        for x in read_json(verify(pin))['cases']}
    cases={x['case_id']:x for x in req['cases']}
    if len(cases)!=len(req['cases']) or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',c) for c in cases):
        raise InvalidArtifact('safe unique source IDs required')
    seen=set();group_ids=set();identities={}
    for g in req['groups']:
        ids=g['members'];rows=[cases[x] for x in ids]
        if (g['protein_id'] in group_ids or len(set(ids))!=len(ids) or seen.intersection(ids)
                or any(x.get('root_case_id',x['case_id'])!=g['protein_id'] for x in rows)):
            raise InvalidArtifact('source/protein group overlap or mismatch')
        if g['kind']=='ten_folds':
            if (len(ids)!=10 or any(sum(x['source_conditioning_metal']==z for x in rows)!=5 for z in ('Ca','La'))
                    or len({x['source_structure']['sha256'] for x in rows})!=10):
                raise InvalidArtifact('exact five distinct Ca plus five La folds required; three-fold union is not supported')
            canon=[x for x in rows if x.get('canonical_coordinate_match')]
            if len(canon)!=1 or canon[0]['case_id']!=g['canonical_case_id'] or canon[0]['source_conditioning_metal']!='La':
                raise InvalidArtifact('one predeclared canonical La member required')
        elif g['kind']=='crystal_singleton':
            if len(ids)!=1 or rows[0]['normalization']!='selected_crystal_chain':
                raise InvalidArtifact('explicit selected crystal singleton required')
        else:raise InvalidArtifact('undeclared source group policy')
        actual=[source_identity(x,req['config']) for x in rows]
        if any(x!=actual[0] for x in actual[1:]):
            raise InvalidArtifact('actual sequence, numbering, roles or assembly differs within protein')
        if g['kind']=='ten_folds':
            canonical=canon[0];historical=historical_sources.get(canonical['case_id'],{})
            if (not historical.get('canonical_coordinate_match') or
                    not matching_source(canonical,historical) or
                    historical.get('root_case_id')!=g['protein_id']):
                raise InvalidArtifact('canonical membership must match an actual archived calibration source; unknown state anchors are unsupported')
        identities[g['protein_id']]=cache_key(actual[0]);seen.update(ids);group_ids.add(g['protein_id'])
    if not seen or seen!=set(cases):raise InvalidArtifact('every source must belong to exactly one declared group')
    return identities


def request(sources, groups, output, archive_path, crystal_sources=None, crystals=()):
    source=read_json(sources);arc=read_json(archive_path);rows=[];decl=[]
    if source['config']!=arc['config']:raise InvalidArtifact('source inventory configuration differs')
    for gid in groups:
        selected=[copy.deepcopy(x) for x in source['cases'] if x.get('root_case_id')==gid]
        canon=[x['case_id'] for x in selected if x.get('canonical_coordinate_match')]
        decl.append({'protein_id':gid,'kind':'ten_folds','members':[x['case_id'] for x in selected],
            'canonical_case_id':canon[0] if len(canon)==1 else None,
            'selection_rule':'all five Ca and five La members declared in supplied source inventory'})
        rows.extend(selected)
    if crystals:
        cs=read_json(crystal_sources)
        if cs['config']!=source['config']:raise InvalidArtifact('crystal configuration differs')
        for cid in crystals:
            x=copy.deepcopy(next(x for x in cs['cases'] if x['case_id']==cid))
            x.update(root_case_id=cid,source_conditioning_metal=None,canonical_coordinate_match=False)
            rows.append(x);decl.append({'protein_id':cid,'kind':'crystal_singleton','members':[cid],
                'canonical_case_id':cid,'selection_rule':'explicit existing crystal selector'})
    req={'protocol_id':PROTOCOL,'profile':PROFILE,'config':source['config'],'cases':rows,'groups':decl,
        'source_inventory':record(sources),'crystal_inventory':record(crystal_sources) if crystal_sources else None,
        'archive':record(archive_path),'reference':arc['reference'],'release':arc['release'],
        'production_changed':False,'new_molecular_calls':0}
    req['source_identities']=check_request(req);candidate.put(output,req)
    return {'request':record(output),'groups':len(decl),'sources':len(rows),'new_molecular_calls':0}


def matching_source(expected, actual):
    if 'source_structure' not in actual:return False
    original=actual.get('original_source_structure',actual['source_structure'])
    return original==expected['source_structure'] and all(actual.get(k)==expected.get(k)
        for k in ('case_id','roles','metal','pqq','assembly','normalization'))


def prepare(request, output, source_mode='reuse-exact', agreement=DEFAULT_PLAN):
    start=time.monotonic();req=read_json(request);identities=check_request(req);arc=read_json(verify(req['archive']))
    if source_mode not in ('reuse-exact','fresh'):raise InvalidArtifact('unsupported source preparation mode')
    out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces/')
    out.mkdir(parents=True,exist_ok=False);source_rows=[];catalog={}
    for pin in arc['source_preparations']:
        old=read_json(verify(pin))
        if old['config']!=req['config']:raise InvalidArtifact('archive source configuration differs')
        for row in old['cases']:catalog.setdefault(row['case_id'],[]).append((row,pin))
    for src in req['cases']:
        matches=[(x,p) for x,p in catalog.get(src['case_id'],[]) if matching_source(src,x.get('source',{}))]
        if source_mode=='reuse-exact' and len(matches)==1:
            row=copy.deepcopy(matches[0][0]);row['source_preparation_reuse']=matches[0][1]
        elif source_mode=='reuse-exact':
            row={'case_id':src['case_id'],'source':src,'status':'unavailable','reason':'exact_source_preparation_not_available'}
        elif src['normalization']=='protenix_generic_PQQ':row=folds.prepare_one(src,req['config'],out/'source')
        else:
            try:row=source_prepare.prepare_case(src,req['config'],out/'source'/src['case_id'])
            except Exception as exc:row={'case_id':src['case_id'],'source':src,'status':'unsupported','reason':str(exc)}
        row.setdefault('source',src);source_rows.append(row)
    write_new(out/'SOURCE_PREPARATION.json',{'request':record(request),'config':req['config'],'cases':source_rows,
        'source_mode':source_mode,'new_molecular_calls':0})
    source_index={x['case_id']:x for x in source_rows};frozen=[]
    for g in req['groups']:
        selected=set()
        for cid in g['members']:
            row=source_index[cid]
            if row['status']=='prepared':
                audit=read_json(verify(row['representations']['context']['preparation']))
                selected.update(context.fragment_id(x) for x in audit['added_fragments'])
        fragments=[context.fragment_record(k) for k in sorted(selected)]
        path=out/'groups'/g['protein_id']/'UNION.json'
        write_new(path,{**g,'fragments':fragments,'source_request':record(request),
            'rule':'union_of_all_supported_declared_source_fragment_identities_before_energies'})
        frozen.append({**g,'union':record(path),'fragments':fragments})
    write_new(out/'GROUP_SELECTION.json',{'request':record(request),'groups':frozen})
    results=[]
    for g in frozen:
        local=[]
        for cid in g['members']:
            original=source_index[cid]
            if original['status']!='prepared':
                row={'case_id':cid,'source':original['source'],'status':'source_preparation_unavailable',
                    'reason':original.get('reason',original['status']),'original_preparation':original}
            else:row=context.prepare_one(original,g,req['config'],out/'union')
            local.append(row)
        canon=next(x for x in local if x['case_id']==g['canonical_case_id'])
        for row in local:
            if row['status']=='prepared_awaiting_group_check':
                if 'state_key' not in canon:row.update(status='canonical_state_unavailable',reason='declared canonical member unsupported')
                elif row['protein_key']!=canon['protein_key'] or row['state_key']!=canon['state_key']:
                    row.update(status='group_state_mismatch',reason='source atom/H-parent/cap/charge/protein identity differs')
                else:row.update(status='prepared',canonical_signature=canon['state_signature'])
            row['group_id']=g['protein_id'];results.append(row)
    pp=out/'PREPARATION.json';write_new(pp,{'request':record(request),'groups':frozen,'cases':results,
        'source_preparation':record(out/'SOURCE_PREPARATION.json'),'config':req['config'],
        'denominator':len(results),'supported':sum(x['status']=='prepared' for x in results),
        'new_molecular_calls':0,'wall_seconds':time.monotonic()-start})
    plan={'protocol_id':PROTOCOL,'profile':PROFILE,'request':record(request),'preparation':record(pp),
        'group_selection':record(out/'GROUP_SELECTION.json'),'agreement':record(agreement),
        'reference':req['reference'],'source_identities':identities,'settings':union.SETTINGS,
        'pool_settings':union.POOL_SETTINGS,'implementation':{n:record(Path(__file__).parent/n) for n in
            ('pqq_union_candidate.py','pqq_adaptive_candidate.py','consistent_context.py','union_adaptive.py',
             'pqq_fast_prepare.py','accommodation_folds.py','nikasha_pool.py')},
        'fresh_scientific_execution_status':'not_implemented_until_numerical_profile_qualification',
        'production_changed':False,'automatic_fallback':False,'new_molecular_calls':0}
    candidate.put(out/'plan.json',plan);return dry_run(out/'plan.json')


def checked(plan):
    p=read_json(plan);req=read_json(verify(p['request']));prep=read_json(verify(p['preparation']))
    if (p['protocol_id']!=PROTOCOL or p['profile']!=PROFILE or p['settings']!=union.SETTINGS
            or p['pool_settings']!=union.POOL_SETTINGS or p['production_changed'] or p['automatic_fallback']):
        raise InvalidArtifact('candidate protocol/profile changed')
    for pin in p['implementation'].values():verify(pin)
    verify(p['agreement']);selection=read_json(verify(p['group_selection']))
    if (check_request(req)!=p['source_identities'] or selection['groups']!=prep['groups'] or
            prep['request']!=p['request'] or p['reference']!=req['reference']):
        raise InvalidArtifact('source/group preparation differs')
    if [x['case_id'] for x in prep['cases']]!=[x['case_id'] for x in req['cases']]:
        raise InvalidArtifact('declared member denominator/order changed')
    for row in prep['cases']:
        if row['status']!='prepared':continue
        g=next(g for g in prep['groups'] if g['protein_id']==row['group_id'])
        if row['union']!=g['union'] or read_json(verify(g['union']))['fragments']!=g['fragments']:
            raise InvalidArtifact('fixed union membership changed')
        signature=read_json(verify(row['state_signature']));canonical=read_json(verify(row['canonical_signature']))
        if signature!=canonical or cache_key(signature)!=row['state_key']:raise InvalidArtifact('prepared group state changed')
        for ep in row['representations']['context']['endpoints'].values():verify(ep['xyz'])
    return p,req,prep


def dry_run(plan):
    p,req,prep=checked(plan);n=prep['supported']
    return {'plan':record(plan),'status':'prepared_for_archive_replay','profile':PROFILE,
        'groups':len(req['groups']),'denominator':prep['denominator'],'supported':n,
        'unavailable':[{'case_id':x['case_id'],'status':x['status'],'reason':x.get('reason')} for x in prep['cases'] if x['status']!='prepared'],
        'future_union_candidate_maximum_calls_before_reuse':{'origin_MACE':2*n,'optimizer_starts':2*n,
            'cross_MACE':2*n,'native_GFN2':12*n,'DFT':0},
        'fresh_original_local_static_requires_separate_components':True,
        'fresh_execution_status':p['fresh_scientific_execution_status'],'new_molecular_calls':0,
        'three_fold_union_compatible':False,'production_changed':False}


def replay(plan, output):
    p,req,prep=checked(plan);arc=read_json(verify(req['archive']));ref=read_json(verify(req['reference']))
    old=read_json(verify(arc['comparison']));original=read_json(verify(arc['original_calibration']))
    baseline={x['case_id']:x['methods']['context_composite'] for x in old['rows']}
    original_band=old['bands']['context_composite'];union_band=old['bands']['context_union']
    aliases={x['actual_union_case_id']:x['case_id'] for x in ref['rows']}
    canonical_original={x['case_id']:x for x in original['rows'] if x['representation']=='context'}
    archived_prepared={x['case_id']:x for pin in arc['union_preparations'] for x in read_json(verify(pin))['cases']}
    archived_local={x['case_id']:x for pin in arc['source_preparations'] for x in read_json(verify(pin))['cases']}
    prepared_local={x['case_id']:x for x in read_json(verify(prep['source_preparation']))['cases']}
    sources={x['case_id']:x for x in req['cases']};available={};manifest_cache={}
    for pin in arc['pool_collections']:
        c=read_json(verify(pin));m=read_json(verify(c['manifest']));sm=read_json(verify(m['source_manifest']))
        if (c['protocol_id']!=union.POOL_PROTOCOL or m['settings']!=ref['settings'] or
                m['model']!=ref['model'] or sm['settings']!=ref['optimizer_settings']):
            raise InvalidArtifact('archive scoring model/settings differ')
        manifest_cache[pin['sha256']]=(m,sm)
        for row in c['cases']:
            if row['case_id'] in available and available[row['case_id']][0]!=row:
                raise InvalidArtifact('conflicting duplicate archived case')
            available[row['case_id']]=(row,pin)
    rows=[]
    for prepared in prep['cases']:
        cid=prepared['case_id'];aid=aliases.get(cid,cid);src=sources[cid]
        local=prepared_local[cid];old_local=archived_local.get(cid)
        raw_match=bool(old_local and matching_source(src,old_local.get('source',{})))
        local_match=raw_match and local['status']=='prepared' and old_local['status']=='prepared'
        if local_match:
            local_match=all(context.reusable_state(local['representations']['context']['endpoints'][z],
                old_local['representations']['context']['endpoints'][z]) for z in ('Ca','La'))
        base=baseline.get(cid)
        if base is None and aid in canonical_original:
            value=canonical_original[aid]['composite_R_model_kcal_mol'];base={'R':value,'decision':decision(value,original_band)}
        row={'case_id':cid,'root_case_id':prepared['group_id'],'source_conditioning_metal':src.get('source_conditioning_metal'),
            'canonical_coordinate_match':src.get('canonical_coordinate_match',False),'source_structure':src['source_structure'],
            'preparation_status':prepared['status'],'source_domain_status':'consumed_reference' if raw_match else 'unvalidated_input_domain',
            'original_static':base if base and local_match else {'R':None,'decision':'unavailable'},
            'original_static_archive_compatible':local_match,'union_static':{'R':None,'decision':'unavailable'},
            'candidate':{v:{'R':None,'decision':'unavailable'} for v in ('operational','mathematical')},
            'status':'unavailable','reason':prepared.get('reason','compatible_archive_not_available')}
        if prepared['status']=='prepared' and aid in available:
            previous=archived_prepared.get(cid)
            if (not previous or not matching_source(src,previous['source']) or
                    prepared['state_key']!=previous.get('state_key') or
                    read_json(verify(prepared['union']))['fragments']!=read_json(verify(previous['union']))['fragments']):
                raise InvalidArtifact('archived raw source, fixed membership or chemical identity differs: '+cid)
            case,pin=available[aid];m,sm=manifest_cache[pin['sha256']]
            origins=[t for t in sm['tasks'] if t['case_id']==aid]
            if len(origins)!=2:raise InvalidArtifact('archived paired origin mapping unavailable')
            for t in origins:
                ep=prepared['representations']['context']['endpoints'][t['metal']]
                if (ep['charge']!=t['charge'] or ep['multiplicity']!=t['multiplicity'] or
                        not context.same_geometry(xyz(verify(ep['xyz'])),xyz(verify(t['xyz'])))):
                    raise InvalidArtifact('archive source/union state or coordinates differ: '+cid)
            origin=pool.score(case['matrix']['Ca']['origin']['components'],case['matrix']['La']['origin']['components'])
            row['union_static']={**origin,'R':origin['composite_R_model_kcal_mol'],
                'decision':decision(origin['composite_R_model_kcal_mol'],union_band)}
            chosen=case['pool']
            if case['status']=='prepared' and pool.choose_rows(case['matrix'],CANDIDATES)!=chosen:
                raise InvalidArtifact('archive matrix/selection algebra differs')
            if chosen['status']=='available':
                for variant in ('operational','mathematical'):
                    values=chosen[variant];value=values['composite_R_model_kcal_mol']
                    row['candidate'][variant]={**values,'R':value,'decision':decision(value,ref['variants'][variant]['bands']),
                        'delta_R_from_union_origin':value-origin['composite_R_model_kcal_mol']}
            row.update(status=chosen['status'],reason=chosen.get('reason'),source_collection=pin,
                matrix=case['matrix'],pool=chosen,candidate_geometry=case['candidates'])
        # Static compatibility is independent of candidate admission. A valid
        # static score remains visible when the candidate is unavailable.
        rows.append(row)
    index={x['case_id']:x for x in rows};aggregates=[]
    for g in req['groups']:
        if g['kind']!='ten_folds':continue
        la=[c for c in g['members'] if sources[c]['source_conditioning_metal']=='La' and c!=g['canonical_case_id']]
        ca=[c for c in g['members'] if sources[c]['source_conditioning_metal']=='Ca']
        members={'La4':la,'Ca5':ca,**{'La3_'+str(i+1):list(v) for i,v in enumerate(itertools.combinations(la,3))}}
        methods={'original_static':original_band,'union_static':union_band,
            'union_adaptive':ref['variants']['operational']['bands']}
        normalized={cid:{'methods':{'original_static':x['original_static'],'union_static':x['union_static'],
            'union_adaptive':x['candidate']['operational']}} for cid,x in index.items()}
        for name,ids in members.items():
            values={method:strict_summary(ids,normalized,method,band,'unknown') for method,band in methods.items()}
            for val in values.values():val.pop('outcome',None)
            aggregates.append({'protein_id':g['protein_id'],'descriptor':name,'members':ids,'methods':values})
        balanced={}
        for method,band in methods.items():
            arms=[next(x for x in aggregates if x['protein_id']==g['protein_id'] and x['descriptor']==a)['methods'][method] for a in ('La4','Ca5')]
            value=None if any(x['R'] is None for x in arms) else sum(x['R'] for x in arms)/2
            balanced[method]={'R':value,'decision':decision(value,band),'missing_members':sorted({v for a in arms for v in a['missing_members']})}
        aggregates.append({'protein_id':g['protein_id'],'descriptor':'balanced','members':la+ca,'methods':balanced})
    result={'protocol_id':PROTOCOL,'profile':PROFILE,'plan':record(plan),'archive':req['archive'],
        'reference':req['reference'],'rows':rows,'aggregates':aggregates,'denominator':len(rows),
        'available':sum(x['status']=='available' for x in rows),'new_molecular_calls':0,'production_changed':False,
        'affinity_probability':None,'aquo_referenced_score':None,'three_fold_union_compatible':False,
        'interpretation':'Exact archived ten-fold preparation policy; correlated structural summaries, not biological or thermal populations.'}
    candidate.put(output,result);return {k:v for k,v in result.items() if k not in ('rows','aggregates')}


def report(result, output):
    d=read_json(result);checked(verify(d['plan']))
    lines=['# Opt-in ten-fold union/adaptive archive replay','',f"Available {d['available']}/{d['denominator']}; zero molecular calls; production unchanged.",'',
        '| Source | Original static R / call | Union static R / call | Union adaptive R / call | Status |','|---|---|---|---|---|']
    for r in d['rows']:
        a=r['candidate']['operational'];b=r['original_static'];c=r['union_static']
        lines.append(f"| {r['case_id']} | {b['R']} / {b['decision']} | {c['R']} / {c['decision']} | {a['R']} / {a['decision']} | {r['status']}: {r.get('reason') or ''} |")
    lines+=['','Each score is Ca−La in model kcal/mol on its own frozen protocol scale. All ten declared folds determine fixed context membership.',
        'Strict La4/Ca5/balanced and all four La triples are retained in the JSON. Missing members remain unavailable.',
        'Three-fold union preparation, fresh scientific execution and cheaper numerical profiles are not qualified by this replay.']
    with Path(output).open('x') as f:f.write('\n'.join(lines)+'\n')
    return {'report':record(output),'new_molecular_calls':0}


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='op',required=True)
    q=sub.add_parser('archive')
    for name in ('preparation','crystals','comparison','reference','output'):q.add_argument('--'+name,required=True)
    q.add_argument('--release',default=str(standard.DEFAULT_RELEASE))
    q=sub.add_parser('request');q.add_argument('--sources',required=True);q.add_argument('--groups',nargs='*',default=[])
    q.add_argument('--archive',dest='archive_path',required=True);q.add_argument('--output',required=True)
    q.add_argument('--crystal-sources');q.add_argument('--crystals',nargs='*',default=[])
    q=sub.add_parser('prepare');q.add_argument('--request',required=True);q.add_argument('--output',required=True)
    q.add_argument('--source-mode',choices=('reuse-exact','fresh'),default='reuse-exact');q.add_argument('--agreement',default=str(DEFAULT_PLAN))
    for op in ('dry-run','replay'):
        q=sub.add_parser(op);q.add_argument('--plan',required=True)
        if op=='replay':q.add_argument('--output',required=True)
    q=sub.add_parser('report');q.add_argument('--result',required=True);q.add_argument('--output',required=True)
    args=vars(parser.parse_args());op=args.pop('op').replace('-','_');print(json.dumps(globals()[op](**args),indent=2))

if __name__=='__main__':main()
