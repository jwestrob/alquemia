"""Finite full100 scope and exact completed-pilot reuse; no molecular execution."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from nikasha_pool import choose_rows
from motion_envelope_pool import POOL_PROTOCOL
from motion_envelope_scalar import PROFILE,recipe
import slsqp_precision

PROTOCOL='Nikasha_motion_envelope4p3_all100_transfer_scope_v1'


def key(row):return row['selection_id'],row['case_id']


def audit_reuse(prepared,case,collection,manifest,search,reference):
    """Same original prepared object and exact saved maps, not fragment-list overlap."""
    if case['source']['union']!=prepared or case['pool']['status']!='available':
        raise InvalidArtifact('completed pilot does not contain this exact prepared source')
    if (manifest['protocol_id']!=POOL_PROTOCOL or manifest['numerical_policy_id']!=PROFILE or
        manifest['model']!=reference['model'] or manifest['settings']!=reference['pool_settings'] or
        search['settings']!=reference['optimizer_settings'] or search['settings']!=slsqp_precision.SETTINGS):
        raise InvalidArtifact('completed candidate protocol/model/settings differ')
    if choose_rows(case['matrix'],[x['id']for x in case['candidates']])!=case['pool']:
        raise InvalidArtifact('completed pool algebra differs')
    checked=0
    for metal in ('Ca','La'):
        task=next(t for t in search['tasks']if(t['case_id'],t['metal'])==(case['case_id'],metal))
        if task['mapping']!=prepared['maps'][metal]:raise InvalidArtifact('exact physical map pin differs')
        verify(task['mapping']);verify(task['xyz'])
        for cell in case['matrix'][metal].values():
            if cell['status']!='complete':raise InvalidArtifact('incomplete saved component')
            mace=read_json(verify(cell['MACE']))
            if mace['energy_eV']!=cell['components']['MACE_eV']:raise InvalidArtifact('native energy receipt differs')
            verify(cell['xyz'])
            for medium,low in cell['low'].items():
                if low['status']!='complete' or low['observed_TolE_hartree']!=1e-10:
                    raise InvalidArtifact('saved scalar profile unavailable')
                lm=read_json(verify(low['manifest']));lt=next(t for t in lm['tasks']if t['task_id']==low['task_id'])
                if (lt['charge'],lt['multiplicity'])!=(task['charge'],task['multiplicity']):
                    raise InvalidArtifact('scalar state differs')
                if verify(lt['input']).read_text()!=recipe(task['charge'],medium,'fresh'):
                    raise InvalidArtifact('scalar recipe differs')
                for field in ('receipt','output'):verify(low[field])
                component='GFN2_'+('vacuum'if medium=='vacuum'else'ALPB')+'_hartree'
                if low['energy_hartree']!=cell['components'][component]:raise InvalidArtifact('scalar component differs')
                checked+=1
        alias=case['aliases']['adaptive_'+metal]
        actual=read_json(verify(alias['proposal_receipt']))
        if actual['status']!='proposal_available' or actual['manifest']!=manifest['source_manifest']:
            raise InvalidArtifact('actual optimizer receipt differs')
    return {'status':'exact_pilot_pool_reuse','collection':record(collection),'case_id':case['case_id'],
        'source_manifest':manifest['source_manifest'],'physical_source_row_exact':True,
        'mapping_pins_exact':True,'scalar_cells_verified':checked,'pool':case['pool']}


def inventory(preparation,prior_audit,pilot,reference,agreement,output):
    p=read_json(preparation);a=read_json(prior_audit);d=read_json(pilot);ref=read_json(reference)
    selection=read_json(verify(p['selection']));pm=read_json(verify(d['manifest']));sm=read_json(verify(pm['source_manifest']))
    if a['preparation']!=record(preparation) or a['selection']!=p['selection']:
        raise InvalidArtifact('prior all-source physical archive audit differs')
    if ref['collection']!=record(pilot)or any(v['status']!='available'for v in ref['variants'].values()):
        raise InvalidArtifact('actual pilot reference unavailable or belongs to another collection')
    verify(ref['implementation']);verify(selection['agreement']);record(agreement)
    triples=p['triples'];primary=[r for r in p['cases']if'primary_triple_member'in r['roles']]
    if (len(triples),sum(t['status']=='prepared'for t in triples),len(primary))!=(100,94,104):
        raise InvalidArtifact('declared100/94/104 ledger differs')
    audits={key(r):r for r in a['rows']};completed={key(c['source']['union']):c for c in d['cases']}
    rows=[];new=[];reused=[]
    for r in primary:
        old=audits[key(r)]
        if old['origin_reuse']or old['pool_reuse']or old['candidate_coordinate_matches_requiring_further_state_audit']:
            raise InvalidArtifact('additional old-coordinate match requires explicit audit')
        if r['status']!='prepared':raise InvalidArtifact('unexpected missing source; retain and revise inventory')
        for pin in r['maps'].values():verify(pin)
        entry={'case_id':r['case_id'],'selection_id':r['selection_id'],'triple_ids':r['triple_ids'],
            'source':r['source'],'state_key':r['state_key'],'atom_count':r['atom_count'],
            'La_charge':r['new_La_charge'],'maps':r['maps'],'prepared':r,'status':'new_pool_required','reuse':None}
        if key(r)in completed:
            entry['reuse']=audit_reuse(r,completed[key(r)],pilot,pm,sm,ref);entry['status']='exact_pilot_pool_reuse';reused.append(entry)
        else:new.append(entry)
        rows.append(entry)
    lookup={key(r)for r in rows}
    for t in triples:
        if t['status']=='prepared'and not all((t['selection_id'],cid)in lookup for cid in t['members']):
            raise InvalidArtifact('prepared triple lacks exact member')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shards=[]
    for i,part in enumerate((new[:51],new[51:])):
        tasks=[]
        for row in part:
            r=row['prepared'];cid=r['case_id']+'__envelope_'+r['selection_id'].rsplit('_',1)[-1]
            for metal,ep in r['representations']['context']['endpoints'].items():
                verify(ep['xyz']);tasks.append({'task_id':cid+'__context__'+metal,'case_id':cid,'source_case_id':r['case_id'],
                    'selection_id':r['selection_id'],'metal':metal,'xyz':ep['xyz'],'charge':ep['charge'],'multiplicity':ep['multiplicity']})
        count=len(part)
        shard={'protocol_id':PROTOCOL,'shard':i,'cases':part,'origin_tasks':tasks,
            'counts':{'sources':count,'origin_MACE':2*count,'searches':2*count,'maximum_cross_MACE':2*count,
                'origin_GFN2':4*count,'maximum_candidate_GFN2':8*count},
            'execution_status':'scope_only_not_a_runnable_scientific_manifest','molecular_calls':0}
        path=out/f'shard_{i}/SCOPE.json';write_new(path,shard);shards.append(record(path))
    counts={'triple_denominator':100,'complete_preparations':94,'prior_unavailable_triples':6,
        'source_membership_pairs':len(rows),'whole_pilot_pool_reuses':len(reused),'new_pools':len(new),
        'new_origin_MACE':2*len(new),'new_searches':2*len(new),'maximum_new_cross_MACE':2*len(new),
        'new_origin_GFN2':4*len(new),'maximum_new_candidate_GFN2':8*len(new),
        'maximum_new_GFN2_total':12*len(new),'new_DFT':0,'new_folds':0}
    result={'protocol_id':PROTOCOL,'preparation':record(preparation),'prior_archive_audit':record(prior_audit),
        'pilot_collection':record(pilot),'reference':record(reference),'agreement':record(agreement),
        'frozen_UTC':datetime.now(timezone.utc).isoformat(),'implementation':record(__file__),
        'scalar_profile':PROFILE,'model':ref['model'],'optimizer_settings':ref['optimizer_settings'],
        'pool_settings':ref['pool_settings'],'triples':triples,'rows':rows,'shards':shards,'counts':counts,
        'additional_old_origin_reuse':0,'molecular_calls_in_inventory':0,'new_reference_fit':False,
        'proposed_resources':{'warm_GPU_per_shard':{'CPU':32,'memory_MiB':200000,'H200':1},
            'scalar_per_shard':{'CPU':32,'memory_MiB':65536,'workers':32,'ranks_per_cell':1}},
        'execution_status':'awaiting_parent_finite_scope_coordination','production_changed':False}
    write_new(out/'INVENTORY.json',result);return {'inventory':record(out/'INVENTORY.json'),'counts':counts,'shards':shards}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('preparation','prior_audit','pilot','reference','agreement','output'):p.add_argument('--'+k.replace('_','-'),required=True)
    print(json.dumps(inventory(**vars(p.parse_args())),indent=2))
