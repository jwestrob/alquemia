"""Exact reuse audit of new-ftol pools for unchanged triple-membership contexts."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from second_shell_context import parent_state
from coordination_preparation_context import geometry
import consistent_context as context
import slsqp_precision as precision
import slsqp_precision_expansion as expansion
import nikasha_pool as pool
from compact_solvation import completed,input_text
import union_triple_preparation as prep


def audit(canonical,reference,agreement,output):
    start=time.monotonic();c=read_json(canonical);m=prep.validate_selection(verify(c['selection']));ref=read_json(reference)
    if ref['optimizer_settings']!=precision.SETTINGS or ref['model']!=m['config']['model']:
        raise InvalidArtifact('actual new-ftol profile/model required')
    saved={};validations=[]
    for pin in ref['collections']:
        d=read_json(verify(pin));mp=verify(d['manifest']);pm=read_json(mp);sm=read_json(verify(pm['source_manifest']))
        module=precision if len(sm['cases'])==4 else expansion
        validations.append(module.validate_pool(mp))
        for row in d['cases']:
            source=row['source'];actual=source['actual_union_case_id']
            if actual in saved:raise InvalidArtifact('ambiguous actual source in new-ftol pools')
            saved[actual]=(row,pin,pm,sm)
    rows=[]
    for prepared in c['cases']:
        cid=prepared['case_id'];r={'case_id':cid,'status':'new_pool_required','reason':'changed_context','pool_collection':None}
        if prepared['tenfold_comparison']['exact_reuse_eligible']:
            if cid not in saved:raise InvalidArtifact('expected unchanged-context new-ftol pool missing: '+cid)
            row,pin,pm,sm=saved[cid];source=row['source'];old=source['union']
            if (old['state_key']!=prepared['state_key'] or old['original_core']!=prepared['original_core'] or
                    sm['settings']!=precision.SETTINGS or pm['settings']!=ref['settings'] or
                    row['pool']['status']!='available' or source['actual_union_case_id']!=cid):
                raise InvalidArtifact('source graph/core/state/settings unavailable: '+cid)
            state=parent_state(prepared['original_core'],m['config']['topology'],require_endpoint_receipts=False)
            pa=read_json(verify(prepared['representations']['context']['preparation']));maps={}
            for z in ('Ca','La'):
                t=next(t for t in sm['tasks'] if (t['case_id'],t['metal'])==(row['case_id'],z))
                ep=prepared['representations']['context']['endpoints'][z]
                if not context.reusable_state(ep,t):raise InvalidArtifact('origin coordinates/charge differ')
                g=geometry(state,pa,xyz(verify(ep['xyz'])),xyz(verify(prepared['original_core']['endpoints'][z]['xyz'])))
                g=json.loads(json.dumps(g));prior=read_json(verify(t['mapping']))
                if g!=prior:raise InvalidArtifact('physical coordinate/cap mapping differs: '+cid+' '+z)
                maps[z]=t['mapping']
                for name,cell in row['matrix'][z].items():
                    native=read_json(verify(cell['MACE']))
                    if cell['components']['MACE_eV']!=native['energy_eV'] or native['status'] not in ('computed','complete'):
                        raise InvalidArtifact('actual native component differs')
                    for medium in ('vacuum','alpb'):
                        oldlow=cell['low'][medium];actual=completed(verify(oldlow['manifest']),oldlow['task_id'])
                        if actual is None or any(oldlow[k]!=actual[k] for k in actual):raise InvalidArtifact('actual low-level receipt differs')
                        if actual['energy_hartree']!=cell['components']['GFN2_'+('vacuum' if medium=='vacuum' else 'ALPB')+'_hartree']:
                            raise InvalidArtifact('component algebra differs')
                        lm=read_json(verify(actual['manifest']));lt=next(t for t in lm.get('all_tasks',lm['tasks']) if t['task_id']==actual['task_id'])
                        if (not context.same_geometry(xyz(verify(cell['xyz'])),xyz(verify(lt['xyz']))) or
                                lt['charge']!=ep['charge'] or lt['multiplicity']!=ep['multiplicity'] or
                                verify(lt['input']).read_text().replace(' MaxIter 500\n','')!=input_text(ep['charge'],ep['multiplicity'],medium,'native')):
                            raise InvalidArtifact('actual candidate low-level physical state/recipe differs')
            if pool.choose_rows(row['matrix'],[x['id'] for x in row['candidates']])!=row['pool']:
                raise InvalidArtifact('saved minimal-pool selection algebra differs')
            r.update(status='exact_new_ftol_pool_reuse',reason=None,pool_collection=pin,source_case_id=row['case_id'],
                source_manifest=pm['source_manifest'],mapping=maps,source=row['source'],pool=row,
                numerical_provenance='actual archived eight-rank scalar receipts; rank1 separately qualified, no relabeling')
        rows.append(r)
    if sum(r['status']=='exact_new_ftol_pool_reuse' for r in rows)!=19 or len(rows)!=28:
        raise InvalidArtifact('expected19 exact and9 changed canonical/control contexts required')
    result={'protocol_id':'Nikasha_three_La_membership_new_ftol_pool_reuse_v1','canonical_preparation':record(canonical),
        'existing_new_ftol_reference':record(reference),'agreement':record(agreement),'rows':rows,
        'denominator':28,'complete_pool_reuses':19,'new_pools':9,'separate_stress_new_pools':1,
        'source_pool_validations':validations,'new_molecular_calls':0,'wall_seconds':time.monotonic()-start,
        'implementation':record(__file__),'new_reference':None}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('rows','source_pool_validations')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('canonical','reference','agreement','output'):p.add_argument('--'+k,required=True)
    print(json.dumps(audit(**vars(p.parse_args())),indent=2))
