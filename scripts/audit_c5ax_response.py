"""Read-only C5AXV8 replica, component and physical-response audit."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from functools import lru_cache
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from mace_site_kinematics import Kinematics
from adaptive_force_diagnostic import project
from nikasha_pool import score

CID='c5axv8-pqq-la_model'
@lru_cache(None)
def load(path,sha):return read_json(verify({'path':path,'sha256':sha}))
def data(pin):return load(pin['path'],pin['sha256'])

def geometry(mapping,q):
    kin=Kinematics(mapping);p=kin.positions_only(q);out={}
    for i,a in enumerate(mapping['source_atom_metadata']):
        if not a:continue
        if (a['resnum'] in (198,275,317,319) and a['atom'] in ('OE1','OE2','OD1','OD2')) or (a['resname']=='PQQ' and a['atom'] in ('O5','O7','N6')):
            key=f"{a['chain']}/{a['resnum']}/{a['resname']}/{a['atom']}";out[key]=float(np.linalg.norm(p[i]-p[kin.metal]))
    return out


def response(row,pin,canonical=False):
    coll=data(pin);case=next(c for c in coll['cases'] if c['case_id']==row['case_id']);poolm=data(coll['manifest']);mp=poolm['source_manifest'];m=data(mp)
    meta=next(c for c in m['cases'] if c['case_id']==row['case_id']);ts={t['metal']:t for t in m['tasks'] if t['case_id']==row['case_id']}
    mapping=data(ts['Ca']['mapping'])['context']
    if mapping!=data(ts['La']['mapping'])['context']:raise InvalidArtifact('paired mappings differ')
    result={'pool_collection':pin,'proposal_manifest':mp,'physical_source_case_id':meta['actual_union_case_id'],
        'state_signature':meta['origin_row']['group_state_signature'],'same_context_atoms':meta['atom_count'],
        'source_protonated':data(ts['Ca']['source_preparation'])['source'],'map':ts['Ca']['mapping'],
        'old_context_atoms':meta['union']['original_atom_count'],'core_parent':meta['union']['original_core']['parent'],
        'roles':data(meta['union']['original_core']['parent']).get('fixed_core',{}).get('requested_roles'),
        'charges':{z:t['charge'] for z,t in ts.items()},'origin_distances_A':geometry(mapping,np.zeros(len(mapping['modes']))),
        'origin_components':score(case['matrix']['Ca']['origin']['components'],case['matrix']['La']['origin']['components']),
        'pool':case['pool'],'endpoints':{},'physical_ids':mapping['physical_ids']}
    for z,t in ts.items():
        ep=(row['endpoints'] if canonical else row['precision_numerical_comparison']['endpoints'])[z];rp=ep['new_receipt'];r=data(rp);stages={}
        for stage in ('origin','proposal'):
            pt=r[stage];f=np.load(verify(pt['forces']),allow_pickle=False);q=np.asarray(pt['full_q']);coords,unit,raw,normed,lengths=project(mapping,q,f)
            if not np.allclose(coords,[a[1:] for a in xyz(verify(pt['coordinate']))],atol=1e-9,rtol=0):raise InvalidArtifact('existing coordinate/J replay differs')
            if not np.allclose(raw[t['active_indices']],pt['gradient_kcal_mol_rad'],atol=1e-8,rtol=0):raise InvalidArtifact('saved active gradient does not replay')
            modes=[{'id':mo['id'],'kind':mo['kind'],'unit':mo['unit'],'active':i in t['active_indices'],
                'raw_gradient_kcal_per_unit':float(raw[i]),'normalized_gradient_kcal_mol_A':float(normed[i]),'physical_norm_A_per_unit':float(lengths[i])} for i,mo in enumerate(mapping['modes'])]
            stages[stage]={'coordinate':pt['coordinate'],'MACE':pt['MACE'],'forces':pt['forces'],'full_q':pt['full_q'],'distances_A':geometry(mapping,q),'modes':modes,
                'metal_gradient_norm_kcal_mol_A':float(np.linalg.norm(normed[:3])),
                'maximum_omitted_angular_load_kcal_mol_A':max(abs(x['normalized_gradient_kcal_mol_A']) for x in modes if x['kind']!='metal_translation' and not x['active'])}
        result['endpoints'][z]={'receipt':rp,'active_mode_ids':t['active_mode_ids'],'selector':t['selector'],'optimizer':r['optimizer'],'boundary_flag':r['boundary_flag'],'final_geometry':r['final_geometry'],
            'force_hamiltonian':'native_vacuum_OMOL','composite_gradient':None,'stages':stages}
    origin=result['origin_components']['composite_R_model_kcal_mol']
    result['accommodation_delta_R_kcal_mol']=case['pool']['operational']['composite_R_model_kcal_mol']-origin
    result['selected_components_work']={z:case['pool']['rows'][z]['work_from_origin_kcal_mol'][case['pool']['rows'][z]['operational_candidate']] for z in ('Ca','La')}
    return result


def analyze(precision225,precision34,exclusions,strict32,output,strict225=None):
    p=read_json(precision225);p34=read_json(precision34);excluded=read_json(exclusions);strict=read_json(strict32)
    rows=[r for r in p['rows'] if r['root_case_id']==CID]
    if len(rows)!=9:raise InvalidArtifact('all nine noncanonical replicas required')
    canonical=next(r for r in p34['rows'] if r['case_id']==CID);allrows=[]
    for r in rows:
        v={'case_id':r['case_id'],'source_conditioning_metal':r['source_conditioning_metal'],'status':r['precision_status'],
            'expected_class':r['expected_class'],'methods':r['methods'],'reason':r['precision_reason'],'role':'noncanonical_consumed_replica'}
        if r['precision_status']=='available':v['response']=response(r,r['precision_collection'])
        else:v['prior_source_diagnosis']=next(x for x in excluded['rows'] if x['case_id']==r['case_id'])
        allrows.append(v)
    can={'case_id':CID,'source_conditioning_metal':'La','status':canonical['status'],'expected_class':canonical['known_class'],'role':'canonical_calibration',
         'precision_scores':canonical['scores'],'response':response(canonical,canonical['collection'],True),
         'strict32':next(r for r in strict['rows'] if r['case_id']==CID)}
    allrows.append(can)
    if len({r['response']['physical_source_case_id'] if 'response' in r else r['case_id'] for r in allrows})!=10:raise InvalidArtifact('canonical/replica source duplicates')
    if len({r['response']['state_signature']['sha256'] for r in allrows if 'response' in r})!=1:raise InvalidArtifact('La union composition differs')
    if strict225 is not None:
        s=read_json(strict225);idx={r['case_id']:r for r in s['rows']};pools={r['case_id']:r for r in s['actual_pools']}
        for r in allrows:
            if r['case_id'] in idx:
                r['strict225']=idx[r['case_id']]
                if r['case_id'] in pools:r['strict_pool']=pools[r['case_id']]
    component_spans={}
    for field in ('native_R_model_kcal_mol','solvation_delta_R_kcal_mol','composite_R_model_kcal_mol'):
        values=[r['response']['origin_components'][field] for r in allrows if 'response' in r];component_spans[field]={'min':min(values),'max':max(values),'range':max(values)-min(values)}
    result={'scope':'all10 declared C5AXV8 sources;5 supportedLa and5offsiteCa; no outcome exclusion',
        'inputs':{k:record(v) for k,v in [('precision225',precision225),('precision34',precision34),('exclusions',exclusions),('strict32',strict32)]}|{'strict225':record(strict225) if strict225 else None},
        'strict225_status':'included' if strict225 else 'not_available_at_this_snapshot','rows':allrows,'bands':p['bands'],'origin_component_spans':component_spans,
        'new_molecular_calls':0,'new_geometries':0,'new_thresholds':0,'force_interpretation':'archived vacuum MACE gradients; not composite-score gradients','implementation':record(__file__)}
    write_new(output,result);return {'output':record(output),'declared_sources':10,'supported':5,'strict225_included':bool(strict225),'new_molecular_calls':0}

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    for k in ('precision225','precision34','exclusions','strict32','output'):a.add_argument('--'+k,required=True,type=Path)
    a.add_argument('--strict225',type=Path);print(json.dumps(analyze(**vars(a.parse_args())),indent=2))
