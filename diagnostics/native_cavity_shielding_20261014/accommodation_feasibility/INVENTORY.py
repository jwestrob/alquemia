#!/usr/bin/env python3
"""Read-only artifact/identity audit; no coordinate arithmetic or scientific imports."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
VERIFIED={}


def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    if isinstance(p,dict):
        if 'path' in p and 'sha256' in p:
            if p['path'] not in VERIFIED:VERIFIED[p['path']]=pin(p['path'])
            assert VERIFIED[p['path']]['sha256']==p['sha256'],p['path']
        else:
            for v in p.values():verify(v)
    elif isinstance(p,list):
        for v in p:verify(v)


def read(p):
    if isinstance(p,dict):verify(p);p=p['path']
    return json.loads(Path(p).read_text())


def main():
    sources={
        'inventory':pin(ROOT/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json'),
        'origin_preparation':pin(ROOT/'workspaces/native_cavity_shielding_20261014/preparation_v1/execution/PREPARATION.json'),
        'origin_native_collection':pin(ROOT/'workspaces/native_cavity_shielding_20261014/native_v1/execution/FINAL_COLLECTION.json'),
        'component_v2_collection':pin(ROOT/'workspaces/ca_coverage_recovery_20261012/two_ca_solvent_q9_v1/FINAL_COLLECTION.json'),
        'native_manifest':pin(ROOT/'workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/manifest.json')}
    assert sources['inventory']['sha256']=='757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343'
    inv=read(sources['inventory']);orig=read(sources['origin_preparation']);components=read(sources['component_v2_collection'])
    assert orig['status']=='prepared' and orig['prepared_endpoints']==4
    component_rows={r['source_case_id']:r for r in components['rows']}
    rows=[];motions=[]
    for source in inv['rows']:
        for k in ('source_preparation','state_signature','maps','prepared_admission','spectator_admission'):verify(source[k])
        cp=read(source['source_preparation']);maps=read(source['maps']['Ca'])['context'];verify(cp['source'])
        assert source['maps']['Ca']['sha256']==source['maps']['La']['sha256']
        local_physical={link[2]:link[0] for link in maps['core_links'] if link[1]=='source'}
        for owner in ('Ca','La'):
            prop=read(source['native_proposals'][owner]);manifest=read(prop['manifest'])
            task,=[t for t in manifest['tasks'] if t['metal']==owner]
            assert prop['status']=='proposal_available' and prop['proposal']['geometry_checks']['pass']
            assert prop['proposal']['full_q'][:3]==[0.0,0.0,0.0]
            moving=sorted({i for k in task['active_indices'] for i in maps['modes'][k]['moving_indices']})
            unmapped=[i for i in moving if i not in local_physical]
            boundary=[]
            for link in maps['core_links']:
                if link[1]=='cap' and (link[2] in moving or link[3] in moving):
                    boundary.append({'local_cap_index':link[0],'retained':maps['source_atom_metadata'][link[2]],
                                     'omitted':maps['source_atom_metadata'][link[3]],'raw_link':link})
            motion={'source_case_id':source['case_id'],'candidate':'adaptive_'+owner,
                    'proposal':source['native_proposals'][owner],'proposal_manifest':prop['manifest'],
                    'active_indices':task['active_indices'],'active_modes':[maps['modes'][i] for i in task['active_indices']],
                    'full_q':prop['proposal']['full_q'],'saved_geometry_admission':prop['proposal']['geometry_checks'],
                    'moving_physical_members':[{'physical_index':i,'local_index':local_physical.get(i),'source':maps['source_atom_metadata'][i]} for i in moving],
                    'moving_members_missing_from_local_block':unmapped,'cap_links_touching_mode_members':boundary,
                    'boundary_flag':prop['boundary_flag'],'physical_positions_count':len(maps['positions_A'])}
            verify(motion['proposal']);verify(motion['proposal_manifest']);motions.append(motion)
        for metal in ('Ca','La'):
            origin,=[r for r in orig['rows'] if r['source_case_id']==source['case_id'] and r['metal']==metal]
            for candidate in ('origin','adaptive_Ca','adaptive_La'):
                row={'source_case_id':source['case_id'],'role':source['role'],'metal':metal,'candidate':candidate,
                     'status':'unavailable','reason':None,'source_preparation':source['source_preparation'],
                     'state_signature':source['state_signature'],'maps':source['maps'],'physical_source':cp['source'],
                     'origin_union':origin['variants']['union'],'origin_cap_removed':origin['variants']['cap_removed']}
                try:
                    cell=source['matrix'][metal][candidate];native=read(cell['native_frozen_ALPB_result']);ni=native['input']
                    comp=component_rows[source['case_id']]['components'][metal+'__'+candidate]
                    assert native['status']=='available' and native['manifest']==sources['native_manifest']
                    assert ni['xyz']==cell['xyz'] and ni['maps']==source['maps'] and ni['state_signature']==source['state_signature']
                    assert ni['source_preparation']==source['source_preparation'] and ni['atom_count']==cell['atom_count']==192
                    assert ni['charge']==cell['charge'] and ni['multiplicity']==cell['multiplicity']==1
                    assert ni['valence_electrons']==cell['GFN2_valence_electrons']==588
                    assert comp['accepted'] and all(c['passed'] and c['informative'] for c in native['checks'])
                    assert all(c['informative'] for c in native['vacuum_precision'])
                    for medium in ('vacuum','water'):
                        e=ni['endpoints'][medium];point=read(e['point']);ex=read(e['execution'])
                        assert e['point']==comp['inputs'][medium] and e['status']=='available' and ex['returncode']==0
                        assert e['executed_xyz']['sha256']==cell['xyz']['sha256']==native[medium]['xyz']['sha256']
                        assert e['charges']['sha256']==native[medium]['charges']['sha256']
                        assert all(point['backends'][b]['status']=='available' for b in ('xtb','tblite'))
                        assert point['backends']['xtb']['solver_evidence']['convergence_declared']
                        assert point['backends']['tblite']['solver_evidence']['library_returned_without_error']
                        verify(e);verify(native[medium])
                    verify(cell['xyz']);verify(cell['MACE']);verify(cell['MACE_request']);verify(comp['verified_output_evidence'])
                    row.update(status='available',xyz=cell['xyz'],charge=cell['charge'],multiplicity=1,atom_count=192,valence_electrons=588,
                               MACE=cell['MACE'],MACE_request=cell['MACE_request'],native_result=cell['native_frozen_ALPB_result'],
                               vacuum_charges=ni['endpoints']['vacuum']['charges'],vacuum_SCC=ni['endpoints']['vacuum'],
                               native_vacuum_replay=native['vacuum'],native_components_Eh=native['vacuum']['components_Eh'],
                               native_diagonal_checks=native['checks'],vacuum_precision=native['vacuum_precision'],
                               component_v2={k:v for k,v in comp.items() if k!='verified_output_evidence'})
                except Exception as exc:row.update(reason=f'{type(exc).__name__}: {exc}')
                rows.append(row)
    nonorigin=[r for r in rows if r['candidate']!='origin']
    result={'status':'complete_inventory','scope':'metadata, identity and hashes only; no numerical coordinate/energy analysis',
            'sources':sources,'declared_local_cells':12,'declared_nonorigin_cells':8,
            'available_nonorigin_cells':sum(r['status']=='available' for r in nonorigin),
            'missing_nonorigin_cells':[r for r in nonorigin if r['status']!='available'],
            'rows':rows,'saved_motion_metadata':motions,'verified_artifacts':list(VERIFIED.values()),
            'new_coordinates':False,'new_scientific_calls':0,'union_accommodation_admission':'not_executed'}
    with (OUT/'INVENTORY.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'inventory':pin(OUT/'INVENTORY.json'),'available_nonorigin_cells':result['available_nonorigin_cells'],
                      'missing_nonorigin_cells':len(result['missing_nonorigin_cells']),
                      'unmapped_motion_members':[r['moving_members_missing_from_local_block'] for r in motions],
                      'verified_artifacts':len(VERIFIED)}))


if __name__=='__main__':main()
