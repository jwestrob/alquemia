#!/usr/bin/env python3
"""Metadata-only inventory: copy existing evidence, never evaluate geometry/energy."""
import argparse, hashlib, json, re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
METALS=('Ca','La');ORDER=('origin','adaptive_Ca','adaptive_La')
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,('changed input',p);return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')

def main(output):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    pack=ROOT/'workspaces/exaa_source_packing_20261014/packing_v1/manifest.json'
    comparison=ROOT/'workspaces/exaa_source_packing_20261014/packing_comparison_v1/execution/RESULT.json'
    scores=ROOT/'workspaces/pqq_inconclusive_review_20261014/scores_v1/execution/RESULT.json'
    dft=ROOT/'workspaces/accommodation_fold_DFT_20260920/run_v1/final_collection.json'
    old_audit=ROOT/'diagnostics/orca_parameter_audit_20261012/q9_response_diagnosis/DFT_INVENTORY.json'
    pm=read(pack);cr=read(comparison);sc=read(scores);rows=[]
    wanted=['q9z4j7-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-'+str(i) for i in (3,0)]
    for cid in wanted:
        request=next(r for r in pm['requests'] if r['case_id']==cid);r=request['record'];p=request['prepared_source']
        s=next(x for x in sc['rows'] if x['record_id']==r['record_id']);c=next(x for x in cr['rows'] if x['case_id']==cid)
        assert s['status']=='available' and c['status']=='available'
        state=read(r['state_signature']);audit=read(r['source_preparation']);assert audit['new_atom_count']==192
        assert state['endpoint_charges']=={'Ca':-3,'La':-2} and state['multiplicities']=={'Ca':1,'La':1} and not state['water_inventory']
        assert s['case_key']==request['original_review_identity']['case_key'];matrix={z:{} for z in METALS}
        for z in METALS:
            for q in ORDER:
                x=r['matrix'][z][q];verify(x['xyz']);mm=read(x['MACE']);v=s['matrix'][z][q];n=read(v['native_result']);t=n['input']
                assert n['status']=='available' and mm['energy_eV']==x['MACE_eV']==v['MACE_eV']
                assert t['xyz']['sha256']==x['xyz']['sha256'] and t['state_signature']==r['state_signature']
                assert t['atom_count']==192 and t['charge']==state['endpoint_charges'][z] and t['multiplicity']==1
                rq=read(mm['request']);assert rq['xyz']['sha256']==x['xyz']['sha256']
                assert n['G_qvac_hartree']==v['G_qvac_hartree'] and n['vacuum']['components_Eh']==v['native_components_qvac_hartree']
                matrix[z][q]={'xyz':x['xyz'],'atom_count':t['atom_count'],'charge':t['charge'],'multiplicity':t['multiplicity'],
                    'GFN2_valence_electrons':t['valence_electrons'],'MACE':x['MACE'],'MACE_eV':x['MACE_eV'],
                    'MACE_request':mm['request'],'MACE_model':mm['model'],'native_frozen_ALPB_result':v['native_result'],
                    'native_frozen_ALPB_hartree':v['G_qvac_hartree'],'native_components_hartree':v['native_components_qvac_hartree'],
                    'self_consistent_GFN2_correction_hartree_diagnostic_only':v['U_self_consistent_hartree'],
                    'GFN2_endpoints':t['endpoints'],'saved_origin_relative_work_kcal_mol':c['small']['works'][z][q],
                    'compatible_DFT':None,'compatible_DFT_wavefunction':None}
        for path in r['maps'].values():verify(path)
        for path in request['native_proposals'].values():verify(path)
        rows.append({'case_id':cid,'role':'target' if cid.endswith('sample-3') else 'control','record_id':r['record_id'],
            'case_key':request['original_review_identity']['case_key'],'source_conditioning_metal':'2Ca_v1','selection_id':r['selection_id'],
            'state_signature':r['state_signature'],'state':state,'maps':r['maps'],'source_preparation':r['source_preparation'],
            'prepared_admission':r['prepared_admission'],'full_source':p['source']['source_structure'],
            'source_geometry_audit':p['source']['source_geometry_audit'],'spectator_admission':p['spectator_admission'],
            'original_core':p['original_core'],'native_proposals':request['native_proposals'],'matrix':matrix,
            'saved_selection':c['small']['selected'],'saved_same_candidate_differential_works':c['small']['same_candidate_differential_work'],
            'provenance':'192atom original context from packing preparation manifest, never the249atom expanded matrix',
            'DFT_status':'no_exact_compatible_completed_reference_identified'})
    known=read(dft);history=[]
    for r in known['rows']:
        if r['biological_group']!='q9z4j7':continue
        eps={}
        for z,e in r['endpoints'].items():
            ip=verify(e['input']);xp=verify(e['xyz']);op=verify(e['output']);verify(e['receipt'])
            lines=ip.read_text().splitlines();geometry=[line.strip() for line in lines if line.strip().lower().startswith('* xyz')]
            eps[z]={'status':e['status'],'input':e['input'],'xyz':e['xyz'],'output':e['output'],'receipt':e['receipt'],
                'XYZ_atom_count_header':int(xp.read_text().splitlines()[0]),'input_geometry_directive':geometry,
                'method_headers':[line for line in lines if line.strip().startswith('!')],
                'same_directory_wavefunctions':[pin(p) for p in sorted(op.parent.glob('*.gbw'))],
                'compatible_with_new_inventory':False,'reason':'Different La-conditioned source and smaller73atom context; current two-Ca192atom states/coordinates are not these endpoints.'}
        history.append({'source_case_id':r['case_id'],'conditioning':r['source_conditioning_metal'],
            'preparation_status':r['preparation_status'],'preparation_reason':r['preparation_reason'],'endpoints':eps})
    result={'protocol':'ExaA_twoCa_saved192atom_donor_response_inventory_v1','status':'metadata_ready','rows':rows,
        'candidate_order':ORDER,'declared_sources':2,'candidate_geometries_per_source':3,'metal_geometry_cells':12,
        'inputs':{'packing_manifest':pin(pack),'saved_comparison':pin(comparison),'saved_scores':pin(scores),
                  'historical_DFT_collection':pin(dft),'prior_bounded_DFT_audit':pin(old_audit),'implementation':pin(__file__)},
        'historical_DFT_source_rows':history,'exact_compatible_DFT_cells':[],'exact_compatible_DFT_wavefunctions':[],
        'search_limits':'Named native DFT collection and prior donor-response audit; original current two-Ca source trees contain no native DFT result. Not an exhaustive search of every external archive.',
        'state_caveats':['Dry fixed oxidized PQQ3minus hypothesis; no new protonation/redox/occupancy states.',
          'Spectator Ca remains in full source outside local context.',
          'GFN2 valence-electron counts are backend metadata, not a DFT all-electron/ECP electron-count prescription.',
          'Old73atom La-conditioned DFT origins and249atom packing endpoints are not substitutes for these exact192atom inputs.'],
        'numerical_verification':'No new geometry/numerical evaluation. Prior real preparation and saved-output gates remain pinned; root must qualify any DFT inputs on Slurm.',
        'new_molecular_calls':0,'submissions':0,'energy_units':{'MACE':'eV','native_ALPB':'hartree','saved_works':'kcal/mol'}}
    write(out/'INVENTORY.json',result);print(json.dumps(pin(out/'INVENTORY.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();main(a.output)
