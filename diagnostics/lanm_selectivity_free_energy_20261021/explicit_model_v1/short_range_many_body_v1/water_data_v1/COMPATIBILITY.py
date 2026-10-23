#!/usr/bin/env python3
"""Saved-only chemical/input compatibility and separate own-basis labels."""
import argparse, json, time
from pathlib import Path
from MATERIALIZE import pin,read,save,verify,geometry

def tokens(inp,drop_field=False):
    ignored={'parameters','polar-print','save-uinduce','response-matrix-export','digits','openmp-threads'}
    if drop_field:ignored.add('external-field')
    lines=[s.split('#')[0].split() for s in verify(inp['key']).read_text().splitlines()]
    return [t for t in lines if t and t[0].lower() not in ignored]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    tick=time.monotonic();cpu=time.process_time();d=read(pin(a.dataset));fixtures=read(d['fixtures'])['fixtures'];forces=read(d['forces']);fields=read(d['fields']);static=read(d['static'])
    checks={};asym={};rhs={};allgeom={}
    for id,p in fixtures.items():
        f=read(p);allgeom[id]=f['geometry'];op=f['native_operator'];A=op['A0_Aminus3']
        asym[id]=max(abs(A[i][j]-A[j][i]) for i in range(len(A)) for j in range(len(A)))
        rhs[id]=max(abs(x-y) for x,y in zip(op['Ed_e_per_A2'],op['Ep_e_per_A2']))
        checks[id]={'water_only_types':set(f['geometry']['types'])<={36,37,500,501},'donor_override_types_403_105_absent':not set(f['geometry']['types'])&{403,105},'old_local_response_zero':f['baseline_identity']['old_h_H_zero'],'no_array_modification':True}
        assert all(checks[id].values())
    allgeom.update({id:s['geometry'] for id,s in static.items()})
    forcejoins={}
    for id,r in forces.items():
        f=read(fixtures[id]);old=r['baseline_model']['provenance']['effective_input'];new=f['inputs']
        same=(old['xyz']['sha256']==new['xyz']['sha256'] and old['parameters']['sha256']==new['parameters']['sha256'] and tokens(old)==tokens(new))
        assert same, ('force baseline mismatch',id)
        forcejoins[id]={'same_effective_model_and_geometry':True,'native_authority':r['baseline_model']['provenance']['native_authority']}
    fieldjoins={}
    for id,r in fields.items():
        f=read(r['baseline_zero_operator']);inp=r['model_finite_field_record']['runtime_input'];zero=f['inputs']
        same=(inp['xyz']['sha256']==zero['xyz']['sha256'] and inp['parameters']['sha256']==zero['parameters']['sha256'] and tokens(inp,True)==tokens(zero))
        assert same, ('field baseline mismatch',id)
        fieldjoins[id]={'same_zero_field_physics_and_geometry':True,'declared_field_only_change':r['model_finite_field_record']['job']['field_au']}
    targets=read(d['inputs']['targets']);qm=read(targets['artifacts']['quantum_full_precision_review']);common=read(qm['collection']);pair=read(common['pair_collection']);cfg=read(pair['configuration']);panel=read(cfg['reference_panel']);cells={x['cell_id']:x for x in panel['cells']};own={}
    for id,r in pair['reference_rows'].items():
        cell=cells[id];parts=id.split('__')
        if len(parts)==4:mid=f'{parts[1]}_central_{parts[2]}__dimer_{parts[3][6:]}'
        else:mid='ion_'+cell['metal']
        idx=cell.get('real_parent_atom_indices',list(range(len(cell['symbols']))))
        xyz=[cell['positions_A'][i] for i in idx]
        assert xyz==allgeom[mid]['xyz_A'],('own basis geometry',id,mid)
        own[id]={'model_fixture_id':mid,'reference_cell':cell,'native_reference':r,'consumed_development':True,'reference_variant':'original matched own-fragment basis; distinct from common-parent-basis label','common_basis_target_replacement':False,'same_real_geometry':True}
    ownpin=save(a.out/'MATCHED_OWN_BASIS.json',{'source':common['pair_collection'],'panel':cfg['reference_panel'],'rows':own,'reference_decomposition':pair['reference_decomposition'],'original_clusters':pair['original_clusters']})
    report={'status':'saved_water_baseline_compatibility_checked','dataset':pin(a.dataset),'implementation':pin(__file__),'fixture_checks':checks,'force_baseline_joins':forcejoins,'field_baseline_joins':fieldjoins,'raw_A_max_asymmetry_Aminus3_by_fixture':asym,'raw_RHS_channel_max_difference_by_fixture':rhs,'raw_arrays_symmetrized':False,'matched_own_basis':ownpin,'missing':{'static_operators':list(static),'baseline_coordinate_derivatives':['dA0/dR','dEd/dR','dEp/dR','permanent moment frame derivatives'],'later_radial14_candidate_operators':'excluded; different fitted baseline and h/H','wet_donor17621':'not represented here; old water parameters incompatible with shared water baseline'},'limits':['This is saved-array/input compatibility, not a new stability or reference-accuracy qualification.','Full forces cannot be predicted by holding baseline operator and fields fixed.','No new SCF references replace the original targets; precision comparison remains separate.'],'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1},'native_calls':0,'quantum_calls':0,'operator_solves':0}
    save(a.out/'COMPATIBILITY.json',report);print(json.dumps({'status':report['status'],'matched_own_basis_states':len(own),'nonzero_raw_asymmetry':sum(x!=0 for x in asym.values()),'cost':report['cost']}))
if __name__=='__main__':main()
