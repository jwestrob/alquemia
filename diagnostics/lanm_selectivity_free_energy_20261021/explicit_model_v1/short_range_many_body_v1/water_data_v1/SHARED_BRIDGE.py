#!/usr/bin/env python3
"""Metadata proof that disjoint donor overrides leave the saved water operator intact."""
import argparse,json,time
from pathlib import Path
from MATERIALIZE import read,pin,verify,save

IGNORE={'parameters','polar-print','save-uinduce','response-matrix-export','digits','openmp-threads'}
def val(s):
    try:return float(s.replace('D','E'))
    except ValueError:return s
def records(p):return [t for line in verify(p).read_text().splitlines() if (t:=line.split('#')[0].split())]
def identity(inp,types):
    allrec=records(inp['parameters'])+records(inp['key']);key=records(inp['key'])
    classes={int(t[1]):int(t[2]) for t in allrec if t[0].lower()=='atom'}
    present_classes={classes[t] for t in types};pairs={}
    for t in allrec:
        k=t[0].lower()
        if k not in ('vdwpr','polpair'):continue
        a,b=map(int,t[1:3]);present=present_classes if k=='vdwpr' else set(types)
        if a in present and b in present:pairs[(k,*sorted((a,b)))]=[val(s) for s in t[3:]]
    body=[[val(s) for s in t] for t in key if t[0].lower() not in IGNORE|{'vdwpr','polpair'}]
    return body,{':'.join(map(str,k)):v for k,v in sorted(pairs.items())},{str(t):classes[t] for t in set(types)}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',type=Path,required=True);ap.add_argument('--shared-baseline',type=Path,required=True);ap.add_argument('--compatibility',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();tick=time.monotonic();cpu=time.process_time()
    data=read(pin(a.dataset));shared=read(pin(a.shared_baseline));verify(shared['key']);verify(shared['parameter_file']);verify(shared['water_authority']);own=read(data['baseline'])
    assert read(shared['water_authority'])==own
    inputs={'key':shared['key'],'parameters':shared['parameter_file']};fixtures=read(data['fixtures'])['fixtures'];static=read(data['static']);checks={}
    for id,r in [(id,read(p)) for id,p in fixtures.items()]+list(static.items()):
        typ=r['geometry']['types'];left=identity(r['inputs'],typ);right=identity(inputs,typ)
        assert r['inputs']['parameters']['sha256']==shared['parameter_file']['sha256']
        assert left==right,('shared union changes applicable physics',id)
        checks[id]={'same_nonpair_key_body':True,'same_parameter_file':True,'same_applicable_pair_records':True,'active_pair_records':left[1],'type_to_class':left[2],'donor_types_absent':not set(typ)&{403,105},'donor_classes_absent':not set(left[2].values())&{403,55}}
        assert checks[id]['donor_types_absent'] and checks[id]['donor_classes_absent']
    compat=read(pin(a.compatibility));assert compat['dataset']==pin(a.dataset)
    bridge={'status':'water_exactly_compatible_with_explicit_disjoint_shared_baseline','implementation':pin(__file__),'water_dataset':pin(a.dataset),'shared_baseline':pin(a.shared_baseline),'shared_key':shared['key'],'checks':checks,'no_saved_input_or_array_mutation':True,'native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    bridgepin=save(a.out/'SHARED_BASELINE_BRIDGE.json',bridge)
    save(a.out/'BUNDLE.json',{'status':'materialized_water_data_complete','schema':'short_range_many_body_water_bundle_v1','dataset':pin(a.dataset),'compatibility':pin(a.compatibility),'matched_own_basis':compat['matched_own_basis'],'canonical_shared_baseline':pin(a.shared_baseline),'shared_baseline_bridge':bridgepin,'reference_lookup':'DATASET.reference_model_joins is authoritative; per-fixture reference_ids is not an embedded target list.','force_claim':'Full force requires baseline derivative hooks; frozen cached operators alone are insufficient.','reference_variants':'Common-parent-basis, own-basis and total E/force records remain distinct; no target substitution.'})
    print(json.dumps({'status':bridge['status'],'compatible_water_states':len(checks),'cost':bridge['cost']}))
if __name__=='__main__':main()
