#!/usr/bin/env python3
"""Stage exact common-baseline native bridge inputs, never invoke a native binary."""
import argparse,hashlib,json,platform,time
from pathlib import Path

IDS=('aquo9_central_La__dimer_01','aquo9_central_Dy__far__metal_two_water','matched_dry_acetate_mono_central_La')
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,('changed source',p);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):
    p=Path(p);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');return pin(p)
def stage(cfg,out):
    tick=time.monotonic();cpu=time.process_time();out.mkdir(parents=True,exist_ok=False)
    assert tuple(cfg['fixtures'])==IDS
    base=read(cfg['shared_baseline']);wa=read(cfg['water_union_admission']);da=read(cfg['donor_union_admission'])
    assert wa['shared_baseline']==cfg['shared_baseline'] and da['shared_baseline']==cfg['shared_baseline']
    assert base['old_h_H_zero'];key=verify(base['key']).read_bytes();prm=verify(base['parameter_file']).read_bytes()
    words={s.split('#')[0].split()[0].lower() for s in key.decode().splitlines() if s.split('#')[0].split()}
    assert not words&{'ewald','pme','a-axis','b-axis','c-axis','external-field','polar-local-response','polar-radial-response','polar-pair-response','exchpol','exchange-polar'}
    entries={}
    for id in IDS:
        f=read(cfg['fixtures'][id]);assert f['compatible'] and f['native_operator'] is not None
        if id in wa['checks']:assert wa['checks'][id]['same_applicable_pair_records']
        else:assert da['checks']['donor/'+id]['compatible']
        xyz=verify(f['inputs']['xyz']).read_bytes();lines=[s.split() for s in xyz.decode().splitlines() if s.strip()];assert int(lines[0][0])==len(f['geometry']['xyz_A'])
        assert [[float(s) for s in row[2:5]] for row in lines[1:]]==f['geometry']['xyz_A']
        folder=out/id;folder.mkdir();(folder/'system.xyz').write_bytes(xyz);(folder/'system.key').write_bytes(key);(folder/'amoeba09.prm').write_bytes(prm)
        entries[id]={'fixture':cfg['fixtures'][id],'geometry':f['geometry'],'baseline':cfg['shared_baseline'],'source_inputs':f['inputs'],'inputs':{'xyz':pin(folder/'system.xyz'),'key':pin(folder/'system.key'),'parameters':pin(folder/'amoeba09.prm')},'native_operator_saved':f['native_operator'],'zero_reference':f['baseline_observables'],'coordinate_identity':'exact original XYZ bytes; common-union key has proved same applicable parameters','native_execution_released':False}
    r={'status':'prepared_exact_three_vacuum_bridge_centers','implementation':pin(__file__),'configuration':cfg,'entries':entries,'phase_declaration':cfg['phase_declaration'],'native_bridge_source':cfg['native_bridge_source'],'build':None,'root_runtime_release':None,'new_native_calls':0,'new_quantum_calls':0,'operator_solves':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'host':platform.node()}}
    save(out/'MANIFEST.json',r);return r
def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=stage(json.loads(a.inputs.read_text()),a.out);print(json.dumps({'status':r['status'],'cost':r['cost']}))
if __name__=='__main__':main()
