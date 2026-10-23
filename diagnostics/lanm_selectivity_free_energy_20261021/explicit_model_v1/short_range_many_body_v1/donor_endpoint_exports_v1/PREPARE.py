"""Copy eight existing La donor endpoints and bind saved targets; no numerical code."""
import argparse, ast, hashlib, json, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
REL = Path('lanm_selectivity_free_energy_20261021/explicit_model_v1/short_range_many_body_v1')
D, W = REPO/'diagnostics'/REL, REPO/'workspaces'/REL
IDS = tuple(f'pbe0__La__{direction}__{shift}__complex' for direction in ('radial','rotation')
            for shift in ('minus_010','minus_005','plus_005','plus_010'))

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p;return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    start,cpu=time.monotonic(),time.process_time()
    source=pin(W/'donor_data_v1/displaced_v1/RESULT.json');assert source['sha256'].startswith('64997f1a')
    center=pin(W/'donor_data_v1/union_v1/kernel_views/matched_dry_acetate_mono_central_La.json');assert center['sha256'].startswith('c6ec944d')
    old=read(source);f=read(center);baseline=f['baseline'];base=read(baseline)
    assert f['compatible'] and base['old_h_H_zero'] and old['baseline']==baseline
    cfgpath=W/'native_force_qualification_v1/scout_v2/CONFIG.json';cfg=json.loads(cfgpath.read_text())
    names=('implementation','protocol','metadata_helper','kernel','combiner','native_binary','build','bridge_ready')
    methods={k:cfg[k] for k in names}
    for p in methods.values():verify(p)
    assert methods['implementation']['sha256']=='6737906cbc088d68068bc35c19bb9d6a4423b00481ca5390e2b5a7f3b5f70668'
    assert methods['native_binary']['sha256']=='68a47932f106d7d40bedaac1ddc034f6a314df182f8db07f318b1c9742f05e6a'
    assert read(methods['bridge_ready'])['root_build']==methods['build']
    assert read(methods['build'])['binary']==methods['native_binary']
    package=Path(__file__).resolve().parent
    for name in ('PREPARE.py','RUN.py'):ast.parse((package/name).read_text())
    a.out.mkdir(parents=True,exist_ok=False);states={};evidence={}
    key=verify(base['key']).read_bytes();prm=verify(base['parameter_file']).read_bytes()
    for id in IDS:
        r=old['rows'][id];g=r['geometry'];c=r['coordinate'];assert r['baseline_compatibility']['compatible']
        assert r['inputs']['key']['sha256']==f['inputs']['key']['sha256'] and r['inputs']['parameters']['sha256']==base['parameter_file']['sha256']
        assert all(g[k]==f['geometry'][k] for k in ('atomic_numbers','atom_ids','metal_atom','ligands','tinker_types','bonds_zero_based'))
        xyz=verify(r['inputs']['xyz']).read_bytes();lines=[s.split() for s in xyz.decode().splitlines() if s.strip()]
        assert int(lines[0][0])==8 and [[float(s) for s in row[2:5]] for row in lines[1:]]==g['xyz_A']
        for p in r['inputs'].values():verify(p)
        verify(r['reference']['source']);verify(c['source_geometry'])
        folder=a.out/id;folder.mkdir();(folder/'system.xyz').write_bytes(xyz);(folder/'system.key').write_bytes(key);(folder/'amoeba09.prm').write_bytes(prm)
        states[id]={'id':id,'case_id':f['id'],'coefficient_mode':'zero','direction':c['direction'],
                    'step_tag':c['shift_label'].split('_')[1],'step_A':abs(c['q_A']),'sign':-1 if c['q_A']<0 else 1,
                    'fixture':center,'geometry':g,'center_geometry':f['geometry'],'parameters':{'pairs':{},'triples':[]},
                    'inputs':{'xyz':pin(folder/'system.xyz'),'key':pin(folder/'system.key'),'parameters':pin(folder/'amoeba09.prm')}}
        evidence[id]={'coordinate':c,'source_inputs':r['inputs'],'baseline_compatibility':r['baseline_compatibility'],
                      'historical_model':r['model'],'original_quantum_target':r['reference'],'consumed_development':True,
                      'original_endpoint_operator':None,'original_endpoint_gradient':None}
    result={'status':'prepared_eight_existing_La_donor_baseline_exports','implementation':pin(__file__),
            'runtime_source':pin(package/'RUN.py'),'slurm_wrapper':pin(package/'RUN.sbatch'),
            'development_plan':pin(D/'water_data_v1/DEVELOPMENT_PLAN.txt'),'source_endpoints':source,
            'shared_baseline':baseline,'center_fixture':center,'methods':methods,'parent_runtime_configuration':pin(cfgpath),
            'phase':read(methods['bridge_ready'])['phase'],'units':cfg['units'],'timeout_seconds':300,
            'ordered_ids':list(IDS),'states':states,'evidence':evidence,
            'execution':{'qos':'standard_expand','workers':4,'allocated_cpus':8,'threads_each_controller':1,'threads_each_native':1,
                         'maximum_native_programs':8,'maximum_dense_baseline_solves':8,'maximum_mutual_operator_applications':192,
                         'new_centers':0,'new_quantum_calls':0,'new_fit_calls':0,'root_release':None},
            'scope':'Fresh baseline exports only at eight literal saved geometries. No new quantum targets, coefficients, center or force-regression admission.',
            'reused_runtime_function':'RUN.run_state only; its 102-state main/gate is not invoked',
            'preparation_cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'numerical_calls':0}}
    save(a.out/'MANIFEST.json',result);print(json.dumps(pin(a.out/'MANIFEST.json')))

if __name__=='__main__':main()
