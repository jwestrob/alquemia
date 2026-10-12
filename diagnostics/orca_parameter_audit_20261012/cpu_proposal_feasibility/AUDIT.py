"""File-only CPU proposal feasibility inventory; never imports a model/optimizer."""
import email, hashlib, json, statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
PINS={}
def pin(path):
    p=Path(path).resolve();d={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};PINS[str(p)]=d;return d
def read(path):
    if isinstance(path,dict):assert pin(path['path'])['sha256']==path['sha256'];path=path['path']
    else:pin(path)
    return json.loads(Path(path).read_text())
def packages(root):
    result={};site=Path(root)/'lib/python3.11/site-packages'
    for name in ['scipy','numpy','torch','mace_torch','ase','e3nn']:
        paths=list(site.glob(name+'-*.dist-info/METADATA'))
        if paths:
            assert len(paths)==1;d=email.message_from_string(paths[0].read_text())
            result[name]={'version':d.get('Version'),'metadata':pin(paths[0])}
    return result
archive=ROOT/'workspaces/motion_envelope_20260923/pilot34_searches_v1'
m=read(archive/'manifest.json')
source_ids=['a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1','q9z4j7-pqq-la_model__conditioned_La__seed-1_sample-0']
fixtures=[];costs=[]
for t in m['tasks']:
    path=archive/'proposals'/t['task_id']/'result.json'
    # This is a finite 68-receipt accounting pass, not an archive crawl.
    r=json.loads(path.read_text())
    calls=len({q['result']['path'] for q in r['requests'] if not q.get('reused')})
    costs.append({'task_id':t['task_id'],'status':r['status'],'new_MACE_calls':calls,
                  'optimizer':r['optimizer'],'wall_seconds':r['wall_seconds'],'boundary_flag':r.get('boundary_flag')})
    if t['source_case_id'] not in source_ids:continue
    pin(path)
    points={}
    for name,point in [('origin',t['origin_reuse']['point']),('proposal',r['proposal'])]:
        for key in ['coordinate','MACE','forces']:
            assert pin(point[key]['path'])['sha256']==point[key]['sha256']
        native=read(point['MACE'])
        actual_model=native.get('model') or read(native['manifest'])['model']
        assert actual_model==m['model']
        assert native['energy_eV']==point['MACE_eV'] and native['forces']==point['forces']
        if 'request' in native:
            request=read(native['request']);assert request['xyz']==point['coordinate']
            assert (request['charge'],request['multiplicity'])==(t['charge'],t['multiplicity'])
        else:
            nm=read(native['manifest']);nt=next(x for x in nm['tasks'] if x['task_id']==native['task_id'])
            assert nt['xyz']==point['coordinate'] and nt['charge']==t['charge']
            assert nt.get('spin_multiplicity',nt.get('multiplicity'))==t['multiplicity']
        points[name]=point
    read(t['mapping'])
    fixtures.append({'task_id':t['task_id'],'source_case_id':t['source_case_id'],'metal':t['metal'],
      'charge':t['charge'],'multiplicity':t['multiplicity'],'atom_count':int(Path(t['xyz']['path']).read_text().splitlines()[0]),
      'mapping':t['mapping'],'active_mode_ids':t['active_mode_ids'],'selector':t['selector'],
      'proposal_receipt':pin(path),'points':points,'historical_cost':costs[-1],
      'final_geometry':r['final_geometry'],'residual':r['residual']})
assert len(fixtures)==4
implementation={}
for name in ['motion_envelope_run.py','slsqp_precision.py','adaptive_completion.py','adaptive_angular_proposals.py',
             'accommodation_proposals.py','accommodation_nonlinear.py','mace_site_kinematics.py','adaptive_force_diagnostic.py',
             'union_adaptive.py','mace_omol.py','context_factorial_mace.py']:
    p=pin(ROOT/'scripts'/name);implementation[name]={'current':p,'historical_snapshot':m['implementation'].get(name),
       'snapshot_bytes_equal':p['sha256']==m['implementation'].get(name,{}).get('sha256')}
parent=Path(m['cpu_python']).parent.parent;worker=Path(m['gpu_python']).parent.parent
optimizer={}
for name in ['wrapper','kernel']:
    p=m['optimizer_software'][name];assert pin(p['path'])['sha256']==p['sha256']
    other=worker/'lib/python3.11/site-packages/scipy/optimize'/Path(p['path']).name
    optimizer[name]={'parent':p,'model_environment':pin(other),'bytes_equal':p['sha256']==pin(other)['sha256']}
energy=read(ROOT/'workspaces/orca_parameter_audit_20261012/context_mace_v1/FINAL.json')
gpu=read(archive/'gpu_summary_1211342.json');execution=read(archive/'EXECUTION_1211342.json')
assert sum(r['new_MACE_calls'] for r in costs)==gpu['new_MACE_calls']==998
result={'status':'feasible_isolated_adapter_required_unexecuted','implementation':pin(__file__),
 'historical_manifest':pin(archive/'manifest.json'),'exact_settings':m['settings'],'model':m['model'],
 'environments':{'optimizer_parent':{'invocation':m['cpu_python'],'packages':packages(parent),'optimizer_version':m['optimizer_software']['version']},
                 'persistent_CPU_worker':{'invocation':m['gpu_python'],'packages':packages(worker),'venv_configuration':pin(worker/'pyvenv.cfg')}},
 'optimizer_binary_comparison':optimizer,'source_implementation':implementation,'force_proposal_fixtures':fixtures,
 'historical_68_proposal_cost':{'new_MACE_calls':998,'range':[min(r['new_MACE_calls'] for r in costs),max(r['new_MACE_calls'] for r in costs)],
      'median_calls':statistics.median(r['new_MACE_calls'] for r in costs),'GPU_summary':pin(archive/'gpu_summary_1211342.json'),
      'execution_receipt':pin(archive/'EXECUTION_1211342.json'),'GPU_elapsed_seconds':gpu['wall_seconds'],'outer_elapsed_seconds':execution['wall_seconds'],
      'rows':costs,'CPU_force_or_optimization_timing_measured':False},
 'energy_only_CPU_evidence':{'receipt':pin(ROOT/'workspaces/orca_parameter_audit_20261012/context_mace_v1/FINAL.json'),
      'new_calls':energy['actual_new_calls'],'model_load_seconds':energy['model_load_seconds'],'elapsed_seconds':energy['wall_seconds'],
      'CPU_seconds':energy['process_cpu_seconds'],'peak_RSS_KiB':energy['peak_host_RSS_KiB'],
      'qualification_scope':'Four GPU-energy origin replays and16 cross-context energies; force/proposal parity not tested.'},
 'minimal_scout_recommendation':{'source_case_id':source_ids[0],'phase1_energy_force_pairs':5,
      'phase1_states':'Ca origin, La origin, Ca saved proposal, La saved proposal, repeat Ca origin after state/geometry switches.',
      'phase2_optimizer_starts':2,'historical_nonzero_model_calls_for_same_two_starts':28,
      'new_call_count_if_same_optimizer_trace':33,'hard_guaranteed_call_count':False,
      'reason':'SLSQP maxiter200 limits iterations, not a fixed nfev; do not silently add maxfun or claim an exact new count.',
      'outer_resource_stop':'If needed, declare a separate task call/wall cap before execution; exhausted cap means unavailable, never converged.',
      'proposed_runtime_allocation':'Original optimizer parent one BLAS thread + persistent CPU worker eight Torch threads, interop1, other BLAS1; charge conservative9 of aggregate16.',
      'gates':{'energy_replay_kcal_mol':0.01,'maximum_Cartesian_force_error_eV_A':0.001,'proposal_energy_difference_kcal_mol':0.001,
        'selector':'Same ordered selected mode ids/criteria/redundancies; numerical diagnostic values separately reported.',
        'proposal':'Existing optimizer success, energy nonincrease, final physical geometry, boundary flags; report actual full-q and coordinate differences without declaring them bit-identical.'},
      'new_DFT_or_solvent_calls':0,'status':'Recommendation only; no scout launched'},
 'new_molecular_calls':0,'new_optimizer_calls':0,'new_model_imports':0,'new_dependencies':0,'pins':list(PINS.values())}
with (OUT/'RESULT.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'fixtures':len(fixtures),'archived_proposal_calls':998,'optimizer_kernel_matches_model_environment':optimizer['kernel']['bytes_equal'],'molecular_calls':0}))
