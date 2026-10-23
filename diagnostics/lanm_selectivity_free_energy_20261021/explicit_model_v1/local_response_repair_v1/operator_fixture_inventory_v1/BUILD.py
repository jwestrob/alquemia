#!/usr/bin/env python3
"""One-shot saved metadata inventory; no native executor or parameter search."""
import hashlib, json, os, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
D = ROOT / 'diagnostics/lanm_selectivity_free_energy_20261021'
W = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021'
OUT = Path(__file__).resolve().parent
START, CPU = time.perf_counter(), time.process_time()
FILES, CACHE = {}, {}

def digest(b): return hashlib.sha256(b).hexdigest()
def pin(p):
    p = Path(p).resolve()
    if str(p) not in CACHE: CACHE[str(p)] = dict(path=str(p), sha256=digest(p.read_bytes()))
    return CACHE[str(p)]
def register(p, label=None):
    q = p if isinstance(p, dict) else pin(p)
    assert pin(q['path']) == q, q['path']
    key = label or 'file_' + q['sha256'][:16]
    # Identical bytes at different historical paths remain distinguishable.
    if key in FILES and FILES[key] != q: key += '_' + digest(q['path'].encode())[:8]
    FILES[key] = q
    return key
def read(p):
    if isinstance(p, dict): register(p); p = p['path']
    return json.loads(Path(p).read_text())
def inputs(q): return {k: register(q[k]) for k in ('xyz','key','parameters')}
def geometry(q):
    register(q['xyz'])
    lines = Path(q['xyz']['path']).read_text().splitlines()
    n = int(lines[0].split()[0]); atoms = []
    for i, line in enumerate(lines[1:n+1], 1):
        s = line.split(); assert int(s[0]) == i
        atoms.append(dict(type=int(s[5]), xyz_A=list(map(float,s[2:5])), bonds=sorted(map(int,s[6:]))))
    assert len(atoms) == n
    return atoms
def key_lines(q, field=None):
    register(q['key']); lines = []
    removed = []
    for raw in Path(q['key']['path']).read_text().splitlines():
        s = raw.split('#',1)[0].split()
        if not s: continue
        s[0] = s[0].lower()
        if s[0] in ('polar-print','save-uinduce'):
            assert len(s)==1; removed.append(s[0]); continue
        if s[0]=='external-field':
            assert field is not None and list(map(float,s[1:]))==field
            removed.append('external-field'); continue
        lines.append(s)
    if field is not None: assert removed.count('external-field')==1
    return lines
def signature(q, field=None):
    register(q['parameters'])
    return digest(json.dumps([geometry(q),key_lines(q,field),q['parameters']['sha256']],sort_keys=True).encode())
def pointer(artifact, path): return dict(artifact=artifact, json_path=path)

target_path = D/'explicit_model_v1/local_response_repair_v1/targets/TARGET_INVENTORY.json'
target = read(target_path); register(target_path,'target_inventory')
data = read(target['artifacts']['development_data']); register(target['artifacts']['development_data'],'development_data')
mr = read(target['artifacts']['model_moments']); register(target['artifacts']['model_moments'],'model_moments')
mm = read(mr['manifest']); register(mr['manifest'],'moment_manifest')
tr_path = W/'explicit_model_v1/collective_mechanics_options_v1/run_v1/model/RESULT.json'
tr = read(tr_path); register(tr_path,'trimer_model'); tm=read(tr['manifest']); register(tr['manifest'],'trimer_manifest')
fr_path = W/'explicit_model_v1/conditional_susceptibility_v1/model/run_v1/execution/RESULT.json'
fr = read(fr_path); register(fr_path,'field_model'); fm=read(fr['manifest']); register(fr['manifest'],'field_manifest')
tc=read(target['artifacts']['trimer_comparison']);register(target['artifacts']['trimer_comparison'],'trimer_comparison')
for k in ('development_objective','body_helper','development_selected_review','dipole_comparison','dipole_arithmetic_review','field_comparison'):
    register(target['artifacts'][k],k)
for k in ('moment_review','trimer_review','analyze','native_units'):
    register(fm['sources'][k],k)
for k,q in mr['config']['physics_sources'].items(): register(q,'baseline_'+k)
assert data['provenance']['candidate_result']==mr['candidate_result']==mm['candidate_result']==tm['candidate_result']
assert data['baseline_parameters']==mm['selected_parameters']==tm['selected_parameters']==fm['selected_parameters']
assert fr['error'] is None and tr['error'] is None and mr['error'] is None

ROWS = {}
def new(cid, category, q, parent, atom_map):
    assert cid not in ROWS
    atoms=geometry(q); assert sum(a['type'] in (500,501) for a in atoms)==1
    r=dict(category=category,metal='La' if atoms[0]['type']==500 else 'Dy', atoms=len(atoms),
           input=inputs(q), geometry_types_bonds_sha256=digest(json.dumps(atoms,sort_keys=True).encode()),
           effective_input_sha256=signature(q),parent=parent,parent_atom_indices=atom_map,
           cached=[], aliases=[], matrix_export='requires_new_instrumented_call')
    ROWS[cid]=r; return r
def moments(r, artifact, cid, saved, q, field='moments'):
    assert saved['status']=='complete' and signature(saved['input'])==signature(q)
    z=saved[field]
    r['cached'].append(dict(**pointer(artifact,['rows',cid,field]),energy_kcal_mol=z['energy_kcal_mol'],
        induced_Debye=pointer(artifact,['rows',cid,field,'induced_Debye']),
        energy_output=register(z['analyze_output']),induced_output=register(z['induced_output'])))

for cid, model in data['baseline_model'].items():
    q=model['provenance']['effective_input']; g=data['rows'][cid]['geometry']
    r=new(cid,'original_development',q,cid,list(range(g['atoms'])))
    assert [a['xyz_A'] for a in geometry(q)]==g['positions_A']
    r['cached'].append(dict(**pointer('development_data',['baseline_model',cid]),
        energy_kcal_mol=model['energy'],forces_available=True,
        native_authority=register(model['provenance']['native_authority'])))
    r['uses']=['original_25_energy_force_controls']
    if cid in mm['rows']:
        z=mm['rows'][cid]; assert signature(z['input'])==signature(q)
        r['aliases'].append(pointer('moment_manifest',['rows',cid]))
        moments(r,'model_moments',cid,mr['rows'][cid],z['input'])
        r['uses'].append('parent_C_P_B_dipoles')

for cid, z in mm['rows'].items():
    if z['kind'] not in ('dimer','ion'): continue
    q=z['input']; parent=z['parent_geometry_id']; r=new(cid,'literal_pair_dimer' if z['kind']=='dimer' else 'bare_ion_control',q,parent,z['parent_atom_indices'])
    src=geometry(z['source_input']); actual=geometry(q); mapping=z['parent_atom_indices']
    assert [[a['type'],a['xyz_A']] for a in actual]==[[src[i]['type'],src[i]['xyz_A']] for i in mapping]
    assert key_lines(q)==key_lines(z['source_input'])
    r['source_manifest']=pointer('moment_manifest',['rows',cid]);r['uses']=['pair_C_P_B_dipoles','trimer_cycles_as_applicable']
    moments(r,'model_moments',cid,mr['rows'][cid],q)

for cid,z in tm['model_rows'].items():
    if z['kind']!='metal_two_water': continue
    q=z['input'];r=new(cid,'near_far_trimer',q,z['parent_geometry_id'],z['parent_atom_indices'])
    src=geometry(z['source_input']); assert [[a['type'],a['xyz_A']] for a in geometry(q)]==[[src[i]['type'],src[i]['xyz_A']] for i in z['parent_atom_indices']]
    assert key_lines(q)==key_lines(z['source_input'])
    r['source_manifest']=pointer('trimer_manifest',['model_rows',cid]);r['uses']=['E3_mu3']
    moments(r,'trimer_model',cid,tr['rows'][cid],q,'parsed_moments')

assert len(ROWS)==54 and len({v['effective_input_sha256'] for v in ROWS.values()})==54
assert len({json.dumps(key_lines({k:FILES[v] for k,v in r['input'].items()})) for r in ROWS.values()})==2
FIELDS={}
for state,z in fm['rows'].items():
    cid=z['id'];assert cid in ROWS and signature(z['input'])==ROWS[cid]['effective_input_sha256']
    if state=='ion_La': continue
    ROWS[cid]['aliases'].append(pointer('field_manifest',['rows',state])); ROWS[cid]['uses'].append('conditional_susceptibility')
    jobs={}
    for jid,j in fm['jobs'].items():
        if j['state']!=state: continue
        assert signature(j['input'],j['field_MV_cm'])==ROWS[cid]['effective_input_sha256']
        actual=fr['rows'][jid];assert actual['status']=='complete' and actual['job']==j
        jobs[jid]=dict(field_au=j['field_au'],input=inputs(j['input']),cached=pointer('field_model',['rows',jid,'moments']))
    FIELDS[state]=dict(fixture=cid,zero=pointer('field_manifest',['rows',state,'baseline']),signed_field_targets=jobs,
                       summary=pointer('field_model',['summary','states',state]))
assert len(FIELDS)==4 and sum(len(x['signed_field_targets']) for x in FIELDS.values())==17 # 16 fields plus prior zero identity

CACHED={}
for artifact, manifest, result, mf in [('model_moments',mm['rows'],mr,'moments'),('trimer_model',tm['model_rows'],tr,'parsed_moments')]:
    for cid,z in manifest.items():
        if any(a['type'] in (500,501) for a in geometry(z['input'])): continue
        assert result['rows'][cid]['status']=='complete'
        CACHED[cid]=dict(source=pointer(artifact,['rows',cid,mf]),input=inputs(z['input']),new_calls=0)
assert len(CACHED)==46
for cycle,z in tc['cycles'].items(): assert set(z['model_terms']) <= set(ROWS)|set(CACHED)

result=dict(status='source_inventory_complete_unreleased',implementation=pin(__file__),artifacts=FILES,
 candidate_result=data['provenance']['candidate_result'],candidate_pin_authority='development_data provenance and selected raw review; no fit trajectory reread',
 selected_parameters=data['baseline_parameters'],counts=dict(original_development=10,literal_pair_dimers=34,near_far_trimers=8,ion_controls=2,total_baseline_exports=54,
 field_state_aliases=4,cached_signed_field_rows=16,cached_field_zero_identity=1,metal_free_cached_rows=46,new_QM_calls=0,new_native_calls=0),
 first_gate='aquo9_central_La__dimer_01',remaining_after_actual_parity=[k for k in ROWS if k!='aquo9_central_La__dimer_01'],rows=ROWS,
 field_aliases=FIELDS,cached_metal_absent=CACHED,trimer_cycle_joins={k:v['model_terms'] for k,v in tc['cycles'].items()},
 scope=dict(selection='complete pre-existing consumed geometry inventory; no outcome-based selection or symmetry reduction',
  deduplication='Exact Cartesian coordinates, atom order/types/bonds, parameter bytes and all active key tokens; comments/whitespace and only polar-print/save-uinduce are ignored. Signed fields are separate cached targets, not extra zero-field exports.',
  instrumentation='Only RESPONSE-MATRIX-EXPORT reporting keyword may be appended; reviewed source/binary and actual parity are future root-owned gates.',
  pending='A, Ed, Ep and independent u/v arrays require new instrumented calls; old .uind is a baseline u target, not proof of separate v or a stored operator.',
  calls='54 prospective instrumented analyze EM calls, first1+remaining53; two snapshots per call are not separate native invocations. Any fresh uninstrumented parity call is additional and owned by the scout protocol.',
  ions='Kept as controls. Proposed oxygen-only correction vanishes on bare ions; cached zero-field ion energies/moments need no candidate reevaluation.',
  forces='A fixed-coordinate operator does not provide coordinate derivatives; original force gates cannot be certified from these exports alone.',
  reserved='All reserved/held and earlier orientation-validation fixtures excluded exactly as frozen target inventory. No new loss, thresholds, fit or molecular execution.'),
 cost=dict(wall_seconds=time.perf_counter()-START,cpu_seconds=time.process_time()-CPU,local_threads=1,new_native_calls=0,verified_file_pins=len(CACHE)))
destination=OUT/'FIXTURE_INVENTORY.json'
assert not destination.exists()
destination.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(inventory=pin(destination),counts=result['counts'],cost=result['cost'])))
