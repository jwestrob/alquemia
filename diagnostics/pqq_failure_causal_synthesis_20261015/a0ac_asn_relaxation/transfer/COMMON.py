"""Metadata utilities for the finite transfer; no scientific imports."""
import copy, hashlib, importlib.util, json, sys
from pathlib import Path
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/mace_omol.py').is_file())
WORK=ROOT/'workspaces/pqq_failure_causal_synthesis_20261015'
HERE=ROOT/'diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer'
METALS=('Ca','La');ORDER3=('origin','adaptive_Ca','adaptive_La');ORDER4=ORDER3+('amide_flip',);ORDER6=ORDER4+('relaxed_flip_Ca','relaxed_flip_La')
ORDERS={'common3':list(ORDER3),'common4':list(ORDER4),'common6':list(ORDER6)}
INVENTORY_SHA='83bc0b0373c4453a71f87eccf98675fe7d6f4b3a1a6ff3708b24ec1a2a75c6c0'
def require(v,why):
    if not v:raise ValueError(why)
def pin(p):
    p=Path(p).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return {'path':str(p),'sha256':h.hexdigest()}
def check(p):require(pin(p['path'])==p,'Changed artifact: '+p['path']);return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def load(name,p):
    s=importlib.util.spec_from_file_location(name,check(p));m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def canonical_gate(collection):
    cp=pin(collection) if not isinstance(collection,dict) else collection;r=read(cp);c=read(r['config']);ep=pin(Path(cp['path']).with_name('EXECUTION.json'));e=read(ep)
    require(e['status']=='complete' and e['config']==r['config'] and r['all25_common6_available'],'Actual complete canonical solvent gate required')
    inv=read(c['calibration_inventory']);require([x['record_id'] for x in r['rows']]==inv['calibration_rule']['designated_record_ids'],'Designated25 membership changed')
    require(r['pool_orders']==ORDERS and all(r['calibrations']['common6'][p]['status']=='available' and r['calibrations']['common6'][p]['available']==25 for p in ('operational','mathematical')),'Frozen common6 canonical calibration gate unavailable')
    require(r['calibrations']['common3']==r['old_calibrations']==read(c['canonical_collection'])['calibrations'],'Original common3 calibration changed')
    return {'collection':cp,'config':r['config'],'execution':ep,'calibrations':r['calibrations'],'old_calibrations':r['old_calibrations']}
def context_key(r):
    value=[r['source_case_id'],r['selection_id'],r['source_preparation']['sha256'],r['state_signature']['sha256'],[r['maps'][z]['sha256'] for z in METALS]]
    return 'context_'+hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()[:24]
def cohort_ids(ledger):
    return {'expanded135':ledger['cohorts']['expanded135'],'current_ca125':[r['record_id'] for r in ledger['cohorts']['current_ca125_annotated_overlay']],'two_ca15':ledger['cohorts']['two_ca15']}
def scope_record(r,entry,artifacts):
    route=entry['source_route'];sm_pin=artifacts[route['original_search_manifest']];sm=read(sm_pin)
    tasks=[sm['tasks'][t['manifest_task_index']] for t in route['tasks']]
    require(len(tasks)==2 and all(t['task_id']==v['task_id'] and t['metal']==v['metal'] for t,v in zip(tasks,route['tasks'])),'Original task index differs')
    parent=artifacts[route['core_parent']];core=read(parent);proton=artifacts[route['protonation']];source=artifacts[route['protonated_source']]
    require(read(proton)['output']==source==read(r['source_preparation'])['source'],'Exact current source route differs')
    role=core['fixed_core']['requested_roles'].get('anchor_asparagine')
    require(role==route['original_role_selector'],'Original role metadata differs')
    normalized=dict(chain=role['chain'],resnum=int(role['resnum']),insertion_code=role.get('icode',role.get('insertion_code','')),resname=role['resname']) if isinstance(role,dict) else None
    bc=read(r['baseline_manifest'])['cases'][r['baseline_manifest_case_index']]
    require(bc['case_id']==r['case_id'],'Original baseline case differs')
    if 'cases' in sm:
        sc,=[v for v in sm['cases'] if v['case_id']==r['case_id']];active=sc['union']['source']['metal'];active_evidence={'manifest':sm_pin,'case_id':sc['case_id'],'field':'union.source.metal'}
        require(sc['union']['original_core']['parent']==parent,'Original core lineage differs')
    else:
        active=bc['native_source_record']['source']['origin_row']['source']['source']['metal'];active_evidence={'manifest':r['baseline_manifest'],'case_index':r['baseline_manifest_case_index'],'field':'native_source_record.source.origin_row.source.source.metal'}
    spectator=bc.get('new_source_admission',{}).get('spectator_admission')
    if spectator:
        sp=read(spectator);require(sp['active_selector']==active and sp['full_protonated_source']==source and sp['status']=='passed' and sp['passed'],'Exact spectator lineage differs')
    require(bool(spectator)==('two_ca15' in entry['membership']),'Two-Ca spectator admission missing or unexpected')
    result={k:r[k] for k in ('record_id','case_id','source_case_id','selection_id','protein','maps','state_signature','source_preparation')}
    result.update(core_parent=parent,protonation=proton,protonated_source=source,search_manifest=sm_pin,original_tasks=tasks,
        role_name='anchor_asparagine',original_role_selector=role,role=normalized,included_role='anchor_asparagine' in core['fixed_core']['included_roles'],
        original_xyz={z:r['matrix'][z]['origin']['xyz'] for z in METALS},six_saved_cells=r['matrix'],active_selector=active,active_selector_evidence=active_evidence,spectator_admission=spectator)
    return result
def inventories():
    ip=pin(HERE.parent/'transfer_inventory/INVENTORY_PLAN.json');require(ip['sha256']==INVENTORY_SHA,'Frozen transfer ledger differs');ledger=read(ip);a=ledger['artifacts'];invp=a[ledger['authority']['inventory']];inv=read(invp)
    cohorts=cohort_ids(ledger);current=set(cohorts['expanded135'])|set(cohorts['current_ca125']);entries={r['record_id']:r for r in ledger['records']};contexts={}
    for r in inv['records']:
        if r['record_id'] not in current:continue
        e=entries[r['record_id']];key=context_key(r)
        if key in contexts:
            prior=contexts[key]['record'];require(all(r[k]==prior[k] for k in ('case_id','source_case_id','selection_id','maps','state_signature','source_preparation','matrix')),'Duplicate context evidence differs')
            contexts[key]['memberships'].append(r['record_id']);continue
        contexts[key]={'context_id':key,'record':r,'memberships':[r['record_id']],'scope':scope_record(r,e,a)}
    memberships=[rid for s in contexts.values() for rid in s['memberships']]
    require(len(current)==len(memberships)==len(set(memberships))==260 and set(memberships)==current and len(contexts)==257,'Current membership/context denominator differs')
    return ip,ledger,invp,inv,list(contexts.values())

def cache_registry(contexts,gate):
    """Metadata aliases only. All producing receipts and exact common3 remain visible."""
    canonical=read(WORK/'asn_canonical_search_GPU_v1/execution/FINAL_COLLECTION.json');cp=canonical['config'];config=read(cp)
    prep_pin=config['preparation'];prep=read(prep_pin);canonical_prep={r['source_case_id']:r for r in prep['rows']}
    a0ac_pin=pin(WORK/'asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json');a0ac=read(a0ac_pin);ac=read(a0ac['config']);ap=read(ac['preparation']);a0ac_prep={r['source_case_id']:r for r in ap['rows']}
    canonical_pin=pin(WORK/'asn_canonical_search_GPU_v1/execution/FINAL_COLLECTION.json')
    canonical_native=read(gate['collection']);native_config=read(gate['config'])
    require(native_config['score_collection']==canonical_pin and native_config['score_config']==cp and canonical_native['config']==gate['config'],'Canonical native/MACE producer lineage differs')
    a0ac_native_pin=pin(WORK/'asn_relaxation_solvent_v1/execution/FINAL_COLLECTION.json');a0ac_native=read(a0ac_native_pin);anc=read(a0ac_native['config']);ane=read(Path(a0ac_native_pin['path']).with_name('EXECUTION.json'))
    require(anc['score_collection']==a0ac_pin and anc['score_config']==a0ac['config'] and ane['config']==a0ac_native['config'] and ane['status']=='complete' and a0ac_native['available_new_cells']==40,'A0AC native/MACE producer lineage differs')
    sources=[(canonical,canonical_pin,prep_pin,canonical_prep,gate['collection'],canonical_native),
             (a0ac,a0ac_pin,ac['preparation'],a0ac_prep,a0ac_native_pin,a0ac_native)]
    for coll,collpin,prepared_pin,prepared,native_pin,native in sources:
        ex=pin(Path(collpin['path']).with_name('EXECUTION.json'));require(read(ex)['status']=='complete' and read(ex)['config']==coll['config'],'Reuse producer incomplete')
        for row in coll['rows']:
            sid=row['source_case_id'];matches=[s for s in contexts if s['record']['source_case_id']==sid and all(s['record']['matrix'][z][q]['MACE']==row['matrix'][z][q]['receipt'] and s['record']['matrix'][z][q]['xyz']==row['matrix'][z][q]['xyz'] for z in METALS for q in ORDER3)]
            if not matches:continue
            require(len(matches)==1 and row['status']=='complete','Ambiguous/incomplete exact executed alias')
            s=matches[0];old=prepared[sid]
            require(all(old['original_row'][k]==s['record'][k] for k in ('maps','state_signature','source_preparation')),'Reused physical source/state differs')
            nr,=[r for r in native['rows'] if r['source_case_id']==sid]
            require(all(nr['evidence'][z][q]['status']=='available' and nr['evidence'][z][q]['MACE']==row['matrix'][z][q]['receipt'] and nr['evidence'][z][q]['xyz']==row['matrix'][z][q]['xyz'] for z in METALS for q in ORDER6),'Native/MACE exact endpoint lineage differs')
            if 'reuse' in s:
                require(all(s['reuse']['executed_result_alias']['matrix'][z][q]['receipt']==row['matrix'][z][q]['receipt'] for z in METALS for q in ORDER6),'Duplicate reuse differs');continue
            s['reuse']={'prepared_row':old,'executed_result_alias':row,'producing_config':coll['config'],'producing_collection':collpin,'producing_execution':ex,'producing_preparation':prepared_pin,
                        'solvent_collection':native_pin,'solvent_row':nr,'repeat_search_permitted':False}
    require(sum('reuse' in s for s in contexts)==34,'Exact25 plus A0AC10 shared-one reuse differs')
    return contexts
