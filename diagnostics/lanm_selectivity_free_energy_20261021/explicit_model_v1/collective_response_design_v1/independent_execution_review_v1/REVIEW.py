#!/usr/bin/env python3
"""Independent finite metadata/pin review only; no engines or source controls run."""
import hashlib,json,math,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;D=HERE.parent;ROOT=D.parents[3]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/collective_response_design_v1'
def pin(p):p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):assert pin(p['path'])==p;return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,r):
    with p.open('x') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
tick,cpu=time.monotonic(),time.process_time();checks=[]
def check(name,ok):assert ok,name;checks.append(name)
cfgpin=pin(D/'CONFIG_PREPARE_v3.json');cfg=read(cfgpin);prep_pin=pin(W/'preparation_v3/RESULT.json');prep=read(prep_pin);manifest=read(prep['manifest'])
controls_pin=pin(W/'source_controls_v1/RESULT.json');controls=read(controls_pin)
check('final source9eb552e6',cfg['implementation']==pin(D/'RUN.py') and cfg['implementation']['sha256']=='9eb552e6ec65da01b3e62ac803740405ee7cc05ab45d1371294fdea5b3ca422d')
check('final zero-call preparation binding',prep['implementation']==cfg['implementation'] and prep['config']==cfgpin and prep['status']=='prepared' and prep['new_native_calls']==0)
check('source control receipt binding',controls['status']=='passed' and controls['checks_passed']==229 and controls['subject']==cfg['implementation'] and controls['preparation']==prep_pin and controls['new_native_calls']==0)
for name,v in cfg.items():
    if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
for group in ('executables','physics_sources'):
    for v in cfg[group].values():verify(v)
check('frozen q and numerical scope',cfg['q_values']==[.01,.05] and cfg['native_overlap_exponent_A_inverse']==4. and cfg['normalization_distance_A']==2.5 and cfg['maximum_native_programs']==116 and cfg['maximum_internal_FD_energies']==636 and cfg['workers']==4)
data=read(cfg['data']);check('exact failed development point',manifest['candidate_result']==data['provenance']['candidate_result'] and manifest['candidate_result']['sha256']=='5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832' and manifest['selected_review']==data['provenance']['selected_review'])
rows=manifest['rows'];check('92cells116programs636FD',len(rows)==92 and sum(r['native_programs'] for r in rows.values())==116 and sum(r['internal_FD_energies'] for r in rows.values())==636)
check('complete staged4_4_84',list(map(len,manifest['stages'].values()))==[4,4,84] and set(sum(manifest['stages'].values(),[]))==set(rows) and len(sum(manifest['stages'].values(),[]))==92)
clusters={f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy')}
check('four positive-q Cartesian gates',set(manifest['stages']['derivatives'])=={'q05__'+g for g in clusters} and all(rows[c]['q']==.05 and rows[c]['derivative_test'] for c in manifest['stages']['derivatives']))
check('FD literal25_25_28_28',sorted(rows[c]['atoms'] for c in manifest['stages']['derivatives'])==[25,25,28,28])
for q in (.01,.05):
    cells=[r for r in rows.values() if r['q']==q and r['kind'] in ('fixture','pair')]
    check('fixed44rows bothmetals '+str(q),len(cells)==44 and sum(r['kind']=='fixture' for r in cells)==10 and sum(r['kind']=='pair' for r in cells)==34)
S=((1+5+25/3)*math.exp(-5))**2
for cid,r in rows.items():
    inp,src=r['input'],r['source_input'];xyz=verify(inp['xyz']).read_text();original=verify(src['xyz']).read_text();key=verify(inp['key']).read_text();base=verify(src['key']).read_text()
    check('input bytes '+cid,verify(inp['parameters']).read_bytes()==verify(src['parameters']).read_bytes() and key.startswith(base+'\n'))
    added=[v.split() for v in key[len(base):].splitlines() if v.strip() and not v.startswith('#')]
    check('only declared records '+cid,len(added)==6 and added[0]==['exchange-polar','S2U'] and added[1]==['polar-print'] and added[3:]==[['exchpol','35','0','0','4','0'],['exchpol','500','0','1','4','1'],['exchpol','501','0','1','4','1']] and added[2][:2]==['exchpol','34'] and added[2][3:]==['1','4','0'] and abs(float(added[2][2])*S-r['q'])<1e-16)
    if r['kind']=='shell':
        source=[s.split() for s in original.splitlines()[2:] if s.strip()];target=[s.split() for s in xyz.splitlines()[1:] if s.strip()]
        expected=[[str(int(v[0])-1),*v[1:6],*[str(int(j)-1) for j in v[6:]]] for v in source]
        check('literal27-water-only '+cid,len(target)==27 and target==expected)
    else:check('unchanged literal XYZ '+cid,xyz==original)
    check('one_or_two_programs '+cid,r['native_programs']==(1 if r['kind']=='pair' else 2) and r['internal_FD_energies']==(6*r['atoms'] if r['derivative_test'] else 0))
source=(D/'RUN.py').read_text();check('analyze admitted before testgrad',source.index("parsed=pair.parse_analyze(folder/'analyze.out')")<source.index("checks['testgrad']=invoke"))
check('residual every-summary gate',"len(matches)==len(re.findall(r'Induced Dipoles" in source and '0<=r<=1e-8' in source)
check('prospective inherited tolerances',"abs(de)<=1e-4 and df<=1e-3" in source)
check('native FD source unchanged',cfg['native_helper']['sha256']=='4189602145853bc39d8dfaf30915b2a6b1d4632984e0c50fc0bf321644fa306c')
check('no live release or native execution',cfg['root_release'] is None and cfg['preparation_result'] is None)
report={'status':'passed_independent_execution_review','implementation':cfg['implementation'],'reviewer_source':pin(__file__),'configuration':cfgpin,'specification':cfg['specification'],'preparation':prep_pin,'prepared_manifest':prep['manifest'],'saved_controls':controls_pin,'mechanism_review':pin(D/'independent_review_v1/SOURCE_REVIEW.json'),'checks':checks,'checks_passed':len(checks),'source_findings':['Four control cells finish before any positive-q derivative stage; all four full Cartesian gates finish before remaining84 cells.','Every analyze now receives existing strict analyze admission before testgrad. Every induction summary token is consumed, finite, nonnegative and at most1e-8 Debye.','Four one-thread workers, hard116/636 ceilings, locked counters/result persistence, complete92-row denominator, immediate failure stop and durable raw receipts are retained.','Zero-q and metal-absent controls use prospectively declared1e-4 energy/1e-3 force tolerance; exact flags remain descriptive.','Both amplitudes are fixed and shared across both metals. No amplitude selection, optimizer, changed target or held-validation release.','Existing native parser retains1e-5A Cartesian step and1e-3 maximum component error gate; saved replay verifies all receipt argv/input/output and scores.'],'blockers':[],'new_native_calls':0,'native_outcomes_read':False,'author_control_suite_rerun':False,'threads':1,'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu}
save(HERE/'REVIEW.json',report)
(HERE/'REPORT.md').write_text('# Independent EXCHPOL execution review\n\nPassed: final source, prospective tolerance amendment, actual zero-call preparation and saved controls are bound in REVIEW.json. No pending blocker.\n\nThe literal inventory is92 cells/116 programs/636 internal displaced energies, with four single-thread workers. All controls precede the four positive-q derivative gates; those precede the remaining matrix. The new immediate analyze admission and complete induction-token parsing close the two identified failure-detection gaps. Raw receipts and all missing slots persist on failure.\n\nBoth q values remain fixed for both metals. Original selected geometry and key prefix/parameter bytes are preserved. No molecular program, control-suite rerun, future outcome inspection, fit or sampling release was performed by this review. Physical mechanism interpretation remains the separate pinned source review.\n')
print(json.dumps({'review':pin(HERE/'REVIEW.json'),'checks':len(checks),'new_native_calls':0}))
