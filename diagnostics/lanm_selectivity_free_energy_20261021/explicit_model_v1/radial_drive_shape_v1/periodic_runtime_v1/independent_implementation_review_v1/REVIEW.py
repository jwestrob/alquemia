"""Read-only source/metadata receipt; imports no author or numerical code."""
import datetime, hashlib, json, time
from pathlib import Path

HERE=Path(__file__).resolve().parent
P=HERE.parent
CHECKS=[]
PINS={}
def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as f: sha=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(path),'sha256':sha}
def verify(p):
    if p['path'] not in PINS:
        assert pin(p['path'])==p,p
        PINS[p['path']]=p
    else: assert PINS[p['path']]==p
    return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def check(label,value):
    assert value,label
    CHECKS.append({'label':label,'passed':True})
def main():
    wall=time.monotonic();cpu=time.process_time()
    ready_pin=pin(P/'IMPLEMENTATION_READY.json')
    check('Exact requested frozen READY',ready_pin['sha256']=='2f24ade0b91add734a5b84f7c260de65ff97981830cd74f53b872e470220de6c')
    r=read(ready_pin);pkg=read(r['package']);native=read(r['native_sources']);controls=read(r['controls'])
    contract=read(r['source_contract']);build=read(r['build_handoff']);template=read(r['required_review_template'])
    for group in ('execution_sources','certificate_sources','implementations'):
        check('Exact '+group+' package/template joins',r[group]==pkg[group]==template[group])
        for p in r[group].values():verify(p)
    for k in ('native_sources','controls'):check('Exact '+k+' joins',r[k]==pkg[k]==template[k])
    for p in pkg['metadata'].values():verify(p)
    for p in native['files'].values():verify(p)
    inverse=read(native['inverse'])
    for k in ('parent_induce','parent_reporter'):verify(inverse[k])
    for p in controls['source_pins'].values():verify(p)
    for k in ('source_generator','test_source','producer','declaration','plan','plan_review'):verify(pkg[k])
    check('64 author controls passed; read only',controls['status']=='passed_periodic_v2_runtime_source_controls' and controls['check_count']==64 and len(controls['checks'])==64 and all(x['passed'] for x in controls['checks']) and controls['error'] is None)
    check('Controls bind exact generated native source',controls['native_manifest']==r['native_sources'])
    check('Control and source-generation work contains no scientific execution',all(v==0 for v in controls['counts'].values()) and all(v==0 for v in native['counts'].values()))
    check('Normal fatal guards retained in generated source record',native['normal_fatal_guards_unchanged'] is True)
    check('Candidate and execution remain unbound',all(r[k] is None for k in ('candidate','composition','build_release','numerical_release')) and all(pkg[k] is None for k in ('candidate','build','runtime_release')))
    decl=read(pkg['declaration']);radial=read(build['radial_source_manifest'])
    check('Build source inventory exactly inherited',build['source_inventory']['radial_archive_files']==radial['source_outputs'] and len(radial['source_outputs'])==615)
    check('Build wrapper objects exactly generated',build['runtime_native_manifest']==r['native_sources'] and build['source_inventory']['isolated_objects']==native['files'])
    check('Radial and identity source authorities unchanged',build['radial_source_manifest']==decl['parents']['radial_native_manifest'] and contract['identity_sources']=={'IDENTITY.py':decl['parents']['identity_identity_validator'],'STATE.py':decl['parents']['identity_parser']})
    for p in contract['identity_sources'].values():verify(p)
    idnative=read(decl['parents']['identity_source_manifest'])
    check('Shared V2 writer/gate bytes exactly admitted',all(native['files'][k]['sha256']==v['sha256'] for k,v in idnative['outputs'].items()))
    old=build['old_baseline_reuse'];oldbuild=read(old['build_result']);oldroot=read(old['root_admission'])
    check('Old cutoff12 executable and archive metadata join actual build',old['executable']==oldbuild['executables']['gradient_old'] and old['library']==oldbuild['library'] and oldroot['result']==old['build_result'])
    check('Old bypass preserves V1 grammar and independent module lineage',old['ticket_schema']=='STAGE2_BOUND_INPUT_V1' and old['state_schema']=='PERIODIC_PARITY_STATE_V1' and '/pme_export_v1/build_v2/' in old['polar_module']['path'])
    check('Matched radial modules; parity excludes qualification induce',build['compile_argv_pattern'][6]=='-I<new-radial-build>/mod' and 'induce_qualification.f.o' not in build['link_objects_before_archive']['parity'] and all('induce_qualification.f.o' in build['link_objects_before_archive'][k] for k in ('initialize','gradient_clone')))
    check('One-worker prospective library target and unchanged optimization flags',build['archive_build_argv'][-4:]==['--target','tinker','--parallel','1'] and '-DTINKER_OPENMM=OFF' in build['archive_configuration_argv'] and build['new_architecture_flags']==[])
    check('Initializer body scope explicit before certificate assembly','Full V2 state body, h/H and identity admission occur in CERTIFY before assembly.' in contract['initializer_scope'])
    check('No old central reuse and exact prospective program accounting',contract['prospective_counts']==decl['prospective_counts'] and '70nativeprograms,33freshcertificates' in contract['prospective_counts']['if_original32_displacement_and36_gradient_structure_is_later_declared'])
    manifest=json.loads((P/'MANIFEST_TEMPLATE.json').read_text())
    check('Manifest is intentionally unbound',manifest['status']=='unbound_periodic_v2_runtime_inputs' and manifest['campaign'] is None and manifest['rows']=={} and manifest['model']['candidate'] is None)
    result=dict(template)
    result.update(ready=ready_pin,package=r['package'],independent_reviewer_implementation=pin(__file__),build_handoff=r['build_handoff'],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=CHECKS,verified_pins=list(PINS.values()),review_scope='Independent source and saved author-control review only; no imported runtime, control rerun, build, native, FFT, assembly, factor or solve.',findings=[],controls_rerun=False,source_only=True,build_admitted=False,numerical_admission=False,candidate=None,limitations=['Initializer completion admits PME/header readback only; full V2 body/h/H is checked by CERTIFY before assembly.','Actual new radial archive, matching modules and three new links remain unbuilt; old absent executable remains its distinct admitted V1 lineage.','Nominal finite stability and sampled parity are ordinary-double checks; e_app remains unknown and no trajectory-wide or physical-accuracy claim follows.','Future finite force collection must bind anchored strain projections and declared printing allowances; no E/F panel is released here.'],counts={'author_control_checks_read':64,'new_control_runs':0,'native_calls':0,'builds':0,'FFT_replays':0,'assemblies':0,'factorizations':0,'solver_calls':0},cost={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1})
    report='Independent periodic V2 runtime implementation source review: PASS\n\nThe phase policy removes central bootstrap circularity and requires a fresh cutoff12 initialized state, assembly and nominal factor before central parity. Displacements reuse campaign-level parity/scout only and retain their own complete exact-state chain. Raw unequal Ed/Ep reach the shared V2 writer; corrected RHS is observed separately. The generated Stage1 action loops and numerical factor/parity routines retain their inherited equations, signs, tolerances and finite limits.\n\nBoth radial slots, all original coefficients, settings and actual state remain explicit. Old absent bypass uses its separate V1 ticket/parser and admitted cutoff12 executable; all new V2 objects require the matching new radial archive/polar module. The eight-object/three-link handoff leaves normal fatal guard bodies intact. No cutoff9 state/CSC/PME/certificate shortcut remains.\n\nRead 64 passed author controls without rerunning them. Candidate, build and runtime releases remain unbound. Header-only initializer scope is explicit; CERTIFY performs full body/h/H checks before assembly. This receipt admits source correspondence only, not a build, molecular execution, uniform application-error bound, trajectory stability or physical accuracy.\n'
    (HERE/'REPORT.txt').write_text(report)
    result['report']=pin(HERE/'REPORT.txt')
    with (HERE/'REVIEW.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'review':pin(HERE/'REVIEW.json'),'checks':len(CHECKS),'pins':len(PINS),'cost':result['cost']}))
if __name__=='__main__':main()
