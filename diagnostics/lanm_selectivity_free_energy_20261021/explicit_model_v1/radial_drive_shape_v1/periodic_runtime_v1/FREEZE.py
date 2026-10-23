"""Freeze candidate-independent source metadata only; imports no scientific code."""
import ast
import datetime
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'diagnostics').is_dir() and (p/'workspaces').is_dir())
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
RD=Path('explicit_model_v1/radial_drive_shape_v1')
J=Path('explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1')
PN=J/'periodic_native_qualification_v1';S2=PN/'stage2_certificate_v1';EX=S2/'execution_v1'
OUT=W/RD/'periodic_runtime_v1/preparation_v1'

def pin(p):
    p=Path(p).resolve()
    with p.open('rb') as f: h=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(p),'sha256':h}
def read(p):return json.loads(Path(p).read_text())
def verify(p):
    assert pin(p['path'])==p,p
    return read(p['path'])
def write(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def literal(path,name):
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in node.targets):return ast.literal_eval(node.value)
    raise ValueError(name)

def main():
    OUT.mkdir(parents=True,exist_ok=False)
    declaration=read(HERE/'ROOT_DECLARATION.json');parents=declaration['parents']
    for p in parents.values():assert pin(p['path'])==p,p
    native_pin=pin(W/RD/'periodic_runtime_v1/native_source_v1/MANIFEST.json');native=verify(native_pin)
    controls_pin=pin(W/RD/'periodic_runtime_v1/controls_v1/RESULT.json');controls=verify(controls_pin)
    assert controls['status']=='passed_periodic_v2_runtime_source_controls'
    assert controls['check_count']==64 and controls['error'] is None
    for p in controls['source_pins'].values():assert pin(p['path'])==p,p
    execution={n:pin(HERE/n) for n in ['COMMON.py','POLICY.py','METADATA.py','NATIVE.py','CERTIFY.py','PARITY_REPLAY.py','GRADIENT_PARSE.py','ROOT_DECLARATION.json']}
    identity={'IDENTITY.py':parents['identity_identity_validator'],'STATE.py':parents['identity_parser']}
    settings=literal(identity['IDENTITY.py']['path'],'SETTINGS')
    old=verify(parents['execution_manifest']);recipe=verify(old['stage1_link_recipe'])
    radial=verify(parents['radial_native_manifest'])
    certificate_sources={n:execution[n] for n in ['COMMON.py','POLICY.py','METADATA.py','CERTIFY.py']}
    certificate_sources.update(identity)
    implementations={'native':execution['NATIVE.py'],'certificate':execution['CERTIFY.py']}
    oldbuild_pin=pin(W/S2/'cutoff12_baseline_v1/build_run_v1/execution/RESULT.json');oldbuild=verify(oldbuild_pin)
    oldroot=pin(W/S2/'cutoff12_baseline_v1/build_run_v1/ROOT_ACTUAL_REVIEW.json')
    assert verify(oldroot)['result']==oldbuild_pin
    oldbaseline={'executable':oldbuild['executables']['gradient_old'],'library':oldbuild['library'],
        'polar_module':next(p for p in recipe['module_files'] if p['path'].endswith('/polar.mod')),
        'gate':parents['cutoff12_gate_actual_parent'],
        'reporter':native['files']['GRADIENT.f90'],'parser':parents['gradient_parser'],
        'root_admission':oldroot,'build_result':oldbuild_pin,'ticket_schema':'STAGE2_BOUND_INPUT_V1','state_schema':'PERIODIC_PARITY_STATE_V1'}
    build={
      'status':'candidate_independent_build_handoff_not_released',
      'radial_source_manifest':parents['radial_native_manifest'],'runtime_native_manifest':native_pin,
      'source_inventory':{'radial_archive_files':radial['source_outputs'],'isolated_objects':native['files'],
                          'normal_source_guards_remain_fatal':True},
      'inherited_recipe_authorities':{
          'archive_builder':pin(D/J/'periodic_stability_definition_v1/pme_export_v1/BUILD_LOGIN_ONE.py'),
          'successful_parent_build':recipe['parent_build'],'wrapper_builder':pin(D/EX/'BUILD.py'),
          'wrapper_link_recipe':old['stage1_link_recipe']},
      'compiler_and_link_dependencies':recipe['compiler_and_link_dependencies'],
      'archive_configuration_argv':['cmake','-S',str(Path(radial['source_root'])/'cmake'),'-B','<fresh-build-dir>',
          '-DTINKER_DIR='+radial['source_root'],'-DCMAKE_BUILD_TYPE=RELWITHDEBINFO',
          '-DCMAKE_Fortran_COMPILER=/usr/bin/f95','-DCMAKE_Fortran_FLAGS_RELWITHDEBINFO=-O2 -g -fbacktrace',
          '-DFFTW_LIB=/usr/lib/x86_64-linux-gnu/libfftw3.a',
          '-DFFTW_THREADED_LIB=/usr/lib/x86_64-linux-gnu/libfftw3_threads.a','-DTINKER_OPENMM=OFF'],
      'archive_build_argv':['cmake','--build','<fresh-build-dir>','--target','tinker','--parallel','1'],
      'archive_target_note':'tinker is the add_library target in the pinned CMakeLists; building only this target avoids unused analyze/testgrad executables. This is a prospective target change from the parent recipe and needs build-adapter review.',
      'wrapper_flags':recipe['flags'],'new_architecture_flags':[],
      'compile_argv_pattern':['/usr/bin/f95',*recipe['flags'],'-I<new-radial-build>/mod','-I<fresh-wrapper-mod>','-J<fresh-wrapper-mod>','-c','<pinned-source>','-o','<fresh-object>'],
      'compile_order':['GATE_V2.f90','STATE_WRITER_V2.f90','STATE_ENTRY_POINTS.f90','induce_qualification.f','INITIALIZE.f90','GRADIENT.f90','PARITY_SETUP.f','PARITY_V2.f90'],
      'link_objects_before_archive':{
          'initialize':['INITIALIZE.f90.o','GATE_V2.f90.o','STATE_WRITER_V2.f90.o','STATE_ENTRY_POINTS.f90.o','induce_qualification.f.o'],
          'gradient_clone':['GRADIENT.f90.o','GATE_V2.f90.o','STATE_WRITER_V2.f90.o','STATE_ENTRY_POINTS.f90.o','induce_qualification.f.o'],
          'parity':['PARITY_V2.f90.o','PARITY_SETUP.f.o','STATE_WRITER_V2.f90.o','STATE_ENTRY_POINTS.f90.o']},
      'link_tail':['<new-radial-build>/libtinker.a',*recipe['argv'][recipe['argv'].index(recipe['library']['path'])+1:-2]],
      'old_baseline_reuse':oldbaseline,
      'archive_rules':['Compile all V2 objects against the new radial polar.mod containing lrpshape. Never mix old and radial module directories.',
          'The whole induce_qualification object precedes the radial archive only for initialize and gradient_clone. The parity link uses normal radial archive and isolated no-solve setup; it must not link the qualification induce object.',
          'Reuse the admitted old cutoff12 executable exactly, recording its old archive/module/gate/reporter/parser chain. No V2 writer or ticket enters that executable.',
          'Pin all new archive/module/object/executable outputs and loader evidence. Do not treat any current source pin as a successful build.'],
      'expected_build_result':{'status':'periodic_v2_build_complete','native_sources':native_pin,
          'executables':{'initialize':None,'parity':None,'gradient_clone':None,'gradient_old':None},
          'radial_library':None,'radial_polar_module':None,'old_baseline':oldbaseline},
      'expected_runtime_outputs':{
          'initialize':['certificate_actual_state.dat (V2)','pme_norm_export_001.dat','stage2_initialization_receipt.dat'],
          'parity':['parity_state.dat (V2)','pme_norm_export_001.dat','original ten-probe action files'],
          'gradient_clone':['stage2_gradient.dat (V2 observation; active full state before solve)','certificate_actual_state.dat for ACTIVE'],
          'gradient_old':['stage2_gradient.dat (original V1 bypass observations)']},
      'node_requirements':'One build worker; BLAS/OMP1; unset LD_PRELOAD/LD_LIBRARY_PATH. Verify compiler/dependencies and actual loader on the allocated node before any later application. Shared paths and absence of march/mtune flags are not node admission.',
      'candidate':None,'build_release':None,'build_result':None,'native_calls':0}
    # Verify the named library target without invoking CMake.
    cmake=Path(radial['source_root'])/'cmake/CMakeLists.txt'
    assert 'add_library (tinker STATIC' in cmake.read_text() or 'add_library(tinker STATIC' in cmake.read_text()
    write(HERE/'BUILD_HANDOFF.json',build)
    contract={'status':'source_only_periodic_v2_runtime_contract','declaration':pin(HERE/'ROOT_DECLARATION.json'),
      'execution_sources':execution,'native_sources':native_pin,'identity_sources':identity,
      'implementations':implementations,'certificate_sources':certificate_sources,'settings':settings,
      'phase_order':declaration['phase_order'],'central_displaced_identity':declaration['central_displaced_identity'],
      'initializer_scope':'One direct no-solve initialization and PME/site/summary replay; V2 header/settings/transform readback only. Full V2 state body, h/H and identity admission occur in CERTIFY before assembly.',
      'preserved':['Raw unequal direct1/direct2; explicit e900/e901 even at zero field or r0; cell and cutoff12 identity.',
          'Original real/reciprocal action routines, factor algorithm, ten probes, numerical gates, caps and no-retry invocation.',
          'Normal radial native backend/stability fatal guards; isolated qualification callsites only.',
          'V1 old absent bypass has independent archive/module/grammar; no legacy central/state/certificate reuse.'],
      'future_bindings_required':['Explicit candidate/composition and literal all-state inputs with explicit polar-eps1e-8 plus both response records.',
          'Reviewed new radial archive and four-role build bundle (old absent executable may be exact reuse).',
          'Fresh central initialize/certificate/root admission, then central operator parity/root admission, then finite central/bypass scout.',
          'Only after scout: each displaced ACTIVE row requires same campaign parity plus its own fresh state/certificate/root admission.',
          'Separate downstream collector configuration must use the declared anchored strain projection and printing allowance; this package does not release or collect a finite E/F panel.'],
      'prospective_counts':declaration['prospective_counts'],'source_controls':controls_pin,
      'limitations':['No actual full 64-cube V2 body parse, build, link or runtime validation performed.',
          'A nominal finite same-state certificate and sampled parity do not establish trajectory-wide stability or physical accuracy.',
          'Candidate, runtime inputs, build and all numerical releases remain null.'],
      'counts':{'native_calls':0,'builds':0,'FFT_replays':0,'assembly_calls':0,'factorizations':0,'solver_calls':0}}
    write(HERE/'SOURCE_CONTRACT.json',contract)
    manifest={'status':'unbound_periodic_v2_runtime_inputs','campaign':None,'state_schema':'PERIODIC_PARITY_STATE_V2',
       'model':{'candidate':None,'composition':None,'coefficients':None,'radial':None,'settings':settings},
       'execution_sources':execution,'native_sources':native_pin,'identity_sources':identity,
       'native_source_admission':parents['radial_port_source_review'],'identity_source_admission':parents['identity_root_review'],
       'helpers':old['helpers'],'nominal_contract':parents['nominal_contract'],
       'environment':None,'native_runtime':None,'loader_command':None,'physics_metadata':None,
       'probe_source':None,'probe_data':None,'rows':{}}
    write(HERE/'MANIFEST_TEMPLATE.json',manifest)
    write(HERE/'ROW_TEMPLATE.json',{'state_id':None,'role':None,'cell':None,
       'source_inputs':{'xyz':None,'key':None,'parameters':None},'inputs':{'xyz':None,'key':None,'parameters':None}})
    for kind in ('native','certificate'):
        cfg={'status':'unbound_periodic_v2_'+kind+'_configuration','implementation':implementations[kind],
             'manifest':None,'build':None,'build_admission':None,'state_id':None,
             'central_parity_admission':None,'central_scout_review':None}
        rel={'status':'unreleased_periodic_v2_'+kind,'implementation':implementations[kind],
             'configuration':None,'output':None,'source_review':None,'state_id':None}
        if kind=='native':
            cfg.update(operation=None,certificate=None,certificate_admission=None,probe_source=None,probe_data=None)
            rel.update(operation=None,maximum_native_programs=1)
        else:
            cfg.update(initializer=None,caps={'assembly':{'wall':180,'AS':8*1024**3,'CPU':180,'BLAS':1},'factor':{'wall':900,'AS':32*1024**3,'CPU':10800,'BLAS':12}},nominal_contract=parents['nominal_contract'],
                       certificate_sources=certificate_sources,authorities={'identity_root_review':parents['identity_root_review'],'radial_source_review':parents['radial_port_source_review']})
            rel.update(maximum_assemblies=1,maximum_numeric_factors=1,maximum_solve_probes=10)
        write(HERE/(kind.upper()+'_CONFIG_TEMPLATE.json'),cfg);write(HERE/(kind.upper()+'_RELEASE_TEMPLATE.json'),rel)
    template={'status':'passed_periodic_v2_runtime_source_review','execution_sources':execution,'native_sources':native_pin,
       'implementations':implementations,'certificate_sources':certificate_sources,'source_contract':pin(HERE/'SOURCE_CONTRACT.json'),
       'controls':controls_pin,'ready':None,'package':None,'independent_reviewer_implementation':None}
    write(HERE/'SOURCE_REVIEW_TEMPLATE.json',template)
    package={'status':'frozen_periodic_v2_runtime_source_package','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
       'producer':pin(__file__),'declaration':pin(HERE/'ROOT_DECLARATION.json'),'plan':declaration['source_plan'],
       'plan_review':declaration['independent_design_review'],'execution_sources':execution,'native_sources':native_pin,
       'implementations':implementations,'certificate_sources':certificate_sources,'helpers':old['helpers'],
       'source_generator':pin(HERE/'GENERATE_NATIVE.py'),'test_source':pin(HERE/'TEST_RUNTIME.py'),'controls':controls_pin,
       'metadata':{n:pin(HERE/n) for n in ['BUILD_HANDOFF.json','SOURCE_CONTRACT.json','MANIFEST_TEMPLATE.json','ROW_TEMPLATE.json',
           'NATIVE_CONFIG_TEMPLATE.json','NATIVE_RELEASE_TEMPLATE.json','CERTIFICATE_CONFIG_TEMPLATE.json','CERTIFICATE_RELEASE_TEMPLATE.json','SOURCE_REVIEW_TEMPLATE.json']},
       'candidate':None,'build':None,'runtime_release':None,'counts':contract['counts'],
       'source_generation_cost':native['cost'],'manufactured_control_cost':controls['cost']}
    write(OUT/'MANIFEST.json',package)
    ready={'status':'ready_for_independent_periodic_v2_runtime_source_review','package':pin(OUT/'MANIFEST.json'),
       'execution_sources':execution,'native_sources':native_pin,'implementations':implementations,
       'certificate_sources':certificate_sources,'source_contract':pin(HERE/'SOURCE_CONTRACT.json'),
       'build_handoff':pin(HERE/'BUILD_HANDOFF.json'),'controls':controls_pin,
       'required_review_template':pin(HERE/'SOURCE_REVIEW_TEMPLATE.json'),
       'required_review_status':'passed_periodic_v2_runtime_source_review','candidate':None,'composition':None,
       'build_release':None,'numerical_release':None,'counts':contract['counts'],
       'note':'Plan READY.json remains unchanged. This implementation handoff creates no run or build authorization.'}
    write(HERE/'IMPLEMENTATION_READY.json',ready)
    print(json.dumps({'ready':pin(HERE/'IMPLEMENTATION_READY.json'),'package':pin(OUT/'MANIFEST.json')}))

if __name__=='__main__':main()
