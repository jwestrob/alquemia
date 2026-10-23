"""Twelve fixed inactive strain endpoints; no build or central invocation."""
import hashlib,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;PARENT=HERE.parent/'execution_v1'
if hashlib.sha256((PARENT/'COMMON.py').read_bytes()).hexdigest()!='b13a4d20fbd9cc0f17db35b72244c2f9cd4d0d263b65466a7bce07e47decdb6e':raise ValueError('Changed runtime helpers')
sys.path.insert(0,str(PARENT))
from COMMON import *
DIRECTIONS=('strain_0','strain_1','strain_2');STEP_LABELS=('1.25e-05','6.25e-06')
IDS=tuple(d+'__'+h+'__'+s for d in DIRECTIONS for h in STEP_LABELS for s in ('minus','plus'))
DECLARATION_SHA='54936c17d7ecca77ebc4aaa56a2461668008ebd75a3e2847341a758709ae135b'

def validate_manifest(m):
 require(m['status']=='prepared_twelve_anchored_strain_refinement_inputs' and tuple(m['rows'])==IDS and m['maximum_native_programs']==12,'Exact12 finite endpoints')
 require(m['root_declaration']['sha256']==DECLARATION_SHA,'Exact root scope');d=read(m['root_declaration'])
 require(d['directions_to_test']==list(DIRECTIONS) and d['half_steps_A']==[1.25e-5,6.25e-6] and d['maximum_new_gradient_programs']==12 and d['new_central_programs']==0,'Declared direction/step/counts')
 authority=read(m['source_inputs']);require(authority['declaration']==m['root_declaration'],'Source authority');base=read(m['baseline_manifest']);parent=load(authority['baseline_binding'],'strain_original_binding');parent.validate_manifest(base)
 build=read(m['build']);require(m['build']==d['parent_build'] and build['manifest']==m['baseline_manifest'] and build['source_manifest']==base['native_sources'] and build['status']=='stage2_build_complete','Original cutoff12 build lineage')
 br=read(m['build_root_review']);require(br['status']=='admitted_cutoff12_baseline_build' and br['result']==m['build'] and br['manifest']==m['baseline_manifest'],'Original root build admission')
 cm=read(d['central_manifest']);require(m['central_gradient']==d['central_gradient']==cm['rows']['central']['baseline_gradient'] and m['central_result']==cm['rows']['central']['baseline_result'],'Exact admitted central reuse')
 central=read(m['central_result']);require(central['status']=='stage2_gradient_complete' and central['inputs']==base['rows']['central']['inputs'] and central['gradient_output']==m['central_gradient'],'Central input/result identity');verify(m['central_gradient'])
 evidence=read(d['evidence']);require(evidence['status']=='reviewed_exact_term_attribution_electrostatic_strain_residual' and read(evidence['result'])['complete'],'Completed attribution prerequisite')
 for k in ('helpers','collector','environment','native_runtime','loader_command','candidate','composition','native_sources','directions','displaced_manifest'):
  require(m[k]==base[k],'Unchanged authority '+k)
 require(m['directions']==d['directions'],'Fixed directions');verify(m['directions'])
 for x in m['execution_sources'].values():verify(x)
 require(m['execution_sources']['RUN.py']==pin(HERE/'RUN.py') and m['execution_sources']['BINDING.py']==pin(__file__),'Actual versioned runtime')
 require(verify(m['execution_sources']['RUN.py']).read_bytes()==verify(authority['baseline_runner']).read_bytes(),'Unchanged native runner')
 require(m['geometry_preparer']==authority['geometry_preparer'],'Unchanged geometry/serialization methods');verify(m['geometry_preparer'])
 require(len({r['output'] for r in m['rows'].values()})==12,'12 unique outputs')
 for sid,row in m['rows'].items():
  require(row['role']=='BYPASS_OLD_ABSENT' and row['inputs']==row['source_inputs'],'Exact inactive runtime inputs; observer suffix already inherited')
  for k in ('key','parameters'):require(row['inputs'][k]==central['inputs'][k],'Same central '+k)
  require(row['direction']==sid.split('__')[0] and row['half_step_A']==float(sid.split('__')[1]) and row['sign']==(-1 if sid.endswith('minus') else 1),'State identity')
  for x in row['inputs'].values():verify(x)
 return base,build

def source_context(config,release,out,implementation,kind,cpus,blas=1):
 cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text())
 require(kind=='native' and cpus==blas==1 and cfg['operation']=='gradient','One gradient only; no build/initializer')
 require(cfg['status']=='bound_anchored_strain_refinement_configuration' and cfg['implementation']==rel['implementation']==pin(implementation),'Exact runtime source')
 require(rel['status']=='released_anchored_strain_refinement' and rel['configuration']==pin(config) and Path(rel['output']).resolve()==Path(out).resolve(),'Exact root release/output')
 m=read(cfg['manifest']);validate_manifest(m);sid=cfg['state_id'];require(sid in IDS and Path(out).resolve()==Path(m['rows'][sid]['output']).resolve(),'Fixed state/output')
 require(cfg['build']==m['build'] and cfg['build_root_review']==m['build_root_review'] and cfg['certificate'] is None and cfg['certificate_admission'] is None,'Unchanged baseline build/no active certificate')
 review=read(rel['source_review']);require(review['status']=='passed_anchored_strain_refinement_source_review' and review['manifest']==cfg['manifest'] and review['implementation']==cfg['implementation'] and review['binding']==pin(__file__),'Frozen source review')
 env=read(m['environment']);require(Path(sys.prefix).resolve()==Path(env['prefix']).resolve(),'Pinned isolated interpreter')
 for x in env['runtime_files']:verify(x)
 require(all(os.environ.get(k)=='1' for k in THREADS),'One-thread environment')
 return cfg,m,rel,allocation(m,1)

def native_build(cfg,m):
 b=read(cfg['build']);require(cfg['build']==m['build'] and b['status']=='stage2_build_complete' and b['manifest']==m['baseline_manifest'] and b['source_manifest']==m['native_sources'],'Exact reused binary lineage');verify(b['executables']['gradient_old']);return b
