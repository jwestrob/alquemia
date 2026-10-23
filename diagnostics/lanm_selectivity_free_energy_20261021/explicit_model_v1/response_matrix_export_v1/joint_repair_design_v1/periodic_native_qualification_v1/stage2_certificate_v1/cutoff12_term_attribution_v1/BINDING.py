"""Five exact baseline states; observer-only source delta and one-call policy."""
import hashlib,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'execution_v1'
COMMON_SHA='b13a4d20fbd9cc0f17db35b72244c2f9cd4d0d263b65466a7bce07e47decdb6e'
if hashlib.sha256((PARENT/'COMMON.py').read_bytes()).hexdigest()!=COMMON_SHA:raise ValueError('Changed runtime helpers')
sys.path.insert(0,str(PARENT))
from COMMON import *
IDS=('central',)+tuple('strain_1__'+h+'__'+s for h in ('1e-04','5e-05') for s in ('minus','plus'))
DECLARATION_SHA='7b8af01e4c8491292802a5c25fbf7efae5f994a490cfec6eaced96316c76b3f1'
def validate_manifest(m):
 require(m['status']=='prepared_five_term_observer_inputs' and tuple(m['rows'])==IDS and m['maximum_native_programs']==5,'Exact five states')
 require(m['root_declaration']['sha256']==DECLARATION_SHA,'Prospective declaration')
 decl=read(m['root_declaration']);require(tuple(decl['states'])==IDS and decl['maximum_new_native_gradient_programs']==5,'Declared finite panel')
 authority=read(m['source_inputs']);require(authority['declaration']==m['root_declaration'],'Preparation authority')
 base=read(decl['baseline_manifest']);load(authority['original_binding'],'original_cutoff_binding').validate_manifest(base)
 collection=read(decl['cutoff12_collection']);review=read(decl['evidence']);require(review['result']==decl['cutoff12_collection'],'Admitted original failed total gates')
 require(collection['complete'] and collection['counts']['admitted_new_gradients']==9,'Complete old nine baseline states')
 build=read(decl['baseline_build']);require(build['status']=='stage2_build_complete' and build['error'] is None and collection['build']==decl['baseline_build'],'Original build authority')
 require(m['parent_build']==decl['baseline_build'],'Exact reused build')
 for name in ('helpers','collector','environment','native_runtime','loader_command','candidate','composition','stage1_link_recipe','directions','displaced_manifest'):
  require(m[name]==base[name],'Unchanged scientific authority '+name)
 for item in m['execution_sources'].values():verify(item)
 require(m['execution_sources']['RUN.py']==pin(HERE/'RUN.py') and m['execution_sources']['BINDING.py']==pin(__file__) and m['observer']==pin(HERE/'OBSERVER.py'),'Actual sources')
 source=read(m['native_sources']);original=read(base['native_sources']);observer=load(m['observer'],'term_source_observer')
 require(source['parent_manifest']==base['native_sources'] and source['original_gradient']==authority['original_gradient'],'Native lineage')
 oldsource=read(authority['original_source_manifest']);require(authority['original_gradient'] in oldsource['source_outputs'],'Original gradient manifest membership')
 oldbuild=read(authority['original_build']);require(oldbuild['manifest']==authority['original_source_manifest'] and oldbuild['library']==read(m['stage1_link_recipe'])['library'],'Unchanged archive build/source')
 old=verify(authority['original_gradient']).read_text();new=verify(source['files']['GRADIENT_OBSERVER.f']).read_text()
 require(new==observer.render(old) and observer.inverse(new)==old,'Only exact insertion observer')
 for name in ('GATE.f90','GRADIENT.f90','STATE_WRITER.f90'):
  require(source['files'][name]==original['files'][name],'Unchanged main/gate/writer source '+name);verify(source['files'][name]);verify(m['reused_objects'][name])
  if name=='GATE.f90':
   cmd=build['commands']['compile_GATE_f90'];require(cmd['returncode']==0 and cmd['argv'][cmd['argv'].index('-o')+1]==m['reused_objects'][name]['path'],'Actual cutoff12 gate object')
  else:require(m['reused_objects'][name]==base['reused_objects'][name],'Original harness object')
 require(len({r['output'] for r in m['rows'].values()})==5,'Five immutable outputs')
 for sid,row in m['rows'].items():
  require(row['role']=='BYPASS_OLD_ABSENT' and row['inputs']==base['rows'][sid]['inputs'] and row['source_inputs']==base['rows'][sid]['source_inputs'],'Exact same inactive inputs')
  for item in row['inputs'].values():verify(item)
  oldrow=read(row['baseline_result']);require(row['baseline_result']==collection['rows'][sid]['result'] and row['baseline_gradient']==collection['rows'][sid]['gradient_output']==oldrow['gradient_output'],'Exact saved output authority')
  require(oldrow['inputs']==row['inputs'] and oldrow['state_id']==sid and oldrow['status']=='stage2_gradient_complete','Baseline row identity');verify(row['baseline_gradient'])
 return base,build

def source_context(config,release,out,implementation,kind,cpus,blas=1):
 cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text())
 require(kind in ('build','native') and cpus==blas==1,'One CPU finite phase')
 require(cfg['status']=='bound_term_observer_'+kind+'_configuration' and cfg['implementation']==rel['implementation']==pin(implementation),'Bound exact phase source')
 require(rel['status']=='released_term_observer_'+kind and rel['configuration']==pin(config) and Path(rel['output']).resolve()==Path(out).resolve(),'Exact root release/output')
 m=read(cfg['manifest']);validate_manifest(m);review=read(rel['source_review'])
 require(review['status']=='passed_cutoff12_term_observer_source_review' and review['manifest']==cfg['manifest'] and review['binding']==pin(__file__) and review['implementations'][kind]==cfg['implementation'],'Exact source review')
 env=read(m['environment']);require(Path(sys.prefix).resolve()==Path(env['prefix']).resolve(),'Pinned interpreter')
 for item in env['runtime_files']:verify(item)
 require(all(os.environ.get(k)=='1' for k in THREADS),'One thread environment')
 if kind=='build':require(Path(out).resolve()==Path(m['build_output']).resolve(),'One immutable build')
 else:
  require(cfg['operation']=='gradient' and cfg['state_id'] in IDS and cfg['certificate'] is None and cfg['certificate_admission'] is None,'Inactive gradient only')
  require(Path(out).resolve()==Path(m['rows'][cfg['state_id']]['output']).resolve(),'Fixed state output')
  native_build(cfg,m);admission=read(cfg['build_root_review'])
  require(admission['status']=='admitted_cutoff12_term_observer_build' and admission['result']==cfg['build'] and admission['manifest']==cfg['manifest'],'Actual root build admission')
 return cfg,m,rel,allocation(m,1)

def native_build(cfg,m):
 b=read(cfg['build'])
 require(b['status']=='stage2_build_complete' and b['manifest']==cfg['manifest'] and b['source_manifest']==m['native_sources'] and b['compile_attempts']==b['link_attempts']==1 and b['native_programs']==0 and b['error'] is None,'One successful observer compile/link')
 require(set(b['executables'])=={'gradient_old'} and b['library']==read(m['stage1_link_recipe'])['library'],'Unchanged native archive')
 verify(b['executables']['gradient_old']);return b
