"""Metadata/source preparation only; no compiler, native or numerical replay."""
from BINDING import *
import difflib

def main():
 tick=time.monotonic();cpu=time.process_time();authority=read(pin(HERE/'INPUTS.json'));decl=read(authority['declaration'])
 base=read(decl['baseline_manifest']);build=read(decl['baseline_build']);collection=read(decl['cutoff12_collection'])
 root=Path(decl['baseline_manifest']['path']).parents[2]/'cutoff12_term_attribution_v1'
 out=root/'preparation_v1';out.mkdir(parents=True,exist_ok=False);source=root/'source_v1';source.mkdir(exist_ok=False)
 observer=load(pin(HERE/'OBSERVER.py'),'term_prepare_observer');old=verify(authority['original_gradient']).read_text();new=observer.render(old)
 (source/'GRADIENT_OBSERVER.f').write_text(new)
 (source/'SOURCE.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='original/gradient.f',tofile='observer/gradient.f')))
 original=read(base['native_sources']);files=dict(original['files']);files['GRADIENT_OBSERVER.f']=pin(source/'GRADIENT_OBSERVER.f')
 native=dict(status='prepared_gradient_observer_source',parent_manifest=base['native_sources'],original_gradient=authority['original_gradient'],generator=pin(HERE/'OBSERVER.py'),files=files,diff=pin(source/'SOURCE.diff'),original_inverse_exact=True)
 save(source/'MANIFEST.json',native)
 m={k:base[k] for k in ('helpers','collector','environment','native_runtime','loader_command','candidate','composition','stage1_link_recipe','directions','displaced_manifest')}
 objects=dict(base['reused_objects']);objects['GATE.f90']=pin(build['commands']['compile_GATE_f90']['argv'][-1])
 m.update(status='prepared_five_term_observer_inputs',producer=pin(__file__),source_inputs=pin(HERE/'INPUTS.json'),root_declaration=authority['declaration'],parent_build=decl['baseline_build'],native_sources=pin(source/'MANIFEST.json'),observer=pin(HERE/'OBSERVER.py'),
  execution_sources={n:pin(HERE/n) for n in ('BINDING.py','BUILD.py','RUN.py','OBSERVER.py')},original_wrapper=authority['original_wrapper'],reused_objects=objects,maximum_native_programs=5,rows={},build_output=str(root/'build_run_v1/execution'),resources=base['resources'])
 for sid in IDS:
  oldrow=base['rows'][sid];r=collection['rows'][sid]
  m['rows'][sid]=dict(role='BYPASS_OLD_ABSENT',inputs=oldrow['inputs'],source_inputs=oldrow['source_inputs'],baseline_result=r['result'],baseline_gradient=r['gradient_output'],output=str(root/'run_v1/tasks'/sid/'execution'))
 validate_manifest(m);save(out/'MANIFEST.json',m)
 result=dict(status='term_observer_preparation_passed',implementation=pin(__file__),manifest=pin(out/'MANIFEST.json'),source_manifest=m['native_sources'],states=list(IDS),source_inverse_exact=True,new_geometries=0,native_calls=0,builds=0,scientific_replays=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 save(out/'RESULT.json',result);print(json.dumps(result))
if __name__=='__main__':main()
