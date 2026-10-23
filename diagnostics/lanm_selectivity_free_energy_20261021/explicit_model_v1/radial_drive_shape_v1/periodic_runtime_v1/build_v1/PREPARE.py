"""Metadata/hash-only build preparation; no compiler, subprocess or candidate read."""
import json
import shutil
from pathlib import Path
import BUILD as B

HERE=Path(__file__).resolve().parent
def main(out):
    out.mkdir(parents=True,exist_ok=False)
    declaration=B.read(B.pin(HERE/'ROOT_DECLARATION.json'));parents=declaration['parents']
    h=B.read(parents['BUILD_HANDOFF.json']);ready=B.read(parents['IMPLEMENTATION_READY.json'])
    oldmanifest=B.read(B.read(B.pin(HERE.parent/'ROOT_DECLARATION.json'))['parents']['execution_manifest'])
    root=B.pin(HERE.parent/'ROOT_SOURCE_REVIEW.json');admission=B.read(root)
    independent=admission['independent_review'];review=B.read(independent)
    B.require(admission['ready']==review['ready']==parents['IMPLEMENTATION_READY.json'],'Exact runtime authority')
    sources={'BUILD.py':B.pin(HERE/'BUILD.py'),'RUN.sbatch':B.pin(HERE/'RUN.sbatch'),'TEST_BUILD.py':B.pin(HERE/'TEST_BUILD.py'),'PREPARE.py':B.pin(__file__),
       'COMMON.py':B.pin(Path(h['inherited_recipe_authorities']['wrapper_builder']['path']).parent/'COMMON.py'),
       'allocation12':oldmanifest['helpers']['slurm12'],'loader':oldmanifest['helpers']['slurm1']}
    tools={name:B.pin(shutil.which(name)) for name in ('cmake','make','ar','ranlib','as','ld','ldd')}
    tools['scontrol']=B.pin('/usr/local/bin/scontrol')
    inputs={'status':'prepared_candidate_independent_v2_build_inputs','implementation':B.pin(HERE/'BUILD.py'),
       'declaration':B.pin(HERE/'ROOT_DECLARATION.json'),'handoff':parents['BUILD_HANDOFF.json'],'runtime_ready':parents['IMPLEMENTATION_READY.json'],
       'runtime_source_review':independent,'runtime_root_admission':root,'sources':sources,'tools':tools,'native_runtime':oldmanifest['native_runtime'],
       'candidate':None,'composition':None,'counts':{'compiler_calls':0,'native_calls':0,'scheduler_calls':0}}
    B.save(out/'INPUTS.json',inputs)
    cfg={'status':'unbound_periodic_v2_build_configuration','implementation':B.pin(HERE/'BUILD.py'),'inputs':B.pin(out/'INPUTS.json'),
       'runtime_source_review':independent,'runtime_root_admission':root,'archive_workers':12,'candidate':None,'composition':None}
    B.save(HERE/'CONFIG_TEMPLATE.json',cfg)
    release={'status':'unreleased_one_periodic_v2_build','implementation':cfg['implementation'],'configuration':None,'output':None,'build_source_review':None,
       'maximum_archive_builds':1,'maximum_wrapper_compiles':8,'maximum_wrapper_links':3,'maximum_native_programs':0}
    B.save(HERE/'RELEASE_TEMPLATE.json',release)
    print(json.dumps({'inputs':B.pin(out/'INPUTS.json'),'configuration_template':B.pin(HERE/'CONFIG_TEMPLATE.json')}))

if __name__=='__main__':
    import sys
    main(Path(sys.argv[1]).resolve())
