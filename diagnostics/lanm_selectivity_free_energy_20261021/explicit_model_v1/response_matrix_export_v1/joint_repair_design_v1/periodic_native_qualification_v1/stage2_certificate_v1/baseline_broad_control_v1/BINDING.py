"""Eight response-disabled inputs; original native run and runtime helpers retained."""
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / 'execution_v1'
COMMON_SHA = 'b13a4d20fbd9cc0f17db35b72244c2f9cd4d0d263b65466a7bce07e47decdb6e'
if hashlib.sha256((PARENT/'COMMON.py').read_bytes()).hexdigest() != COMMON_SHA:
    raise ValueError('Changed original runtime helper')
sys.path.insert(0, str(PARENT))
from COMMON import *

DECLARATION_SHA = 'e09055bd32bbc417ce457d52388cc719b2404e4696be650f3e5ff0dc4a18e25c'
IDS = tuple(direction+'__'+step+'__'+sign
            for direction in ('whole_box','strain_1')
            for step in ('1e-04','5e-05') for sign in ('minus','plus'))

def disabled_key(raw):
    """Delete exactly the two original activation lines; preserve every other byte."""
    kept=[];seen=[]
    for line in raw.splitlines(keepends=True):
        words=line.split()
        if words and words[0].lower()==b'polar-local-response':
            require(len(words)==7 and words[1] in (b'900',b'901') and words[2]==b'402',
                    'Exact two water response records')
            seen.append(words[1])
        else:
            kept.append(line)
    require(sorted(seen)==[b'900',b'901'],'Two unique activation records required')
    return b''.join(kept)

def validate_manifest(m):
    require(m['status']=='prepared_eight_response_disabled_baseline_inputs','Prepared eight-slot manifest')
    require(tuple(m['rows'])==IDS and m['maximum_native_programs']==8,'Exact eight slots/cap')
    require(m['root_declaration']['sha256']==DECLARATION_SHA,'Prospective root declaration')
    declaration=read(m['root_declaration'])
    require(declaration['directions']==['whole_box','strain_1'] and declaration['maximum_new_native_programs']==8
            and declaration['half_steps_A']==[1e-4,5e-5] and declaration['signs']==[-1,1], 'Declared directions/steps')
    parent=read(m['parent_manifest']);build=read(m['build'])
    require(parent['status']=='prepared_stage2_execution_inputs_no_numerical_calls','Original manifest')
    require(build['status']=='stage2_build_complete' and build['manifest']==m['parent_manifest']
            and build['source_manifest']==parent['native_sources']==m['native_sources'],'Original admitted build bridge')
    for item in m['execution_sources'].values():verify(item)
    require(m['execution_sources']['RUN.py']==pin(HERE/'RUN.py') and m['execution_sources']['BINDING.py']==pin(__file__),
            'Actual adapter sources')
    require(m['parent_execution_sources']==parent['execution_sources'],'Original execution source inventory')
    for item in m['parent_execution_sources'].values():verify(item)
    for key in ('helpers','collector','environment','native_runtime','loader_command','candidate','composition',
                'native_sources','native_source_review','stage1_actual_admission','directions','displaced_manifest'):
        require(m[key]==parent[key],'Unchanged authority '+key)
    for item in m['helpers'].values():verify(item)
    verify(m['collector'])
    central=read(m['central_analytic_reuse']);scout=read(m['central_scout_review'])
    require(m['central_analytic_reuse']==declaration['central_analytic_reuse'] and
            central['state_id']=='BYPASS_OLD_ABSENT' and central['role']=='BYPASS_OLD_ABSENT' and
            central['status']=='stage2_gradient_complete' and central['build']==m['build'] and
            central['manifest']==m['parent_manifest'],'Exact admitted old baseline center')
    require(scout['status']=='admitted_stage2_central_and_bypass_scout' and scout['build']==m['build'] and
            scout['manifest']==m['parent_manifest'] and scout['all_three_bypass_tokens_equal'] is True,
            'Successful original central bypass')
    failure=read(declaration['evidence'])
    require(failure['result']==m['active_full_result'],'Preserved failed active panel')
    verify(m['active_full_result'])
    for sid,row in m['rows'].items():
        old=parent['rows'][sid]
        require(row['role']=='BYPASS_OLD_ABSENT' and row['parent_active_inputs']==old['inputs'],
                'Inactive old executable for exact former active slot')
        for label in ('xyz','parameters'):
            require(row['inputs'][label]==row['source_inputs'][label]==old['inputs'][label],
                    'Literal unchanged '+label)
            verify(row['inputs'][label])
        require(verify(row['inputs']['key']).read_bytes()==disabled_key(verify(old['inputs']['key']).read_bytes()),
                'Runtime key changes only two activation lines')
        require(verify(row['source_inputs']['key']).read_bytes()==disabled_key(verify(old['source_inputs']['key']).read_bytes()),
                'Prepared key changes only two activation lines')
        require(verify(row['inputs']['key']).read_bytes()==verify(row['source_inputs']['key']).read_bytes()+b'\nPME-NORM-EXPORT 1\n',
                'Original observer-only suffix')
    return parent,build

def source_context(config,release,out,implementation,kind,cpus,blas=1):
    cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text())
    require(kind=='native' and cpus==blas==1,'One native CPU; no build/initializer')
    require(cfg['status']=='bound_stage2_baseline_control_configuration' and cfg['operation']=='gradient',
            'Bound baseline gradient only')
    require(cfg['implementation']==rel['implementation']==pin(implementation),'Actual adapter source')
    require(rel['status']=='released_stage2_baseline_control' and rel['configuration']==pin(config),'Exact root release')
    require(Path(rel['output']).resolve()==Path(out).resolve(),'Released immutable output')
    m=read(cfg['manifest']);parent,build=validate_manifest(m);sid=cfg['state_id']
    require(sid in IDS and cfg['build']==m['build'] and cfg['central_scout_review']==m['central_scout_review'],
            'Exact row/build/scout')
    require(cfg['certificate'] is None and cfg['certificate_admission'] is None,'Inactive branch has no active certificate')
    require(Path(out).resolve()==Path(m['rows'][sid]['output']).resolve(),'One fixed output per original slot')
    review=read(rel['source_review'])
    require(review['status']=='passed_stage2_baseline_broad_control_source_review' and review['manifest']==cfg['manifest']
            and review['implementation']==cfg['implementation'] and review['binding']==pin(__file__),'Exact source review')
    env=read(m['environment']);require(Path(sys.prefix).resolve()==Path(env['prefix']).resolve(),'Pinned isolated interpreter')
    for item in env['runtime_files']:verify(item)
    require(all(os.environ.get(k)=='1' for k in THREADS),'One-thread environment')
    return cfg,m,rel,allocation(m,1)

def native_build(cfg,manifest):
    build=read(cfg['build'])
    require(cfg['build']==manifest['build'] and build['status']=='stage2_build_complete'
            and build['manifest']==manifest['parent_manifest'] and build['source_manifest']==manifest['native_sources'],
            'Explicit unchanged old-build lineage')
    verify(build['executables']['gradient_old'])
    return build
