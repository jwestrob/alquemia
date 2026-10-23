"""One fixed cutoff12 baseline build and nine separately released gradient states."""
import hashlib
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'execution_v1'
COMMON_SHA='b13a4d20fbd9cc0f17db35b72244c2f9cd4d0d263b65466a7bce07e47decdb6e'
if hashlib.sha256((PARENT/'COMMON.py').read_bytes()).hexdigest()!=COMMON_SHA:
    raise ValueError('Changed original runtime helpers')
sys.path.insert(0,str(PARENT))
from COMMON import *

DECLARATION_SHA='e7246ed093197b200fc096ecbbdb8e4cc03b57db0603fa06d74b1e2a074ae36a'
IDS=('central',)+tuple(d+'__'+h+'__'+s for d in ('whole_box','strain_1') for h in ('1e-04','5e-05') for s in ('minus','plus'))

def cutoff_key(raw):
    lines=raw.splitlines(keepends=True);count=0;out=[]
    for line in lines:
        words=line.split()
        require(not words or words[0].lower()!=b'polar-local-response','Response must remain absent')
        if words and words[0].lower()==b'ewald-cutoff':
            require(words==[b'ewald-cutoff',b'9.0'],'Exact old cutoff line')
            require(line.count(b'9.0')==1,'Unambiguous cutoff token')
            line=line.replace(b'9.0',b'12.0');count+=1
        out.append(line)
    require(count==1,'Exactly one Ewald cutoff record')
    return b''.join(out)

def validate_manifest(m):
    require(m['status']=='prepared_nine_cutoff12_baseline_inputs' and tuple(m['rows'])==IDS and m['maximum_native_programs']==9,'Exact nine states')
    require(m['root_declaration']['sha256']==DECLARATION_SHA,'Prospective cutoff declaration')
    decl=read(m['root_declaration']);require(tuple(decl['states'])==IDS and decl['maximum_new_native_gradient_programs']==9,'Declared finite panel')
    base=read(m['baseline_manifest']);central=read(decl['central_baseline_result']);parent=read(m['parent_manifest']);oldbuild=read(m['parent_build'])
    require(decl['baseline_inputs_manifest']==m['baseline_manifest'] and oldbuild['status']=='stage2_build_complete'
            and oldbuild['manifest']==m['parent_manifest'] and oldbuild['source_manifest']==parent['native_sources'],'Admitted original source/build')
    require(read(decl['evidence'])['result']==m['baseline_failed_collection'],'Preserved cutoff9 failure')
    verify(m['baseline_failed_collection'])
    for item in m['execution_sources'].values():verify(item)
    require(m['execution_sources']['RUN.py']==pin(HERE/'RUN.py') and m['execution_sources']['BINDING.py']==pin(__file__),'Actual runtime sources')
    for key in ('helpers','collector','environment','native_runtime','loader_command','candidate','composition','stage1_link_recipe','directions','displaced_manifest'):
        require(m[key]==parent[key],'Unchanged inherited authority '+key)
    for item in m['helpers'].values():verify(item)
    verify(m['collector'])
    native=read(m['native_sources']);original=read(parent['native_sources'])
    old=verify(original['files']['GATE.f90']).read_bytes();new=verify(native['files']['GATE.f90']).read_bytes()
    require(old.count(b'ewaldcut/=9d0')==1 and new==old.replace(b'ewaldcut/=9d0',b'ewaldcut/=12d0'),'Only native harness cutoff guard changes')
    for name in ('GRADIENT.f90','STATE_WRITER.f90'):
        require(native['files'][name]==original['files'][name],'Unchanged helper source '+name)
        verify(native['files'][name]);obj=m['reused_objects'][name];verify(obj)
        command=oldbuild['commands']['compile_'+name.replace('.','_')]
        require(command['returncode']==0 and command['argv'][command['argv'].index('-o')+1]==obj['path']
                and command['argv'][command['argv'].index('-c')+1]==native['files'][name]['path'],'Actual old object/source join')
    require(len({r['output'] for r in m['rows'].values()})==9,'Nine immutable outputs')
    for sid,row in m['rows'].items():
        oldrow=central if sid=='central' else base['rows'][sid]
        require(row['role']=='BYPASS_OLD_ABSENT' and row['cutoff9_inputs']==oldrow['inputs'],'Exact inactive old input authority')
        for label in ('xyz','parameters'):
            require(row['source_inputs'][label]==row['inputs'][label]==oldrow['inputs'][label],'Literal unchanged '+label)
            verify(row['inputs'][label])
        require(verify(row['inputs']['key']).read_bytes()==cutoff_key(verify(oldrow['inputs']['key']).read_bytes()),'Runtime cutoff-only delta')
        require(verify(row['source_inputs']['key']).read_bytes()==cutoff_key(verify(oldrow['source_inputs']['key']).read_bytes()),'Prepared cutoff-only delta')
        require(verify(row['inputs']['key']).read_bytes()==verify(row['source_inputs']['key']).read_bytes()+b'\nPME-NORM-EXPORT 1\n','Unchanged observer suffix')
    return parent,oldbuild

def source_context(config,release,out,implementation,kind,cpus,blas=1):
    cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text())
    require(kind in ('build','native') and cpus==blas==1,'One CPU finite phase')
    require(cfg['status']=='bound_cutoff12_'+kind+'_configuration' and cfg['implementation']==rel['implementation']==pin(implementation),'Bound exact phase source')
    require(rel['status']=='released_cutoff12_'+kind and rel['configuration']==pin(config) and Path(rel['output']).resolve()==Path(out).resolve(),'Exact root output/release')
    m=read(cfg['manifest']);validate_manifest(m)
    review=read(rel['source_review'])
    require(review['status']=='passed_cutoff12_baseline_source_review' and review['manifest']==cfg['manifest']
            and review['binding']==pin(__file__) and review['implementations'][kind]==cfg['implementation'],'Exact source review')
    env=read(m['environment']);require(Path(sys.prefix).resolve()==Path(env['prefix']).resolve(),'Pinned isolated interpreter')
    for item in env['runtime_files']:verify(item)
    require(all(os.environ.get(k)=='1' for k in THREADS),'One-thread environment')
    if kind=='build':
        require(Path(out).resolve()==Path(m['build_output']).resolve(),'One immutable build output')
    else:
        require(cfg['operation']=='gradient' and cfg['state_id'] in IDS and cfg['certificate'] is None and cfg['certificate_admission'] is None,'Inactive gradient only')
        require(Path(out).resolve()==Path(m['rows'][cfg['state_id']]['output']).resolve(),'Fixed state output')
        build=native_build(cfg,m);admission=read(cfg['build_root_review'])
        require(admission['status']=='admitted_cutoff12_baseline_build' and admission['result']==cfg['build']
                and admission['manifest']==cfg['manifest'],'Actual separately reviewed build')
    return cfg,m,rel,allocation(m,1)

def native_build(cfg,m):
    b=read(cfg['build'])
    require(b['status']=='stage2_build_complete' and b['manifest']==cfg['manifest'] and b['source_manifest']==m['native_sources']
            and b['compile_attempts']==b['link_attempts']==1 and b['native_programs']==0 and b['error'] is None,'Exact successful one-compile/one-link build')
    require(set(b['executables'])=={'gradient_old'} and b['library']==read(m['stage1_link_recipe'])['library'],'Unchanged old model archive')
    verify(b['executables']['gradient_old']);return b
