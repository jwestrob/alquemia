"""Isolated OpenMP build for the existing native dense source-point loop."""
from pathlib import Path
import argparse
import difflib
import shutil
import tarfile
from affordable_common import InvalidArtifact,read_json,record,verify,write_new


def patch(text):
    start=text.index('subroutine build_phi_dense(')
    end=text.index('end subroutine build_phi_dense',start)
    body=text[start:end]
    before='    do icav = 1, ncav\n'
    after='\n    deallocate(vscales, vscales_rel, v4pi2lp1, stat=info)'
    if body.count(before)!=1 or body.count(after)!=1:
        raise InvalidArtifact('native dense source loop changed')
    body=body.replace(before,'    !$omp parallel do default(shared) private(icav,im,c,v) schedule(static)\n'+before)
    body=body.replace(after,'    !$omp end parallel do\n'+after)
    return text[:start]+body+text[end:]


def prepare(parent_manifest,plan,output):
    parent=read_json(parent_manifest)
    if parent['protocol']!='isolated_pyddx_0p9p0_portable_cmake_build_v2':
        raise InvalidArtifact('original portable stock build required')
    for pin in parent['pins']:verify(pin)
    archive=next(pin for pin in parent['pins'] if pin['path']==parent['source_archive'])
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);source=root/'source';source.mkdir()
    with tarfile.open(verify(archive)) as f:f.extractall(source,filter='data')
    tree=source/'pyddx-0.9.0';target=tree/'src/ddx_multipolar_solutes.f90';before=target.read_text();original=record(target)
    after=patch(before);target.write_text(after)
    (root/'source_parallel.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='upstream/ddx_multipolar_solutes.f90',tofile='parallel/ddx_multipolar_solutes.f90')))
    unchanged={}
    with tarfile.open(verify(archive)) as f:
        for p in sorted((tree/'src').iterdir()):
            if not p.is_file() or p==target:continue
            native=f.extractfile('pyddx-0.9.0/src/'+p.name).read()
            if native!=p.read_bytes():raise InvalidArtifact('unrelated native source changed')
            unchanged[p.name]=record(p)
    patched=root/'pyddx-0.9.0-dense-source-openmp.tar.gz'
    with tarfile.open(patched,'w:gz') as f:f.add(tree,arcname='pyddx-0.9.0')
    shutil.copyfile(__file__,root/Path(__file__).name)
    wrapper=Path(__file__).resolve().parents[1]/'diagnostics/mace_omol_20260917/run_ddx_build_v2.sbatch'
    m={**parent,'protocol':'isolated_pyddx_dense_source_openmp_build_v1','source_archive':str(patched),
       'parent_build_manifest':record(parent_manifest),'plan':record(plan),'builder':record(root/Path(__file__).name),
       'patch':record(root/'source_parallel.patch'),'original_native_source':original,'parallel_native_source':record(target),
       'unchanged_source_files':unchanged,'wrapper':record(wrapper),'scientific_calculations':0}
    # The original hash is historical, not a reference to the now-patched path.
    m['original_native_source']={'archive':archive,'member':'pyddx-0.9.0/src/ddx_multipolar_solutes.f90','sha256':original['sha256']}
    m['pins']=[*parent['pins'],record(patched),m['plan'],m['builder'],m['patch'],m['parallel_native_source'],m['wrapper'],*unchanged.values()]
    write_new(root/'manifest.json',m);return m


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('parent-build-manifest','plan','output'):p.add_argument('--'+k,required=True)
    a=p.parse_args();m=prepare(a.parent_build_manifest,a.plan,a.output)
    print({'protocol':m['protocol'],'unchanged_native_source_files':len(m['unchanged_source_files']),'scientific_calculations':0})
