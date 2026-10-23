"""One separately released isolated Stage2 build, no native program execution."""
import argparse
from pathlib import Path
import resource
import time
from COMMON import *

def run(config,release,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),'commands':{},'native_programs':0,'compile_attempts':0,'link_attempts':0,'error':None}
    save(out/'RESULT.json',result)
    try:
        cfg,m,rel,alloc=source_context(config,release,out,__file__,'build',1)
        result.update(manifest=cfg['manifest'],source_manifest=m['native_sources'],allocation=alloc)
        src=read(m['native_sources'])['files'];recipe=read(m['stage1_link_recipe'])
        for p in [recipe['library'],*recipe['module_files'],*recipe['compiler_and_link_dependencies']]:verify(p)
        compiler=recipe['argv'][0];flags=['-static-libgcc','-O2','-g','-fbacktrace','-fopenmp']
        require(recipe['argv'][1:6]==flags,'original compiler flags')
        include=next(x for x in recipe['argv'] if x.startswith('-I'))
        tail=recipe['argv'][recipe['argv'].index(recipe['library']['path']):-2]
        work=out/'build';work.mkdir();mods=work/'mod';mods.mkdir()
        objects={}
        for name in ('GATE.f90','STATE_WRITER.f90','induce_qualification.f','GRADIENT.f90','INITIALIZE.f90'):
            source=verify(src[name]);obj=work/(name+'.o');objects[name]=obj
            argv=[compiler,*flags,include,'-I'+str(mods),'-J'+str(mods),'-c',str(source),'-o',str(obj)]
            result['compile_attempts']+=1;require(result['compile_attempts']<=5,'fixed five compile commands')
            invoke(result,out,'compile_'+name.replace('.','_'),argv,work,180,3*1024**3,180)
        exe={}
        for name,main,clone in [('initialize','INITIALIZE.f90',True),('gradient_clone','GRADIENT.f90',True),('gradient_old','GRADIENT.f90',False)]:
            executable=work/name
            selected=[objects[main],objects['GATE.f90'],objects['STATE_WRITER.f90']]
            if clone:selected.append(objects['induce_qualification.f'])
            argv=[compiler,*flags,*map(str,selected),*tail,'-o',str(executable)]
            result['link_attempts']+=1;require(result['link_attempts']<=3,'fixed three links')
            invoke(result,out,'link_'+name,argv,work,180,3*1024**3,180)
            loader=load(m['helpers']['slurm1'],'stage2_loader')
            checked=subprocess.run([m['loader_command']['path'],str(executable)],check=True,capture_output=True,text=True,timeout=30)
            require(not checked.stderr.strip(),'loader stderr')
            runtime=loader.native_loader(checked.stdout,m['native_runtime'])
            exe[name]=pin(executable);result.setdefault('native_runtime',{})[name]=runtime
        result.update(status='stage2_build_complete',executables=exe,library=recipe['library'],compiler=pin(compiler))
    except BaseException as exc:result.update(status='stage2_build_failed',error=repr(exc))
    result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,'controller_peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1,'scheduler_cost_separate':True}
    save(out/'RESULT.json',result);return 0 if result['status']=='stage2_build_complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();raise SystemExit(run(a.config,a.release,a.out))
