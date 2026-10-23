"""One GATE compile plus one old-library link; no native program execution."""
import argparse
from BINDING import *

def run(config,release,out):
    out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'started','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),'commands':{},
            'native_programs':0,'compile_attempts':0,'link_attempts':0,'error':None}
    save(out/'RESULT.json',result)
    try:
        cfg,m,rel,alloc=source_context(config,release,out,__file__,'build',1)
        result.update(manifest=cfg['manifest'],source_manifest=m['native_sources'],allocation=alloc,reused_objects=m['reused_objects'],parent_build=m['parent_build'])
        src=read(m['native_sources'])['files'];recipe=read(m['stage1_link_recipe'])
        for item in [recipe['library'],*recipe['module_files'],*recipe['compiler_and_link_dependencies']]:verify(item)
        compiler=recipe['argv'][0];flags=['-static-libgcc','-O2','-g','-fbacktrace','-fopenmp']
        require(recipe['argv'][1:6]==flags,'Original compiler flags')
        include=next(x for x in recipe['argv'] if x.startswith('-I'));tail=recipe['argv'][recipe['argv'].index(recipe['library']['path']):-2]
        work=out/'build';work.mkdir();mods=work/'mod';mods.mkdir();obj=work/'GATE.f90.o'
        argv=[compiler,*flags,include,'-I'+str(mods),'-J'+str(mods),'-c',str(verify(src['GATE.f90'])),'-o',str(obj)]
        result['compile_attempts']=1;invoke(result,out,'compile_GATE_f90',argv,work,180,3*1024**3,180)
        main=verify(m['reused_objects']['GRADIENT.f90']);writer=verify(m['reused_objects']['STATE_WRITER.f90']);exe=work/'gradient_old'
        argv=[compiler,*flags,str(main),str(obj),str(writer),*tail,'-o',str(exe)]
        result['link_attempts']=1;invoke(result,out,'link_gradient_old',argv,work,180,3*1024**3,180)
        loader=load(m['helpers']['slurm1'],'cutoff12_loader')
        ld=subprocess.run([m['loader_command']['path'],str(exe)],check=True,capture_output=True,text=True,timeout=30)
        require(not ld.stderr.strip(),'Loader stderr')
        result.update(status='stage2_build_complete',executables={'gradient_old':pin(exe)},library=recipe['library'],compiler=pin(compiler),
                      native_runtime={'gradient_old':loader.native_loader(ld.stdout,m['native_runtime'])})
    except BaseException as exc:result.update(status='stage2_build_failed',error=repr(exc))
    result['cost']={'wall_seconds':time.monotonic()-tick,'controller_CPU_seconds':time.process_time()-cpu,
                    'controller_peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1,'scheduler_cost_separate':True}
    save(out/'RESULT.json',result);return 0 if result['status']=='stage2_build_complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();raise SystemExit(run(a.config,a.release,a.out))
