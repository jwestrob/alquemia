#!/usr/bin/env python3
"""Single released link and OFF/ON pair. No retry or native fallback."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import socket
import subprocess
import sys
import time


def pin(path):
    path=Path(path).resolve()
    return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def verify(item):
    if pin(item['path'])!=item:raise ValueError('Changed pin: '+item['path'])
    return Path(item['path'])


def main(config_path,out_path):
    start,cpu=time.monotonic(),time.process_time()
    out=Path(out_path).resolve();out.mkdir(parents=True,exist_ok=False)
    r=dict(status='failed',error=None,configuration=pin(config_path),implementation=pin(__file__),
           execution_mode='login',job_id=None,host=socket.gethostname(),PID=os.getpid(),
           started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
           native_program_invocations=0,link_invocations=0,
           cells={m:dict(status='not_attempted') for m in ['OFF','ON']})
    def save():
        (out/'RESULT.json').write_text(json.dumps(r,indent=2)+'\n')
    def command(argv,cwd,stem,kind):
        rec=dict(status='starting',argv=list(map(str,argv)),cwd=str(cwd),wall_limit_seconds=120)
        r[stem]=rec;save();wall=time.monotonic();usage=resource.getrusage(resource.RUSAGE_CHILDREN)
        proc=None
        try:
            with (out/(stem+'.stdout')).open('w') as stdout,(out/(stem+'.stderr')).open('w') as stderr:
                proc=subprocess.Popen(rec['argv'],cwd=cwd,stdout=stdout,stderr=stderr,start_new_session=True)
                rec.update(status='running',PID=proc.pid)
                if kind=='native':r['native_program_invocations']+=1
                if kind=='link':r['link_invocations']+=1
                save()
                try:rec['returncode']=proc.wait(timeout=120)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL);proc.wait();rec['returncode']=proc.returncode
                    raise RuntimeError(stem+' exceeded its declared 120 s wall limit')
            if rec['returncode']!=0:raise RuntimeError(stem+' returned '+str(rec['returncode']))
            rec['status']='complete'
        except Exception as error:
            rec.update(status='failed',error=repr(error));raise
        finally:
            end=resource.getrusage(resource.RUSAGE_CHILDREN)
            rec.update(wall_seconds=time.monotonic()-wall,child_CPU_seconds=end.ru_utime+end.ru_stime-usage.ru_utime-usage.ru_stime,
                       stdout=pin(out/(stem+'.stdout')) if (out/(stem+'.stdout')).exists() else None,
                       stderr=pin(out/(stem+'.stderr')) if (out/(stem+'.stderr')).exists() else None)
            save()
    try:
        cfg=json.loads(Path(config_path).read_text())
        for k in ['root_release','root_source_review','ready','manifest','harness','replay','build_result','library','source_manifest']:
            verify(cfg[k])
        for item in cfg['modules']:verify(item)
        source_review=json.loads(Path(cfg['root_source_review']['path']).read_text())
        if source_review['status']!='passed_initialized_PME_qualification_source_review' or any(source_review[k]!=cfg[k] for k in ['ready','manifest','harness','replay']):
            raise ValueError('Source review does not bind this qualification')
        release=json.loads(Path(cfg['root_release']['path']).read_text())
        if release['native_program_cap']!=2 or release['link_cap']!=1 or release['configuration_sources']!=cfg['configuration_sources']:
            raise ValueError('Release mismatch')
        if os.environ.get('SLURM_JOB_ID') or os.environ.get('SLURM_STEP_ID') or socket.gethostname()!='biotite':
            raise ValueError('Wrong execution envelope')
        if any(os.environ.get(k)!='1' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']):
            raise ValueError('One-thread envelope required')
        build=json.loads(Path(cfg['build_result']['path']).read_text())
        if build['status']!='complete' or build['library']!=cfg['library'] or build['manifest']!=cfg['source_manifest']:
            raise ValueError('Build not admitted')
        m=json.loads(Path(cfg['manifest']['path']).read_text())
        if m['native_program_cap']!=2 or m['gradient_call_cap']!=2 or m['internal_reciprocal_applications']!=54:
            raise ValueError('Changed finite inventory')
        r.update(root_release=cfg['root_release'],build_result=cfg['build_result'],source_manifest=cfg['source_manifest'],
                 library=cfg['library'],modules=cfg['modules'],manifest=cfg['manifest'],root_source_review=cfg['root_source_review'])
        exe=out/'pme_harness'
        build_dir=Path(cfg['library']['path']).parent
        link=['/usr/bin/f95','-static-libgcc','-O2','-g','-fbacktrace','-fopenmp','-I'+str(build_dir/'mod'),cfg['harness']['path'],
              cfg['library']['path'],'/usr/lib/x86_64-linux-gnu/libfftw3.a','/usr/lib/x86_64-linux-gnu/libfftw3_threads.a',
              '/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so','/usr/lib/x86_64-linux-gnu/libpthread.a','-o',str(exe)]
        command(link,out,'link','link')
        r['harness_executable']=pin(exe);save()
        for mode in ['OFF','ON']:
            cell=out/mode;cell.mkdir()
            for key,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]:
                shutil.copy2(verify(m['cells'][mode][key]),cell/name)
            r['cells'][mode].update(status='prepared',inputs={name:pin(cell/name) for name in ['system.xyz','system.key','amoeba09.prm']});save()
            try:
                verify(r['harness_executable'])
                command([exe,'system.xyz',mode],cell,mode,'native')
                if not(cell/'harness_ef.dat').read_text().endswith('END 6\n'):raise ValueError('Incomplete EF')
                if mode=='ON' and (cell/'harness_apply.dat').read_text().splitlines()[-1].split()!=['END','54']:
                    raise ValueError('Incomplete applications')
                r['cells'][mode]['status']='complete'
            except Exception as error:
                r['cells'][mode].update(status='failed',error=repr(error));raise
            finally:save()
        command([sys.executable,cfg['replay']['path'],cfg['manifest']['path'],out/'OFF',out/'ON',out/'replay'],out,'replay','saved')
        r['saved_replay']=pin(out/'replay/RESULT.json')
        rr=json.loads((out/'replay/RESULT.json').read_text())
        if rr['status']!='passed_initialized_PME_export_numerical_qualification':raise ValueError('Replay not passed')
        r.update(status='initialized_PME_export_qualification_complete',gradient_calls_completed=2,internal_reciprocal_applications_completed=54)
    except Exception as error:r['error']=repr(error)
    finally:
        usage=resource.getrusage(resource.RUSAGE_CHILDREN)
        r.update(wall_seconds=time.monotonic()-start,controller_CPU_seconds=time.process_time()-cpu,
                 child_CPU_seconds=usage.ru_utime+usage.ru_stime,peak_child_RSS_KiB=usage.ru_maxrss,
                 scheduler_allocated_CPU_seconds=0,maximum_native_workers=1)
        # Record completed partial evidence without inventing unprinted work.
        r['observed_complete_gradient_files']=sum((out/m/'harness_ef.dat').is_file() and (out/m/'harness_ef.dat').read_text().endswith('END 6\n') for m in ['OFF','ON'])
        save()
    print(json.dumps({k:r[k] for k in ['status','error','native_program_invocations','link_invocations','wall_seconds']}))
    return 0 if r['status']=='initialized_PME_export_qualification_complete' else 1


if __name__=='__main__':raise SystemExit(main(*sys.argv[1:]))
