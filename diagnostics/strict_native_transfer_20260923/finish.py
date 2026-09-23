"""Only terminal collection/comparison/reporting for the four declared jobs."""
import argparse,datetime,fcntl,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new

def finish(run):
    run=Path(run).resolve();sub=read_json(run/'SUBMISSION.json');jobs=[x['job_id'] for x in sub['jobs']]
    if len(jobs)!=4:raise RuntimeError('exact four submitted jobs required')
    with (run/'finish.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (run/'FINISH.json').exists():return read_json(run/'FINISH.json')
        impl=read_json(run/'shard_0/manifest.json')['implementation']['strict_native_transfer.py'];verify(impl)
        reporter=record(Path(__file__).with_name('report_transfer.py'))
        start=time.time();recovered=[]
        write_new(run/'OBSERVER.json',{'pid':os.getpid(),'started_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'submission':record(run/'SUBMISSION.json'),'implementation':record(__file__),'scorer':impl,'reporter':reporter,
            'scope':'wait exact4 jobs; collect missing terminal records once; compare/report/test only; no molecular calls or retries'})
        while True:
            q=subprocess.run(['squeue','-h','-j',','.join(jobs),'-o','%i'],capture_output=True,text=True)
            if q.returncode:raise RuntimeError('scheduler query failed: '+q.stderr)
            live=set(q.stdout.split())
            if not live:break
            time.sleep(45)
        commands=[]
        def execute(args,name):
            logfile=run/name
            with logfile.open('x') as f:r=subprocess.run(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
            commands.append({'command':args,'returncode':r.returncode,'log':record(logfile)})
            if r.returncode:raise RuntimeError('report-only operation failed: '+name)
        try:
            for i in range(4):
                dest=run/f'shard_{i}'/'COLLECTION.json'
                if not dest.exists():
                    execute([sys.executable,str(verify(impl)),'collect','--manifest',str(dest.parent/'manifest.json'),'--output',str(dest)],f'RECOVERY_COLLECT_{i}.log');recovered.append(i)
            target=run/'COMPARISON.json'
            if not target.exists():execute([sys.executable,str(verify(impl)),'compare','--inventory',str(run/'INVENTORY.json'),'--collections',*[str(run/f'shard_{i}'/'COLLECTION.json') for i in range(4)],'--output',str(target)],'FINAL_COMPARE.log')
            report=Path(__file__).with_name('REPORT.md')
            if not (run/'SUMMARY.json').exists():execute([sys.executable,str(verify(reporter)),'--comparison',str(target),'--submission',str(run/'SUBMISSION.json'),'--output',str(report)],'FINAL_REPORT.log')
            execute([sys.executable,'-m','unittest','discover','-s','tests','-p','test_strict_native_transfer.py','-v'],'FINAL_TESTS.log')
            result={'status':'complete','comparison':record(target),'summary':record(run/'SUMMARY.json'),'costs':record(run/'COSTS.json'),'report':record(report),'tests':record(run/'FINAL_TESTS.log')}
        except Exception as exc:result={'status':'reporting_failed','reason':str(exc)}
        result.update(commands=commands,recovered_collection_shards=recovered,new_molecular_calls=0,observer_wall_seconds=time.time()-start)
        write_new(run/'FINISH.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',required=True,type=Path)
    print(json.dumps(finish(**vars(p.parse_args())),indent=2),flush=True)
