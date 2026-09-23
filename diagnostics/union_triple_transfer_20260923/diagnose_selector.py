"""Read-only projection replay on the CPU host that rejected exact metadata."""
import json,os,sys,time
from pathlib import Path
from affordable_common import read_json,record,verify,write_new
import union_triple_transfer_run as runner


def differences(a,b,path=''):
    out=[]
    if type(a)!=type(b):return [{'path':path,'expected':str(type(a)),'actual':str(type(b)),'kind':'type'}]
    if isinstance(a,dict):
        for k in sorted(a.keys()|b.keys()):
            if k not in a or k not in b:out.append({'path':path+'.'+k,'kind':'key'})
            else:out.extend(differences(a[k],b[k],path+'.'+k))
    elif isinstance(a,list):
        if len(a)!=len(b):out.append({'path':path,'kind':'length'})
        else:
            for i,(x,y) in enumerate(zip(a,b)):out.extend(differences(x,y,path+'['+str(i)+']'))
    elif a!=b:out.append({'path':path,'expected':a,'actual':b,'kind':'float' if isinstance(a,float) else 'discrete',
                        'absolute_difference':abs(a-b) if isinstance(a,float) else None})
    return out


start=time.monotonic();result=[]
for path in sys.argv[2:]:
    m=read_json(path)
    for c in m['cases']:
        ts,choice=runner.build_tasks(c['origin_row'],m['model'])
        old=[t for t in m['tasks'] if t['case_id']==c['case_id']]
        result.append({'case_id':c['case_id'],'manifest':record(path),
            'same_selected_IDs':[x['id'] for x in choice['selected']]==[x['id'] for x in c['selection']['selected']],
            'selection_differences':differences(c['selection'],choice),'task_differences':differences(old,ts)})
write_new(sys.argv[1],{'host':os.uname().nodename,'job_id':os.environ.get('SLURM_JOB_ID'),'rows':result,
    'new_molecular_calls':0,'wall_seconds':time.monotonic()-start,'implementation':record(__file__)})
