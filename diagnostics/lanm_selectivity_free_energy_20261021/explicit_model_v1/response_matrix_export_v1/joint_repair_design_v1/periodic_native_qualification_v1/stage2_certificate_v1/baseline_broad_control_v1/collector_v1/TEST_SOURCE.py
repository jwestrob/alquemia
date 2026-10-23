"""Small manufactured arithmetic/accounting controls; no actual endpoint reads."""
import ast
from COLLECT import *

def test(out):
    tick=time.monotonic();cpu=time.process_time();checks=[]
    def check(name,ok):require(ok,name);checks.append(name)
    def rejects(name,fn):
        try:fn()
        except ValueError:checks.append(name);return
        raise ValueError('Did not reject '+name)
    inp=read(pin(HERE/'INPUTS.json'));mod=parent(inp['parent_collector'])
    check('Frozen original numerical gate',mod.GATE==Decimal('.001') and mod.STEPS==(Decimal('.0001'),Decimal('.00005')))
    def synthetic(slope,polslope):
        rows={}
        for h in mod.STEPS:
            rows[h]=tuple({'energy_token':str(Decimal(20)+Decimal(sign)*h*Decimal(slope)),
                           'polarization_token':str(Decimal(5)+Decimal(sign)*h*Decimal(polslope))} for sign in (-1,1))
        return rows
    baseline=synthetic(2,1);active=synthetic(5,4)
    slopes=component_slopes(mod,active)
    check('Total/polar/nonpolar components',all(x['components']['total']['FD']==5 and x['components']['polarization']['FD']==4 and x['components']['nonpolar']['FD']==1 for x in slopes['steps']))
    inc=increment(mod,active,baseline,Decimal(5),Decimal(2),Decimal('1e-12'),Decimal('2e-12'))
    check('Exact fixed-direction increment',all(x['FD_increment']==x['analytic_increment']==3 and x['absolute_residual']==0 for x in inc['steps']))
    check('Increment descriptive only','passed' not in inc and all('passed' not in x for x in inc['steps']))
    check('Printing uncertainties retained',all(x['printing_only_bound']>0 for x in inc['steps']))
    result=mod.slope_checks('synthetic',Decimal(2),Decimal(0),{h:tuple(x['energy_token'] for x in pair) for h,pair in baseline.items()})
    check('Original slope check reused',result['passed'] and result['two_step']['passed'])
    subs=[{'job_id':str(100+i)} for i in range(8)]
    raw='\n'.join(str(100+i)+'|test|COMPLETED|0:0|3|1|00:02.0' for i in range(8))
    saved=scheduler(raw+'\n100.batch|batch|COMPLETED|0:0|3|1|00:02.0',subs)
    check('Parent allocations counted once',saved['allocated_CPU_seconds']==24 and len(saved['records'])==9)
    rejects('Missing denominator rejected',lambda:scheduler('\n'.join(raw.splitlines()[:-1]),subs))
    rejects('Duplicate parent rejected',lambda:scheduler(raw+'\n'+raw.splitlines()[0],subs))
    rejects('Unrelated job rejected',lambda:scheduler(raw+'\n999|x|COMPLETED|0:0|3|1|00:02.0',subs))
    rejects('Failed job preserved',lambda:scheduler(raw.replace('100|test|COMPLETED','100|test|FAILED'),subs))
    rejects('Changed CPU count rejected',lambda:scheduler(raw.replace('|3|1|','|3|2|'),subs))
    tree=ast.parse((HERE/'COLLECT.py').read_text())
    check('No execution imports',not any(isinstance(x,(ast.Import,ast.ImportFrom)) and any(n.name.split('.')[0] in ('subprocess','scipy') for n in x.names) for x in ast.walk(tree)))
    check('Exact eight declared states',len(IDS)==8 and len(set(IDS))==8)
    out.mkdir(parents=True,exist_ok=False)
    save(out/'RESULT.json',{'status':'passed_focused_baseline_collector_controls','implementation':pin(__file__),'collector':pin(HERE/'COLLECT.py'),
         'inputs':pin(HERE/'INPUTS.json'),'checks':checks,'count':len(checks),'actual_endpoint_reads':0,'native_calls':0,'FFT_calls':0,'factors':0,'scheduler_queries':0,
         'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}})
    print(json.dumps({'result':pin(out/'RESULT.json'),'passed':len(checks)}))

if __name__=='__main__':test(Path(sys.argv[1]))
