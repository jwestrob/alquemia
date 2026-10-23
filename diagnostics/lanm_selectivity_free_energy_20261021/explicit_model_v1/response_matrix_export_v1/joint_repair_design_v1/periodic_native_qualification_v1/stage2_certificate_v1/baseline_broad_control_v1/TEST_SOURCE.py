"""Finite source/key guards only; never invokes native, scheduler or scientific helpers."""
import ast
import copy
from BINDING import *

def test(manifest,out):
    tick=time.monotonic();cpu=time.process_time();checks=[]
    def check(label,condition):require(condition,label);checks.append(label)
    def rejects(label,fn):
        try:fn()
        except ValueError:checks.append(label);return
        raise ValueError('Did not reject '+label)
    raw=b'# retained\nPOLAR-LOCAL-RESPONSE 900 402 .1 .2 .3 .4\ncutoff 9\npolar-local-response 901 402 .1 .2 .3 .4\n'
    check('Only exact activation lines removed',disabled_key(raw)==b'# retained\ncutoff 9\n')
    rejects('Missing activation rejected',lambda:disabled_key(b'cutoff 9\n'))
    rejects('Duplicate activation rejected',lambda:disabled_key(raw.replace(b'901',b'900')))
    rejects('Wrong namespace rejected',lambda:disabled_key(raw.replace(b'901',b'902')))
    rejects('Wrong water role rejected',lambda:disabled_key(raw.replace(b'402',b'403')))
    old=(PARENT/'NATIVE.py').read_text();new=(HERE/'RUN.py').read_text()
    check('Exact single import delta',new==old.replace('from COMMON import *','from BINDING import *'))
    functions=lambda text:{x.name:ast.dump(x,include_attributes=False) for x in ast.parse(text).body if isinstance(x,ast.FunctionDef)}
    check('All native functions AST identical',functions(old)==functions(new))
    m=read(pin(manifest));validate_manifest(m);checks.append('Actual eight metadata bindings and deletion joins')
    wrong=copy.deepcopy(m);wrong['rows'].pop(IDS[-1]);rejects('Missing slot rejected',lambda:validate_manifest(wrong))
    wrong=copy.deepcopy(m);wrong['rows'][IDS[0]]['role']='ACTIVE';rejects('Active role rejected',lambda:validate_manifest(wrong))
    wrong=copy.deepcopy(m);wrong['rows'][IDS[0]]['inputs']['key']=wrong['rows'][IDS[0]]['parent_active_inputs']['key']
    rejects('Unremoved response rejected',lambda:validate_manifest(wrong))
    check('Original invoke helper retained',invoke.__module__=='COMMON' and pin(PARENT/'COMMON.py')['sha256']==COMMON_SHA)
    check('All fixed outputs distinct',len({r['output'] for r in m['rows'].values()})==8)
    out.mkdir(parents=True,exist_ok=False)
    save(out/'RESULT.json',{'status':'passed_focused_baseline_source_controls','implementation':pin(__file__),
         'manifest':pin(manifest),'checks':checks,'count':len(checks),'source':m['execution_sources'],
         'native_calls':0,'scheduler_calls':0,'initializers':0,'factors':0,
         'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}})
    print(json.dumps({'result':pin(out/'RESULT.json'),'passed':len(checks),'native_calls':0}))

if __name__=='__main__':test(Path(sys.argv[1]),Path(sys.argv[2]))
