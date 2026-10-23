"""Small key-rejection and source-identity checks; no compiler/native invocation."""
import ast
import argparse
from BINDING import *

def run(out,manifest):
    start=time.monotonic();cpu=time.process_time();checks=[]
    def check(name,ok):require(ok,name);checks.append(name)
    raw=b'parameters test.prm\newald-cutoff 9.0\nvdw-cutoff 9.0\newald-alpha 0.4\n'
    check('only Ewald cutoff token changes',cutoff_key(raw)==raw.replace(b'ewald-cutoff 9.0',b'ewald-cutoff 12.0'))
    check('line endings preserved',cutoff_key(raw.replace(b'\n',b'\r\n')).count(b'\r\n')==4)
    for name,bad in [('active response',raw+b'polar-local-response 900 402 1 2 3 4\n'),('duplicate cutoff',raw+b'ewald-cutoff 9.0\n'),('wrong cutoff',raw.replace(b'9.0\nvdw',b'12.0\nvdw')),('missing cutoff',b'parameters x\n')]:
        try:cutoff_key(bad)
        except ValueError:check('reject '+name,True)
        else:raise AssertionError(name)
    m=read(pin(manifest));parent,oldbuild=validate_manifest(m)
    check('nine metadata chains validate',len(m['rows'])==9)
    original=verify(parent['execution_sources']['NATIVE.py']).read_text();actual=(HERE/'RUN.py').read_text()
    fun=lambda t:[ast.dump(x,include_attributes=False) for x in ast.parse(t).body if isinstance(x,ast.FunctionDef)]
    check('all native functions AST identical',fun(original)==fun(actual))
    tree=ast.parse((HERE/'BUILD.py').read_text());calls=[x for x in ast.walk(tree) if isinstance(x,ast.Call) and isinstance(x.func,ast.Name) and x.func.id=='invoke']
    check('one compile and one link invocation',len(calls)==2 and {x.args[2].value for x in calls}=={'compile_GATE_f90','link_gradient_old'})
    tree=ast.parse((HERE/'BIND.py').read_text());check('binder has no execution calls',not any(isinstance(x,ast.Call) and ((isinstance(x.func,ast.Name) and x.func.id in ('invoke','exec','eval')) or (isinstance(x.func,ast.Attribute) and x.func.attr in ('Popen','run','system'))) for x in ast.walk(tree)))
    check('nine original XYZ and parameter identities',all(r['inputs'][k]==r['cutoff9_inputs'][k] for r in m['rows'].values() for k in ('xyz','parameters')))
    check('no initializer or factor in declared counts',m['counts']['new_initializers']==m['counts']['new_factors']==0)
    result={'status':'passed_cutoff12_source_controls','implementation':pin(__file__),'manifest':pin(manifest),'checks':checks,'check_count':len(checks),
            'compile_calls':0,'native_calls':0,'factor_calls':0,'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    out.mkdir(parents=True,exist_ok=False);save(out/'RESULT.json',result);print(json.dumps({'result':pin(out/'RESULT.json'),'checks':len(checks)}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);a=p.parse_args();run(a.out,a.manifest)
