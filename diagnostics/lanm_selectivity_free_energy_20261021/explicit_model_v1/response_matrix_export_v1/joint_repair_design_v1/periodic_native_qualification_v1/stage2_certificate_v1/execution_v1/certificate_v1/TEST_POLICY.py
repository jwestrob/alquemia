"""Finite artificial/AST controls. No actual geometry/backend/solver invocation."""
import argparse
import ast
import copy
import tempfile
from unittest.mock import patch
from pathlib import Path
import CERTIFY as C

def run(out):
    tick=C.time.monotonic();cpu=C.time.process_time();checks=[];here=Path(__file__).resolve().parent
    def check(name,ok):
        C.require(ok,name);checks.append(name)
    original=here.parents[2]/'stage1_source_v1/REPLAY.py'
    old=ast.parse(original.read_text());new=ast.parse((here/'STATE.py').read_text())
    def node(tree,name):return next(x for x in tree.body if isinstance(x,(ast.FunctionDef,ast.ClassDef)) and x.name==name)
    for name in ['Rows','floats','compare','read_current']:check(name+' AST unchanged',ast.dump(node(old,name))==ast.dump(node(new,name)))
    expected=ast.get_source_segment(original.read_text(),node(old,'compact_response')).replace('def compact_response(sites,coefficients,gamma,pdamp):','def compact_response(sites,coefficients,gamma,pdamp,cell):').replace('d-=49.324*np.rint(d/49.324);r=np.linalg.norm(d,axis=1)','cell=np.asarray(cell,dtype=float)\n    if cell.shape!=(3,) or not np.isfinite(cell).all() or min(cell)<44:raise ValueError("Declared orthogonal cell")\n    d-=cell*np.rint(d/cell);r=np.linalg.norm(d,axis=1)')
    check('compact response changes explicit cell only',ast.dump(ast.parse(expected).body[0])==ast.dump(node(new,'compact_response')))
    state=C.load(C.pin(here/'STATE.py'),'tiny_stage2_state');np=state.np
    sites=[{'xyz':[0.,0.,0.],'alpha':1.,'type':900},{'xyz':[2.,0.,0.],'alpha':1.,'type':402},{'xyz':[0.,2.,0.],'alpha':1.,'type':402}]
    gamma=np.array([0.,.4,.4]);damp=np.ones(3);coef=[.1,.2,.3,.4]
    h,b=state.compact_response(sites,coef,gamma,damp,[49.324]*3)
    shifted=copy.deepcopy(sites);shifted[1]['xyz'][0]+=49.324
    hs,bs=state.compact_response(shifted,coef,gamma,damp,[49.324]*3)
    check('manufactured whole-cell image invariance',np.allclose(h,hs,rtol=0,atol=1e-12) and np.allclose(b,bs,rtol=0,atol=1e-12))
    hz,bz=state.compact_response(sites,[0.]*4,gamma,damp,[49.324]*3)
    check('manufactured zero response',np.array_equal(hz,np.zeros((3,3))) and np.array_equal(bz,np.zeros(3)))
    check('manufactured collective metal-only H',b[0]>0 and np.array_equal(b[1:],np.zeros(2)))
    rejected=False
    try:state.compact_response(sites,coef,gamma,damp,[40.]*3)
    except ValueError:rejected=True
    check('unsupported cell rejected',rejected)
    # A failed first child must leave no factor launch and preserve phase result.
    with tempfile.TemporaryDirectory() as td:
        p=Path(td);config=p/'config.json';release=p/'release.json';C.save(config,{});C.save(release,{})
        cfg={'state_id':'tiny','manifest':None,'build':None,'initializer':None}
        m={'candidate':None,'native_sources':None,'rows':{'tiny':{'source_inputs':{},'inputs':{}}}}
        init={'state_output':None,'pme_export':None}
        launched=[]
        def failed_invoke(r,o,label,*args):
            launched.append(label);r['commands'][label]={'status':'failed'};(o/label).mkdir();C.save(o/label/'RESULT.json',{'status':'manufactured_failure','error':'retained'});raise RuntimeError('manufactured first phase failure')
        with patch.object(C,'context',return_value=(cfg,m,{}, {},init)),patch.object(C,'invoke',side_effect=failed_invoke):
            code=C.supervise(config,release,p/'out')
        r=C.json.loads((p/'out/RESULT.json').read_text())
        check('first phase failure stops factor launch',code==2 and launched==['assembly'] and 'factor' not in r['commands'])
        check('failed phase pin and reason preserved',r['status']=='stage2_state_certificate_failed' and r['assembly']['sha256']==C.pin(p/'out/assembly/RESULT.json')['sha256'] and r['error'] is not None)
    source=(here/'CERTIFY.py').read_text();tree=ast.parse(source)
    calls=[n.func.attr for n in ast.walk(node(tree,'factor')) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)]
    check('one factor and one original sampled audit call',calls.count('factor_once')==1 and calls.count('audit_factor')==1)
    check('no new eigensolver or optimizer',not any(x in source for x in ['eigvalsh','minimize(','linprog(','cholesky(']))
    check('writable index recovery preserved',"np.array(idx,copy=True,order='C')" in source and "np.array(ptr,copy=True,order='C')" in source)
    check('phase limits fixed',C.CAPS=={'assembly':{'wall':180,'AS':8*1024**3,'CPU':180,'BLAS':1},'factor':{'wall':900,'AS':32*1024**3,'CPU':10800,'BLAS':12}})
    check('failed factor discarded before residual solves',source.index("if error is not None or ww:")<source.index('backend.audit_factor') and "del f;return" in source)
    check('no fixed old graph count reused','16371009' not in source and '1815015' not in source)
    check('fixed nominal helper reused','nominal_helper.nominal_parameters(a,a)' in source and "m['nominal_contract']['sha256']==nominal_helper.NOMINAL_SHA" in source)
    check('ten solve cap retained',"result['solve_attempts']<=10" in source and "backend.ten_probes(N,(0,1,2))" in source)
    check('all helper source pins checked',"review['certificate_sources']==cfg['certificate_sources']" in source)
    out.mkdir(parents=True,exist_ok=False)
    result={'status':'passed_stage2_certificate_source_controls','implementation':C.pin(__file__),'sources':{n:C.pin(here/n) for n in ['CERTIFY.py','STATE.py','BIND.py','RUN.sbatch']},'state_source_authority':C.pin(original),'checks':checks,'count':len(checks),'actual':{'wall_seconds':C.time.monotonic()-tick,'CPU_seconds':C.time.process_time()-cpu,'numerical_threads':1,'largest_artificial_response_sites':3,'native_calls':0,'assembly_calls':0,'symbolic_calls':0,'numeric_factor_calls':0,'solve_calls':0}}
    C.save(out/'RESULT.json',result);print(C.json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.out)
