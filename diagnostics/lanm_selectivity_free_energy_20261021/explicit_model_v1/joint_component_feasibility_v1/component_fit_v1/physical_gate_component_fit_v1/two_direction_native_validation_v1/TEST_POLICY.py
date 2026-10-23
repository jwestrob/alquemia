"""Changed-wrapper controls only, with native/solver paths forbidden."""
import json
from pathlib import Path
import socket
import sys
import threading
import time
from types import SimpleNamespace
import RUN

def main(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time();checks=[];calls=[]
    def check(name,ok):checks.append({'name':name,'passed':bool(ok)});assert ok,name
    def reject(name,fn):
        try:fn()
        except ValueError:check(name,True)
        else:check(name,False)
    def forbidden(*a,**k):calls.append(True);raise AssertionError('No native or solver in policy controls')
    cfg=json.loads((RUN.HERE/'CONFIG_TEMPLATE.json').read_text());old=RUN.module(cfg['original_adapter'],'two_point_policy_parent')
    class Fake:
        _parse=forbidden
        _assemble=forbidden
    base=SimpleNamespace(Evaluator=Fake);parent=old.evaluator_class(base);child=RUN.evaluator_class(old,base)
    for name in ('_cell','_invoke','_import_saved','required_programs','_parse','_assemble'):check('inherited_'+name,getattr(parent,name).__code__ is getattr(child,name).__code__)
    obj=child.__new__(child);obj.mode='run';obj.stop_event=threading.Event();obj.batch={'hostname':socket.gethostname()};obj._gated=True;obj.ledger={'native_program_invocations':112,'vector_keys':[]}
    reject('cap112_before_native',obj._gate);obj.ledger['native_program_invocations']=111;obj._gate();check('111_boundary_allowed',True)
    obj.mode='prepare';reject('preparation_no_native',obj._gate);obj.mode='run';obj.stop_event.set();reject('stop_prevents_next_invocation',obj._gate);obj.stop_event.clear()
    obj.allowed={};reject('unknown_vector_rejected',lambda:obj.evaluate([0.]*12));reject('cannot_add_third_vector',lambda:obj.authorize([0.]*12,'trial',1))
    v=[0.]*12;key=RUN.digest(v);obj.allowed={key:{'parameters':v}};obj.ledger['vector_keys']=[key];reject('no_duplicate_point',lambda:obj.evaluate(v));obj.ledger['vector_keys']=['a','b'];reject('two_point_limit',lambda:obj.evaluate(v))
    obj.ledger.update(vector_keys=[],native_program_invocations=57);obj.expected_programs={key:56};obj.required_programs=lambda v:56;reject('insufficient_whole_point_budget',lambda:obj.evaluate(v))
    chi={'La_dimer':1.,'La_trimer':3.,'Dy_dimer':4.,'Dy_trimer':7.};C={'aquo9_central_La':[0.,0.,2.],'aquo9_central_Dy':[0.,0.,3.]}
    metrics={'q':.2,'component_errors_kcal_mol':{'P8':2.,'B8':-3.,'P9':4.,'B9':-5.},'all25_original_physical_gates':False}
    answer={'response':{'chi_states_a0_cubed':chi,'C_dipoles_Debye':C}};pred={'component_errors_kcal_mol':{k:0. for k in metrics['component_errors_kcal_mol']},'actual_component_q':.1}
    targets={'chi_a0_cubed':{k:0. for k in chi},'C_Debye':{'aquo9_central_La':[0.,0.,3.],'aquo9_central_Dy':[0.,0.,2.]}}
    summary=RUN.summarize(answer,pred,SimpleNamespace(assess=lambda *a:metrics),None,targets)
    check('physical_failure_preserved_without_model_admission',summary['metrics']['all25_original_physical_gates'] is False and summary['model_admitted'] is False and summary['selected_candidate'] is None)
    check('signed_electrical_contrast',summary['W9_C_Dy_minus_La_Debye']=={'model':[0.,0.,1.],'reference':[0.,0.,-1.]})
    check('conditional_increment',summary['chi_conditional_increments_a0_cubed']['Dy_minus_La']=={'model':1.,'reference':0.,'error':1.})
    check('zero_native_solver_calls',not calls)
    result={'status':'passed_two_direction_policy_controls','implementation':RUN.pin(__file__),'runner':RUN.pin(RUN.__file__),'checks':checks,'native_calls':0,'solver_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'result':RUN.pin(out/'RESULT.json'),'checks':len(checks)}))
if __name__=='__main__':main(sys.argv[1])
