"""Manufactured failure-preservation checks only; no real outcome reads."""
import copy
import json
import math
from pathlib import Path
import sys
import time
import REVIEW as R

def main(out):
    tick,cpu=time.perf_counter(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False);checks=[]
    def check(name,ok):checks.append({'name':name,'passed':bool(ok)});assert ok,name
    def reject(name,fn):
        try:fn()
        except ValueError:check(name,True)
        else:check(name,False)
    geos=('water_minus','water_central','water_plus','aquo8_central','aquo9_central')
    model={};data={};old={'per_fixture':{k:{} for k in ('La','Dy','paired')},'energy_controls':{'balanced_exchange':{},'individual_radial_work':{},'paired_radial_work':{}}}
    for geo in geos:
        for metal in ('La','Dy'):
            key=geo+'_'+metal;forces=[[0.,0.,0.] for _ in range(4)]
            if key=='water_minus_Dy':forces[0][0]=6.
            model[key]={'energy':2. if key=='water_plus_Dy' else 0.,'forces':forces}
            data[key]={'geometry':{'positions_A':[[0.,0.,0.],[2.,0.,0.],[2.,1.,0.],[2.,0.,1.]]},'bodies':[{'category':'ion','indices':[0]},{'category':'water','indices':[1,2,3]}],'reference':{'energy':0.,'forces':[[0.,0.,0.] for _ in range(4)]}}
        for kind in ('La','Dy','paired'):
            fail=geo=='water_minus' and kind in ('Dy','paired')
            old['per_fixture'][kind][geo]={'raw_body_components_kcal_mol_A':[6. if fail else 0.]+[0.]*8,'rigid_RMS_kcal_mol_A':math.sqrt(6.) if fail else 0.,'maximum_body_component_error_kcal_mol_A':6. if fail else 0.,'passed':not fail}
    for geo in ('water_minus','water_central','water_plus','aquo9_central'):old['energy_controls']['balanced_exchange'][geo]={'error_kcal_mol':2. if geo=='water_plus' else 0.,'passed':geo!='water_plus'}
    for metal in ('La','Dy'):
        for sign in ('minus','plus'):old['energy_controls']['individual_radial_work'][metal+'_'+sign]={'error_kcal_mol':2. if (metal,sign)==('Dy','plus') else 0.,'passed':(metal,sign)!=('Dy','plus')}
    for sign in ('minus','plus'):old['energy_controls']['paired_radial_work'][sign]={'error_kcal_mol':2. if sign=='plus' else 0.,'passed':sign!='plus'}
    v={'model':model,'score':{'old_EF':old,'original_EF_checks_passed':False}}
    gates=R.independent25(v,data)
    check('complete25_with5_failures',len(gates)==25 and sum(not x['passed'] for x in gates.values())==5)
    check('body_failure_preserved',not gates['force/Dy/water_minus']['passed'])
    check('paired_force_failure_preserved',not gates['force/paired/water_minus']['passed'])
    check('energy_failure_preserved',not gates['energy/paired_radial_work/plus']['passed'])
    corrupt=copy.deepcopy(v);corrupt['score']['original_EF_checks_passed']=True;reject('false_global_pass_rejected',lambda:R.independent25(corrupt,data))
    corrupt=copy.deepcopy(v);corrupt['score']['old_EF']['per_fixture']['Dy']['water_minus']['passed']=True;reject('false_body_pass_rejected',lambda:R.independent25(corrupt,data))
    corrupt=copy.deepcopy(v);corrupt['model']['water_plus_Dy']['energy']=2.1;reject('wrong_saved_cycle_rejected',lambda:R.independent25(corrupt,data))
    reject('flag_change_not_arithmetic_tolerance',lambda:R.compare({'passed':False},{'passed':True},'manufactured'))
    reject('nonfinite_rejected',lambda:R.near(float('nan'),1.,'manufactured'))
    R.compare({'value':1.},{'value':1.+1e-12},'inherited arithmetic');check('inherited_arithmetic_only',True)
    result={'status':'passed_two_direction_saved_reviewer_manufactured_controls','implementation':R.pin(__file__),'reviewer':R.pin(R.__file__),'checks':checks,'real_outcomes_read':0,'native_calls':0,'solver_calls':0,'cost':{'wall_seconds':time.perf_counter()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'result':R.pin(out/'RESULT.json'),'checks':len(checks)}))
if __name__=='__main__':main(sys.argv[1])
