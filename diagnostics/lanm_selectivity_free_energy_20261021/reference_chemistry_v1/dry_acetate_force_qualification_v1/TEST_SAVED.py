"""Narrow saved/synthetic delta controls; no subprocess or native invocation."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('dry_test',HERE/'REFERENCE.py');a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a)
OUT=a.WS/'reference_chemistry_v1/dry_acetate_force_qualification_v1/saved_controls_v1'
PANEL=a.WS/'reference_chemistry_v1/dry_acetate_force_qualification_v1/preparation_v1/PANEL.json'

def main():
    start=time.monotonic();cpu=time.process_time();OUT.mkdir(parents=True,exist_ok=False);p=a.load(PANEL);p['_path']=str(PANEL);checks=[]
    def check(ok,name):
        if not ok:raise AssertionError(name)
        checks.append(name)
    def rejects(fn,name):
        try:fn()
        except (ValueError,KeyError,TypeError):checks.append(name);return
        raise AssertionError('Accepted '+name)
    # Clearly synthetic analytic linear/quadratic energies; no new quantum data.
    q=copy.deepcopy(p);slopes={('La','radial'):5.,('Dy','radial'):8.,('La','rotation'):3.,('Dy','rotation'):-2.};rows={}
    for m in a.METALS:
        for d in a.DIRECTIONS:q['central'][m]['analytic'][d]=slopes[m,d]
    for c in p['cells']:
        s=slopes[c['metal'],c['direction']];x=c['q_A'];e=(s*x+4*x*x)/a.b.HA
        rows[c['cell_id']]={'energy_hartree':e,'electronic_energy_hartree':.8*e,'D4_energy_hartree':.2*e}
    scored=a.score(q,rows);check(scored['declared_checks']==18 and scored['available_checks']==18 and scored['passed_checks']==18,'All18 synthetic controls retained/pass')
    for d,expected in [('radial',3.),('rotation',-5.)]:
        for step in ('0.01','0.005'):check(abs(scored['directions'][d]['paired']['derivatives_kcal_mol_A'][step]['PBE0_D4']-expected)<1e-12,'Paired sign '+d+' '+step)
    shifted=copy.deepcopy(rows);shifted[a.cid('La','radial','plus_005')]['energy_hartree']+=1e-5
    ss=a.score(q,shifted);check(not ss['primary_checks']['radial/La/analytic_0.005']['passed'] and not ss['primary_checks']['radial/paired/analytic_0.005']['passed'] and not ss['all_numerical_checks_passed'],'Single energy defect cannot hide in paired/full gate')
    missing=copy.deepcopy(rows);del missing[a.cid('Dy','rotation','minus_010')];ss=a.score(q,missing)
    check(ss['declared_checks']==18 and ss['directions']['rotation']['Dy'] is None and ss['directions']['rotation']['paired'] is None and not ss['all_numerical_checks_passed'],'Missing endpoints preserve unavailable comparisons')
    for bad in (float('nan'),float('inf')):
        corrupt=copy.deepcopy(rows);corrupt[a.ALL[0]]['energy_hartree']=bad;rejects(lambda:a.score(q,corrupt),'Nonfinite score rejected '+str(bad))
    # Transform/projection duality on actual coordinates with synthetic forces.
    g=p['central']['La']['original_cell']['geometry'];forces=[[math.sin(i+1),math.cos(i+2),.1*i] for i in range(8)];analytic=a.projection(forces,g)
    for d in a.DIRECTIONS:
        h=1e-5; energies=[]
        for t in (-h,h):energies.append(-sum(forces[i][k]*v for i,x in enumerate(a.transform(g,d,t)) for k,v in enumerate(x)))
        derivative=(energies[1]-energies[0])/(2*h);check(abs(derivative-analytic[d])<1e-8,'Synthetic force/tangent energy derivative '+d)
    text=Path(p['central']['La']['output']['path']).read_text();xc=a.xc_metadata(text)
    check(xc['PostSCFGGA']=='off' and any('Rotationally invariant' in s for s in xc['settings']) and len(xc['counts'])==4,'Actual center XC/grid convention retained')
    altered=text.replace('Gradients option       PostSCFGGA      .... off','Gradients option       PostSCFGGA      .... on')
    check(a.xc_metadata(altered)['PostSCFGGA']=='on','Changed PostSCFGGA reported without false equality')
    rejects(lambda:a.xc_metadata(text.replace('PostSCFGGA','MissingFlag')),'Missing PostSCFGGA rejected')
    # Synthetic source/stage receipts exercise release guards only; never run().
    review=OUT/'SYNTHETIC_SOURCE_REVIEW.json';review.write_text(json.dumps({'status':'passed_dry_acetate_force_source_review','implementation':a.b.pin(a.__file__),'panel':a.b.pin(PANEL)}))
    release={'status':'released_dry_acetate_force_scout','stage':'scout','implementation':a.b.pin(a.__file__),'panel':a.b.pin(PANEL),'ranks':4,'maximum_total_native_calls':16,'maximum_stage_calls':4,'cell_ids':list(a.SCOUT),'task_output_directories':{c:str(OUT/'SYNTHETIC_NEVER_RUN'/c) for c in a.SCOUT},'independent_source_review':a.b.pin(review)}
    def gate(r,c):
        f=OUT/'SYNTHETIC_RELEASE.json';f.write_text(json.dumps(r));return a.release_gate(p,f,c)
    check(gate(release,a.SCOUT[0])['stage']=='scout','Exact four-cell scout gate')
    rejects(lambda:gate(release,a.REMAINING[0]),'Scout cannot release remaining12')
    for key,value in [('maximum_stage_calls',5),('maximum_total_native_calls',17),('ranks',8),('implementation',{'path':a.__file__,'sha256':'0'*64})]:
        bad=copy.deepcopy(release);bad[key]=value;rejects(lambda:gate(bad,a.SCOUT[0]),'Release mismatch '+key)
    rem=copy.deepcopy(release);rem.update(status='released_dry_acetate_force_remaining',stage='remaining',maximum_stage_calls=12,cell_ids=list(a.REMAINING),task_output_directories={c:str(OUT/'SYNTHETIC_NEVER_RUN'/c) for c in a.REMAINING})
    rejects(lambda:gate(rem,a.REMAINING[0]),'Remaining release needs actual scout')
    sr={'status':'dry_acetate_force_scout_complete','stage':'scout','adapter':a.b.pin(a.__file__),'panel':a.b.pin(PANEL),'admitted_cells':4,'stage_numerical_checks_passed':True};sf=OUT/'SYNTHETIC_SCOUT.json';sf.write_text(json.dumps(sr));rem['scout_result']=a.b.pin(sf)
    rr=OUT/'SYNTHETIC_ROOT_REVIEW.json';rr.write_text(json.dumps({'status':'passed_dry_acetate_force_scout','result':rem['scout_result'],'implementation':a.b.pin(a.__file__),'panel':a.b.pin(PANEL)}));rem['scout_review']=a.b.pin(rr)
    check(gate(rem,a.REMAINING[0])['stage']=='remaining','Exact twelve-cell continuation gate')
    sr['stage_numerical_checks_passed']=False;sf.write_text(json.dumps(sr));rem['scout_result']=a.b.pin(sf);rejects(lambda:gate(rem,a.REMAINING[0]),'Failed scout prevents continuation')
    # Execute only the saved collector with sixteen absent paths.
    inventory=OUT/'SYNTHETIC_MISSING_INVENTORY.json';inventory.write_text(json.dumps(dict.fromkeys(a.ALL)))
    a.collect(PANEL,inventory,'full',OUT/'missing_collection')
    result=json.loads((OUT/'missing_collection/RESULT.json').read_text())
    check(len(result['slots'])==16 and result['admitted_cells']==0 and result['scoring']['declared_checks']==18 and result['scoring']['available_checks']==0 and not result['stage_numerical_checks_passed'],'Actual missing collector keeps16slots/18checks')
    result={'status':'passed_saved_delta_controls','implementation':a.b.pin(a.__file__),'test_source':a.b.pin(__file__),'panel':a.b.pin(PANEL),'checks':checks,'passed':len(checks),'new_native_calls':0,'synthetic_fixtures_explicit':True,'threads':1,'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}}
    pin=a.b.write(OUT/'RESULT.json',result);print(json.dumps({'result':pin,'passed':len(checks),'cost':result['cost']},indent=2))

if __name__=='__main__':main()
