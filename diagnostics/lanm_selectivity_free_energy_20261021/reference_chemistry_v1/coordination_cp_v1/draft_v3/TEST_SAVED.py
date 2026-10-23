#!/usr/bin/env python3
"""One targeted saved-output/deck/release control set; zero engine calls."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('coordination_cp_test',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
b=m.b

def main(panel,out):
    tick,cpu=time.monotonic(),time.process_time();p=m.load(panel);out.mkdir(parents=True,exist_ok=False);checks=[];sources=[]
    def check(ok,label): b.check(ok,label);checks.append(label)
    def rejects(fn,label):
        try: fn()
        except (ValueError,KeyError,TypeError): checks.append(label);return
        raise ValueError('Failed rejection: '+label)
    for c in p['cells']:
        g=p['geometries'][c['geometry_id']];text=m.deck(p,c,4,3000)
        check(text==b.verify(c['template_input']).read_text(),'exact_deck:'+c['cell_id'])
        check([list(map(float,l.split()[1:])) for l in text.split('* xyz ',1)[1].splitlines()[1:-1]]==g['positions_A'],'coordinates:'+c['cell_id'])
    # Exercise the reused N-center parser on actual 25/28-center outputs.
    for gid,s in p['saved_complexes'].items():
        g=p['geometries'][gid];c=m.make_cell(g,m.ROLES[0]);c.update(role='complex',d4=True,ghost_mask=[False]*len(g['symbols']),expected_electrons=8+10*g['water_count'],expected_effective_nuclear_charges=[11]+[8,1,1]*g['water_count'])
        r=m.parse(b.verify(s['output']).read_text(),c,0)
        check(r['status']=='complete' and r['energy_hartree']==s['final_energy_hartree'],'whole_complex_raw:'+gid)
        check(len(r['final_COSX_correction_hartree'])==1 and len(r['total_energy_after_final_integration_hartree'])==1,'native_final_stage:'+gid)
        sources.append(s['output'])
    # Reuse all four actual DFT water ghost roles; do not calculate them again.
    path=m.dft.WS/'dft_water_cp_v1/prepared_v1/PANEL.json'
    wp=json.loads(path.read_text())
    samples=[]
    for c0 in wp['cells']:
        if c0['d4']:continue
        c=copy.deepcopy(c0);c['method']='pbe0'
        q=m.dft.WS/'dft_water_cp_v1/run_v1/tasks'/c['cell_id']/'execution/RESULT.json';old=json.loads(q.read_text())
        text=b.verify(old['output']).read_text();r=m.parse(text,c,old['returncode'])
        check(r['status']=='complete' and r['energy_hartree']==old['energy_hartree'],'saved_ghost_raw:'+c['cell_id']);sources.extend((b.pin(q),old['output']));samples.append((text,c))
    text,c=samples[0]
    for field,value in [('expected_electrons',c['expected_electrons']+2),('charge',c['charge']+1),('expected_basis_functions',c['expected_basis_functions']+1),('expected_ecp_core',999),('d4',True),('expected_effective_nuclear_charges',[99]*4)]:
        bad=copy.deepcopy(c);bad[field]=value
        check(m.parse(text,bad,0)['status']!='complete','reject_wrong_'+field)
    check(m.parse(text.replace('ORCA TERMINATED NORMALLY','ORCA STOPPED'),c,0)['status']!='complete','reject_missing_termination')
    check(m.parse(text,c,1)['status']!='complete','reject_native_nonzero')
    import re
    nonfinite=re.sub(r'FINAL SINGLE POINT ENERGY\s+'+m.cp.FLOAT,'FINAL SINGLE POINT ENERGY nan',text)
    check(m.parse(nonfinite,c,0)['status']!='complete','reject_nonfinite_energy')
    # Release admission alone: never invoke run or a native program.
    rel=dict(status='approved_coordination_cp_eight_cells',implementation=p['adapter'],panel=b.pin(panel),cell_ids=[c['cell_id'] for c in p['cells']],maximum_native_calls=8,ranks=4,
        preparation=b.pin(Path(panel).parent/'PREFLIGHT.json'),independent_review=b.pin(HERE/'ROOT_DECLARATION.md'),execution_root=str((out/'nonexecuted').resolve()))
    cid=p['cells'][0]['cell_id'];dest=Path(rel['execution_root'])/'tasks'/cid/'execution'
    for key,val in [('status','unreleased'),('maximum_native_calls',9),('ranks',8),('cell_ids',rel['cell_ids'][:-1])]:
        bad=copy.deepcopy(rel);bad[key]=val;f=out/('rejected_'+key+'.json');b.write(f,bad)
        rejects(lambda:m.gate(p,panel,f,cid,dest),'reject_release_'+key)
    good=out/'synthetic_gate_only.json';b.write(good,rel)
    rejects(lambda:m.gate(p,panel,good,cid,out/'wrong'),'reject_wrong_execution_path')
    dest.mkdir(parents=True);rejects(lambda:m.gate(p,panel,good,cid,dest),'reject_existing_attempt')
    result=dict(status='saved_controls_passed',implementation=b.pin(m.__file__),test_source=b.pin(__file__),panel=b.pin(panel),checks=checks,count=len(checks),source_pins=sources,
        local_threads=1,new_native_calls=0,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,command=sys.argv)
    print(json.dumps(b.write(out/'RESULT.json',result)))

if __name__=='__main__':main(sys.argv[1],Path(sys.argv[2]))
