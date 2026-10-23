#!/usr/bin/env python3
"""Targeted extraction, existing parser and release controls; no native calls."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('pair_saved_test',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
b=m.b

def main(path,out):
    tick,cpu=time.monotonic(),time.process_time();p=m.load(path);out.mkdir(parents=True,exist_ok=False);checks=[];sources=[]
    def check(ok,label):b.check(ok,label);checks.append(label)
    def rejects(fn,label):
        try:fn()
        except (ValueError,KeyError,TypeError):checks.append(label);return
        raise ValueError('Failed rejection: '+label)
    per_parent={gid:[] for gid in m.GIDS}
    for c in p['cells']:
        par=p['parent_geometries'][c['parent_geometry_id']];g=c['geometry'];sy,x=b.xyz_read(g['xyz']['path']);text=m.deck(p,c,4,3000)
        check(sy==[par['symbols'][i] for i in c['parent_atom_indices']] and x==[par['positions_A'][i] for i in c['parent_atom_indices']] and text==b.verify(c['template_input']).read_text(),'literal_atom_subset_and_deck:'+c['cell_id'])
        if c['role']=='dimer':per_parent[c['parent_geometry_id']].append(c['water_index'])
    for gid,ix in per_parent.items():check(ix==list(range(1,p['parent_geometries'][gid]['water_count']+1)),'complete_literal_sum:'+gid)
    # The qualified parser must still admit both actual real-atom four-center outputs.
    old=b.load_panel(p['parent_reference_panel']['path'])
    for metal in ('La','Dy'):
        g=next(g for g in old['geometries'] if g['geometry_id']=='water_central_'+metal)
        c=copy.deepcopy(next(c for c in p['cells'] if c['role']=='dimer' and c['metal']==metal));c.update(symbols=g['symbols'],positions_A=g['positions_A'])
        rp=m.helper.cp.ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/pilot_v1/tasks'/('ecp_pbe0__water_central_'+metal)/'execution/RESULT.json'
        r=json.loads(rp.read_text());actual=m.parse(b.verify(r['output']).read_text(),c,r['returncode'])
        check(actual['status']=='complete' and actual['energy_hartree']==r['energy_hartree'],'prior_real_dimer_native_parser:'+metal);sources.extend((b.pin(rp),r['output']))
    # Isolate the new zero-dispersion guard from native syntax (which awaits execution).
    ion=next(c for c in p['cells'] if c['role']=='bare_ion');original=m.helper.parse
    try:
        m.helper.parse=lambda *args:dict(status='complete',checks={'upstream':True},D4_energy_hartree=0.)
        check(m.parse('',ion,0)['status']=='complete','isolated_D4_zero_guard_accepts')
        m.helper.parse=lambda *args:dict(status='complete',checks={'upstream':True},D4_energy_hartree=.001)
        check(m.parse('',ion,0)['status']=='failed_admission','isolated_D4_nonzero_guard_rejects')
    finally:m.helper.parse=original
    mutations=[]
    bad=copy.deepcopy(p);bad['cells'][0]['geometry']['positions_A'][1][0]+=.001;mutations.append(('changed_subset_position',bad))
    bad=copy.deepcopy(p);bad['cells'][0]['ghost_mask'][1]=True;mutations.append(('introduced_ghost',bad))
    bad=copy.deepcopy(p);bad['cells']=bad['cells'][:-1];mutations.append(('missing_own_ion',bad))
    bad=copy.deepcopy(p);bad['cells'][-1]['expected_ecp_core']=0;mutations.append(('changed_ion_core',bad))
    for label,bad in mutations:
        f=out/(label+'.json');b.write(f,bad);rejects(lambda:m.load(f),'reject_'+label)
    rel=dict(status='approved_matched_pair_36_cells',implementation=p['adapter'],panel=b.pin(path),cell_ids=[c['cell_id'] for c in p['cells']],maximum_native_calls=36,ranks=4,
        preparation=b.pin(Path(path).parent/'PREFLIGHT.json'),independent_review=b.pin(HERE/'ROOT_DECLARATION.md'),execution_root=str((out/'nonexecuted').resolve()))
    cid=p['cells'][0]['cell_id'];dest=Path(rel['execution_root'])/'tasks'/cid/'execution'
    for key,val in [('status','unreleased'),('maximum_native_calls',37),('cell_ids',rel['cell_ids'][:-1]),('ranks',8)]:
        bad=copy.deepcopy(rel);bad[key]=val;f=out/('reject_release_'+key+'.json');b.write(f,bad);rejects(lambda:m.gate(p,path,f,cid,dest),'reject_release_'+key)
    f=out/'synthetic_gate_only.json';b.write(f,rel);rejects(lambda:m.gate(p,path,f,cid,out/'wrong'),'reject_wrong_attempt_path')
    dest.mkdir(parents=True);rejects(lambda:m.gate(p,path,f,cid,dest),'reject_existing_attempt')
    r=dict(status='saved_controls_passed',implementation=b.pin(m.__file__),test_source=b.pin(__file__),panel=b.pin(path),checks=checks,count=len(checks),source_pins=sources,new_native_calls=0,local_threads=1,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,command=sys.argv,native_own_ion_admission_pending=True)
    print(json.dumps(b.write(out/'RESULT.json',r)))

if __name__=='__main__':main(sys.argv[1],Path(sys.argv[2]))
