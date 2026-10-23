#!/usr/bin/env python3
"""One focused saved metadata/control set. No chemistry or preparation rerun."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('trimer_preparation',HERE/'PREPARE.py')
p=importlib.util.module_from_spec(s);s.loader.exec_module(p)


def main():
    tick=time.monotonic(); cpu=time.process_time(); checks=[]
    cfg=json.loads((HERE/'CONFIG_PREPARE.json').read_text())
    receipt=p.read(p.pin(p.W/'explicit_model_v1/collective_mechanics_options_v1/preparation_v1/RESULT.json'))
    m=p.read(receipt['manifest']); common=p.module(cfg['common_adapter'],'trimer_control_common')
    panel=p.read(cfg['common_panel']); model=p.read(cfg['model_result']); mm=p.read(model['manifest'])
    def check(ok,label):p.require(ok,label);checks.append(label)
    def rejects(fn,label):
        try:fn()
        except (ValueError,KeyError):checks.append(label);return
        raise AssertionError(label)
    expected={8:{'near':[1,5],'far':[1,6]},9:{'near':[1,5],'far':[1,3]}}
    check(len(m['reference_cells'])==len(m['model_rows'])==16 and len(m['cycles'])==8,'complete fixed denominators')
    for gid,g in panel['geometries'].items():
        choice=p.selection(g)
        check({k:v['water_indices'] for k,v in choice.items()}==expected[g['water_count']],'geometry-selected indices '+gid)
        check(choice==m['selections'][gid],'exact selection replay '+gid)
        bad=copy.deepcopy(g);bad['positions_A'][0][0]=float('nan')
        rejects(lambda:p.selection(bad),'nonfinite geometry rejects '+gid)
    qids={c['cell_id'] for c in m['reference_cells']}
    for c in m['reference_cells']:
        text=p.verify(c['template_input']).read_text()
        check(text==common.deck(panel,c,4,3000),'exact qualified deck rendering '+c['cell_id'])
        check(len(c['real_parent_atom_indices'])==(7 if c['role']=='metal_two_water' else 6),'real-site mask '+c['cell_id'])
        check(c['expected_electrons']==(28 if c['charge']==3 else 20),'electron count '+c['cell_id'])
    for mid,row in m['model_rows'].items():
        parent=mm['rows'][row['parent_geometry_id']]
        xyz=[r.split() for r in p.verify(row['input']['xyz']).read_text().splitlines()[1:]]
        check(len(xyz)==row['atoms'] and [int(r[5]) for r in xyz]==row['types'],'typed literal model '+mid)
        check(all(row['input'][k]['sha256']==parent['input'][k]['sha256'] for k in ('key','parameters')),'unchanged model terms '+mid)
        check(row['frame_signatures']==[parent['frame_signatures'][a] for a in row['parent_atom_indices']],'exact local frames '+mid)
    for cid,cycle in m['cycles'].items():
        q=cycle['reference_terms'];t=cycle['model_terms'];gid=cycle['parent_geometry_id']
        check(len(q)==len(t)==7 and sum(q.values())==sum(t.values())==1,'seven-term signs '+cid)
        check(sum(k in qids for k in q)==2 and all(k in qids or k in m['reused_reference'] for k in q),'new/reused reference joins '+cid)
        masks={k:(m['model_rows'][k] if k in m['model_rows'] else mm['rows'][k])['parent_atom_indices'] for k in t}
        check(p.proof(masks,t,mm['rows'][gid]['atoms'])==cycle['permanent_and_pair_cancellation'],'cancellation replay '+cid)
        wrong=dict(t);wrong[next(iter(wrong))]*=-1
        rejects(lambda:p.proof(masks,wrong,mm['rows'][gid]['atoms']),'wrong cycle sign rejects '+cid)
        incomplete=dict(masks);incomplete.pop(next(iter(incomplete)))
        rejects(lambda:p.proof(incomplete,t,mm['rows'][gid]['atoms']),'missing fragment rejects '+cid)
    g=panel['geometries']['aquo8_central_La']; near=p.selection(g)['near']
    rejects(lambda:p.definition(common,g,near,'near','unsupported'),'unsupported fragment role rejects')
    bad=dict(near,water_indices=[5,1]);rejects(lambda:p.definition(common,g,bad,'near','two_water'),'unordered atom selection rejects')
    out=p.W/'explicit_model_v1/collective_mechanics_options_v1/controls_v1';out.mkdir(exist_ok=False)
    r={'status':'passed_focused_saved_metadata_controls','implementation':p.pin(__file__),
       'preparation_source':cfg['implementation'],'manifest':receipt['manifest'],'checks_passed':len(checks),
       'checks':checks,'new_native_calls':0,'preparation_repeated':False,
       'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
    print(json.dumps({'result':p.write(out/'RESULT.json',r),'checks':len(checks)}))


if __name__=='__main__':main()
