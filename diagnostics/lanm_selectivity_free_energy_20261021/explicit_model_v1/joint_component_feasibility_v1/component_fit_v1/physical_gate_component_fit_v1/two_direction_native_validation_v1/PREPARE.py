#!/usr/bin/env python3
"""Two fixed points: admitted cache metadata inventory, zero molecular calls."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read(p):require(pin(p['path'])==p,'Immutable input');return json.loads(Path(p['path']).read_text())
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def main(out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out);require(not out.exists(),'Fresh inventory')
    inp=json.loads((HERE/'INPUTS.json').read_text());s=inp['sources']
    for p in s.values():require(pin(p['path'])==p,'Exact input '+p['path'])
    plan=read(s['plan']);panel=read(plan['panel_result']);affine=read(plan['affine_result']);root=read(s['panel_root_review']);audit=read(root['saved_review'])
    require(root['status']=='admitted_finite_final_local_attribution' and root['result']==plan['panel_result'] and root['affine_result']==plan['affine_result'] and root['next_plan']==s['plan'],'Exact root panel admission/plan')
    require(audit['status']=='passed_final_local_attribution_saved_review' and audit['result']==plan['panel_result'] and audit['affine_result']==plan['affine_result'],'Actual panel saved admission')
    fit=read(s['fit_result']);old=read(s['reuse_inventory']);domain=read(s['domain']);manifest=read(s['manifest']);center=read(s['selected'])
    require(root['center']==s['selected'] and digest(center['parameters'])==plan['center_key'],'Exact final center')
    require(panel['status']=='completed_final_local_attribution_panel' and panel['error'] is None and len(panel['vectors'])==24,'Terminal24 derivative panel')
    hp=s['inventory_helper'];sp=importlib.util.spec_from_file_location('admitted_metadata_helper',hp['path']);helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)
    native=helper.pure_functions(s['evaluator'],('require','checked','effective'))
    registry=dict(old['registry'])
    for source,value in ((s['fit_result'],fit),(plan['panel_result'],panel)):
        for key,row in value['vectors'].items():
            require(row['status']=='complete' and digest(row['parameters'])==key,'Admitted complete metadata')
            registry[key]={'parameters':row['parameters'],'artifact':row['artifact'],'kind':'vector','source':source}
    order=[plan['center_key'],*sorted(k for k in registry if k!=plan['center_key'])];index={m:{} for m in ('La','Dy')}
    for key in order:
        require(digest(registry[key]['parameters'])==key,'Exact metadata vector key')
        for metal in index:index[metal].setdefault(digest(native['effective'](registry[key]['parameters'],metal,domain)),key)
    cells={m:[cid for cid,r in manifest['rows'].items() if r['metal']==m] for m in index}
    costs={m:sum(1+int(manifest['rows'][cid]['needs_gradient']) for cid in ids) for m,ids in cells.items()};require(costs=={'La':28,'Dy':28},'Unchanged46/56 panel')
    rows=[];missing={m:set() for m in index}
    for row in plan['rows']:
        pred=next(p for p in affine['cases'][row['case']]['predictions'] if p['fraction']==row['fraction'])
        require(pred['parameters']==row['parameters'] and pred['key']==row['key']==digest(row['parameters']),'Exact fixed half witness')
        require(pred['literal_original_box_pass'] and pred['retained_exact_sign_feasible'],'Diagnostic witness retained constraints/domain')
        matches={}
        for metal in index:
            effective=native['effective'](row['parameters'],metal,domain);key=digest(effective);match=index[metal].get(key)
            if match is None:missing[metal].add(key)
            matches[metal]={'effective_values':effective,'effective_key':key,'matching_saved_key':match,'cells':cells[metal],'potential_new_programs':costs[metal] if match is None else 0}
        rows.append(dict(row,metal_cache_matches=matches,exact_full_vector_match=registry.get(row['key']),potential_new_programs=sum(x['potential_new_programs'] for x in matches.values())))
    require([r['id'] for r in rows]==['B_half','D_half'] and all(r['potential_new_programs']<=r['expected_new_programs'] for r in rows),'Only declared two points and budgets')
    total=sum(len(v)*costs[m] for m,v in missing.items());require(total<=112,'112 maximum')
    out.mkdir(parents=True,exist_ok=False)
    result={'status':'prepared_two_direction_metadata_no_execution','implementation':pin(__file__),'inputs':pin(HERE/'INPUTS.json'),'sources':s,'plan':s['plan'],
        'center_key':plan['center_key'],'selected':s['selected'],'rows':rows,'registry':registry,'saved_metal_index':index,
        'unique_metadata_vectors':len(registry),'historical_vector_bodies_read':0,'mutable_ledger_reads':0,'potential_new_programs':total,
        'native_calls':0,'solver_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (out/'INVENTORY.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'inventory':pin(out/'INVENTORY.json'),'potential_new_programs':total,'metadata_vectors':len(registry)}))
if __name__=='__main__':main(sys.argv[1])
