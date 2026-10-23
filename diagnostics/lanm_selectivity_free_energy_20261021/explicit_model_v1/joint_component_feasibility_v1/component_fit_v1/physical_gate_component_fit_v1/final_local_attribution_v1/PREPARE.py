#!/usr/bin/env python3
"""Exact final-probe cache inventory only. Never import an evaluator or solver."""
import ast
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent

def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def read(authority):
    assert pin(authority['path'])==authority
    return json.loads(Path(authority['path']).read_text())

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def pure_functions(authority,names):
    assert pin(authority['path'])==authority
    tree=ast.parse(Path(authority['path']).read_text())
    selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert {n.name for n in selected}==set(names)
    env={'np':np,'math':math,'json':json,'hashlib':hashlib}
    exec(compile(ast.Module(body=selected,type_ignores=[]),authority['path'],'exec'),env)
    return env

def main():
    tick,cpu=time.monotonic(),time.process_time()
    assert len(sys.argv)==2
    out=Path(sys.argv[1]).resolve()
    assert not out.exists()
    inp=json.loads((HERE/'INPUTS.json').read_text());s=inp['sources']
    for p in s.values():assert pin(p['path'])==p
    fit,audit=read(s['fit_result']),read(s['fit_review'])
    assert fit['status']=='bounded_physical_gate_component_fit_complete' and fit['error'] is None
    assert audit['status']=='passed_physical_gate_component_fit_actual_review' and audit['result']==s['fit_result']
    assert audit['selected']==fit['selected']==s['selected'] and audit['selected_native_admission_passed']
    assert fit['selected_key']==inp['center_key'] and fit['selected_metrics']['all25_original_physical_gates']
    center,domain,manifest=read(s['selected']),read(s['domain']),read(s['manifest'])
    old=read(s['reuse_inventory']);prop=read(s['last_proposal'])
    assert old['status']=='prepared_immutable_exact_reuse_index'
    assert len(manifest['rows'])==46 and len(manifest['static'])==48
    core=pure_functions(s['core'],('require','token','checked','probe_points'))
    native=pure_functions(s['evaluator'],('require','checked','effective'))
    assert core['token'](center['parameters'])==inp['center_key']==digest(center['parameters'])
    points=core['probe_points'](center['parameters'],domain)
    assert len(points)==24 and len({digest(v) for v in points})==24
    # Only terminal-index metadata are joined here. Historical cell bodies/ledgers
    # and their derived LAPACK fields are neither reopened nor recomputed.
    registry=dict(old['registry'])
    for key,row in fit['vectors'].items():
        assert row['status']=='complete' and digest(row['parameters'])==key
        registry[key]={'parameters':row['parameters'],'artifact':row['artifact'],'kind':'vector','source':s['fit_result']}
    for key,row in registry.items():assert digest(row['parameters'])==key
    effective={metal:{} for metal in ('La','Dy')}
    # Prefer the already admitted selected body; other exact matches are metadata
    # promises only until a future executor verifies their immutable cell records.
    order=[inp['center_key'],*sorted(k for k in registry if k!=inp['center_key'])]
    for key in order:
        for metal in effective:
            eff=digest(native['effective'](registry[key]['parameters'],metal,domain))
            effective[metal].setdefault(eff,key)
    cells={m:[cid for cid,row in manifest['rows'].items() if row['metal']==m] for m in effective}
    costs={m:sum(1+int(manifest['rows'][cid]['needs_gradient']) for cid in ids) for m,ids in cells.items()}
    assert costs=={'La':28,'Dy':28}
    rows=[];new_effective={m:set() for m in effective};spans=[]
    for i,p in enumerate(domain['parameters']):
        width=p['upper']-p['lower'];minus,plus=points[2*i:2*i+2]
        spans.append({'parameter':p['name'],'minus':minus[i],'center':center['parameters'][i],'plus':plus[i],
            'minus_normalized_distance':(center['parameters'][i]-minus[i])/width,
            'plus_normalized_distance':(plus[i]-center['parameters'][i])/width,
            'actual_normalized_span':(plus[i]-minus[i])/width,
            'clamped_minus':minus[i]==p['lower'],'clamped_plus':plus[i]==p['upper']})
        for sign,v in [(-1,minus),(1,plus)]:
            key=digest(v);matches={};missing=0
            for metal in effective:
                values=native['effective'](v,metal,domain);token=digest(values);match=effective[metal].get(token)
                if match is None:
                    new_effective[metal].add(token);missing+=costs[metal]
                matches[metal]={'effective_values_a_b_c_d_R_epsilon':values,'effective_key':token,
                    'matching_saved_vector':None if match is None else {'key':match,**registry[match]},
                    'cells':cells[metal],'potential_new_programs':costs[metal] if match is None else 0}
            rows.append({'ordinal':len(rows)+1,'parameter':p['name'],'sign':sign,'parameters':v,'key':key,
                'exact_full_vector_match':registry.get(key),'metal_cache_matches':matches,'potential_new_programs':missing,
                'future_native_result':None})
    result={'status':'prepared_final_local_probe_inventory_no_execution',
        'implementation':pin(__file__),'inputs':pin(HERE/'INPUTS.json'),'sources':s,
        'center_key':inp['center_key'],'selected':s['selected'],'center_parameters':center['parameters'],
        'cache_scope':{'historical_registry_vectors':len(old['registry']),'terminal_vectors':len(fit['vectors']),
            'unique_joined_vectors':len(registry),'historical_vector_bodies_read':0,'mutable_ledger_reads':0,
            'admission':'Inherited exact terminal/reuse indexes. Any future imported cell still requires original cache-identity and raw/provenance checks; no derived LAPACK-byte comparison.'},
        'probe_rows':rows,'spans':spans,'full_vector_matches':sum(r['exact_full_vector_match'] is not None for r in rows),
        'new_effective_metal_states':{m:len(v) for m,v in new_effective.items()},
        'potential_new_programs_after_exact_metadata_reuse':sum(len(v)*costs[m] for m,v in new_effective.items()),
        'potential_programs_before_reuse':24*56,'per_metal_programs':costs,
        'last_saved_Jacobian':{'proposal':s['last_proposal'],'origin_key':prop['origin_key'],
             'at_final_selected':prop['origin_key']==inp['center_key'],'classification':'off-center step7-selected-origin evidence only'},
        'future_bindings':inp['future_bindings'],'native_calls':0,'solver_calls':0,'reference_calls':0,
        'array_package_version':np.__version__,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    assert not result['last_saved_Jacobian']['at_final_selected']
    assert result['potential_new_programs_after_exact_metadata_reuse']<=1344
    out.mkdir(parents=True,exist_ok=False)
    (out/'INVENTORY.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'inventory':pin(out/'INVENTORY.json'),'counts':{k:result[k] for k in ('full_vector_matches','new_effective_metal_states','potential_new_programs_after_exact_metadata_reuse')},'spans':spans,'cost':result['cost']},indent=2))

if __name__=='__main__':main()
