"""External solver-specific correction; original molecular artifacts stay immutable."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
TOLG = r'^\s*Orbital Gradient\s+TolG\s+\.{3,}\s*('+NUMBER+r')'
TRAH = r'^\s*Converg\.\s+threshold\s+\(grad\. norm\)\s+\.{3,}\s*('+NUMBER+r')'
SUMMARY = r'^\s*Last Orbital Gradient\s+\.{3}\s*'+NUMBER+r'\s+Tolerance\s*:\s*('+NUMBER+r')'

def pin(path):
    path=Path(path).resolve(); h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return {'path':str(path),'sha256':h.hexdigest()}

def check(record):
    if pin(record['path'])['sha256']!=record['sha256']: raise ValueError('Hash mismatch: '+record['path'])
    return Path(record['path'])

def write(path,data):
    with Path(path).open('x') as f: json.dump(data,f,indent=2); f.write('\n')

def original(manifest):
    record=manifest['implementation']['REFERENCE.py']['snapshot']; path=check(record)
    sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location('_pede_original_reference',path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

def actual_trah_threshold(text):
    evidence={key:re.findall(pattern,text,re.M) for key,pattern in [('active_TolG',TOLG),('TRAH_gradient_norm',TRAH),('optional_summary_tolerance',SUMMARY)]}
    if not evidence['active_TolG'] or not evidence['TRAH_gradient_norm']:
        raise ValueError('Missing printed active TolG or TRAH gradient-norm threshold')
    values=[float(x.replace('D','E').replace('d','e')) for row in evidence.values() for x in row]
    if any(not math.isfinite(v) or v<=0 for v in values) or len(set(values))!=1:
        raise ValueError('Conflicting/nonphysical printed gradient thresholds')
    return {'threshold':values[0], 'evidence':evidence,
            'source_lines':[line for line in text.splitlines() if any(re.match(p,line) for p in (TOLG,TRAH,SUMMARY))]}

def gradient_criterion(text,residual):
    evidence=actual_trah_threshold(text)
    evidence['observed_residual']=residual
    evidence['passed']=math.isfinite(residual) and 0<=residual<=evidence['threshold']
    return evidence

def collect(module,manifest,manifest_path,preflight):
    result=module.collect(manifest,manifest_path,preflight)
    for row in result['cells']:
        if row.get('final_solver')!='TRAH' or 'checks' not in row: continue
        row['previous_TRAH_check']=row['checks']['final_TRAH_gradient']
        try:
            output=check(row['output']); evidence=gradient_criterion(output.read_text(errors='replace'),row['final_trah_macro']['orbital_residual'])
            row['final_TRAH_criterion']=evidence
            row['checks']['final_TRAH_gradient']=evidence['passed']
        except Exception as exc:
            row['checks']['final_TRAH_gradient']=False
            row['final_TRAH_criterion_error']=f'{type(exc).__name__}: {exc}'
        row['status']='complete_diagnostic' if all(row['checks'].values()) else 'failed_collection'
    result['schema']='pede_fixed_origin_dft_collection.v2'
    result['complete_cells']=sum(r['status']=='complete_diagnostic' for r in result['cells'])
    result['collector_correction']='Read actual active TRAH/TolG thresholds; optional summary may omit Orbital Gradient; no changed criterion'
    return result
