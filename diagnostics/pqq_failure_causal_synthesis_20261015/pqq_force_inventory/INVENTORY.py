#!/usr/bin/env python3
"""Read exact receipt/maps metadata and hash files; never load coordinates/forces."""
import hashlib
import json
from pathlib import Path
import re

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'scripts/context_factorial_gfn2.py').is_file())
HERE=Path(__file__).resolve().parent

def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}
def check(p):
    if pin(p['path'])!=p:raise ValueError('Changed artifact: '+p['path'])
    return Path(p['path'])
def read(p):return json.loads((check(p) if isinstance(p,dict) else Path(p)).read_text())

summary_path=ROOT/'diagnostics/exaa_donor_response_20261014/FINAL_RESULT.json'
summary=read(summary_path);manifest=read(summary['manifest'])
source_path=ROOT/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json';sources=read(source_path)
schema_path=ROOT/'workspaces/pede_solvent_response_20261013/cofactor_review_v2/manifest.json';schema=read(schema_path)
dftevidence={r['task_id']:r for r in summary['endpoint_diagnostics']};rows=[];contexts=[]
for source in manifest['sources']:
    archived,=[r for r in sources['rows'] if r['case_id']==source['case_id']]
    cp=read(source['source_preparation']);mapping=read(source['maps']['Ca']);core=read(archived['original_core']['parent'])
    for p in source['maps'].values():check(p)
    check(source['state_signature'])
    fragment,=[f for f in core['qm_fragments'] if f['role']=='pqq_cofactor']
    atoms=[r for r in cp['mapping']['source_to_qm'] if r.get('fragment')==fragment['id']]
    named={r['name']:r['qm_index'] for r in atoms}
    if len(named)!=27 or fragment['microstate_id']!='pqq_ox_3minus_v1':raise ValueError('PQQ inventory/state differs')
    moves=[m['id'] for m in mapping['context']['modes'] if set(m['moving_indices'])&set(named.values())]
    contexts.append({'sample':source['sample'],'case_id':source['case_id'],'maps':source['maps'],
      'state_signature':source['state_signature'],'source_preparation':source['source_preparation'],
      'original_core_preparation':archived['original_core']['parent'],'microstate_id':fragment['microstate_id'],
      'formal_charge':fragment['formal_charge'],'PQQ_fragment':fragment['id'],'PQQ_indices_zero_based':named,
      'PQQ_atom_identity':[{k:r[k] for k in ('name','element','qm_index','kind','source') if k in r} for r in atoms],
      'PQQ_hydrogen_parents':{r['name']:r['parent_atom'] for r in fragment['atom_records'] if r['element']=='H'},
      'existing_modes_moving_PQQ':moves})
for t in manifest['tasks']:
    d=dftevidence[t['task_id']];out=check(d['output']);check(d['execution']);request=read(t['saved_MACE'])
    for p in (t['input'],t['xyz'],t['source_xyz'],request['forces']):check(p)
    q=read(request['request']);check(q['xyz'])
    if q['xyz']['sha256']!=t['source_xyz']['sha256'] or t['xyz']['sha256']!=t['source_xyz']['sha256']:
        raise ValueError('DFT/MACE exact coordinate identity differs')
    if request['status']!='complete' or q['charge']!=t['charge'] or q['multiplicity']!=t['multiplicity']:
        raise ValueError('MACE endpoint state/receipt differs')
    directory=out.parent
    gradfiles=[pin(p) for p in sorted(directory.iterdir()) if p.is_file() and (p.suffix.lower() in ('.engrad','.grad') or 'gradient' in p.name.lower())]
    header=bool(re.search(r'^\s*(?:CARTESIAN GRADIENT|ORCA SCF GRADIENT)\s*$',out.read_text(),re.M))
    input_text=check(t['input']).read_text()
    gbw=directory/'endpoint.runtime.gbw'
    rows.append({'task_id':t['task_id'],'sample':t['sample'],'candidate':t['candidate'],'metal':t['metal'],
      'charge':t['charge'],'multiplicity':t['multiplicity'],'xyz':t['xyz'],'source_xyz':t['source_xyz'],
      'DFT':{'status':'force_unavailable' if not gradfiles and not header else 'requires_review',
        'method_input':t['input'],'first_input_line':input_text.splitlines()[0],
        'output':d['output'],'execution':d['execution'],'gradient_artifacts':gradfiles,'Cartesian_gradient_header_present':header,
        'GBW':{'path':str(gbw),'exists':gbw.is_file(),'bytes':gbw.stat().st_size if gbw.exists() else None},
        'reason':'Completed single-point energy calculation; no analytic nuclear gradient requested or retained',
        'density_response_qualification':d['density_response_qualification']},
      'MACE':{'status':'saved_force_artifact_present_not_numerically_loaded','receipt':t['saved_MACE'],
        'request':request['request'],'forces':request['forces'],'model':request['model'],'state':request['state'],
        'force_definition':request['force_definition'],'force_array_shape_check':'deferred to Slurm',
        'exact_DFT_coordinate_bytes_matched':True}})
result={'protocol':'ExaA_PQQ_saved_force_inventory_v1','implementation':pin(__file__),
  'DFT_final_result':pin(summary_path),'DFT_manifest':summary['manifest'],'source_inventory':pin(source_path),
  'cofactor_schema_authority':pin(schema_path),'name_map_generic_to_CCD':schema['name_map'],
  'cofactor_graph':schema['schema_bonds'],'quinone_bonds_generic':[['C5','O1'],['C6','O2']],
  'quinone_bonds_CCD':[['C5','O5'],['C4','O4']],'contexts':contexts,'rows':rows,
  'declared_endpoints':8,'DFT_forces_available':sum(r['DFT']['status']!='force_unavailable' for r in rows),
  'MACE_force_files_present':len(rows),'new_molecular_calls':0,'coordinate_or_force_arrays_loaded':False,
  'conclusion':'Saved MACE projection feasible; DFT PQQ force comparison unavailable without a separately declared gradient calculation.'}
with (HERE/'INVENTORY.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps(pin(HERE/'INVENTORY.json')))
