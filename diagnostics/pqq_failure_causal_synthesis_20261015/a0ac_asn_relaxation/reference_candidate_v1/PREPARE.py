#!/usr/bin/env python3
"""Optional metadata-only binding; does not import or execute scientific code."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / 'workspaces/pqq_failure_causal_synthesis_20261015'
SOURCE = 'a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4'
PROTOCOL = 'A0AC_Ca4_selected_old_new_native_r2SCAN3c_vacuum_v1'
PINS = {
    'source_collection': (BASE/'asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json', '2ed05f5a241f4f1eccbb2a680428a4b032a51449df56ad323e138b83cf32fd82'),
    'source_preparation_collection': (BASE/'asn_relaxation_v2/preparation/PREPARATION.json', '6405de39b9d3ffa88aea7f99f895a685ddcd8ec5b1b37529766b72c365d9d57d'),
    'native_parent': (BASE/'mmol_reference_v1/manifest.json', '10c23bbfcdfacb6a5ff3ea5dbaa6b8661f6ef8cc79681289b4d83b1252b58dc8'),
    'mmol_comparison': (BASE/'mmol_comparison_v1/execution/RESULT.json', 'cf0216391de0680b50b746c9799edd5365c23cc6ee1b9855e39844173903ab61'),
    'selection_collection': (BASE/'asn_relaxation_solvent_v1/execution/FINAL_COLLECTION.json', '64668c7ec758da42d34b039855bdc370b39d06815db47a00dfaf69c41a69368d'),
}
SELECTED = {'Ca': ('adaptive_Ca', 'relaxed_flip_Ca'), 'La': ('adaptive_La', 'relaxed_flip_La')}

def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def check(record):
    if pin(record['path']) != record: raise ValueError('Changed input: '+record['path'])
    return Path(record['path'])

def read(record): return json.loads(check(record).read_text())

def write(path, value):
    with path.open('x') as f: json.dump(value, f, indent=2, allow_nan=False); f.write('\n')

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    pins = {k: {'path': str(v[0].resolve()), 'sha256': v[1]} for k, v in PINS.items()}
    collection, prep, old, fixture = (read(pins[k]) for k in ('source_collection','source_preparation_collection','native_parent','mmol_comparison'))
    read(pins['selection_collection'])
    row, = [r for r in collection['rows'] if r['source_case_id'] == SOURCE]
    prepared, = [r for r in prep['rows'] if r['source_case_id'] == SOURCE]
    original = prepared['original_row']
    if collection['preparation'] != pins['source_preparation_collection'] or row['status'] != 'complete' or prepared['status'] != 'prepared':
        raise ValueError('Actual complete target source/preparation required')
    if fixture['manifest'] != pins['native_parent'] or fixture['available_cells'] != 6:
        raise ValueError('Exact completed MMOL fixture required')
    source = {k: original[k] for k in ('case_id','source_case_id','case_key','record_id','selection_id','maps','state_signature','source_preparation','protonated_source')}
    source.update(sample=4, role='target', atom_count=204, cap_count=14, mapping=prepared['mapping'],
                  prepared_admission=prepared['admission'], contact=prepared['contact'],
                  old_searches=prepared['old_searches'], new_proposals=row['proposals'])
    for key in ('state_signature','source_preparation','protonated_source','mapping','prepared_admission'):
        check(source[key])
    for key in ('origin','identities','source'): check(source['contact'][key])
    w = a.output.resolve(); w.mkdir(parents=True, exist_ok=False)
    impl = w/'implementation'; impl.mkdir()
    names = ('run_orca_task_manifest.py','render_orca_runtime_input.py','_base_render_orca_runtime_input.py',
             'result_protocol.py','OLD_REFERENCE.py','COLLECT_V2.py','EXAA_REFERENCE.py','COUNTED_COLLECT.py')
    copies = {name: check(old['implementation'][name]['snapshot']) for name in names}
    copies['MMOL_REFERENCE.py'] = check(old['implementation']['REFERENCE.py']['snapshot'])
    copies.update({name: Path(__file__).with_name(name) for name in ('PREPARE.py','REFERENCE.py','RUN.sbatch','PLAN.md','READY.md')})
    implementation = {}
    for name, src in copies.items():
        dst = impl/name; shutil.copyfile(src, dst)
        implementation[name] = {'source': pin(src), 'snapshot': pin(dst)}
    tasks = []
    for z, candidates in SELECTED.items():
        for phase, candidate in zip(('old','new'), candidates):
            cell = row['matrix'][z][candidate]; mr = read(cell['receipt'])
            if cell['status'] != 'complete' or (cell['atom_count'],cell['charge'],cell['multiplicity']) != (204,-4 if z=='Ca' else -3,1):
                raise ValueError('Wrong endpoint state')
            request = mr['request']; read(request); check(cell['forces'])
            tid = f'sample4_{z}_{candidate}_vacuum'; td = w/'tasks'/tid; td.mkdir(parents=True)
            shutil.copyfile(check(cell['xyz']), td/'core.xyz')
            (td/'endpoint.inp').write_text(f"! r2SCAN-3c NoAutostart DefGrid3 TightSCF\n* xyzfile {cell['charge']} 1 core.xyz\n")
            tasks.append({'task_id':tid,'sample':4,'candidate':candidate,'phase':phase,'metal':z,'medium':'vacuum',
                'atom_count':204,'charge':cell['charge'],'multiplicity':1,'source_xyz':cell['xyz'],
                'input':pin(td/'endpoint.inp'),'xyz':pin(td/'core.xyz'),'output_path':str(td/'endpoint.out'),
                'saved_MACE':cell['receipt'],'saved_MACE_eV':cell['energy_eV'],'saved_MACE_request':request,
                'saved_forces':cell['forces'],'proposal':prepared['old_searches'][z]['result'] if phase=='old' else row['proposals'][z],
                'full_source_admission':cell.get('full_source_admission')})
    for r in old['runtime_files'] + old['collection_fixtures'] + fixture['collections']: check(r)
    m = {'schema':'A0AC_Ca4_selected_native_reference.v1','protocol_id':PROTOCOL,
         'created_utc':dt.datetime.now(dt.timezone.utc).isoformat(), **pins,
         'sources':[source],'tasks':tasks,'implementation':implementation,'orca':old['orca'],
         'python':pin('/usr/bin/python3'),'runtime_files':old['runtime_files'],'expected_version':'6.1.1',
         'expected_basis':'def2-mTZVPP','expected_auxiliary_basis':'def2-mTZVPP/J','La_ECP_core_electrons':46,
         'execution_policy':{'task_runner':pin(impl/'run_orca_task_manifest.py'),
            'runtime_renderer':pin(impl/'render_orca_runtime_input.py'),'workers':1,'mpi_ranks':12,
            'memory':'unchanged allocation-aware renderer and reserve; no override'},
         'parser_fixtures':old['parser_fixtures'],'fixture_manifest':old['fixture_manifest'],
         'collection_fixtures':old['collection_fixtures'],'fixture_join':old['fixture_join'],
         'mmol_collection_fixtures':fixture['collections'],'selected_candidates':SELECTED,
         'automatic_retry':False,'execution_decision_required':True,
         'plan':pin(impl/'PLAN.md'),'interpretive_resolution_kcal_mol':1.0,
         'scientific_scope':'four fixed selected old/new energies; distinct metal-specific paths; no pool minimum, solvent or classifier'}
    write(w/'manifest.json', m)
    print(json.dumps(pin(w/'manifest.json')))

if __name__ == '__main__': main()
