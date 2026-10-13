#!/usr/bin/env python3
"""Slurm-only finite reference adapter; unchanged native task runner underneath."""
import argparse
from collections import Counter
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import traceback

NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'

def pin(path):
    path = Path(path).resolve(); h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return {'path': str(path), 'sha256': h.hexdigest()}

def check(rec):
    if pin(rec['path'])['sha256'] != rec['sha256']: raise ValueError('SHA mismatch: '+rec['path'])
    return Path(rec['path'])

def write(path, data):
    with Path(path).open('x') as f: json.dump(data, f, indent=2); f.write('\n')

def readxyz(path):
    lines = Path(path).read_text().splitlines(); n = int(lines[0]); rows = [r.split() for r in lines[2:] if r.strip()]
    if len(rows) != n or n != 192: raise ValueError('Expected exact192 atoms')
    for row in rows:
        if len(row) != 4 or not all(math.isfinite(float(v)) for v in row[1:]): raise ValueError('Invalid XYZ')
    return rows

def preflight(m, mp):
    import render_orca_runtime_input as runtime
    from run_orca_task_manifest import load_manifest_tasks
    check(m['source_collection']); check(m['orca'])
    for record in m['runtime_files']: check(record)
    for record in m['implementation'].values(): check(record['snapshot'])
    load_manifest_tasks(mp)
    policy = runtime.allocation_policy(workers=1, ranks=12)
    if policy['workers'] != 1 or policy['ranks_per_worker'] != 12: raise ValueError('Need actual12 MPI allocation')
    expected = {f'sample{s}_{z}_{v}' for s in (1,3) for z in ('Ca','La') for v in ('water','vacuum')}
    if {t['task_id'] for t in m['tasks']} != expected or len(m['tasks']) != 8: raise ValueError('Wrong frozen matrix')
    atoms = {}; states = []
    for source in m['sources']:
        for key in ('prepared_admission','source_preparation','state_signature'): check(source[key])
        a = json.loads(Path(source['prepared_admission']['path']).read_text())
        c = json.loads(Path(source['source_preparation']['path']).read_text())
        state = json.loads(Path(source['state_signature']['path']).read_text())
        if a['atom_count'] != 192 or a['cap_count'] != 14 or state['water_inventory'] != []: raise ValueError('Changed physical inventory')
        mapping = c['mapping']['source_to_qm']; indices = [r['qm_index'] for r in mapping]
        if sorted(indices) != list(range(1,192)): raise ValueError('Mapping does not cover all191 nonmetals')
        pair = {}
        for z in ('Ca','La'): pair[z] = readxyz(check(source['endpoints'][z]))
        if pair['Ca'][1:] != pair['La'][1:] or pair['Ca'][0][1:] != pair['La'][0][1:]: raise ValueError('Paired coordinates differ')
        if pair['Ca'][0][0] != 'Ca' or pair['La'][0][0] != 'La': raise ValueError('Wrong native metal')
        for z, rows in pair.items(): atoms[(source['sample'],z)] = rows
        states.append(source['state_signature']['sha256'])
    if len(set(states)) != 1: raise ValueError('Physical states differ')
    atomic_numbers = {'H':1,'C':6,'N':7,'O':8,'S':16,'Ca':20,'La':57}
    records = []
    for task in m['tasks']:
        inp = check(task['input']); xyz = check(task['xyz']); check(task['source_xyz'])
        if task['xyz']['sha256'] != task['source_xyz']['sha256']: raise ValueError('Source byte identity failed')
        header = '! r2SCAN-3c NoAutostart DefGrid3 TightSCF'+(' CPCM(Water)' if task['medium']=='water' else '')
        q = -3 if task['metal']=='Ca' else -2
        if task['charge'] != q or task['multiplicity'] != 1 or inp.read_text() != header+f'\n* xyzfile {q} 1 core.xyz\n': raise ValueError('Changed scientific policy')
        rows = atoms[(task['sample'],task['metal'])]; composition = dict(Counter(r[0] for r in rows))
        full = sum(atomic_numbers[r[0]] for r in rows)-q
        effective = full-(46 if task['metal']=='La' else 0)
        if effective % 2: raise ValueError('Singlet electron parity inconsistent')
        records.append({'task_id':task['task_id'],'composition':composition,'all_electrons':full,'explicit_electrons':effective})
    mpi = shutil.which('mpirun')
    observed = {'mpi':pin(mpi), 'mpi_version':subprocess.check_output([mpi,'--version'],text=True),
                'hostname':socket.gethostname(), 'affinity':sorted(os.sched_getaffinity(0)),
                'environment':{k:v for k,v in os.environ.items() if k.startswith(('SLURM_','OMP','OMPI','MKL','OPENBLAS')) or k in ('PATH','LD_LIBRARY_PATH','LD_PRELOAD')}}
    observed['ldd'] = {name:subprocess.check_output(['ldd',str(Path(m['orca']['path']).parent/name)],text=True)
                       for name in ('orca','orca_leanscf_mpi')}
    libraries = set()
    for output in observed['ldd'].values():
        libraries.update(re.findall(r'(/\S+)\s+\(0x',output))
    observed['libraries'] = [pin(x) for x in sorted(libraries)]
    old = check(m['parser_fixtures'][0]).read_text(); trah = check(m['parser_fixtures'][1]).read_text()
    fixture_checks = {'native_Mulliken139':population(old,'MULLIKEN ATOMIC CHARGES',139) is not None,
                      'native_Loewdin139':population(old,'LOEWDIN ATOMIC CHARGES',139) is not None,
                      'actual_Trah_macro':bool(trah_macro(trah)) and numeric(trah_macro(trah)[-1][2])==5.528534e-7}
    if not all(fixture_checks.values()): raise ValueError('Actual archived parser fixture failed: '+str(fixture_checks))
    return {'status':'passed','manifest':pin(mp),'memory_policy':policy,'electronic_inventory':records,'runtime_observed':observed,
            'archived_parser_checks':fixture_checks}

def numeric(x): return float(x.replace('D','E').replace('d','e'))

def trah_macro(text):
    return re.findall(r'^\s*(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')[^\n]*\((?:NR|AH|TRAH)\s+MAcro\)',text,re.M)

def population(text, title, n):
    if title not in text: return None
    section = text.rsplit(title,1)[1]; rows = []
    for line in section.splitlines()[1:]:
        hit = re.match(r'^\s*(\d+)\s+([A-Za-z]+)\s*:\s*('+NUMBER+r')',line)
        if hit: rows.append({'index':int(hit[1]),'element':hit[2],'charge':numeric(hit[3])})
        elif rows: break
    return rows if [r['index'] for r in rows] == list(range(n)) else None

def collect(m, mp, pf):
    inventory = {r['task_id']:r for r in pf['electronic_inventory']}; rows = []
    for t in m['tasks']:
        row = {k:t[k] for k in ('task_id','sample','metal','medium','charge','multiplicity')}
        out = Path(t['output_path']); receipt = out.with_name(out.name+'.execution.json')
        if not out.exists() or not receipt.exists(): row['status']='not_completed'; rows.append(row); continue
        try:
            r = json.loads(receipt.read_text()); text = out.read_text(errors='replace'); row['execution']=pin(receipt); row['output']=pin(out)
            energy = re.findall(r'FINAL SINGLE POINT ENERGY\s+('+NUMBER+r')',text)
            nel = re.findall(r'Number of Electrons\s+NEL\s+\.+\s+(\d+)',text)
            row['energy_hartree'] = numeric(energy[-1]) if energy else None
            row['printed_explicit_electrons'] = int(nel[-1]) if nel else None
            row['scf_summary'] = [{'criterion':a.strip(),'value':numeric(b),'tolerance':numeric(c)} for a,b,c in
                re.findall(r'^\s*Last (.*?)\s+\.{3}\s*('+NUMBER+r')\s+Tolerance\s*:\s*('+NUMBER+r')',text,re.M)]
            macro = trah_macro(text)
            row['final_solver'] = 'TRAH' if macro else 'DIIS/SOSCF (see raw iteration table)'
            row['final_trah_macro'] = {'iteration':int(macro[-1][0]),'energy':numeric(macro[-1][1]),'orbital_residual':numeric(macro[-1][2])} if macro else None
            row['scf_policy_lines'] = [x.strip() for x in text.splitlines() if any(s in x for s in ('ConvCheckMode','TolE ','TolG ','TolErr ','ConvForced'))]
            row['raw_solver_lines'] = [x for x in text.splitlines() if re.match(r'^\s*\d+\s+',x) and ('MAcro' in x or re.match(r'^\s*\d+\s+-\d+\.\d+\s+',x))]
            row['populations'] = {kind:population(text,title,192) for kind,title in [('Mulliken','MULLIKEN ATOMIC CHARGES'),('Loewdin','LOEWDIN ATOMIC CHARGES')]}
            source = next(s for s in m['sources'] if s['sample']==t['sample'])
            mapping = json.loads(Path(source['source_preparation']['path']).read_text())['mapping']['source_to_qm']
            groups = {'metal':[0], 'PQQ':[x['qm_index'] for x in mapping if x['kind']=='opaque_cofactor'],
                      'caps':[x['qm_index'] for x in mapping if x['kind']=='sigma_link_H'],
                      'protein':[x['qm_index'] for x in mapping if x['kind']=='source']}
            row['group_indices'] = groups; row['group_charges'] = {}
            for kind, values in row['populations'].items():
                if values is not None: row['group_charges'][kind] = {g:sum(values[i]['charge'] for i in indices) for g,indices in groups.items()}
            row['CPCM_dielectric_hartree'] = [numeric(x) for x in re.findall(r'CPCM Dielectric\s*:\s*('+NUMBER+r')',text)]
            row['dipole_lines'] = [x.strip() for x in text.splitlines() if 'dipole' in x.lower() or 'Debye' in x]
            row['basis_ecp_lines'] = [x.strip() for x in text.splitlines() if any(s in x for s in ('def2-mTZVPP','replacing 46','Def2-ECP','ECP PARAMETER'))]
            row['wall_memory_lines'] = [x.strip() for x in text.splitlines() if 'TOTAL RUN TIME' in x or 'Maximum memory' in x]
            retained = [p for p in out.parent.glob('endpoint.runtime*') if p.is_file() and (p.suffix in ('.gbw','.densities','.json') or 'property' in p.name)]
            row['retained_artifacts'] = [pin(p) for p in sorted(retained)]
            checks = {'native_receipt_success':r['returncode']==0 and r['normal_termination'] and r['scf_converged'],
                'manifest_identity':r['manifest']['sha256']==pin(mp)['sha256'],
                'output_identity':r['artifacts']['output']['sha256']==pin(out)['sha256'],
                'version':r['orca_version']=='6.1.1', 'energy_present':row['energy_hartree'] is not None,
                'electron_count':row['printed_explicit_electrons']==inventory[t['task_id']]['explicit_electrons'],
                'populations_complete':all(v is not None for v in row['populations'].values()),
                'gbw_retained':any(p.suffix=='.gbw' for p in retained),
                'native_basis':'def2-mTZVPP' in text,
                'La_ECP':t['metal']!='La' or bool(re.search(r'replacing\s+46\s+core',text))}
            if macro:
                tol = next((v['tolerance'] for v in row['scf_summary'] if v['criterion']=='Orbital Gradient'),None)
                checks['final_TRAH_gradient'] = tol is not None and row['final_trah_macro']['orbital_residual'] <= tol
            row['checks']=checks; row['status']='complete_diagnostic' if all(checks.values()) else 'failed_collection'
            row['density_response_qualification']='unestablished; solver summaries retained, no same-state repeat or derivative qualification'
        except Exception as e: row['status']='collection_error'; row['error']=f'{type(e).__name__}: {e}'
        rows.append(row)
    return {'schema':'pede_fixed_origin_dft_collection.v1','manifest':pin(mp),'cells':rows,
            'complete_cells':sum(r['status']=='complete_diagnostic' for r in rows),
            'declared_cells':8,'scientific_scope':'fixed-origin solvent response only; no fit, labels, forces or production promotion'}

def main():
    p=argparse.ArgumentParser(); p.add_argument('action',choices=('preflight','execute','collect')); p.add_argument('manifest',type=Path)
    p.add_argument('--manifest-sha',required=True); p.add_argument('tasks',nargs='*'); a=p.parse_intermixed_args()
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('All preflight, execution and collection require real Slurm')
    mp=a.manifest.resolve(); check({'path':str(mp),'sha256':a.manifest_sha}); m=json.loads(mp.read_text())
    if m['implementation']['REFERENCE.py']['snapshot'] != pin(__file__): raise ValueError('Adapter identity mismatch')
    w=mp.parent; attempt=w/'attempts'/f"{os.environ['SLURM_JOB_ID']}_{a.action}"; attempt.mkdir(parents=True,exist_ok=False)
    status={'status':'started','manifest':pin(mp),'action':a.action,'tasks':a.tasks,'started_utc':dt.datetime.now(dt.timezone.utc).isoformat()}
    try:
        with (w/'RUN.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            pf=preflight(m,mp); write(attempt/'PREFLIGHT.json',pf)
            if a.action=='execute':
                known={t['task_id'] for t in m['tasks']}
                if not a.tasks or not set(a.tasks)<=known or len(a.tasks)!=len(set(a.tasks)): raise ValueError('Explicit unique frozen task IDs required')
                from run_orca_task_manifest import run_manifest
                for tid in a.tasks:
                    run_manifest(mp,orca_path=Path(m['orca']['path']),workers=1,nprocs=12,selected_tasks={tid},
                                 expected_runner_sha256=m['execution_policy']['task_runner']['sha256'])
                    result=collect(m,mp,pf); write(attempt/(tid+'_COLLECTION.json'),result)
                    cell=next(r for r in result['cells'] if r['task_id']==tid)
                    if cell['status']!='complete_diagnostic': raise RuntimeError('Stop before next cell: '+tid+' '+cell['status'])
            if a.action!='preflight': write(attempt/'COLLECTION.json',collect(m,mp,pf))
            status['status']='complete'
    except BaseException as e:
        status['status']='failed'; status['error']=f'{type(e).__name__}: {e}'; status['traceback']=traceback.format_exc()
        raise
    finally:
        status['finished_utc']=dt.datetime.now(dt.timezone.utc).isoformat(); write(attempt/'EXECUTION.json',status)

if __name__=='__main__': main()
