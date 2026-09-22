"""Three fixed real-source native GFN2 iteration-limit checks; no score overlay."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil

from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new
from affordable_workflow import dry_run, execute as execute_manifest
from compact_solvation import completed, diagnostics, input_text

PROTOCOL = 'adaptive_pool_native_GFN2_MaxIter500_diagnostic_v1'
IDS = ('4MAE__La__at_adaptive_Ca__vacuum', '1H4I__La__origin__vacuum',
       '4MAE__La__origin__vacuum')


def extended_input(task):
    original = verify(task['input']).read_text()
    if original != input_text(task['charge'], task['multiplicity'], 'vacuum', 'native'):
        raise InvalidArtifact('source differs from unchanged primary native vacuum recipe')
    return original.replace('%scf\n', '%scf\n MaxIter 500\n')


def prepare(failed_manifest, pilot_manifest, agreement, output):
    failed_manifest = Path(failed_manifest).resolve()
    pilot_manifest = Path(pilot_manifest).resolve()
    fm = read_json(failed_manifest); pm = read_json(pilot_manifest)
    failed = next(t for t in fm['tasks'] if t['task_id'] == IDS[0])
    if completed(failed_manifest, failed['task_id']) is not None:
        raise InvalidArtifact('selected failure has a successful endpoint')
    failure_output = Path(failed['output_path'])
    if 'SCF NOT CONVERGED AFTER 125 CYCLES' not in failure_output.read_text():
        raise InvalidArtifact('expected actual 125-cycle failure not present')
    sources = [(IDS[0], failed_manifest, failed, None)]
    for case, tid in zip(('1H4I', '4MAE'), IDS[1:]):
        c = next(c for c in pm['cases'] if c['case_id'] == case)
        endpoint = c['matrix']['La']['origin']['low']['vacuum']
        mp = verify(endpoint['manifest']); sm = read_json(mp)
        source = next(t for t in sm.get('all_tasks', sm['tasks']) if t['task_id'] == endpoint['task_id'])
        if completed(mp, source['task_id']) != endpoint:
            raise InvalidArtifact('exact successful control receipt not reproducible')
        sources.append((tid, mp, source, endpoint))
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = out/'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        dst = impl/p.name; shutil.copyfile(p, dst); pins[p.name] = record(dst)
    tasks = []
    for tid, mp, source, endpoint in sources:
        d = out/'tasks'/tid; d.mkdir(parents=True)
        xp = d/'core.xyz'; shutil.copyfile(verify(source['xyz']), xp)
        ip = d/'endpoint.inp'; ip.write_text(extended_input(source))
        op = Path(source['output_path'])
        tasks.append({'task_id': tid, 'case_id': tid.split('__')[0], 'case': tid.split('__')[0],
                      'metal': 'La', 'medium': 'vacuum', 'charge': source['charge'],
                      'multiplicity': source['multiplicity'], 'xyz': record(xp), 'input': record(ip),
                      'output_path': str(d/'endpoint.out'), 'source_manifest': record(mp),
                      'source_task': source, 'source_endpoint': endpoint, 'source_output': record(op),
                      'source_receipt': record(str(op)+'.execution.json'),
                      'source_parameters': record(op.parent/'endpoint.runtime.xtb.json')})
    m = {'protocol_id': PROTOCOL, 'agreement': record(agreement), 'pilot_manifest': record(pilot_manifest),
         'orca': fm['orca'], 'implementation': pins, 'tasks': tasks,
         'execution_resources': {'mpi_ranks': 8, 'concurrent_tasks': 3},
         'execution_policy': {'task_runner': pins['run_orca_task_manifest.py'],
                              'runtime_renderer': pins['render_orca_runtime_input.py']},
         'control_tolerance_kcal_mol': 0.05, 'score_overlay': None,
         'new_high_level_calls': 0, 'production_changed': False}
    write_new(out/'manifest.json', m)
    return validate(out/'manifest.json')


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or tuple(t['task_id'] for t in m['tasks']) != IDS:
        raise InvalidArtifact('fixed diagnostic membership changed')
    if m['control_tolerance_kcal_mol'] != .05 or m['execution_resources'] != {'mpi_ranks':8,'concurrent_tasks':3}:
        raise InvalidArtifact('frozen tolerance/resources changed')
    verify(m['pilot_manifest'])
    for pin in m['implementation'].values(): verify(pin)
    for t in m['tasks']:
        parent = read_json(verify(t['source_manifest']))
        source = next(s for s in parent.get('all_tasks', parent['tasks']) if s['task_id']==t['source_task']['task_id'])
        if source != t['source_task'] or (t['charge'], t['multiplicity']) != (source['charge'], source['multiplicity']):
            raise InvalidArtifact('source task/state changed')
        if verify(t['xyz']).read_bytes() != verify(source['xyz']).read_bytes():
            raise InvalidArtifact('coordinates changed')
        if verify(t['input']).read_text() != extended_input(source):
            raise InvalidArtifact('input changes beyond MaxIter500')
        for k in ('source_output','source_receipt','source_parameters'): verify(t[k])
        if t['source_endpoint'] and completed(verify(t['source_manifest']),source['task_id']) != t['source_endpoint']:
            raise InvalidArtifact('source control completion changed')
    return dry_run(manifest)


def observed(text):
    patterns = {'MaxIter': r'Maximum # iterations\s+MaxIter\s+\.{4}\s+(\d+)',
                'TolE': r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',
                'charge': r'Total Charge\s+Charge\s+\.{4}\s+(-?\d+)',
                'multiplicity': r'Multiplicity\s+Mult\s+\.{4}\s+(\d+)',
                'electrons': r'Number of Electrons\s+NEL\s+\.{4}\s+(\d+)',
                'dimension': r'Basis Dimension\s+Dim\s+\.{4}\s+(\d+)',
                'initial_guess': r'^INITIAL GUESS:\s*(\S+)'}
    d = {k: re.findall(p,text,re.M) for k,p in patterns.items()}
    if any(len(v)!=1 for v in d.values()): raise InvalidArtifact('missing/ambiguous native SCF settings')
    d = {k: v[0] for k,v in d.items()}
    d['special_mixer'] = 'INFO: Using special xTB SCF mixer' in text
    return d


def collect(manifest, output):
    validate(manifest); m = read_json(manifest); rows = []
    for t in m['tasks']:
        op = Path(t['output_path']); rp = Path(str(op)+'.execution.json')
        row = {'task_id':t['task_id'], 'status':'unavailable', 'energy_hartree':None,
               'source_endpoint':t['source_endpoint'], 'difference_kcal_mol':None,
               'control_pass':None, 'reason':None}
        text = op.read_text() if op.exists() else ''
        cycles = re.findall(r'SCF (?:NOT )?CONVERGED AFTER\s+(\d+) CYCLES',text)
        row['SCF_cycles'] = int(cycles[-1]) if cycles else None
        row['artifacts'] = [record(p) for p in (op,rp) if p.exists()]
        row['execution_receipt'] = read_json(rp) if rp.exists() else None
        pin = completed(manifest,t['task_id'])
        if pin:
            try:
                audit = diagnostics(pin,t); settings = observed(text)
                original = observed(verify(t['source_output']).read_text())
                if settings.pop('MaxIter') != '500' or original.pop('MaxIter') != '125' or settings != original:
                    raise InvalidArtifact('actual SCF/state settings differ beyond iteration limit')
                if settings['initial_guess'] != 'SAD' or not settings['special_mixer']:
                    raise InvalidArtifact('intended fresh native driver not active')
                if read_json(verify(audit['parameter_export'])) != read_json(verify(t['source_parameters'])):
                    raise InvalidArtifact('native parameter exports changed')
                row.update(status='complete', **pin, **audit, observed_settings=settings,
                           observed_MaxIter=500, parameters_identical=True)
                if t['source_endpoint']:
                    delta = (pin['energy_hartree']-t['source_endpoint']['energy_hartree'])*HA_TO_KCAL
                    row.update(difference_kcal_mol=delta, control_pass=abs(delta)<=.05)
            except (InvalidArtifact, OSError, KeyError, ValueError) as exc:
                row.update(status='unsupported',reason=str(exc))
        elif op.exists(): row['reason']='native endpoint not converged or receipt invalid'
        rows.append(row)
    result = {'protocol_id':PROTOCOL, 'manifest':record(manifest), 'rows':rows,
              'complete':sum(r['status']=='complete' for r in rows), 'denominator':3,
              'both_controls_pass':all(r['control_pass'] is True for r in rows[1:]),
              'formerly_failed_cell_complete':rows[0]['status']=='complete',
              'recovery_overlay_applied':False, 'production_changed':False}
    write_new(output,result); return result


def execute(manifest):
    validate(manifest)
    return execute_manifest(manifest)


def main():
    p=argparse.ArgumentParser(description=__doc__); sp=p.add_subparsers(dest='op',required=True)
    q=sp.add_parser('prepare')
    for k in ('failed_manifest','pilot_manifest','agreement','output'): q.add_argument('--'+k.replace('_','-'),required=True)
    for name in ('validate','execute','collect'):
        q=sp.add_parser(name); q.add_argument('--manifest',required=True)
        if name=='collect': q.add_argument('--output',required=True)
    args=vars(p.parse_args()); op=args.pop('op'); r=globals()[op](**args)
    print(json.dumps({k:v for k,v in r.items() if k!='rows'},indent=2))


if __name__=='__main__': main()
