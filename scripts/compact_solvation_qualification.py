"""Fixed-state SCF qualification of the existing compact solvent correction."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import shutil

from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new
from affordable_workflow import dry_run, execute as execute_manifest
from compact_solvation import completed, diagnostics, METHOD

PROTOCOL = 'compact_GFN2_same_state_numerical_qualification_v1'
CASES = ('1H4I', 'q9z4j7-pqq-la_model', '4MAE', 'a0acd6b9f2-pqq-la_model', '1F6S', '1GLG')
VARIANTS = {'native_tight_fresh': CASES, 'ordinary_tight_seeded': CASES[:2],
            'ordinary_tight_gbw': CASES[:2], 'ordinary_explicit_gbw': CASES[:2]}
SCOPES = {'primary_checks': ('native_tight_fresh','ordinary_tight_seeded'),
          'orbital_restart_fix': ('ordinary_tight_gbw',),
          'explicit_guess_fix': ('ordinary_explicit_gbw',)}


def recipe(charge, medium, variant):
    if variant not in VARIANTS or medium not in ('vacuum', 'alpb'):
        raise InvalidArtifact('unsupported numerical recipe')
    seeded = variant != 'native_tight_fresh'
    gbw = variant in ('ordinary_tight_gbw','ordinary_explicit_gbw')
    explicit = variant == 'ordinary_explicit_gbw'
    return ('! Native-GFN2-xTB' + ('' if explicit else ' MORead NoAutostart' if gbw else '' if seeded else ' NoAutostart')
            + (' ALPB(Water)' if medium == 'alpb' else '')
            + ('\n%moinp "seed_source.gbw"' if gbw and not explicit else '')
            + '\n%maxcore 2000\n%method\n WriteXTBParam true\n ReadXTBParam false\nend\n'
            + '%scf\n SmearTemp 300\n Convergence Tight\n UseXTBMixer '
            + ('false' if seeded else 'true')
            + ('\n Guess MORead\n MOInp "seed_source.gbw"' if explicit else '') + '\nend\n'
            + f'* xyzfile {charge} 1 core.xyz\n')


def prepare(primary_collection, agreement, output, stage='primary_checks'):
    if stage not in SCOPES: raise InvalidArtifact('undeclared qualification stage')
    primary_collection = Path(primary_collection).resolve()
    col = read_json(primary_collection)
    if col['status'] != 'complete':
        raise InvalidArtifact('completed primary panel required')
    parent = read_json(verify(col['manifest']))
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = out/'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        dst = impl/p.name; shutil.copyfile(p, dst); pins[p.name] = record(dst)
    tasks = []
    for variant in SCOPES[stage]:
        cases = VARIANTS[variant]
        for case in cases:
            for metal in ('Ca', 'La'):
                row = next(r for r in col['rows'] if (r['case_id'], r['representation'], r['metal']) == (case, 'context', metal))
                for medium in ('vacuum', 'alpb'):
                    endpoint = row['endpoints'][medium]
                    original = next(t for t in parent['all_tasks'] if (t['case_id'], t['representation'], t['metal'], t['medium']) == (case, 'context', metal, medium))
                    tid = f'{case}__{metal}__{medium}__{variant}'
                    d = out/'tasks'/tid; d.mkdir(parents=True)
                    xp = d/'core.xyz'; shutil.copyfile(verify(original['xyz']), xp)
                    ip = d/'endpoint.inp'; ip.write_text(recipe(original['charge'], medium, variant))
                    seed = None
                    if variant != 'native_tight_fresh':
                        kind = 'gbw' if variant in ('ordinary_tight_gbw','ordinary_explicit_gbw') else 'xtbw'
                        source = verify(endpoint['output']).parent/('endpoint.runtime.'+kind)
                        dest = d/('seed_source.'+kind); shutil.copyfile(source, dest)
                        seed = {'source': record(source), 'immutable_copy': record(dest),
                                'kind':kind, 'runtime_path': None if kind=='gbw' else str(d/'endpoint.runtime.xtbw')}
                    tasks.append({'task_id': tid, 'case_id': case, 'case': case, 'representation': 'context',
                                  'metal': metal, 'medium': medium, 'variant': variant, 'charge': original['charge'],
                                  'multiplicity': 1, 'input': record(ip), 'xyz': record(xp),
                                  'output_path': str(d/'endpoint.out'), 'primary_task': original,
                                  'primary_endpoint': endpoint, 'restart': seed})
    m = {'protocol_id': PROTOCOL, 'stage': stage, 'method_id': METHOD, 'agreement': record(agreement),
         'primary_collection': record(primary_collection), 'primary_manifest': col['manifest'],
         'orca': parent['orca'], 'implementation': pins, 'tasks': tasks,
         'execution_resources': {'mpi_ranks': 8, 'concurrent_tasks': 8},
         'execution_policy': {'task_runner': pins['run_orca_task_manifest.py'],
                              'runtime_renderer': pins['render_orca_runtime_input.py']},
         'tolerances': {'endpoint_transfer_kcal_mol': .10, 'score_correction_kcal_mol': .20},
         'baseline_changed': False, 'scientific_state_changed': False, 'new_DFT_MACE_calls': 0}
    write_new(out/'manifest.json', m)
    return validate(out/'manifest.json')


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or m['method_id'] != METHOD:
        raise InvalidArtifact('qualification method changed')
    col = read_json(verify(m['primary_collection'])); parent = read_json(verify(m['primary_manifest']))
    verify(m['agreement']); verify(m['orca'])
    for pin in m['implementation'].values(): verify(pin)
    stage = m.get('stage','primary_checks')
    if stage not in SCOPES: raise InvalidArtifact('undeclared qualification stage')
    expected = {(v, c, z, s) for v in SCOPES[stage] for c in VARIANTS[v] for z in ('Ca','La') for s in ('vacuum','alpb')}
    actual = {(t['variant'], t['case_id'], t['metal'], t['medium']) for t in m['tasks']}
    if actual != expected or len(m['tasks']) != len(expected):
        raise InvalidArtifact('fixed numerical scope changed')
    for t in m['tasks']:
        original = next(o for o in parent['all_tasks'] if o['task_id'] == t['primary_task']['task_id'])
        row = next(r for r in col['rows'] if (r['case_id'],r['representation'],r['metal']) == (t['case_id'],'context',t['metal']))
        endpoint = row['endpoints'][t['medium']]
        if original != t['primary_task'] or endpoint != t['primary_endpoint']:
            raise InvalidArtifact('primary source changed')
        if t['charge'] != original['charge'] or t['multiplicity'] != original['multiplicity']:
            raise InvalidArtifact('electronic state changed')
        if verify(t['xyz']).read_bytes() != verify(original['xyz']).read_bytes():
            raise InvalidArtifact('coordinates changed')
        if verify(t['input']).read_text() != recipe(t['charge'], t['medium'], t['variant']):
            raise InvalidArtifact('numerical recipe changed')
        for name in ('output', 'receipt', 'manifest'): verify(endpoint[name])
        if t['variant'] != 'native_tight_fresh':
            if not t['restart'] or verify(t['restart']['source']).read_bytes() != verify(t['restart']['immutable_copy']).read_bytes():
                raise InvalidArtifact('restart source changed')
        elif t['restart'] is not None:
            raise InvalidArtifact('unrecorded restart')
    dry_run(manifest)
    return {'status': 'validated', 'manifest': record(manifest), 'tasks': len(m['tasks'])}


def execute(manifest):
    validate(manifest); m = read_json(manifest)
    # The immutable source survives ORCA rewriting its active restart file.
    for t in m['tasks']:
        if not t['restart'] or t['restart'].get('kind')=='gbw' or Path(t['output_path']).exists(): continue
        r = t['restart']; runtime = Path(r['runtime_path'])
        if runtime.exists():
            if runtime.read_bytes() != verify(r['immutable_copy']).read_bytes():
                raise InvalidArtifact('unverified preexisting restart')
        else: shutil.copyfile(verify(r['immutable_copy']), runtime)
    return execute_manifest(manifest)


def restart_diagnostics(text, variant):
    guesses = re.findall(r'^INITIAL GUESS:\s*(\S+)',text,re.M)
    guess = guesses[-1] if guesses else None
    ignored = 'MOInp will be ignored' in text
    requested = variant != 'native_tight_fresh'
    orbital = variant in ('ordinary_tight_gbw','ordinary_explicit_gbw')
    return {'observed_initial_guess':guess, 'ignored_moinp_warning':ignored,
            'restart_status':'not_requested' if not requested else
                'orbital_restart_confirmed' if orbital and guess=='MOREAD' and not ignored else
                'requested_restart_not_confirmed'}


def collect(manifest, output):
    validate(manifest); m = read_json(manifest); endpoints = []; comparisons = []
    for t in m['tasks']:
        pin = completed(manifest, t['task_id']); op = Path(t['output_path'])
        row = {k: t[k] for k in ('task_id','variant','case_id','metal','medium')}
        row.update(status='unavailable', energy_hartree=None, reason=None, primary=t['primary_endpoint'],
                   artifacts=[record(p) for p in (op,Path(str(op)+'.execution.json')) if p.exists()])
        text = op.read_text() if op.exists() else ''
        row.update(restart_diagnostics(text,t['variant']))
        if pin:
            try:
                audit = diagnostics(pin,t)
                if t['variant']!='native_tight_fresh' and row['restart_status']!='orbital_restart_confirmed':
                    raise InvalidArtifact('requested restart not confirmed in actual output')
                if audit['charge_sanity_status'] != 'pass':
                    raise InvalidArtifact('atomic charge outside declared diagnostic range')
                tol = re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text)
                observed = float(tol[-1]) if len(tol)==1 else None
                params = read_json(verify(audit['parameter_export']))
                if params != read_json(verify(t['primary_endpoint']['parameter_export'])):
                    raise InvalidArtifact('GFN2 parameters changed')
                mixer = 'INFO: Using special xTB SCF mixer' in text
                if mixer != (t['variant']=='native_tight_fresh'):
                    raise InvalidArtifact('actual solver differs')
                if observed is None or observed > 1e-8:
                    raise InvalidArtifact('requested tighter convergence not confirmed')
                cycles = re.findall(r'SCF CONVERGED AFTER\s+(\d+) CYCLES',text)
                restart_lines = [line.strip() for line in text.splitlines() if re.search(r'xtbw|restart|reading.*(?:charge|multipol)',line,re.I)]
                row.update(status='complete', **pin, **audit, observed_TolE_hartree=observed,
                           cycles=int(cycles[-1]) if cycles else None, restart_output_lines=restart_lines,
                           endpoint_difference_kcal_mol=(pin['energy_hartree']-t['primary_endpoint']['energy_hartree'])*HA_TO_KCAL)
            except (InvalidArtifact, OSError, KeyError, IndexError) as exc:
                row.update(status='unsupported', reason=str(exc))
        endpoints.append(row)
    by = {(r['variant'],r['case_id'],r['metal'],r['medium']):r for r in endpoints}
    variants = SCOPES[m.get('stage','primary_checks')]
    for variant in variants:
        cases = VARIANTS[variant]
        for case in cases:
            eps = {z:{s:by[(variant,case,z,s)] for s in ('vacuum','alpb')} for z in ('Ca','La')}
            good = all(r['status']=='complete' for v in eps.values() for r in v.values())
            row = {'variant':variant,'case_id':case,'status':'unavailable','transfer_differences_kcal_mol':None,'score_correction_difference_kcal_mol':None}
            if good:
                diffs = {z:((v['alpb']['energy_hartree']-v['vacuum']['energy_hartree'])-(v['alpb']['primary']['energy_hartree']-v['vacuum']['primary']['energy_hartree']))*HA_TO_KCAL for z,v in eps.items()}
                delta = diffs['Ca']-diffs['La']
                row.update(status='pass' if max(abs(x) for x in diffs.values())<=.10 and abs(delta)<=.20 else 'fail',
                           transfer_differences_kcal_mol=diffs,score_correction_difference_kcal_mol=delta)
            comparisons.append(row)
    result = {'protocol_id':PROTOCOL,'manifest':record(manifest),'collector_implementation':record(__file__),
              'endpoints':endpoints,'comparisons':comparisons,
              'complete_endpoints':sum(r['status']=='complete' for r in endpoints),'endpoint_denominator':len(m['tasks']),
              'variant_status':{v:'pass' if all(r['status']=='pass' for r in comparisons if r['variant']==v) else 'failed_or_unavailable' for v in variants},
              'seed_restart_interpretation':'actual restart output must be reviewed; seeded checks only concern same-basin agreement',
              'original_unseeded_failures_preserved':True,'baseline_changed':False,'calibration_changed':False}
    write_new(output,result); return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare')
    for k in ('primary_collection','agreement','output'):q.add_argument('--'+k.replace('_','-'),required=True)
    q.add_argument('--stage',choices=tuple(SCOPES),default='primary_checks')
    for name in ('validate','execute','collect'):
        q=s.add_parser(name);q.add_argument('--manifest',required=True)
        if name=='collect':q.add_argument('--output',required=True)
    a=vars(p.parse_args());op=a.pop('operation');r=globals()[op](**a)
    print(json.dumps({k:v for k,v in r.items() if k not in ('endpoints','comparisons')},indent=2))

if __name__=='__main__':main()
