"""Finite continuation/collection of the nine declared LanM occupancy states.

Uses the immutable scoring implementation and existing ORCA runner. No additional
structures, chemical states, retries, thresholds, optimization rounds or models.
"""
import argparse
from email.message import EmailMessage
from email.utils import formatdate, make_msgid
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
W = ROOT/'workspaces/lanm_global_occupancy_20260923'
IMPL = W/'scoring_v2/implementation'
sys.path.insert(0, str(IMPL))
import lanm_global_occupancy as score
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from compact_solvation import completed, diagnostics
from structure_informed_starts import scf_details

PARENT = W/'scoring_v2/manifest.json'
FIRST_NATIVE = W/'native_feasibility_v2/manifest.json'
FIRST_MACE = W/'mace_feasibility_v2/result.json'
ALL_MACE = W/'mace_merged_v2/result.json'
FINAL_NATIVE = W/'native_matrix_v2/manifest.json'
DIAG = ROOT/'diagnostics/lanm_global_occupancy_20260923'


def gate():
    m = score.validate(PARENT)
    r = read_json(FIRST_MACE)
    if r['manifest'] != record(PARENT) or r['capability_only']:
        raise InvalidArtifact('first complete-search result required')
    score.validate_mace_rows(m, r)
    if r['selected_states'] != ['Hans_8DQ2__EF12'] or r['rows'][0]['status'] != 'complete':
        raise InvalidArtifact('first whole-chain MACE matrix is unavailable')
    qualification = read_json(W/'mace_feasibility_v2/Hans_8DQ2__EF12/adapter_qualification.json')
    if not qualification['passed']:
        raise InvalidArtifact('actual full-chain adapter check failed')
    nm = read_json(FIRST_NATIVE)
    if nm['parent_manifest'] != record(PARENT) or len(nm['tasks']) != 4:
        raise InvalidArtifact('first native denominator differs')
    checks = []
    for t in nm['tasks']:
        pin = completed(FIRST_NATIVE, t['task_id'])
        if pin is None:
            raise InvalidArtifact('first native pair is unavailable: '+t['task_id'])
        audit = diagnostics(pin, t)
        scf = scf_details(verify(pin['output']).read_text())
        if (audit['charge_sanity_status'] != 'pass' or scf['energy']['tolerance'] != 1e-10
                or not scf['native_mixer_observed']):
            raise InvalidArtifact('first native numerical/state check failed')
        checks.append({'task_id': t['task_id'], 'actual': pin, 'scf': scf,
                       'charge_sum_e': audit['charge_sum_e']})
    return {'status': 'whole_chain_feasibility_pass', 'first_MACE': record(FIRST_MACE),
            'first_native': record(FIRST_NATIVE), 'checks': checks, 'adapter_qualification': qualification}


def propose_rest():
    passed = gate()
    write_new(W/'CONTINUATION_GATE.json', passed)
    m = score.validate(PARENT)
    selected = [read_json(verify(p))['state_id'] for p in m['states']
                if read_json(verify(p))['state_id'] != 'Hans_8DQ2__EF12']
    if len(selected) != 8:
        raise InvalidArtifact('exact remaining eight states required')
    score.run_mace(PARENT, W/'mace_remaining_v2', selected=selected)
    first, rest = read_json(FIRST_MACE), read_json(W/'mace_remaining_v2/result.json')
    result = {'protocol_id': score.PROTOCOL, 'manifest': record(PARENT),
        'capability_only': False, 'selected_states': first['selected_states']+rest['selected_states'],
        'rows': first['rows']+rest['rows'], 'constituent_results': [record(FIRST_MACE), record(W/'mace_remaining_v2/result.json')],
        'molecular_calls_attempted': first['molecular_calls_attempted']+rest['molecular_calls_attempted'],
        'wall_seconds': first['wall_seconds']+rest['wall_seconds']}
    score.validate_mace_rows(m, result)
    write_new(ALL_MACE, result)
    score.native_prepare(PARENT, FINAL_NATIVE.parent, ALL_MACE, ranks=8, workers=4, reuse=[FIRST_NATIVE])
    # Submission is of this already authorized finite matrix only.
    command = ['sbatch', '--parsable', '--output='+str(DIAG/'native_matrix_%j.out'),
               '--error='+str(DIAG/'native_matrix_%j.err'), str(DIAG/'run_native_matrix.sbatch')]
    response = subprocess.run(command, text=True, capture_output=True)
    write_new(W/'NATIVE_SUBMISSION.json', {'command': command, 'returncode': response.returncode,
              'stdout': response.stdout, 'stderr': response.stderr, 'manifest': record(FINAL_NATIVE)})
    if response.returncode:
        raise InvalidArtifact('native matrix submission failed; prepared matrix retained')
    print(response.stdout.strip())


def collect_and_report():
    mr = read_json(ALL_MACE); nm = read_json(FINAL_NATIVE)
    expected = {(r['state_id'], c['candidate'], c['metal']): c
                for r in mr['rows'] for c in r['cells']}
    task_pairs = [(t, FINAL_NATIVE) for t in nm['tasks']]
    for r in nm['reused']:
        mp = verify(r['actual']['manifest']); old = read_json(mp)
        t = next(t for t in old['tasks'] if t['task_id']==r['actual']['task_id'])
        task_pairs.append((t, mp))
    for t, mp in task_pairs:
        c = expected[t['state_id'], t['candidate'], t['metal']]
        if xyz(verify(t['xyz'])) != xyz(verify(c['xyz'])) or t['charge'] != c['charge']:
            raise InvalidArtifact('native/MACE physical coordinates or state differ')
    summary = score.collect(PARENT, ALL_MACE, [FINAL_NATIVE], W/'RESULT_v2.json')
    result = read_json(W/'RESULT_v2.json')
    lines = ['# Whole-chain conditional LanM pilot', '',
        'Complete monomers with fixed source waters; La/Dy competition on a common',
        'origin/La-proposal/Dy-proposal pool within each source and occupancy.',
        'The existing PQQ workflow is unchanged.', '',
        f"Available complete state pools: **{summary['available']}/{summary['states']}**.",
        f"Native failed/unavailable scalar cells: {summary['native_failure_count']}.", '',
        'Positive Hans-minus-Mex exchange contrast denotes stronger relative La',
        'preference for Hans under the same conditional occupancy. These values',
        'are model energy descriptors, not binding free energies or occupancy populations.', '',
        '|Hans source / occupied sites|Static relative La preference|Accommodated relative La preference|',
        '|---|---:|---:|']
    def fmt(v): return 'unavailable' if v is None else f'{v:+.3f}'
    for row in summary['comparisons']:
        lines.append(f"|{row['Hans_state']}|{fmt(row['static_relative_La_preference_kcal_mol'])}|{fmt(row['accommodated_relative_La_preference_kcal_mol'])}|")
    lines += ['', 'Units: model kcal/mol for exchange of the declared full ion inventory.', '',
        '## Structural response', '', '|State|La proposal|Dy proposal|Selected La / Dy|', '|---|---|---|---|']
    for r in result['rows']:
        proposals = r.get('proposals', {})
        def motion(z):
            p = proposals.get(z, {})
            g = p.get('geometry', {})
            value = f"{g['RMS_atom_displacement_A']:.3f} Å RMS" if 'RMS_atom_displacement_A' in g else p.get('reason', 'unavailable')
            return p.get('status', 'unavailable')+'; '+value
        selections = r.get('selections') or {}
        choice = ' / '.join(selections.get(z, {}).get('selected', 'unavailable') for z in score.METALS)
        lines.append(f"|{r['state_id']}|{motion('La')}|{motion('Dy')}|{choice}|")
    lines += ['', '## Limits', '',
        '- A small bound-state monomer pilot. No folding ensemble or Hans dimerization equilibrium.',
        '- Full atoms are retained; native MACE still has a finite local interaction range.',
        '- Dy maximum-spin multiplicities 11/21 are declared hypotheses; other couplings and spin–orbit effects were not tested.',
        '- Sources have different retained water inventories. No absolute cross-source energy pooling occurred.',
        '- Bounded vacuum proposals are selected using actual native GFN2 solvent endpoints. They are not composite stationary minima.',
        '- Actual H-coordinate radial repair precedes all matched calculations; proton inventories are fixed.',
        '- This is consumed development evidence for three structures, not a validated within-series classifier.', '',
        'Primary experimental context: [Mattocks et al., Nature 2023](https://www.nature.com/articles/s41586-023-05945-5).',
        'The known protein-level preference includes oligomerization and cooperative response.', '',
        f"Actual output: `{W/'RESULT_v2.json'}`.",
        f"Finite matrix and scalar receipts: `{FINAL_NATIVE.parent}`.", '']
    report = DIAG/'PILOT_REPORT.md'
    with report.open('x') as f: f.write('\n'.join(lines))
    vault = Path('/home/jwestrob/jwestrob/obsidian-vault/agent-captures/2026-09-23_Nikasha-whole-chain-LanM-occupancy-pilot.md')
    with vault.open('x') as f: f.write('\n'.join(lines))
    write_new(W/'REPORT_RECEIPT.json', {'result': record(W/'RESULT_v2.json'), 'report': record(report), 'vault': record(vault)})
    return summary, report


def notify(summary, report):
    msg = EmailMessage()
    msg['To'] = 'jacobwestroberts@gmail.com'; msg['From'] = 'jwestrob@biotite.berkeley.edu'
    msg['Subject'] = 'Nikasha: first whole-protein LanM occupancy pilot is back'
    msg['Date'] = formatdate(localtime=False); msg['Message-ID'] = make_msgid()
    body = ["Hey Jacob,", '',
        'The small full-chain La/Dy pilot has finished. We tested two Hans crystal structures and one Mex comparator, each with EF1+EF2, EF2+EF3, and four ions.', '',
        'The purpose was to let the protein accommodate either metal and see whether that improves the source-dependent selectivity seen with separate pockets.', '',
        f"{summary['available']} of {summary['states']} complete shared geometry pools are available; {summary['native_failure_count']} native scalar cells failed or remain unavailable.", '',
        'Matched Hans-minus-Mex exchange scores follow. Positive means greater relative La preference for Hans; static -> accommodated, in model kcal/mol:', '']
    for r in summary['comparisons']:
        def fmt(v): return 'unavailable' if v is None else f'{v:+.2f}'
        body.append(r['Hans_state']+': '+fmt(r['static_relative_La_preference_kcal_mol'])+' -> '+fmt(r['accommodated_relative_La_preference_kcal_mol']))
    body += ['', 'These are conditional monomer scores, not binding affinities or a prediction that two/four ions is the equilibrium occupancy. Dimerization, broader unfolding and alternate Dy spin couplings are not in this pilot. The PQQ workflow is unchanged.', '',
        'The detailed report is attached and copied to the vault. This is an automated collection of actual completed outputs; interpretation of which direction to pursue comes next.', '', 'Enjoy the walk!']
    msg.set_content('\n'.join(body))
    msg.add_attachment(report.read_bytes(), maintype='text', subtype='markdown', filename='LanM-pilot-report.md')
    eml = DIAG/'pilot_result_email.eml'
    with eml.open('xb') as f: f.write(msg.as_bytes())
    response = subprocess.run(['/usr/sbin/sendmail', '-t'], input=msg.as_bytes(), capture_output=True)
    write_new(DIAG/'pilot_result_email_receipt.json', {'to': msg['To'], 'subject': msg['Subject'],
        'message_id': msg['Message-ID'], 'returncode': response.returncode,
        'stdout': response.stdout.decode(errors='replace'), 'stderr': response.stderr.decode(errors='replace'),
        'status': 'accepted_by_local_mail_relay' if response.returncode==0 else 'relay_failed',
        'delivery_confirmation': None, 'email': record(eml)})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('operation', choices=['gate', 'propose-rest', 'collect-notify'])
    a = p.parse_args()
    automation = read_json(W/'AUTOMATION.json')
    for pin in automation['files'].values():
        verify(pin)
    if automation['files']['controller'] != record(__file__):
        raise InvalidArtifact('execute pinned finite controller')
    if a.operation == 'gate': print(json.dumps(gate(), indent=2))
    elif a.operation == 'propose-rest': propose_rest()
    else:
        summary, report = collect_and_report(); notify(summary, report)
        print(json.dumps(summary, indent=2))


if __name__ == '__main__': main()
