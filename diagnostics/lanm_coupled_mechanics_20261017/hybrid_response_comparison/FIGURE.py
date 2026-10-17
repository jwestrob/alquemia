#!/usr/bin/env python3
"""Render only a complete saved comparison; editable SVG plus PNG/PDF."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import traceback

PROTOCOL = 'lanm_hybrid_model_and_response_saved_comparison_v1'
SOURCES = ('7CCO', '7CCN')
STAGES = ('Module2_whole', 'hybrid_same_Module2_selection', 'hybrid_parent_pool', 'hybrid_response_pool')
STATES = (('7CCO', 'La'), ('7CCO', 'Dy'), ('7CCN', 'La'), ('7CCN', 'Dy'))
LABELS = ('Legacy whole-LBT\n(Module 2)', 'Hybrid model\nfixed selections',
          'Hybrid model\nreselected parent pool', 'Hybrid model\nresponse pool')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        return {'path': str(p), 'sha256': hashlib.file_digest(f, 'sha256').hexdigest()}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed input: ' + record['path'])
    return Path(record['path'])


def load_complete(path):
    result = json.loads(path.read_text())
    require(result['protocol'] == PROTOCOL and result['status'] == 'complete', 'Complete scientific comparison required')
    config = json.loads(verify(result['config']).read_text())
    verify(result['implementation'])
    verify(config['implementation'])
    require(result['implementation']['sha256'] == config['implementation']['sha256']
            and config['protocol'] == PROTOCOL, 'Comparison source/config identity differs')
    require(result['declared_sources'] == list(SOURCES) and result['declared_stages'] == list(STAGES)
            and result['declared_response_states'] == [list(s) for s in STATES], 'Declared figure denominator differs')
    require([r['source_id'] for r in result['rows']] == list(SOURCES), 'Both ordered LBT sources required')
    values = {}
    for row in result['rows']:
        require(row['status'] == 'available' and [s['stage'] for s in row['stages']] == list(STAGES),
                'Incomplete source/stage sequence')
        ys = [s['R_kcal_mol'] for s in row['stages']]
        require(all(type(y) in (int, float) and math.isfinite(y) for y in ys), 'Missing/nonfinite exchange value')
        values[row['source_id']] = ys
    slots = result['response_terminal_slots']
    require([(s['source_id'], s['metal']) for s in slots] == list(STATES)
            and all(s['status'] == 'admitted' and type(s['converged']) is bool
                    and s['full_force']['converged'] == s['converged'] for s in slots),
            'Four admitted terminal convergence flags required')
    interpretation = result['interpretation']
    require(interpretation['binding_free_energy'] is False and interpretation['qualified_preference_prediction'] is False
            and interpretation['cross_source_pool'] is False and interpretation['cross_region_absolute_energy_subtraction'] is False,
            'Conditional interpretation changed')
    return values, slots


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--result', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Rendering requires Slurm')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tick = time.monotonic()
    receipt = {'status': 'unavailable', 'implementation': pin(__file__), 'input': pin(a.result),
               'slurm_job_id': os.environ['SLURM_JOB_ID'], 'new_molecular_calls': 0, 'artifacts': {}}
    try:
        values, slots = load_complete(a.result)
        # Imported only after complete saved-data admission, inside Slurm.
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                             'svg.fonttype': 'none', 'svg.hashsalt': PROTOCOL,
                             'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
        fig, ax = plt.subplots(figsize=(9.2, 5.5))
        fig.subplots_adjust(left=.115, right=.965, bottom=.31, top=.80)
        colors = ('#2864A5', '#B96A28')
        legends = ('7CCO (La-conditioned)', '7CCN (Lu-conditioned)')
        for index, sid in enumerate(SOURCES):
            line, = ax.plot(range(4), values[sid], marker=('o', 's')[index], markersize=6,
                            linewidth=1.8, color=colors[index], label=legends[index])
            line.set_gid('source_' + sid)
            for x, y in enumerate(values[sid]):
                annotation = ax.annotate(f'{y:.2f}', (x, y), xytext=(0, 10 if index == 0 else -16),
                                         textcoords='offset points', ha='center', fontsize=9, color=colors[index])
                annotation.set_gid(sid + '_' + STAGES[x] + '_value')
        ax.axhline(0, color='#747474', linewidth=.8, linestyle='--', zorder=0)
        ax.grid(axis='y', color='#DDDDDD', linewidth=.6)
        ax.set_axisbelow(True)
        ax.set_xticks(range(4), LABELS)
        ax.set_xlim(-.2, 3.2)
        ax.margins(y=.22)
        ax.set_ylabel('Conditional exchange R (kcal/mol)')
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1.02), ncol=2, frameon=False, borderaxespad=0)
        fig.suptitle('Conditional La/Dy exchange across model and motion stages', x=.115, ha='left', y=.965, fontsize=13)
        n_unconverged = sum(not s['converged'] for s in slots)
        footer = ('R = (Mex Dy − La) − (LBT Dy − La); Mex broad pool stays fixed.\n'
                  'Baseline hybrid region; both LBT sources retained. Signs do not label correctness.\n'
                  f'Final response: {n_unconverged}/4 searches nonconverged. Finite energies, not binding affinities.')
        fig.text(.115, .065, footer, fontsize=9, va='bottom', linespacing=1.5)
        for fmt in ('svg', 'png', 'pdf'):
            target = out / ('conditional_exchange.' + fmt)
            metadata = {'Creator': 'FIGURE.py'}
            if fmt == 'svg':
                metadata.update(Date=None, Description='Input RESULT SHA256: ' + receipt['input']['sha256'])
            elif fmt == 'pdf':
                metadata.update(CreationDate=None, ModDate=None)
            fig.savefig(target, format=fmt, dpi=300, metadata=metadata)
            receipt['artifacts'][fmt] = pin(target)
        plt.close(fig)
        verify(receipt['input'])
        receipt.update(status='complete', matplotlib_version=matplotlib.__version__,
                       plotted_values_kcal_mol=values, stages=list(STAGES),
                       nonconverged_response_searches=n_unconverged,
                       terminal_flags=[{k: s[k] for k in ('source_id', 'metal', 'converged', 'termination')} for s in slots],
                       editable_SVG_text=True, pending_values_used=False)
    except Exception as exc:
        receipt.update(failure=repr(exc), traceback=traceback.format_exc())
    receipt['elapsed_seconds'] = time.monotonic() - tick
    with (out / 'RENDER.json').open('x') as f:
        json.dump(receipt, f, indent=2, allow_nan=False)
        f.write('\n')
    print(json.dumps({'status': receipt['status'], 'receipt': pin(out / 'RENDER.json')}))
    return 0 if receipt['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
