"""Render archived outcome counts only; no fitting or molecular calculations."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

METHODS = ('DFT', 'context_composite')
LABELS = ('Preserved DFT', 'Released context\ncomposite')
OUTCOMES = ('correct', 'wrong', 'inconclusive', 'unavailable')
COLORS = ('#367F82', '#CC694A', '#E7BC60', '#CCD1D9')


def tally(rows, method):
    values = Counter(r['methods'][method]['outcome'] for r in rows)
    if set(values) - set(OUTCOMES): raise ValueError('Unknown archived outcome')
    return {key: values[key] for key in OUTCOMES}


def render(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    data = json.loads(source.read_text()); rows, triples = data['rows'], data['triples']
    if (len(rows), len(triples)) != (225, 100): raise ValueError('Declared population changed')
    groups = Counter(r['root_case_id'] for r in rows)
    triple_groups = Counter(r['root_case_id'] for r in triples)
    if len(groups) != 25 or set(groups.values()) != {9} or set(triple_groups.values()) != {4} or set(groups) != set(triple_groups):
        raise ValueError('The figure requires the same 25 groups and all four triples per group')
    matched = [r for r in rows if all(r['methods'][m]['R'] is not None and
               r['methods'][m]['outcome'] != 'unavailable' for m in METHODS)]
    excluded = [r['case_id'] for r in rows if r not in matched]
    counts = {'individual_common': {m: tally(matched, m) for m in METHODS},
              'strict_triples': {m: tally(triples, m) for m in METHODS}}
    # This is an integrity check against the already published counts, not a new statistic.
    for m in METHODS:
        if counts['strict_triples'][m] != {k: data['counts']['all100_triples'][m][k] for k in OUTCOMES}:
            raise ValueError('Triple count replay differs')
    output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig = plt.figure(figsize=(9.0, 5.35), facecolor='white')
    grid = fig.add_gridspec(2, 2, width_ratios=(3.7, 2.15), left=.215, right=.985,
                           bottom=.21, top=.77, wspace=.07, hspace=.90)
    fig.text(.025, .956, 'PQQ discrimination across structural repeats', fontsize=14, weight='bold')
    fig.text(.025, .91, 'Fixed protocol-specific decision bands • 25 previously used reference proteins', fontsize=9.5, color='#444444')
    legend_labels = ('Correct', 'Wrong', 'Inconclusive', 'Outside coverage / unavailable*')
    fig.legend([Patch(facecolor=c) for c in COLORS], legend_labels, loc='upper left',
               bbox_to_anchor=(.015, .888), ncol=4, frameon=False, fontsize=9,
               handlelength=1.3, columnspacing=1.5)
    panels = [('A', 'Individual folds', len(rows), 'individual_common'),
              ('B', 'Three-fold medians (La-conditioned)', len(triples), 'strict_triples')]
    for row_idx, (letter, title, total, key) in enumerate(panels):
        ax, count_ax = fig.add_subplot(grid[row_idx, 0]), fig.add_subplot(grid[row_idx, 1])
        ax.set_title(f'{letter}   {title}', loc='left', fontsize=10.5, weight='bold', pad=12)
        for i, method in enumerate(METHODS):
            y = 1 - i; left = 0
            c = counts[key][method]
            values = [c[o] for o in OUTCOMES]
            if key == 'individual_common': values[-1] = len(excluded)
            if sum(values) != total: raise ValueError('Bar denominator changed')
            for value, color in zip(values, COLORS):
                ax.barh(y, value, left=left, height=.45, color=color, linewidth=0)
                left += value
            for j, value in enumerate(values):
                count_ax.text((j+.5)/4, (y+.55)/2.1, str(value),
                              transform=count_ax.transAxes, ha='center', va='center', fontsize=10)
        ax.set_yticks([1, 0], LABELS); ax.tick_params(axis='y', length=0, pad=10)
        ax.set_ylim(-.55, 1.55); ax.set_xlim(0, total)
        ax.set_xticks([0, 50, 100, 150, 200, 225] if row_idx == 0 else [0, 25, 50, 75, 100])
        ax.set_xlabel('Declared sources' if row_idx == 0 else 'Declared three-fold subsets', fontsize=8.5)
        ax.spines['left'].set_visible(False); ax.spines['bottom'].set_color('#999999')
        ax.tick_params(axis='x', colors='#555555', labelsize=8)
        count_ax.axis('off')
        for j, label in enumerate(('Correct', 'Wrong', 'Inconcl.', 'Other*')):
            count_ax.text((j+.5)/4, 1.04, label, transform=count_ax.transAxes,
                          ha='center', va='center', fontsize=8, color='#555555')
    fig.text(.025, .088, f'* A: identical {len(matched)}/225-source coverage; {len(excluded)} sources excluded from both matched bars.', fontsize=8.7)
    unavailable = [counts['strict_triples'][m]['unavailable'] for m in METHODS]
    missing_note = f'{unavailable[0]} unavailable for each method' if len(set(unavailable)) == 1 else f'DFT {unavailable[0]} and context {unavailable[1]} unavailable'
    fig.text(.025, .052, '* B: all 100 predeclared subsets retained; each median requires all three members. ' + missing_note + '.', fontsize=8.7)
    fig.text(.025, .018, 'Structural repeats and overlapping subsets are correlated, not independent biological observations. No probability or affinity claim.', fontsize=8.3, color='#444444')
    stem = output / 'structural_comparison'
    fig.savefig(stem.with_suffix('.svg'))
    fig.savefig(stem.with_suffix('.pdf'))
    fig.savefig(stem.with_suffix('.png'), dpi=180)
    plt.close(fig)
    receipt = {'source': {'path': str(source), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()},
               'implementation': {'path': str(Path(__file__).resolve()),
                                  'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
               'matplotlib_version': matplotlib.__version__, 'methods': list(METHODS),
               'groups': 25, 'individual_declared': 225, 'individual_common': len(matched),
               'excluded_source_ids': excluded, 'triple_declared': 100, 'counts': counts,
               'new_molecular_calls': 0, 'new_statistics': False,
               'outputs': [str(stem.with_suffix(ext)) for ext in ('.svg', '.pdf', '.png')]}
    (output / 'figure_data.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True); parser.add_argument('--output', required=True)
    args = parser.parse_args(); result = render(args.source, args.output)
    print(json.dumps({'matched_sources': result['individual_common'], 'counts': result['counts'], 'outputs': result['outputs']}))
