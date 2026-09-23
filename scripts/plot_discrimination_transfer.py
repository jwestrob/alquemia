"""Plot an existing frozen-reference transfer ledger; no fitting or molecular work."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from affordable_common import InvalidArtifact, read_json, record, write_new


def render(comparison, methods, output, xlim=None):
    ledger = read_json(comparison)
    labels = dict(item.split('=', 1) for item in methods)
    if len(labels) != len(methods) or not labels:
        raise InvalidArtifact('Provide distinct method=title arguments')
    rows = ledger['rows']
    if len({r['case_id'] for r in rows}) != len(rows):
        raise InvalidArtifact('Duplicate source rows')
    groups = {}
    for row in rows:
        group = row['root_case_id']
        if group in groups and groups[group] != row['expected_class']:
            raise InvalidArtifact('Contradictory protein labels')
        if row['expected_class'] not in ('Ca', 'La'):
            raise InvalidArtifact('Unsupported reference label')
        groups[group] = row['expected_class']
    order = sorted(groups, key=lambda k: (groups[k], k))
    ypos = {k: i for i, k in enumerate(order)}
    plotted, counts, bands = [], {}, {}
    for method in labels:
        band = ledger['bands'][method]
        lo, hi = band['Ca_max'], band['La_min']
        if not np.isfinite([lo, hi]).all() or lo >= hi:
            raise InvalidArtifact('Separated frozen decision bands required')
        center = (lo + hi) / 2
        bands[method] = {'Ca_max': lo, 'La_min': hi, 'display_center': center}
        counts[method] = dict(Counter(r['methods'][method]['outcome'] for r in rows))
        for r in rows:
            value = r['methods'][method]
            score = value['R']
            actual = ('unavailable' if score is None else
                      'correct' if (score <= lo and r['expected_class'] == 'Ca') or
                                   (score >= hi and r['expected_class'] == 'La') else
                      'inconclusive' if lo < score < hi else 'wrong')
            if value['outcome'] != actual:
                raise InvalidArtifact('Outcome does not match frozen bands: ' + r['case_id'])
            plotted.append({'method': method, 'case_id': r['case_id'],
                'protein': r['root_case_id'], 'expected_class': r['expected_class'],
                'source_conditioning_metal': r['source_conditioning_metal'],
                'R_model_kcal_mol': score,
                'R_minus_band_midpoint_model_kcal_mol': None if score is None else score-center,
                'decision': value['decision'], 'outcome': actual})
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    with (out/'source_values.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(plotted[0]))
        writer.writeheader()
        writer.writerows(plotted)
    plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42,
                         'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    single = len(labels) == 1
    fig, axes = plt.subplots(1, len(labels), figsize=(7.2, 11.5) if single else (4.1*len(labels)+1.5, 10.7),
                             sharey=True, squeeze=False)
    colors = {'Ca': '#2878ad', 'La': '#8e4b9c'}
    allx = [r['R_minus_band_midpoint_model_kcal_mol'] for r in plotted
            if r['R_model_kcal_mol'] is not None]
    xmin, xmax = min(allx+[0]), max(allx+[0])
    pad = max(2, (xmax-xmin)*.06)
    limits = [xmin-pad, xmax+pad] if xlim is None else list(xlim)
    if len(limits) != 2 or not np.isfinite(limits).all() or limits[0] >= limits[1]:
        raise InvalidArtifact('Two finite ascending display limits required')
    if xmin < limits[0] or xmax > limits[1]:
        raise InvalidArtifact('Display limits would hide a plotted value')
    for ax, method in zip(axes[0], labels):
        b = bands[method]
        lower, upper = b['Ca_max']-b['display_center'], b['La_min']-b['display_center']
        ax.axvspan(lower, upper, color='#eeeeee', zorder=0)
        ax.axvline(lower, color='#999999', lw=.7)
        ax.axvline(upper, color='#999999', lw=.7)
        for group in order:
            local = [r for r in plotted if r['method'] == method and r['protein'] == group]
            valid = [r for r in local if r['R_model_kcal_mol'] is not None]
            xs = [r['R_minus_band_midpoint_model_kcal_mol'] for r in valid]
            if xs:
                ax.plot([min(xs), max(xs)], [ypos[group]]*2, color=colors[groups[group]],
                        lw=.65, alpha=.35, zorder=1)
            for r in valid:
                x = r['R_minus_band_midpoint_model_kcal_mol']
                y = ypos[group]+(-.1 if r['source_conditioning_metal'] == 'Ca' else .1)
                ax.scatter(x, y, marker='o' if r['source_conditioning_metal'] == 'Ca' else '^',
                           s=22, color=colors[groups[group]], alpha=.75, linewidths=.25, zorder=2)
                if r['outcome'] == 'wrong':
                    ax.scatter(x, y, marker='x', s=60, color='#c22c27', lw=1.2, zorder=4)
                elif r['outcome'] == 'inconclusive':
                    ax.scatter(x, y, marker='o', s=65, facecolors='none',
                               edgecolors='#dd9300', linewidths=1.1, zorder=3)
            balanced = [p for p in ledger['pools']
                        if p['root_case_id'] == group and p['pool'] == 'balanced']
            if len(balanced) != 1:
                raise InvalidArtifact('One existing strict balanced aggregate per protein required')
            aggregate = balanced[0]['methods'][method]['R']
            if aggregate is not None:
                ax.scatter(aggregate-b['display_center'], ypos[group], marker='D',
                           s=14, facecolors='none', edgecolors='black', linewidths=.9, zorder=5)
            ax.text(1.02, ypos[group], f'{len(valid)}/{len(local)}',
                    transform=ax.get_yaxis_transform(), color='#666666', fontsize=7, va='center')
        c = counts[method]
        ax.set_title(labels[method]+'\n'+
            f"{c.get('correct',0)} correct · {c.get('wrong',0)} wrong\n"
            f"{c.get('inconclusive',0)} inconclusive · {c.get('unavailable',0)} unavailable", fontsize=10)
        ax.set_xlim(*limits)
        ax.set_xlabel('Contrast minus frozen band midpoint\n(model kcal/mol)')
        ax.grid(axis='y', alpha=.15, lw=.5)
        ax.set_yticks(range(len(order)), [x.removesuffix('-pqq-la_model').upper() for x in order])
    axes[0,0].set_ylim(len(order)-.35, -.65)
    for tick, group in zip(axes[0,0].get_yticklabels(), order):
        tick.set_color(colors[groups[group]])
    legend = [Line2D([], [], color=colors[z], marker='s', ls='', label=f'Known {z} class') for z in ('Ca','La')]
    legend += [Line2D([], [], color='#555555', marker=m, ls='', label=t) for m,t in
               [('o','Ca-conditioned source'),('^','La-conditioned source')]]
    legend += [Line2D([], [], color='black', marker='D', fillstyle='none', ls='', label='Strict balanced aggregate'),
               Line2D([], [], color='#c22c27', marker='x', ls='', label='Wrong'),
               Line2D([], [], color='#dd9300', marker='o', fillstyle='none', ls='', label='Inconclusive')]
    fig.legend(handles=legend, loc='lower center', ncol=2 if single else 4,
               bbox_to_anchor=(.5,.07 if single else .035), fontsize=8)
    fig.suptitle('Nikasha: source-structure discrimination', fontsize=13 if single else 14, y=.99)
    caption = (f'{len(rows)} consumed source structures in {len(groups)} protein groups; repeats are not independent biology.\n'
               'Grey: each method’s own inconclusive band. Right labels: available/declared sources. '
               'Only a display offset is subtracted; no rescaling or refitting.')
    if single:
        caption = (f'{len(rows)} consumed source structures in {len(groups)} protein groups.\n'
                   'Repeats are not independent biology. Grey: inconclusive band.\n'
                   'Right labels: available/declared sources.\n'
                   'Only a display offset is subtracted; no rescaling or refitting.')
    fig.text(.5,.014 if single else .012, caption, ha='center', fontsize=8)
    fig.tight_layout(rect=(0,.185 if single else .105,.99,.97), w_pad=2.5)
    for ext in ('svg','pdf','png'):
        fig.savefig(out/f'structural_transfer.{ext}', dpi=180)
    plt.close(fig)
    write_new(out/'receipt.json', {'input': record(comparison), 'implementation': record(__file__),
        'methods': labels, 'frozen_bands': bands, 'counts': counts,
        'declared_sources': len(rows), 'protein_groups': len(groups),
        'display_xlim': limits,
        'new_molecular_calls': 0, 'new_thresholds': False,
        'aggregate_policy': 'Existing complete La4/Ca5 medians, then equal mean; no new aggregation.',
        'outputs': [record(p) for p in sorted(out.iterdir())]})
    return counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', required=True)
    parser.add_argument('--methods', required=True, nargs='+', help='Existing ledger key=display title')
    parser.add_argument('--output', required=True)
    parser.add_argument('--xlim', nargs=2, type=float, help='Display limits only; may not hide plotted values')
    print(json.dumps(render(**vars(parser.parse_args())), indent=2))
