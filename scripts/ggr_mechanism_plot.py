#!/usr/bin/env python3
"""Plot approved GGR comparisons from pinned real collections; never run science.

Missing results stay unavailable. Existing output directories are immutable.
All plotted values and source records are exported alongside PNG/PDF figures.
"""
from __future__ import annotations

import argparse
import csv
import math
import shutil
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from affordable_common import HA_TO_KCAL, InvalidArtifact, energy, read_json, record, verify, write_new
from ggr_workflow import checked_rows

V2 = 'generic_vertical_exchange_native_r2scan3c_v2'
V3 = 'generic_peptide_amide_vertical_native_r2scan3c_v3'
NMA = 'generic_peptide_alpha_caps_native_r2scan3c_dev_v1'
CONNECTED = 'ggr_connected_segment_native_r2scan3c_dev_v1'
COLORS = {'formamide': '#4477AA', 'nma': '#EE7733', 'analytic': '#4477AA',
          'pass': '#228833', 'failed': '#CC3311'}


def _source(row):
    if row.get('source_structure_id'):
        return row['source_structure_id'].upper()
    # These are the explicit case IDs in the two historical collections.
    return {'ggr_1glg_GGR': '1GLG', 'alacta_1f6s_v2_strong_site': '1F6S',
            'alacta_6ip9_v2_strong_site': '6IP9'}.get(row['case'])


def _select(rows, source, protocol):
    matches = [r for r in rows if _source(r) == source and r['protocol_id'] == protocol]
    if len(matches) > 1:
        raise InvalidArtifact(f'ambiguous source/protocol result: {source}, {protocol}')
    return matches[0] if matches else None


def _point(row):
    if row is None:
        return {'status': 'not_supplied', 'value': None}
    result = {'status': row['status'], 'case': row['case'], 'protocol_id': row['protocol_id'],
              'collection': row['_collection'], 'value': None}
    if row['status'] == 'complete':
        value = row['score'].get('S_kcal_mol')
        if value is None:
            result['status'] = 'reporting_reference_unavailable'
        elif not math.isfinite(value):
            raise InvalidArtifact('nonfinite reported score')
        else:
            result['value'] = value
            result['R_hartree'] = row['score']['R_hartree']
    else:
        result['reasons'] = {m: e.get('reason', e['status']) for m, e in row.get('endpoints', {}).items()}
    return result


def _difference(alpha, ggr):
    if alpha is None or ggr is None:
        return {'status': 'not_supplied', 'value': None}
    result = {'alpha_case': alpha['case'], 'GGR_case': ggr['case'], 'status': 'unavailable', 'value': None,
              'alpha_collection': alpha['_collection'], 'GGR_collection': ggr['_collection']}
    if alpha['protocol_id'] != ggr['protocol_id']:
        raise InvalidArtifact('relative ordering requires a common preparation protocol')
    if alpha['status'] == ggr['status'] == 'complete':
        result.update(status='complete', value=(alpha['score']['R_hartree']-ggr['score']['R_hartree'])*HA_TO_KCAL)
    return result


def _save(fig, output, stem, title, records):
    fig.savefig(output / f'{stem}.png', dpi=220, bbox_inches='tight')
    fig.savefig(output / f'{stem}.pdf', bbox_inches='tight', metadata={'Title': title})
    plt.close(fig)
    records[stem] = {suffix: record(output / f'{stem}.{suffix}') for suffix in ('png', 'pdf')}


def _style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', color='#DDDDDD', linewidth=.6)
    ax.set_axisbelow(True)


def ladder(rows, output, records):
    labels = ['Historical v2\nformaldehyde', 'Repaired v3\nformamide', 'Extended amide\nNMA', 'Connected\nnative segment']
    points = [_point(_select(rows, '1GLG', p)) for p in (V2, V3, NMA, CONNECTED)]
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for i, point in enumerate(points):
        if point['value'] is None:
            ax.text(i, .5, point['status'].replace('_', '\n'), transform=ax.get_xaxis_transform(),
                    ha='center', color='#777777', fontsize=9)
        else:
            ax.scatter(i, point['value'], s=65, color='#777777' if i == 0 else '#4477AA', zorder=3)
            ax.annotate(f'{point["value"]:.2f}', (i, point['value']), xytext=(0, 10),
                        textcoords='offset points', ha='center')
    ax.set_xticks(range(4), labels)
    ax.set_xlim(-.4, 3.4)
    ax.margins(y=.2)
    ax.set_ylabel('S in existing aquo reporting gauge (kcal/mol)')
    ax.set_title('GGR 1GLG: representation sensitivity')
    _style(ax)
    fig.text(.5, .015, 'Larger S is more La-like; no calibrated cutoff for these generic models.\n'
             'Frozen source geometry and microstate; representation changes local interactions and CPCM cavity.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .085, 1, 1))
    _save(fig, output, 'ggr_representation_ladder', 'GGR representation sensitivity', records)
    return [{'label': label.replace('\n', ' '), **point} for label, point in zip(labels, points)]


def source_comparison(rows, output, records):
    sources = ['1GLG', '2FW0', '2FVY']
    points = {name: [_point(_select(rows, source, protocol)) for source in sources]
              for name, protocol in [('formamide', V3), ('nma', NMA)]}
    fig, (ax, table_ax) = plt.subplots(2, 1, figsize=(8.2, 6), gridspec_kw={'height_ratios': [3, 1.3]})
    for name, offset, marker in [('formamide', -.07, 'o'), ('nma', .07, 'D')]:
        valid = [(i, p['value']) for i, p in enumerate(points[name]) if p['value'] is not None]
        if valid:
            ax.plot([i+offset for i, _ in valid], [v for _, v in valid], marker=marker,
                    color=COLORS[name], label='Repaired formamide' if name == 'formamide' else 'Extended amide (NMA)',
                    linestyle='none', markersize=7)
    ax.set_xticks(range(3), ['1GLG\ngalactose-bound', '2FW0\nsugar-free/open', '2FVY\nglucose-bound/closed'])
    ax.set_xlim(-.4, 2.4)
    ax.set_ylabel('S, existing reporting gauge (kcal/mol)')
    ax.set_title('GGR source structures under two representations')
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, fontsize=9)
    _style(ax)
    table_ax.axis('off')
    cells = [[source] + [f'{points[name][i]["value"]:.3f}' if points[name][i]['value'] is not None
                        else points[name][i]['status'].replace('_', ' ')
                        for name in ('formamide', 'nma')] for i, source in enumerate(sources)]
    table = table_ax.table(cellText=cells, colLabels=['Source', 'Formamide S', 'NMA S'], loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)
    fig.text(.5, .012, 'One biological observation; different crystal conditions. This is structural robustness,\n'
             'not a causal sugar-only effect. Missing results have no plotted value.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .055, 1, 1))
    _save(fig, output, 'ggr_source_representation', 'GGR source and representation comparison', records)
    return [{'source': source, 'representation': name, **points[name][i]}
            for i, source in enumerate(sources) for name in ('formamide', 'nma')]


def alpha_comparison(rows, output, records, comparison):
    points = []
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for name, protocol, offset, marker in [('formamide', V3, -.08, 'o'), ('nma', NMA, .08, 'D')]:
        ggr = _select(rows, '1GLG', protocol)
        for i, source in enumerate(('1F6S', '6IP9')):
            point = {'source': source, 'representation': name,
                     **_difference(_select(rows, source, protocol), ggr)}
            if comparison and name == 'nma' and point['value'] is not None:
                matches = [p for p in comparison['alpha_minus_GGR'] if p['case'] == point['alpha_case']]
                if len(matches) != 1 or not math.isclose(matches[0]['alpha_minus_GGR_R_kcal_mol'], point['value'], abs_tol=1e-8):
                    raise InvalidArtifact('supplied comparison disagrees with plotted matched contrast')
            points.append(point)
            if point['value'] is None:
                ax.text(i+offset, .5, point['status'].replace('_', '\n'), transform=ax.get_xaxis_transform(),
                        ha='center', fontsize=8, color='#777777')
            else:
                ax.scatter(i+offset, point['value'], s=65, marker=marker, color=COLORS[name], zorder=3)
                ax.annotate(f'{point["value"]:.2f}', (i+offset, point['value']),
                            xytext=(-18 if name == 'formamide' else 18, -16 if name == 'formamide' else 10),
                            textcoords='offset points', ha='center', fontsize=9)
    ax.axhline(0, color='#888888', linewidth=.8)
    ax.set_xticks([0, 1], ['Alpha 1F6S\nCa-derived geometry', 'Alpha 6IP9\nLa-derived geometry'])
    ax.set_xlim(-.5, 1.5)
    ax.margins(y=.15)
    ax.set_ylabel('R(alpha) − R(GGR 1GLG), kcal/mol')
    ax.set_title('Matched alpha-lactalbumin versus GGR contrasts')
    ax.legend(handles=[Line2D([], [], marker='o', linestyle='', color=COLORS['formamide'], label='Repaired formamide'),
                       Line2D([], [], marker='D', linestyle='', color=COLORS['nma'], label='Extended amide (NMA)')],
              frameon=False, fontsize=9)
    _style(ax)
    fig.text(.5, .015, 'Negative: alpha is less La-like than GGR on that protocol. The common aquo offset cancels.\n'
             'Two alpha structures represent one condition-qualified biological observation.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .085, 1, 1))
    _save(fig, output, 'alpha_minus_ggr', 'Matched alpha minus GGR contrasts', records)
    return points


def _read_sensitivity(path):
    data = read_json(path)
    if data.get('schema_version') != 'alquemia.ggr_sensitivity_comparison.v1':
        raise InvalidArtifact('unsupported sensitivity collection schema')
    verify(data['manifest'])
    for endpoint in data['energies'].values():
        if endpoint['status'] == 'complete':
            verify(endpoint['receipt'])
            if energy(verify(endpoint['output'])) != endpoint['energy_hartree']:
                raise InvalidArtifact('sensitivity endpoint energy differs from recorded real output')
    return data


def gradient_checks(paths, output, records):
    blocks = {}
    sources = []
    for path in paths:
        data = _read_sensitivity(path)
        sources.append({'collection': record(path), 'normal_tight_bridge': data['normal_tight_bridge']})
        for row in data['comparisons']:
            key = row['representation'], row['coordinate'], row['phase']
            if key in blocks:
                raise InvalidArtifact(f'duplicate sensitivity block: {key}')
            blocks[key] = row
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.8))
    values = []
    for ri, representation in enumerate(('extended', 'connected')):
        for ci, coordinate in enumerate(('metal', 'peptide')):
            ax = axes[ri, ci]
            nominal = blocks.get((representation, coordinate, 'nominal'))
            half = blocks.get((representation, coordinate, 'half'))
            active = [r for r in (nominal, half) if r is not None]
            plotted, labels = 0, []
            for row in active:
                if 'checks' not in row:
                    values.append({k: row.get(k) for k in ('representation', 'coordinate', 'phase', 'status',
                                                          'reason', 'units', 'gradient_units', 'amplitude_units')})
                    continue
                for metal in ('La', 'Ca', 'R'):
                    check = row['checks'][metal]
                    pred, odd, tol = (check[k] for k in ('predicted_odd_energy_kcal_mol', 'odd_energy_kcal_mol', 'tolerance_kcal_mol'))
                    if not all(math.isfinite(v) for v in (pred, odd, tol)) or tol <= 0:
                        raise InvalidArtifact('invalid numerical consistency result')
                    label = f'{metal} ({"h" if row["phase"] == "nominal" else "h/2"})'
                    ax.errorbar(pred, plotted, xerr=tol, color=COLORS['analytic'], fmt='D', markersize=5,
                                capsize=3, elinewidth=1.2, zorder=2)
                    ax.plot([pred, odd], [plotted, plotted], color='#BBBBBB', linewidth=1)
                    ax.scatter(odd, plotted, color=COLORS['pass'] if check['pass'] else COLORS['failed'],
                               marker='o' if check['pass'] else 'x', s=38, zorder=3)
                    values.append({'representation': representation, 'coordinate': coordinate, 'phase': row['phase'],
                                   'endpoint_or_contrast': metal, 'status': row['status'], 'units': row['units'],
                                   'gradient_units': row.get('gradient_units', row['units']),
                                   'amplitude_units': row.get('amplitude_units',
                                                             {'metal': 'A', 'peptide': 'radian'}[coordinate]),
                                   **check})
                    labels.append(label)
                    plotted += 1
            if not plotted:
                status = '; '.join(f'{r["phase"]}: {r["status"]}' for r in active) or 'Sensitivity collection not supplied'
                ax.text(.5, .5, status, transform=ax.transAxes, ha='center', va='center', wrap=True, fontsize=10)
                ax.set_yticks([])
                ax.set_xticks([])
            else:
                ax.set_yticks(range(plotted), labels)
                ax.set_ylim(plotted-.5, -.5)
            if nominal and not half and nominal.get('half_step_eligible'):
                ax.text(.5, .01, 'Half-step triggered; result unavailable', transform=ax.transAxes,
                        ha='center', color=COLORS['failed'], fontsize=8)
            title = ('Extended amide' if representation == 'extended' else 'Connected peptide')
            title += ': ' + ('metal motion' if coordinate == 'metal' else 'peptide crankshaft')
            ax.set_title(title, fontsize=11)
            ax.set_xlabel('Odd energy change (kcal/mol)')
            _style(ax)
    fig.legend(handles=[Line2D([], [], marker='D', color=COLORS['analytic'], linestyle='', label='Analytic h × gradient'),
                        Line2D([], [], marker='o', color=COLORS['pass'], linestyle='', label='Actual odd energy: pass'),
                        Line2D([], [], marker='x', color=COLORS['failed'], linestyle='', label='Actual odd energy: failed check')],
               loc='upper center', ncol=3, frameon=False, fontsize=9)
    fig.text(.5, .008, 'Intervals show predeclared acceptance tolerances, not statistical uncertainty. h = 0.02 Å or 1°; optional h/2 shown separately.\n'
             'This checks physical derivatives of the frozen CPCM model; no relaxation or entropy correction is validated.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .065, 1, .95))
    _save(fig, output, 'physical_gradient_checks', 'Physical gradient and directional energy checks', records)
    return {'checks': values, 'sources_and_normal_tight_bridges': sources}


def plot(collections, historical, comparison_path, sensitivity, output):
    output = Path(output)
    if output.exists():
        raise InvalidArtifact(f'refusing to overwrite plot directory: {output}')
    if len({str(p.resolve()) for p in collections+historical}) != len(collections+historical):
        raise InvalidArtifact('duplicate collection argument')
    rows, unsupported = [], []
    for path in collections+historical:
        data = checked_rows(path)
        rows.extend({**r, '_collection': record(path)} for r in data['rows'])
        unsupported.extend(read_json(verify(data['manifest'])).get('unsupported_preparations', []))
    comparison = read_json(comparison_path) if comparison_path is not None else None
    if comparison is not None:
        for rec in comparison['collections']+comparison['historical_collections']:
            verify(rec)
        if {record(p)['sha256'] for p in collections} != {r['sha256'] for r in comparison['collections']}:
            raise InvalidArtifact('comparison must describe the supplied current collections exactly')
    output.mkdir(parents=True)
    preserved_script = output / Path(__file__).name
    shutil.copyfile(__file__, preserved_script)
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    artifacts = {}
    result = {'schema_version': 'alquemia.ggr_mechanism_plot.v1', 'implementation': record(preserved_script),
              'historical_live_implementation_at_render': record(__file__),
              'collections': [record(p) for p in collections], 'historical_collections': [record(p) for p in historical],
              'comparison': record(comparison_path) if comparison_path else None,
              'sensitivity_collections': [record(p) for p in sensitivity],
              'scientific_evaluations_executed': 0,
              'unsupported_preparations': unsupported,
              'gradient_units': 'units and gradient_units describe projected_gradient (kcal/mol/A or kcal/mol/radian); '
                                'amplitude_units separately describes amplitude (A or radian)',
              'interpretation': 'All cases are method development; no decision cutoff, fitted model or new scientific comparison.'}
    result['representation_ladder'] = ladder(rows, output, artifacts)
    result['source_by_representation'] = source_comparison(rows, output, artifacts)
    result['alpha_minus_GGR'] = alpha_comparison(rows, output, artifacts, comparison)
    result['physical_gradient_checks'] = gradient_checks(sensitivity, output, artifacts)
    result['figures'] = artifacts
    write_new(output / 'plot_manifest.json', result)
    table_rows = [{'figure': 'representation_ladder', **r} for r in result['representation_ladder']]
    table_rows += [{'figure': 'source_by_representation', **r} for r in result['source_by_representation']]
    table_rows += [{'figure': 'alpha_minus_GGR', **r} for r in result['alpha_minus_GGR']]
    with (output / 'score_plot_values.tsv').open('x', newline='') as stream:
        fields = ['figure', 'label', 'source', 'representation', 'case', 'alpha_case', 'GGR_case', 'status', 'value']
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(table_rows)
    with (output / 'physical_gradient_values.tsv').open('x', newline='') as stream:
        fields = ['representation', 'coordinate', 'phase', 'endpoint_or_contrast', 'status', 'units',
                  'gradient_units', 'amplitude_units', 'amplitude',
                  'projected_gradient', 'predicted_odd_energy_kcal_mol', 'odd_energy_kcal_mol', 'residual_kcal_mol',
                  'residual_per_coordinate_unit', 'tolerance_kcal_mol', 'pass', 'reason']
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(result['physical_gradient_checks']['checks'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--collection', action='append', type=Path, required=True)
    parser.add_argument('--historical-collection', action='append', type=Path, default=[])
    parser.add_argument('--comparison', type=Path)
    parser.add_argument('--sensitivity', action='append', type=Path, default=[])
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    result = plot(args.collection, args.historical_collection, args.comparison, args.sensitivity, args.output_dir)
    print(f'{len(result["figures"])} PNG/PDF figures exported; zero scientific evaluations executed')


if __name__ == '__main__':
    main()
