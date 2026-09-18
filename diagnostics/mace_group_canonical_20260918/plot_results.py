"""Plot actual canonical report values; no fitting, rescoring or hidden offset."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report', required=True); p.add_argument('--output', required=True)
    args = p.parse_args(); d = json.loads(Path(args.report).read_text())
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    rows = sorted(d['scores'].items(), key=lambda item: (item[1]['evaluation_role'] != 'calibration',
                                                       item[1]['expected_class'], item[0]))
    fig, axes = plt.subplots(1, 2, figsize=(12, 11), sharey=True)
    colors = {'Ca': '#327ba5', 'La': '#c86b29'}
    for y, (name, r) in enumerate(rows):
        marker = 'o' if r['evaluation_role'] == 'calibration' else 'D'
        for ax, value in zip(axes, (r['baseline']['published_S_kcal_mol'], r['R_model_kcal'])):
            if value is not None:
                ax.scatter(value, y, color=colors[r['expected_class']], marker=marker, s=36, zorder=3)
            else:
                ax.text(.03, y, 'unavailable', transform=ax.get_yaxis_transform(), va='center', color='#777777', fontsize=8)
    axes[0].set_yticks(range(len(rows)), [name.replace('-pqq-la_model', '') for name, r in rows], fontsize=8)
    axes[0].invert_yaxis()
    # Read released baseline bands from the original pinned inventory.
    manifest = json.loads(Path(d['manifest']['path']).read_text())
    prep = json.loads(Path(manifest['preparation']['path']).read_text())
    cfg = json.loads(Path(prep['config']['path']).read_text())
    inv_pin = cfg['inventory']
    if pin(inv_pin['path']) != inv_pin:
        raise ValueError('baseline inventory hash changed')
    inv = json.loads(Path(inv_pin['path']).read_text()); bands = inv['baseline_bands']
    axes[0].axvline(bands['Ca']['S_kcal_mol_max'], color=colors['Ca'], lw=1, ls='--')
    axes[0].axvline(bands['La']['S_kcal_mol_min'], color=colors['La'], lw=1, ls='--')
    c = d['calibration']
    if c['status'] == 'available_for_research':
        axes[1].axvline(c['Ca_max_R_model_kcal'], color=colors['Ca'], lw=1, ls='--')
        axes[1].axvline(c['La_min_R_model_kcal'], color=colors['La'], lw=1, ls='--')
    axes[0].set(title='Preserved fixed-core DFT baseline', xlabel='Published S (kcal/mol)')
    axes[1].set(title='Unchanged full-protein typed-group MACE', xlabel='Raw R = E(Ca) − E(La) (model kcal)')
    axes[1].ticklabel_format(axis='x', style='plain', useOffset=False)
    for ax in axes:
        ax.axhline(24.5, color='#aaaaaa', lw=.8)
        ax.grid(axis='x', alpha=.18); ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
    gap = c['observed_class_gap_model_kcal']
    subtitle = 'Calibration pending' if gap is None else f'Calibration gap: {gap:.3f} model kcal; {c["status"].replace("_", " ")}'
    fig.suptitle('PQQ class calibration and retrospective crystal transfer', fontsize=15, y=.99)
    fig.text(.5, .96, subtitle, ha='center', fontsize=10)
    legend = [Line2D([], [], marker='o', ls='', color=colors[label], label=label + '-class evidence') for label in ('Ca', 'La')]
    legend += [Line2D([], [], marker='D', ls='', color='black', label='Retrospective crystal transfer')]
    fig.legend(handles=legend, loc='lower center', bbox_to_anchor=(.5, .062), ncol=3, frameon=False)
    fig.text(.025, .043, '25 consumed calibration cases; crystals 1H4I/4MAE repeat calibration sequences. 1KB0 full-protein source is unsupported.', fontsize=9)
    fig.text(.025, .022, 'Different protocols and score scales. Dashed lines are each protocol’s own bands. Prior MACE grouping failures remain.', fontsize=9)
    fig.tight_layout(rect=(0, .10, 1, .94))
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out / ('comparison.' + ext), dpi=170)
    src = out / 'plot_results.py'; src.write_bytes(Path(__file__).read_bytes())
    provenance = {'report': pin(args.report), 'implementation': pin(src), 'baseline_inventory': inv_pin,
                  'matplotlib_version': matplotlib.__version__, 'display_offset': None,
                  'score_formula_changed': False, 'molecular_evaluations': 0}
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    plt.close(fig)


if __name__ == '__main__':
    main()
