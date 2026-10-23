#!/usr/bin/env python3
"""Plot terminal six-parameter water evidence without new molecular evaluations."""
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ['MPLBACKEND'] = 'Agg'
import matplotlib.pyplot as plt
import numpy as np


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--result', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    result = json.loads(Path(args.result).read_text())
    assert result['status'] == 'bounded_six_parameter_water_complete'
    assert result['error'] is None and result['candidate_frozen']
    assert pin(result['config']['path']) == result['config']
    config = json.loads(Path(result['config']['path']).read_text())
    assert pin(config['training_result']['path']) == config['training_result']
    data = json.loads(Path(config['training_result']['path']).read_text())
    assert data['complete'] == 8
    assert pin(result['prior_result']['path']) == result['prior_result']
    prior = json.loads(Path(result['prior_result']['path']).read_text())
    assert result['initial']['parameters'] == prior['selected']['parameters'] + [0.24429636467571708, 0.37332330078584813]
    for key, row in result['initial']['model'].items():
        assert all(row[k] == prior['selected']['model'][key][k] for k in ('energy', 'forces'))
    models = {
        'PBE0-D4': {k: r['reference'] for k, r in data['rows'].items()},
        'Original parameters': {k: r['original_baseline'] for k, r in data['rows'].items()},
        'Four-parameter selection': result['initial']['model'],
        'Six-parameter selection': result['selected']['model'],
    }
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex='col', layout='constrained')
    styles = [('black', 'o'), ('0.55', 'x'), ('#e69f00', 's'), ('#0072b2', 'D')]
    saved = []
    for col, metal in enumerate(('La', 'Dy')):
        ids = ['water_' + s + '_' + metal for s in ('minus', 'central', 'plus')]
        distances, directions = [], []
        for gid in ids:
            geometry = data['rows'][gid]['geometry']
            pos = np.array(geometry['positions_A'])
            direction = np.array(geometry['radial_direction'])
            assert np.isclose(np.linalg.norm(direction), 1)
            assert np.dot(pos[0] - pos[1], direction) > 0
            distances.append(float(np.linalg.norm(pos[0] - pos[1])))
            directions.append(direction)
        for (name, model), (color, marker) in zip(models.items(), styles):
            center_energy = model[ids[1]]['energy']
            energies = [model[k]['energy'] - center_energy for k in ids]
            derivatives = [float(-np.dot(model[k]['forces'][0], d)) for k, d in zip(ids, directions)]
            axes[0, col].plot(distances, energies, color=color, marker=marker, label=name)
            axes[1, col].plot(distances, derivatives, color=color, marker=marker, label=name)
            saved.append({'metal': metal, 'model': name, 'metal_O_distance_A': distances,
                          'relative_energy_kcal_mol': energies,
                          'outward_dE_dR_kcal_mol_A': derivatives})
        axes[0, col].set_title(metal)
        axes[1, col].set_xlabel('Metal–oxygen distance (Å)')
    for ax in axes.flat:
        ax.axhline(0, color='0.8', lw=.7)
        ax.grid(alpha=.15)
    axes[0, 0].set_ylabel('Energy relative to center (kcal/mol)')
    axes[1, 0].set_ylabel('Outward energy derivative (kcal/mol/Å)')
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Water model repair: three consumed configurations')
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    fig.savefig(out / 'water_response.png', dpi=180)
    fig.savefig(out / 'water_response.pdf')
    plt.close(fig)
    (out / 'DATA.json').write_text(json.dumps({
        'inputs': [pin(args.result), result['config'], config['training_result'], result['prior_result']],
        'implementation': pin(__file__), 'series': saved,
        'development_checks_passed': result['development_checks_passed'],
        'interpretation': 'Three fixed configurations joined for readability; not continuous curves or independent validation. Derivative displaces only the metal outward with water fixed. Aquo8 and full rigid-body controls are reported separately.',
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
