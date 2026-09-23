"""Export an existing common-pool energy response and source-mapped distances.

No molecular evaluations, fitting, coordination assignment or population model.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from accommodation_nonlinear import relative_components


def render(collection, case, output, display_radius_A=3.5):
    data = read_json(collection)
    matches = [r for r in data['cases'] if r['case_id'] == case]
    if len(matches) != 1 or matches[0]['pool']['status'] != 'available':
        raise InvalidArtifact('One completed common-pool source is required')
    row = matches[0]
    candidates = ['origin', 'adaptive_Ca', 'adaptive_La']
    if [c['id'] for c in row['candidates']] != candidates:
        raise InvalidArtifact('Expected the actual common three-candidate pool')
    coordinates, pins, receipts, tasks = {}, {}, {}, {}
    for c in row['candidates']:
        pins[c['id']] = c['xyz']
        coordinates[c['id']] = xyz(verify(c['xyz']))
    for metal in ('Ca', 'La'):
        pin = row['aliases']['adaptive_' + metal]['proposal_receipt']
        receipt = read_json(verify(pin))
        manifest = read_json(verify(receipt['manifest']))
        task = next(t for t in manifest['tasks'] if t['task_id'] == receipt['task_id'])
        if receipt['status'] != 'proposal_available' or not receipt['final_geometry']['physical_feasible']:
            raise InvalidArtifact('Admitted physical proposal required')
        receipts[metal], tasks[metal] = receipt, task
    if tasks['Ca']['source_preparation'] != tasks['La']['source_preparation']:
        raise InvalidArtifact('Paired source mapping differs')
    preparation_pin = tasks['Ca']['source_preparation']
    prep = read_json(verify(preparation_pin))
    source_rows = prep['mapping']['source_to_qm']
    arrays = {k: np.array([a[1:] for a in v], dtype=float) for k, v in coordinates.items()}
    origin = coordinates['origin']
    if origin[0][0] not in ('Ca', 'La'):
        raise InvalidArtifact('Declared metal must occupy index zero')
    for name, atoms in coordinates.items():
        if len(atoms) != len(origin) or [a[0] for a in atoms[1:]] != [a[0] for a in origin[1:]]:
            raise InvalidArtifact('Physical composition differs between candidates')
        if not np.array_equal(arrays[name][0], arrays['origin'][0]):
            raise InvalidArtifact('This angular illustration requires a fixed metal')
    distances = []
    for atom in source_rows:
        source = atom.get('source')
        if not source or source['element'] not in ('N', 'O'):
            continue
        i = atom['qm_index']
        values = {k: float(np.linalg.norm(a[i] - a[0])) for k, a in arrays.items()}
        if values['origin'] > display_radius_A:
            continue
        label = f"{source['resname']}{source['resnum']}{source['insertion_code']} {source['atom']}"
        distances.append({'qm_index': i, 'chain': source['chain'], 'atom_label': label, **values})
    if not distances:
        raise InvalidArtifact('No mapped heteroatoms within the declared display radius')
    work_rows = []
    for metal in ('Ca', 'La'):
        matrix = row['matrix'][metal]
        for candidate in candidates:
            if matrix[candidate]['status'] != 'complete':
                raise InvalidArtifact('Every cross-scored matrix cell must be complete')
            work = relative_components(matrix[candidate]['components'], matrix['origin']['components'])
            old = row['pool']['rows'][metal]['work_from_origin_kcal_mol'][candidate]
            if any(abs(work[k] - old[k]) > 1e-8 for k in work):
                raise InvalidArtifact('Raw component algebra differs from archived work')
            work_rows.append({'metal': metal, 'candidate': candidate,
                'selected': candidate == row['pool']['rows'][metal]['operational_candidate'], **work})
    selected = {r['metal']: r for r in work_rows if r['selected']}
    delta = selected['Ca']['composite_kcal_mol'] - selected['La']['composite_kcal_mol']
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    for name, rows in [('distances_A.csv', distances), ('relative_energies_kcal_mol.csv', work_rows)]:
        with (out / name).open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.8), gridspec_kw={'width_ratios': [1, 1.65]})
    colors = {'origin': '#777777', 'adaptive_Ca': '#2878ad', 'adaptive_La': '#8e4b9c'}
    labels = {'origin': 'Origin', 'adaptive_Ca': 'Ca proposal', 'adaptive_La': 'La proposal'}
    for i, metal in enumerate(('Ca', 'La')):
        for j, candidate in enumerate(candidates):
            value = next(r for r in work_rows if r['metal'] == metal and r['candidate'] == candidate)
            x = i + (j-1)*.23
            axes[0].bar(x, value['composite_kcal_mol'], width=.2, color=colors[candidate],
                        label=labels[candidate] if i == 0 else None,
                        edgecolor='black' if value['selected'] else 'none', linewidth=1.4)
            axes[0].text(x, value['composite_kcal_mol']-.3, f"{value['composite_kcal_mol']:.2f}",
                         ha='center', va='top', fontsize=9)
    axes[0].axhline(0, color='#555555', linewidth=.7)
    axes[0].set_xticks([0, 1], ['Ca endpoint', 'La endpoint'])
    axes[0].set_ylim(min(r['composite_kcal_mol'] for r in work_rows)-2, 1.7)
    axes[0].set_ylabel('Composite energy change from own origin\n(model kcal/mol)')
    axes[0].set_title('A  Both metals compete over the same geometries', loc='left', fontsize=10)
    x = np.arange(len(distances))
    for j, candidate in enumerate(candidates):
        axes[1].scatter(x+(j-1)*.14, [r[candidate] for r in distances],
                        color=colors[candidate], s=35, label=labels[candidate], zorder=3)
    for i, d in enumerate(distances):
        axes[1].plot([i]*3, [d[c] for c in candidates], color='#bbbbbb', linewidth=.6)
    axes[1].set_xticks(x, [r['atom_label'].replace(' ', '\n') for r in distances], fontsize=8)
    axes[1].set_ylabel('Metal–heteroatom distance (Å)')
    axes[1].set_title('B  Physical response at the actual source atoms', loc='left', fontsize=10)
    axes[1].grid(axis='y', alpha=.2)
    handles, legend = axes[1].get_legend_handles_labels()
    fig.legend(handles, legend, loc='lower center', ncol=3, bbox_to_anchor=(.5, .17))
    fig.suptitle(case.replace('-pqq-la_model__conditioned_', ': ')
                 .replace('__seed-1_sample-', '-conditioned fold '), fontsize=13)
    fig.text(.5, .105, f"Selected Ca work − selected La work = {delta:+.3f} model kcal/mol; "
             'positive shifts the contrast toward La.', ha='center', fontsize=10)
    flagged = ', '.join(z for z in ('Ca','La') if receipts[z]['boundary_flag']) or 'none'
    fig.text(.5, .035, f"Consumed development example; boundary-limited proposal: {flagged}. "
             'Black outlines mark selected energies.\n'
             f'All source-mapped N/O atoms within {display_radius_A:g} Å at origin are shown; '
             'this is a display rule, not a coordination assignment.\n'
             'Fixed composition and pocket membership; finite energy-selected candidates, '
             'not equilibrium populations or binding free energies.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0,.23,1,.94), w_pad=2.3)
    for ext in ('svg', 'pdf', 'png'):
        fig.savefig(out / ('accommodation_example.' + ext), dpi=180)
    plt.close(fig)
    write_new(out/'receipt.json', {'input': record(collection), 'case_id': case,
        'implementation': record(__file__), 'preparation': preparation_pin, 'coordinate_pins': pins,
        'display_radius_A': display_radius_A, 'source_atom_count': len(distances),
        'active_mode_ids': tasks['Ca']['active_mode_ids'],
        'proposals': {z: row['aliases']['adaptive_'+z]['proposal_receipt'] for z in ('Ca','La')},
        'selected_work_difference_model_kcal_mol': delta, 'new_molecular_calls': 0,
        'new_thresholds': False, 'outputs': [record(p) for p in sorted(out.iterdir())]})
    return {'output': str(out), 'displayed_source_atoms': len(distances), 'delta_R': delta}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--collection', required=True)
    parser.add_argument('--case', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--display-radius-A', dest='display_radius_A', type=float, default=3.5)
    print(render(**vars(parser.parse_args())))
