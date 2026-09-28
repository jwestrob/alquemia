"""Collect only admitted repaired origins; missing cells never become zero."""
import argparse
import json
import sys
import math
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import HA_TO_KCAL, InvalidArtifact, read_json, record, verify, write_new
from metal_environment_force_assembly import prepare_mapping


def compare(config):
    declaration = read_json(config)
    if declaration.get('primary_declaration'):
        verify(declaration['primary_declaration'])
    rows = []
    recipes = []
    seen = set()
    source_inputs = {}
    for entry in declaration['cells']:
        key = (entry['source'], entry['metal'])
        if key in seen:
            raise InvalidArtifact('duplicate source/metal cell')
        seen.add(key)
        manifest_path = Path(entry['manifest'])
        m = read_json(manifest_path)
        tasks = [t for t in m['tasks'] if t['task_id'] == entry['task_id']]
        if len(tasks) != 1:
            raise InvalidArtifact('declared task missing or duplicated')
        t = tasks[0]
        if m['source_id'] != entry['source_id'] or t['metal'] != entry['metal'] or t['configuration'] != 'A':
            raise InvalidArtifact('source/state/origin differs')
        verify(m['inputs']); verify(t['xyz']); verify(t['pointcharges'])
        previous = source_inputs.setdefault(entry['source'], m['inputs']['sha256'])
        if previous != m['inputs']['sha256']:
            raise InvalidArtifact('paired source preparation differs')
        recipes.append((m['method'], m['orca']['sha256'],
                        tuple((z, k, m['assets'][z][k]['sha256']) for z in ('La', 'Dy') for k in ('basis', 'aux'))))
        row = dict(source=entry['source'], source_id=m['source_id'], metal=t['metal'],
                   task_id=t['task_id'], manifest=record(manifest_path), inputs=m['inputs'],
                   electronic_state=t['electronic_state'], initial_guess=m.get('initial_guess', 'PModel'),
                   solver=m.get('solver', 'default'), status='unavailable', energy_hartree=None,
                   reason='collection_not_available')
        row.update(classical_status='unavailable', classical_kcal_mol=None)
        classical_path=Path(entry['classical']) if entry.get('classical') else None
        if classical_path is not None and classical_path.exists():
            classical=read_json(classical_path)
            if classical['inputs'] != m['inputs'] or classical['metal'] != t['metal']:
                raise InvalidArtifact('classical source/state differs')
            if classical.get('QM_MM_Coulomb_included') is not False or classical.get('C4_induction_included') is not False:
                raise InvalidArtifact('classical electrostatic exclusions differ')
            verify(classical['ledger'])
            physical=prepare_mapping(verify(m['inputs']), 'A', t['metal'])['physical_atoms']
            indices={sid:i for i,sid in enumerate(classical['source_ids'])}
            if len(indices)!=len(classical['source_ids']) or set(indices)!={a['id'] for a in physical}:
                raise InvalidArtifact('classical physical membership differs')
            coords=np.asarray(classical['coordinates_A'])
            if coords.shape!=(len(physical),3) or not np.allclose(coords[[indices[a['id']] for a in physical]], [a['xyz_A'] for a in physical], atol=1e-8,rtol=0):
                raise InvalidArtifact('classical coordinates differ')
            total=classical['total_kcal_mol']
            if not math.isfinite(total) or abs(sum(classical['components_kcal_mol'].values())-total)>1e-8:
                raise InvalidArtifact('classical energy components inconsistent')
            row.update(classical_status='complete', classical_kcal_mol=total, classical=record(classical_path),
                       classical_components_kcal_mol=classical['components_kcal_mol'])
        path = Path(entry['collection'])
        if path.exists():
            c = read_json(path)
            if c['manifest'] != record(manifest_path) or c['source_id'] != m['source_id']:
                raise InvalidArtifact('collection manifest/source mismatch')
            r = c['rows'].get(t['task_id'])
            if r is None:
                raise InvalidArtifact('collection omits declared task')
            row.update(collection=record(path), status=r['status'], reason=r.get('reason'))
            if r['status'] == 'complete':
                if r['electronic_state'] != t['electronic_state']:
                    raise InvalidArtifact('collected state mismatch')
                for name in ('output', 'engrad', 'pointcharge_gradient', 'receipt'):
                    verify(r[name])
                row['energy_hartree'] = r['energy_hartree']
        rows.append(row)
    expected = {(s, z) for s in declaration['sources'] for z in ('La', 'Dy')}
    if seen != expected or len(set(recipes)) != 1:
        raise InvalidArtifact('declared denominator or electronic recipe differs')
    contrasts = {}
    for source in declaration['sources']:
        pair = {r['metal']: r for r in rows if r['source'] == source}
        contrasts[source] = ((pair['Dy']['energy_hartree'] - pair['La']['energy_hartree']) * HA_TO_KCAL
                             if all(r['status'] == 'complete' for r in pair.values()) else None)
    differences = {}
    for hans in declaration['hans_sources']:
        mex = declaration['mex_source']
        differences[hans] = (contrasts[hans] - contrasts[mex]
                             if contrasts[hans] is not None and contrasts[mex] is not None else None)
    classical_contrasts={}
    for source in declaration['sources']:
        pair={r['metal']:r for r in rows if r['source']==source}
        classical_contrasts[source]=(pair['Dy']['classical_kcal_mol']-pair['La']['classical_kcal_mol']
                                     if all(r['classical_status']=='complete' for r in pair.values()) else None)
    finite_contrasts={s:contrasts[s]+classical_contrasts[s]
                      if contrasts[s] is not None and classical_contrasts[s] is not None else None
                      for s in declaration['sources']}
    finite_differences={s:finite_contrasts[s]-finite_contrasts[declaration['mex_source']]
                        if finite_contrasts[s] is not None and finite_contrasts[declaration['mex_source']] is not None else None
                        for s in declaration['hans_sources']}
    return dict(declaration=record(config), implementation=record(__file__), rows=rows,
                analysis_role=declaration.get('analysis_role','primary'),
                primary_declaration=declaration.get('primary_declaration'),
                replaced_attempts=declaration.get('replaced_attempts',[]),
                complete_cells=sum(r['status'] == 'complete' for r in rows), declared_cells=len(expected),
                raw_Dy_minus_La_kcal_mol=contrasts, Hans_minus_Mex_exchange_kcal_mol=differences,
                classical_Dy_minus_La_kcal_mol=classical_contrasts,
                finite_electronic_plus_classical_Dy_minus_La_kcal_mol=finite_contrasts,
                finite_Hans_minus_Mex_exchange_kcal_mol=finite_differences,
                sign_convention='Positive exchange means conditional greater relative La preference in Hans',
                energy_scope='finite embedded electronic component only', affinity=None, classification=None,
                finite_energy_scope='embedded electronic plus explicitly matched classical terms; no bulk solvent',
                limitations=declaration['limitations'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = compare(args.config)
    write_new(args.output, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))
