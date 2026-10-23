"""Select the predeclared water subset from the unchanged full Stage A scorer."""

def evaluate(data, model, body, anchors, stage):
    water = tuple(q for q in stage.CONFIGURATIONS if q.startswith('validation_aquo8_'))
    water_cells = {q + '_' + m for q in water for m in stage.METALS}
    stage.require(len(water) == 2 and len(water_cells) == 4, 'Water denominator differs')
    stage.require(set(data['rows']) == set(model) == set(stage.CELLS), 'All16 slots required')
    for gid in set(stage.CELLS) - water_cells:
        stage.require(model[gid] is None and data['rows'][gid].get('reference') is None,
                      'Donor outcomes are outside this water-only collection: ' + gid)
    full = stage.evaluate(data, model, body, anchors)
    groups = {
        'balanced_exchange': {q: full['balanced_exchange'][q] for q in water},
        'orientation_work': {m: full['orientation_work']['aquo8_' + m]
                             for m in ('La', 'Dy', 'paired')},
        'rigid_response': {q + '_' + m: full['per_fixture'][m][q]
                           for q in water for m in ('La', 'Dy', 'paired')},
    }
    counts = {k: {'declared': len(v),
                  'available': sum(x['status'] == 'complete' for x in v.values()),
                  'passed': sum(x['passed'] is True for x in v.values()),
                  'failed': sum(x['passed'] is False for x in v.values())}
              for k, v in groups.items()}
    comparisons = [v for group in groups.values() for v in group.values()]
    complete = all(full['cells'][gid]['status'] == 'complete' for gid in water_cells)
    complete = complete and all(v['status'] == 'complete' for v in comparisons)
    stage.require(len(comparisons) == 11 and full['status'] == 'partial'
                  and full['validation_passed'] is False, 'Full validation cannot pass here')
    return {'status': 'water_subset_complete' if complete else 'water_subset_partial',
            'water_subset_passed': complete and all(v['passed'] is True for v in comparisons),
            'water_failure_observed': any(v['passed'] is False for v in comparisons),
            'counts': counts, 'comparisons': groups, 'full_stage_A': full,
            'water_cells': 4, 'water_comparisons': 11, 'original_reserved_cells': 24,
            'full16_validation_passed': False, 'protein_sampling_released': False,
            'physical_accuracy_qualified': False, 'production_released': False}
