"""Printed energy bookkeeping for the approved GGR representation comparisons."""
import argparse
import importlib.util
from pathlib import Path

from affordable_common import HA_TO_KCAL, InvalidArtifact, record, write_new
from ggr_workflow import checked_rows


def audit(collections, historical, output):
    root = Path(__file__).resolve().parents[1]
    source = root / 'diagnostics/benchmark_set_20260915/audit_ggr.py'
    spec = importlib.util.spec_from_file_location('archived_ggr_audit', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    current = [r for p in collections for r in checked_rows(p)['rows']]
    previous = [r for p in historical for r in checked_rows(p)['rows']]
    needed = {r.get('historical_case') for r in current}
    selected = [r for r in previous if r['case'] in needed and
                (r['protocol_id'] == 'generic_peptide_amide_vertical_native_r2scan3c_v3'
                 or r['case'] == 'ggr_1glg_GGR')]
    rows = []
    for row in selected + current:
        item = {k: row.get(k) for k in ('case', 'lane', 'protocol_id', 'status',
                                        'source_structure_id', 'score')}
        item['endpoints'] = None
        item['contrast_components_kcal_mol'] = None
        if row['status'] == 'complete':
            endpoints = {m: module.components(e) for m, e in row['endpoints'].items()}
            item['endpoints'] = endpoints
            item['contrast_components_kcal_mol'] = {
                k: (endpoints['Ca']['hartree'][k] - endpoints['La']['hartree'][k]) * HA_TO_KCAL
                for k in endpoints['Ca']['hartree']}
        rows.append(item)
    differences = []

    def difference(before, after):
        valid = before['status'] == after['status'] == 'complete'
        values = {k: after['contrast_components_kcal_mol'][k] - v
                  for k, v in before['contrast_components_kcal_mol'].items()} if valid else None
        if valid and abs(values['final'] - (after['score']['S_kcal_mol'] - before['score']['S_kcal_mol'])) > 1e-7:
            raise InvalidArtifact('component contrast difference algebra')
        differences.append({'from_case': before['case'], 'from_protocol': before['protocol_id'],
                            'to_case': after['case'], 'to_protocol': after['protocol_id'],
                            'status': 'complete' if valid else 'unavailable',
                            'delta_R_components_kcal_mol': values})

    for new in current:
        matches = [r for r in rows if r['case'] == new.get('historical_case') and
                   r['protocol_id'] == 'generic_peptide_amide_vertical_native_r2scan3c_v3']
        if matches:
            if len(matches) != 1:
                raise InvalidArtifact('ambiguous historical component comparison')
            after = next(r for r in rows if r['case'] == new['case'])
            difference(matches[0], after)
    for oldname, newname in [('ggr_1glg_nma', 'ggr_1glg_connected'),
                             ('ggr_2fw0_formamide', 'ggr_2fw0_nma'),
                             ('ggr_2fvy_formamide', 'ggr_2fvy_nma')]:
        before, after = ([r for r in rows if r['case'] == name] for name in (oldname, newname))
        if before and after:
            if len(before) != 1 or len(after) != 1:
                raise InvalidArtifact('ambiguous matched representation comparison')
            difference(before[0], after[0])
    result = {'schema_version': 'alquemia.ggr_component_audit.v1',
              'implementation': record(__file__), 'reused_component_parser': record(source),
              'collections': [record(p) for p in collections],
              'historical_collections': [record(p) for p in historical],
              'hartree_to_kcal_mol': HA_TO_KCAL, 'rows': rows, 'differences': differences,
              'interpretation': 'Printed CPCM and SCF-minus-CPCM are bookkeeping on polarized densities, not independent solvent/vacuum calculations or a causal decomposition. Nuclear geometry and fragment membership change the density and cavity together. No components are removed from the score.'}
    write_new(output, result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--collection', type=Path, action='append', required=True)
    parser.add_argument('--historical-collection', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.collection, args.historical_collection, args.output)
