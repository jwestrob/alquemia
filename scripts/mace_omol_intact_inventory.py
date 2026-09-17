"""Inventory intact canonical sources and formal charges; no energy evaluation."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import resource
import time
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from mace_global_prepare import protein, FF, STANDARD

INVENTORY_SHA = 'f51a8641fa50a34087fca9ba76d0e5ca960ac85c1ea5da65a3fa19bed9d5cb21'


def audit(source_inventory, agreement, output):
    start = time.monotonic()
    source = record(source_inventory)
    if source['sha256'] != INVENTORY_SHA:
        raise InvalidArtifact('requires the declared canonical inventory')
    data = read_json(source_inventory)
    if len(data['rows']) != 28 or data['calibration_count'] != 25:
        raise InvalidArtifact('canonical inventory count changed')
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for original in data['rows']:
        name = original['case_id']
        result = {key: original[key] for key in
                  ('case_id', 'expected_class', 'evaluation_role', 'evidence_stratum',
                   'sequence_accession_group', 'prospectively_blind')}
        try:
            pm_ref = original['preparation_provenance']['protonation_manifest']
            pm = read_json(verify(pm_ref))
            source_pdb = pm['output']
            verify(source_pdb)
            if pm.get('detected_missing_residues'):
                raise InvalidArtifact('missing residues require explicit preparation scope')
            if original['explicit_water_inventory']:
                raise InvalidArtifact('unexpected canonical explicit waters')
            if original['cofactor']['formal_charge'] != -3:
                raise InvalidArtifact('unexpected PQQ formal charge')
            prep_row = {'case_id': 'PQQ_' + name if name in ('1H4I', '4MAE') else name,
                        'source_structure': source_pdb,
                        'selected_site': original['assembly']['selected_site']}
            atoms, charge, details = protein(prep_row)
            excluded = [r for r in details['source_residue_inventory'] if not r['selected_protein']]
            selected = [r for r in details['source_residue_inventory'] if r['selected_protein']]
            if any(r['resname'] not in STANDARD for r in selected):
                raise InvalidArtifact('nonstandard residue entered protein topology')
            endpoint_charges = {'La': charge, 'Ca': charge - 1}
            result.update(status='protein_template_supported', protein_atoms=len(atoms),
                          prospective_total_atoms=len(atoms) + 28,
                          protein_charge_e=charge, cofactor_charge_e=-3,
                          endpoint_charges=endpoint_charges,
                          within_reported_training_charge_range=all(-10 <= q <= 10 for q in endpoint_charges.values()),
                          source_protonation_manifest=pm_ref, source_protonated_pdb=source_pdb,
                          selected_chain='A', protein_residue_count=len(selected),
                          excluded_residues=excluded, terminal_additions=details['terminal_additions'],
                          hydrogen_length_changes=len(details['protein_H_moves']),
                          heavy_coordinates_preserved=details['source_heavy_coordinates_preserved'],
                          forcefield_charge_sum_e=details['ff_charge_sum_e'],
                          final_prepared_coordinates=None, energy=None)
        except Exception as exc:
            result.update(status='unsupported', reason=str(exc), exception_type=type(exc).__name__)
        rows.append(result)
        print(name, result['status'], result.get('endpoint_charges'), result.get('reason', ''), flush=True)
    result = {'status':'complete', 'scope':'preparation_readiness_only', 'source_inventory':source,
              'agreement':record(agreement), 'forcefield':record(FF), 'rows':rows,
              'new_energy_evaluations':0, 'new_protonations':0,
              'implementation':{p.name:record(p) for p in Path(__file__).parent.glob('*.py')},
              'wall_seconds':time.monotonic()-start,
              'process_CPU_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime}
    write_new(out/'result.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('source-inventory', 'agreement', 'output'):
        parser.add_argument('--'+key, required=True)
    audit(**vars(parser.parse_args()))
