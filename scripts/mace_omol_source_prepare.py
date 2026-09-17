"""Build audited intact inputs from explicit archived sources, without new folds."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_global_prepare import STANDARD, prepare_case

POLICY = 'omol_source_backed_chain_A_ff19sb_H_explicit_exclusions_v1'
EXCLUDABLE = {'HOH', 'GAL', 'GLC', 'BGC', 'GOL', 'CIT', 'MLA', 'NA', 'ACT', 'CO2'}
METALS = {'CA', 'LA', 'CE', 'PR', 'ND', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'Y'}


def inspect_source(row, exclusions):
    """Reject unsupported polymer/cofactor chemistry before filtered FF matching."""
    import gemmi
    verify(row['source_structure']); core = read_json(verify(row['preparation_manifest']))
    pm_pin = row.get('source_protonation_manifest', row.get('protonation_manifest'))
    if not pm_pin:
        raise InvalidArtifact('archived protonation provenance required')
    pm = read_json(verify(pm_pin))
    if pm['output'] != row['source_structure'] or pm['ph'] != 7.0 or pm['add_missing_residues']:
        raise InvalidArtifact('source differs from the declared archived pH7 preparation')
    if pm.get('repaired_missing_atom_count', 0):
        raise InvalidArtifact('source contains reconstructed nonterminal heavy atoms')
    raw = row.get('raw_source', pm['source']); structure = gemmi.read_structure(str(verify(raw)))
    structure.setup_entities()
    residues = [r for c in structure[0] if c.name == 'A' for r in c]
    if not residues:
        raise InvalidArtifact('raw source lacks chain A')
    observed = []; omitted = []; metals = []; metal_coordinates = []
    pqq = bool(core.get('pqq'))
    allowed = set(exclusions)
    if len(allowed) != len(exclusions) or not allowed <= EXCLUDABLE:
        raise InvalidArtifact('undeclared/unsupported excluded heterogen chemistry')
    for r in residues:
        item = {'resname': r.name, 'resid': str(r.seqid), 'atoms': len(r), 'entity_type': r.entity_type.name}
        observed.append(item)
        if r.entity_type == gemmi.EntityType.Polymer and r.name not in STANDARD:
            raise InvalidArtifact('unsupported source polymer residue: '+r.name+str(r.seqid))
        if r.name in STANDARD:
            continue
        if r.name in METALS:
            if len(r) != 1:
                raise InvalidArtifact('unsupported multi-atom metal residue')
            metals.append(item); metal_coordinates.append(list(r[0].pos)); continue
        if r.name == 'PQQ' and pqq:
            continue
        if r.name not in allowed:
            raise InvalidArtifact('source heterogen lacks an explicit supported exclusion: '+r.name)
        omitted.append(item)
    if len(metals) != 1:
        raise InvalidArtifact('single metal-bearing chain required; retain other metals through a separate policy')
    selected = row['selected_site']['xyz_A']
    if any(abs(x-y)>1e-6 for x,y in zip(metal_coordinates[0], selected)):
        raise InvalidArtifact('selected site does not match the actual source metal')
    if {r['resname'] for r in omitted} != allowed:
        raise InvalidArtifact('exclusion inventory differs from actual source')
    for bond in structure.connections:
        if bond.type == gemmi.ConnectionType.Covale and (bond.partner1.chain_name == 'A' or bond.partner2.chain_name == 'A'):
            if bond.partner1.res_id.name not in STANDARD or bond.partner2.res_id.name not in STANDARD:
                raise InvalidArtifact('unsupported covalent nonstandard source connection: '+bond.name)
    la, ca = (xyz(verify(row['endpoints'][m]['xyz'])) for m in ('La', 'Ca'))
    if not la or la[0][0] != 'La' or ca[0][0] != 'Ca' or la[0][1:] != ca[0][1:] or la[1:] != ca[1:]:
        raise InvalidArtifact('source core pair has different physical coordinates or membership')
    if any(abs(x-y)>1e-6 for x,y in zip(la[0][1:], selected)):
        raise InvalidArtifact('source core metal differs from selected physical site')
    if pqq and not row['case_id'].startswith('PQQ'):
        raise InvalidArtifact('PQQ source requires the existing explicit PQQ preparation route')
    return {'raw_source': raw, 'protonation': pm_pin, 'source_residues': observed,
            'excluded_heterogens': omitted, 'selected_metal': metals[0],
            'missing_residues_reported': pm.get('detected_missing_residues'),
            'scope': 'observed chain A; exclusions are declared model simplifications, not complete holo assembly'}


def validate_bridge(preparation):
    b = preparation['source_bridge']
    for k in ('agreement', 'implementation'):
        verify(b[k])
    for pin in b['implementation_files'].values():
        verify(pin)
    row = read_json(verify(b['source_row'])); evidence = read_json(verify(b['evidence']))
    exclusions = read_json(verify(b['exclusions']))
    if (preparation['policy_id'] != POLICY or preparation['source_audit_row'] != row
            or preparation['evidence'] != evidence or preparation['evidence_use'] != b['evidence_use']
            or len(preparation['explicit_waters']) != b['expected_water_count']
            or b['source_audit'] != inspect_source(row, exclusions)):
        raise InvalidArtifact('source-backed preparation differs from its declared inputs')


@cached_file_checks
def prepare(source_row, evidence, exclusions, expected_waters, evidence_use, agreement, output):
    row = read_json(source_row); label = read_json(evidence); omitted = read_json(exclusions)
    audit = inspect_source(row, omitted)
    if expected_waters < 0:
        raise InvalidArtifact('expected water count must be nonnegative')
    # Output coordinates/chemistry use the same established implementation.
    ref = prepare_case(row, output, allow_terminal_completion=True,
                       expected_water_count=expected_waters, evidence=label,
                       policy_id=POLICY, check_peptide_connectivity=True)
    path = verify(ref); p = read_json(path)
    from mace_global_benchmark import snapshot
    names = tuple(f.name for f in Path(__file__).parent.glob('mace_omol*.py')) + (
        'mace_canonical_run.py', 'mace_canonical.py', 'mace_curvature.py',
        'affordable_response.py', 'mace_mechanics_run.py', 'mace_short_engine.py', 'mace_file_checks.py')
    implementation = snapshot(path.parent, names)
    p['evidence_use'] = evidence_use
    p['source_bridge'] = {'source_row': record(source_row), 'evidence': record(evidence),
                          'exclusions': record(exclusions), 'expected_water_count': expected_waters,
                          'evidence_use': evidence_use, 'agreement': record(agreement),
                          'implementation': implementation[Path(__file__).name],
                          'implementation_files': implementation, 'source_audit': audit}
    # Keep the emitted physical preparation immutable; the public bridge record
    # is a separate file with additional provenance, referring to the same XYZs.
    final = path.parent/'source_preparation.json'; write_new(final, p)
    from mace_omol_prepared import audit_preparation
    checked = audit_preparation(final)
    write_new(path.parent/'source_audit.json', checked)
    return {'status': 'prepared', 'preparation': record(final), 'audit': checked}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('source-row', 'evidence', 'exclusions', 'agreement', 'output'):
        parser.add_argument('--'+key, required=True)
    parser.add_argument('--expected-waters', type=int, required=True)
    parser.add_argument('--evidence-use', required=True,
                        choices=('consumed_method_development', 'threshold_calibration', 'prospective_evaluation'))
    print(json.dumps(prepare(**vars(parser.parse_args())), indent=2))
