"""Report every declared multisite hybrid endpoint without an absolute decision band."""
from __future__ import annotations
import argparse
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_density_comparisons import PROTOCOL,validate_spec


def report(collection,output):
    r=read_json(collection);m=read_json(verify(r['manifest']))
    if r['protocol']!=PROTOCOL or not r['complete']:raise InvalidArtifact('complete actual grouped-panel collection required')
    spec,_=validate_spec(m['sources']['comparison_spec'],m['case_ids'],m['software'])
    primary=r['variants']['primary'];prep=read_json(verify(read_json(verify(read_json(verify(m['sources']['boundary']))['manifest']))['preparation']))
    cases={}
    for name,c in primary['cases'].items():
        mapping=read_json(verify(prep['cases'][name]));physical=read_json(verify(mapping['normalized_global_preparation']))
        state=read_json(verify(mapping['state']))
        cases[name]=dict(**c,source_mapping=prep['cases'][name],evidence_use=physical['evidence_use'],
            selected_metal_identity=physical['selected_metal_identity_alias'],background_metals=physical['background_metals'],
            explicit_water_union=physical['explicit_waters'],ligand_ledger=state['ligand_ledger'],
            physical_charge_by_endpoint={k:v['charge'] for k,v in physical['endpoints'].items()},
            endpoint_environment_correction_kcal={k:e['environment_correction_kcal'] for k,e in c['endpoints'].items()},
            radius_R_changes_kcal={v:r['variants'][v]['cases'][name]['R_kcal']-c['R_kcal'] for v in ('radius_minus','radius_plus')})
    result=dict(protocol=PROTOCOL,scientific_protocol=m['scientific_protocol'],collection=record(collection),
        plan=m['plan'],comparison_spec=m['sources']['comparison_spec'],complete=True,numerical_pass=r['numerical_pass'],
        radius_comparison_pass=r['radius_sensitivity_pass'],supporting_parvalbumin_ordering_pass=primary['ordering_pass'],
        supporting_comparisons=primary['contrasts'],ordered_vectors=primary['ordered_vectors'],cases=cases,
        aequorin_site_classifications=None,calibrated_class=None,reference=None,baseline_changed=False,
        prior_GGR_failure='retained; not superseded by this family extension',biological_groups_added=2,
        evidence_note='parvalbumin supporting cross-study; aequorin protein-level weak direction, site-unresolved',
        scope='consumed method development; no threshold fitting or prospective validation')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);write_new(root/'result.json',result)
    lines=['# Frozen multisite hybrid: complete ordered results','',
        'No absolute calibration or aequorin site classifications. The earlier GGR2FW0 failure remains.', '',
        f'Numerical checks: **{r["numerical_pass"]}**. Supporting parvalbumin ordering: **{primary["ordering_pass"]}**.', '',
        '| Site | Raw R (kcal/mol) |','|---|---:|']
    for group,ids in spec['ordered_vectors'].items():
        for name in ids:lines.append(f'| {name} | {cases[name]["R_kcal"]:.9f} |')
    lines+=['','The common raw metal-energy offset is retained; zero is not a decision boundary.', '',
        '| Supporting comparison | Difference (kcal/mol) | >0.02 |','|---|---:|---|']
    for c in primary['contrasts']:lines.append(f'| {c["left"]} − {c["right"]} | {c["difference_kcal"]:.9f} | {c["pass_"]} |')
    lines+=['','These six comparisons are two sites against three GGR structures, not six independent biological observations.',
        'Aequorin remains the ordered EF1/EF3/EF4 vector. Its global assay cannot label individual sites.',
        'All components, background ions, waters, charge ledgers, radius effects and execution pins are retained in result.json.']
    (root/'REPORT.md').write_text('\n'.join(lines)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--collection',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=report(a.collection,a.output);print({k:r[k] for k in ('complete','numerical_pass','supporting_parvalbumin_ordering_pass')})
