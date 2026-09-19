#!/usr/bin/env python3
"""Reuse the frozen crystal fixed-core preparer for deposited 8GY2 chain A.

No energies, coordinate optimization, invented sites or method changes.
The selector operation only inventories the already identified source site.
"""
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
CAL = ROOT / 'diagnostics/pqq_pmdh_fixed_core_calibration_20260914'
spec = importlib.util.spec_from_file_location('frozen_crystal_preparer', CAL / 'reserved_crystal_holdout/prepare_holdouts.py')
h = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = h
spec.loader.exec_module(h)


def selectors(source, destination):
    source = source.resolve()
    st = h.gemmi.read_structure(str(source)); model = st[0]; chain = model['A']
    metal_selector = h.parse_atom_selector('A:CA803/CA')
    pqq_selector = h.parse_residue_selector('A:PQQ802')
    metal = h.resolve_atom(model, metal_selector)
    row = dict(holdout_role='external_association_control', pdb_id='8GY2',
        biological_class='Ca/membrane_ADH_structural_functional_association',
        source_path=str(source), source_sha256=h.sha256_file(source), source_format='mmCIF',
        resolution_A=str(st.resolution), model='1', selected_chain='A',
        chain_selection_basis='Only catalytic ADH subunit A in deposited heterotrimer; source site identified by PQQ802 and Ca803, independently of scores.',
        metal_selector=metal_selector.label(), metal_element=metal.element.name,
        metal_occupancy=str(metal.occ), metal_bfactor_A2=str(metal.b_iso),
        metal_xyz_A=','.join(str(v) for v in (metal.pos.x,metal.pos.y,metal.pos.z)),
        pqq_selector=pqq_selector.label(), pqq_schema='pdb_ccd_pqq_v1',
        glu_selector='A:GLU215', asn_selector='A:ASN297', catalytic_asp_selector='A:ASP342',
        plus2_selector='A:THR344', include_plus2='false', cationic_partner_selector='A:LYS369',
        core_altlocs='none', noncore_direct_ligand='none', noncore_direct_ligand_metal_A='NA',
        dry_policy='dry_fixed_core_v3: retain selected metal, complete PQQ and homologous protein fragments; exclude source waters, other chains and noncore HEC cofactors without replacement; do not fill vacancies',
        existing_protonated_path='NA',existing_protonated_sha256='NA',protonated_metal_selector='NA',protonated_pqq_selector='NA',
        notes='Ca-associated functional control, not direct La/Ca affinity or proof of La exclusion. Whole-protein MACE preparation is unavailable because heme chemistry is not covered by its frozen preparation.')
    for role, field in h.ROLE_FIELDS.items():
        selector=h.parse_residue_selector(row[field]); residue=h.resolve_residue(model,selector)
        if role in h.ROLE_DISTANCE_FIELDS:
            key=h.fixed_key(selector,role)
            distances=h.fixed.donor_distances(residue,key,metal.pos)
            row[h.ROLE_DISTANCE_FIELDS[role]]=str(distances[0]['distance_A']) if distances else 'NA'
    a=h.resolve_residue(model,h.parse_residue_selector(row['catalytic_asp_selector']))
    b=h.resolve_residue(model,h.parse_residue_selector(row['cationic_partner_selector']))
    partner=h.fixed.asp_partner_geometry(a,b)
    row['cat_asp_partner_contact']=row['catalytic_asp_selector']+'/'+partner['asp_atom']+'--'+row['cationic_partner_selector']+'/'+partner['partner_atom']
    row['cat_asp_partner_A']=str(partner['heavy_atom_distance_A'])
    contacts=h.typed_source_contacts(model,metal_selector,metal)
    row['fixed_typed_CN_3p1']=str(len(contacts))
    row['typed_donor_ledger_3p1_A']=';'.join(x['selector'].label()+'='+str(x['distance_A']) for x in contacts)
    row.update({k: 'NA' if v is None else str(v) for k,v in h.water_inventory(model,'A',metal).items()})
    row['source_altloc_residues']=';'.join(h.altloc_inventory(model)) or 'none'
    raw=h.audit_raw_source(row)
    destination.mkdir(parents=True,exist_ok=False)
    with (destination/'target_spec.tsv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t',lineterminator='\n');w.writeheader();w.writerow(row)
    audit={'source':h.file_record(source),'selector_audit_passed':True,'typed_donors':raw['donor_records'],
           'waters':raw['water_record'],'source_metal':raw['source_metal'],
           'interpretation':'Ca-associated PQQ functional/structural control; no direct La/Ca affinity label; no energy calculation.'}
    (destination/'SELECTOR_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
    return row


def prepare(selector_dir, output):
    selector_dir=selector_dir.resolve();output=output.resolve()
    rows=list(csv.DictReader((selector_dir/'target_spec.tsv').open(),delimiter='\t'));assert len(rows)==1
    pins_path=CAL/'implementation_pins.json';pins=h.fixed.verify_pins(pins_path)
    h.verify_runtime(pins)
    historical=ROOT/'diagnostics/pqq_q46444_1kb0_external_validation_20260915/prepared/external_preparation.json'
    gate=json.loads(historical.read_text())['calibration_gate']
    h.verify_file_record(gate['result'],'published calibration result')
    output.mkdir(parents=True,exist_ok=False)
    newpins={'protocol_id':h.PROTOCOL_ID,'holdout_spec':h.file_record(selector_dir/'target_spec.tsv'),
             'holdout_selector_audit':h.file_record(selector_dir/'SELECTOR_AUDIT.json'),
             'holdout_implementation':{'prepare_holdouts':h.file_record(Path(h.__file__))},
             'wrapper':h.file_record(Path(__file__)),'calibration_implementation_pins':h.file_record(pins_path),
             'agreement':h.file_record(Path(__file__).with_name('AGREEMENT.md'))}
    pp=output/'implementation_pins.json';pp.write_text(json.dumps(newpins,indent=2)+'\n')
    start=time.monotonic();cpu=time.process_time()
    try:
        target=h.prepare_one(rows[0],1,output,gate,pp,newpins,pins_path,pins)
        result={'status':'prepared_unscored','target':target,'protocol_id':h.PROTOCOL_ID,'energy_evaluations':0}
    except Exception as e:
        result={'status':'preparation_failed','failure_type':type(e).__name__,'failure_reason':str(e),'energy_evaluations':0}
    result.update(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu)
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if result['status']!='prepared_unscored':raise SystemExit(1)


def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('selectors');s.add_argument('--source',type=Path,required=True);s.add_argument('--output',type=Path,required=True)
    s=sub.add_parser('prepare');s.add_argument('--selectors',type=Path,required=True);s.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.command=='selectors':selectors(args.source,args.output)
    else:prepare(args.selectors,args.output)

if __name__=='__main__':main()
