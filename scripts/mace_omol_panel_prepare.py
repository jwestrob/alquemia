"""Prepare the declared intact canonical panel from archived protonated sources."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import resource
import time
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_global_prepare import prepare_case, FF
from mace_hybrid import check_atoms
from mace_omol_intact import preparations
from mace_omol_intact_inventory import INVENTORY_SHA

POLICY = 'omol_canonical_chain_A_ff19sb_H_terminal_OXT_v1'
SETTINGS = {'chain':'A', 'protonation':'archived', 'protein_H':'ff19SB_radial',
            'terminal_OXT':'last_residue_C_CA_O_connectivity_reflection_1.25A',
            'PQQ_charge':-3, 'explicit_waters':0, 'charge_clipping':False}
ACTIVE_POLICY = 'omol_canonical_chain_A_ff19sb_H_peptide_connectivity_v2'
ACTIVE_SETTINGS = {**SETTINGS, 'peptide_C_N_limits_A':[1.0,1.8]}


def source_rows(path):
    if record(path)['sha256'] != INVENTORY_SHA:
        raise InvalidArtifact('canonical source inventory changed')
    source = read_json(path)
    if len(source['rows']) != 28 or sum(r['evaluation_role']=='calibration' for r in source['rows']) != 25:
        raise InvalidArtifact('canonical source roles/count changed')
    return source['rows']


def prepare(source_inventory, physical_preparation, agreement, output):
    start = time.monotonic()
    sources = source_rows(source_inventory)
    existing = preparations(physical_preparation)
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for source in sources:
        name = source['case_id']
        row = {k:source[k] for k in ('case_id', 'expected_class', 'evaluation_role',
                                    'sequence_accession_group', 'evidence_stratum', 'prospectively_blind')}
        try:
            if source['explicit_water_inventory'] or source['cofactor']['formal_charge'] != -3:
                raise InvalidArtifact('canonical water/cofactor state differs')
            if name in ('1H4I','4MAE'):
                ref = existing['PQQ_'+name]['source']
                row['preparation_reused'] = True
            else:
                pm_ref = source['preparation_provenance']['protonation_manifest']
                pm = read_json(verify(pm_ref))
                verify(pm['output'])
                if pm.get('detected_missing_residues'):
                    raise InvalidArtifact('missing residues outside preparation scope')
                prep_row = {'case_id':'PQQ_'+name, 'source_structure':pm['output'],
                            'preparation_manifest':source['source_manifest'],
                            'selected_site':source['assembly']['selected_site'],
                            'endpoints':source['endpoints'], 'protonation_manifest':pm_ref}
                evidence = {'evidence_stratum':source['evidence_stratum'],
                            'direction':source['expected_class'], 'group':source['sequence_accession_group']}
                ref = prepare_case(prep_row, out/name, allow_terminal_completion=True,
                                   expected_water_count=0, evidence=evidence, policy_id=ACTIVE_POLICY,
                                   check_peptide_connectivity=True)
                row['preparation_reused'] = False
            prepared = read_json(verify(ref))
            charges = {m:prepared['endpoints'][m]['charge'] for m in ('La','Ca')}
            row.update(status='prepared', preparation=ref, atoms=len(prepared['physical_atoms']),
                       endpoint_charges=charges,
                       outside_reported_training_charge_range=any(not -10<=q<=10 for q in charges.values()),
                       outside_reported_training_size_range=len(prepared['physical_atoms'])>350)
        except Exception as exc:
            row.update(status='unsupported', reason=str(exc), exception_type=type(exc).__name__,
                       score=None)
        rows.append(row)
        print(name, row['status'], row.get('endpoint_charges'), row.get('reason',''), flush=True)
    result = {'schema_version':'alquemia.omol_intact_panel_preparation.v1', 'policy_id':ACTIVE_POLICY,
              'settings':ACTIVE_SETTINGS, 'status':'complete' if all(r['status']=='prepared' for r in rows) else 'incomplete',
              'source_inventory':record(source_inventory), 'physical_preparation':record(physical_preparation),
              'agreement':record(agreement), 'rows':rows, 'forcefield':record(FF),
              'implementation':{p.name:record(p) for p in Path(__file__).parent.glob('*.py')},
              'new_energy_evaluations':0, 'new_protonations':0,
              'wall_seconds':time.monotonic()-start,
              'process_CPU_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime}
    write_new(out/'preparation_manifest.json', result)
    validate(out/'preparation_manifest.json')
    return result


def validate(path):
    data = read_json(path)
    policies={POLICY:SETTINGS, ACTIVE_POLICY:ACTIVE_SETTINGS}
    if (data['schema_version']!='alquemia.omol_intact_panel_preparation.v1'
            or data['policy_id'] not in policies or data['settings']!=policies.get(data['policy_id']) or data['forcefield']!=record(FF)):
        raise InvalidArtifact('intact panel preparation policy changed')
    sources = source_rows(verify(data['source_inventory']))
    existing = preparations(verify(data['physical_preparation']))
    for pin in [data['agreement'], *data['implementation'].values()]:verify(pin)
    if [r['case_id'] for r in data['rows']] != [r['case_id'] for r in sources]:
        raise InvalidArtifact('intact panel inventory changed')
    for row, source in zip(data['rows'], sources):
        for key in ('expected_class','evaluation_role','sequence_accession_group','evidence_stratum','prospectively_blind'):
            if row[key]!=source[key]:raise InvalidArtifact('source evidence/role changed')
        if row['status']=='unsupported':
            if not row.get('reason') or row.get('score') is not None:
                raise InvalidArtifact('unsupported preparation must retain reason and unavailable score')
            continue
        if row['status']!='prepared':raise InvalidArtifact('invalid preparation status')
        p = read_json(verify(row['preparation']))
        name = row['case_id']
        if name in ('1H4I','4MAE'):
            if row['preparation']!=existing['PQQ_'+name]['source'] or row['preparation_reused'] is not True:
                raise InvalidArtifact('existing full crystal preparation was changed')
        else:
            pm = read_json(verify(source['preparation_provenance']['protonation_manifest']))
            expected_evidence={'evidence_stratum':source['evidence_stratum'],
                               'direction':source['expected_class'],'group':source['sequence_accession_group']}
            if (p['policy_id']!=data['policy_id'] or p['case_id']!='PQQ_'+name or row['preparation_reused'] is not False
                    or p['source']!=pm['output'] or p['source_preparation']!=source['source_manifest']
                    or p['evidence']!=expected_evidence):
                raise InvalidArtifact('canonical whole-chain source/policy/evidence changed')
        for key in ('source','source_preparation','forcefield','water_forcefield'):verify(p[key])
        atoms = p['physical_atoms']
        if (p['assembly']!='deposited_chain_A_with_declared_cofactor_and_site_waters'
                or p['explicit_waters'] or p['cofactor_charge_e']!=-3
                or len(atoms)!=row['atoms'] or len({a['id'] for a in atoms})!=len(atoms)
                or [i for i,a in enumerate(atoms) if a['element']=='M']!=[len(atoms)-1]
                or any('cap' in a['kind'] for a in atoms)):
            raise InvalidArtifact('intact physical inventory changed')
        for metal in ('La','Ca'):
            endpoint = p['endpoints'][metal]
            actual = xyz(verify(endpoint['xyz']))
            expected = [(metal if a['element']=='M' else a['element'],*a['xyz_A']) for a in atoms]
            q = p['protein_charge_e'] + p['cofactor_charge_e'] + (3 if metal=='La' else 2)
            if (actual!=expected or endpoint['charge']!=q or endpoint['spin_multiplicity']!=1
                    or endpoint['state']!=check_atoms(actual,q) or row['endpoint_charges'][metal]!=q):
                raise InvalidArtifact('paired coordinates, charge or electron state changed')
        if (row['outside_reported_training_charge_range'] != any(not -10<=q<=10 for q in row['endpoint_charges'].values())
                or row['outside_reported_training_size_range'] != (len(atoms)>350)):
            raise InvalidArtifact('training-domain flag changed')
    if (data['status']=='complete') != all(r['status']=='prepared' for r in data['rows']):
        raise InvalidArtifact('preparation completeness is inconsistent')
    if data['policy_id']==ACTIVE_POLICY:
        from mace_omol_backbone_audit import inspect
        audit=inspect(data)
        if any(a['status']!='pass' and r['status']=='prepared' for a,r in zip(audit['rows'],data['rows'])):
            raise InvalidArtifact('new preparation claims unsupported peptide geometry is prepared')
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('source-inventory','physical-preparation','agreement','output'):
        parser.add_argument('--'+key, required=True)
    prepare(**vars(parser.parse_args()))
