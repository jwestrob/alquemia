"""Opt-in APBS electrostatic transfer of one frozen QM charge distribution.

The six charging calculations use identical grid dimensions, spacings and
centres. Every reaction-field difference cancels discretization self energies.
No second CPCM solvation term is added. See the diagnostic MODEL.md.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import re
import subprocess
import time

import numpy as np

from affordable_common import InvalidArtifact, read_json, write_new, record, verify, cache_key, energy, xyz

PROTOCOL = 'frozen_mbis_apbs_transfer_descriptor_v1'
COULOMB_KCAL_A = 332.063713299
KJ_PER_KCAL = 4.184


def physical_boundary_key(atoms):
    # Match the actual PQR serialization, avoiding binary rounding differences
    # from PDB/OpenMM conversions while preserving the physical cavity.
    return cache_key(sorted((a['id'],*[round(v,10) for v in a['xyz_A']],round(a['radius_A'],6)) for a in atoms))


def mbis_charges(output, coordinates, total_charge):
    energy(output)  # must be a successful endpoint, never a partial charge table
    text=Path(output).read_text()
    if 'MBIS ANALYSIS' not in text:
        raise InvalidArtifact('MBIS charges unavailable; no Mulliken/formal-charge fallback')
    section=text.rsplit('MBIS ANALYSIS',1)[1]
    rows=re.findall(r'^\s*(\d+)\s+([A-Z][a-z]?)\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s*$',section,re.M)
    atoms=xyz(coordinates)
    if len(rows)!=len(atoms) or [int(r[0]) for r in rows]!=list(range(len(atoms))):
        raise InvalidArtifact('incomplete/ambiguous MBIS charge table')
    if [r[1] for r in rows]!=[a[0] for a in atoms]:
        raise InvalidArtifact('MBIS atom ordering mismatch')
    charges=[float(r[2]) for r in rows]
    if not all(math.isfinite(q) for q in charges) or abs(sum(charges)-total_charge)>1e-4:
        raise InvalidArtifact('MBIS charge closure failed; no ECP renormalization is inferred')
    return {'model':'ORCA_MBIS_monopoles', 'charge_e':charges,'sum_e':sum(charges),
            'source_output':record(output),'source_xyz':record(coordinates),
            'ecp_convention':'reported_net_atomic_charge_no_posthoc_shift',
            'electrostatic_quality_status':'requires_independent_ESP_check'}


def transfer_components(charging, direct_coulomb_kcal):
    """APBS charging terms are kJ/mol; output is kcal/mol, converted once."""
    required={'target_total','target_environment','reference_core',
              'homogeneous_total','homogeneous_environment','homogeneous_core'}
    if set(charging)!=required or not all(math.isfinite(v) for v in charging.values()) or not math.isfinite(direct_coulomb_kcal):
        raise InvalidArtifact('missing/nonfinite charging components')
    rf_total=(charging['target_total']-charging['homogeneous_total'])/KJ_PER_KCAL
    rf_env=(charging['target_environment']-charging['homogeneous_environment'])/KJ_PER_KCAL
    rf_ref=(charging['reference_core']-charging['homogeneous_core'])/KJ_PER_KCAL
    grid_cross=(charging['homogeneous_total']-charging['homogeneous_environment']-charging['homogeneous_core'])/KJ_PER_KCAL
    return {'direct_core_environment_kcal_mol':direct_coulomb_kcal,
            'reaction_field_total_kcal_mol':rf_total,'reaction_field_environment_kcal_mol':rf_env,
            'reaction_field_reference_core_kcal_mol':rf_ref,
            'reaction_field_transfer_kcal_mol':rf_total-rf_env-rf_ref,
            'homogeneous_grid_cross_kcal_mol':grid_cross,
            'grid_vs_analytic_cross_error_kcal_mol':grid_cross-direct_coulomb_kcal,
            'delta_U_kcal_mol':direct_coulomb_kcal+rf_total-rf_env-rf_ref}


def validate_state(s):
    for key in ('source','assembly','microstate','explicit_waters','boundary_mapping','charge_quality',
                'physical_atoms','core_atoms','environment_atoms','core_total_charge_e',
                'expected_environment_charge_e','settings'):
        if key not in s:
            raise InvalidArtifact(f'missing environmental state field: {key}')
    verify(s['source'])
    verify(s['boundary_mapping'])
    settings=s['settings']
    for key in ('solute_dielectric','solvent_dielectric','salt_molar','temperature_K','grid_spacing_A','grid_dimensions',
                'grid_center_A','surface','probe_radius_A','charge_discretization','boundary_condition','radii_policy'):
        if key not in settings:
            raise InvalidArtifact(f'missing solver setting: {key}')
    if settings['solute_dielectric']!=1 or settings['salt_molar']!=0:
        raise InvalidArtifact('v1 supports only dielectric-1 interior, zero salt')
    if settings['surface']!='mol' or settings['charge_discretization']!='spl2' or settings['boundary_condition']!='mdh':
        raise InvalidArtifact('unsupported numerical model')
    if len(settings['grid_dimensions'])!=3 or any(n<33 or (n-1)%32 for n in settings['grid_dimensions']):
        raise InvalidArtifact('APBS dimensions must be 32k+1')
    if any(not math.isfinite(settings[k]) or settings[k]<=0 for k in ('grid_spacing_A','solvent_dielectric','temperature_K','probe_radius_A')) or len(settings['grid_center_A'])!=3 or not all(math.isfinite(v) for v in settings['grid_center_A']):
        raise InvalidArtifact('invalid grid')
    for field in ('physical_atoms','core_atoms','environment_atoms'):
        ids=set()
        for a in s[field]:
            if a['id'] in ids or len(a['xyz_A'])!=3 or not all(math.isfinite(v) for v in a['xyz_A']) or not math.isfinite(a['radius_A']) or a['radius_A']<=0 or not math.isfinite(a.get('charge_e',float('nan'))):
                raise InvalidArtifact(f'invalid/duplicate atom in {field}')
            ids.add(a['id'])
    q=sum(a['charge_e'] for a in s['core_atoms']); env=sum(a['charge_e'] for a in s['environment_atoms'])
    if not all(math.isfinite(s[k]) for k in ('core_total_charge_e','expected_environment_charge_e')) or abs(q-s['core_total_charge_e'])>1e-4 or abs(env-s['expected_environment_charge_e'])>1e-6:
        raise InvalidArtifact('endpoint/environment charge closure failure')
    core={a['id'] for a in s['core_atoms']}; environment={a['id'] for a in s['environment_atoms']}
    if core & environment:
        raise InvalidArtifact('force-field charge on QM source atom')
    for a in s['core_atoms']:
        for b in s['environment_atoms']:
            if np.linalg.norm(np.asarray(a['xyz_A'])-b['xyz_A'])<1.0:
                raise InvalidArtifact('environment charge too close to core/cap')
    if s['charge_quality']['status']!='passed':
        raise InvalidArtifact('charge electrostatic quality is unvalidated')
    verify(s['charge_quality']['receipt'])


def pqr_lines(atoms):
    return '\n'.join(f"ATOM {i+1} X MOL 1 {a['xyz_A'][0]:.10f} {a['xyz_A'][1]:.10f} {a['xyz_A'][2]:.10f} {a['charge_e']:.10f} {a['radius_A']:.6f}"
                     for i,a in enumerate(atoms))+'\n'


def prepare(state_path, output):
    s=read_json(state_path); validate_state(s)
    output=Path(output)
    if output.exists(): raise InvalidArtifact(f'refusing existing output {output}')
    # Physical cavity atoms remain present with zero charge when represented in QM.
    physical={a['id']:dict(a,charge_e=0.0) for a in s['physical_atoms']}
    target=dict(physical)
    for a in s['environment_atoms']:
        if a['id'] not in physical or a['xyz_A']!=physical[a['id']]['xyz_A']:
            raise InvalidArtifact('unmapped environmental atom')
        target[a['id']]=a
    env=list(target.values())
    for a in s['core_atoms']:
        if a['id'] not in physical:
            # A synthetic cap may not enlarge the physical protein boundary.
            if a.get('kind')!='cap' or not any(np.linalg.norm(np.asarray(a['xyz_A'])-b['xyz_A'])+a['radius_A']<=b['radius_A']+1e-6 for b in physical.values()):
                raise InvalidArtifact('core/cofactor/cap absent from physical boundary')
        elif not np.allclose(a['xyz_A'],physical[a['id']]['xyz_A'],atol=1e-9,rtol=0) or a['radius_A']!=physical[a['id']]['radius_A']:
            raise InvalidArtifact('QM atom changes physical cavity coordinates/radius')
        target[a['id']]=a
    total=list(target.values()); core=s['core_atoms']
    # Zero-charge cap spheres are included in env too, so its cavity is exactly
    # the total cavity; the containment check above keeps the physical boundary fixed.
    env_ids={a['id'] for a in env}
    env.extend(dict(a,charge_e=0.0) for a in core if a['id'] not in env_ids)
    output.mkdir(parents=True)
    files={}
    for label,atoms in [('total',total),('environment',env),('core',core)]:
        p=output/f'{label}.pqr'; p.write_text(pqr_lines(atoms)); files[label]=record(p)
    c=s['settings']; blocks=['read',' mol pqr total.pqr',' mol pqr environment.pqr',' mol pqr core.pqr','end']
    dimensions=' '.join(str(x) for x in c['grid_dimensions'])
    center=' '.join(str(x) for x in c['grid_center_A'])
    spacing=c['grid_spacing_A']
    labels=['target_total','target_environment','reference_core','homogeneous_total','homogeneous_environment','homogeneous_core']
    for i,label in enumerate(labels):
        solvent=c['solvent_dielectric'] if i<3 else c['solute_dielectric']
        blocks.extend([f'elec name {label}',' mg-manual',f' dime {dimensions}',f' grid {spacing} {spacing} {spacing}',
                       f' gcent {center}',f' mol {i%3+1}',' lpbe',' bcfl mdh',f" pdie {c['solute_dielectric']}",
                       f' sdie {solvent}',' chgm spl2',' srfm mol',f" srad {c['probe_radius_A']}",
                       ' swin 0.3',' sdens 10.0',f" temp {c['temperature_K']}",' calcenergy total',' calcforce no','end'])
    blocks.extend(f'print elecEnergy {i+1} end' for i in range(6)); blocks.append('quit')
    inp=output/'transfer.in'; inp.write_text('\n'.join(blocks)+'\n')
    direct=0.0
    for a in core:
        for b in s['environment_atoms']:
            direct+=COULOMB_KCAL_A*a['charge_e']*b['charge_e']/np.linalg.norm(np.asarray(a['xyz_A'])-b['xyz_A'])
    manifest={'protocol_id':PROTOCOL,'state':record(state_path),'cache_key':cache_key({'state':s,'protocol':PROTOCOL,'implementation':record(__file__)}),
              'cache_scope':'preparation_only; executable identity required for result reuse',
              'implementation':record(__file__),'input':record(inp),'pqr':files,'labels':labels,
              'direct_coulomb_kcal_mol':float(direct),'status':'prepared_not_executed',
              'settings':c, 'physical_boundary_hash':physical_boundary_key(s['physical_atoms'])}
    write_new(output/'apbs_manifest.json',manifest)
    return manifest


def collect(manifest_path, output_path):
    m=read_json(manifest_path)
    verify(m['state']); verify(m['input'])
    for r in m['pqr'].values(): verify(r)
    text=Path(output_path).read_text()
    values=re.findall(r'Global net ELEC energy\s*=\s*([-+\d.eE]+)\s+kJ/mol',text)
    if len(values)!=6 or 'Thanks for using APBS' not in text:
        raise InvalidArtifact('APBS incomplete/unrecognized output; no correction available')
    charging=dict(zip(m['labels'],map(float,values)))
    components=transfer_components(charging,m['direct_coulomb_kcal_mol'])
    return {'protocol_id':PROTOCOL,'components':components,'charging_energies_kJ_mol':charging,'source_manifest':record(manifest_path),
            'physical_boundary_hash':m['physical_boundary_hash'],
            'output':record(output_path),'physical_validation_status':'not_yet_validated',
            'decision':'uncalibrated_protocol'}


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare'); q.add_argument('--state',type=Path,required=True); q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('collect'); q.add_argument('--manifest',type=Path,required=True); q.add_argument('--apbs-output',type=Path,required=True); q.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.op=='prepare': result=prepare(a.state,a.output)
    else: result=collect(a.manifest,a.apbs_output); write_new(a.output,result)
    print(result.get('status',result.get('physical_validation_status')))


if __name__=='__main__': main()
