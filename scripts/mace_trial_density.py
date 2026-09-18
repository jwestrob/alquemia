"""Audit saved field-responsive trial densities for a new, fixed GK functional."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,energy,read_json,record,verify,write_new,xyz
from density_embedding import parse_potential
from mace_density_gk_hybrid import validate as validate_parent
from mace_responsive_charges import validate as validate_charges

PROTOCOL='saved_responsive_trial_density_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1'
PINS=dict(quantum='6e36ccf621e20172b07eb4673c9ede04047c2fe19f614d4b79080b7ec0729975',
          charges='fe1ab39ef2ef2ef45f43b49c57a5bce23561fa25ee6c060dfbca43b27df72145',
          prior_report='a8f0a304a3786803f3fd98e39dba40ecee93c6e7cc16eb95176d977818a2b126')
TOL=dict(intrinsic_variational_kcal=.05,algebra_kcal=1e-7,coordinate_A=1e-11)


def audit(parent,quantum,charges,prior_report,plan,output):
    start=time.monotonic();cpu=time.process_time()
    for name,path in dict(quantum=quantum,charges=charges,prior_report=prior_report).items():
        if record(path)['sha256']!=PINS[name]:raise InvalidArtifact('declared responsive archive differs: '+name)
    p=read_json(parent);pm=validate_parent(verify(p['manifest']))
    if not p['complete'] or not p['numerical_pass']:raise InvalidArtifact('numerically qualified parent required')
    q=read_json(quantum);cr=read_json(charges);old=read_json(prior_report)
    validate_charges(verify(cr['manifest']))
    cm=read_json(verify(cr['manifest']));qm=read_json(verify(q['manifest']))
    if cm['quantum']!=record(quantum) or cr['status']!='complete' or not cr['projection_gate_pass']:
        raise InvalidArtifact('responsive charges have incompatible or failed quantum source')
    terms=read_json(verify(pm['terms']));bm=read_json(verify(pm['sources']['boundary_manifest']));rows={}
    for ct in cm['tasks']:
        tid=ct['task_id'];case=ct['case_id'];metal=ct['metal'];row=cr['rows'][tid];qr=q['rows'][tid]
        qt=next(t for t in qm['tasks'] if t['task_id']==tid)
        bt=next(t for t in bm['tasks'] if t['case_id']==case and t['metal']==metal and t['mode']=='source')
        meta=read_json(verify(bt['boundary']));proj=read_json(verify(ct['projection']))
        if xyz(verify(ct['xyz']))!=xyz(verify(bt['source_QM_xyz'])) or ct['state']!=meta['state']:
            raise InvalidArtifact('responsive nuclei or physical state differ')
        if proj['physical_ids']!=read_json(verify(bt['source_projection']))['physical_ids']:
            raise InvalidArtifact('responsive source support differs')
        if ct['charge']!=bt['QM_formal_charge'] or qt['multiplicity']!=1:
            raise InvalidArtifact('responsive electronic state differs')
        if ct['source_receipt']!=qr['receipt'] or ct['source_output']!=qr['output']:
            raise InvalidArtifact('responsive output receipt mismatch')
        for pin in ct['files'].values():verify(pin)
        for pin in ct['source_wavefunctions'].values():verify(pin)
        raw=energy(verify(qr['output']))
        if raw!=qr['energy_hartree']:raise InvalidArtifact('raw embedded energy mismatch')
        group=read_json(verify(ct['probe_groups']));points=np.loadtxt(verify(ct['points']),skiprows=1)
        phi=parse_potential(verify(row['potential']),points)
        indices=np.array(group['environment_indices']);weights=read_json(verify(group['environment_weights']))
        pc=np.loadtxt(verify(qt['pointcharges']),skiprows=1)
        if pc.shape!=(len(indices),4) or not np.array_equal(pc[:,0],weights['weights_e']):
            raise InvalidArtifact('old generating field charge/order changed')
        if not np.allclose(pc[:,1:],points[indices]*BOHR_TO_A,rtol=0,atol=TOL['coordinate_A']):
            raise InvalidArtifact('old generating field coordinates changed')
        direct=HA_TO_KCAL*math.fsum(pc[:,0]*phi[indices]);intrinsic=raw-direct/HA_TO_KCAL
        prior=old['cases'][case]['endpoints'][metal]
        errors=dict(old_field=direct-row['direct_coupling_kcal_mol']['exact'],
            intrinsic=intrinsic*HA_TO_KCAL-prior['intrinsic_core_kcal_mol'])
        if max(abs(v) for v in errors.values())>TOL['algebra_kcal']:
            raise InvalidArtifact('intrinsic energy subtraction does not reproduce archive')
        cost=(intrinsic-terms[tid]['DFT_vacuum_hartree'])*HA_TO_KCAL
        rows[tid]=dict(case_id=case,metal=metal,source_task=ct,quantum_task=qt,
            embedded_energy_hartree=raw,old_field_coupling_kcal=direct,intrinsic_core_hartree=intrinsic,
            vacuum_core_hartree=terms[tid]['DFT_vacuum_hartree'],intrinsic_polarization_cost_kcal=cost,
            variational_pass=cost>=-TOL['intrinsic_variational_kcal'],algebra_errors_kcal=errors,
            core_coordinates_exact=True,physical_state_exact=True,wavefunctions_available=True,
            short_context_kcal=terms[tid]['short_context_kcal'],short_receipts=terms[tid],
            projected_charge_e=row['projected_charge_e'],source_charge_receipt=row['execution_receipt'],
            existing_center_potential=row['potential'],existing_center_points=ct['points'],
            existing_environment_indices=group['environment_indices'],boundary=bt['boundary'])
    complete=len(rows)==8
    result=dict(protocol=PROTOCOL,parent=record(parent),parent_manifest=p['manifest'],quantum=record(quantum),
        charges=record(charges),prior_report=record(prior_report),plan=record(plan),implementation=record(__file__),
        complete=complete,gates_pass=complete and all(r['variational_pass'] for r in rows.values()),rows=rows,
        tolerances=TOL,new_DFT_calls=0,new_charge_fits=0,new_potential_utilities=0,new_native_solves=0,
        numerical_score=None,calibrated_class=None,baseline_changed=False,
        receipt=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu))
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    a=s.add_parser('audit')
    for name in ('parent','quantum','charges','prior-report','plan','output'):a.add_argument('--'+name,required=True)
    a=p.parse_args();r=audit(a.parent,a.quantum,a.charges,a.prior_report,a.plan,a.output)
    print(json.dumps({k:r[k] for k in ('complete','gates_pass','new_DFT_calls','new_potential_utilities')}))


if __name__=='__main__':main()
