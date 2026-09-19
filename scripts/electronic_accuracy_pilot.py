#!/usr/bin/env python3
"""Prepare a finite real-core correlated electronic diagnostic with existing runners."""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import shutil
import re

from affordable_common import InvalidArtifact, record, verify, read_json, write_new, xyz, energy

PROTOCOL = 'real_core_dlpno_ccsdt1_tightpno_cpcm_ptes_tz_v1'


def electronic_input(charge, metal):
    basis = '  NewGTO Ca "cc-pwCVTZ" end\n' if metal == 'Ca' else '  NewECP La "def2-ECP" end\n'
    return ('! DLPNO-CCSD(T1) TightPNO def2-TZVPPD def2/J RIJCOSX TightSCF DefGrid3 CPCM(Water) NoAutostart\n'
        '%maxcore 3000\n%basis\n  AuxC "AutoAux"\n' + basis + 'end\n'
        '%method\n  NewNCore Ca 10 end\n  NewNCore La 46 end\nend\n'
        '%cpcm\n  CPCMccm 2\nend\n' + f'* xyzfile {charge} 1 core.xyz\n')


def sources(root):
    hp=root/'workspaces/hydration_network_20260918/core_transfer_v1/manifest.json'
    hm=read_json(hp)
    hc=read_json(hp.parent/'collection_1201867.json')
    bp=root/'workspaces/baseline_benchmark_20260915/run_v1/manifest.json'
    bm=read_json(bp)
    rows=[]
    for t in bm['tasks']:
        if t['case']=='ggr_1glg_GGR' and t['lane']=='repaired':
            e=hm['ggr_reused']['endpoints'][t['metal']]
            rows.append(dict(case='GGR_1GLG',metal=t['metal'],xyz=t['xyz'],charge=t['charge'],
                source_manifest=t['source_manifest'],baseline=e,source_input=t['input'],source_tasks=record(bp)))
    for t in hm['tasks']:
        if t['case']!='1F6S':continue
        metal=t['metal']
        old=t['original_endpoints'][metal]
        # The archived output receipt records its immutable original input.
        rp=read_json(verify(old['receipt']))
        for suffix,xp,endpoint in [('original',t['original_xyz'],old),
                ('water_prepared',t['xyz'],next(r['result'] for r in hc['rows'] if r['case']=='1F6S' and r['metal']==metal))]:
            rows.append(dict(case='ALPHA_1F6S_'+suffix,metal=metal,xyz=xp,charge=t['charge'],
                source_manifest=t['parent'],baseline=endpoint,
                source_input=t['input'] if suffix=='water_prepared' else rp['artifacts']['template_input'],source_tasks=record(hp)))
    return rows,hm


def prepare(root, agreement, output):
    root=Path(root).resolve();out=Path(output).resolve()
    rows,hm=sources(root)
    if len(rows)!=6:raise InvalidArtifact('six source endpoints required')
    # Validate all original artifacts before writing any task.
    for r in rows:
        atoms=xyz(verify(r['xyz']))
        if sum(a[0] in ('Ca','La') for a in atoms)!=1 or atoms[0][0]!=r['metal']:
            raise InvalidArtifact('unexpected metal/atom ordering')
        measured=energy(verify(r['baseline']['output']))
        if abs(measured-r['baseline']['energy_hartree'])>1e-9:raise InvalidArtifact('baseline energy mismatch')
        verify(r['source_manifest']);verify(r['source_input'])
    out.mkdir(parents=True,exist_ok=False)
    pp=out/'preparation.json'
    write_new(pp,dict(protocol_id=PROTOCOL,agreement=record(agreement),sources=rows,
        baseline_changed=False,implementation=record(__file__)))
    tasks=[]
    for r in rows:
        name=r['case']+'_'+r['metal'];td=out/name;td.mkdir()
        xp=td/'core.xyz';shutil.copyfile(verify(r['xyz']),xp)
        ip=td/'endpoint.inp';ip.write_text(electronic_input(r['charge'],r['metal']))
        tasks.append(dict(task_id=name,case=r['case'],metal=r['metal'],charge=r['charge'],multiplicity=1,
            input=record(ip),xyz=record(xp),output_path=str(td/'endpoint.out'),source_xyz=r['xyz'],
            source_input=r['source_input'],source_manifest=r['source_manifest'],task_type='single_point',
            protocol_id=PROTOCOL))
    mp=out/'manifest.json'
    write_new(mp,dict(schema_version='alquemia.electronic_accuracy.v1',protocol_id=PROTOCOL,
        agreement=record(agreement),preparation=record(pp),orca=hm['orca'],
        execution_policy=hm['execution_policy'],execution_resources=dict(mpi_ranks=16,concurrent_tasks=4),
        tasks=tasks,baseline_changed=False))
    from affordable_workflow import dry_run
    return dry_run(mp)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','agreement','output'):p.add_argument('--'+name,required=True)
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
