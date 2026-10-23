#!/usr/bin/env python3
"""Prepare exactly two consumed-fixture programs; never compile or invoke native."""
import hashlib
import json
from pathlib import Path
import shutil
import time


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    start,cpu=time.monotonic(),time.process_time()
    d=Path(__file__).resolve().parent
    root=Path.cwd()
    w=Path(str(d).replace('/diagnostics/','/workspaces/',1))/'preparation_v1'
    source=root/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/sampling_cpu_refactor_v1/qualification_login_preparation_v1/MANIFEST.json'
    authority=json.loads(source.read_text())
    row=authority['rows']['aquo9_central_La__near__two_water']
    assert row['atoms']==6 and row['types']==[36,37,37,36,37,37]
    for item in row['input'].values():
        assert pin(item['path'])==item,item
    old=Path(row['input']['key']['path']).read_text()
    forbidden={'a-axis','b-axis','c-axis','alpha','beta','gamma','ewald','ewald-boundary','pme-grid',
               'ewald-alpha','pewald-alpha','pme-order','ppme-order','fft-package','external-field',
               'polar-local-response','polar-donor-response','response-matrix-export','pme-norm-export','dewald'}
    assert not any(line.split()[0].lower() in forbidden for line in old.splitlines() if line.split() and not line.lstrip().startswith('#'))
    addition='''\n# Isolated initialized-PME export qualification; conducting native default.
a-axis 44.0
b-axis 44.0
c-axis 44.0
alpha 90.0
beta 90.0
gamma 90.0
ewald
ewald-cutoff 12.0
vdw-cutoff 12.0
ewald-alpha 0.3
pewald-alpha 0.3
pme-grid 32 32 32
pme-order 5
ppme-order 5
fft-package FFTW
polarization MUTUAL
'''
    w.mkdir(parents=True,exist_ok=False)
    cells={}
    for mode in ('OFF','ON'):
        out=w/mode
        out.mkdir()
        shutil.copy2(row['input']['xyz']['path'],out/'system.xyz')
        shutil.copy2(row['input']['parameters']['path'],out/'amoeba09.prm')
        (out/'system.key').write_text(old+addition+('PME-NORM-EXPORT 1\n' if mode=='ON' else ''))
        cells[mode]={k:pin(out/name) for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]}
    assert Path(cells['ON']['key']['path']).read_text()==Path(cells['OFF']['key']['path']).read_text()+'PME-NORM-EXPORT 1\n'
    manifest=dict(status='prepared_no_link_or_native_calls',implementation=pin(__file__),specification=pin(d/'SPEC.md'),
                  source_fixture_manifest=pin(source),source_fixture_id=row['id'],source_input=row['input'],
                  source_fixture_charge=row['charge'],harness=pin(d/'PME_HARNESS.f90'),snapshot_reader=pin(d.parent/'SNAPSHOT.py'),
                  cells=cells,native_program_cap=2,gradient_call_cap=2,internal_reciprocal_applications=54,
                  internal_displaced_energies=0,workers=1,root_release=None,build_result=None,harness_executable=None,
                  numerical_tolerances=dict(norm_replay_relative_to_max_one=1e-12,operator_relative_to_max_one=1e-11),
                  cost=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,threads=1),
                  native_calls=0,link_invocations=0)
    (w/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(pin(w/'MANIFEST.json')))


if __name__=='__main__':main()
