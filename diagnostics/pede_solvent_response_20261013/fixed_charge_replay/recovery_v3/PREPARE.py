#!/usr/bin/env python3
"""Metadata-only correction of the exact XYZ-reader conversion module."""
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[4]
OLD=ROOT/'workspaces/pede_solvent_response_20261013/fixed_charge_replay_v2'
NEW=ROOT/'workspaces/pede_solvent_response_20261013/fixed_charge_replay_v3'
UP=ROOT/'workspaces/pede_solvent_response_20261013/fixed_charge_replay_v1/upstream'


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    mp=OLD/'manifest.json'
    if pin(mp)['sha256']!='d1392c78cf3c0446a73ad4337db105f9ab371d10d2932b0586023ba223820689':
        raise ValueError('parent manifest differs')
    original=json.loads(mp.read_text());m=json.loads(mp.read_text())
    failed=OLD/'execution_v1/FINAL_COLLECTION.json'
    result=json.loads(failed.read_text())
    if result['status']!='failed' or result['completed_component_evaluations']!=4:
        raise ValueError('unexpected parent execution')
    NEW.mkdir(exist_ok=False);impl=NEW/'implementation';impl.mkdir()
    for name,record in original['implementation'].items():
        source=Path(record['path'])
        if pin(source)!=record:
            raise ValueError('parent implementation changed')
        if name=='REPLAY.f90':
            source=Path(__file__).parent/name
        shutil.copyfile(source,impl/name);m['implementation'][name]=pin(impl/name)
    shutil.copyfile(__file__,impl/'PREPARE_COORDINATE_RECOVERY.py')
    m['implementation']['PREPARE_COORDINATE_RECOVERY.py']=pin(impl/'PREPARE_COORDINATE_RECOVERY.py')
    m['coordinate_reader_revision']=3
    m['technical_recovery']={'parent_manifest':pin(mp),'failed_parent_collection':pin(failed),
        'failed_parent_qualification':pin(OLD/'execution_v1/QUALIFICATION.json'),
        'fix':'input XYZ conversion now imports mctc_io_convert, matching actual xTB read_structure/read_xyz chain',
        'source_reader_chain':[pin(UP/path) for path in (
            'xtb/src/io/reader.f90','mctc-lib/src/mctc/io/read/xyz.f90',
            'mctc-lib/src/mctc/io/convert.f90','xtb/src/mctc/convert.f90')],
        'native_parameter_constants_changed':False,'tolerances_changed':False,
        'prior_inadmissible_diagonals':4,'prior_vacuum_charge_evaluations':0}
    for key in ('cases','tolerances','native_model','build','required_setup_lines','HA_TO_KCAL','fixed_contrast_definition'):
        if m[key]!=original[key]:
            raise ValueError('unexpected source/settings change: '+key)
    with (NEW/'manifest.json').open('x') as f:
        json.dump(m,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(pin(NEW/'manifest.json')))


if __name__=='__main__':
    main()
