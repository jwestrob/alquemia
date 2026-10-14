#!/usr/bin/env python3
"""Metadata-only successor for failed CMake integration; no scientific execution."""
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[4]
OLD = ROOT/'workspaces/pede_solvent_response_20261013/fixed_charge_replay_v1'
NEW = ROOT/'workspaces/pede_solvent_response_20261013/fixed_charge_replay_v2'


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    mp=OLD/'manifest.json'
    if pin(mp)['sha256']!='2469b38c0d227df8806aa48763608c50672eb1660a14532679c65185c214b5ca':
        raise ValueError('parent manifest differs')
    original=json.loads(mp.read_text())
    m=json.loads(mp.read_text())
    failed=OLD/'execution_v1/FINAL_COLLECTION.json'
    result=json.loads(failed.read_text())
    if result['status']!='failed' or result['completed_component_evaluations']!=0:
        raise ValueError('parent is not the declared zero-evaluation failure')
    NEW.mkdir(exist_ok=False)
    impl=NEW/'implementation';impl.mkdir()
    for name,record in original['implementation'].items():
        source=Path(record['path'])
        if pin(source)!=record:
            raise ValueError('parent implementation changed')
        if name=='CMakeLists.txt':
            source=Path(__file__).parent/name
        shutil.copyfile(source,impl/name)
        m['implementation'][name]=pin(impl/name)
    shutil.copyfile(__file__,impl/'PREPARE_RECOVERY.py')
    m['implementation']['PREPARE_RECOVERY.py']=pin(impl/'PREPARE_RECOVERY.py')
    m['build_integration_revision']=2
    m['technical_recovery']={'parent_manifest':pin(mp),'failed_parent_collection':pin(failed),
        'failed_configure_stderr':pin(OLD/'execution_v1/configure.err'),
        'fix':'replace mctc-lib target ALIAS with upstream imported INTERFACE integration',
        'upstream_pattern':pin(OLD/'upstream/xtb/cmake/modules/xtb-utils.cmake'),
        'scientific_settings_changed':False,'native_evaluations_repeated':0}
    for key in ('cases','tolerances','native_model','build','required_setup_lines','HA_TO_KCAL','fixed_contrast_definition'):
        if m[key]!=original[key]:
            raise ValueError('unexpected scientific/build input change: '+key)
    with (NEW/'manifest.json').open('x') as f:
        json.dump(m,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(pin(NEW/'manifest.json')))


if __name__=='__main__':
    main()
