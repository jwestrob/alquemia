"""Produce an isolated explicit-operator initializer; never edit the library."""
import difflib
import hashlib
import json
from pathlib import Path


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def derive(source):
    start = source.index('      subroutine periodicresponsesetup (handled)\n')
    stop = source.index('      subroutine periodicresponsebackendguard\n', start)
    old = source[start:stop]
    assert old.count('      call periodicresponsebackendguard\n') == 1
    new = old.replace('subroutine periodicresponsesetup (handled)',
                      'subroutine parityresponsesetup (handled)', 1)
    new = new.replace('c     This source has no admission of a real periodic execution.\n'
                      '      call periodicresponsebackendguard\n',
                      'c     Standalone explicit-operator initialization only.\n'
                      'c     Normal setup and both fatal guards remain unchanged.\n')
    assert new != old and 'call periodicresponsebackendguard' not in new
    return old, new


def main(authorities, out):
    auth = json.loads(Path(authorities).read_text())
    manifest_pin = auth['parent_source']
    assert pin(manifest_pin['path']) == manifest_pin
    manifest = json.loads(Path(manifest_pin['path']).read_text())
    source_pin = next(x for x in manifest['source_outputs']
                      if x['path'].endswith('/source/induce.f'))
    assert pin(source_pin['path']) == source_pin
    old, new = derive(Path(source_pin['path']).read_text())
    out = Path(out)
    for name, text in [('PARITY_SETUP.f', new), ('SETUP_AUTHORITY.f', old),
                       ('SETUP.diff', ''.join(difflib.unified_diff(
                           old.splitlines(True), new.splitlines(True),
                           fromfile='immutable/periodicresponsesetup',
                           tofile='isolated/parityresponsesetup')))]:
        with (out/name).open('x') as stream:
            stream.write(text)
    record = dict(status='derived_operator_only_setup_source',
                  implementation=pin(__file__), parent_manifest=manifest_pin,
                  parent_induce=source_pin, setup=pin(out/'PARITY_SETUP.f'),
                  authority=pin(out/'SETUP_AUTHORITY.f'), diff=pin(out/'SETUP.diff'),
                  normal_source_modified=False, native_calls=0, builds=0)
    (out/'SETUP_PROVENANCE.json').write_text(json.dumps(record, indent=2)+'\n')


if __name__ == '__main__':
    import sys
    main(*sys.argv[1:])
