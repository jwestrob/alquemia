#!/usr/bin/env python3
"""Create an immutable instrumented baseline source copy; no compiler/native."""
import difflib
import hashlib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/tinker_source'
OUT = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/source_v1'


def pin(path):
    return {'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def prepare():
    assert not OUT.exists(), 'Preserve source versions'
    helper = (HERE/'EXPORT_HELPERS.f').read_text()
    assert all(len(s)<=72 for s in helper.splitlines() if s and s[0] not in ('c','C','*')), 'Fixed-form code longer than72 columns'
    original = (BASE/'source/induce.f').read_text()
    start = original.index('      subroutine induce\n'); end = original.index('      subroutine induce0a\n')
    first = original[start:end]
    assert first.count('      logical header\n') == 1
    first = first.replace('      logical header\n','      logical header,matrix_export\n')
    marker = 'c     choose the method for computation of induced dipoles\nc\n'
    assert first.count(marker) == 1
    first = first.replace(marker,marker+'      call responsematrixflag (matrix_export)\nc\n')
    marker = 'c     update the lists of previous induced dipole values\nc\n'
    assert first.count(marker) == 1
    first = first.replace(marker,'c     Reporting hook preserves original state before continuation.\nc\n      if (matrix_export)  call exportresponsematrix\nc\n'+marker)
    changed = original[:start]+first+original[end:]+helper
    OUT.mkdir(parents=True)
    shutil.copytree(BASE/'source',OUT/'source'); shutil.copytree(BASE/'cmake',OUT/'cmake')
    (OUT/'source/induce.f').write_text(changed)
    diff = ''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),'a/source/induce.f','b/source/induce.f'))
    (OUT/'EXPORT.patch').write_text(diff)
    changes = [p.name for p in (OUT/'source').iterdir() if p.is_file() and p.read_bytes() != (BASE/'source'/p.name).read_bytes()]
    assert changes == ['induce.f'],changes
    result = {'status':'baseline_export_source_prepared_uncompiled','implementation':pin(Path(__file__)),
              'helper':pin(HERE/'EXPORT_HELPERS.f'),'source_root':str(OUT),'baseline_root':str(BASE),'changed_files':changes,
              'baseline_inputs':[pin(p) for p in sorted((BASE/'source').iterdir()) if p.is_file()],
              'source_outputs':[pin(p) for p in sorted((OUT/'source').iterdir()) if p.is_file()],
              'cmake':pin(OUT/'cmake/CMakeLists.txt'),'patch':pin(OUT/'EXPORT.patch'),
              'source_checks':{'only_induce_changed':True,'helper_fixed_form_within72':True,'original_induce0a_and_all_kernels_unchanged':changed[changed.index('      subroutine induce0a\n'):].startswith(original[end:])},
              'compiler_calls':0,'new_native_calls':0,'physics_changed':False}
    (OUT/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'manifest':pin(OUT/'MANIFEST.json'),'patch':result['patch']}))


if __name__ == '__main__':
    prepare()
