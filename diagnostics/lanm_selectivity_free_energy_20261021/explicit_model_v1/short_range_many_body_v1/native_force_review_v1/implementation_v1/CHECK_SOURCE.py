"""Independent metadata/inverse-source review only; no imported scientific code."""
import ast
import hashlib
import json
from pathlib import Path
import re
import time

ROOT = Path(__file__).resolve().parents[4]
SRMB = ROOT / 'explicit_model_v1/short_range_many_body_v1'
PACKAGE = SRMB / 'native_force_bridge_v1'
OUT = Path(__file__).resolve().parent

def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def run():
    start, cpu = time.monotonic(), time.process_time()
    checks, seen = [], {}
    def check(name, ok):
        if not ok:
            raise AssertionError(name)
        checks.append(name)
    def verify_tree(value):
        if isinstance(value, dict):
            if set(('path', 'sha256')).issubset(value):
                actual = pin(value['path'])
                check('pin:' + value['path'], actual['sha256'] == value['sha256'])
                seen[actual['path']] = actual
            else:
                for child in value.values():
                    verify_tree(child)
        elif isinstance(value, list):
            for child in value:
                verify_tree(child)
    manifest = json.loads((PACKAGE/'generated_v1/MANIFEST.json').read_text())
    recipe = json.loads((PACKAGE/'BUILD_RECIPE.json').read_text())
    source_map = json.loads((PACKAGE/'SOURCE_MAP.json').read_text())
    for value in (manifest, recipe, source_map):
        verify_tree(value)
    for name in ('GENERATE.py', 'COMBINE.py'):
        ast.parse((PACKAGE/name).read_text())
        check('Python AST:' + name, True)
    def parent_routine(path, name):
        lines = Path(path).read_text().splitlines(keepends=True)
        first = next(i for i, line in enumerate(lines) if re.match(r'^      subroutine '+name+r'(?:\s|\()', line))
        last = next(i for i in range(first+1, len(lines)) if lines[i].strip() == 'end')
        return ''.join(lines[first:last+1])
    p = (PACKAGE/'generated_v1/srmb_epolar1a.f').read_text()
    g = (PACKAGE/'generated_v1/srmb_native_gradient.f').read_text()
    inverse_p = p.replace('      subroutine srmb_epolar1a\n', '      subroutine epolar1a\n').replace(
        'c     Supplied fixed uind/uinp; initialized udir/udirp. No SCF.\n', '      call induce\n')
    inverse_g = g.replace('      subroutine srmb_native_gradient (energy,derivs)\n', '      subroutine gradient (energy,derivs)\n').replace(
        '      if (use_polar)  call srmb_epolar1a\n', '      if (use_polar)  call epolar1\n')
    check('exact inverse polarization body', inverse_p == parent_routine(manifest['parents']['epolar1.f']['path'], 'epolar1a'))
    check('exact inverse total gradient body', inverse_g == parent_routine(manifest['parents']['gradient.f']['path'], 'gradient'))
    check('no direct induce in private clones', not re.search(r'^\s*call\s+induce\b', p+g, re.M|re.I))
    bridge = (PACKAGE/'BRIDGE.f90').read_text()
    check('no induce/moments/energy direct bridge call', not re.search(r'^\s*call\s+(induce\w*|moments|energy)\b', bridge, re.M|re.I))
    check('one direct-field callsite', len(re.findall(r'^\s*call dfield0a\(', bridge, re.M)) == 1)
    check('one mutual unit-column callsite', len(re.findall(r'^\s*call ufield0a\(', bridge, re.M)) == 1)
    check('one private gradient callsite', len(re.findall(r'^\s*call srmb_native_gradient\(', bridge, re.M)) == 1)
    imported = set(re.findall(r'^\s*use\s+(\w+)\s*$', bridge+p+g, re.M|re.I)) - {'srmb_force_state'}
    check('all imported native modules explicitly pinned', imported == set(recipe['modules']))
    commands = recipe['commands']
    check('exactly three compile commands then one link', len(commands) == 4 and all('-c' in row for row in commands[:3]) and '-c' not in commands[3])
    objects = ['{fresh_build_dir}/BRIDGE.o', '{fresh_build_dir}/srmb_epolar1a.o', '{fresh_build_dir}/srmb_native_gradient.o']
    check('explicit private objects before original archive', all(commands[3].index(obj) < commands[3].index(recipe['archive']['path']) for obj in objects))
    check('no compiler or native execution by generator receipt', manifest['compiler_calls'] == manifest['native_calls'] == manifest['numerical_tests'] == 0)
    result = {'status': 'passed_independent_native_force_source_metadata_checks',
              'checks': checks, 'check_count': len(checks), 'unique_pins_verified': len(seen),
              'verified_pins': list(seen.values()), 'implementation': pin(__file__),
              'scientific_code_imported': False, 'compiler_calls': 0, 'native_calls': 0,
              'solves': 0, 'numerical_replay': False,
              'cost': {'wall_seconds': time.monotonic()-start, 'cpu_seconds': time.process_time()-cpu}}
    with (OUT/'SOURCE_CHECKS.json').open('x') as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({k:result[k] for k in ('status','check_count','unique_pins_verified','cost')}))

if __name__ == '__main__':
    run()
