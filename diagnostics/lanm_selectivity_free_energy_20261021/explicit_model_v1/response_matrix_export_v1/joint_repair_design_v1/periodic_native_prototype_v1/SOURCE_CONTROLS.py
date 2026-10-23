#!/usr/bin/env python3
"""Narrow source joins and artificial expression controls; never invokes native."""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import re
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
OUT = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/periodic_native_prototype_v1'


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def routines(text):
    return {m[1]: m[0] for m in re.finditer(
        r'(?m)^      subroutine (\w+)[^\n]*\n.*?^      end\n', text, re.S)}


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def native_expression_kernel(r):
    # Literal scalar operation order of the new Fortran helper, not its execution.
    if r >= 21:
        return 0., 0.
    fzero = 1 / (((1+5+25/3)*math.exp(-5))**2)
    t = 2*r; p = 1+t+t*t/3
    f = fzero*(p*math.exp(-t))**2
    fp = -(4/3)*f*t*(1+t)/p
    if r <= 20:
        return f, fp
    u = r-20
    sw = (1-u)**3*(1+3*u+6*u*u)
    swp = -30*u*u*(1-u)**2
    return f*sw, fp*sw+f*swp


def main():
    tick = time.monotonic(); cpu = time.process_time(); checks = []
    def check(name, ok, details=None):
        checks.append({'name': name, 'passed': bool(ok), 'details': details})
        if not ok:
            raise AssertionError(name)
    manifest_path = OUT/'source_v1/MANIFEST.json'
    manifest = json.loads(manifest_path.read_text())
    parent = json.loads(Path(manifest['parent_manifest']['path']).read_text())
    new = Path(manifest['source_root'])/'source'
    old = Path(parent['source_root'])/'source'
    directory = OUT/'source_controls_v1'
    directory.mkdir(exist_ok=False)
    report = {'status': 'started', 'implementation': pin(__file__),
              'source_manifest': pin(manifest_path), 'checks': checks,
              'new_native_calls': 0, 'builds': 0, 'actual_box_evaluations': 0}
    try:
        check('all_parent_source_pins_exact',
              all(pin(p['path']) == p for p in parent['source_outputs']))
        check('all_prepared_source_pins_exact',
              all(pin(p['path']) == p for p in manifest['source_outputs']))
        check('all_declared_input_pins_exact',
              all(pin(p['path']) == p for p in manifest['inputs'].values()))
        changed = sorted(p.name for p in old.iterdir() if p.is_file()
                         and p.read_bytes() != (new/p.name).read_bytes())
        check('only_four_declared_source_files_changed',
              changed == ['epolar1.f', 'epolar3.f', 'induce.f', 'polar.f'])
        ni = (new/'induce.f').read_text(); oi = (old/'induce.f').read_text()
        nr, oldr = routines(ni), routines(oi)
        check('unchanged_original_induce_routines', all(
            nr[k] == v for k, v in oldr.items()
            if k not in ('localresponsesetup', 'localresponsegradient', 'induce0a')),
            {'count': len(oldr)-3})
        prelude = ('      lrperiodic = .false.\n'
                   '      call periodicresponsesetup(periodic_handled)\n'
                   '      if (periodic_handled) return\n')
        restored = nr['localresponsesetup'].replace(prelude, '').replace(
            'sites(28),periodic_handled', 'sites(28)')
        check('vacuum_setup_exact_after_namespace_predispatch',
              restored == oldr['localresponsesetup'])
        dispatch = ('      if (lrperiodic) then\n'
                    '         call periodicresponsegradient\n'
                    '         return\n      end if\n')
        check('vacuum_gradient_exact_after_periodic_dispatch',
              nr['localresponsegradient'].replace(dispatch, '') ==
              oldr['localresponsegradient'])
        certificate = ('c     The source-only periodic path has no implemented stability backend.\n'
                       '      if (lrperiodic) call periodicresponsestability\n')
        check('original_solver_exact_except_explicit_periodic_certificate',
              nr['induce0a'].replace(certificate, '') == oldr['induce0a'])
        check('certificate_after_direct_preparation_before_mutual_iteration',
              nr['induce0a'].index('call dfield0 (field,fieldp)') <
              nr['induce0a'].index('call periodicresponsestability') <
              nr['induce0a'].index('udir(j,i) = polarity(i) * field(j,i)'))
        for name in ('periodicresponsebackendguard', 'periodicresponsestability'):
            r = nr[name]
            check(name+'_unconditional_fatal', 'call fatal' in r and
                  not re.search(r'(?m)^\s+if\b', r) and 'UNIMPLEMENTED' in r)
        setup = nr['periodicresponsesetup']
        check('zero_bypass_precedes_native_modes_and_backend_guard',
              setup.index('if (.not.lractive) return') <
              setup.index('bad=.not.use_bounds') <
              setup.index('call periodicresponsebackendguard'))
        check('backend_guard_after_geometry', setup.index('call periodicresponsegeometry') <
              setup.index('call periodicresponsebackendguard'))
        for token in ('ia.ne.900.and.ia.ne.901', 'ib.ne.402', 'type(i).eq.403',
                      'type(i).eq.405', 'type(i).eq.414', 'nwat.gt.3980',
                      'nna.ne.7 .or. ncl.ne.10', 'min(xbox,ybox,zbox).lt.44.0d0',
                      "boundary.ne.'TINFOIL'", '.not.orthogonal',
                      'npole.ne.n', 'polarity(i).lt.1.0d-8', 'nuse.ne.n',
                      'abs(c)+abs(d).gt.0.5d0*dom+1.0d-15',
                      "'RESPONSE-MATRIX-EXPORT'", "'POLAR-DONOR-RESPONSE'",
                      'i12(1,hj).ne.i', 'all(i12(1:2,oj).ne.i)',
                      'i12(1,i).eq.i12(2,i)'):
            check('guard_'+token, token in setup)
        pg = nr['periodicresponsegeometry']; dg = nr['periodicresponsegradient']
        check('periodic_work_arrays_dynamic', 'allocatable' in pg and
              'allocatable' in dg and '(28)' not in pg+dg and '(0:28)' not in pg+dg)
        check('nearest_image_and_compact_kernel_once',
              pg.count('vec(:,i)-sides*anint(vec(:,i)/sides)') == 1 and
              pg.count('call periodicresponsekernel(r,ff(i),fr(i))') == 1)
        check('switch_product_rule_before_field_multiply',
              nr['periodicresponsekernel'].index('fr=fr*sw+f*swr') <
              nr['periodicresponsekernel'].index('f=f*sw'))
        check('both_common_field_hooks_unchanged', nr['dfield0'] == oldr['dfield0'] and
              nr['ufield0'] == oldr['ufield0'])
        check('all_native_frame_and_moment_sources_identical', all(
            (new/f).read_bytes() == (old/f).read_bytes()
            for f in ('empole1.f', 'torque.f', 'moments.f', 'epolar.f')))
        ep1n = routines((new/'epolar1.f').read_text()); ep1o = routines((old/'epolar1.f').read_text())
        check('all_gradient_subkernels_unchanged_except_old_hook_removal', all(
            ep1n[k] == v.replace('      if (lractive) call localresponsegradient\n', '')
            for k, v in ep1o.items() if k != 'epolar1'))
        check('single_common_gradient_hook',
              (new/'epolar1.f').read_text().count('call localresponsegradient') == 1 and
              'call localresponsegradient' in ep1n['epolar1'])
        ep3n = routines((new/'epolar3.f').read_text()); ep3o = routines((old/'epolar3.f').read_text())
        correction = ep3n['localresponseanalysis'].split('c     Reduced stationary')[1]
        correction = 'c     Reduced stationary'+correction[:correction.index('      return\n')]
        check('analysis_subkernels_unchanged_except_old_hook_removal', all(
            ep3n[k] == v.replace(correction, '') for k, v in ep3o.items() if k != 'epolar3'))
        check('single_common_pairwise_only_analysis_hook',
              (new/'epolar3.f').read_text().count('call localresponseanalysis') == 1 and
              'if (lractive.and.pairwise)' in ep3n['epolar3'])
        added = [line for line in Path(manifest['patch']['path']).read_text().splitlines()
                 if line.startswith('+') and not line.startswith('+++')]
        long = [line for line in added if len(line[1:]) > 72 and not line[1:].startswith('c')]
        check('added_fixed_form_statements_fit_72_columns', not long, long)
        src = Path(manifest['inputs']['compact_implementation']['path'])
        spec = importlib.util.spec_from_file_location('compact_frozen', src)
        compact = importlib.util.module_from_spec(spec); spec.loader.exec_module(compact)
        for radius in (0., 1e-6, 2.5, 8., 20., 20.000001, 20.25, 20.5, 20.75, 20.999999, 21., 22.):
            a = native_expression_kernel(radius); b = compact.radial(radius)[:2]
            check('radial_expression_'+str(radius), all(
                abs(x-y) <= 2e-14*max(abs(x), abs(y), 1e-300)
                for x, y in zip(a, b)), {'Fortran_expression': a, 'prototype': b})
        # Unequal dual channels and one switch-region water exercise both passes.
        vectors = [[2.4, .3, -.2], [-1.7, 2.2, .4], [20.5, .1, -.2]]
        coef = dict(a=.1, b=-.08, c=.1, d=-.03, kappa=.7, alpha=1.1)
        response = compact.response(vectors, **coef)
        u = [[.2, -.1, .3], [-.3, .4, .2], [.1, .6, -.2]]
        v = [[-.1, .3, .4], [.2, -.2, .5], [.4, -.1, .3]]
        um, vm = [.4, -.2, .1], [-.1, .3, .2]
        cross = dot(um, vm); Khalf = .5*332.0637133
        zi=[]; term=[]; ph=[]; pr=[]; radii=[]
        for vec, uu, vv, f, q in zip(vectors, u, v, response['f'], response['q']):
            r=math.sqrt(dot(vec, vec)); radii.append(r)
            V, JV=compact.regularized_v(vec, coef['kappa'])
            phi=dot(V, vec)/dot(vec, vec); ph.append(phi)
            phip_over_r=(dot(vec, [sum(JV[k][ell]*vec[ell] for ell in range(3))
                                  for k in range(3)])/dot(vec, vec)-phi)/dot(vec, vec)
            pr.append(phip_over_r)
            z=3*dot(V, [x+y for x,y in zip(uu,vv)]);zi.append(z)
            term.append(coef['b']*f*z/(1+q)**2)
        prefix=[0.]
        for value in term: prefix.append(prefix[-1]+value)
        suffix=[0.]*(len(term)+1)
        for i in range(len(term)-1,-1,-1): suffix[i]=suffix[i+1]+term[i]
        gradients=[]
        for j,vec in enumerate(vectors):
            w=[x+y for x,y in zip(u[j],v[j])]
            pref=coef['a']+coef['b']*response['q'][j]/(1+response['q'][j])
            ci=pref*zi[j]+prefix[j]+suffix[j+1]
            native=[Khalf*(cross*response['dbeta'][j][k]-ci*response['df'][j][k]
                    -3*response['f'][j]*pref*(ph[j]*w[k]+pr[j]*dot(vec,w)*vec[k]))
                    for k in range(3)]
            direct=[Khalf*(cross*response['dbeta'][j][k]-sum(
                    (u[i][l]+v[i][l])*response['dh'][i][j][l][k]
                    for i in range(3) for l in range(3))) for k in range(3)]
            gradients.append(native)
            check('unequal_dual_two_pass_contraction_water_'+str(j),
                  all(abs(x-y) <= 2e-12*max(1.,abs(y)) for x,y in zip(native,direct)),
                  {'two_pass': native, 'direct_Jacobian': direct})
        reaction=[-sum(g[k] for g in gradients) for k in range(3)]
        check('equal_opposite_metal_reaction', all(
              abs(sum(g[k] for g in gradients)+reaction[k]) == 0 for k in range(3)))
        strain=[[.1,.2,-.1],[.2,-.3,.05],[-.1,.05,.2]]
        vir=[[sum(.5*(r[k]*g[l]+r[l]*g[k]) for r,g in zip(vectors,gradients))
              for l in range(3)] for k in range(3)]
        contracted=sum(vir[k][l]*strain[k][l] for k in range(3) for l in range(3))
        explicit=sum(dot(g,[dot(row,r) for row in strain]) for r,g in zip(vectors,gradients))
        check('symmetric_image_virial_matches_affine_chain',
              abs(contracted-explicit) <= 2e-12*max(1.,abs(explicit)))
        report['status']='passed_periodic_prototype_source_controls'
        report['limitations']=['No compiler or Fortran execution',
            'Artificial expression checks are not native numerical qualification',
            'No PME backend or scalable stability admission',
            'No actual candidate or box evaluated']
    except Exception as exc:
        report['status']='failed_periodic_prototype_source_controls'
        report['error']=repr(exc)
    report['cost']={'wall_seconds':time.monotonic()-tick,
                    'CPU_seconds':time.process_time()-cpu,'threads':1}
    report['check_count']=len(checks)
    (directory/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'result':pin(directory/'RESULT.json'),'status':report['status'],
                      'checks':len(checks),'cost':report['cost']}))
    if report['status'].startswith('failed'):
        raise SystemExit(1)


if __name__=='__main__': main()
