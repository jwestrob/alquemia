#!/usr/bin/env python3
"""Guard-only private source revision. No compiler or numerical calls."""
import copy, difflib, hashlib, json, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.with_name('native_force_bridge_v1')
ROOT = next(p for p in HERE.parents if (p/'workspaces').is_dir())
WORK = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
PARENT = WORK/'response_matrix_export_v1/source_v2/source'
BUILD = WORK/'short_range_many_body_v1/native_force_bridge_v1/build_v1'

def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p): return json.loads(Path(p).read_text())
def save(p,x): Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def statements(p):
    return [s.strip().lower() for s in Path(p).read_text().splitlines()
            if s.strip() and s[0] not in 'cC*!']

def main():
    t=time.monotonic();cpu=time.process_time()
    old=read(OLD/'READY.json');recipe=read(OLD/'BUILD_RECIPE.json')
    assert pin(OLD/'BRIDGE.f90')==old['sources']['BRIDGE.f90']
    assert recipe['archive']['sha256']=='37135742acaf53548f472966b1aa940d22fdb2ab72dcb4913eb9bb2524e838e2'
    assert pin(recipe['archive']['path'])==recipe['archive']
    expected_extra=['subroutine extra1','use atoms','use deriv','use energi','implicit none','integer i',
                    'ex = 0.0d0','do i = 1, n','dex(1,i) = 0.0d0','dex(2,i) = 0.0d0',
                    'dex(3,i) = 0.0d0','end do','return','end']
    assert statements(PARENT/'extra1.f')==expected_extra
    assert statements(PARENT/'kextra.f')==['subroutine kextra','implicit none','return','end']
    assert pin(PARENT/'extra1.f')['sha256']=='28ee5bf85df7749f85d31fac38a9e462c460628eacb2289fcd949c503c1bdef8'
    assert '      use_extra = .true.' in (PARENT/'field.f').read_text()
    assert '      if (use_extra)  call extra1' in (PARENT/'gradient.f').read_text()
    source=(OLD/'BRIDGE.f90').read_text()
    first="    call insist(n>=1.and.n<=28.and.npole==n,'finite 1..28 all-atom response scope')"
    observe="""    ! Preserve all native switches; report them even if a later guard rejects.
    write(*,'(a,5(1x,l1))') 'SRMB_GUARD_CHGTRN_RXNFLD_ORBIT_EXTRA_METAL', &
                          use_chgtrn,use_rxnfld,use_orbit,use_extra,use_metal
    write(*,'(a)') 'SRMB_EXTRA_POLICY SOURCE_BOUND_ZERO_ONLY_EXTRA1_28ee5bf8'
    flush(6)
"""
    old_guard='use_chgtrn.or.use_rxnfld.or.use_orbit.or.use_extra.or.use_metal'
    new_guard='use_chgtrn.or.use_rxnfld.or.use_orbit.or.use_metal'
    rationale="""    ! field.f defaults use_extra=T. Exact linked extra1.f only zeros ex/dex;
    ! kextra.f is empty. The source/build contract pins this archive and noop.
    ! Keep use_extra and its original gradient dispatch, never turn physics off.
"""
    after='    gradient_calls=gradient_calls+1\n'
    post="""    call insist(ex==0d0.and.all(dex(:,1:n)==0d0),'source-bound no-op extra1 changed energy or gradient')
    write(*,'(a)') 'SRMB_EXTRA_ZERO_ENERGY_GRADIENT T'
    flush(6)
"""
    assert source.count(first)==source.count(old_guard)==source.count(after)==1
    new=source.replace(first,observe+first).replace('    call insist(.not.('+old_guard, rationale+'    call insist(.not.('+new_guard).replace(after,after+post)
    restored=new.replace(observe,'').replace(rationale,'').replace(new_guard,old_guard).replace(post,'')
    assert restored==source
    (HERE/'BRIDGE.f90').write_text(new)
    (HERE/'BRIDGE.diff').write_text(''.join(difflib.unified_diff(source.splitlines(True),new.splitlines(True),fromfile='native_force_bridge_v1/BRIDGE.f90',tofile='native_force_bridge_v2/BRIDGE.f90')))
    r=copy.deepcopy(recipe)
    r['commands']=[copy.deepcopy(recipe['commands'][0]),copy.deepcopy(recipe['commands'][-1])]
    r['commands'][0]=[str(HERE/'BRIDGE.f90') if a==str(OLD/'BRIDGE.f90') else a for a in r['commands'][0]]
    reused=[]
    for name in ('srmb_epolar1a.o','srmb_native_gradient.o'):
        obj=pin(BUILD/'build'/name);reused.append(obj)
        r['commands'][1]=[obj['path'] if a=='{fresh_build_dir}/'+name else a for a in r['commands'][1]]
    r.update(status='source_only_guard_recovery_build_recipe',scope='one BRIDGE compile and one link; reuse two immutable original private objects; no native invocation',
             source_files=[pin(HERE/'BRIDGE.f90')],reused_objects=reused,reused_object_build=pin(BUILD/'RESULT.json'),
             inherited_recipe=pin(OLD/'BUILD_RECIPE.json'),noop_sources={n:pin(PARENT/n) for n in ('extra1.f','kextra.f')},
             max_build_workers=1,actual_compiler_calls=0)
    save(HERE/'BUILD_RECIPE.json',r)
    evidence={'status':'source_guard_recovery_prepared','generator':pin(__file__),'parent_ready':pin(OLD/'READY.json'),
              'failed_scout':pin(WORK/'short_range_many_body_v1/native_force_qualification_v1/scout_v1/execution/RESULT.json'),
              'native_authorities':{n:pin(PARENT/n) for n in ('field.f','kextra.f','extra1.f','kchgtrn.f','mechanic.f','gradient.f')},
              'bridge':pin(HERE/'BRIDGE.f90'),'diff':pin(HERE/'BRIDGE.diff'),'recipe':pin(HERE/'BUILD_RECIPE.json'),
              'checks':{'parent_pin':True,'extra1_exact_zero_only_executable_body':True,'kextra_exact_empty_body':True,
                        'original_extra_dispatch_retained':True,'all_other_bridge_text_exact_after_delta_removal':True,
                        'native_switches_never_assigned_by_delta':True},
              'native_calls':0,'compiler_calls':0,'numerical_tests':0,
              'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    save(HERE/'SOURCE_CHECKS.json',evidence)
    template=copy.deepcopy(old)
    template.update(status='source_ready_for_native_force_qualification',root_build=None,parent_ready=pin(OLD/'READY.json'),
                    guard_recovery=pin(HERE/'SOURCE_CHECKS.json'),private_objects_reused=True,source_generator_ran_once=True)
    template['sources']['BRIDGE.f90']=pin(HERE/'BRIDGE.f90')
    template['sources']['BUILD_RECIPE.json']=pin(HERE/'BUILD_RECIPE.json')
    template['sources']['PREPARE_SOURCE.py']=pin(__file__)
    template['sources']['BRIDGE.diff']=pin(HERE/'BRIDGE.diff')
    save(HERE/'RUNTIME_READY_TEMPLATE.json',template)
    ready={'status':'guard_only_recovery_source_ready_no_build_or_native','bridge':pin(HERE/'BRIDGE.f90'),
           'recipe':pin(HERE/'BUILD_RECIPE.json'),'source_checks':pin(HERE/'SOURCE_CHECKS.json'),
           'runtime_ready_template':pin(HERE/'RUNTIME_READY_TEMPLATE.json'),'parent_ready':pin(OLD/'READY.json'),
           'phase':old['phase'],'kernel':old['kernel'],'combiner':old['sources']['COMBINE.py'],
           'native_calls':0,'compiler_calls':0}
    save(HERE/'READY.json',ready)
    print(json.dumps({'ready':pin(HERE/'READY.json'),'bridge':ready['bridge'],'recipe':ready['recipe'],'cost':evidence['cost']},indent=2))

if __name__=='__main__':main()
