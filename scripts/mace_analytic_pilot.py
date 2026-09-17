"""Versioned analytic-multipole pilot using the existing MACE executor."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_hybrid import TASK_IDS,EV_TO_KCAL,rotation
from mace_rotation import rotate_vectors

SCHEMA='alquemia.mace_analytic.v1'
PROTOCOL='mace_polar_1m_analytic_multipole_vacuum_r2scan3c_pilot_v1'
LARGE_PROTOCOL='mace_polar_1l_analytic_multipole_vacuum_r2scan3c_pilot_v1'


def prepare(source_collection,rotation_manifest,kernel_tests,agreement,output,target_software=None,model_variant='medium'):
    from mace_hybrid import dry_run
    source=read_json(source_collection);sm=read_json(verify(source['manifest']))
    rotated=read_json(rotation_manifest)
    dry_run(verify(source['manifest']));dry_run(rotation_manifest)
    if source['status']!='complete':raise InvalidArtifact('completed finite-displacement reference required')
    m=copy.deepcopy(sm);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={}
    for name in ('mace_hybrid.py','mace_rotation.py','mace_analytic.py','mace_analytic_pilot.py',
                 'mace_blocked.py','mace_local_memory.py','mace_realspace_compat.py','affordable_common.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name);pins[name]=record(impl/name)
    tasks=[copy.deepcopy(next(t for t in sm['tasks'] if t['task_id']==name)) for name in TASK_IDS]
    for name in TASK_IDS:
        t=copy.deepcopy(next(t for t in rotated['tasks'] if t['task_id']==name))
        t['task_id']=name+'_rotate';tasks.append(t)
    for variant in ('primary','rotate'):
        tasks.extend(copy.deepcopy(t) for t in sm['tasks'] if t['kind']=='full' and t['variant']==variant)
    if model_variant not in ('medium','large'):
        raise InvalidArtifact('unknown checkpoint variant')
    if model_variant=='large':
        if not target_software or source['protocol_id']!=PROTOCOL:
            raise InvalidArtifact('large comparison requires analytic medium reference and pinned target software')
        target=read_json(target_software)
        m['software']=record(target_software);m['model']['checkpoint']=target['checkpoint']
    elif target_software:
        raise InvalidArtifact('unexpected checkpoint override for medium pilot')
    m.update(schema_version=SCHEMA,protocol_id=LARGE_PROTOCOL if model_variant=='large' else PROTOCOL,
             model_variant=model_variant,implementation=pins,tasks=tasks,
             agreement=record(agreement),finite_reference=record(source_collection),rotation_manifest=record(rotation_manifest),
             kernel_test_receipt=record(kernel_tests),run_inventory={'distinct_mace_energy_force_calls':12,'new_DFT_endpoints':0})
    m.pop('kernel_validation_reference',None);m.pop('kernel_validation_tolerances',None)
    m['model']['pair_kernel']={'kernel_id':'graph044_analytic_multipoles_blocked_v1','tile_atoms':256,
                              'edge_tile':2048 if model_variant=='large' else 4096,
                              'node_tile':128 if model_variant=='large' else 256,'core_edge_tile':128,'core_node_tile':17}
    for t in tasks:
        t.pop('cache_key');t['cache_key']=cache_key({'task':t,'model':m['model'],'software':m['software'],'implementation':pins})
    write_new(out/'manifest.json',m);return dry_run(out/'manifest.json')


def validate(m):
    source=read_json(verify(m['finite_reference']));sm=read_json(verify(source['manifest']))
    rotated=read_json(verify(m['rotation_manifest']))
    tests=read_json(verify(m['kernel_test_receipt']))
    for key in ('kernel','test_source','log','checkpoint','source_collection'):verify(tests[key])
    if tests.get('base_test_source'):verify(tests['base_test_source'])
    if tests['status']!='pass' or tests['kernel']['sha256']!=m['implementation']['mace_analytic.py']['sha256']:
        raise InvalidArtifact('kernel verification is missing or stale')
    variant=m.get('model_variant','medium')
    protocol={'medium':PROTOCOL,'large':LARGE_PROTOCOL}.get(variant)
    if tests['checkpoint']!=m['model']['checkpoint'] or m['protocol_id']!=protocol:
        raise InvalidArtifact('kernel test/checkpoint or protocol mismatch')
    old=copy.deepcopy(sm['model']);new=copy.deepcopy(m['model'])
    old.pop('pair_kernel');new.pop('pair_kernel')
    if variant=='large':
        if source['protocol_id']!=PROTOCOL:raise InvalidArtifact('large requires analytic medium comparator')
        old.pop('checkpoint');new.pop('checkpoint')
        a,b=read_json(verify(sm['software'])),read_json(verify(m['software']))
        if any(a[k]!=b[k] for k in ('python','requirements','backend_source_inventory','package_pins','graph_git_commit')) or b['checkpoint']!=m['model']['checkpoint']:
            raise InvalidArtifact('large pilot changes software beyond checkpoint')
    elif m['software']!=sm['software']:
        raise InvalidArtifact('medium pilot changes software')
    if old!=new or m['tolerances']!=sm['tolerances']:
        raise InvalidArtifact('unapproved physical model/software/tolerance change')
    settings=m['model']['pair_kernel']
    if settings['kernel_id']!='graph044_analytic_multipoles_blocked_v1' or settings['tile_atoms']<1:
        raise InvalidArtifact('unsupported analytic kernel')
    names=TASK_IDS+[n+'_rotate' for n in TASK_IDS]+[f'1h4i_full_{metal}_{v}' for v in ('primary','rotate') for metal in ('La','Ca')]
    if [t['task_id'] for t in m['tasks']]!=names:
        raise InvalidArtifact('analytic task inventory changed')
    for t in m['tasks']:
        if t['kind']=='core' and t['variant']=='rotate':
            expected=copy.deepcopy(next(x for x in rotated['tasks'] if x['task_id']==t['task_id'].removesuffix('_rotate')))
            expected['task_id']=t['task_id']
        else:expected=next(x for x in sm['tasks'] if x['task_id']==t['task_id'])
        if {k:v for k,v in t.items() if k!='cache_key'}!={k:v for k,v in expected.items() if k!='cache_key'}:
            raise InvalidArtifact('analytic pilot changes source state/coordinates')


def add_analysis(c,m):
    rows=c['rows'];core_checks=[];pairs=[];changes=[]
    for name in TASK_IDS:
        a,b=rows[name],rows[name+'_rotate']
        if a['status']!='computed' or b['status']!='computed':continue
        de=(b['energy_eV']-a['energy_eV'])*EV_TO_KCAL
        df=float(np.max(np.abs(np.load(verify(b['forces']))@rotation()-np.load(verify(a['forces'])))))
        dd=float(np.max(np.abs(np.load(verify(b['density_coefficients']))-rotate_vectors(np.load(verify(a['density_coefficients'])),rotation()))))
        core_checks.append({'task_id':name,'energy_delta_kcal_mol':de,'force_max_delta_eV_A':df,
                            'density_covariance_max_abs':dd,'pass':abs(de)<=m['tolerances']['energy_kcal_mol'] and df<=m['tolerances']['force_max_eV_A']})
    for size in (33,36):
        la,ca=f'1h4i_qm{size}_La',f'1h4i_qm{size}_Ca'
        if all(rows[n]['status']=='computed' for n in (la,ca,la+'_rotate',ca+'_rotate')):
            delta=((rows[ca+'_rotate']['energy_eV']-rows[la+'_rotate']['energy_eV'])-(rows[ca]['energy_eV']-rows[la]['energy_eV']))*EV_TO_KCAL
            pairs.append({'partition':f'qm{size}','paired_R_delta_kcal_mol':delta,'pass':abs(delta)<=m['tolerances']['energy_kcal_mol']})
    old=read_json(verify(m['finite_reference']))
    for name in TASK_IDS+['1h4i_full_La_primary','1h4i_full_Ca_primary']:
        a,b=old['rows'][name],rows[name]
        if b['status']!='computed':continue
        energy_key='large_minus_medium_energy_kcal_mol' if m.get('model_variant')=='large' else 'analytic_minus_finite_energy_kcal_mol'
        changes.append({'task_id':name,energy_key:(b['energy_eV']-a['energy_eV'])*EV_TO_KCAL,
                        'force_max_change_eV_A':float(np.max(np.abs(np.load(verify(b['forces']))-np.load(verify(a['forces']))))),
                        'density_max_change':float(np.max(np.abs(np.load(verify(b['density_coefficients']))-np.load(verify(a['density_coefficients'])))))})
    c['analytic_comparison']={'core_rotation_checks':core_checks,'core_paired_rotation_checks':pairs,
       'core_gate':len(core_checks)==4 and len(pairs)==2 and all(r['pass'] for r in core_checks+pairs),
       'method_changes':changes,'finite_reference':m['finite_reference']}
    if m.get('model_variant')=='large':
        c['analytic_comparison']['comparison_label']='analytic_large_minus_analytic_medium'
        c['analytic_comparison']['reference_field_note']='finite_reference is the legacy schema key; this pinned reference uses analytic medium'
    if c.get('core_partition'):
        key='partition_change_vs_medium_kcal_mol' if m.get('model_variant')=='large' else 'partition_change_vs_finite_kcal_mol'
        c['analytic_comparison'][key]=c['core_partition']['hybrid_partition_shift_kcal_mol']-old['core_partition']['hybrid_partition_shift_kcal_mol']


def core_gate(manifest):
    from mace_hybrid import collect
    c=collect(manifest)
    return {'status':'pass' if c['analytic_comparison']['core_gate'] else 'failed',
            'core_checks':c['analytic_comparison']['core_rotation_checks']}


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare')
    for key in ('source-collection','rotation-manifest','kernel-tests','agreement','output'):q.add_argument('--'+key,required=True)
    q.add_argument('--target-software');q.add_argument('--model-variant',choices=('medium','large'),default='medium')
    q=sub.add_parser('core-gate');q.add_argument('--manifest',required=True);q.add_argument('--output')
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.source_collection,a.rotation_manifest,a.kernel_tests,a.agreement,a.output,a.target_software,a.model_variant)
    else:
        r=core_gate(a.manifest)
        if a.output:write_new(a.output,r)
    print(json.dumps(r,indent=2))
    if r['status']=='failed':raise SystemExit(1)

if __name__=='__main__':main()
