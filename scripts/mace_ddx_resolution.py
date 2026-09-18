"""Separate basis and integration errors on an archived real conductor pair."""
from pathlib import Path
import argparse
import json
import shutil
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_ddx_cpcm import preflight as parent_preflight
from mace_ddx_recovery import execute_group

PROTOCOL='fixed_source_ddCPCM_basis_quadrature_separation_v1'
GRIDS=dict(integration18=dict(lmax=18,n_lebedev=2030),integration24=dict(lmax=24,n_lebedev=3470),
           basis30=dict(lmax=30,n_lebedev=3470),integration30=dict(lmax=30,n_lebedev=5810))


def prepare(parent,collection,plan,output):
    p=parent_preflight(parent); c=read_json(collection)
    if p['protocol']!='fixed_source_full_protein_ddCPCM_component_refinement_v2' or c['manifest']!=record(parent) or not c['complete']:
        raise InvalidArtifact('matching actual refined conductor collection required')
    root=Path(output).resolve(); root.mkdir(parents=True,exist_ok=False); impl=root/'implementation'; impl.mkdir()
    for name,pin in p['implementation'].items(): shutil.copyfile(verify(pin),impl/name)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    m=dict(protocol=PROTOCOL,parent=record(parent),collection=record(collection),plan=record(plan),
           implementation={f.name:record(f) for f in impl.glob('*.py')},
           **{k:p[k] for k in ('model','energy_prefactor','python','module','tolerances')},grids=GRIDS,
           requested_new_forward_solves=8,new_DFT_calls=0,new_MACE_calls=0,groups=[])
    source=next(g['source'] for g in p['groups'] if g['group_id']=='GGR_2FW0_primary_primary')
    for grid in GRIDS:
        g=dict(group_id='GGR_2FW0_'+grid,case_id='GGR_2FW0',variant='primary',grid=grid,
               source=source,states=['Ca','La'],new_states=['Ca','La'],reuse={})
        g['cache_key']=cache_key(dict(group=g,model=m['model'],grid=GRIDS[grid],implementation=m['implementation'],plan=m['plan']))
        m['groups'].append(g)
    m['cache_key']=cache_key(m);write_new(root/'manifest.json',m);return preflight(root/'manifest.json')


def preflight(path):
    m=read_json(path); p=parent_preflight(verify(m['parent'])); c=read_json(verify(m['collection']))
    if m['cache_key']!=cache_key({k:v for k,v in m.items() if k!='cache_key'}) or m['protocol']!=PROTOCOL or m['grids']!=GRIDS:
        raise InvalidArtifact('resolution model or cache differs')
    for k in ('model','energy_prefactor','python','module','tolerances'):
        if m[k]!=p[k]: raise InvalidArtifact('physical/numerical method changed: '+k)
    for pin in (m['plan'],*m['implementation'].values()): verify(pin)
    if c['manifest']!=m['parent'] or not c['complete'] or len(m['groups'])!=4 or m['requested_new_forward_solves']!=8:
        raise InvalidArtifact('resolution experiment inventory differs')
    source=next(g['source'] for g in p['groups'] if g['group_id']=='GGR_2FW0_primary_primary')
    for g,grid in zip(m['groups'],GRIDS):
        expected=dict(group_id='GGR_2FW0_'+grid,case_id='GGR_2FW0',variant='primary',grid=grid,
                      source=source,states=['Ca','La'],new_states=['Ca','La'],reuse={})
        expected['cache_key']=cache_key(dict(group=expected,model=m['model'],grid=GRIDS[grid],implementation=m['implementation'],plan=m['plan']))
        if g!=expected: raise InvalidArtifact('resolution source/group changed')
    return m


def collect(path,output):
    m=preflight(path); root=Path(path).resolve().parent; groups={}; checks=[]; prior=read_json(verify(m['collection']))
    for g in m['groups']:
        pin=record(root/'groups'/g['group_id']/'attempt_0001/result.json'); r=read_json(verify(pin))
        if r['cache_key']!=g['cache_key'] or r['manifest']!=record(path): raise InvalidArtifact('resolution execution mismatch')
        verify(r['arrays']);groups[g['grid']]=dict(receipt=pin,**r)
        checks.append(dict(name=g['grid']+'_reciprocity',error=r['reciprocity_error_kcal'],tolerance=.05,
                           pass_=r['reciprocity_error_kcal'] is not None and r['reciprocity_error_kcal']<=.05))
        for state,row in r['rows'].items():
            expected={**m['model'],**GRIDS[g['grid']]}
            parameters=all(row.get('actual_parameters',{}).get(k)==v for k,v in expected.items() if k not in ('model','enable_force'))
            ok=row['status']=='computed' and parameters and row.get('contraction_pass') and row.get('passive_energy_pass')
            checks.append(dict(name=g['grid']+'_'+state+'_native_state',pass_=bool(ok)))
    comparisons={}
    for name,new,old in (
        ('integration18',groups['integration18'],prior['groups']['GGR_2FW0_primary_primary']),
        ('integration24',groups['integration24'],prior['groups']['GGR_2FW0_primary_refined']),
        ('basis24_to_30',groups['basis30'],groups['integration24']),
        ('integration30',groups['integration30'],groups['basis30'])):
        values={}
        for metal in ('Ca','La'):
            a=new['rows'][metal].get('energy_kcal_mol');b=old['rows'][metal].get('energy_kcal_mol')
            values[metal]=None if a is None or b is None else a-b
        values['R']=None if None in values.values() else values['Ca']-values['La']
        for term,value in values.items(): checks.append(dict(name=name+'_'+term,error=value,tolerance=.1,pass_=value is not None and abs(value)<=.1))
        comparisons[name]=values
    r=dict(protocol=PROTOCOL,manifest=record(path),groups=groups,comparisons_kcal=comparisons,checks=checks,
           complete=all(g['status']=='complete' for g in groups.values()),checks_pass=all(c['pass_'] for c in checks),
           full_model_qualified=False,new_score=None,new_DFT_calls=0,new_MACE_calls=0,baseline_changed=False)
    write_new(output,r);return r


def execute(path,output):
    m=preflight(path);root=Path(path).resolve().parent
    for g in m['groups']:execute_group(m,path,g['group_id'],root/('group_'+g['group_id']+'.json'))
    return collect(path,output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare')
    for k in ('parent','collection','plan','output'):p.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','collect'):
        p=sub.add_parser(op);p.add_argument('--manifest',required=True)
        if op!='dry-run':p.add_argument('--output',required=True)
    a=parser.parse_args()
    if a.command=='prepare':r=prepare(a.parent,a.collection,a.plan,a.output)
    elif a.command=='dry-run':r=preflight(a.manifest)
    elif a.command=='execute':r=execute(a.manifest,a.output)
    else:r=collect(a.manifest,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('implementation','groups')}))
