"""Common physical GK source/reference boundary for declared real inputs."""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path
import shutil
import time
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_density_frameworks import validate as frameworks
from mace_density_gk_boundary import boundary,additional_key,verify_native,comparable,METALS,TOL
from mace_tinker_framework_solver import native_call

PROTOCOL='source_graph_common_Ca2018_GK_boundary_v1'


def key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol','mode','frameworks','preparation','software','charges','implementation','tolerances')}))


def prepare(framework_manifest,software,output,charges=None):
    start=time.monotonic();cpu=time.process_time();fm=frameworks(framework_manifest)
    p=read_json(verify(fm['preparation']));sw=read_json(software)
    if sw['returncode'] or sw['library']!=read_json(verify(fm['software']))['library']:
        raise InvalidArtifact('qualified common-boundary backend required')
    cr=read_json(charges) if charges else None
    if cr is not None and (cr['status']!='complete' or not cr['projection_gate_pass']):
        raise InvalidArtifact('actual qualified projected density charges required')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    paths={k:verify(v) for k,v in fm['implementation'].items()}
    for n in ('mace_density_boundary_panel.py','mace_density_gk_boundary.py','mace_density_multipoles.py',
              'mace_native_field_accounting.py','mace_qm_field.py','mace_density_frameworks.py'):
        paths[n]=Path(__file__).with_name(n)
    for n,src in paths.items():shutil.copyfile(src,impl/n)
    m=dict(protocol=PROTOCOL,mode='source' if charges else 'zero_source',frameworks=record(framework_manifest),
        preparation=fm['preparation'],software=record(software),charges=record(charges) if charges else None,
        implementation={p.name:record(p) for p in impl.glob('*.py')},tolerances=TOL,tasks=[],cases={},
        new_energy_calls=0,new_response_solves=0,baseline_changed=False)
    for case,pin in p['cases'].items():
        c=read_json(verify(pin));state=read_json(verify(c['state']))
        old=next(t for t in fm['tasks'] if t['case_id']==case)
        mapping=read_json(verify(old['mapping']));base=read_json(verify(old['parameters']))
        meta=boundary(state,mapping,base);physical={a['id']:a for a in state['physical_atoms']};ids=meta['physical_ids']
        meta.update(state=c['state'],base_parameters=old['parameters'],base_mapping=old['mapping'],
            positions_A=[physical[i]['xyz_A'] for i in ids],assembly=state['assembly'],microstate=state['microstate'],
            explicit_waters=state['explicit_waters'],evidence=state['evidence'],evidence_use=state['evidence_use'])
        d=root/'cases'/case;d.mkdir(parents=True);write_new(d/'boundary.json',meta);m['cases'][case]=record(d/'boundary.json')
        template=verify(old['xyz']).read_text().splitlines()
        for metal in ('Ca','La'):
            formal=c['endpoints'][metal]['charge'];srcq={};source_receipt=None
            if cr is not None:
                row=cr['rows'][case+'_'+metal];proj=read_json(verify(row['projection']))
                if set(proj['physical_ids'])!=set(meta['source_support_ids']):raise InvalidArtifact('source support differs')
                if any(pos!=physical[i]['xyz_A'] for i,pos in zip(proj['physical_ids'],proj['coordinates_A'])):
                    raise InvalidArtifact('source projected coordinates differ')
                actual=math.fsum(row['projected_charge_e'])
                if abs(actual-formal)>TOL['QM_charge_e']:raise InvalidArtifact('source charge closure fails')
                srcq=dict(zip(proj['physical_ids'],row['projected_charge_e']));source_receipt=row['execution_receipt']
            tid=case+'_'+metal+'_'+m['mode'];d=root/'tasks'/tid;d.mkdir(parents=True);n=len(ids)
            q=[meta['environment_charge_by_id'][i]+srcq.get(i,0.) for i in ids]
            lines=[f'{n} {tid}; unchanged common physical GK boundary',*template[1:]]
            lines.append(f'{n} {metal} '+' '.join(format(v,'.17g') for v in physical['metal']['xyz_A'])+f' {METALS[metal]["type"]}')
            (d/'framework.xyz').write_text('\n'.join(lines)+'\n')
            opts=['parameters '+str(verify(fm['native_parameters'])),'multipoleterm ONLY','polarizeterm','solvateterm',
                  'solvate GK','dielectric 1.0','gk-radius SOLUTE','gkc 2.455','polarization MUTUAL','polar-iter 100','polar-eps 1e-7']
            (d/'framework.key').write_text('\n'.join(opts+additional_key())+'\n')
            frozen=[i+1 for i,pid in enumerate(ids) if pid in meta['source_support_ids']]
            (d/'framework.freeze').write_text(str(len(frozen))+'\n'+'\n'.join(map(str,frozen))+'\n')
            (d/'framework.charges').write_text(str(n)+'\n'+''.join(f'{i} {v:.17g}\n' for i,v in enumerate(q,1)))
            t=dict(task_id=tid,case_id=case,metal=metal,mode=m['mode'],boundary=m['cases'][case],
                charge_e=q,frozen_indices=frozen,QM_formal_charge=formal,QM_actual_charge=math.fsum(srcq.values()) if cr else None,
                source_charge_execution=source_receipt,source_QM_xyz=c['endpoints'][metal]['xyz'],
                xyz=record(d/'framework.xyz'),key=record(d/'framework.key'),mask=record(d/'framework.freeze'),overrides=record(d/'framework.charges'))
            t['cache_key']=key(t,m);receipt=native_call(sw['executable'],t,str(verify(t['overrides'])),d/'native.log',1)
            r=dict(status='failed',task_id=tid,cache_key=t['cache_key'],receipt=receipt)
            try:
                if receipt['returncode']:raise InvalidArtifact('native boundary initialization failed')
                r.update(verify_native(verify(receipt['log']).read_text(),t,meta,base))
            except Exception as exc:r['failure_reason']=str(exc)
            write_new(d/'result.json',r);m['tasks'].append(t)
            if r['status']=='failed':
                write_new(root/'failure.json',dict(result=record(d/'result.json'),attempted_initializations=len(m['tasks'])))
                raise InvalidArtifact(r['failure_reason'])
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m);return collect(root/'manifest.json',root/'result.json')


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['tolerances']!=TOL or m['mode'] not in ('source','zero_source'):
        raise InvalidArtifact('boundary configuration differs')
    for p in [m['frameworks'],m['preparation'],m['software'],*m['implementation'].values(),*m['cases'].values()]:verify(p)
    if m['charges']:verify(m['charges'])
    if len(m['tasks'])!=2*len(m['cases']):raise InvalidArtifact('paired boundary inventory differs')
    for t in m['tasks']:
        if key(t,m)!=t['cache_key']:raise InvalidArtifact('boundary cache differs')
        for k in ('boundary','source_QM_xyz','xyz','key','mask','overrides'):verify(t[k])
        if t['source_charge_execution']:verify(t['source_charge_execution'])
    return m


def collect(path,output):
    m=validate(path);root=Path(path).resolve().parent;rows={};checks=[];results={}
    for t in m['tasks']:
        rp=root/'tasks'/t['task_id']/'result.json';r=read_json(rp);meta=read_json(verify(t['boundary']))
        if r['cache_key']!=t['cache_key'] or r['status']!='native_density_boundary_prepared':raise InvalidArtifact('boundary result unavailable')
        verify_native(verify(r['receipt']['log']).read_text(),t,meta,read_json(verify(meta['base_parameters'])))
        rows[t['task_id']]=dict(result=record(rp),status=r['status']);results[t['task_id']]=r
    for case in m['cases']:
        a,b=[results[case+'_'+e+'_'+m['mode']] for e in ('Ca','La')]
        if m['mode']=='zero_source':checks.append(dict(name=case+'_environment_identity',pass_=comparable(a)==comparable(b)))
        else:
            x,y=[copy.deepcopy(r['parameters']) for r in (a,b)]
            for p in (x,y):
                p['atoms'][-1]['atomic_number']=0
                for atom in p['atoms']:atom['charge_e']=0.
            checks.append(dict(name=case+'_paired_cavity',pass_=x==y))
    r=dict(protocol=PROTOCOL,mode=m['mode'],manifest=record(path),cases=m['cases'],tasks=rows,checks=checks,
        complete=True,gates_pass=all(c['pass_'] for c in checks),actual_native_initializations=len(rows),
        native_wall_seconds=sum(r['receipt']['wall_seconds'] for r in results.values()),
        native_CPU_seconds=sum(r['receipt']['child_CPU_seconds'] for r in results.values()),
        new_energy_calls=0,new_response_solves=0,numerical_score=None,baseline_changed=False)
    write_new(output,r);return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare')
    for n in ('frameworks','software','output'):s.add_argument('--'+n,required=True)
    s.add_argument('--charges')
    s=sub.add_parser('audit');s.add_argument('--manifest',required=True)
    a=p.parse_args()
    if a.command=='prepare':
        r=prepare(a.frameworks,a.software,a.output,a.charges);print(json.dumps({k:r[k] for k in ('mode','complete','gates_pass','actual_native_initializations')}))
    else:print(json.dumps(dict(status='pass',tasks=len(validate(a.manifest)['tasks']))))
