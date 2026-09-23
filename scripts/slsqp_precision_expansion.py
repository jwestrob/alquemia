"""Fixed34-source expansion of the original qualified SLSQP precision policy."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import scipy
import slsqp_precision as pilot
import union_adaptive as original
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz,paired,HA_TO_KCAL
from mace_hybrid import EV_TO_KCAL
from mace_site_kinematics import Kinematics
from adaptive_force_diagnostic import project,preview

FAILED=('q92wy9-pqq-la_model__conditioned_Ca__seed-1_sample-2','q60ar6-pqq-la_model__conditioned_La__seed-1_sample-0')
NONCANONICAL=('a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1','a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-3',
              'a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4','a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4')


def population(reference):
    ref=read_json(reference)
    canonical=[{k:r[k] for k in ('case_id','known_class','role')} for r in ref['rows']]
    if len(canonical)!=28 or sum(r['role']=='calibration' for r in canonical)!=25:raise InvalidArtifact('canonical denominator differs')
    rows=canonical+[{'case_id':c,'known_class':'La','role':'consumed_noncanonical_development'} for c in NONCANONICAL]
    rows += [{'case_id':c,'known_class':'La' if c.startswith('q92') else 'Ca','role':'consumed_optimizer_limit_failure'} for c in FAILED]
    if len({r['case_id'] for r in rows})!=34:raise InvalidArtifact('source duplicate')
    return rows


def prepare(canonical,pilot_source,canonical_pool,pilot_pool,failed_sources,reuse,reference,agreement,output):
    if len(failed_sources)!=2:raise InvalidArtifact('two exact failed-source manifests required')
    sources=[canonical,pilot_source,*failed_sources];parents=[read_json(p) for p in sources];all_rows=population(reference)
    oldpools=[canonical_pool,pilot_pool,None,None];reused=read_json(reuse);rpm=read_json(verify(reused['manifest']));rsm=read_json(verify(rpm['source_manifest']))
    if rsm['settings']!=pilot.SETTINGS or rpm['settings']!=original.POOL_SETTINGS or reused['available']!=4 or set(rsm['declared_case_ids'])!=set(pilot.CASES):raise InvalidArtifact('exact successful precision pilot required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];cases=[];records=[]
    for row in all_rows:
        cid=row['case_id'];i=next((i for i,p in enumerate(parents) if any(c['case_id']==cid for c in p['cases'])),None)
        if i is None:raise InvalidArtifact('declared actual source absent: '+cid)
        parent=parents[i];case=next(c for c in parent['cases'] if c['case_id']==cid)
        if parent['protocol_id']!=original.PROTOCOL or parent['settings']!=pilot.original_completion.SETTINGS:raise InvalidArtifact('old protocol differs')
        receipts={z:record(Path(sources[i]).parent/'proposals'/(cid+'__'+z)/'result.json') for z in ('Ca','La')}
        statuses={z:read_json(verify(p))['status'] for z,p in receipts.items()}
        if cid in FAILED and statuses!={'Ca':'proposal_available','La':'unavailable'}:raise InvalidArtifact('declared old failure differs')
        row={**row,'parent':record(sources[i]),'old_pool':record(oldpools[i]) if oldpools[i] else None,'old_receipts':receipts,'reused_precision_pilot':cid in pilot.CASES}
        records.append(row)
        if row['reused_precision_pilot']:
            if next(c for c in rsm['cases'] if c['case_id']==cid)!=case:raise InvalidArtifact('reused source context differs')
            continue
        cases.append(case);tasks.extend(next(t for t in parent['tasks'] if (t['case_id'],t['metal'])==(cid,z)) for z in ('Ca','La'))
    p=parents[0];m={k:p[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','resources','optimizer_software')}
    m.update(protocol_id=pilot.PROTOCOL,settings=pilot.SETTINGS,gates=pilot.GATES,agreement=record(agreement),reference=record(reference),
      all_sources=records,sources=[r for r in records if not r['reused_precision_pilot']],reuse=record(reuse),cases=cases,tasks=tasks,
      declared_case_ids=[c['case_id'] for c in cases],population='numerical_expansion34_new30',declared_total=34,maximum_optimizer_starts=60,
      maximum_cross_MACE_calls=60,maximum_GFN2_calls=240,new_origin_calls=0,new_DFT_calls=0,production_changed=False,
      implementation=original.snapshot(out/'implementation'))
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest):
    m=read_json(manifest);expected_rows=population(verify(m['reference']));ids=[r['case_id'] for r in expected_rows if r['case_id'] not in pilot.CASES]
    if (m['protocol_id']!=pilot.PROTOCOL or m['settings']!=pilot.SETTINGS or m['gates']!=pilot.GATES or m['declared_case_ids']!=ids or
        m['population']!='numerical_expansion34_new30' or m['declared_total']!=34 or len(m['cases'])!=30 or len(m['tasks'])!=60 or
        m['maximum_optimizer_starts']!=60 or m['maximum_cross_MACE_calls']!=60 or m['maximum_GFN2_calls']!=240):raise InvalidArtifact('frozen expansion scope differs')
    for k in ('agreement','reference','reuse','software','orca','cpu_executable','gpu_executable'):verify(m[k])
    for p in m['implementation'].values():verify(p)
    for k in ('wrapper','kernel'):verify(m['optimizer_software'][k])
    if m['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer software changed')
    if {k for k in pilot.SETTINGS if pilot.SETTINGS[k]!=pilot.original_completion.SETTINGS[k]}!={'optimizer_ftol'}:raise InvalidArtifact('additional optimizer change')
    for row,want in zip(m['all_sources'],expected_rows):
        if any(row[k]!=v for k,v in want.items()):raise InvalidArtifact('reference/transfer membership differs')
        for z,pin in row['old_receipts'].items():
            old=read_json(verify(pin))
            if old['task_id']!=row['case_id']+'__'+z or old['manifest']!=row['parent']:raise InvalidArtifact('old proposal receipt differs')
    expected=[]
    for cid,source,c in zip(ids,m['sources'],m['cases']):
        if source['case_id']!=cid or c['case_id']!=cid:raise InvalidArtifact('source order differs')
        parent=read_json(verify(source['parent']))
        if parent['settings']!=pilot.original_completion.SETTINGS or next(x for x in parent['cases'] if x['case_id']==cid)!=c:raise InvalidArtifact('old source changed')
        if source['old_pool']:
            old=read_json(verify(source['old_pool']));om=read_json(verify(old['manifest']))
            if om['source_manifest']!=source['parent']:raise InvalidArtifact('old source pool differs')
        elif cid not in FAILED:raise InvalidArtifact('undeclared missing old pool')
        for k in ('model','software','orca','cpu_python','gpu_python','resources'):
            if m[k]!=parent[k]:raise InvalidArtifact('electronic/runtime method changed')
        data=[];pair=[]
        for z in ('Ca','La'):
            t=next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(cid,z));expected.append(t['task_id']);pair.append(t)
            if t!=next(q for q in parent['tasks'] if q['task_id']==t['task_id']):raise InvalidArtifact('task changed')
            ep=c['origin_row']['native_endpoints'][z];native,forces=original.origin(ep,m['model']);kin=Kinematics(read_json(verify(t['mapping']))['context'])
            coords,v,raw,normed,_=project(kin.data,np.zeros(len(kin.modes)),forces);data.append((kin,v,normed));point=t['origin_reuse']['point']
            if not np.array_equal(raw[t['active_indices']],point['gradient_kcal_mol_rad']) or native['energy_eV']!=t['q0']['components']['MACE_eV']:raise InvalidArtifact('q0 force or energy changed')
            if not np.allclose(coords,[a[1:] for a in xyz(verify(t['xyz']))],atol=1e-12,rtol=0):raise InvalidArtifact('q0 map differs')
            pilot.angular.final_geometry(kin,t,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
        paired(verify(pair[1]['xyz']),verify(pair[0]['xyz']),pair[1]['charge'],pair[0]['charge'])
        if read_json(verify(pair[0]['mapping']))!=read_json(verify(pair[1]['mapping'])):raise InvalidArtifact('paired mapping changed')
        choice=preview(data[0][0].modes,data[0][1],data[0][2],data[1][2])
        if choice!=c['selection'] or any(t['selector']!=choice for t in pair):raise InvalidArtifact('paired selector changed')
    if [t['task_id'] for t in m['tasks']]!=expected:raise InvalidArtifact('task order changed')
    return {'status':'validated','manifest':record(manifest),'total_sources':34,'reused_pools':4,'new_searches':60,'cross_MACE_maximum':60,'GFN2_maximum':240,'new_q0_calls':0,'new_DFT_calls':0}


def engine():
    p=pilot.private('slsqp_precision');p.validate=validate;return p.engine()
def execute(manifest):return engine()[0].execute(manifest)
def collect(manifest,output):return engine()[0].collect(manifest,output)
def prepare_pool(proposals,agreement,output):return engine()[0].prepare_pool(proposals,agreement,output)
def validate_pool(manifest):return engine()[0].validate_pool(manifest)
def execute_mace(manifest):return engine()[1].execute_mace(manifest)
def collect_pool(manifest,output):
    validate_pool(manifest);return engine()[1].collect(manifest,output)


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    specs={'prepare':('canonical','pilot_source','canonical_pool','pilot_pool','failed_sources','reuse','reference','agreement','output'),
       'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'prepare_pool':('proposals','agreement','output'),
       'validate_pool':('manifest',),'execute_mace':('manifest',),'collect_pool':('manifest','output')}
    for op,args in specs.items():
        q=sub.add_parser(op)
        for arg in args:q.add_argument('--'+arg.replace('_','-'),type=Path,required=True,**({'nargs':'+'} if arg=='failed_sources' else {}))
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
if __name__=='__main__':main()
