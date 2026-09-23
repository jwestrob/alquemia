"""Finite declared100-triple transfer: exact pool audit and existing-runner inputs."""
from __future__ import annotations
import argparse
import copy
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz,paired
from second_shell_context import parent_state
from coordination_preparation_context import geometry
from mace_site_kinematics import Kinematics
from adaptive_force_diagnostic import project,preview
from adaptive_origin_recovery import snapshot
import consistent_context as context
import compact_solvation as solvent
import nikasha_pool as pool
import slsqp_precision as precision
from slsqp_precision_transfer import qualification
import union_adaptive as union
import union_triple_preparation as prep

PROTOCOL='Nikasha_three_La_membership_all100_transfer_v1'
REFERENCE_ID='Nikasha_three_La_membership_ftol1e8_canonical25_v1'
CANDIDATES=['origin','adaptive_Ca','adaptive_La']


def identity(c):return c['case_id']+'__union_'+c['selection_id'].rsplit('__',1)[1]


def mapped(c,config):
    state=parent_state(c['original_core'],config['topology'],require_endpoint_receipts=False)
    audit=read_json(verify(c['representations']['context']['preparation']));maps={}
    for z in ('Ca','La'):
        ep=c['representations']['context']['endpoints'][z]
        maps[z]=geometry(state,audit,xyz(verify(ep['xyz'])),xyz(verify(c['original_core']['endpoints'][z]['xyz'])))
    maps=json.loads(json.dumps(maps))
    if maps['Ca']!=maps['La']:raise InvalidArtifact('paired physical maps differ')
    e=c['representations']['context']['endpoints'];paired(verify(e['La']['xyz']),verify(e['Ca']['xyz']),e['La']['charge'],e['Ca']['charge'])
    return maps


def mapping_equivalence(current,archived):
    """Existing coordinate-copy tolerance; every discrete graph/mode stays exact."""
    a=copy.deepcopy(current);b=copy.deepcopy(archived);maximum=0.
    for section in ('core','context'):
        for key in ('positions_A','core_positions_A'):
            if key not in a[section] and key not in b[section]:continue
            if key not in a[section] or key not in b[section]:raise InvalidArtifact('mapping coordinate fields differ')
            x=np.asarray(a[section].pop(key));y=np.asarray(b[section].pop(key))
            if x.shape!=y.shape or not np.allclose(x,y,rtol=0,atol=1e-12):raise InvalidArtifact('physical map coordinates differ')
            maximum=max(maximum,float(np.max(np.abs(x-y))))
    if a!=b:raise InvalidArtifact('physical graph, cap rule or mode mapping differs')
    return maximum


def audit_pool(c,config,comparison_row,reference):
    """Read actual receipts; no coordinate or numerical-profile substitution."""
    cid=c['case_id'];pin=comparison_row['precision_collection'];d=read_json(verify(pin))
    pm=read_json(verify(d['manifest']));sm=read_json(verify(pm['source_manifest']))
    row=next(r for r in d['cases'] if r['source'].get('actual_union_case_id',r['case_id'])==cid)
    old=row['source']['union']
    if (row['pool']['status']!='available' or sm['settings']!=precision.SETTINGS or sm['model']!=config['model'] or
        sm['software']!=config['software'] or pm['settings']!=reference['settings'] or
        old['state_key']!=c['state_key'] or old['original_core']!=c['original_core']):
        raise InvalidArtifact('actual pool profile/model/source/state differs: '+cid)
    maps=mapped(c,config);pins={};projection=[];low_manifests={};coordinate_delta={}
    for z in ('Ca','La'):
        t=next(t for t in sm['tasks'] if (t['case_id'],t['metal'])==(row['case_id'],z))
        ep=c['representations']['context']['endpoints'][z]
        archived_map=read_json(verify(t['mapping']))
        if not context.reusable_state(ep,t):
            raise InvalidArtifact('new context physical map or coordinates differ: '+cid)
        coordinate_delta[z]=mapping_equivalence(maps[z],archived_map)
        pins[z]=t['mapping'];native,forces=union.origin(row['source']['origin_row']['native_endpoints'][z],config['model'])
        kin=Kinematics(archived_map['context']);_,v,raw,normed,_=project(kin.data,np.zeros(len(kin.modes)),forces)
        if not np.array_equal(raw[t['active_indices']],t['origin_reuse']['point']['gradient_kcal_mol_rad']):
            raise InvalidArtifact('actual q0 force projection differs')
        projection.append((kin,v,normed,t))
        for cell in row['matrix'][z].values():
            actual_native=read_json(verify(cell['MACE']))
            if actual_native['status'] not in ('computed','complete') or actual_native['energy_eV']!=cell['components']['MACE_eV']:
                raise InvalidArtifact('actual MACE component differs')
            for medium in ('vacuum','alpb'):
                low=cell['low'][medium];mp=verify(low['manifest']);actual=solvent.completed(mp,low['task_id'])
                if actual is None or any(low[k]!=actual[k] for k in actual):raise InvalidArtifact('actual GFN2 receipt differs')
                if str(mp) not in low_manifests:low_manifests[str(mp)]=read_json(mp)
                lm=low_manifests[str(mp)];lt=next(t for t in lm.get('all_tasks',lm['tasks']) if t['task_id']==low['task_id'])
                component='GFN2_'+('vacuum' if medium=='vacuum' else 'ALPB')+'_hartree'
                if (actual['energy_hartree']!=cell['components'][component] or lt['charge']!=ep['charge'] or lt['multiplicity']!=ep['multiplicity'] or
                    not context.same_geometry(xyz(verify(lt['xyz'])),xyz(verify(cell['xyz']))) or
                    verify(lt['input']).read_text().replace(' MaxIter 500\n','')!=solvent.input_text(ep['charge'],ep['multiplicity'],medium,'native')):
                    raise InvalidArtifact('actual candidate geometry/state/recipe differs')
    choice=preview(projection[0][0].modes,projection[0][1],projection[0][2],projection[1][2])
    if any(p[3]['selector']!=choice for p in projection) or row['source']['selection']!=choice:
        raise InvalidArtifact('actual paired four-mode selector differs')
    if pool.choose_rows(row['matrix'],CANDIDATES)!=row['pool']:raise InvalidArtifact('actual pool algebra differs')
    return {'status':'exact_precision_pool_reuse','collection':pin,'source_manifest':pm['source_manifest'],'source_case_id':row['case_id'],
        'mapping':pins,'maximum_map_coordinate_copy_difference_A':coordinate_delta,'pool':row,
        'numerical_provenance':'actual original receipts retain their one/eight-rank allocation; rank1 separately qualified'}


def audit(preparation,origin_reuse,precision_comparison,reference,agreement,output,workers=4):
    start=time.monotonic();p=read_json(preparation);sel=prep.validate_selection(verify(p['selection']));r=read_json(origin_reuse)
    comparison=read_json(precision_comparison);ref=read_json(reference)
    if (p['triple_denominator'],p['prepared_triples'],p['unique_source_selection_pairs'])!=(100,94,125):raise InvalidArtifact('frozen100/94/125 population differs')
    if r['preparation']!=record(preparation):raise InvalidArtifact('origin audit preparation differs')
    if (ref['reference_id']!=REFERENCE_ID or ref['optimizer_settings']!=precision.SETTINGS or ref['canonical_denominator']!=25 or
        ref['crystals_or_stress_used_for_fit'] or any(v['status']!='available' for v in ref['variants'].values())):
        raise InvalidArtifact('complete independent threefold canonical25 reference required')
    out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces/')
    out.mkdir(parents=True,exist_ok=False)
    cmp={x['case_id']:x for x in comparison['rows']};orig={(x['selection_id'],x['case_id']):x for x in r['rows']}
    cases=[c for c in p['cases'] if not c['is_separate_stress_probe']];rows=[];futures={}
    if len(cases)!=125 or len({identity(c) for c in cases})!=125:raise InvalidArtifact('exact125 unique primary pairs required')
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for c in cases:
            if c['status']!='prepared':raise InvalidArtifact('unexpected new preparation failure')
            if c['tenfold_comparison']['exact_reuse_eligible']:
                futures[executor.submit(audit_pool,c,sel['config'],cmp[c['case_id']],ref)]=c
        results={}
        for future in as_completed(futures):
            c=futures[future];results[identity(c)]=future.result()
            print(json.dumps({'pair_id':identity(c),'status':'exact_precision_pool_reuse'}),flush=True)
    for c in cases:
        pid=identity(c);reuse=orig[c['selection_id'],c['case_id']]
        row={'pair_id':pid,'case_id':c['case_id'],'selection_id':c['selection_id'],'source':c,'origin_reuse':reuse,
            'pool_reuse':results.get(pid),'status':'exact_precision_pool_reuse' if pid in results else 'new_pool_required'}
        if pid not in results:
            maps=mapped(c,sel['config']);pins={}
            for z in ('Ca','La'):
                dest=out/'maps'/(pid+'__'+z+'.json');write_new(dest,maps[z]);pins[z]=record(dest)
            row['mapping']=pins
        rows.append(row)
    counts={'triple_denominator':100,'complete_prepared_triples':94,'old_unavailable_triples':6,'source_membership_pairs':125,
        'full_pool_reuses':len(results),'new_source_pools':len(rows)-len(results),
        'new_origin_MACE':sum(len(x['origin_reuse']['missing_native_endpoints']) for x in rows if not x['pool_reuse']),
        'new_origin_GFN2':sum(len(x['origin_reuse']['missing_solvent_cells']) for x in rows if not x['pool_reuse'])}
    counts.update(new_searches=2*counts['new_source_pools'],maximum_cross_MACE=2*counts['new_source_pools'],maximum_candidate_GFN2=8*counts['new_source_pools'])
    if (counts['full_pool_reuses'],counts['new_source_pools'],counts['new_origin_MACE'],counts['new_origin_GFN2'])!=(70,55,38,76):
        raise InvalidArtifact('expected finite work changed; record and review before execution: '+json.dumps(counts))
    result={'protocol_id':PROTOCOL,'preparation':record(preparation),'origin_reuse':record(origin_reuse),
        'precision_comparison':record(precision_comparison),'reference':record(reference),'agreement':record(agreement),
        'settings':precision.SETTINGS,'pool_settings':ref['settings'],'config':sel['config'],'triples':p['triples'],'rows':rows,'counts':counts,
        'original_context_proposal_reuse':'unavailable: archived original-context adaptive_completion profile uses old optimizer_ftol; only exact origins reused',
        'original_context_profile':record(Path(__file__).parent/'adaptive_completion.py'),
        'frozen_UTC':datetime.now(timezone.utc).isoformat(),'new_molecular_calls':0,'production_changed':False,
        'wall_seconds':time.monotonic()-start,'workers':workers,'implementation':record(__file__)}
    write_new(out/'AUDIT.json',result);return {k:v for k,v in result.items() if k not in ('rows','triples','config')}


def prepare_origins(audit,template,rank_qualification,output):
    a=read_json(audit);base=read_json(template);rank=qualification(rank_qualification);cfg=a['config']
    if a['protocol_id']!=PROTOCOL or a['counts']['new_origin_MACE']!=38 or a['counts']['new_origin_GFN2']!=76:
        raise InvalidArtifact('audited primary finite counts required')
    if base['model']!=cfg['model'] or base['software']!=cfg['software']:raise InvalidArtifact('actual model/runtime differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=snapshot(Path(__file__).parent,out/'implementation')
    inputs={'protocol_id':PROTOCOL,'audit':record(audit),'template':record(template),'rank_qualification':rank,
        'agreement':a['agreement'],'reference':a['reference'],'config':cfg,'implementation':impl,'new_molecular_calls':0}
    ip=out/'INPUTS.json';write_new(ip,inputs)
    mm={'protocol_id':'mace_omol_0_100m_vacuum_descriptor_v1','stage':PROTOCOL,'preparation':record(ip),
        'agreement':a['agreement'],'software':cfg['software'],'model':cfg['model'],'implementation':impl,'tasks':[],'reused':{},'verification_references':{}}
    tasks=[]
    for row in a['rows']:
        if row['pool_reuse']:continue
        for z in row['origin_reuse']['missing_native_endpoints']:
            ep=row['source']['representations']['context']['endpoints'][z];cid=row['pair_id']
            t={'task_id':cid+'__context__'+z,'case_id':cid,'source_case_id':row['case_id'],'selection_id':row['selection_id'],
                'representation':'context','metal':z,'metal_index':0,'kind':'core','xyz':ep['xyz'],'charge':ep['charge'],
                'spin_multiplicity':ep['multiplicity'],'energy_component':'MACE_OMOL_total_vacuum_energy'}
            t['cache_key']=cache_key({'task':t,'model':cfg['model'],'software':cfg['software'],'implementation':impl});mm['tasks'].append(t)
        for missing in row['origin_reuse']['missing_solvent_cells']:
            z,medium=missing.split('_',1)
            ep=row['source']['representations']['context']['endpoints'][z];cid=row['pair_id'];tid=cid+'__context__'+z+'__'+medium
            td=out/'solvent'/'tasks'/tid;td.mkdir(parents=True);xp=td/'core.xyz';shutil.copyfile(verify(ep['xyz']),xp)
            body=solvent.input_text(ep['charge'],ep['multiplicity'],medium,'native').replace('%scf\n','%scf\n MaxIter 500\n')
            inp=td/'endpoint.inp';inp.write_text(body)
            t={'task_id':tid,'case_id':cid,'case':cid,'source_case_id':row['case_id'],'selection_id':row['selection_id'],
                'metal':z,'medium':medium,'solver':'native','representation':'context','charge':ep['charge'],
                'multiplicity':ep['multiplicity'],'xyz':record(xp),'input':record(inp),'output_path':str(td/'endpoint.out')}
            t['scientific_key']=cache_key({'xyz_sha256':t['xyz']['sha256'],'charge':ep['charge'],'multiplicity':ep['multiplicity'],
                'medium':medium,'solver':'native','method':solvent.METHOD,'input_body':body,'orca':base['orca']});tasks.append(t)
    sm={'protocol_id':solvent.PROTOCOL,'method_id':solvent.METHOD,'transfer_protocol_id':PROTOCOL,'inputs':record(ip),
        'agreement':a['agreement'],'orca':base['orca'],'implementation':impl,'rank_qualification':rank,'tasks':tasks,'all_tasks':tasks,'reused':{},
        'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},'GFN2_maxiter':500,'allocated_cpus':32,'allocated_memory_MiB':65536,
        'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']}}
    write_new(out/'mace'/'manifest.json',mm);write_new(out/'solvent'/'manifest.json',sm)
    preflight=validate_origins(out/'solvent'/'manifest.json');write_new(out/'PREFLIGHT.json',preflight)
    ready={'audit':record(audit),'origin_MACE_manifest':record(out/'mace'/'manifest.json'),'origin_GFN2_manifest':record(out/'solvent'/'manifest.json'),
        'counts':a['counts'],'new_molecular_calls':0,'submission_authorized':False,'cpu_python':base['cpu_python'],'gpu_python':base['gpu_python']}
    write_new(out/'READY.json',ready);return ready


def validate_origins(manifest):
    from affordable_workflow import dry_run
    m=read_json(manifest);inputs=read_json(verify(m['inputs']));a=read_json(verify(inputs['audit']))
    qualification(verify(m['rank_qualification']));rows={r['pair_id']:r for r in a['rows']}
    expected={(r['pair_id'],z,s) for r in a['rows'] if not r['pool_reuse']
              for z,s in (x.split('_',1) for x in r['origin_reuse']['missing_solvent_cells'])}
    if (m['transfer_protocol_id']!=PROTOCOL or len(m['tasks'])!=76 or m['all_tasks']!=m['tasks'] or
        {(t['case_id'],t['metal'],t['medium']) for t in m['tasks']}!=expected or
        m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}):
        raise InvalidArtifact('finite76 original-state scalar tasks required')
    for pin in m['implementation'].values():verify(pin)
    for t in m['tasks']:
        r=rows[t['case_id']];ep=r['source']['representations']['context']['endpoints'][t['metal']]
        if (r['case_id'],r['selection_id'])!=(t['source_case_id'],t['selection_id']) or not context.reusable_state(ep,t):
            raise InvalidArtifact('source/membership/state differs')
        body=solvent.input_text(t['charge'],t['multiplicity'],t['medium'],'native').replace('%scf\n','%scf\n MaxIter 500\n')
        if verify(t['input']).read_text()!=body:raise InvalidArtifact('native scalar recipe differs')
    return dry_run(manifest)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    for op,fields in {'audit':('preparation','origin_reuse','precision_comparison','reference','agreement','output'),
                      'prepare_origins':('audit','template','rank_qualification','output'),
                      'validate_origins':('manifest',)}.items():
        q=s.add_parser(op)
        for f in fields:q.add_argument('--'+f.replace('_','-'),required=True)
        if op=='audit':q.add_argument('--workers',type=int,default=4)
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
