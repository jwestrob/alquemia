#!/usr/bin/env python3
"""Six-row preparation recovery; all other terminal rows retain their exact receipts."""
import argparse,os,shutil,time,traceback
from COMMON import *
LOCAL=Path(__file__).resolve().parent
PROTOCOL='uniform_terminal_amide_transfer_preparation_recovery_v2'
PARENT_SHA='5498bcfa3d11efb5e84a063231fabf60e28499314dfcbd81d5861f2fe8ba472d'
CE='context_716ff8132b6b6536a6ce0dc9'
CAP={'context_a2a0096c29ff3a1c6f246e82','context_922a5cee4869dd437a14d56f','context_8e39d5d61eeabb1222d871f1','context_fabc02c52734b6bff84b9728','context_f18dd62cf49978f15c75f7f1'}
FAILED={CE}|CAP

def parent_gate(p,mp,execution):
    require(p['status']==execution['status']=='complete' and p['manifest']==mp==execution['manifest'] and execution['job_id']=='5205','Actual terminal5205 preparation required')
    require(len(p['rows'])==257 and sum(len(r['memberships']) for r in p['rows'])==260 and len({r['context_id'] for r in p['rows']})==257,'Original population differs')
    require(p['prepared_contexts']==251 and p['reused_contexts']==34 and {r['context_id'] for r in p['rows'] if r['status']!='prepared'}==FAILED,'Exact six original failures required')
    for r in p['rows']:
        if r['context_id'] in FAILED:require(r['reason']==('ValueError: Unsupported source element/alternate location' if r['context_id']==CE else 'ValueError: Unsupported selected moving cap boundary'),'Original failure class differs')

def freeze(previous,output):
    pp=pin(previous);require(pp['sha256']==PARENT_SHA,'Actual completed5205 collection differs');p=read(pp);mp=p['manifest'];old=read(mp);ep=p['execution'];parent_gate(p,mp,read(ep))
    rowpins={r['context_id']:pin(Path(pp['path']).parent/r['context_id']/'PREPARATION.json') for r in p['rows']}
    require(all(read(rowpins[r['context_id']])==r for r in p['rows']),'Original per-row receipts differ')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir();files={}
    for name in ('COMMON.py','RECOVER.py','SOURCE.py','GEOMETRY.py','REBASE.py','CAP_CHECK.py','RUN.sbatch','PLAN.md'):
        dest=out/'implementation'/name;shutil.copyfile(LOCAL/name,dest);files[name]=pin(dest)
    c=copy.deepcopy(old);c.update(protocol=PROTOCOL,implementation=files,output=str(out/'execution'),
       previous_preparation=pp,previous_manifest=mp,previous_execution=ep,previous_row_receipts=rowpins,
       retry_context_ids=[r['context_id'] for r in p['rows'] if r['context_id'] in FAILED],exact_prior_row_reuse=251,
       new_energy_calls=0,new_optimizations=0)
    write(out/'manifest.json',c);print(json.dumps(pin(out/'manifest.json')))

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports')
    mp={'path':str(Path(path).resolve()),'sha256':sha};c=read(mp)
    require(c['protocol']==PROTOCOL and c['implementation']['RECOVER.py']==pin(__file__) and c['python']==pin(sys.executable),'Frozen runner/interpreter differs')
    for p in list(c['implementation'].values())+list(c['kernels'].values()):check(p)
    require(canonical_gate(c['canonical_gate']['collection'])==c['canonical_gate'],'Canonical gate changed')
    parent=read(c['previous_preparation']);parent_gate(parent,c['previous_manifest'],read(c['previous_execution']));oldc=read(c['previous_manifest'])
    require(c['previous_preparation']['sha256']==PARENT_SHA and c['sources']==oldc['sources'] and c['settings']==oldc['settings'] and c['kernels']==oldc['kernels'],'Population or original kernels changed')
    require(set(c['retry_context_ids'])==FAILED and c['exact_prior_row_reuse']==251,'Finite recovery scope changed')
    rows=copy.deepcopy(parent['rows']);sources={s['context_id']:s for s in c['sources']};original={r['context_id']:r for r in parent['rows']}
    for r in rows:require(read(c['previous_row_receipts'][r['context_id']])==r,'Saved source row changed')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    receipt={'status':'failed','manifest':mp,'job_id':os.environ['SLURM_JOB_ID'],'previous_preparation':c['previous_preparation'],'new_energy_calls':0,'new_optimizations':0}
    try:
        base=load('_recovery_source_parser',c['source_parser']);receipt['placement']=base.allocation()
        old=load('_recovery_original_flip',c['source_preparation_helper']);generic=load('_recovery_generic_flip',c['generic_contact_helper'])
        graph=load('_recovery_graph',c['kernels']['graph_inventory']);rebase=load('_recovery_rebase',c['implementation']['REBASE.py']);cap=load('_recovery_caps',c['implementation']['CAP_CHECK.py'])
        flip=load('_recovery_flip',c['implementation']['GEOMETRY.py']);source=load('_recovery_declared_source',c['implementation']['SOURCE.py'])
        for k in ('same_res','identity','base_row','contact_check'):setattr(flip,k,getattr(generic,k))
        sys.path.insert(0,c['source_snapshot'])
        import numpy as np
        from mace_site_kinematics import Kinematics
        import adaptive_angular_proposals as angular
        for module in list(sys.modules.values()):
            f=getattr(module,'__file__',None)
            if f and Path(f).resolve().parent==Path(c['source_snapshot']).resolve():require(Path(f).name in c['source_modules'] and pin(f)==c['source_modules'][Path(f).name],'Imported source differs')
        for k in ('maximum_physical_heavy_displacement_A','angular_bound_radian','final_physical_tolerance_A','angular_tolerance_radian','physical_boundary_margin_A','angular_boundary_margin_radian'):require(angular.SETTINGS[k]==c['settings'][k],'Original physical bound differs')
        rs=next(s for s in c['sources'] if 'reuse' in s and not s['scope']['spectator_admission']);sc=rs['scope'];facade=source.DeclaredSource(old,sc,base)
        require(facade.physical_source(check(sc['protonated_source']),base)==old.physical_source(check(sc['protonated_source']),base),'Real supported source parser regression differs')
        g=graph.inspect(sc);require(g['status']=='identity_and_graph_supported','Real regression graph unavailable');inp={k:sc[k] for k in ('source_preparation','state_signature','protonated_source','protonation','maps')};inp['xyz']=sc['original_xyz']
        rd=out/'one_ion_regression';rd.mkdir();rr=flip.prepare({'record':sc,'graph':g,'inputs':inp},rd,c,np,base,Kinematics,facade)
        reg=generic.regression(rr,rs['reuse']['prepared_row']['original_row']);write(rd/'REGRESSION.json',reg);receipt['real_one_ion_regression']=pin(rd/'REGRESSION.json')
        for i,previous in enumerate(parent['rows']):
            cid=previous['context_id']
            if cid not in FAILED:continue
            d=out/cid;d.mkdir();s=sources[cid];sc=s['scope']
            try:
                g=graph.inspect(sc);write(d/'GRAPH.json',g);require(g['status']=='identity_and_graph_supported' and g['role']['resname']=='ASN','Original neutral anchor Asn graph unsupported')
                facade=source.DeclaredSource(old,sc,base);write(d/'SOURCE_IDENTITY.json',facade.evidence)
                require((facade.active['element']=='Ce')==(cid==CE),'Declared source-extension scope changed')
                inp={k:sc[k] for k in ('source_preparation','state_signature','protonated_source','protonation','maps')};inp['xyz']=sc['original_xyz']
                fd=d/'flip';fd.mkdir();fr=flip.prepare({'record':sc,'graph':g,'inputs':inp},fd,c,np,base,Kinematics,facade)
                rd=d/'rebase';rd.mkdir();row=rebase.prepare_row({'row':fr,'scope':sc,'contact_identities':pin(fd/'CONTACT_IDENTITIES.json')},rd,c,np,Kinematics,angular,facade,base,generic)
                require(all(v['physical_feasible'] for v in read(row['admission'])['q0_geometry_checks'].values()),'Rebased q0 physical admission failed')
                row['context_id']=cid
                checks=cap.check_row(row,np,Kinematics,angular,generic);write(d/'CAP_AXIS_CHECK.json',checks)
                require(checks['status']==('passed' if cid in CAP else 'not_applicable'),'Actual cap-axis qualification failed')
                rows[i]={**previous,**row,'scope':sc,'context_id':cid,'memberships':s['memberships'],'graph':pin(d/'GRAPH.json'),'reuse_only':False,
                   'source_identity_admission':pin(d/'SOURCE_IDENTITY.json'),'cap_axis_check':pin(d/'CAP_AXIS_CHECK.json'),'recovery_from':c['previous_row_receipts'][cid]}
                rows[i].pop('traceback',None)
            except Exception as e:rows[i]={**previous,'status':'unavailable','reason':f'{type(e).__name__}: {e}','traceback':traceback.format_exc(),'recovery_from':c['previous_row_receipts'][cid]}
            write(d/'PREPARATION.json',rows[i])
        require(all(r==original[r['context_id']] for r in rows if r['context_id'] not in FAILED),'A non-retried prior row changed')
        receipt.update(status='complete',exact_prior_row_reuse=251,retried_context_ids=c['retry_context_ids'])
    except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc())
    finally:
        receipt.update(elapsed_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu);write(out/'EXECUTION.json',receipt)
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'manifest':mp,'status':receipt['status'],'rows':rows,'execution':pin(out/'EXECUTION.json'),
          'declared_contexts':257,'declared_memberships':260,'prepared_contexts':sum(r['status']=='prepared' for r in rows),'reused_contexts':sum(r.get('reuse_only',False) for r in rows),
          'cohorts':c['cohorts'],'ledger':c['ledger'],'previous_preparation':c['previous_preparation'],'prior_row_receipts':c['previous_row_receipts'],
          'new_energy_calls':0,'new_optimizations':0,'all_memberships_retained':True})
    return 0 if receipt['status']=='complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);f=s.add_parser('freeze');f.add_argument('--previous-preparation',required=True);f.add_argument('--output',required=True)
    e=s.add_parser('execute');e.add_argument('--manifest',required=True);e.add_argument('--manifest-sha',required=True);a=p.parse_args()
    if a.action=='freeze':freeze(a.previous_preparation,a.output)
    else:sys.exit(execute(a.manifest,a.manifest_sha))
