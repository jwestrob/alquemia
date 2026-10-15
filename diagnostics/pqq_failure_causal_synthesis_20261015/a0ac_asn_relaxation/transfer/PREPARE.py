#!/usr/bin/env python3
"""Finite transfer metadata binding and Slurm-only graph/flip/rebase preparation."""
import argparse, os, shutil, time, traceback
from COMMON import *
PROTOCOL='uniform_terminal_amide_transfer_preparation_v1'

def freeze(canonical_collection,output):
    gate=canonical_gate(canonical_collection)
    ip,ledger,invp,inv,sources=inventories();cache_registry(sources,gate)
    templatep=pin(WORK/'asn_canonical_relaxation_v1/manifest.json');template=read(templatep)
    gm=read(template['geometry_manifest']);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);(out/'implementation').mkdir()
    files={}
    for name in ('COMMON.py','PREPARE.py','GEOMETRY.py','PREPARE.sbatch','PLAN.md'):
        dest=out/'implementation'/name;shutil.copyfile(HERE/name,dest);files[name]=pin(dest)
    kernels={k:ledger['artifacts'][ledger['helpers'][k]['artifact']] for k in ('graph_inventory','fixed_flip','rebase','search','gpu_bridge')}
    kernels['conditional_solvent']=read(gate['config'])['implementation']
    c={k:template[k] for k in ('prior_config','source_preparation_helper','source_parser','generic_contact_helper','python','source_snapshot','source_modules','worker_manifest','settings','technical_recovery')}
    c.update(protocol=PROTOCOL,implementation=files,kernels=kernels,canonical_gate=gate,template=templatep,geometry_template=template['geometry_manifest'],
       ledger=ip,inventory=invp,cohorts=ledger['cohorts'],historical_records=ledger['records'],sources=sources,declared_contexts=257,declared_memberships=260,
       reuse_contexts=34,maximum_new_contexts=223,output=str(out/'execution'),new_energy_calls=0,new_optimizations=0)
    write(out/'manifest.json',c);print(json.dumps(pin(out/'manifest.json')))

def blank(s):
    r=s['record'];scope=s['scope'];cp=read(r['source_preparation']);st=read(r['state_signature'])
    return {**{k:r[k] for k in ('source_case_id','case_id','record_id','selection_id','protein','expected_class','biological_group','source_conditioning_metal')},
       'context_id':s['context_id'],'memberships':s['memberships'],'scope':scope,'atom_count':cp['new_atom_count'],'endpoint_charges':st['endpoint_charges'],'water_inventory':st['water_inventory'],
       'role':scope['role'],'original_four_modes':{t['metal']:t['active_mode_ids'] for t in scope['original_tasks']},'terminal_mode_selected_by_metal':None,
       'status':'unavailable','reason':'not prepared','new_energy_calls':0,'new_optimizations':0}

def cached(s):
    r=blank(s);reuse=s['reuse'];prior=reuse['prepared_row'];alias=reuse['executed_result_alias'];gate=read(prior['admission'])
    require(prior['status']=='prepared' and gate['pass'] and alias['status']=='complete','Existing prepared/executed alias unavailable')
    require(alias['prepared_mapping']==prior['mapping'] and alias['admission']==prior.get('prepared_state_alias',prior['admission']),'Executed prepared-state alias differs')
    for p in (prior['mapping'],prior['contact']['origin'],prior['contact']['identities'],prior['contact']['source'],reuse['producing_config'],reuse['producing_collection'],reuse['producing_execution'],reuse['producing_preparation']):check(p)
    require(read(reuse['producing_execution'])['status']=='complete','Actual cached producer failed')
    templates=prior.get('task_templates',prior.get('tasks'))
    require(templates is not None and {t['metal'] for t in templates}==set(METALS),'Actual reused task pair missing')
    modes={z:s['scope']['original_tasks'][next(i for i,t in enumerate(s['scope']['original_tasks']) if t['metal']==z)]['active_mode_ids'] for z in METALS}
    for t in templates:
        require(t['mapping']==prior['mapping'] and t['contact']==prior['contact'] and t['active_mode_ids']==modes[t['metal']],'Actual reuse mapping/mode differs')
    r.update(status='prepared',reason=None,mapping=prior['mapping'],contact=prior['contact'],admission=alias['admission'],
       task_templates=templates,original_row=prior['original_row'],reuse=reuse,executed_result_alias=alias,
       terminal_mode_selected_by_metal=prior.get('terminal_mode_selected_by_metal',{z:any(m.endswith('/chi2') and '/'+str(s['scope']['role']['resnum'])+'/' in m for m in modes[z]) for z in METALS}),
       reuse_only=True,repeat_search_permitted=False)
    for k in ('producing_config','producing_collection','producing_execution','producing_preparation'):r[k]=reuse[k]
    return r

def execute(path,sha):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Actual Slurm required before scientific imports')
    mp={'path':str(Path(path).resolve()),'sha256':sha};c=read(mp);require(c['protocol']==PROTOCOL and c['implementation']['PREPARE.py']==pin(__file__) and c['python']==pin(sys.executable),'Frozen runner/interpreter differs')
    require(canonical_gate(c['canonical_gate']['collection'])==c['canonical_gate'],'Canonical release evidence changed')
    for p in list(c['implementation'].values())+list(c['kernels'].values()):check(p)
    require(len(c['sources'])==257 and sum(len(s['memberships']) for s in c['sources'])==260,'Finite membership changed')
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time();rows=[blank(s) for s in c['sources']]
    receipt={'status':'failed','manifest':mp,'job_id':os.environ['SLURM_JOB_ID'],'new_energy_calls':0,'new_optimizations':0}
    try:
        base=load('_transfer_source_parser',c['source_parser']);receipt['placement']=base.allocation()
        old=load('_transfer_original_flip',c['source_preparation_helper']);generic=load('_transfer_generic_flip',c['generic_contact_helper'])
        graph=load('_transfer_named_graph',c['kernels']['graph_inventory']);rebase=load('_transfer_rebase',c['kernels']['rebase']);flip=load('_transfer_identity_flip',c['implementation']['GEOMETRY.py'])
        for k in ('same_res','identity','base_row','contact_check'):setattr(flip,k,getattr(generic,k))
        sys.path.insert(0,c['source_snapshot'])
        import numpy as np
        from mace_site_kinematics import Kinematics
        import adaptive_angular_proposals as angular
        for module in list(sys.modules.values()):
            f=getattr(module,'__file__',None)
            if f and Path(f).resolve().parent==Path(c['source_snapshot']).resolve():require(Path(f).name in c['source_modules'] and pin(f)==c['source_modules'][Path(f).name],'Imported source differs')
        for k in ('maximum_physical_heavy_displacement_A','angular_bound_radian','final_physical_tolerance_A','angular_tolerance_radian','physical_boundary_margin_A','angular_boundary_margin_radian'):require(angular.SETTINGS[k]==c['settings'][k],'Original physical bound differs')
        # Actual one-ion source exercises the only changed geometry routing.
        rs=next(s for s in c['sources'] if 'reuse' in s and not s['scope']['spectator_admission'])
        sc=rs['scope'];g=graph.inspect(sc);require(g['status']=='identity_and_graph_supported','Real regression graph unavailable')
        inputs={k:sc[k] for k in ('source_preparation','state_signature','protonated_source','protonation','maps')};inputs['xyz']=sc['original_xyz']
        rd=out/'one_ion_regression';rd.mkdir();rr=flip.prepare({'record':sc,'graph':g,'inputs':inputs},rd,c,np,base,Kinematics,old)
        regression=generic.regression(rr,rs['reuse']['prepared_row']['original_row']);write(rd/'REGRESSION.json',regression);receipt['real_one_ion_regression']=pin(rd/'REGRESSION.json')
        for i,s in enumerate(c['sources']):
            d=out/s['context_id'];d.mkdir()
            try:
                if 'reuse' in s:rows[i]=cached(s)
                else:
                    g=graph.inspect(s['scope']);write(d/'GRAPH.json',g);require(g['status']=='identity_and_graph_supported' and g['role']['resname']=='ASN','Original neutral anchor Asn graph unsupported')
                    sc=s['scope'];inputs={k:sc[k] for k in ('source_preparation','state_signature','protonated_source','protonation','maps')};inputs['xyz']=sc['original_xyz']
                    fd=d/'flip';fd.mkdir();fr=flip.prepare({'record':sc,'graph':g,'inputs':inputs},fd,c,np,base,Kinematics,old)
                    rd=d/'rebase';rd.mkdir();rs={'row':fr,'scope':sc,'contact_identities':pin(fd/'CONTACT_IDENTITIES.json')}
                    row=rebase.prepare_row(rs,rd,c,np,Kinematics,angular,old,base,generic)
                    require(all(v['physical_feasible'] for v in read(row['admission'])['q0_geometry_checks'].values()),'Rebased q0 physical admission failed')
                    rows[i]={**rows[i],**row,'scope':sc,'context_id':s['context_id'],'memberships':s['memberships'],'graph':pin(d/'GRAPH.json'),'reuse_only':False}
                write(d/'PREPARATION.json',rows[i])
            except Exception as e:
                rows[i].update(status='unavailable',reason=f'{type(e).__name__}: {e}',traceback=traceback.format_exc())
                if not (d/'PREPARATION.json').exists():write(d/'PREPARATION.json',rows[i])
        receipt['status']='complete'  # terminal inventory, not universal geometry admission
    except BaseException as e:receipt.update(reason=str(e),traceback=traceback.format_exc())
    finally:
        receipt.update(elapsed_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu)
        write(out/'EXECUTION.json',receipt)
        write(out/'PREPARATION.json',{'protocol':PROTOCOL,'manifest':mp,'status':receipt['status'],'rows':rows,'execution':pin(out/'EXECUTION.json'),'declared_contexts':257,'declared_memberships':260,
           'prepared_contexts':sum(r['status']=='prepared' for r in rows),'reused_contexts':sum(r.get('reuse_only',False) for r in rows),'cohorts':c['cohorts'],'ledger':c['ledger'],
           'new_energy_calls':0,'new_optimizations':0,'all_memberships_retained':True})
    return 0 if receipt['status']=='complete' else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);f=s.add_parser('freeze');f.add_argument('--canonical-collection',required=True);f.add_argument('--output',required=True)
    e=s.add_parser('execute');e.add_argument('--manifest',required=True);e.add_argument('--manifest-sha',required=True);a=p.parse_args()
    if a.action=='freeze':freeze(a.canonical_collection,a.output)
    else:sys.exit(execute(a.manifest,a.manifest_sha))
