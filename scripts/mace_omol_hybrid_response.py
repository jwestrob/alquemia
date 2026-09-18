"""Bounded metal response of the matched-H vacuum DFT/masked-MACE hybrid."""
from __future__ import annotations
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,BOHR_TO_A,read_json,record,verify,write_new,xyz,cache_key,paired
from affordable_response import read_engrad
from mace_hybrid import EV_TO_KCAL,accepted_attempt,check_atoms,write_xyz
from mace_file_checks import cached_file_checks
from mace_omol_gradient_run import CONFIG,check_parent,model as gradient_model
from mace_omol_ablation import ADAPTER as MASK,COMPONENT
from mace_omol_ablation_run import SEMANTICS
from mace_metal_response import grid,shift,matrix
from mace_bounded_response import ball
POLICY='matched_normalized_H_vacuum_hybrid_bounded_metal_response_v1'
STAGE='matched_hybrid_response'
CASES=('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9')


@cached_file_checks
def sources(preparation,static_report,ggr_gradients,full_check=False):
    p=read_json(preparation);r=read_json(static_report);gr=read_json(ggr_gradients)
    if set(p['cases'])!=set(CASES) or r['status']!='complete' or not r['numerical_gate_pass'] or gr['status']!='complete' or not gr['numerical_gate_pass']:
        raise InvalidArtifact('complete matched static hybrid and gradient qualification required')
    qm=read_json(verify(r['quantum']['manifest']));cm=read_json(verify(r['MACE']['manifest']));gm=read_json(verify(gr['manifest']))
    if qm['preparation']!=record(preparation) or cm['preparation']!=record(preparation):raise InvalidArtifact('static preparation differs')
    if full_check:
        from mace_omol_matched_h import collect_quantum,collect_mace
        if collect_quantum(verify(r['quantum']['manifest']))['rows']!=r['quantum']['rows'] or collect_mace(verify(r['MACE']['manifest']))['rows']!=r['MACE']['rows']:
            raise InvalidArtifact('static archive differs from actual completed receipts')
    states={};fulls={};reused={}
    for name,ref in p['cases'].items():
        c=read_json(verify(ref));gid=c['global_id'];states[name]={}
        for metal,e in c['endpoints'].items():
            key=name+'_'+metal;d=r['quantum']['rows'][key];raw=read_engrad(verify(d['engrad']))
            for v in read_json(verify(d['receipt']))['artifacts'].values():verify(v)
            if abs(raw['energy_Ha']-d['energy_hartree'])>1e-8 or not np.allclose(raw['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A,d['gradient_kcal_mol_per_A'],atol=1e-10,rtol=1e-10):
                raise InvalidArtifact('actual vacuum gradient differs')
            qr=next(t for t in qm['tasks'] if t['task_id']==key)
            if qr['xyz']!=e['xyz'] and qr['xyz']['sha256']!=e['xyz']['sha256']:raise InvalidArtifact('vacuum coordinates differ')
            core={'xyz':e['xyz'],'charge':e['charge'],'spin_multiplicity':1,'state':e['state'],'metal_index':0,
                  'archived_center_energy_eV':r['MACE']['rows'][key]['energy_eV'],'case_mapping':ref}
            w=p['whole_reuses'][gid][metal];parent=read_json(verify(read_json(verify(w['collection']))['manifest']))
            wt=next(t for t in parent['tasks'] if t['task_id']==w['selected']['task_id'])
            whole={k:copy.deepcopy(wt[k]) for k in ('xyz','charge','spin_multiplicity','state','metal_index','preparation','assembly','microstate','explicit_waters','evidence') if k in wt}
            whole['archived_center_energy_eV']=w['selected']['result']['energy_eV']
            if whole['xyz']!=w['selected']['input'] or whole['xyz']!=read_json(verify(c['normalized_global_preparation']))['endpoints'][metal]['xyz']:
                raise InvalidArtifact('normalized full physical coordinates differ')
            if gid in fulls and metal in fulls[gid] and fulls[gid][metal]!=whole:raise InvalidArtifact('GGR physical system differs between partitions')
            fulls.setdefault(gid,{})[metal]=whole
            states[name][metal]={'core':core,'full':whole,'global_id':gid,'DFT_source':d,'DFT_gradient_kcal_mol_A':d['gradient_kcal_mol_per_A'][0]}
            if gid=='GGR_1GLG':
                gt=next(t for t in gm['tasks'] if t['metal']==metal and t['gradient_variant']=='center' and t['derivative_backend']=='checkpointed')
                result=gr['rows'][gt['task_id']];mp=verify(gr['manifest'])
                if gt['source_xyz']!=whole['xyz'] or gt['charge']!=whole['charge'] or not np.allclose([a[1:] for a in xyz(verify(gt['xyz']))],[a[1:] for a in xyz(verify(whole['xyz']))],atol=1e-12,rtol=0):
                    raise InvalidArtifact('GGR analytic center is not the exact normalized full state')
                if not any(accepted_attempt(d,gt,mp)==result for d in (mp.parent/'execution'/gt['task_id']).glob('attempt_*')):
                    raise InvalidArtifact('GGR gradient lacks a valid execution receipt')
                reused[gid+'_'+metal]={'result':result,'manifest':gr['manifest'],'task_id':gt['task_id'],'metal_index':whole['metal_index']}
    return p,r,gr,gm,states,fulls,reused


def task(source,name,metal,kind,point,u,output):
    rows,index=shift(xyz(verify(source['xyz'])),metal,u)
    if index!=source['metal_index']:raise InvalidArtifact('selected source metal index differs')
    tid=kind+'__'+name+'_'+metal+'__'+point;xp=Path(output)/(tid+'.xyz');write_xyz(xp,rows)
    return {**copy.deepcopy(source),'task_id':tid,'case_id':name,'metal':metal,'kind':kind,'variant':'primary',
            'xyz':record(xp),'source_xyz':source['xyz'],'point':point,'displacement_A':list(u),
            'state':check_atoms(rows,source['charge']),'energy_only':point!='center',
            'descriptor_gradient_experiment':POLICY,'derivative_backend':'checkpointed',
            'energy_component':COMPONENT,'charge_feature_adapter':MASK,'output_semantics':SEMANTICS,
            'capture_native_readout':True}


@cached_file_checks
def prepare(preparation,static_report,ggr_gradients,agreement,output):
    from mace_omol import common,seal
    start=time.monotonic();cpu=time.process_time();p,r,gr,gm,states,fulls,reused=sources(preparation,static_report,ggr_gradients,True)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);jobs=[]
    groups={'centers':[(v['core'],name,metal,'core','center',[0.,0.,0.]) for name,row in states.items() for metal,v in row.items()]}
    groups['centers'] += [(s,gid,metal,'full','center',[0.,0.,0.]) for gid,row in fulls.items() if gid!='GGR_1GLG' for metal,s in row.items()]
    for gid,row in fulls.items():groups['grid_'+gid]=[(s,gid,metal,'full',point,u) for metal,s in row.items() for point,u in grid().items() if point!='center']
    for group,items in groups.items():
        _,d,m=common(verify(p['inventory']),verify(p['software']),agreement,out/group,STAGE)
        m.update(protocol_id=POLICY,model=gradient_model(verify(p['software'])),gradient_settings=CONFIG,
                 core_report=gm['core_report'],source_preparation=record(preparation),static_report=record(static_report),
                 ggr_gradients=record(ggr_gradients),group=group,output_semantics=SEMANTICS,
                 tasks=[task(s,name,metal,kind,point,u,d) for s,name,metal,kind,point,u in items])
        check_parent(m);seal(d,m);jobs.append(record(d/'manifest.json'))
    result={'protocol_id':POLICY,'agreement':record(agreement),'source_preparation':record(preparation),
            'static_report':record(static_report),'ggr_gradients':record(ggr_gradients),'jobs':jobs,
            'reused_GGR_whole_gradients':reused,'new_gradient_centers':12,'new_whole_grid_calls':216,
            'preparation_wall_seconds':time.monotonic()-start,'preparation_CPU_seconds':time.process_time()-cpu}
    write_new(out/'preparation.json',result);return {'status':'prepared','jobs':jobs,'new_MACE_calls':228}


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);p,r,gr,gm,states,fulls,reused=sources(verify(m['source_preparation']),verify(m['static_report']),verify(m['ggr_gradients']))
    if (m['protocol_id']!=POLICY or m['stage']!=STAGE or m['model']!=gradient_model(verify(p['software'])) or
            m['software']!=p['software'] or m['gradient_settings']!=CONFIG or m['core_report']!=gm['core_report']):
        raise InvalidArtifact('matched hybrid response method/settings differ')
    check_parent(m)
    for ref in [m['agreement'],*m['implementation'].values()]:verify(ref)
    if m['group']=='centers':
        expected={'core__'+name+'_'+metal+'__center':(v['core'],name,metal,'core','center',[0.,0.,0.]) for name,row in states.items() for metal,v in row.items()}
        expected.update({'full__'+gid+'_'+metal+'__center':(s,gid,metal,'full','center',[0.,0.,0.]) for gid,row in fulls.items() if gid!='GGR_1GLG' for metal,s in row.items()})
    elif m['group'].startswith('grid_') and m['group'][5:] in fulls:
        gid=m['group'][5:];expected={'full__'+gid+'_'+metal+'__'+point:(s,gid,metal,'full',point,u) for metal,s in fulls[gid].items() for point,u in grid().items() if point!='center'}
    else:raise InvalidArtifact('undeclared task group')
    if len(m['tasks'])!=len(expected) or {t['task_id'] for t in m['tasks']}!=set(expected):raise InvalidArtifact('finite response task inventory differs')
    for t in m['tasks']:
        s,name,metal,kind,point,u=expected[t['task_id']]
        fields={**s,'case_id':name,'metal':metal,'kind':kind,'variant':'primary','source_xyz':s['xyz'],'point':point,
                'displacement_A':u,'energy_only':point!='center','descriptor_gradient_experiment':POLICY,'derivative_backend':'checkpointed',
                'energy_component':COMPONENT,'charge_feature_adapter':MASK,'output_semantics':SEMANTICS,'capture_native_readout':True}
        for k,v in fields.items():
            if k!='xyz' and t.get(k)!=v:raise InvalidArtifact('task physical or method field differs: '+k)
        wanted,_=shift(xyz(verify(s['xyz'])),metal,u);actual=xyz(verify(t['xyz']))
        if [x[0] for x in actual]!=[x[0] for x in wanted] or not np.allclose([x[1:] for x in actual],[x[1:] for x in wanted],atol=5e-10,rtol=0):raise InvalidArtifact('source displacement differs')
        if check_atoms(actual,t['charge'])!=t['state']:raise InvalidArtifact('charge/electron state differs')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('response scientific cache differs')
    return {'status':'pass','manifest':record(manifest),'tasks':len(m['tasks'])}


@cached_file_checks
def collect(manifest):
    mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[];checks=[]
    for t in m['tasks']:
        candidates=[]
        for d in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(d,t,mp);attempts.append({'task_id':t['task_id'],'attempt':str(d),'accepted':r is not None})
            if r is not None:candidates.append(r)
        r=candidates[-1] if candidates else {'status':'unavailable','energy_eV':None};rows[t['task_id']]=r
        if r['status']=='computed' and t['point']=='center':
            error=(r['energy_eV']-t['archived_center_energy_eV'])*EV_TO_KCAL
            checks.append({'task_id':t['task_id'],'name':'center_scalar_agreement','error_kcal':error,'tolerance':.01,'pass':abs(error)<=.01})
    return {'status':'complete' if all(r['status']=='computed' for r in rows.values()) else 'incomplete',
            'manifest':record(mp),'rows':rows,'attempts':attempts,'checks':checks,
            'numerical_checks_pass':all(c['pass'] for c in checks),'physical_force_claim':False}


@cached_file_checks
def assess(prepared,collections,output):
    prep=read_json(prepared);p,r,gr,gm,states,fulls,reused=sources(verify(prep['source_preparation']),verify(prep['static_report']),verify(prep['ggr_gradients']),True)
    rows={};pins=[];checks=[]
    for path in read_json(collections)['collections']:
        c=read_json(path);mp=verify(c['manifest']);validate(mp);actual=collect(mp)
        if c!=actual or c['status']!='complete' or not c['numerical_checks_pass']:raise InvalidArtifact('actual initial calculations incomplete/failed')
        if c['manifest'] not in prep['jobs'] or set(rows)&set(c['rows']):raise InvalidArtifact('wrong or duplicated initial calculation')
        rows.update(c['rows']);pins.append(record(path));checks.extend(c['checks'])
    if len(pins)!=4 or len(rows)!=228:raise InvalidArtifact('incomplete initial inventory')
    for gidmetal,item in reused.items():rows['full__'+gidmetal+'__center']=item['result']
    matrices={};odd_checks={}
    for gid,ends in fulls.items():
        matrices[gid]={};odd_checks[gid]={}
        for metal,s in ends.items():
            center=rows['full__'+gid+'_'+metal+'__center'];values={point:(rows['full__'+gid+'_'+metal+'__'+point]['energy_eV']-center['energy_eV'])*EV_TO_KCAL for point in grid()}
            matrices[gid][metal]={scale:matrix(values,scale) for scale in ('coarse','fine')}
            gradient=np.load(verify(center['gradient']))[s['metal_index']]*EV_TO_KCAL;oc=[]
            for scale,h in [('coarse',.02),('fine',.01)]:
                for i,axis in enumerate('xyz'):
                    odd=(values[scale+'_'+axis+'p']-values[scale+'_'+axis+'m'])/2;pred=h*gradient[i];tol=max(.01,.01*abs(pred))
                    oc.append({'scale':scale,'axis':axis,'odd_kcal':odd,'gradient_prediction_kcal':float(pred),'tolerance_kcal':float(tol),'error_kcal':float(odd-pred),'pass':bool(abs(odd-pred)<=tol)})
            odd_checks[gid][metal]=oc
    result_rows={}
    for name,ends in states.items():
        result_rows[name]={}
        for metal,s in ends.items():
            gid=s['global_id'];gcore=np.load(verify(rows['core__'+name+'_'+metal+'__center']['gradient']))[0]*EV_TO_KCAL
            gfull=np.load(verify(rows['full__'+gid+'_'+metal+'__center']['gradient']))[s['full']['metal_index']]*EV_TO_KCAL
            g=np.array(s['DFT_gradient_kcal_mol_A'])+gfull-gcore;k=matrices[gid][metal]
            fine=ball(k['fine'],g);coarse=ball(k['coarse'],g);pred=copy.deepcopy(fine)
            if fine['status']=='solved' and coarse['status']=='solved':
                u=np.array(fine['displacement_A']);error=float(.5*u@(k['fine']-k['coarse'])@u);delta=float(np.linalg.norm(u-coarse['displacement_A']))
                pred.update(refinement_error_kcal_mol=error,refinement_displacement_A=delta,
                    status='eligible_for_native_check' if abs(error)<=.02 and delta<=.01 and all(c['pass'] for c in odd_checks[gid][metal]) else 'numerical_response_not_qualified')
            result_rows[name][metal]={'prediction':pred,'gradient_kcal_mol_A':g.tolist(),'DFT_gradient_kcal_mol_A':s['DFT_gradient_kcal_mol_A'],
                                     'context_gradient_kcal_mol_A':(gfull-gcore).tolist(),'matrices_kcal_mol_A2':{key:value.tolist() for key,value in k.items()}}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);source=out/'assessment_source.py';source.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete','protocol_id':POLICY,'prepared':record(prepared),'collections':pins,'implementation':record(source),
            'rows':result_rows,'odd_checks':odd_checks,'center_checks':checks,'predictions_are_not_validated_scores':True,
            'baseline_changed':False,'qualified_scores':None,'calibrated_class':None}
    write_new(out/'result.json',result);return {'status':'complete','predictions':{name:{m:v['prediction'] for m,v in ends.items()} for name,ends in result_rows.items()}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for name,keys in [('prepare',('preparation','static-report','ggr-gradients','agreement','output')),('validate',('manifest',)),('collect',('manifest',)),('assess',('prepared','collections','output'))]:
        q=sub.add_parser(name)
        for key in keys:q.add_argument('--'+key,required=True)
    args=vars(parser.parse_args());command=args.pop('command');print(json.dumps(globals()[command](**args),indent=2))
