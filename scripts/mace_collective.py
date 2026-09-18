"""Opt-in whole-protein simultaneous substitutions in the existing masked MACE model."""
from __future__ import annotations
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import Z,EV_TO_KCAL,rotation,write_xyz,accepted_attempt
from mace_omol_locality import endpoint
from mace_omol_ablation_run import PROTOCOL as PARENT,descriptor_model

POLICY='masked_omol_whole_protein_collective_metal_substitution_v1'
STAGE='collective_metals'
COUNTS={'A0A7':6,'HEW5':8,'RTX':8,'PARV_4CPV':2,'AEQ_1SL8':3}


def joint_state(atoms,charge,indices,metal):
    ids=set(indices)
    if (metal not in ('Ca','La') or not ids or len(ids)!=len(indices) or charge!=int(charge)
        or any(not isinstance(i,int) or not 0<=i<len(atoms) for i in ids)
        or ids!={i for i,a in enumerate(atoms) if a[0] in ('Ca','La')}
        or any(atoms[i][0]!=metal for i in ids) or any(a[0] not in Z for a in atoms)):
        raise InvalidArtifact('explicit collective metal inventory differs')
    electrons=sum(Z[a[0]] for a in atoms)-charge
    if electrons%2 or not np.isfinite([a[1:] for a in atoms]).all():raise InvalidArtifact('collective singlet parity/coordinate failure')
    return {'atoms':len(atoms),'all_electron_count':electrons,'multiplicity':1}


@cached_file_checks
def sources(config):
    cfg=read_json(config);kh=read_json(verify(cfg['khoury']));mu=read_json(verify(cfg['multisite']));groups={};parent=None;reference=None
    if kh['status']!='complete' or not kh['numerical_gate_pass'] or mu['status']!='complete' or not mu['numerical_gate_pass']:raise InvalidArtifact('complete actual parent studies required')
    records={g:[(x['site'],{**x['result'],'preparation':x['preparation']}) for x in kh['rows'] if x['domain']==g] for g in ('A0A7','HEW5','RTX')}
    records.update({g:[(n,read_json(verify(mu['scores'][n]['report']))) for n in mu['site_order'] if n.startswith(g+'_')] for g in ('PARV_4CPV','AEQ_1SL8')})
    for group,items in records.items():
        if len(items)!=COUNTS[group]:raise InvalidArtifact('fixed site count changed')
        values=[];physical=None;ca_values=[];first=None
        for label,r in items:
            if r['protocol_id']!=PARENT or not r['numerical_gate_pass']:raise InvalidArtifact('parent descriptor protocol differs')
            if reference is None:reference=r['factorization_reference']
            if reference!=r['factorization_reference']:raise InvalidArtifact('atomic model reference differs')
            p=read_json(verify(r['preparation']));verify(p['source']);verify(p['source_preparation'])
            atoms=p['physical_atoms'];byid={a['id']:('Ca' if a['kind'] in ('background_metal','selected_metal') else a['element'],*a['xyz_A']) for a in atoms}
            if len(byid)!=len(atoms):raise InvalidArtifact('duplicate physical identity')
            if physical is None:physical=byid
            if physical!=byid:raise InvalidArtifact('site substitutions have different physical systems')
            for metal in ('Ca','La'):
                t,_=endpoint(r['rows'][metal+'_bound_primary']);m=read_json(verify(r['rows'][metal+'_bound_primary']['manifest']))
                if parent is None:parent=m
                if m['model']!=parent['model'] or m['software']!=parent['software']:raise InvalidArtifact('actual scalar method differs across proteins')
                if t['preparation']!=r['preparation']:raise InvalidArtifact('endpoint/preparation differs')
                expected=[(metal if a['kind']=='selected_metal' else a['element'],*a['xyz_A']) for a in atoms]
                if xyz(verify(t['xyz']))!=expected:raise InvalidArtifact('source-mapped endpoint differs')
                if metal=='Ca':ca_values.append(r['rows']['Ca_bound_primary']['energy_eV'])
            values.append({'site':label,'single_La_energy_eV':r['rows']['La_bound_primary']['energy_eV'],'old_score_model_kcal':r['R_mask_model_kcal'],'source_preparation':r['preparation']})
            if first is None:
                first={'preparation':r['preparation'],'Ca':r['rows']['Ca_bound_primary'],'task':copy.deepcopy(next(t for t in read_json(verify(r['manifest']))['tasks'] if t['task_id']=='Ca_bound_primary'))}
        if max(ca_values)-min(ca_values)>.01/EV_TO_KCAL:raise InvalidArtifact('archived all-Ca permutation gate failed')
        prep=read_json(verify(first['preparation']));ids=[i for i,a in enumerate(prep['physical_atoms']) if a['kind'] in ('selected_metal','background_metal')]
        if len(ids)!=COUNTS[group]:raise InvalidArtifact('mapped metal inventory differs')
        if first['task']['charge']!=prep['protein_charge_e']+prep['cofactor_charge_e']+2*len(ids):
            raise InvalidArtifact('archived all-Ca formal charge closure differs')
        groups[group]={**first,'metal_indices':ids,'metal_atom_ids':[prep['physical_atoms'][i]['id'] for i in ids],
                       'ordered_sites':values,'metal_count':len(ids),'evidence':prep['evidence'],'assembly':prep['assembly'],
                       'explicit_waters':prep['explicit_waters'],'protein_charge_e':prep['protein_charge_e'],'source_microstate':prep['microstate']}
    ref=read_json(verify(reference))
    if ref['status']!='pass' or ref['conversion']!=EV_TO_KCAL or not all(c['pass'] for c in ref['checks']):raise InvalidArtifact('unsupported isolated-atom reference')
    return groups,parent,reference,kh['ggr_scores']


def make_task(group,g,variant,out):
    t=copy.deepcopy(g['task']);t.pop('cache_key');atoms=xyz(verify(t['xyz']));indices=g['metal_indices'];metal='Ca' if variant=='Ca_replay' else 'La'
    coords=[(metal if i in indices else a[0],*a[1:]) for i,a in enumerate(atoms)];source=t['xyz'];selected=t['metal_index']
    permutation=list(range(len(coords)))
    if variant=='La_rotate':coords=[(a[0],*(np.array(a[1:])@rotation().T)) for a in coords]
    if variant=='La_permute':
        permutation.reverse();coords=[coords[i] for i in permutation];indices=[len(atoms)-1-i for i in indices];selected=len(atoms)-1-selected
    charge=t['charge']+(0 if metal=='Ca' else g['metal_count']);tid=group+'_'+variant;xp=out/(tid+'.xyz');write_xyz(xp,coords)
    t.update(task_id=tid,case_id=group,metal=metal,metal_index=selected,charge=charge,xyz=record(xp),source_xyz=source,
             state=joint_state(coords,charge,indices,metal),variant=variant,position='all_source_sites_occupied',joint_metal_indices=indices,
             joint_metal_atom_ids=g['metal_atom_ids'],joint_protocol=POLICY,source_atom_permutation=permutation,
             source_microstate=g['source_microstate'],microstate='fixed_source_protonation; all_recorded_sites_'+metal+('2plus' if metal=='Ca' else '3plus'))
    return t


def task_specs():return [(g,'La_primary') for g in COUNTS]+[('RTX',v) for v in ('Ca_replay','La_rotate','La_permute')]


@cached_file_checks
def prepare(config,agreement,output):
    from mace_omol import common,seal
    start=time.monotonic();cpu=time.process_time();groups,parent,ref,ggr=sources(config)
    _,out,m=common(verify(parent['inventory']),verify(parent['software']),agreement,output,STAGE)
    m.update(protocol_id=POLICY,model=descriptor_model(verify(parent['software'])),config=record(config),groups=groups,factorization_reference=ref,ggr_scores=ggr,
             tasks=[make_task(g,groups[g],v,out) for g,v in task_specs()],preparation_wall_seconds=time.monotonic()-start,preparation_CPU_seconds=time.process_time()-cpu)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    m=read_json(manifest);groups,parent,ref,ggr=sources(verify(m['config']))
    if (m['protocol_id']!=POLICY or m['stage']!=STAGE or m['groups']!=groups or m['model']!=descriptor_model(verify(parent['software']))
        or m['software']!=parent['software'] or m['factorization_reference']!=ref or m['ggr_scores']!=ggr or len(m['tasks'])!=8):raise InvalidArtifact('collective source/model inventory changed')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    for t,(group,variant) in zip(m['tasks'],task_specs()):
        g=groups[group];src=g['task'];old=xyz(verify(src['xyz']));actual=xyz(verify(t['xyz']));metal='Ca' if variant=='Ca_replay' else 'La';q=src['charge']+(0 if metal=='Ca' else g['metal_count'])
        want=[(metal if i in g['metal_indices'] else a[0],*a[1:]) for i,a in enumerate(old)]
        perm=list(range(len(want)));indices=g['metal_indices'];selected=src['metal_index']
        if variant=='La_rotate':want=[(a[0],*(np.array(a[1:])@rotation().T)) for a in want]
        if variant=='La_permute':perm.reverse();want=[want[i] for i in perm];indices=[len(old)-1-i for i in indices];selected=len(old)-1-selected
        if [a[0] for a in want]!=[a[0] for a in actual] or not np.allclose([a[1:] for a in want],[a[1:] for a in actual],atol=1e-12,rtol=0):raise InvalidArtifact('collective source geometry differs')
        fields={**{k:v for k,v in src.items() if k not in ('cache_key','xyz','task_id','case_id','metal','metal_index','charge','state','variant','position','source_xyz','microstate')},
                'task_id':group+'_'+variant,'case_id':group,'metal':metal,'metal_index':selected,'charge':q,'state':joint_state(actual,q,indices,metal),'variant':variant,'position':'all_source_sites_occupied',
                'source_xyz':src['xyz'],'joint_metal_indices':indices,'joint_metal_atom_ids':g['metal_atom_ids'],'joint_protocol':POLICY,'source_atom_permutation':perm,
                'source_microstate':g['source_microstate'],'microstate':'fixed_source_protonation; all_recorded_sites_'+metal+('2plus' if metal=='Ca' else '3plus')}
        if {k:v for k,v in t.items() if k not in ('xyz','cache_key')}!=fields:raise InvalidArtifact('collective task settings differ')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):raise InvalidArtifact('collective scientific cache differs')
    return {'status':'pass','tasks':8,'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};attempts=[]
    for t in m['tasks']:
        good=[]
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r=accepted_attempt(a,t,mp);attempts.append({'task_id':t['task_id'],'path':str(a),'accepted':r is not None,'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
            if r is not None:good.append(r)
        rows[t['task_id']]=good[-1] if good else {'status':'unavailable','energy_eV':None}
    return {'manifest':record(mp),'rows':rows,'attempts':attempts,'status':'complete' if all(v['status']=='computed' for v in rows.values()) else 'incomplete'}


@cached_file_checks
def report(manifest,output):
    start=time.monotonic();c=collect(manifest);m=read_json(manifest)
    if c['status']!='complete':raise InvalidArtifact('collective endpoints incomplete')
    rows=c['rows'];checks=[];scores={};ref=read_json(verify(m['factorization_reference']));atom=ref['Ca_minus_La_disconnected_atom_model_eV']
    def check(name,error,tol=.01):checks.append({'name':name,'error_model_kcal':float(error),'tolerance':tol,'pass':bool(abs(error)<=tol)})
    kh=read_json(verify(read_json(verify(m['config']))['khoury']))
    ggr=read_json(verify(kh['ggr_sources'][0]))['scores']['GGR_1GLG']
    for metal in ('Ca','La'):endpoint(ggr['endpoints'][metal]['bound'])
    one_ca=ggr['endpoints']['Ca']['bound']['energy_eV'];one_la=ggr['endpoints']['La']['bound']['energy_eV']
    check('actual_GGR_n1_score_identity',(one_ca-one_la-atom)*EV_TO_KCAL-ggr['R_mask_model_kcal'],1e-6)
    check('actual_GGR_n1_nonadditivity',one_la+0*one_ca-one_la,1e-6)
    for name,r in rows.items():check(name+'_readout',r['native_readout']['component_sum_error_kcal_mol'])
    for variant in ('La_rotate','La_permute'):check('RTX_'+variant,(rows['RTX_'+variant]['energy_eV']-rows['RTX_La_primary']['energy_eV'])*EV_TO_KCAL)
    check('RTX_all_Ca_replay',(rows['RTX_Ca_replay']['energy_eV']-m['groups']['RTX']['Ca']['energy_eV'])*EV_TO_KCAL)
    for group,g in m['groups'].items():
        ca=g['Ca']['energy_eV'];la=rows[group+'_La_primary']['energy_eV'];n=g['metal_count'];single=sum(s['single_La_energy_eV'] for s in g['ordered_sites'])
        joint=((ca-la)/n-atom)*EV_TO_KCAL;mean=((n*ca-single)/n-atom)*EV_TO_KCAL;k=la+(n-1)*ca-single;old=float(np.mean([s['old_score_model_kcal'] for s in g['ordered_sites']]))
        check(group+'_single_mean_replay',mean-old,1e-6);check(group+'_nonadditivity_identity',joint-mean+k/n*EV_TO_KCAL,1e-6)
        scores[group]={'site_count':n,'all_Ca_eV':ca,'all_La_eV':la,'sum_single_La_eV':single,'atom_reference_eV':atom,'joint_per_metal_model_kcal':joint,
                       'single_substitution_mean_model_kcal':mean,'change_model_kcal':joint-mean,'nonadditivity_eV':k,'ordered_single_site_scores':[s['old_score_model_kcal'] for s in g['ordered_sites']],
                       'evidence':g['evidence'],'calibrated_class':None}
    contrasts=[];parv=[]
    for group in ('A0A7','HEW5','RTX','PARV_4CPV'):
        for gid,value in m['ggr_scores'].items():
            s=scores[group];row={'positive_group':group,'GGR':gid,'old_mean_margin_model_kcal':s['single_substitution_mean_model_kcal']-value,
                               'joint_margin_model_kcal':s['joint_per_metal_model_kcal']-value,'pass':s['joint_per_metal_model_kcal']-value>.02}
            (parv if group=='PARV_4CPV' else contrasts).append(row)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);sp=out/'report_source.py';sp.write_bytes(Path(__file__).read_bytes())
    result={'status':'complete','protocol_id':POLICY,'collection':c,'scores':scores,'checks':checks,'numerical_gate_pass':all(x['pass'] for x in checks),
            'domain_comparisons':contrasts,'domain_pass_count':sum(x['pass'] for x in contrasts),'domain_denominator':9,'domain_biological_groups':3,'GGR_biological_groups':1,
            'parvalbumin_supporting_comparisons':parv,'aequorin_site_labels':None,'assay_cooperativity_reproduced':False,'affinity_or_free_energy':None,'baseline_changed':False,
            'implementation':record(sp),'new_MACE_calls':8,'new_DFT_calls':0,'wall_seconds':time.monotonic()-start}
    write_new(out/'result.json',result);return {k:result[k] for k in ('status','numerical_gate_pass','domain_pass_count','domain_denominator','scores','domain_comparisons','parvalbumin_supporting_comparisons')}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    for op,keys in [('prepare',('config','agreement','output')),('validate',('manifest',)),('collect',('manifest',)),('report',('manifest','output'))]:
        p=sub.add_parser(op)
        for key in keys:p.add_argument('--'+key,required=True)
    a=vars(parser.parse_args());op=a.pop('command');print(json.dumps(globals()[op](**a),indent=2))
