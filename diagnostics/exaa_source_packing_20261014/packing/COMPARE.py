#!/usr/bin/env python3
"""Compare fixed-pool packing works using saved outputs; Slurm only."""
import argparse, hashlib, json, math, os, resource, socket, time, traceback
from pathlib import Path

HERE=Path(__file__).resolve().parent
ORDER=('origin','adaptive_Ca','adaptive_La'); METALS=('Ca','La')
TERMS=('Gelec','Ghb','Gsasa','Gshift')
EV=23.06054783061903; HA=627.509474; TOL=1e-8

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,('changed input',p);return Path(p['path'])
def read(p):return json.loads((verify(p) if isinstance(p,dict) else Path(p)).read_text())
def write(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def close(a,b,residuals):
    assert math.isfinite(a) and math.isfinite(b)
    d=abs(a-b);residuals.append(d);assert d<=TOL,('algebra mismatch',a,b,d)
def difference(a,b):return {k:a[k]-b[k] for k in a}

def freeze(a):
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False)
    prep=read(a.preparation);assert prep['status']=='prepared' and len(prep['rows'])==2
    old=read(a.old_scores);new=Path(a.new_collection).resolve();cfg=new.parent.parent/'CONFIG.json'
    c={'implementation':pin(__file__),'plan':pin(HERE/'COMPARE_PLAN.md'),
       'preparation':pin(a.preparation),'old_scores':pin(a.old_scores),'old_review':old['input'],
       'new_collection_path':str(new),'new_config':pin(cfg),'output':str(out/'execution'),
       'denominator':2,'new_energy_calls':0,'algebra_tolerance_kcal_mol':TOL,
       'policy':{'candidate_order':ORDER,'origin_tolerance_kcal_mol':.1,'EV_TO_KCAL':EV,'HA_TO_KCAL':HA}}
    write(out/'CONFIG.json',c);print(json.dumps(pin(out/'CONFIG.json')))

def cell(saved,expected,state,old=False):
    """Resolve exact native + MACE evidence; return only internal raw values."""
    if saved is None:raise ValueError('required matrix cell unavailable')
    native_pin=saved['native_result'] if old else saved['solvent'];native=read(native_pin)
    assert native['status']=='available';t=native['input']
    assert t['xyz']['sha256']==expected['xyz']['sha256'];verify(t['xyz']);verify(expected['xyz'])
    assert t['state_signature']==state;verify(state)
    assert t['charge']==expected['charge'] and t['multiplicity']==expected['multiplicity']
    mp=expected['MACE'] if old else saved['MACE'];mm=read(mp);rq=read(mm['request'])
    assert mm['status'] in ('complete','computed')
    assert rq['xyz']['sha256']==expected['xyz']['sha256'];verify(rq['xyz'])
    assert rq['charge']==expected['charge'] and rq['multiplicity']==expected['multiplicity']
    assert mm['energy_eV']==saved['MACE_eV']
    if old:assert mm['energy_eV']==expected['MACE_eV']
    assert any(u['MACE']==mp and u['MACE_eV']==mm['energy_eV'] for u in t['uses'])
    assert native['G_qvac_hartree']==saved['G_qvac_hartree']
    assert native['vacuum']['charge_source']=='vacuum'
    assert native['vacuum']['G_Eh']==native['G_qvac_hartree']
    v=native['vacuum']['components_Eh'];assert set(v)==set(TERMS)
    assert abs(sum(v.values())-native['G_qvac_hartree'])<=1e-12
    if old:assert v==saved['native_components_qvac_hartree']
    return {'MACE_eV':mm['energy_eV'],'solvent_Eh':native['G_qvac_hartree'],'components_Eh':v,
            'evidence':{'MACE':mp,'request':mm['request'],'native':native_pin,'xyz':expected['xyz'],'state':state}}

def works(matrix):
    out={z:{} for z in METALS}
    for z in METALS:
        for q in ORDER:
            v,o=matrix[z].get(q),matrix[z].get('origin')
            if v is None or o is None:out[z][q]=None;continue
            dm=(v['MACE_eV']-o['MACE_eV'])*EV;ds=(v['solvent_Eh']-o['solvent_Eh'])*HA
            out[z][q]={'MACE_kcal_mol':dm,'solvent_kcal_mol':ds,'total_kcal_mol':dm+ds,
                      **{k+'_kcal_mol':(v['components_Eh'][k]-o['components_Eh'][k])*HA for k in TERMS}}
    complete=all(out[z][q] is not None for z in METALS for q in ORDER)
    selected=None
    if complete:
        selected={}
        for policy in ('operational','mathematical'):
            s={}
            for z in METALS:
                best=min(ORDER,key=lambda q:out[z][q]['total_kcal_mol'])
                q=best if policy=='mathematical' or out[z][best]['total_kcal_mol']<-.1 else 'origin'
                s[z]={'candidate':q,'work':out[z][q]}
            s['differential_work']=difference(s['Ca']['work'],s['La']['work']);selected[policy]=s
    paired={q:difference(out['Ca'][q],out['La'][q]) if all(out[z][q] is not None for z in METALS) else None for q in ORDER}
    return {'status':'available' if complete else 'partial_or_unavailable','works':out,'selected':selected,'same_candidate_differential_work':paired}

def run(path):
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_JOB_NODELIST'),'actual Slurm required'
    assert int(os.environ['SLURM_CPUS_ON_NODE'])==1
    c=read(path);assert c['implementation']==pin(__file__);read(c['old_review']);verify(c['plan'])
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    receipt={'status':'failed','config':pin(path),'job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'new_molecular_calls':0}
    try:
        prep=read(c['preparation']);old=read(c['old_scores']);review=read(c['old_review']);nc=read(c['new_config'])
        assert nc['input']==c['preparation'] and len(nc['cells'])==12
        np=Path(c['new_collection_path']);new=read(np) if np.is_file() else None
        if new:assert new['config']==c['new_config'] and new['denominator']==2
        oldrows={r['record_id']:r for r in old['rows']};refs={r['record_id']:r for r in review['cases']}
        newrows={r['case_id']:r for r in new['rows']} if new else {}
        assert not set(newrows)-{p['case_id'] for p in prep['rows']}
        rows=[];residuals=[]
        for p in prep['rows']:
            rid=p['original_review_identity']['record_id'];o=oldrows[rid];ref=refs[rid]
            assert o['case_key']==ref['case_key']==p['original_review_identity']['case_key']
            assert o['source_case_id']==p['case_id'] and ref['record']==p['original_record']
            assert o['status']=='available';n=newrows.get(p['case_id']);fail=[];matrices={}
            for context in ('small','expanded'):
                mat={z:{} for z in METALS}
                for z in METALS:
                    for q in ORDER:
                        try:
                            if context=='small':
                                e=dict(p['old_matrix'][z][q],charge=o['state']['endpoint_charges'][z],multiplicity=o['state']['multiplicities'][z])
                                v=cell(o['matrix'][z][q],e,p['original_record']['state_signature'],True)
                            else:
                                if n is None:raise ValueError('new collection/source missing')
                                e=p['matrix'][z][q];v=cell(n['matrix'][z][q],e,p['state_signature'])
                            mat[z][q]=v
                        except Exception as exc:
                            mat[z][q]=None;fail.append({'context':context,'metal':z,'candidate':q,'reason':type(exc).__name__+': '+str(exc)})
                matrices[context]=mat
            small=works(matrices['small']);large=works(matrices['expanded'])
            if small['status']=='available':
                for z in METALS:
                    for q in ORDER:close(small['works'][z][q]['total_kcal_mol'],o['works'][z][q],residuals)
                for policy in ('operational','mathematical'):
                    s=small['selected'][policy];expected=ref['score']['pool'][policy]
                    for z in METALS:
                        assert s[z]['candidate']==expected[z]['candidate']
                        close(s[z]['work']['total_kcal_mol'],expected[z]['work_kcal_mol'],residuals)
                    close(s['differential_work']['total_kcal_mol'],expected['R_model_kcal_mol']-ref['score']['pool']['origin']['R_model_kcal_mol'],residuals)
                assert {z:small['selected']['operational'][z]['candidate'] for z in METALS}==o['geometry_choices']
                close(small['selected']['operational']['differential_work']['total_kcal_mol'],o['accommodation_R_shift'],residuals)
            if large['status']=='available':
                assert n['status']=='available'
                for z in METALS:
                    for q in ORDER:
                        for k,v in n['works'][z][q].items():close(large['works'][z][q][k],v,residuals)
                for policy in ('operational','mathematical'):
                    for z in METALS:assert large['selected'][policy][z]['candidate']==n['selected'][policy][z]['candidate']
                    close(large['selected'][policy]['differential_work']['total_kcal_mol'],n['selected'][policy]['differential_work_kcal_mol'],residuals)
            changes={z:{q:difference(large['works'][z][q],small['works'][z][q]) if large['works'][z][q] is not None and small['works'][z][q] is not None else None for q in ORDER} for z in METALS}
            paired={q:difference(large['same_candidate_differential_work'][q],small['same_candidate_differential_work'][q]) if large['same_candidate_differential_work'][q] is not None and small['same_candidate_differential_work'][q] is not None else None for q in ORDER}
            selections=None
            if small['selected'] and large['selected']:
                selections={policy:{'per_metal':{z:difference(large['selected'][policy][z]['work'],small['selected'][policy][z]['work']) for z in METALS},
                    'differential_work':difference(large['selected'][policy]['differential_work'],small['selected'][policy]['differential_work'])} for policy in ('operational','mathematical')}
            evidence={cx:{z:{q:v['evidence'] if v else None for q,v in m[z].items()} for z in METALS} for cx,m in matrices.items()}
            rows.append({'case_id':p['case_id'],'original_review_identity':p['original_review_identity'],'role':p['role'],
               'status':'available' if not fail else 'partial_or_unavailable','failures':fail,'new_collection_failures':n.get('failures',[]) if n else [],
               'small':small,'expanded':large,'expanded_minus_small_works':changes,'expanded_minus_small_same_candidate_differential_work':paired,
               'expanded_minus_small_selected_works':selections,'evidence':evidence})
        result={'denominator':2,'rows':rows,'config':pin(path),'new_collection':pin(np) if np.is_file() else None,
          'status':'complete' if all(r['status']=='available' for r in rows) else 'partial_or_unavailable',
          'actual_saved_algebra_checks':{'tolerance_kcal_mol':TOL,'count':len(residuals),'maximum_residual_kcal_mol':max(residuals,default=None)},
          'absolute_crosssize_energy_comparisons':False,'classifications':None,'new_molecular_calls':0,'production_changed':False}
        write(out/'RESULT.json',result)
        lines=['# Fixed-pool packing accommodation comparison','','All values below are works relative to each context/metal origin, in kcal/mol. No calibrated class calls.','',
               '| Source | Metal | Candidate | Small total | Expanded total | Expanded − small |','|---|---|---|---:|---:|---:|']
        fmt=lambda v:'unavailable' if v is None else f'{v["total_kcal_mol"]:.6f}'
        for r in rows:
            for z in METALS:
                for q in ORDER:lines.append(f'| {r["case_id"].rsplit("sample-",1)[-1]} | {z} | {q} | {fmt(r["small"]["works"][z][q])} | {fmt(r["expanded"]["works"][z][q])} | {fmt(r["expanded_minus_small_works"][z][q])} |')
        lines+=['','Full component, selection-policy and unavailable-cell details are in RESULT.json.']
        (out/'REPORT.md').write_text('\n'.join(lines)+'\n');receipt.update(status=result['status'],result=pin(out/'RESULT.json'))
    except BaseException as exc:receipt.update(reason=str(exc),traceback=traceback.format_exc());raise
    finally:
        u=resource.getrusage(resource.RUSAGE_SELF);receipt.update(wall_seconds=time.monotonic()-start,CPU_seconds=u.ru_utime+u.ru_stime,peak_RSS_KiB=u.ru_maxrss)
        write(out/'EXECUTION.json',receipt)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    f=s.add_parser('freeze');f.add_argument('--new-collection',required=True);f.add_argument('--old-scores',required=True);f.add_argument('--preparation',required=True);f.add_argument('--output',required=True)
    f=s.add_parser('run');f.add_argument('--config',required=True);a=p.parse_args()
    freeze(a) if a.action=='freeze' else run(a.config)
