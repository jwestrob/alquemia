"""Saved-only independent dry-body replay; no producer arithmetic or native calls."""
import hashlib
import json
import math
import mmap
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[6]
D = ROOT / 'diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/force_direction_diagnostic_v1'
W = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/force_direction_diagnostic_v1'

def pin(p):
    p=Path(p).resolve(); h=hashlib.sha256()
    with p.open('rb') as f:
        for c in iter(lambda:f.read(1048576),b''): h.update(c)
    return {'path':str(p),'sha256':h.hexdigest()}

def selected_only(p):
    # Balanced JSON scan avoids loading the historical fit trajectory.
    with Path(p).open('rb') as f, mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as m:
        key=b'\n  "selected": '; i=m.find(key)
        assert i>=0
        start=i+len(key); assert m[start]==ord('{')
        depth=0; string=False; escape=False
        for j in range(start,len(m)):
            c=m[j]
            if string:
                if escape: escape=False
                elif c==92: escape=True
                elif c==34: string=False
            elif c==34: string=True
            elif c in (123,91): depth+=1
            elif c in (125,93):
                depth-=1
                if depth==0: return json.loads(m[start:j+1])
    raise RuntimeError('Selected object absent/incomplete')

def dot(a,b): return math.fsum(x*y for x,y in zip(a,b))
def minus(a,b): return [x-y for x,y in zip(a,b)]
def body_sum(x,f,ids):
    center=[math.fsum(x[i][k] for i in ids)/len(ids) for k in range(3)]
    net=[math.fsum(f[i][k] for i in ids) for k in range(3)]
    torque=[math.fsum((x[i][(k+1)%3]-center[(k+1)%3])*f[i][(k+2)%3]
                     -(x[i][(k+2)%3]-center[(k+2)%3])*f[i][(k+1)%3] for i in ids) for k in range(3)]
    return net,torque

def main():
    start=time.monotonic(); cpu=time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    p=W/'run_v1/RESULT.json'; result=json.loads(p.read_text())
    assert result['status']=='saved_force_direction_decomposition_complete'
    for key in ('implementation','declaration','fit','training'): assert pin(result[key]['path'])==result[key]
    selected=selected_only(result['fit']['path']); data=json.loads(Path(result['training']['path']).read_text())
    assert selected['status']=='complete' and len(result['fixtures'])==36
    assert len({(r['kind'],r['fixture']) for r in result['fixtures']})==36
    rows=[]; errors=[]; fraction_error=[]
    for item in result['fixtures']:
        if item['environment']!='dry': continue
        name=item['fixture']; kind=item['kind']; row=data['rows'][name+'_La']
        assert row['geometry']['positions_A']==data['rows'][name+'_Dy']['geometry']['positions_A']
        assert sorted(b['category'] for b in row['bodies'])==['ion','ligand']
        axis=row['geometry']['radial_direction']; assert abs(dot(axis,axis)-1)<1e-10
        # Sum physical model/reference responses first; subtract afterwards.
        contributions=dict(radial_translation=0.,transverse_translation=0.,rotation=0.)
        raw=[]; signed=[]
        for body in sorted(row['bodies'],key=lambda b:b['category']):
            all_res={}
            for metal in ('La','Dy'):
                mf=selected['model'][name+'_'+metal]['forces']; rf=data['rows'][name+'_'+metal]['reference']['forces']
                mn,mt=body_sum(row['geometry']['positions_A'],mf,body['indices'])
                rn,rt=body_sum(row['geometry']['positions_A'],rf,body['indices'])
                all_res[metal]=(minus(mn,rn),minus(mt,rt))
            if kind=='paired': net,torque=(minus(all_res['Dy'][k],all_res['La'][k]) for k in (0,1))
            else: net,torque=all_res[kind]
            radial=dot(net,axis); trans=minus(net,[radial*q for q in axis])
            w=1/6 if body['category']=='ion' else 1/12
            contributions['radial_translation']+=w*radial**2
            contributions['transverse_translation']+=w*dot(trans,trans)
            raw.extend(net)
            if body['category']=='ligand':
                contributions['rotation']+=w*dot(torque,torque);raw.extend(torque)
                signed=[radial,math.sqrt(dot(trans,trans)),math.sqrt(dot(torque,torque))]
        saved=selected['score']['per_fixture'][kind][name]
        errors.extend(abs(a-b) for a,b in zip(raw,saved['raw_body_components_kcal_mol_A']))
        errors.append(abs(math.fsum(contributions.values())-saved['rigid_RMS_kcal_mol_A']**2))
        errors.extend(abs(v-item['RMS_squared_contributions'][k]) for k,v in contributions.items())
        assert item['saved_passed']==saved['passed']==(saved['rigid_RMS_kcal_mol_A']<=2 and max(map(abs,raw))<=5)
        rows.append({'kind':kind,'fixture':name,'contributions':contributions,'passed':saved['passed'],'ligand_radial_transverse_torque':signed})
    assert len(rows)==9 and sum(not r['passed'] for r in rows)==8 and max(errors)<1e-8
    aggregates={}
    for kind in ('La','Dy','paired'):
        rr=[r for r in rows if r['kind']==kind]; c={k:math.fsum(r['contributions'][k] for r in rr)/3 for k in rows[0]['contributions']}
        frac={k:v/math.fsum(c.values()) for k,v in c.items()}
        for k,v in frac.items(): fraction_error.append(abs(v-result['aggregates'][kind+'/dry']['fraction_of_squared_error'][k]))
        aggregates[kind]={'fractions':frac,'failed':sum(not r['passed'] for r in rr),'pooled_RMS':math.sqrt(math.fsum(c.values()))}
    assert max(fraction_error)<1e-12
    out=W/'independent_review_v1';out.mkdir(exist_ok=False)
    review={'status':'passed_saved_force_direction_review','implementation':pin(__file__),'result':pin(p),'producer':result['implementation'],'declaration':result['declaration'],'fit':result['fit'],'training':result['training'],'source_review':'Full36 production loop preserves existing category/body/translation/rotation weights and component signs. Independent arithmetic limited to all9 dry fixtures, not a rerun of full native admission.','dry_rows':rows,'dry_aggregates':aggregates,'max_independent_numeric_difference':max(errors),'max_fraction_difference':max(fraction_error),'interpretation':'Eight of nine dry donor force fixtures fail original gates; La dry squared error is mainly rotational and paired dry squared error mainly radial. A response with exact no-water bypass and water-O-only short-range changes cannot alter these fixed dry donor forces. This localizes a model/reference discrepancy, not a unique donor-physics cause; dry reference directional derivative qualification remains separate.','new_native_calls':0,'held_outcomes_read':False,'threads':1,'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}}
    (out/'REVIEW.json').write_text(json.dumps(review,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'review':pin(out/'REVIEW.json'),'max_difference':max(errors),'dry_aggregates':aggregates,'cost':review['cost']},indent=2))

if __name__=='__main__':main()
