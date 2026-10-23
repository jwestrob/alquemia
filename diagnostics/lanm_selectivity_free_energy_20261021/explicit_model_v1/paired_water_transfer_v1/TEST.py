#!/usr/bin/env python3
"""Saved-native parser and physical arithmetic controls. Zero native calls."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('transfer_test_run',HERE/'RUN.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r)

def main(data_path,out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);data=r.loaded(data_path);cfg=data['config'];body=r.module(cfg['body_helper'],'transfer_test_body');score=r.module(r.pin(HERE/'SCORE.py'),'transfer_test_score');native=r.module(cfg['native_helper'],'transfer_test_native')
    records=[]
    def check(ok,name):
        assert ok,name;records.append({'case':name,'passed':True})
    refs={k:copy.deepcopy(v['reference'])for k,v in data['rows'].items()};base={k:copy.deepcopy(v['baseline'])for k,v in data['rows'].items()}
    a=score.evaluate(data,refs,body);check(a['all_development_checks_passed'],'Reference self-target zero residual')
    check(len(a['radial_work'])==24 and len(a['balanced_exchange'])==12 and len(a['matched_water_addition'])==3 and len(a['rigid_response'])==48,'Every individual/paired finite denominator')
    b=score.evaluate(data,base,body);check(all(x['absolute_error_change_kcal_mol']==0 and not x['deteriorated']for k in ('radial_work','balanced_exchange','matched_water_addition')for x in b[k]),'Unchanged baseline energy comparison')
    check(all(x['RMS_error_change']==0 and x['maximum_error_change']==0 and not x['deteriorated']for x in b['rigid_response']),'Unchanged baseline rigid comparison')
    shifted=copy.deepcopy(refs)
    for k,v in shifted.items():v['energy']+=111 if k.endswith('_La') else 222
    c=score.evaluate(data,shifted,body);check(all(abs(x['errors']['candidate']['signed_kcal_mol'])<1e-8 for k in ('radial_work','balanced_exchange','matched_water_addition')for x in c[k]),'Constant element offsets cancel every balanced cycle')
    changed=copy.deepcopy(refs);changed['water_plus_Dy']['energy']+=2;c=score.evaluate(data,changed,body)
    check(abs(next(x for x in c['radial_work']if x['id']=='water_plus_Dy_minus_La')['errors']['candidate']['signed_kcal_mol']-2)<1e-8,'Paired radial work retains Dy-minus-La sign')
    check(not c['development_groups']['water_radial']['candidate_passed'] and c['development_groups']['hydrated_donor']['candidate_passed'],'Individual failure not averaged across families')
    changed=copy.deepcopy(refs);gid='hydrated_acetamide_central_La';row=data['rows'][gid];ids=next(b['indices']for b in row['bodies']if b['category']=='ligand');i,j=ids[:2];v=[row['geometry']['positions_A'][i][k]-row['geometry']['positions_A'][j][k]for k in range(3)]
    for k in range(3):changed[gid]['forces'][i][k]+=v[k];changed[gid]['forces'][j][k]-=v[k]
    c=score.evaluate(data,changed,body);x=next(x for x in c['rigid_response']if x['geometry']=='hydrated_acetamide_central' and x['kind']=='La')
    check(x['models']['candidate']['rigid_RMS_kcal_mol_A']<1e-12 and x['models']['candidate']['atom_vector_RMS_kcal_mol_A']>0,'Internal central force cancels rigid response but remains diagnostic')
    for family in ('water','hydrated_acetamide'):
        synthetic=copy.deepcopy(refs);p=data['paths'][family+'_La'];ids=p['moving_atoms'];u=p['direction']
        for side in ('minus','central','plus'):
            k=family+'_'+side+'_La';synthetic[k]['forces']=[[0.,0.,0.]for _ in synthetic[k]['forces']]
            for i in ids:synthetic[k]['forces'][i]=[-v/len(ids)for v in u]
        c=score.evaluate(data,synthetic,body);x=next(v for v in c['coarse_radial_diagnostics']if v['family']==family and v['metal']=='La')
        check(all(abs(v-1)<1e-12 for v in x['models']['candidate']['slopes_minus_central_plus_kcal_mol_A']),'Outward coordinate derivative sign '+family)
    check(len(r.d.CENTERS)+2*(len(r.d.RADIAL)+len(r.d.HYDRATED))==48 and sum(6*data['rows'][k]['geometry']['atoms']for k in r.d.FD)==1380,'48 program and1380 internal FD budgets')
    for gid in r.d.CELLS:
        key,render=r.d.render(cfg,data['rows'][gid],list(r.d.e.INITIAL));check(render['original_key_reused'] and key==r.verify(data['rows'][gid]['input']['key']).read_text(),'Original parameter exact replay '+gid)
    baseline=data['rows']['water_central_La']['baseline']['native'];gpath=r.verify(baseline['gradient_output']);g=r.gradient(gpath,4,True,native)
    check(g['fd_passed'] and g['fd_total_separation_A']==1e-5 and g['fd_each_side_displacement_A']==5e-6,'Actual native Cartesian output and FD convention')
    text=gpath.read_text()
    variants={'missing_numerical_atom':re.sub(r'^\s*Numer\s+4\s+.*$','',text,flags=re.M),'polarization_failure':text+'\nInduced dipoles not converged\n','nonfinite':text+'\nNaN\n','missing_energy':text.replace('Total Potential Energy','Missing Energy')}
    for label,content in variants.items():
        p=out/(label+'.synthetic.out');p.write_text(content)
        try:r.gradient(p,4,True,native)
        except (AssertionError,ValueError):check(True,'Reject '+label)
        else:raise AssertionError('Accepted '+label)
    p=out/'bad_fd.synthetic.out';lines=[]
    for line in text.splitlines():
        if re.match(r'^\s*Numer\s+1\s+',line):
            fields=line.split();fields[2]=str(float(fields[2])+.1);line=' '.join(fields)
        lines.append(line)
    p.write_text('\n'.join(lines)+'\n');check(not r.gradient(p,4,True,native)['fd_passed'],'FD failure remains failed under original tolerance')
    result={'status':'passed','implementation':r.pin(__file__),'run_source':r.pin(HERE/'RUN.py'),'data_source':r.pin(HERE/'DATA.py'),'score_source':r.pin(HERE/'SCORE.py'),'prepared_data':r.pin(data_path),'records':records,'new_native_calls':0,'reserved_validation_opened':False,'elapsed_seconds':time.monotonic()-start}
    r.save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'checks':len(records),'result':r.pin(out/'RESULT.json')}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.data,a.out)
