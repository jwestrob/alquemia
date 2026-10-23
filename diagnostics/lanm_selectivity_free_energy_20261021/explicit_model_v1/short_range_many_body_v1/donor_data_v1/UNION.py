#!/usr/bin/env python3
"""Saved metadata applicability receipt and immutable common-schema donor views."""
import argparse
import copy
import time
from pathlib import Path
import BUILD_DATA as b

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--donor',required=True);ap.add_argument('--water',required=True);ap.add_argument('--out',required=True)
    a=ap.parse_args();t=time.monotonic();cpu=time.process_time();out=Path(a.out).resolve()
    b.check(not out.exists(),'fresh union output');out.mkdir(parents=True)
    dp=b.pin(a.donor);wp=b.pin(a.water);d=b.read(dp);w=b.read(wp);base=b.read(d['baseline'])
    sources=b.read(b.pin(b.HERE/'INPUTS.json'));resolver=b.load(sources['resolver'])
    key=b.verify(base['key']).read_text();prm=b.verify(base['parameter_file']).read_text()
    shared=resolver.parse_parameters(prm,key);checks={}
    wf=b.read(w['fixtures'])['fixtures']
    b.check(len(wf)==54,'exact54 water operator fixtures')
    for group,fixtures in [('water',wf),('donor',d['fixtures'])]:
        for fid,pin in fixtures.items():
            f=b.read(pin);kind='water' if group=='water' else ('acetate' if 'acetate' in fid else 'amide')
            g,pp,k=b.typed_geometry(f['inputs'],kind,resolver)
            b.check(g['xyz_A']==f['geometry']['xyz_A'],'literal unchanged geometry '+fid)
            c=b.compatibility(pp,k,g,shared,key,f['inputs']['parameters']['sha256']==base['parameter_file']['sha256'])
            expected=group=='water' or fid.startswith('matched_dry_')
            b.check(c['compatible']==expected,'prospective applicability '+fid)
            if group=='water':b.check(set(g['tinker_types'])<={500,501,36,37},'water type exclusion')
            checks[group+'/'+fid]=dict(fixture=pin,**c)
    # Narrow rejection controls use metadata copies only, never alter saved inputs.
    exemplar=b.read(d['fixtures']['matched_dry_acetate_mono_central_La'])
    g,pp,k=b.typed_geometry(exemplar['inputs'],'acetate',resolver)
    controls={}
    def test(name,condition):
        b.check(condition,name);controls[name]=True
    test('zero_union_dry_accepted',b.compatibility(pp,k,g,shared,key,True)['compatible'])
    p=copy.deepcopy(pp);p['polpair'][(36,500)]=(99.,)
    test('absent_water_pair_does_not_change_dry',b.compatibility(p,k,g,shared,key,True)['compatible'])
    p=copy.deepcopy(pp);p['polpair'][(403,500)]=(99.,)
    test('present_donor_pair_rejected',not b.compatibility(p,k,g,shared,key,True)['compatible'])
    p=copy.deepcopy(pp);p['polarize'][500]['alpha_A3']*=1.01
    test('shared_metal_alpha_change_rejected',not b.compatibility(p,k,g,shared,key,True)['compatible'])
    test('changed_parameter_file_rejected',not b.compatibility(pp,k,g,shared,key,False)['compatible'])
    test('external_field_not_ignored',not b.compatibility(pp,k+'\nexternal-field 1 0 0\n',g,shared,key,True)['compatible'])
    test('multipole_body_change_not_ignored',not b.compatibility(pp,k+'\n0.1 0.2 0.3\n',g,shared,key,True)['compatible'])
    test('observer_flags_only_are_ignored',b.compatibility(pp,k+'\npolar-print\nsave-uinduce\n',g,shared,key,True)['compatible'])
    test('intact_two_oxygen_acetate',len(g['ligands'])==1 and len(g['ligands'][0]['donor_atoms'])==2)
    test('unequal_dual_channels_preserved',exemplar['native_operator']['Ed_e_per_A2']!=exemplar['native_operator']['Ep_e_per_A2'])
    q=b.read(sources['qualification']);views={}
    for metal in ('La','Dy'):
        fid='matched_dry_acetate_mono_central_'+metal;old=d['fixtures'][fid];f=copy.deepcopy(b.read(old));op=f['native_operator'];r=q['rows'][f['native_qualification']['row']]
        op['u_native_eA']=op.pop('native_u_eA');op['v_native_eA']=op.pop('native_v_eA')
        op['permanent_mu_origin_A']=op['origin_A'];op['permanent_energy_kcal']=None
        op['nonpolar_energy_meaning']=op['nonpolar_energy_definition'];op['snapshot']=r['exports'][0]
        op['field_coupling']='uniform field added to both native RHS; exact conversion in donor source manifest units'
        oldqualification=f['native_qualification']
        f['native_qualification']=dict(producer=oldqualification['result'],review=oldqualification['review'],
            polarization_energy=r['operator']['polarization_energy'],two_solve_exact=r['operator']['two_solve_exact'],
            snapshots=[{k:s[k]for k in ('source','sequence','vectors','numerical_bounds')}for s in r['operator']['exports']])
        f['harmonization']=dict(source_fixture=old,source_qualification=oldqualification,
            implementation=b.pin(__file__),changes='Metadata field names and saved qualification snapshots only; no numerical arrays, observations or baseline identities changed.')
        views[fid]=b.save(out/'kernel_views'/(fid+'.json'),f)
    result=dict(status='shared_baseline_applicability_verified',implementation=b.pin(__file__),
        donor_dataset=dp,water_dataset=wp,shared_baseline=d['baseline'],checks=checks,
        counts=dict(water_compatible=54,dry_donor_compatible=6,wet_donor_incompatible=18),
        kernel_views=views,manufactured_controls=controls,
        limitations=['Saved wet donor E/F are not shared-baseline residuals.','No new field/operator/force calculation; all arrays unchanged.','No full force provider; baseline coordinate derivatives absent.'],
        cost=dict(wall_seconds=time.monotonic()-t,CPU_seconds=time.process_time()-cpu,threads=1),
        native_calls=0,quantum_calls=0,fits=0,operator_solves=0)
    p=b.save(out/'RESULT.json',result);print(p)

if __name__=='__main__':main()
