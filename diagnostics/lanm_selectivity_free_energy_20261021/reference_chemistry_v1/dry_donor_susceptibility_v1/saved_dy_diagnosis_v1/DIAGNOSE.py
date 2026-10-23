#!/usr/bin/env python3
"""Five saved Dy outputs: stage attribution, not corrected target admission."""
import argparse, hashlib, json, math, os, re, time
from decimal import Decimal, getcontext
from pathlib import Path
getcontext().prec=50
D=Decimal
TAGS=('zero','minus_h','minus_half','plus_half','plus_h')

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def checked(q):
    data=Path(q['path']).read_bytes();assert hashlib.sha256(data).hexdigest()==q['sha256'];return data.decode()

def bound(s):return D(5)*D(10)**(D(s).as_tuple().exponent-1)

def plain(x):
    if isinstance(x,D):return {'decimal':str(x),'float':float(x)}
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,list):return [plain(v) for v in x]
    return x

def main():
    p=argparse.ArgumentParser();p.add_argument('--collection',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    tick,cpu=time.monotonic(),time.process_time();a.out.mkdir(parents=True,exist_ok=False)
    collection=json.loads(a.collection.read_text());assert collection['status']=='dry_donor_susceptibility_reference_unqualified' and collection['complete']==10 and collection['response']['declared_checks']==16
    direction=[D(str(x)) for x in collection['response']['field_direction_lab']];rows={};raw={}
    for tag in TAGS:
        cid='pbe0__Dy__dry_mono_central__Ft_'+tag;q=collection['rows'][cid]
        assert q['status']=='complete' and q['returncode']==0
        text=checked(q['output']);lines=text.splitlines();inp=checked(q['input']);prop=checked(q['properties']['property_file'])
        assert all(t in prop for group in q['properties']['raw_tokens']['dipoles'].values() for t in group)
        assert q['properties']['raw_tokens']['electronic'] in prop
        header=[(i+1,l) for i,l in enumerate(lines) if 'ITER' in l and 'Energy' in l]
        end=[i for i,l in enumerate(lines) if 'SCF CONVERGED AFTER' in l];assert len(end)==1;end=end[0]
        itrows=re.findall(r'^\s*(\d+)\s+(-?\d+\.\d{10,})\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*$', '\n'.join(lines[:end]),re.M)
        assert itrows;last=itrows[-1];epre=D(last[1]);posttok=q['properties']['raw_tokens']['electronic'];epost=D(posttok)
        native_e=re.findall(r'^\s*Total Energy\s*:\s*([-\d.]+)\s+Eh',text,re.M)[-1];assert abs(epost-D(native_e))<=bound(posttok)+bound(native_e)
        corr=re.findall(r'Exchange energy change after final integration\s*:\s*([-\d.]+)\s+Eh',text);assert len(corr)==1
        delta=epost-epre;assert abs(delta-D(corr[0]))<=bound(last[1])+bound(posttok)+bound(corr[0])
        mu=sum(t*D(x) for t,x in zip(direction,q['properties']['raw_tokens']['dipoles']['total']))
        d4tok=q['properties']['raw_tokens']['D4'];finaltok=q['properties']['raw_tokens']['final'];closure=D(finaltok)-epost-D(d4tok)
        assert abs(closure)<=D(str(q['properties']['energy_stage_closure_bound_hartree']))
        normalized='\n'.join(l for l in inp.splitlines() if not re.match(r'^\s*EField\s',l,re.I))
        orb=[];at=text.find('ORBITAL ENERGIES')
        for match in re.finditer(r'^\s*(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$',text[at:],re.M):
            index,occ,eh,ev=match.groups();orb.append({'index':int(index),'occupation':float(occ),'hartree':float(eh),'eV':float(ev)})
        occ=[x for x in orb if x['occupation']>0];virt=[x for x in orb if x['occupation']==0]
        orbital_edges={'HOMO':max(occ,key=lambda x:x['index']),'LUMO':min(virt,key=lambda x:x['index'])} if occ and virt else None
        row={'cell_id':cid,'result':q['result'],'output':q['output'],'property':q['properties']['property_file'],'input':q['input'],
             'normalized_input_sha256':hashlib.sha256(normalized.encode()).hexdigest(),'field_au':q['cell']['field_au'],
             'charge':q['cell']['charge'],'electrons':q['cell']['expected_electrons'],'multiplicity':q['cell']['multiplicity'],'physical_parent_multiplicity':q['cell']['physical_multiplicity_of_parent'],
             'property_origin_Bohr':q['properties']['origin_Bohr'],'scf_iterations':int(last[0]),'SCF_iteration_headers':header,'last_iteration_row':last,
             'energy_pre_final_hartree':epre,'energy_final_electronic_hartree':epost,'final_integration_delta_hartree':delta,'printed_delta_hartree':corr[0],
             'final_total_hartree':D(finaltok),'D4_hartree':D(d4tok),'electronic_final_D4_closure_hartree':closure,'projected_mu_au':mu,
             'printed_scf_criteria':q['parsed']['printed_scf_criteria'],'native_stop':q['parsed']['native_scf_stop_signals'],
             'native_mode':q['parsed']['native_scf_mode'],'grid_point_counts_XC_COSX1_2_3':list(map(int,re.findall(r'Total number of grid points\s*\.\.\.\s*(\d+)',text))),
             'COSX_settings_signature':q['parsed']['native_energy_stages']['COSX_grid_settings_signature'],'orbital_edges':orbital_edges,
             'warnings':[l.strip() for l in lines if re.search(r'WARNING|Warning:|ConvCheckMode',l)],
             'native_elapsed_seconds':q['elapsed_seconds'],'allocated_CPU_seconds_proxy':q['allocated_elapsed_CPU_seconds_proxy']}
        rows[tag]=row;raw[tag]={'terminal_scf_and_correction':{'first_line':max(0,end-18)+1,'text':'\n'.join(lines[max(0,end-18):end+15])},'criteria_and_stage_evidence':q['parsed']['native_energy_stages']['selected_native_lines']}
    assert len({x['normalized_input_sha256'] for x in rows.values()})==1
    assert len({str(x['grid_point_counts_XC_COSX1_2_3']) for x in rows.values()})==1
    assert len({str(x['COSX_settings_signature']) for x in rows.values()})==1
    assert len({x['D4_hartree'] for x in rows.values()})==1
    assert all(x['property_origin_Bohr']==[0.,0.,0.] for x in rows.values())
    calculations={};zero=rows['zero'];mu0=zero['projected_mu_au']
    for step,h in [('h',D('.001')),('half',D('.0005'))]:
        minus,plus=rows['minus_'+step],rows['plus_'+step]
        oddpost=-(plus['energy_final_electronic_hartree']-minus['energy_final_electronic_hartree'])/(2*h)
        oddpre=-(plus['energy_pre_final_hartree']-minus['energy_pre_final_hartree'])/(2*h)
        correction=oddpost-oddpre;even=plus['projected_mu_au']+minus['projected_mu_au']-2*mu0;beta=even/h**2
        chi={name:-(plus[key]+minus[key]-2*zero[key])/h**2 for name,key in [('post','energy_final_electronic_hartree'),('pre','energy_pre_final_hartree')]}
        check=collection['response']['per_metal']['Dy']['estimates'][step]
        assert abs(float(oddpost)-check['odd_energy_moment'])<3e-10
        calculations[step]={'field_au':h,'mu0':mu0,'odd_final':oddpost,'odd_pre_final_diagnostic':oddpre,'odd_final_minus_mu0':oddpost-mu0,
            'odd_pre_final_minus_mu0_diagnostic':oddpre-mu0,'odd_final_integration_contribution':correction,'fraction_of_odd_mismatch_from_final_integration':correction/(oddpost-mu0),
            'printed_correction_odd_bound':(bound(plus['printed_delta_hartree'])+bound(minus['printed_delta_hartree']))/(2*h),
            'chi_energy_final':chi['post'],'chi_energy_pre_final_diagnostic':chi['pre'],'chi_final_integration_contribution':chi['post']-chi['pre'],
            'chi_dipole':(plus['projected_mu_au']-minus['projected_mu_au'])/(2*h),'even_mu_numerator':even,'beta_from_even_mu_diagnostic':beta,
            'leading_cubic_energy_odd_shift_predicted_from_mu':even/6,'odd_mismatch_after_leading_mu_cubic_diagnostic':oddpost-mu0-even/6,
            'original_odd_gate_limit':D('.0001'),'original_odd_gate_passes':abs(oddpost-mu0)<=D('.0001')}
    scale={'odd_mismatch_h_over_half':calculations['h']['odd_final_minus_mu0']/calculations['half']['odd_final_minus_mu0'],
           'even_moment_h_over_half':calculations['h']['even_mu_numerator']/calculations['half']['even_mu_numerator'],
           'expected_ratio_if_leading_smooth_cubic_term_dominates':4,'meaning':'Smooth leading hyperpolarizability alone does not reproduce observed odd residual scaling; higher-order response and numerical residuals are not uniquely separated.'}
    result={'status':'saved_dy_stage_diagnosis_complete_original_panel_remains_unqualified','implementation':pin(__file__),'collection':pin(a.collection),'rows':rows,'calculations':calculations,'truncation_scaling':scale,
        'La_control_summary_reuse_only':collection['response']['per_metal']['La']['estimates'],
        'native_inputs_equal_except_field':True,'same_reported_XC_and_COSX_grids':True,'D4_field_independent':True,
        'energy_stage_join_error_found':False,'original_checks':collection['response']['checks'],'original_numerical_qualification':False,
        'pre_final_energy_is_diagnostic_only':True,'corrected_target_admitted':False,'unique_cause_identified':False,'new_native_calls':0,'saved_raw_outputs_read':5,
        'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (a.out/'RESULT.json').write_text(json.dumps(plain(result),indent=2)+'\n');(a.out/'RAW_EXCERPTS.json').write_text(json.dumps(raw,indent=2)+'\n')
    print(json.dumps({'result':pin(a.out/'RESULT.json'),'calculations':plain(calculations),'scaling':plain(scale)}))

if __name__=='__main__':main()
