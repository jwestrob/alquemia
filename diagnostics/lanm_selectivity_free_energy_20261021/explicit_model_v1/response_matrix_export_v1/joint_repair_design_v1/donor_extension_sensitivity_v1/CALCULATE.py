"""Saved fixed-probe arithmetic only: no producer imports, scoring, fit or native calls."""
import hashlib
import json
import math
import os
import re
import time
from decimal import Decimal, getcontext
from pathlib import Path
import numpy as np

getcontext().prec = 45
T0, C0 = time.monotonic(), time.process_time()
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
OUT = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/donor_extension_sensitivity_v1/actual_v1'
PINS = {}

def pin(p):
    p = Path(p).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def verify(p):
    assert pin(p['path']) == p, p
    PINS[p['path']] = p['sha256']
    return Path(p['path'])

def read(p): return json.loads(verify(p).read_text())
def txt(p): return verify(p).read_text()
def dec(s): return Decimal(str(s).replace('D','E'))
def ar(dim, *values): return 128*dim*np.finfo(float).eps*max([1.]+[float(np.max(np.abs(x))) for x in values])
def dot(a,b): return math.fsum(float(x)*float(y) for x,y in zip(a,b))

for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:
    assert os.environ.get(key) == '1', key
decl = read(pin(HERE/'DECLARATION.json'))
r = read(decl['qualified121_result']); review = read(decl['qualified121_independent_review'])
m = read(decl['qualified121_manifest']); src = read(decl['design_sources'])
assert r['status']=='donor_extension_numerical_qualification_complete' and r['numerical_checks_passed'] and r['error'] is None
assert r['native_program_invocations']==121 and r['internal_FD_energies_completed']==576
assert review['qualification_passed'] is True
assert r['manifest']==decl['qualified121_manifest']
rp=read(src['donor_authorities']['reference_panel'])
rr=read(src['donor_authorities']['reference_derivative_result'])
inv=read(src['donor_inventory']); model40=read(src['model40_result']); la5=read(src['La_quantum_field_collection'])
assert rr['scoring']['all_numerical_checks_passed'] and rr['admitted_cells']==16
assert model40['summary']['numerical_checks_passed'] and model40['complete_cells']==40
assert la5['complete']==5 and la5['numerical_checks_passed']
for name in ['model40_actual_review','La_quantum_field_actual_review']:
    read(src[name])
read(src['donor_authorities']['reference_actual_review'])
B=m['units']['bohr_A']; Ha=m['units']['hartree_kcal_mol']; K=m['units']['electric']; t=np.array(m['field_direction'])
observables=['energy_kcal_mol','radial_dE_dq_kcal_mol_A','rotation_dE_dq_kcal_mol_A','projected_mu_ea0','chi_a0_cubed']
used_rows={}

def record(cid):
    v=r['rows'][cid]; assert v['status']=='complete' and v['numerical_checks_passed']
    verify(v['output']);verify(v['start'])
    for p in v['runtime_input'].values():verify(p)
    used_rows[cid]=dict(output=v['output'],start=v['start'],input=v['runtime_input'],profile=v['job']['profile'],cell=v['job']['cell'])
    return v

def energy(v):
    tokens=re.findall(r'^\s*Total Potential Energy\s*:\s*([-+\d.EeDd]+)\s+Kcal/mole',txt(v['output']),re.M)
    assert len(tokens)==1
    e=dec(tokens[0]); assert float(e)==v['value']['energy']
    return e

def full_gradient(v):
    text=txt(v['output']); full=[z.split()[2:] for z in text.splitlines() if z.startswith('DONOR_RESPONSE_TOTAL_GRADIENT')]
    if full:
        assert len(full)==8
        g=np.array([[float(dec(q)) for q in z] for z in full]); assert np.array_equal(g,np.array(v['response']['gradient']))
        precision='native full central gradient (17 decimal digits)'
    else:
        lines=[z.split() for z in text.splitlines() if z.split() and z.split()[0]=='Anlyt']
        rows=[z for z in lines if len(z)==6 and z[1].isdigit()]
        assert [int(z[1]) for z in rows]==list(range(1,9))
        g=np.array([[float(dec(q)) for q in z[2:5]] for z in rows]);assert np.array_equal(g,np.array(v['value']['gradient']))
        precision='native baseline eight-decimal gradient'
    return g,precision

def moment(v):
    text=txt(v['output']); perm=[z.split()[1:] for z in text.splitlines() if z.startswith('DONOR_RESPONSE_PERMANENT_LAB_EA')]
    assert len(perm)==2
    p=[dec(q) for q in perm[-1]]
    dual=[z.split() for z in text.splitlines() if z.startswith('LOCAL_RESPONSE_INDUCED')]
    if dual:
        dual=dual[-8:];assert [int(z[1]) for z in dual]==list(range(1,9))
        uv=[[dec(q) for q in z[2:]] for z in dual]
    else:
        lines=[z.split() for z in txt(v['exports'][1]).splitlines() if z.startswith('VECTOR')]
        assert len(lines)==24
        uv=[]
        for i in range(8):uv.append([dec(lines[3*i+j][4]) for j in range(3)]+[dec(lines[3*i+j][5]) for j in range(3)])
    mu=[(p[j]+sum((z[j]+z[j+3])/2 for z in uv))/dec(B) for j in range(3)]
    mv=np.array([float(q) for q in mu]); assert np.max(abs(mv-v['mu_lab_au']))<=ar(24,mv)
    return dot(t,mv),[str(q) for q in mu]

def chi(v):
    lines=[z.split() for z in txt(v['exports'][1]).splitlines() if z.startswith('A ')]
    assert len(lines)==24 and [int(z[1]) for z in lines]==list(range(1,25))
    A=np.array([[float(dec(q)) for q in z[2:]] for z in lines]); z=np.tile(t,8)
    value=Ha/(K*B*B)*dot(z,np.linalg.solve(A,z))
    return value

central={}; targets={}; finite={}; rotation={}
for metal in ['La','Dy']:
    fixture=m['rows'][metal+'_central']; frame=fixture['rotation_frame'];xyz=np.array(fixture['positions_A']);inds=frame['ligand_indices']
    assert frame==rp['central'][metal]['analytic']['frame']
    radial=np.zeros((8,3));radial[inds]=frame['radial_unit']
    rt=np.zeros((8,3));rt[inds]=np.cross(frame['rotation_unit'],xyz[inds]-frame['centroid_A'])/frame['rotation_length_A']
    target=dict(radial_dE_dq_kcal_mol_A=rp['central'][metal]['analytic']['radial'],rotation_dE_dq_kcal_mol_A=rp['central'][metal]['analytic']['rotation'],projected_mu_ea0=dot(t,inv['rows'][metal+'_central']['total_dipole_lab_au']))
    target['chi_a0_cubed']=None
    if metal=='La':
        target['chi_a0_cubed']={step:{key:la5['response']['per_metal']['La']['estimates'][step][key] for key in ['chi_energy','chi_moment']} for step in ['h','half']}
    targets[metal]=target;central[metal]={}
    for profile in ['baseline','drive','stiffness','combined']:
        mode='absent' if profile=='baseline' else profile
        ac=f'{mode}__{metal}_central__analyze' if profile=='baseline' else f'{mode}__{metal}_central__zero__analyze'
        gc=f'absent__{metal}_central__testgrad' if profile=='baseline' else f'{mode}__{metal}_central__fd100__testgrad'
        a=record(ac);g=record(gc); grad,prec=full_gradient(g);e=energy(a);mu,mutokens=moment(a);ch=chi(a)
        assert a['job']['source_input']==fixture['input'] and g['job']['source_input']==fixture['input']
        if profile!='baseline':
            expected=r['operator']['states'][metal+'_'+profile]['matrix_chi_a0_cubed']; assert abs(ch-expected)<=ar(24,ch,expected)
        values=dict(energy_kcal_mol=float(e),radial_dE_dq_kcal_mol_A=float(np.sum(grad*radial)),rotation_dE_dq_kcal_mol_A=float(np.sum(grad*rt)),projected_mu_ea0=mu,chi_a0_cubed=ch)
        errors={k:values[k]-target[k] for k in ['radial_dE_dq_kcal_mol_A','rotation_dE_dq_kcal_mol_A','projected_mu_ea0']}
        if metal=='La':errors['chi_a0_cubed']={step:{key:ch-v for key,v in est.items()} for step,est in target['chi_a0_cubed'].items()}
        central[metal][profile]=dict(values=values,model_minus_reference=errors,energy_decimal=str(e),mu_lab_decimal_au=mutokens,gradient_precision=prec,analyze_row=ac,gradient_row=gc,field_estimates=None if profile=='baseline' else r['field']['states'][metal+'_'+profile]['estimates'])
    q=central[metal];base=q['baseline']['values']
    finite[metal]={'changes':{p:{k:q[p]['values'][k]-base[k] for k in observables} for p in ['drive','stiffness','combined']},'separate_secants_per_unit_amplitude':{p:{k:(q[p]['values'][k]-base[k])/.1 for k in observables} for p in ['drive','stiffness']},'mixed_finite_contrast':{k:q['combined']['values'][k]-q['drive']['values'][k]-q['stiffness']['values'][k]+base[k] for k in observables}}
    rotation[metal]={}
    for step in [.01,.005]:
        stepkey=str(step);tag='010' if step==.01 else '005'; endpoints={}
        for profile in ['baseline','combined']:
            es={};mus={}
            for sign in ['minus','plus']:
                cid=f'pbe0__{metal}__rotation__{sign}_{tag}__complex';mode='absent' if profile=='baseline' else 'combined'
                v=record(f'{mode}__{cid}__analyze');es[sign]=energy(v)
                if profile=='combined':mus[sign]=moment(v)[0]
            slope=float((es['plus']-es['minus'])/(2*dec(step)))
            endpoints[profile]=dict(energy_slope_kcal_mol_A=slope,model_minus_reference_energy_slope=slope-rr['scoring']['directions']['rotation'][metal]['derivatives_kcal_mol_A'][stepkey]['PBE0_D4'])
            if profile=='combined':
                ds=(mus['plus']-mus['minus'])/(2*step);qds=inv['saved_quantum_dipole_geometry_derivatives'][metal][stepkey]['tangent_projection_ea0_per_A']
                endpoints[profile].update(dmu_dq_ea0_per_A=ds,model_minus_reference_dmu_dq=ds-qds)
            else:endpoints[profile]['dmu_dq_inherited_both_field_steps']=model40['summary']['rotation_dipole_derivatives'][metal][stepkey]['field_steps']
        rotation[metal][stepkey]=dict(reference_energy_slope_kcal_mol_A=rr['scoring']['directions']['rotation'][metal]['derivatives_kcal_mol_A'][stepkey]['PBE0_D4'],reference_dmu_dq_ea0_per_A=inv['saved_quantum_dipole_geometry_derivatives'][metal][stepkey]['tangent_projection_ea0_per_A'],profiles=endpoints)
paired={p:{k:central['Dy'][p]['values'][k]-central['La'][p]['values'][k] for k in observables} for p in ['baseline','drive','stiffness','combined']}
paired_finite={p:{k:finite['Dy'][p][k]-finite['La'][p][k] for k in observables} for p in ['mixed_finite_contrast']}
result=dict(status='saved_fixed_donor_sensitivity_complete',implementation=pin(__file__),declaration=pin(HERE/'DECLARATION.json'),input_result=decl['qualified121_result'],input_review=decl['qualified121_independent_review'],targets=targets,central=central,finite_changes=finite,paired_Dy_minus_La_model_values=paired,paired_mixed_contrast=paired_finite,combined_rotation=rotation,used_native_rows=used_rows,verified_sources=[dict(path=k,sha256=v) for k,v in PINS.items()],physical_accuracy_qualified=False,parameters_selected=None,fit_performed=False,rank_or_local_identifiability_claim=False,Dy_and_paired_chi_reference_qualified=False,new_native_calls=0,new_parameter_evaluations=0,limits=decl['limits'],cost=dict(wall_seconds=time.monotonic()-T0,CPU_seconds=time.process_time()-C0,threads=1))
OUT.mkdir(parents=True,exist_ok=False)
(OUT/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(result=pin(OUT/'RESULT.json'),central={metal:{p:x['values'] for p,x in q.items()} for metal,q in central.items()},finite=finite,cost=result['cost'])))
