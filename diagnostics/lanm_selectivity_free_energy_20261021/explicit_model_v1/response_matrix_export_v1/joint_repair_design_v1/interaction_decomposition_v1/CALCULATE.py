#!/usr/bin/env python3
"""Direct saved C/P/B arithmetic; no engine, fitter, or validation reader."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
PARENTS=tuple(f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy'))
HA=627.509474
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read(p):require(pin(p['path'])==p,'Changed source '+p['path']);return json.loads(Path(p['path']).read_text())
def save(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def subtract(a,b):return {k:a[k]-b[k] for k in ('C','P','B')}
def termsum(c,e):return math.fsum(n*e[k] for k,n in c.items())
def assess(model,reference):
    e=subtract(model,reference)
    return dict(model_minus_reference_kcal_mol=e,closure_residual=e['C']-math.fsum((e['P'],e['B'])),
                opposite_signed_pair_collective_errors=e['P']*e['B']<0,
                sum_absolute_component_errors_kcal_mol=abs(e['P'])+abs(e['B']),
                absolute_total_error_kcal_mol=abs(e['C']))

def calculate(out):
    tick,cpu=time.monotonic(),time.process_time()
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'One numerical thread')
    src=json.loads((HERE/'INPUTS.json').read_text());r={k:read(p) for k,p in src.items() if not k.endswith('_declaration')}
    for k in ('common_basis_declaration','matched_pair_declaration'):require(pin(src[k]['path'])==src[k],'Declaration pin')
    selected=r['selected'];em=r['evaluator_manifest'];target=r['targets'];common=r['common_basis'];pair=r['matched_pair'];panel=r['common_basis_panel'];root=r['selected_root_review'];reg=r['regression']
    require(reg['status']=='consumed_trimer_regression_complete' and reg['selected']==src['selected'] and r['regression_review']['result']==src['regression'] and r['regression_review']['status']=='passed_consumed_trimer_arithmetic_review','Admitted consumed-regression candidate join')
    require(root['status']=='admitted_selected_candidate_for_consumed_regression_only' and root['selected']==src['selected'] and root['all_selected_energy_moment_subsets_exact'] is True,'Exact root-admitted selected subset')
    raw_review=read(root['selected_raw_review']);require(raw_review['selected']==src['selected'] and len(raw_review['rows'])==46,'Selected46 raw authority')
    require(selected['status']=='complete' and len(selected['cells'])==46 and all(v['status']=='complete' for v in selected['cells'].values()),'Full selected46 denominator')
    require(r['evaluator_configuration']['manifest']==src['evaluator_manifest'] and em['sources']['target_inventory']==src['targets'],'Static/evaluator/target joins')
    require(common['status']=='common_basis_collection_complete' and common['new_complete']==68 and common['reused_complete']==12 and common['panel']==src['common_basis_panel'] and common['pair_collection']==src['matched_pair'],'Full common-basis collection')
    require(r['common_basis_review']['status']=='passed_independent_actual_common_basis_review' and r['common_basis_review']['collector_result']==src['common_basis'],'Admitted common-basis raw authority')
    require(pair['status']=='matched_pair_collection_complete' and pair['reference_complete']==36 and pair['model_complete']==34,'Original D4 pair denominator')
    static={k:v['energy_kcal_mol'] for k,v in em['static'].items()}
    energies={'baseline':dict(static,**{k:v['energy'] for k,v in em['baseline_cells'].items()}),
              'selected':dict(static,**{k:v['energy'] for k,v in selected['cells'].items()})}
    group=target['groups']['pair_and_collective_dipoles'];require(set(group['parents'])==set(PARENTS) and len(group['pair_cycles'])==34,'Original34 literal pairs/four parents')
    parents={};maxclosure=0.;d4parts={}
    for g in PARENTS:
        n=int(g[4]);metal=g.rsplit('_',1)[1];pairs=[x for x in group['pair_cycles'] if x['parent']==g]
        require(len(pairs)==n and {x['water_index'] for x in pairs}==set(range(1,n+1)),'Literal parent pair denominator')
        xyz=em['rows'][g]['positions_A'];require(xyz==panel['geometries'][g]['positions_A'],'Exact full-parent model/reference geometry')
        C={g:1,'ion_'+metal:-1,g+'__shell':-1};P={}
        for p in pairs:
            for k,c in p['model_terms'].items():P[k]=P.get(k,0)+c
            dimer=g+f"__dimer_{p['water_index']:02d}"
            rc=next(c for c in panel['cells'] if c['parent_geometry_id']==g and c['water_index']==p['water_index'] and c['role']=='metal_water')
            require(em['rows'][dimer]['positions_A']==[xyz[i] for i in rc['real_parent_atom_indices']],'Exact extracted dimer coordinates')
        require(P['ion_'+metal]==-n,'N-fold ion subtraction')
        B=dict(C)
        for k,c in P.items():B[k]=B.get(k,0)-c
        values={}
        for label,e in energies.items():
            require(all(k in e and math.isfinite(e[k]) for k in set(C)|set(P)),'Complete finite model terms')
            cc,pp=termsum(C,e),termsum(P,e);bb=cc-pp;expanded=termsum(B,e);residual=bb-expanded
            require(abs(residual)<=1e-9,'Model literal-cycle arithmetic closure');maxclosure=max(maxclosure,abs(residual));values[label]=dict(C=cc,P=pp,B=bb)
        electronic=common['decomposition']['per_parent_electronic'][g]
        ref=dict(C=electronic['cluster_interaction_hartree']*HA,P=electronic['pair_sum_hartree']*HA,B=electronic['cooperative_remainder_hartree']*HA)
        parents[g]=dict(water_count=n,model_terms=dict(C=C,P=P,B=B),baseline_kcal_mol=values['baseline'],selected_kcal_mol=values['selected'],selected_minus_baseline_kcal_mol=subtract(values['selected'],values['baseline']),
                        reference_electronic_only_kcal_mol=ref,reference_electronic_plus_D4_kcal_mol=None,total_energy_errors=None,
                        unavailable_reason='No saved real-shell D4 authority; electronic-only reference is not a total-energy error target.')
        complex_d4=common['reused_terms'][g]['D4_energy_hartree'];ion_d4=pair['reference_rows']['pbe0__'+metal+'__bare_ion']['D4_energy_hartree']
        dimer_d4=[pair['reference_rows'][f'pbe0__aquo{n}__{metal}__water_{i:02d}']['D4_energy_hartree'] for i in range(1,n+1)]
        require(all(math.isfinite(x) for x in (complex_d4,ion_d4,*dimer_d4)),'Finite saved original D4 terms')
        d4parts[g]=dict(complex=complex_d4,bare_ion=ion_d4,dimers=dimer_d4)
    paired={}
    for n in (8,9):
        la,dy=(f'aquo{n}_central_{m}' for m in ('La','Dy'))
        require(panel['geometries'][la]['positions_A']==panel['geometries'][dy]['positions_A'] and panel['geometries'][la]['symbols'][1:]==panel['geometries'][dy]['symbols'][1:],'Own-shell/water dispersion cancellation geometry')
        electronic=subtract(parents[dy]['reference_electronic_only_kcal_mol'],parents[la]['reference_electronic_only_kcal_mol'])
        # Only paired D4 is evaluated: no unknown shell/water value is set to zero.
        C4=math.fsum((d4parts[dy]['complex'],-d4parts[la]['complex'],-d4parts[dy]['bare_ion'],d4parts[la]['bare_ion']))*HA
        P4=math.fsum([*d4parts[dy]['dimers'],*[-x for x in d4parts[la]['dimers']],-n*d4parts[dy]['bare_ion'],n*d4parts[la]['bare_ion']])*HA
        d4=dict(C=C4,P=P4,B=C4-P4);ref={k:electronic[k]+d4[k] for k in electronic}
        values={label:subtract(parents[dy][label+'_kcal_mol'],parents[la][label+'_kcal_mol']) for label in energies}
        paired['W'+str(n)]=dict(reference_electronic_kcal_mol=electronic,reference_D4_kcal_mol=d4,reference_total_kcal_mol=ref,
                               **{label+'_kcal_mol':v for label,v in values.items()},errors={label:assess(v,ref) for label,v in values.items()})
    coordination={k:subtract(paired['W9'][k],paired['W8'][k]) for k in ('reference_electronic_kcal_mol','reference_D4_kcal_mol','reference_total_kcal_mol','baseline_kcal_mol','selected_kcal_mol')}
    coordination['errors']={label:assess(coordination[label+'_kcal_mol'],coordination['reference_total_kcal_mol']) for label in energies}
    mapping=dict(C='cluster_contrast',P='pair_contrast',B='cooperative_remainder')
    reference_replay={k:coordination['reference_total_kcal_mol'][k]-common['decomposition']['electronic_plus_original_D4_kcal_mol'][v] for k,v in mapping.items()}
    baseline_replay={k:coordination['baseline_kcal_mol'][k]-common['decomposition']['unchanged_model_kcal_mol'][v] for k,v in mapping.items()}
    require(max(abs(x) for x in reference_replay.values())<=1e-9*HA,'Existing common-basis coordination reference replay')
    require(max(abs(x) for x in baseline_replay.values())<=80*.5e-8+1e-10,'Existing printed baseline coordination replay')
    result=dict(status='complete_saved_pair_collective_diagnostic',implementation=pin(__file__),declaration=pin(HERE/'DECLARATION.md'),inputs=pin(HERE/'INPUTS.json'),sources=src,
                definition='C=E(M+shell)-E(M)-E(shell); P=sum(E(M+Wi)-E(M)-E(Wi)); B=C-P',hartree_to_kcal_mol=HA,
                parents=parents,paired_Dy_minus_La=paired,coordination_W9_minus_W8=coordination,source_D4_terms_hartree=d4parts,
                existing_reference_replay_residual_kcal_mol=reference_replay,existing_baseline_replay_residual_kcal_mol=baseline_replay,maximum_model_cycle_closure_kcal_mol=maxclosure,
                independent_reference_SCF_statistics=r['common_basis_review']['scf_statistics'],original_EF_checks_passed=reg['original_EF_checks_passed'],
                per_parent_total_reference_available=False,physical_threshold_added=False,physical_accuracy_qualified=False,new_native_calls=0,new_parameter_evaluations=0,fresh_outcomes_read=False,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1,allocated_scheduler_CPU_seconds=0))
    out=Path(out);out.mkdir(parents=True,exist_ok=False);save(out/'RESULT.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();r=calculate(a.out)
    print(json.dumps(dict(status=r['status'],paired=r['paired_Dy_minus_La'],coordination=r['coordination_W9_minus_W8'],cost=r['cost']),indent=2))
