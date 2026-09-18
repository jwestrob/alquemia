"""Compare all frozen GGR structural controls without fitting a decision band."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_density_gk_hybrid import TRIAL_PROTOCOL,TOL
from mace_density_panel import PROTOCOL


def compare(panel,output):
    r=read_json(panel);m=read_json(verify(r['manifest']))
    if r['protocol']!=PROTOCOL or not r['complete']:raise InvalidArtifact('complete executed expansion required')
    parent=read_json(verify(m['sources']['reference']))
    if parent['protocol']!=TRIAL_PROTOCOL:raise InvalidArtifact('frozen scientific reference differs')
    new=r['variants']['primary'];old=parent['variants']['primary'];cases={**old['cases'],**new['cases']}
    order=['GGR_extended','GGR_2FW0','GGR_2FVY']
    if set(new['cases'])!={'GGR_2FW0','GGR_2FVY'}:raise InvalidArtifact('declared transfer inventory differs')
    contrasts=[]
    for a in ('ALPHA_1F6S','ALPHA_6IP9'):
        for g in order:
            parts={k:cases[a]['components_R_kcal'][k]-cases[g]['components_R_kcal'][k] for k in cases[a]['components_R_kcal']}
            value=cases[a]['R_kcal']-cases[g]['R_kcal']
            if abs(math.fsum(parts.values())-value)>TOL['algebra_kcal']:raise InvalidArtifact('contrast component algebra differs')
            contrasts.append(dict(alpha=a,GGR=g,difference_kcal=value,components_kcal=parts,pass_=value>TOL['ordering_kcal'],
                origin='previous_development' if g=='GGR_extended' else 'new_structure_transfer'))
    values=[cases[g]['R_kcal'] for g in order]
    result=dict(protocol=PROTOCOL,scientific_protocol=TRIAL_PROTOCOL,panel=record(panel),reference=m['sources']['reference'],
        scientific_plan=m['plan'],preparation=read_json(verify(read_json(verify(m['sources']['boundary']))['manifest']))['preparation'],
        new_numerical_pass=r['numerical_pass'],new_radius_sensitivity_pass=r['radius_sensitivity_pass'],
        ordering_pass=new['ordering_pass'],new_contrasts=new['contrasts'],all_six_contrasts=contrasts,
        all_six_ordering_pass=all(c['pass_'] for c in contrasts),GGR_order=order,GGR_R_kcal=values,GGR_range_kcal=max(values)-min(values),
        prior_partition_kcal=old['partition_kcal'],prior_partition_pass=old['partition_pass'],new_partition_test=None,
        biological_groups=2,interpretation='consumed structural robustness; no independent biological or blind validation',
        evidence_strata={'GGR':'direct same-assay La/Ca affinity direction','ALPHA':'qualified competition/cross-study affinity direction'},
        calibrated_class=None,aqueous_affinity_score=None,baseline_changed=False,cases=cases)
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);write_new(root/'result.json',result)
    lines=['# Frozen hybrid: GGR structure transfer','',
        f'New ordering gate: **{"PASS" if result["ordering_pass"] else "FAIL"}**. Numerical checks: **{r["numerical_pass"]}**; radius checks: **{r["radius_sensitivity_pass"]}**.',
        '', 'Same frozen model; two consumed biological groups. No calibrated absolute score or production change.', '',
        '| Comparison | Alpha minus GGR (kcal/mol) | Fixed >0.02 gate |','|---|---:|---|']
    for c in contrasts:lines.append(f'| {c["alpha"]} − {c["GGR"]} | {c["difference_kcal"]:.9f} | {"pass" if c["pass_"] else "fail"} |')
    lines+=['',f'GGR three-structure range: {result["GGR_range_kcal"]:.9f} kcal/mol. All structures retained.',
        f'Prior 1GLG connected-minus-extended partition: {old["partition_kcal"]:.9f} kcal/mol; no new partition test.',
        '', 'Full components, receipts, source states and exclusions are linked by the machine-readable result.']
    (root/'REPORT.md').write_text('\n'.join(lines)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--panel',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=compare(a.panel,a.output);print({k:r[k] for k in ('ordering_pass','all_six_ordering_pass','GGR_range_kcal')})
