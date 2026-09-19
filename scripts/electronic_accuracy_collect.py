#!/usr/bin/env python3
"""Collect completed real-core CC endpoints without borrowing a baseline reference."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
from affordable_common import InvalidArtifact,HA_TO_KCAL,record,verify,read_json,write_new,digest,energy,xyz

NUMBER=r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'


def compare(energies):
    """Balanced intersite exchange contrast; common method-specific offset cancels."""
    rs={case:es['Ca']-es['La'] for case,es in energies.items() if set(es)=={'Ca','La'}}
    ds={}
    for case,r in rs.items():
        if case.startswith('ALPHA_') and 'GGR_1GLG' in rs:
            ds[case]=(r-rs['GGR_1GLG'])*HA_TO_KCAL
    a='ALPHA_1F6S_original';b='ALPHA_1F6S_water_prepared'
    return {'R_hartree':rs,'alpha_minus_GGR_kcal_mol':ds,
        'water_preparation_delta_R_kcal_mol':None if a not in rs or b not in rs else (rs[b]-rs[a])*HA_TO_KCAL,
        'aquo_reference':None,'calibrated_classification':None}


def cc_output(path):
    """Require the actual CC total and solvent scheme, retaining diagnostic lines."""
    text=Path(path).read_text();total=energy(path)
    # Only completed correlated output can satisfy this collector. In particular
    # a native DFT baseline output is never a successful substitute.
    totals=re.findall(r'E\(CCSD\(T(?:1)?\)\)\s+\.*\s*('+NUMBER+')',text)
    if not totals:raise InvalidArtifact('missing correlated CCSD(T) total')
    cc_total=float(totals[-1].replace('D','E'))
    if abs(cc_total-total)>2e-8:raise InvalidArtifact('final/CC total disagreement')
    if not re.search(r'CPCM scheme\s+\.*\s*PTES',text):
        raise InvalidArtifact('actual PTES coupled-cluster scheme not demonstrated')
    if not re.search(r'Perturbative triple excitations\s+\.*\s*ON',text):
        raise InvalidArtifact('triples not enabled in actual correlated calculation')
    if re.search(r'(?:CCSD|MDCI).*NOT CONVERGED|COUPLED CLUSTER.*NOT CONVERGED',text,re.I):
        raise InvalidArtifact('correlation nonconvergence')
    # Preserve actual printed terms, including any not known before this new
    # executable's first run. These are diagnostics, not additive corrections.
    labels=('E(0)','E(CORR)','E(TOT)','E(CCSD)','E(CCSD(T))','E(CCSD(T1))','Triples Correction',
            'Final correlation energy','C-PCM corr. term','T1 diagnostic',
            'Norm of the','CPCM scheme','Frozen core','frozen core','Number of basis functions',
            'Number of Electrons','Number of electrons','Number of correlated electrons',
            'Total Charge','Multiplicity','ECP Def2-ECP','PNO occupation number cut-off',
            'Smallest eigenvalue','Number of eigenvalues below threshold','Iterative triples')
    lines=[s.strip() for s in text.splitlines() if any(label in s for label in labels)]
    components={}
    for line in lines:
        match=re.search(r'^(.+?)\s+(?:\.\.\.+|:)\s*('+NUMBER+r')\s*$',line)
        if match:components[match[1].strip()]=float(match[2].replace('D','E'))
    return {'total_hartree':total,'printed_components':components,'diagnostic_lines':lines,
        'CPCM_accounting':'included_in_total_not_added_again','output':record(path)}


def state_accounting(text,task):
    atoms=xyz(verify(task['xyz']));z={'H':1,'C':6,'N':7,'O':8,'Ca':20,'La':57}
    if any(a[0] not in z for a in atoms):raise InvalidArtifact('unsupported source element for core audit')
    metal=task['metal'];ecp=46 if metal=='La' else 0
    expected_ne=sum(z[a[0]] for a in atoms)-ecp-task['charge']
    expected_fc=2*sum(a[0] in ('C','N','O') for a in atoms)+(10 if metal=='Ca' else 0)
    ne=re.findall(r'Number of electrons\s+\.*\s*(\d+)',text)
    fc=re.findall(r'Frozen core treatment\s+\.*\s*chemical core \((\d+) el\)',text)
    nc=re.findall(r'Number of correlated electrons\s+\.*\s*(\d+)',text)
    if not ne or not fc or not nc:raise InvalidArtifact('actual electron/core accounting missing')
    if (int(ne[-1]),int(fc[-1]),int(nc[-1]))!=(expected_ne,expected_fc,expected_ne-expected_fc):
        raise InvalidArtifact('actual electron/core accounting differs from intended policy')
    if metal=='La' and 'ECP Def2-ECP (replacing 46 core electrons' not in text:
        raise InvalidArtifact('native La46 ECP not demonstrated')
    return {'explicit_electrons':expected_ne,'frozen_explicit_electrons':expected_fc,
        'correlated_electrons':expected_ne-expected_fc,'ECP_electrons':ecp}


def collect(manifest):
    from run_orca_task_manifest import load_manifest_tasks,_completed_attempt_is_valid
    mp=Path(manifest).resolve();m,tasks=load_manifest_tasks(mp)
    prep=read_json(verify(m['preparation']));rows=[];es={};bs={}
    source={(s['case'],s['metal']):s for s in prep['sources']}
    metadata={t['task_id']:t for t in m['tasks']}
    for t in tasks:
        original=metadata[t['task_id']]
        case=original['case'];metal=original['metal'];s=source[(case,metal)]
        bs.setdefault(case,{})[metal]=energy(verify(s['baseline']['output']))
        op=t['output'];rp=Path(str(op)+'.execution.json')
        r={'case':case,'metal':metal,'task_id':t['task_id'],'status':'not_run','result':None,'receipt':None}
        if op.exists():
            r['status']='incomplete';r['output_path']=str(op)
            if rp.exists():
                r['receipt']=record(rp);r['output']=record(op)
                try:
                    if not _completed_attempt_is_valid(rp,op,manifest_sha256=digest(mp),task=t,
                        runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
                        raise InvalidArtifact('completed execution receipt invalid')
                    r['result']=cc_output(op)
                    r['state_accounting']=state_accounting(op.read_text(),original)
                    r['status']='complete';r['receipt']=record(rp)
                    es.setdefault(case,{})[metal]=r['result']['total_hartree']
                except (ValueError,KeyError,OSError) as exc:
                    r['status']='failed_or_unsupported';r['reason']=str(exc)
                    r['result']=None
        rows.append(r)
    ep=mp.parent/'budget_events.jsonl'
    events=[json.loads(s) for s in ep.read_text().splitlines()] if ep.exists() else []
    return {'protocol_id':m['protocol_id'],'manifest':record(mp),'rows':rows,
        'complete_endpoints':sum(r['status']=='complete' for r in rows),'total_endpoints':len(rows),
        'baseline':compare(bs),'correlated':compare(es),'baseline_changed':False,
        'cost_events':[e for e in events if 'allocated_core_seconds' in e],
        'interpretation':'combined electronic-recipe sensitivity; no unique error attribution or new affinity bands'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();result=collect(a.manifest);write_new(a.output,result)
    print(json.dumps({k:result[k] for k in ('complete_endpoints','total_endpoints','baseline','correlated')},indent=2))
