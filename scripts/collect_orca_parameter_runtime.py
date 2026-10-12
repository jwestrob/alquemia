#!/usr/bin/env python3
"""Read-only numerical collection of the completed, finite 60 native controls."""
from pathlib import Path
import csv, json, shutil
import orca_parameter_runtime_qualification as b

W=b.W

def main():
    manifests=[W/x for x in ['manifest_La8.json','manifest_aquo50.json','manifest_Ca2.json']]
    jobs=[j for p in manifests for j in b.read(p)['jobs']]
    rows={j['id']:b.read(W/'results'/j['id']/'receipt.json') for j in jobs}
    assert len(rows)==60 and len(list((W/'results').glob('*/receipt.json')))==60
    for p in manifests:
        for j in b.read(p)['jobs']:
            r=rows[j['id']]
            assert r['manifest_sha256']==b.sha(p) and r['inputs_after_match']
            for name,entry in j['files'].items():
                assert b.sha(W/'results'/j['id']/name)==entry['sha256']
            for entry in r['artifacts'].values():
                assert b.sha(entry['path'])==entry['sha256']
    la_ids=[j['id'] for j in b.read(manifests[0])['jobs']]
    ca_ids=[j['id'] for j in b.read(manifests[2])['jobs']]
    source_tables={n:b.checked_source('la_followup/native_reciprocal_controls/'+n) for n in ['SUMMARY.json','endpoint_comparison.tsv','all_atom_charge_comparison.tsv']}
    endpoints={r['job']:r for r in csv.DictReader(source_tables['endpoint_comparison.tsv'].open(),delimiter='\t')}
    charges={}
    for r in csv.DictReader(source_tables['all_atom_charge_comparison.tsv'].open(),delimiter='\t'):
        charges.setdefault(r['job'],[]).append((int(r['atom_index_0based']),float(r['upstream_charge_e'])))
    la_upstream=[]
    for jid in la_ids+ca_ids:
        original=rows[jid]['job'].get('baseline_job_id',jid)
        ordered=sorted(charges[original]);assert [i for i,_ in ordered]==list(range(222))
        ref={'energy_hartree':float(endpoints[original]['upstream_energy_hartree']),'charges_e':[q for _,q in ordered],'accepted':True}
        la_upstream.append({'job':jid,'comparison':b.compare(rows[jid]['output'],ref),
                            'reference_status':'accepted published-parameter upstream values preserved in the checked historical comparison tables; no new upstream calls'})
    interim=b.read(W/'RESULTS.json');pairs=interim['aquo_pairs']
    runtime=b.read(W/'OBSERVED_RUNTIME_FINAL.json')
    for e in runtime['observed_files']:assert b.sha(e['path'])==e['sha256']
    pre={e['path']:e['sha256'] for e in b.read(W/'RUNTIME_BINDING.json')['files']}
    assert all(e['sha256']==pre[e['path']] for e in runtime['observed_files'] if e['path'] in pre)
    # Check only owned historical PIDs; no other sessions are enumerated or signalled.
    live=[]
    for r in rows.values():
        for p in r['observed_processes']:
            cmd=Path('/proc')/str(p['pid'])/'cmdline'
            if cmd.exists():
                try:
                    if str(W).encode() in cmd.read_bytes():live.append(p['pid'])
                except FileNotFoundError:pass
    assert not live
    executions={p.stem:b.read(p) for p in W.glob('EXECUTION_*.json')}
    known_cpu=sum(e.get('reaped_child_user_seconds',0)+e.get('reaped_child_system_seconds',0) for e in executions.values())
    failures=[{'job':jid,'process':r['process'],'failed_output_checks':[k for k,v in r['output']['checks'].items() if not v],
               'failed_property_checks':[k for k,v in r['output']['property_checks'].items() if not v],
               'convergence':r['output']['convergence']} for jid,r in rows.items() if not r['passed']]
    profiles=[]
    for jid,r in rows.items():
        if r['job']['variant'] not in ['corrected','corrected_file','corrected_tight']:continue
        name='control.xtb.json' if 'control.xtb.json' in r['job']['files'] else 'parameters.xtb.json'
        profiles.append({'job':jid,'parameter':r['job']['files'][name],'state':r['job']['expected'],
                         'medium':r['job']['medium'],'native_accepted':r['passed'],'receipt':b.pin(W/'results'/jid/'receipt.json')})
    result={**interim,'schema':'native_runtime_qualification_results.v2','updated_utc':b.now(),'status':'complete','maximum_authorized_native_endpoints':60,
            'La8':{'complete':8,'accepted_and_archived_parity':sum(rows[j]['passed'] and rows[j]['comparison_to_archived']['accepted_parity'] for j in la_ids)},
            'Ca_corrected_full_file2':b.read(W/'EXECUTION_Ca2.json'),'La_Ca_upstream_comparisons':la_upstream,
            'aquo_corrected_summary':{'cells':14,'accepted_against_both':sum(x['corrected_vs_xtb']['accepted_parity'] and x['corrected_vs_tblite']['accepted_parity'] for x in pairs),
                'maximum_energy_residual_Eh':max(abs(x[k]['energy_delta_hartree']) for x in pairs for k in ['corrected_vs_xtb','corrected_vs_tblite']),
                'maximum_all_atom_charge_residual_e':max(x[k]['max_atom_charge_delta_e'] for x in pairs for k in ['corrected_vs_xtb','corrected_vs_tblite']),
                'roundtrip_raw_pass':sum(x['roundtrip_vs_default']['numerical_pass'] for x in pairs),
                'roundtrip_accepted_pass':sum(x['roundtrip_vs_default']['accepted_parity'] for x in pairs)},
            'failures':failures,'manifests':[b.pin(p) for p in manifests],'plans':[b.pin(W/x) for x in ['PLAN.json','PLAN_EXTENSION_Ca2.json']],
            'source_tables':{k:b.pin(v) for k,v in source_tables.items()},'observed_runtime':b.pin(W/'OBSERVED_RUNTIME_FINAL.json'),
            'additional_runtime_files_observed_after_initial_binding':[e for e in runtime['observed_files'] if e['path'] not in pre],
            'corrected_import_evidence':profiles,'owned_processes_remaining':live,'cost':{
                'all60_summed_endpoint_wall_seconds':sum(r['process']['elapsed_seconds'] for r in rows.values()),
                'measured_reaped_child_CPU_seconds_for_aquo50_and_Ca2':known_cpu,
                'La8_CPU_seconds':None,'La8_CPU_status':'not measured by original stage wrapper; actual per-endpoint wall times retained',
                'maximum_process_group_RSS_bytes':max(r['process']['peak_process_group_rss_bytes'] for r in rows.values()),
                'max_observed_threads_per_process':max(r['max_process_threads'] for r in rows.values()),
                'maximum_concurrent_numerical_children':4,'GPU_seconds':0,'Slurm_jobs':0,'new_xtb_or_tblite_calls':0},
            'metadata_caveats':['Ca2 receipt job.parameter_source remains null and allowed_input_changes empty from copied baseline; actual input/parameter byte hashes, stdout import lines and separate declared extension manifest establish the real import. Historical execution receipts are unchanged.',
                'Initial runtime binding covered 10 files. Actual process monitoring added orca_guess and dash; final12-file observed identity is pinned after execution.',
                'Corrected aquo profiles change only the targeted metal four arrays; they are distinct full-file identities and do not jointly qualify mixed-metal combinations.'],
            'scope_limit':'Exact Linux ORCA 6.1.1 implementation controls, supplied geometries and printed-state conventions only. Spin energetics, ligand transfer, affinities, biological ranking and other builds remain unqualified.'}
    for p in [W/'RESULTS.json',b.REPORT/'RESULTS.json']:b.write(p,result,replace=True)
    for name in ['RUNTIME_BINDING.json','OBSERVED_RUNTIME_FINAL.json','manifest_La8.json','manifest_aquo50.json','manifest_Ca2.json','PLAN.json','PLAN_EXTENSION_Ca2.json']:
        shutil.copyfile(W/name,b.REPORT/name)
    b.write(b.REPORT/'APPROVAL_INPUTS.json',{'source':'Actual runtime identity and corrected-import evidence; guard owns policy/approval generation.','runtime':b.pin(W/'RUNTIME_BINDING.json'),'observed_runtime':b.pin(W/'OBSERVED_RUNTIME_FINAL.json'),'imports':profiles},replace=True)
    summary=result['aquo_corrected_summary'];cost=result['cost']
    lines=['# Linux native-GFN2 parameter runtime qualification — 12 October 2026','',
        'All **60 declared native controls completed** on the pinned Linux ORCA 6.1.1 runtime. **52 passed every original state, process, population and printed SCF criterion**. Eight default/roundtrip endpoints retain their density failures. No retry, tolerance change, production edit or independent-solver call was made.','',
        '## Results','',
        '| Control | Actual result |','|---|---|',
        '| Original Ca/La eight-cell replication | 8/8 accepted; energies and all 222 atomic charges reproduce archived Linux controls |',
        '| Corrected full La JSON imported into Ca, both media | 2/2 accepted; energies and all 222 charges exactly equal current defaults |',
        '| Cu1/Cu2/Mn2/Fe2/Fe3/Co2/Ni2 corrected aquo, both media | 14/14 accepted against both archived xTB 6.7.1 and tblite 0.7.0 |',
        '| Unchanged aquo export/import roundtrip | 14/14 raw numerical matches; 10/14 pairs accepted under all density gates |',
        '| Ca2/Mg2/Sr2/Zn2 default aquo | All 8 native endpoints admitted; 7/8 upstream parity passes; Zn water charge parity fails |','',
        f'Corrected aquo maximum residuals are {summary["maximum_energy_residual_Eh"]:.9g} Eh and {summary["maximum_all_atom_charge_residual_e"]:.9g} e, below the frozen 1e-6 Eh/1e-5 e thresholds. Unchanged roundtrip thresholds remain 1e-8 Eh/1e-7 e. Comparisons reuse accepted archived upstream outputs; no Linux xTB/tblite execution is implied.','',
        'Default and unchanged-roundtrip Cu2 vacuum/water, Fe2 water and Fe3 water fail printed density limits despite normal termination. Their exact roundtrip agreement does not admit them. Corrected_tight runs pass the original density criteria. Zinc water retains a maximum all-atom charge residual 8.74284557e-5 e; its small energy residual does not rescue parity.','',
        'The original La correction is exactly three arrays with declared shells and all other fields unchanged. Fresh corrected La reproduces the archived result: the motivating Ca−La solvent correction changes by **+22.600061589 kcal/mol**. This example-specific correction does not restore old calibration or establish protein prediction.','',
        '## Protocol, identity and evidence','',
        'The immutable original plan declared 58 calls. Before execution, PLAN_EXTENSION_Ca2.json declared two unchanged-Ca full-file import controls and the new ceiling 60. Stage1 used the exact existing eight fixtures. The aquo 50 geometry and native input byte hashes match the synced source artifacts. The source inputs establish TolE 1e-10 for default/roundtrip and 1e-12 for corrected_tight, overriding the task prose ambiguity before stage2. Other solver/state/temperature settings and all gates were retained.','',
        'RUNTIME_BINDING.json pins the actual driver, startup/SCF helpers, shared libraries and environment. OBSERVED_RUNTIME_FINAL.json adds actually observed orca_guess and dash, for 12 runtime files. No preregistration of those two additional observed helpers is claimed. All initial pins still match; the separate guard binds the final full identity. APPROVAL_INPUTS.json supplies each real corrected import, exact full-file SHA, medium and printed-state convention. La/Ca use full-file SHA a4f83ca72d58d6e80f634e61bb8e99b7af0f26c38cc4a9295e85215ecf72d05c.','',
        'RESULTS.json preserves raw versus accepted comparisons, all 60 receipt pins, original failures, source tables, finite manifests and costs. The actual full-precision property populations are compared, including charge sums and atom indices. Each imported file is byte checked before/after and stdout records the import. Ca2 receipts retain harmless copied null parameter_source/empty allowed_input_changes metadata; actual hashes, import lines and the separate declared extension establish the real change. No historical receipt was rewritten.','',
        '## Cost and limits','',
        f'Actual summed endpoint wall time: {cost["all60_summed_endpoint_wall_seconds"]:.3f} s. Measured reaped child CPU time for the 50 aquo plus 2 Ca controls: {known_cpu:.3f} s. CPU time for the first 8 was not separately captured; their measured wall times remain in receipts. Peak owned process-group RSS: {cost["maximum_process_group_RSS_bytes"]} bytes. At most 4 concurrent single-thread numerical children; every observed process had 1 thread. Zero Slurm jobs/GPU/new independent-solver calls; no owned processes remain.','',
        'These controls qualify this runtime and explicit parameter/state combinations for further guarded research. Distinct target-metal full JSON profiles do not qualify mixed-metal unions. Spin-independent GFN2 does not establish physical spin-state energetics. Protein affinities, discrimination, accommodation/calibration and untested builds remain separate questions.','']
    (b.REPORT/'REPORT.md').write_text('\n'.join(lines))
    (W/'CHECKPOINT.md').write_text('# Completed native runtime qualification\n\n60/60 declared calls complete;52 accepted,8 density-gate failures preserved.14/14 corrected aquo comparisons pass both archived upstream engines;Znwater charge parity fails. No owned processes remain. Full report: diagnostics/orca_parameter_audit_20261012/runtime/REPORT.md. Source identities, commands, costs, failures and next-step limitations: RESULTS.json and manifests. No rerun authorized by this checkpoint.\n')
    b.status('complete',[],completed=60,accepted=52,next_step='Independent review and guarded dispatch dry-run; no molecular reruns')
    print(json.dumps({'completed':60,'accepted':52,'aquo':summary,'cost':cost}),flush=True)

if __name__=='__main__':main()
