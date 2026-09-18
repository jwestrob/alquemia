# Frozen responsive-density hybrid fails one GGR structural control

**Completed; ordering gate fails.** Both new structures pass all18 numerical,
identity, rotation and radius checks. 2FVY retains the expected alpha-over-GGR
ordering; 2FW0 reverses both comparisons. The earlier development pass remains
real but does not transfer across all three GGR structures. Baseline/default
and every previous scientific record remain unchanged.

| Alpha minus GGR, kcal/mol | 1GLG extended (previous) | 2FW0 (new) | 2FVY (new) |
|---|---:|---:|---:|
| Alpha1F6S | +10.127062434 | −14.764159019 | +8.927577642 |
| Alpha6IP9 | +9.252282954 | −15.638938500 | +8.052798162 |

All four new differences were required to exceed+0.02; two fail. The complete
six-comparison display retains every structure, without averaging or selection.
GGR spans24.891221454kcal across1GLG/2FW0/2FVY. The earlier1GLG connected-minus-
extended partition remains−0.541233142kcal; no new partition test was performed.
Two consumed biological groups remain two groups. No absolute calibration,
aquo-reference decision or prospective validation is claimed.

## What this teaches us

For2FW0 minus2FVY, the full contrast changes+23.691736661kcal:

| Component of the difference | kcal/mol |
|---|---:|
| Intrinsic responsive quantum core | +0.242251473 |
| Direct density–protein coupling | +13.421787726 |
| Permanent GK transfer | +12.600209216 |
| Induced environment transfer | −3.206575910 |
| MACE short full-minus-core | +0.634064156 |

The environmental terms account for+22.815421032kcal. Thus this new discrepancy
is not principally a change in the intrinsic quantum or learned short-context
contrast. It is also far above measured numerical error. This localizes the
problem; it does **not** prove which environmental approximation is responsible,
or exclude a real structure/state effect. The densities respond toff19SB,
then are evaluated with AMOEBA/GK; they remain non-self-consistent trial densities.
GK still uses projected CHELPG sources. Source exclusions and crystal differences
remain confounded, including missing sugar/ions and reported2FVY Glu149 damage.
The sugar-free/open2FW0 failure cannot be discarded because2FVY agrees.

Do not promote this candidate or tune its threshold to these structures.
Retain baseline and pursue a specific accuracy/representation question with
new declared versions. Multisite preparation support proceeds independently;
it is not evidence that this failed candidate will improve on other families.

## Executed science and physical checks

Four native r2SCAN-3c embedded endpoints, eight exact MACE short forwards,
four native CHELPG fits, four49-point/ESP potential queries, two AMOEBA
frameworks/two native parameter reads, eight source/reference initializations,
27static energies/24field queries/30response solves actually ran. No scientific
retry was needed. All fixed chargefit/projection/field gates pass. Source and
reference have identical physical boundaries and external permanent states.

Full model: maxresponse-refinement environment error6.141e−8kcal, maxrigid
component error5.230e−12kcal, maxrelative radius effect0.000343146kcal, both
vacuum identity errors0. All18checks pass. Numerical credibility is supported
for this implementation; predictive robustness fails; cost is practical for
this small pilot but production cost relative to a matched baseline is unmeasured.

Two initial quantum *preflight* failures occurred before any DFT (cross-host
roundoff diagnostic equality and missing snapshot import). Actual endpoints
were later recollected using the existing qualified analytic-gradient parser;
ORCA's conditional numerical-gradient warning had triggered the old parser.
No coordinates, energies or quantum job were changed to repair parsing.
Observationv1 failed frozen import preflight before utilities. The final native
preparation initially rejected identical XYZ bytes at different copied paths;
verified byte equality fixed that before preparation/native execution.

New panel tests: all7distinct tests pass across preflight and final actual-output
runs (6pass102.894s before execution; final integration1pass15.474s). Existing
hybrid6pass33.062s. Earlier input/short/observation15 and legacy14 tests remain
recorded in COMPONENT_STATUS. Legacy trial recollection reproduces every score
and decision exactly; four rigid dipole contractions differ<=3.553e−15 across
hosts. The regression permits1e−12 only for those contractions; remaining fields
must match exactly. No final integration test is fabricated or left unrun.

## Measured execution cost

| Job / phase | Wall s | Allocated core-s | Actual CPU s | GPU s |
|---|---:|---:|---:|---:|
|1201149/1201152, failed preflights |6|384|4.346487|0|
|1201154,4DFT |232|14848|12744.437828|0|
|1201158,8MACE short |112|1792|127.215067|112|
|1201161,4fits+4queries |502|2008|1961.252|0|
|1201162,native controls/collection |409|26176|2812.966|0|
| **Recorded allocation total** | **1261** | **45208** | **17650.217382** | **112** |

Wall is the sum of job durations, not elapsed project time. This includes the
full development checks on two structures, not an ordinary production pair.
Native kernels used100.439151208s; native child processes152.151604574wall/
2566.599330CPU. Primary six native calls cost16.415333wall/274.729563CPU for2FVY,
16.386246wall/274.158CPU for2FW0, excluding runner/preparation/DFT/MACE/utilities.
Native job peakRSS1516432KiB; density utility peakRSS687624KiB. Neither usesGPU.
Six successful local preparation stages additionally recorded207.301142wall/
186.765411parentCPU seconds; native initialization children are separate in
receipts. Unmetered shell/debug work is not represented as zero cost.

## Implementation, provenance and operations

Scientific expression unchanged from5f7d164:
`saved_responsive_trial_density_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1`.
New runner protocol:
`declared_source_graph_responsive_density_GK_POLAR_panel_v1`.
It reuses the native solver, task receipts/recovery and collector equations.
The configurable source/input/short/boundary adapters retain old defaults.
No baseline cache satisfies this challenger. Missing results remain unavailable.

Pre-output scope: [expansion plan](TRIAL_DENSITY_GK_EXPANSION_PLAN.md) and
[query scope](TRIAL_DENSITY_GK_EXPANSION_QUERY_SCOPE.md). Compact actual pins,
unrounded comparisons/components/checks/costs:
[RESULT](TRIAL_DENSITY_GK_EXPANSION_RESULT.json). Full outputs live under
`workspaces/mace_omol_20260917/trial_gk_expansion_*`.

From repository root, inspect the immutable executed manifest without rerunning:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/trial_gk_expansion_native_v1/implementation/mace_density_panel.py dry-run --manifest workspaces/mace_omol_20260917/trial_gk_expansion_native_v1/manifest.json
```

Rebuild a comparison into a fresh output directory:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_density_panel_report.py --panel workspaces/mace_omol_20260917/trial_gk_expansion_native_v1/collection_job_1201162.json --output workspaces/mace_omol_20260917/trial_gk_expansion_comparison_review
```

`mace_density_panel.py prepare --config PATH --output NEW_DIRECTORY`, `dry-run`,
`execute`, and `collect` expose the manifested operations. Existing concrete
config is `trial_gk_expansion_native_config_v1/config.json`; original57tasks
are complete and should not be resubmitted. Execution uses the recorded64CPU
wrapper `run_trial_gk_expansion_native.sbatch`. Explicit `--retry-failed` exists
for technical recovery; it does not authorize changing frozen science.
