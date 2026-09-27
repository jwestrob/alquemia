# Response transfer preparation: 4MAE ready; 1F6S representation choice explicit

**4MAE's six-cell inputs are ready from existing matched artifacts, without
building a new force-field system.** The consumed non-PQQ identity is frozen as
alpha-lactalbumin 1F6S chain-A strong site. A compatible historical 52-atom
representation exists for 1F6S, but it is not the initially inventoried 40-atom
repair and has not been silently substituted. No molecular calculation or job
was launched by this preparation task; root owns subsequent execution.

## 4MAE: exact canonical dry core and matched parent

Manifest: `workspaces/metal_environment_response_20260926/transfer_preparation/4MAE_v1/INPUTS.json`.
It directly supports the existing `metal_environment_reference.py prepare`
interface, with Ca/La endpoint dictionaries and A/B point-charge dictionaries.

The canonical core has 80 atoms: metal, full PQQ, Glu172, Asn256, Asp299, Asp301
and Arg326 side chains. Exact source/carve atom records and original XYZ bytes
are retained. Core charges are Ca −3 / La −2, multiplicity 1. All-electron counts
are 378 and 414; native La ECP46 gives 368 explicit electrons. No biological label
or released decision threshold is used by this preparation.

The already serialized scaffold parent inventory uses the exact same protonated
source, SHA256 `92b0bc51200e5407cb78d5dfd7150e3f319cec60dca7785a336b97340e0c0d9c`.
Its 8,826 source/archived-completion atom coordinates and ff19SB charges are reused;
no new System or force-field template is instantiated. The pre-existing terminal
OXT completion at Glu577 is explicitly pinned. Every original source hydrogen
coordinate is retained; the older whole-protein radial H normalization is not
imported. The conditional parent protein charge is −3 e.

The five side-chain boundaries use the same original policy as the 1H4I scout:
remove source atoms represented in QM and their adjacent CA charge, then restore
each retained residue's prescribed charge by equal increments to its backbone
N and C. Five CB/CA cap maps and physical Jacobians are saved. The resulting
field has **8,774 charges, sum −1 e**, equal for both metals. PQQ and metal have
no MM force-field charges or invented parameters.

Same predeclared proximity rule selects **Thr154 HG1**: its oxygen is 2.52049 Å
from the nearest core heavy atom. B rotates that actual H by +10° about CB→OG1,
moving it 0.19317 Å while preserving the inherited 1.18042 Å O–H length and
C–O–H angle. Only field index 2407 changes coordinates. The closest H-heavy
nonbonded distance after perturbation is 1.90039 Å, above the 1.2 Å screen.
Source hydrogen geometry remains a limitation, not an equilibrated structure.

**Conditional-source limitation:** this exact canonical dry preparation excludes
all source waters and the directly coordinating 15P603 crystallization adduct,
without replacement. Those prior exclusions remain explicit in `INPUTS.json`;
this experiment neither restores them nor claims to model the complete original
experimental coordination state. Source assembly, protonation and exclusions
are fixed within the comparison.

Four real-fixture checks pass (`CHECK_4MAE.json`): paired coordinates/state/electron
parity; only the selected H moves with unchanged charge list; exact field-file
serialization and maps; and preparation/dry-run of all six cells through the
existing runner at 6×57 MPI. No molecular integration test ran here. Numerical
qualification remains unestablished, and no full-hybrid or biological validation
claim follows from these preparation checks.

## 1F6S: frozen biological identity, two distinct archived representations

Selection basis recorded before new energies: a compact repaired amide region,
two explicitly mapped coordinating waters, and complete archived chain-A termini.
This is a consumed strong-site affinity fixture, not a fresh blind label. It was
not selected by response sign or classification. `NONPQQ_INVENTORY.json` pins all
paths and the checks below.

The original repaired **40-atom** core is in:
`workspaces/benchmark_set_20260915/prepared/alacta_1f6s_v2/strong_site/amide_v3/repair_manifest.json`.
Its source includes Asp82/Asp87/Asp88 side chains, whole peptide-amide donors for
Lys79/Asp84, and waters211/212. Its Ca charge is −1 and La charge 0. The general
peptide boundaries require their actual graph/cap treatment; the PQQ side-chain
N/C rule cannot simply be copied onto them.

The whole-protein preparation summary initially inspected stores no individual
force-field charges. A deeper, existing archive **does** contain an ff19SB atomic
charge background and native point-charge inputs:
`workspaces/mace_omol_20260917/full_boundary_GB_v1/states/ALPHA_1F6S/state.json`
and `responsive_quantum_v2/ALPHA_1F6S_{Ca,La}/`.
This corrects the initial incomplete inventory; these are not AMOEBA charges.

That existing native representation has **52 QM atoms and 1,880 external charges
summing to −4 e**, Ca/La core charge −1/0. Both metal fields are byte-identical.
Each point-charge file exactly regenerates from the archived physical-coordinate
and weight map. The source retains normalized protein/water H geometry and the
separate NMA source preparation. Its general graph boundary removes the mapped
projection support and redistributes residue-local charge onto actual bonded
exterior recipients, which can include H atoms. The full ledger is retained.
It is a different approximation from the PQQ side-chain-only boundary policy.

Thus the **blocker for using the initially selected 40-atom input** is the absence
of a currently established matching field/cap transformation for that particular
preparation. It is not absence of all parameterized 1F6S representations. A viable
next preparation route is to explicitly adopt the existing 52-atom normalized
representation after root reviews its physical cap/source map and boundary
policy, preserving its H coordinates and water inventory throughout. That needs
no force-field rebuild and must not be called a replay of the 40-atom repair.

For orientation planning only, the same proximity rule on the 52-atom archive
selects Thr86 hydroxyl (O-to-core distance4.64009 Å). No B configuration, new
endpoint manifest, molecular evaluation or submission has been produced for it.
No method was switched after seeing a new response, and no reserved outcomes
were opened.

## Reproduction and restart

From the repository root, reproduce 4MAE into a new directory:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/metal_environment_response_20260926/transfer_preparation/prepare_4mae.py \
  --repository "$PWD" \
  --output workspaces/metal_environment_response_20260926/transfer_preparation/4MAE_reproduction
```

`check_4mae.py` takes `--repository`, `--inputs` and a new `--output` JSON path.
The current `4MAE_v1` inputs and implementation pins are frozen; do not rewrite
while root's separately owned calculations use them. Root supplied the execution
plan and owns collection. This task has zero new molecular calls, zero allocated
CPU/GPU work and no measured process-memory record. Shared production, runners
and existing results remain unchanged.
