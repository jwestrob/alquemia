# Whole-LBT fixed-deformation model completion

17 October 2026 UTC. Root declares this finite adjunct to the active coupled
mechanics plan before execution. The independent whole-peptide DFT scout13051
continues unchanged. Its admission/cost gate still controls further DFT.

## Question and scope

What deformation work does the existing whole-MACE plus fixed-charge ALPB model
assign to the exact before/after geometry pair already frozen for the independent
whole-peptide reference? Complete the missing cheap model components so that an
admitted reference can be compared with both mechanical descriptions promptly.
This is a model comparison; neither method's accuracy follows from this test.

Use only the four cells in `../lbt_whole_reference/INVENTORY.json`, SHA
`48d9f40f1bc7f5f927131d5f8e1b2ad3d1fa12187f7b8c321b7cb5375efe4baa`.
Both after cells use the same actual13017 La terminal coordinates, regardless of
the energy sign. Preserve all256 atoms and IDs, target index255, physical charge
-2, free NH3+/COO- termini, five deprotonated acidic groups and physical La/Dy
multiplicities1/6. The Dy MACE spin input is not the effective singlet used by the
separate frozen-f DFT calculation. The admitted80-step endpoint is nonconverged.

| Cell | Whole-MACE | Fixed-FF ALPB | Baseline hybrid |
| --- | --- | --- | --- |
| before La | Reuse energy and forces | Reuse energy and gradient | Reuse existing result |
| before Dy | Reuse energy; forces unavailable | Reuse energy and gradient | Reuse existing result |
| after La | One new energy/force request | Reuse energy and gradient | Reuse existing result |
| after Dy | One new energy/force request | One new energy/gradient request | Unavailable |

The source-only cache audit found no compatible after whole-MACE or after-Dy
solvent receipt in the current mechanics/coupled-response workspaces. The
before-Dy whole-MACE result is energy-only; its missing forces remain explicit.
The existing172-atom hybrid MACE result cannot substitute for whole256 MACE.
Do not replace this after geometry with the separately optimized Dy endpoint.

The hard ceiling is **two whole-MACE requests and one fixed-charge ALPB request**.
There are zero DFT, GFN, MM, optimization, source preparation or new-geometry
requests. The closed failed hybrid-response branch remains closed. Its missing
after-Dy fixed-geometry result is reported unavailable, not silently filled.

## Immutable model and execution

Reuse the actual before-point native MACE helper, checkpoint, model head,
float64 precision, multiplicity encoding and qualified software versions. Bind
their hashes in the execution config. Use one warm GPU worker on an actual
Slurm GPU allocation. Reuse the original native whole-LBT fixed-FF ALPB binary,
its linked libraries and exact256 charges with formal+3 target. Preserve native
La/Dy identity and fixed-charge solvent recipe; no SCC or parameter fallback.

All numerical work, including saved replay and preflight, runs in Slurm. Use
ordinary shared allocation rules and actual granted resources, no memory/time
or scheduler-priority overrides. Root owns submissions and completion wakes.
Outputs are versioned and fresh; preserve partial results/failures and actual
attempt counts. No retry or continuation is included in this declaration.

## Gates and collection

Preflight makes zero molecular calls. It must verify source/hash joins, exact
physical order/coordinates, whole charge and metal spin inputs, actual13049
four-state preparation and frozen directions, two saved whole-MACE energy/state
receipts, the available before-La force array, and all three existing ALPB raw
receipts/gradients. Do not require a nonexistent before-Dy MACE force array or
manufacture it. Reuse the established raw-native admission and unit conventions.

Run requires the actual passed preflight and unchanged scientific signatures.
Only after-La/after-Dy MACE and after-Dy ALPB are executable. Admit every new
native result against its exact physical input, model, charge/multiplicity,
finite energy and full256-atom force/gradient. Verify native ALPB scalar and
gradient provenance by the same parser as the archived three cells. Preserve
failed or unavailable slots without computing a partial-survivor differential.

For each metal report native-unit differences before conversion:
`W_M = E_M(after) - E_M(before)`, then `delta_W = W_Dy - W_La`.
Show MACE vacuum, fixed-FF ALPB and their sum separately. Show the existing
La hybrid deformation work for context and label the Dy hybrid work missing.
No absolute-energy subtraction between different Hamiltonians is meaningful.

Retain full new forces, all available force tables and projections along the
two directions frozen by13049. Whole-MACE energy coverage can reach4/4; force
coverage is at most3/4 under this finite scope. Keep that distinction in outputs.
Use the existing1kcal/mol response reporting resolution; it is not an accuracy
gate, fitted threshold or statistical error bar. Record calls, native execution
time, allocation CPU/GPU time and exact input/output provenance.

## Interpretation and next decision

A different sign or magnitude between the two model descriptions localizes
Hamiltonian dependence. Agreement does not qualify either model. The independent
DFT scout and any subsequently admitted finite reference are required to assess
physical support. No sign-based source selection, force-constant fit, stiffness
rescaling, preferred-metal classification, production change or new structural
search follows from these three model calls. Original state, vacuum-solvent and
terminal-chemistry limitations remain in force.
