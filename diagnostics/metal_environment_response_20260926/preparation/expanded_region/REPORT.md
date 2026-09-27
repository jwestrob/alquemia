# Thr159 expanded-region preparation — ready, not scored

**The same real A/B configurations are now prepared with Thr159 in the electronic
region. No electronic calculation, optimization or scheduler job was run.**
The preparation does not resolve the ongoing native SCF qualification or the
missing complete-hybrid interactions.

Current input manifest:
`workspaces/metal_environment_response_20260926/preparation/expanded_region_v1/INPUTS.json`.
The endpoint keys are `Ca_A`, `La_A`, `Ca_B`, `La_B`; each records XYZ, field,
charge, multiplicity, electron counts and source map. All files are hash-pinned.

## Physical question and exact scope

Thr159 was selected as the nearest complete exterior hydroxyl before reading the
response results. Its addition asks whether the electronic environmental response
survives moving that real group across the representation boundary. No region was
chosen for a favorable sign. It is a partition diagnostic: the changing HG1 is
now a quantum atom, so this is not a second fixed-core environmental test.

The exact original protonated 1H4I source and archived A/B preparation were reused.
No newer scaffold source, protonation, H normalization, force-field convention or
water inventory was substituted. All 54 old core coordinates and mappings remain
unchanged. Source-graph traversal after cutting the actual CB–CA bond identifies
exactly CB, HB, OG1, HG1, CG2, HG21, HG22 and HG23 as the complete Thr side chain.
These eight source atoms retain their actual coordinates; one 1.09 Å H link along
CB→CA completes the local fragment. New core size is 63 atoms.

All original physical protein atoms remain accounted for in the pinned source
and maps. The CA coordinate remains a physical boundary endpoint, with link
forces projected onto CB and CA. CA's direct point charge is removed according
to the old boundary prescription; it is not a physical deletion or a freely
moving synthetic cap. The other three original cap maps remain unchanged.

A and B differ only at **QM index 57, Thr159 HG1**. This is exactly the earlier
+10° orientation change, including its inherited 1.18537 Å O–H length. The expanded
external point-charge files are byte-identical between A and B. No full-system
minimum, equilibrium state or hydrated-ion reference is claimed.

## Charge/state closure and original boundary policy

The exact original `affordable_state.py::prepare_skeleton` policy removes charges
on every promoted source atom and adjacent CA, then divides the residue-local
charge correction equally between retained backbone N and C. This is the same
policy used for Glu177, Asn261 and Asp303, not a newly fitted convention.

Thr159's full ff19SB charge is zero within rounding. The retained MM part initially
sums to −0.0137 e; N and C each receive **+0.00685 e**. The promoted Thr fragment
has formal charge zero, so the final retained residue charge is zero. The field
contains 9,078 charges and remains −6 e. The ledger records every removed source
index, recipient and increment.

| Endpoint | QM charge | Multiplicity | All-electron count | r2SCAN-3c explicit electrons |
|---|---:|---:|---:|---:|
| Ca A/B | −3 | 1 | 316 | 316 |
| La A/B | −2 | 1 | 352 | 306, with 46-electron La ECP |

Counts follow actual elements and the unchanged original formal charges. Native
runtime output must still verify these settings before accepting any result.
The total conditional-system charges remain Ca −9 and La −8.

## What this does not establish

Charge closure and exact coordinate transfer do not guarantee partition
consistency. Eight force-field monopoles are replaced by responsive electrons
and nuclei, a synthetic cap is added, and the CA charge is redistributed. Those
are actual changes to the approximate Hamiltonian. The new calculations can
quantify their combined effect but cannot identify its unique cause from a score
jump alone. The original force-field H geometry remains imperfect.

The energy ledger still lacks defensible complete metal/PQQ–exterior cross
repulsion/dispersion and a complete covalent-boundary mechanical definition.
These inputs support an explicitly limited electronic-component test; they do
not authorize a complete QM/MM affinity, relaxation correction or classifier.
No thresholds, labels, production paths or library outcomes were modified.

## Checks and resumption

Preparation assertions verify the actual bonded component, source parameter
charges, original core retention, complete old environment atom accounting,
source-coordinate retention, state parity and unchanged total charges. Four
additional real-input checks pass: paired Ca/La coordinates, only HG1 moving,
new physical cap Jacobians, and byte-identical expanded A/B fields. See
`CHECK_RESULT.json`. Coordinate finite differences are an algebra check,
not molecular force qualification. No scientific integration test ran.

Reproduce to a new directory from the repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/metal_environment_response_20260926/preparation/expanded_region/prepare.py \
  --repository "$PWD" \
  --inputs workspaces/metal_environment_response_20260926/preparation/scout_v3/INPUTS.json \
  --output workspaces/metal_environment_response_20260926/preparation/expanded_region_reproduction
```

Writers refuse overwriting an existing output directory. Preparation took about
nine local wall seconds; process CPU/RSS were not instrumented. Zero molecular
calls and zero allocated compute. Root owns whether/when an electronic test is
justified after the ongoing native SCF diagnosis; this task submitted no jobs.
