# Nikasha PLM candidate: technical supplement

## Protocol identities

| Layer | Pinned identity |
|---|---|
| Operations package | `Nikasha_PLM_three_source_envelope_operations_v1` |
| Source preparation | `fixed_core_PQQ_source_to_complete_context_v1` |
| Envelope preparation | `Nikasha_explicit_three_La_motion_envelope4p3_preparation_v1` |
| Executor interface | `Nikasha_explicit_three_La_envelope4p3_strict_native_execution_v2` |
| Reference | `Nikasha_motion_envelope4p3_adaptive_strict_native_canonical25_v1` |

The new operations layer runs the existing source preparer before the existing
envelope planner; it introduces no new score, Hamiltonian or decision bands.
The previously executed two-protein molecular integration uses interface v1;
v2 removes its fixed two-protein scope and preserves the scientific recipe.
Fresh preparation and v2 execution preflight have been checked on those same
real sources. There was no new molecular endpoint calculation in this packaging
task. The historical integration and new preparation receipts remain separate.

## Preparation and motion

- Exactly three declared La-conditioned sources per protein, exact source
  assembly, matching sequence/numbering/role definitions. Different proteins or
  folds never compete by absolute total energy.
- Frozen seeded standard-protein protonation at pH 7, one CPU thread; retained
  source heavy atoms. Missing heavy-atom repair is unsupported. PQQ remains
  `C14H3N2O8`, formal charge −3. The original dry core and water inventory remain
  fixed. The generic amide repair and contextual-water policies are separately
  versioned; this dry PQQ protocol does not acquire variable water occupancy.
- Fixed union of complete polar fragments selected within 4.3 Å across the three
  sources, using the existing source atom/bond/cap graph. This is the established
  3.5 Å contact distance plus the 0.8 Å admitted displacement. It is not a
  converged global electrostatic boundary or a whole-protein mechanical model.
- Source-connected whole-donor torsions; synthetic caps follow their mappings.
  Force ranking uses heavy-atom displacement normalization and removes redundant
  directions. Both metals use the same selected active subspace. PQQ and metal
  position are fixed in this branch.
- Up to four angular coordinates; ±0.8 rad bounds and maximum source-heavy
  displacement 0.8 Å. SLSQP, maximum 200 iterations. Shifted vacuum-MACE objective
  scaled by `0.03674932217934773`; operative `ftol=1e-8` is on that
  Hartree-equivalent objective, not an unscaled eV tolerance.
- Existing contact/chirality/final-admissibility checks remain active. A search
  touching a bound is reported; it is not declared an unconstrained minimum.

## Energy and selection

For each metal M and candidate geometry q:

```text
E_M(q) = E_MACE-OMOL,vac,M(q)
       + E_native-GFN2,ALPB(water),M(q) − E_native-GFN2,vac,M(q)

Q = {original, Ca-proposed, La-proposed}
R_math = min(q in Q) E_Ca(q) − min(q in Q) E_La(q)
```

Each row scores the same candidate set. The operational variant retains the
origin for improvements below 0.1 model kcal/mol; otherwise it selects the
lower-energy candidate according to the existing tie policy. Both variants and
each candidate's native/solvent/composite work are retained. Any required failed
matrix cell makes the pool unavailable. Candidate counts are not populations.

MACE eV is converted once by `23.06054783061903` kcal/mol per eV; native GFN2
Hartree once by `627.509474` kcal/mol per Hartree. Solvent correction to the
paired contrast is `delta_solv,Ca − delta_solv,La`. The baseline MACE term is
vacuum, so this does not add a second CPCM energy to a DFT endpoint. This remains
a mixed-model descriptor and supplies no water chemical potential or entropy.

The protein score is the median of all three operational contrasts. The complete
member range and mathematical variant remain visible. Larger R is more La-like.
The canonical-only operational bands are `Ca_max = −405463.5325440529` and
`La_min = −405456.4611478038` model kcal/mol; the interval is inconclusive.
Reference SHA256:
`4a1d3ae6e83b5b6ace408f729c39390d44f673554e338df932fd38fc971f67d9`.
Historical DFT aquo-referenced S and raw composite R have different gauges.

## Executables and numerical settings

- MACE checkpoint: `MACE-omol-0-extra-large-1024.model`, head `omol`, float64,
  nonperiodic, model cutoff 6 Å; physical total charge and multiplicity inputs,
  no added dispersion term. SHA256:
  `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`.
- Recorded model environment: `mace-torch 0.3.16`, `torch 2.8.0`,
  `graph-longrange 0.4.4`; environment/checkpoint files pinned in the release.
- ORCA 6.1.1 native GFN2-xTB, installed native parameters and mixer, fresh
  `NoAutostart`, `TolE 1e-10` Hartree, `SmearTemp 300`, `MaxIter 500`, one rank
  per scalar task. Native electronic smearing is not a configurational entropy
  correction. No silent standalone-xTB substitution.
- Closed-shell Ca/La source endpoints, charge differing by +1 on substitution;
  atom order/coordinates/maps and formal-charge accounting are checked.
- Runtime profile: one H200, 32 CPUs, 200000 MiB host RAM; scalar calls execute
  independently within the allocation. MACE calculator stays warm within a batch.
  All exact configuration and qualification files are pinned in
  `params/pqq_plm_envelope_v1.json`; each execution snapshots its implementation.

## Evidence denominators and limitations

| Comparison | Completed finding | Interpretation |
|---|---|---|
| Tenfold-context candidate, 225 alternative structures / 25 groups | 207 correct, 0 wrong, 1 inconclusive, 17 unavailable | Strongest structure-transfer result; different context policy from routine three-source route |
| Same-context/same-solver static, same 208 supported structures | 199 correct, 1 wrong, 8 inconclusive | Accommodation adds useful information under this development comparison |
| Practical three-source primary, 100 correlated triples | 91 correct, 9 unavailable | Six prior preparation exclusions plus three triples sharing one failed scalar cell |
| Separate qualified two-start recovery sensitivity | 94 correct, 6 unavailable | No new band/geometry; not an automatic retry policy |
| Practical individual source/context pairs after recovery | 100 correct, 4 inconclusive of 104 | Do not advertise the tenfold single-ambiguity claim for this route |
| Two unlabeled PLM proteins / six sources | Both median calls Ca-supported; smaller ranges | Physical response / application feasibility, not independent accuracy labels |

The 25 original calibration inputs and three consumed crystal controls retain
their known classes under the appropriate reference. These are development and
transfer observations, not newly blind experiments. Three-source ranges increase
in 52 of 91 primary matched triples relative to the same-context origin; universal
spread reduction is not supported. The separate Ca-conditioned A0A3 sample3 probe
becomes inconclusive. Those limits remain alongside median reference fidelity.

The earlier actual two-protein integration used 183 MACE evaluations and 72
native GFN2 calls in 251 allocated seconds on one H200/32 CPUs. Source handling
took another 16.955 seconds and reused archived protonation. The new packaging
task separately measured fresh source/envelope preparation; its receipt is in
the delivery report. Do not combine these separate runs into a newly measured
end-to-end throughput or a matched DFT speedup.

Archived native DFT donor-displacement checks support limited native-MACE
response directions and approximate magnitudes. The solvent-added response
disagreement persists after tighter scalar convergence, and neither force
consistency of the full composite nor a complete binding free energy is claimed.

## Biological join and archive

The original 176-protein XoxF export retains 200 source-gene mappings and its
existing phylogenetic/transcript records. It is not automatically the final
all-PQQ cohort. Extend that export from real crosswalks when additional families
are added. Preserve original DFT fields, unavailable expression and the existing
normalization definitions; the candidate join creates no new expression matrix.

Keep source manifests, preparation ledgers, execution/result JSON, all failed
endpoint receipts, canonical reference, exact software/model pins, and the
final table receipt with the supplement archive. The SOP provides runnable
commands and immutable output conventions. Full development history remains
available through the linked reports; it is not necessary to reproduce it in
the main text.
