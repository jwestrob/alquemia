# Exact 1H4I additive component ledger — preparation only

The proposed finite additive bookkeeping is implementable with the recovered
parameters. This is **not yet a qualified hybrid force model** and, with no
solvent, is **not a solvent-consistent whole-protein relaxation model**. No
optimization should follow this preparation. Root's `../LEDGER_PLAN.md` governs
any later component evaluations.

## Defined candidate energy

E = E_embedded_QM(real QM + mapped caps; qtilde_MM, DoEQ=false)
  + E_FF_bonded(terms with at least one real MM atom)
  + E_LJ(MM,MM) + E_Coulomb(qtilde_MM,qtilde_MM)
  + E_LJ(real QM,MM).

All QM–QM classical terms are omitted. No independent cap FF particles, cap LJ,
classical QM–MM Coulomb, C4 induction, second solvation term, or bulk-ion offset
is added. The embedded QM term includes its existing r2SCAN-3c terms; no second
QM dispersion term is added. MM-only and cross bonded terms use the original
OpenMM force objects/parameters, including actual CMAP tables in pinned XML.

[ORCA's additive QM/MM documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/multiscalesimulations/qmmm-molecules.html#additive-qmmm)
assigns MM internal and boundary bonded/LJ interactions to the force field.
That supports the decomposition, **not equivalence of our precise link/charge
scheme to native ORCA**. This ledger is an explicitly named candidate.
The capped electronic energy still approximates the covalent boundary;
retaining real cross bonded terms does not prove an exact cap cancellation.

## Concrete accounting decisions

- 9,141 physical atoms: 9,090 MM and 51 real QM. Three mapped caps make the
  electronic input 54 atoms; caps are not extra mechanical particles.
- The three omitted point-charge boundary CA atoms remain real MM atoms with
  their native LJ and bonded interactions, but qtilde=0. Removing those atoms
  from mechanics would sever the scaffold and omit repulsion.
- Retain 9,225 bonds, 16,655 angles, 23,666 torsions and 597 CMAP terms. Omit
  20 bonds, 27 angles and 24 torsions supported wholly within real QM.
- All 49,903 native protein exceptions are classified: 49,765 MM–MM,
  67 cross, 71 wholly QM. Cross LJ exceptions retain exact native pair sigma
  and epsilon. Wholly QM exceptions are omitted with their classical pairs.
- MM–MM exception **charge products must change** after redistribution. The
  builder derives bond distances from actual native bonds, verifies original
  products against the XML's native 1–4 scale, then forms scale*qtilde_i*qtilde_j.
  Original exception products would be inconsistent with the new MM Hamiltonian.
- MM electrostatic exceptions are not blindly applied to quantum embedding.
  The electronic embedding already uses its declared redistributed charges;
  an additional classical cross Coulomb term would double-count that interaction.
- PQQ and metal have no covalent links to MM, so all their real QM–MM pairs use
  nonbonded LJ. PQQ charge placeholders are never read as charges.
- Finite nonperiodic NoCutoff convention; no PME, reaction field or dispersion
  tail correction. Native exception parameters replace generic pairs.
  [OpenMM specifies Lorentz–Berthelot mixing and exception replacement](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#lennard-jones-interaction).
  Amber Rmin/2 converts to sigma_nm = 0.2*(Rmin/2)_A / 2^(1/6);
  epsilon_kJ/mol = 4.184*epsilon_kcal/mol.

## Explicit unqualified metal parameter choice

Installed `frcmod.ions234lm_126_tip3p`:

| Metal | Rmin/2 Å | epsilon kcal/mol | File provenance |
|---|---:|---:|---|
| Ca2+ | 1.649 | 0.10592870 | CM set, Li et al. JCTC 2013 9 2733 |
| La3+ | 1.718 | 0.15060822 | IOD set, Li et al. JPCB 2015 119 883 |

These are the declared pure 12–6 candidates, not C4-stripped 12–6–4 fits.
Their hydration parameterization does not establish embedded-QM transferability.
The full file is hash-pinned in `LEDGER.json`. No parameter was chosen by a
computed classification. PQQ types/coordinates and protein terms remain those
of the exact preparation already reported.

## Artifacts and implementation

`scripts/metal_environment_mechanics.py` prepares the ledger and rejects an
unsupported source size, altered mapping, unsupported exception graph/scaling,
PQQ physical charges or mismatched coordinates. It does not create an OpenMM
Context or evaluate any energies or forces.

Workspace: `workspaces/metal_environment_response_20260926/hybrid_preparation_v1/ledger_v1/`
contains `LEDGER.json`, `particles.json`, `exceptions.json`, `bonded_terms.json`,
and `caps.json`. Core indices, source identities, cap Jacobians, native force
indices and source-file pins remain available for physical force projection.
The cap Jacobians are origin values: moved configurations must regenerate them
from the actual link geometry, retaining the archived serialization convention.

The first preparation process was stopped by its owner before writing artifacts:
an unnecessarily quadratic omitted-term count was replaced by direct support
classification. This was code preparation only; no molecular work was repeated.

## Proposed finite cheap checks — not executed

After root freezes the execution plan, instantiate separate OpenMM force groups
for native retained bonded terms, MM LJ, MM Coulomb and cross LJ. Prefer native
bond/angle/torsion/CMAP objects; use declared NoCutoff pair groups and explicit
exception replacement. No new physics engine is needed.

1. Count every real pair/term once; check all 67 cross exceptions individually,
   and all changed MM charge products against native scaling.
2. On exact real A and B, evaluate components separately; report identical
   metal-independent MM terms and the metal-LJ contribution separately. These
   are diagnostics, not improved response labels.
3. Repeat and jointly rigid-transform the full real system; compare every
   component energy and transformed force. Do not hide an electronic gate
   failure beneath small classical residuals.
4. Directional finite differences in the already declared MM hydroxyl and
   mapped boundary directions, plus one real metal direction, check each
   classical force group. This qualifies derivatives of the declared model,
   not biological or metal-parameter accuracy.
5. Combine gradients only after cap chain-rule mapping and actual MM embedding
   gradients are included. Original electronic rotational failures remain
   visible; a classical success cannot retrospectively pass those gates.

No solvent has been supplied by this exercise. Later relaxation requires solvent
in the same force-generating energy and a separately qualified model; do not
attach the old GB scalar correction or call finite dry mechanics that model.

## Reproduction and completed tests

From the repository root, use a new output directory:

```bash
python scripts/metal_environment_mechanics.py \
  --source-manifest workspaces/metal_environment_response_20260926/preparation/scout_v3/INPUTS.json \
  --protein-result workspaces/metal_environment_response_20260926/hybrid_preparation_v1/protein_v2/RESULT.json \
  --pqq-cross-lj workspaces/metal_environment_response_20260926/hybrid_preparation_v1/pqq_types_v1/PQQ_CROSS_LJ.json \
  --metal-parameters /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/dat/leap/parm/frcmod.ions234lm_126_tip3p \
  --output workspaces/metal_environment_response_20260926/hybrid_preparation_v1/ledger_reproduction
python -m pytest -q tests/test_metal_environment_mechanics.py
```

Four real-artifact parameter-accounting tests pass (4.27 s). These verify physical
partition/caps, changed exception products, retained bonded support, absent QM
classical charges and rejection of a mislabeled metal parameter family. They
are **not** molecular energy/force qualification tests. Zero molecular calls.
