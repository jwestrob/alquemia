# Same-metal orbital initialization helper — 28 September 2026

Completed frozen-f Hans orbitals can now be staged as an explicit initial guess for a nearby geometry. This does not reuse the source energy or qualify the new endpoint; new SCF convergence, gradient admission and execution receipts remain necessary. No QM/Slurm evaluation or shared-runner modification was performed.

## Interface

`scripts/metal_environment_orbital_seed.py` exposes:

- `validate_seed(source_manifest, collection, source_task_id, target_manifest, target_task_id)`
- `stage_seed(source_manifest, collection, source_task_id, target_manifest, target_task_id, directory)`

The first returns initialization-only provenance. The second creates a **new** directory containing `initial.gbw` and `SEED.json`. It verifies source and copied hashes, refuses overwrite, and rechecks output/receipt/manifest/collection pins after copying. Failed staging is not a reusable seed; no `SEED.json` means no admitted seed.

Require an accepted source collection row and a normal, converged, zero-return execution receipt whose immutable output matches the accepted energy. GBW access occurs only after these checks. Reject absent/empty GBW or modification after execution finished. Historical receipts did not hash GBW; the helper pins the retained completed file now rather than claiming a contemporaneous GBW pin.

Require the same executable pin, complete method string, same-metal basis/ECP and auxiliary hashes, exact electronic-state metadata, ordered elements/source IDs/cap identities, and boundary mapping. Only geometric fields/Jacobians are excluded from mapping-identity comparison. Geometry and point-charge changes are recorded explicitly. The preparation/runner still must admit the target Hamiltonian and geometry; this helper is not a replacement for that validation.

Frozen-f Dy uses an effective singlet valence representation and retains physical sextet metadata. This helper rejects changing that representation, state, or metal. It is not the earlier ECP28 explicit-4f sextet path.

## Runner integration contract

Retain `NoAutostart` in the method. In the target SCF block replace the initial-guess directive with:

```
Guess MORead
MOInp "initial.gbw"
```

Copy the pinned `initial.gbw` into the actual execution directory before launching; pin it in execution/cache provenance. Do not append a conflicting second Guess directive. A changed initial guess is a new execution attempt, not permission to reuse a successful target energy cache entry. There is no automatic fresh-guess fallback. The new target must independently converge and pass existing scientific gates.

This syntax was checked against the [ORCA 6.1 initial-guess manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html), retrieved 28 September. No seeded ORCA run has yet qualified convergence, speedup, or preservation of the converged electronic solution.

## Real-artifact verification

`python -m pytest -q tests/test_metal_environment_orbital_seed.py`: **9 passed** (6.70 s).

Fixtures: `workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2/{manifest.json,FINAL_COLLECTION.json}` and their completed GBW/output/receipts. Both La_A→La_B and Dy_A→Dy_B stage successfully; max movement is 0.0829489053 Å and the environmental point-charge file is unchanged. Malformed copies reject cross-metal transfer, altered physical state, method, basis/ECP, executable, source identity, and unavailable source status. These are interface tests, not new molecular evaluations.

Manifest SHA256: `708ba36d326f1435f1fe6d2a09e455fce0d52b2a6083bab9e1dc3d7e31627925`.
La_A GBW SHA256: `baa621eb0b665ee8973c7b77f2ac5bfb6b53162257dfbe550da6f0e59658a945`.
Dy_A GBW SHA256: `49d05fadfa678efbbaf7c75a57cee0620f567e8e879921a825c515ff446eaca9`.

Runnable validation-only example (repository root):

```bash
python scripts/metal_environment_orbital_seed.py \
  --source-manifest workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2/manifest.json \
  --collection workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2/FINAL_COLLECTION.json \
  --source-task-id Dy_A \
  --target-manifest workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2/manifest.json \
  --target-task-id Dy_B
```

The B target above is already complete: this example validates compatibility only, and does not authorize rerunning it. Root integrates with a separately declared future displacement scout.
