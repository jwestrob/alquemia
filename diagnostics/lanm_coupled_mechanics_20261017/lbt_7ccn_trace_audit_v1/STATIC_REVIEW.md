# Independent static review — saved trace audit

**Verdict: pass for the declared single saved-only Slurm analysis.** No material source, provenance, arithmetic-order or scope blocker was found. Actual replay remains untested until that analysis runs; no additional preflight is requested.

## Frozen bundle

| Artifact | SHA-256 |
| --- | --- |
| `AUDIT.py` | `eeb09ad815e2a9ca76df4329c08b16b91ed71c5e61e07bf6778c4b5dc44127ab` |
| `SPEC.md` | `d8c519e4a5719e25a292fe765a4efc89749a26136ce2a605140aa88439ec1dcb` |
| `CONFIG.json` | `6c9b4f90768001c69d0c48a9a2da7a0b49f6636e2ddef6f6307408111c7f360c` |
| `READY.json` | `45ed68c15292ba8d847a7126a5156893219040331ddf3570a657bb5946a6d262` |
| Root science declaration | `9151692784e2e340e6967ef3d2baccb38a07208b84b168492baf82a40b88d01c` |

Review checked source/AST, literal saved JSON identities and174 selected file/hash pairs. All81 configured point/checkpoint pairs and80 request links match actual15458. No scientific module was imported, numerical replay executed, molecular call made or job submitted. Only this review file was written.

## Admission and zero-call scope

Pinned completion helper `551d964b.prepare` receives the actual run configuration and preserves the original source/state/model/solvent, previous cache and derivative-admission lineage. It must return with no initialized model, no requests and zero native costs. Before trace processing, both evaluator entry points are replaced with explicit failure tripwires. No optimizer search, alternate trajectory or native engine invocation is introduced; final zero-call assertions remain.

The denominator is fixed at all81 states92–172 and all80 transitions, with four predeclared20-transition windows plus the overall window. Every checkpoint retains the original initial point, settings and exact historical request/evaluated-point prefixes. The final checkpoint's numerical state equals CHECKPOINT_172; only the extra final no-trial trace row differs, and the code separately checks that complete trace. Partial failures preserve unavailable slots and replay diagnostics.

## Replay arithmetic and outputs

- The code follows the original operation order: rigid projection, original L-BFGS helper, direction projection, descent fallback, maximum-atom trust scaling, accepted multiplier/coordinates, native Armijo work, projected secant pair, curvature/history admission and truncation, then trust update.
- Each transition starts from its own archived checkpoint, avoiding accumulation of replay roundoff over an artificial trajectory. Next-checkpoint history/trust, force metrics, rigid removal and accepted coordinates are checked against saved evidence. H0 scale, clipping, fallback, curvature and actual/predicted decrease are diagnostics, not modified optimizer controls.
- Original gradient signs and conversion association are preserved. Full component-gradient closure and native endpoint work remain exact gates. The actual native solvent keys are `Gelec`, `Ghb`, `Gsasa`, `Gshift`, matching the implementation. Window work comes from native endpoint differences; per-step sums and solvent parts receive separate closure checks.
- Coordinate replay uses absolute1e−10Å with no relative allowance; other derived values use1e−10 absolute plus1e−10 relative; work closure uses absolute1e−10 kcal/mol. These are frozen recomputation allowances, not force-accuracy bounds or relaxed scientific convergence criteria. Failed comparisons are recorded before admission stops.
- Attribution preserves all256 IDs and signed forces, groups residues by their physical residue identity, and keeps La separate. Component/residue convergence booleans from the generic metric helper do not define global model convergence. Endpoint gradient projections are local linear diagnostics, not integrated work or a finite-difference test; optimizer secants are not a physical Hessian.

## Runtime and interpretation

The actual staged launcher matches READY: original MACE Python environment, qualified `node-224-2t-8gpu-1`/`gpu_h200`, one CPU, no GPU request, all thread controls1, no memory/time override and no requeue. It executes the full supplied argv. Its literal `standard` QoS header requires root's declared `--qos=standard_expand` submission override.

The report keeps physical accuracy and convergence distinct from replay success. Budget exhaustion, solvent-dominated lowering, force oscillation or trust scaling alone do not prove stagnation or a defective interaction. The fixed audit can identify tested bookkeeping mismatches and characterize optimizer behavior; it cannot independently validate terminal forces or authorize another optimization/reference. No source changes or additional scientific tests are requested by this review.
