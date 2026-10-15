# Asn290 terminal-amide preparation: actual result

Slurm **4614 completed all 10 sources / 20 Ca–La cells** under the frozen [rotamer plan](../ASN_ROTAMER_PLAN.md). All five Ca-conditioned and five La-conditioned A0ACD6B9F2 sources remain in the denominator, with all 60 original shared-pool cells retained as baseline evidence. Preparation made **zero energy calls**; admission does not establish a favorable candidate or improved classification.

## Construction and admission

Each source uses its actual CB–CG axis to rotate OD1, ND2, HD21 and HD22 together by exactly 180°. The emitted CG and every other coordinate line remain unchanged. The paired Ca/La endpoints retain the dry 204-atom membership, original maps/protons, Ca −4 / La −3 singlet states, all 14 caps and their fixed anchors. This is a different terminal-amide conformer, not an O/N name swap.

All sources passed the frozen covalent, rigid-fragment, source-identity and paired-coordinate checks. The target Ca-conditioned sample 4 has maximum heavy displacement **2.317802402662702 Å**: this declared discrete construction is outside the old 0.8 Å search domain. The old displacement gate was not claimed to pass, and no continuous rotation path was admitted.

Each contact check used all **9,529 physical source atoms**, including H, in an exact-identity union of **9,546 centers** with the local context and separate synthetic caps. Moving source atoms were not duplicated as stale exterior centers. There were no newly introduced contacts below the inherited severe-overlap limits: **1.0 Å for heavy pairs / 0.55 Å when either atom is H**, using `new < limit AND origin >= limit`. These are extreme-overlap gates, not a low-strain criterion.

## Contact limits that matter for scoring

The ten saved nearest-omitted-source summaries place the post-flip closest heavy pairs at **4.4756–5.2747 Å**, and H-involving pairs at **2.7412–3.8994 Å**. For target Ca-conditioned sample 4, OD1–omitted Gly289 CA changes **6.147382532427929 → 5.232481319812116 Å**; HD21–omitted Gly268 HA3 changes **6.029637965914703 → 3.604474525771083 Å**. The full source supplies geometry checks only; its energy is not added to the fixed 204-atom calculation.

Some much shorter H contacts are **inside the scored context**. For example, Ca-conditioned sample 2 introduces Asn290 HD21–its backbone H at **0.9610772821283347 Å**, versus **4.701552190500495 Å** at origin. Such an unrelaxed H–H contact passes the 0.55 Å gate but may carry a substantial energy penalty. Conversely, the target's nearest H entry is the **unchanged covalent ND2–HD21 bond, 1.1827345433359084 Å**. The saved nearest-contact fields include bonded pairs: the roughly 0.96–1.18 Å local entries must not all be called clashes or newly introduced contacts. No geometry or threshold was adjusted after these observations.

## Evidence, cost and next stage

- [Actual preparation](../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_preparation_v1/execution/PREPARATION.json), SHA256 `508e7594aea7ad6198c62ce931d84a5f19d72b71070edcdbc782787f8cd9ca30`, pins every per-source admission receipt and endpoint.
- [Execution receipt](../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_preparation_v1/execution/EXECUTION.json): one allocated CPU on `node-64-768g-3`; **5.166875166003592 s** process wall time, **4.195693803 CPU-s**, peak RSS **93,748 KiB**. These process measurements are not the scheduler's full billed duration.
- Frozen manifest SHA256 `ae4184400f6908039ebb1eaa2121181e92510738d88366d5175b9ecc9df624fb`; preparation implementation SHA256 `45edb9a955782b67b07f8627e42e9e05a4f92c12fb2112fa1a620bcb053e111c`.

Root submitted the separate 20-cell MACE stage as **4617**, with its completion wake armed. At this report's saved-output check, `asn_rotamer_scoring_v1/execution/FINAL_COLLECTION.json` was absent. Energetic preference, common-pool changes and the conditional solvent-stage decision therefore remain unreported here. Production and calibration are unchanged.
