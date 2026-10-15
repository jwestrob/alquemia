# Flipped-basin runner: independent static review

Reviewed `RUN.py` SHA256 `82f78a22e1e6f342bc4c3a0d1443310a2607b867bf4417ec8c07ed093c58656b` and `RUN.sbatch`. No remaining launch blocker identified. This is source review, not executed numerical qualification; root owns the actual Ca-conditioned sample0 scout.

- The starting energies and forces come from the actual fixed-flip receipts. Before model startup, the runner checks the exact starting coordinates and reprojects the saved forces through the rebased Jacobian. It retains the original four mode identities/order and the original private `completion.optimize` implementation, settings, objective units and SLSQP runtime pins.
- The derived evaluator adds the declared full-source contact admission before delegating to the original cached evaluator. Final accepted proposals require the unchanged optimizer/physical gates and their actual full-source admission. Own-metal final results are reused; opposite-metal calls retain the same physical coordinates, atom order, charge and multiplicity policy.
- Both successful searches and both cross-metal results are needed for the complete common6 pool. All ten source rows remain visible, including unexecuted or failed rows. Common3/4/6 use the existing chooser, original-first ties, separate mathematical/operational results and strict −0.1 kcal/mol operational rule. The continuation flag uses the absolute mathematical common4→6 paired change, without class bands.
- The initial stage contains only the declared first source. Remaining sources require the actual completed scout collection. Exclusive output creation, original interrupted-search handling and final execution/failure receipts preserve missing work rather than silently restarting it.

Two concrete issues were corrected before binding: guard evidence now lives outside the evaluator-owned cache directory, avoiding its exclusive-directory-creation collision; preflight now verifies the actually imported SciPy version and SLSQP wrapper/kernel hashes rather than only copied metadata.

No imports, numerical tests, model calls or submissions were performed for this review. The finite scout must establish actual runtime behavior and cost; source review does not establish an equilibrium rotamer, composite minimum or predictive improvement.
