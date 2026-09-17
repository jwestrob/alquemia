# Saved-output audit, declared after observing checkpoint disagreement

Standing pilot authorization applies; no new model or DFT evaluation. The large
pilot completed, passed rotation/charge checks, and its raw full-system contrast
differs from medium by about 671 kcal/mol. This is a consumed development case.

Audit exactly the completed analytic medium and large primary La/Ca pairs on
the same 9,141-atom 1H4I geometry. Verify all source identities and coordinates.
Compare the recorded energy components and their remaining additive reference
term, without assigning a unique physical cause to a learned component.

Export every physical atom's predicted endpoint charges, force magnitudes and
direct-model grad(E_Ca-E_La) = F_La-F_Ca, in eV/Angstrom. Summarize charge response
and squared gradient norms in fixed metal-distance shells [0,5), [5,10), [10,20),
[20,40), [40,infinity) Angstrom. Also report the ten largest direct contrast
gradient norms and ten largest checkpoint force changes, with nearest-atom
distances/identities. These descriptive selections have no decision threshold.
They are not independent validation or evidence of a specific error mechanism.

The direct MACE gradients are not gradients of the subtractive DFT/MACE hybrid;
matching DFT core gradients would also be required. No relaxation estimate,
uncertainty covariance, entropy, new label or numerical score is added. Save
unrounded arrays and a traceable all-atom table under workspaces/.
