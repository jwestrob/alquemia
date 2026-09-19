# Native minimum checks, frozen before these outputs

The optional interior, stable cheap minima selected in StageB receive native
DFT EnGrad. Energy error must be ≤max(0.05kcal/mol,10% of the actual DFT energy
change). Maximum physical-coordinate gradient at the native endpoint must be
≤0.4kcal/mol per Å/radian (the energy floor0.02kcal over the0.05 test step).
Stationarity remains necessary, not sufficient, for harmonic basin free energies.
This is a diagnostic acceptance rule, not a fitted biological classifier.
Negative/soft modes are retained, and minimum trust-boundary failures do not
qualify by lowering their gradients or increasing the allowed radius afterward.
The original radial and coupled local-direction criteria remain unchanged.
