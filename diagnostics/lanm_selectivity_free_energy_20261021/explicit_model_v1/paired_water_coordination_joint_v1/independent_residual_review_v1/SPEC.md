# Independent terminal residual interpretation

Read only completed fit17676, selected-ten raw review17683 and residual summary17686.
Check their exact result joins and reconstruct the last normalized Jacobian from
the twelve stored derivative residuals and actual parameter spans. Independently
use an SVD projection at the frozen rcond1e-6; compare the saved rank, direction,
loss and origin/selected identity. Inspect finite caps, parameter/trust contacts,
and the actual complete final trials. No native call, objective reevaluation,
candidate change, fitting, fresh validation outcome, or extension is authorized.

The question is whether the failure is explained by a simple iteration/bound
limit or by residual components outside the local six-parameter response span.
This is a finite-difference, weighted engineering-residual diagnosis. It cannot
prove a global parameter-family limitation, model accuracy, or a unique omitted
physical term. Use one local numerical thread and retain actual arithmetic costs.
