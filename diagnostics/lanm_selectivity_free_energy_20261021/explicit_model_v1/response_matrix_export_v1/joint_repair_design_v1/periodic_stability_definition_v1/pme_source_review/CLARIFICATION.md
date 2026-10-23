# Intended discrete operator versus numerical application

This clarification preserves the original source report and its pins. Its unresolved FFT/application-error statement must not be read as withholding the conditional mathematical PSD proof.

Define C from the stored real interpolation/conversion coefficients and Q from the stored multiplier. If scatter and gather represent C and its transpose, Q is nonnegative and conjugate-even, and F denotes the exact unnormalized DFT, the intended discrete operator `C^T F* Q F C` is mathematically real symmetric positive semidefinite. That statement does not depend on a forward-error proof of the floating-point FFT implementation. Outward bounds on its stored-weight norm and on the other intended operator terms can certify that mathematical matrix.

Finite arithmetic used to apply the operator, solve the induction equations, or evaluate energies and derivatives is a separate numerical accuracy question. It must be controlled for a reliable implementation, but should not be conflated with the mathematical PSD statement. Actual non-even stored Q and unequal stored scatter/gather weights are different: they change the intended operator and must be rejected or explicitly included in its analysis. The source-only review implemented neither a correction nor a periodic certificate.
