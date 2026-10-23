# Nine baseline states at a 12 Å Ewald cutoff

This package implements the prospective ROOT_DECLARATION exactly. It retains the failed active panel and the failed eight-endpoint, 9 Å baseline controls. It prepares one new center and the same eight whole-box/strain-1 endpoints, with local response absent.

Only the EWALD-CUTOFF key token changes from 9.0 to 12.0. Literal XYZ, parameters, all other key bytes, steps, alpha 0.4, 64³ mesh, interpolation order 5, and original six derivative gates remain unchanged. Runtime keys retain the exact original observer suffix.

The qualification harness has one source constant change: `ewaldcut/=9d0` becomes `ewaldcut/=12d0`. BUILD.py compiles that GATE source once and links one old-model executable, reusing the admitted 20239 GRADIENT and STATE_WRITER objects. Their actual source/compile/object joins are verified. Compiler, flags, module inputs, libraries and loader checks remain pinned. No native model kernel changes. The original objects' public gate interface is unchanged.

RUN.py differs from the original NATIVE.py only in its context import. All run, invocation, input staging, parsing, process limits and counter behavior remain unchanged. The context restricts execution to nine fixed inactive rows and gradient-only operation. Every row has one attempt, a fresh immutable output, 1 CPU, 3 GiB address-space cap and 3600 s CPU/wall limits. Build phases each have 180 s CPU/wall limits. The original one-node, one-task standard_expand wrapper supplies default scheduler memory/time and no-requeue.

The build needs a root source review and exact release. After its successful result and separate root build admission, BIND.py writes the nine native configurations and unbound release templates. It launches nothing. Each native call still needs its exact root release. No initializer, factor, automatic retry, extra center reuse, or candidate qualification is authorized here.

The separate saved collector will evaluate the six original two-direction finite-difference/two-step checks using the newly computed center. It will compare with the admitted 9 Å results descriptively and preserve failures. Changed cutoff is a numerical diagnostic, not proof of PME convergence, a custom-response qualification, physical accuracy, or dynamics readiness.
