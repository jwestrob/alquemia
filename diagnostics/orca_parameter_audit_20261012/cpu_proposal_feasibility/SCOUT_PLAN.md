# Declared CPU force/proposal parity scout

12 October 2026; declared before model loading or scientific calls. Root explicitly
allocates eight total active numerical threads under Jacob's temporary <=16
login-thread authorization. Parent uses one BLAS thread while the model child
waits; child uses eight Torch threads while the parent blocks. No concurrent
parent/child numerical work, other model workers, GPU, Slurm variables, package
installation, production edits, solvent or DFT calls.

Question: can the unchanged expanded-context native-vacuum OMOL force and
four-angle SLSQP proposal method run with a persistent CPU calculator and the
original optimizer binary? Energy-only CPU agreement is prior evidence only.

Exact consumed fixture is A0A3F2YLY8, Ca-conditioned seed1/sample1, envelope
e05e200a13d024a5:202 atoms, chargesCa−3/La−2, multiplicity1. Actual archived Ca
proposal reaches0.8Å boundary; La is interior at0.678478Å. All XYZ/map/force/
receipt pins come from pilot34_searches_v1/manifest.json and its actual results.

1. Five fresh energy/analytic-force calls, one loaded CPU model: Ca origin, La
   origin, Ca saved proposal, La saved proposal, Ca origin repeated after all
   state/geometry switches. Verify model/head/state/order and finite forces.
   Energy difference <=0.01kcal/mol and maximum Cartesian force component error
   <=0.001eV/Å, the existing scripts/mace_omol.py TOL precedent. Also compare
   repeat to first CPU result under those same gates. Project both fresh origin
   forces through the unchanged Kinematics/project/preview code and require
   the same ordered four mode IDs, criteria and redundancy decisions. Record
   numerical load differences without misusing the historical1e-12 copy-only
   tolerance as a new device-parity threshold.
2. Only if every phase1 check passes, exactly two SLSQP starts on these archived
   tasks, q=0, one Ca and one La. Reuse the now-qualified exact archived origins
   as the historical CompletionProposal implementation does; all nonzero trial
   energy/forces are fresh CPU calls and final cache reuse remains explicit.
   Preserve SciPy1.17.1 original parent kernel, four selected angular modes,
   analytic gradients/cap chain rules, shifted Hartree-equivalent objective,
   ftol1e-8, maxiter200, angular bounds±0.8rad and heavy-displacement0.8Å.
   Preserve all existing final geometry/energy/admission/boundary checks.
   Compare final energy to archived proposal with the existing precision
   proposal-energy gate0.001kcal/mol; record coordinate/full-q/gradient differences
   and require the same boundary/admission states. Do not require identical
   iteration trajectories or fabricate coordinate identity.

Finite scope is five replays plus two unchanged optimizer starts. **No new call
or wall-time ceiling**, fallback, retry, optimizer tolerance change, new start,
outcome-based selection or107-case expansion. The archived starts used15+13
new trial calls; this is a historical observation, not a promised CPU call count
or speed prediction. Record actual calls, failures, both process CPU times,
wall time, load cost, peakRSS and all force/results as fresh immutable files.

Transport is the only new execution component: original lanm_qmmm optimizer
parent, existing MACE venv CPU worker, evaluate(request_pin)->(result,pin).
Parent and model environments have different SLSQP binary hashes; therefore do
not move the optimizer into the model environment. Unchanged historical source
snapshot modules are imported directly. Own signals terminate/reap only the
owned worker process. Failure receipts and unfinished stages remain explicit.

Root owns the later scientific decision and any larger campaign.
