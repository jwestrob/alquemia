# MACE local curvature screen

Exact archived physical GGR displacements; zero new DFT. Medium primary, large sensitivity.
Candidate curvature is MACE vacuum + frozen-monopole OBC-II, with charges reevaluated at each geometry.
No combined analytic gradient, relaxation correction, stable-basin or biological-accuracy claim.

| Model | Representation | Motion | Metal | DFT even | MACE+GB even | Difference | Curvature gate |
|---|---|---|---|---:|---:|---:|---|
| medium | extended | metal | La | 0.112904 | 0.111455 | -0.001449 | True |
| medium | extended | metal | Ca | 0.049819 | 0.050272 | 0.000453 | True |
| medium | extended | metal | R | -0.063085 | -0.061183 | 0.001902 | True |
| medium | extended | peptide | La | -0.010722 | -0.011348 | -0.000626 | True |
| medium | extended | peptide | Ca | 0.006093 | 0.004070 | -0.002023 | True |
| medium | extended | peptide | R | 0.016815 | 0.015418 | -0.001397 | True |
| medium | connected | metal | La | 0.115941 | 0.113810 | -0.002131 | True |
| medium | connected | metal | Ca | 0.050864 | 0.051942 | 0.001078 | True |
| medium | connected | metal | R | -0.065077 | -0.061868 | 0.003209 | True |
| medium | connected | peptide | La | -0.007299 | -0.006614 | 0.000686 | True |
| medium | connected | peptide | Ca | 0.009256 | 0.007783 | -0.001474 | True |
| medium | connected | peptide | R | 0.016556 | 0.014396 | -0.002160 | True |
| large | extended | metal | La | 0.112904 | 0.112398 | -0.000506 | True |
| large | extended | metal | Ca | 0.049819 | 0.050401 | 0.000583 | True |
| large | extended | metal | R | -0.063085 | -0.061996 | 0.001089 | True |
| large | extended | peptide | La | -0.010722 | -0.011250 | -0.000528 | True |
| large | extended | peptide | Ca | 0.006093 | 0.004443 | -0.001650 | True |
| large | extended | peptide | R | 0.016815 | 0.015692 | -0.001123 | True |
| large | connected | metal | La | 0.115941 | 0.113494 | -0.002447 | True |
| large | connected | metal | Ca | 0.050864 | 0.052322 | 0.001458 | True |
| large | connected | metal | R | -0.065077 | -0.061172 | 0.003905 | True |
| large | connected | peptide | La | -0.007299 | -0.006892 | 0.000407 | True |
| large | connected | peptide | Ca | 0.009256 | 0.007928 | -0.001328 | True |
| large | connected | peptide | R | 0.016556 | 0.014820 | -0.001735 | True |

Even energies are kcal/mol. Negative curvature is retained, not clamped.
Peptide-path curvature includes the curved coordinate and is not a Cartesian eigenvalue.

```json
{
  "medium": {
    "vacuum_gradient_pass": true,
    "curvature_pass": true,
    "anchored_prediction_pass": true
  },
  "large": {
    "vacuum_gradient_pass": true,
    "curvature_pass": true,
    "anchored_prediction_pass": true
  }
}
```

All projections, secant curvatures, anchored residuals and vacuum/solvent/short components are in result.json.
A passing direction does not validate a coupled scaffold model. response_model_not_validated.

## Interpretation and measured cost

Both checkpoint screens pass. Largest DFT-anchored displacement-energy error:
medium0.003766078kcal/mol; large0.004461813. Largest even-term error:
medium0.003208894; large0.003904629. These are two physical directions at one
small amplitude in two representations of ONE GGR source, not independent
proteins or validation of a full stiffness matrix. In particular, the absolute
0.005kcal tolerance floor permits some relative errors above25% for small
peptide terms. Both negative La peptide curvatures remain negative.

First derivatives differ more than curvatures: extended-model Ca metal derivative
is +1.633484kcal/mol/Angstrom for DFT, versus−0.945421 for medium MACE+GB.
Consequently this supports using DFT energy/gradient anchors; it does not support
unqualified minimization with MACE's own forces. The total MACE+GB gradient has
not been implemented: finite energy differences include charge response, whereas
fixed-charge GB forces do not. No environment-corrected DFT gradient is claimed.

All40MACE+40GB calls completed,zero failed inference attempts and zero newDFT.
Jobs1200731/1200732 took156/176seconds; GB1200733/1200734 took24/21seconds.
Total377allocatedGPU seconds,6032allocatedcore-seconds,473.063actualCPU seconds.
The jobs overlapped, so this sum is not turnaround time. Model inference itself
totalled26.613082/45.962509seconds for medium/large; most batch time is per-task
model/process startup. Peak GPU allocated495767552/947634688bytes; reserved
601882624/1103101952bytes. Preparation/engineering time was not fully profiled.
One attempted large-collection read occurred before its file existed; it failed
before preparation or inference and was retried after completion.

**Continue the structural-response route.** Next validate coupled physical
coordinates and transfer beyond GGR before enabling any relaxation correction.
The direct MACE affinity descriptor remains rejected on the declared development
ordering tests; baseline/default unchanged. A curvature-screen success does not
complete the active discriminator goal.
