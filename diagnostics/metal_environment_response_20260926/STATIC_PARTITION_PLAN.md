# Retrospective static paired-contrast partition diagnostic

No new molecular calls. Existing small/expanded A/B outputs already inspected.
Compute D_X=[E_Ca(X)-E_La(X)]_expanded-[E_Ca(X)-E_La(X)]_small for X=A,B.
This subtracts metal-paired contrasts, not unbalanced absolute region energies.
NeutralThr added identically to bothmetals; matched Hamiltonian and fixedphysical
source, commonnetcharges. Metal-independent internal/environment-only terms cancel
within eachcontrast if identical; missing metal-dependent LJ and boundary electronic
representation do not. Thus D is electronic-component partition sensitivity only,
not fullhybridscore, affinity or calibrationtransfer. Unknownexpanded numerical
uncertainty and originalfailedrigidgate remain. No labels/thresholds fit.
This asks a different question from the response(B-A)partitionchange: stablelocal
response does not automatically imply stable staticCa-minusLa contrast.
