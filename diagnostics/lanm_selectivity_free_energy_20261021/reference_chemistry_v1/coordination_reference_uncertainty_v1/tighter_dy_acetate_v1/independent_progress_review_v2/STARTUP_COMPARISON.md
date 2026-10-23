# Exact-role startup comparison

The completed Tight control **18382** used four MPI ranks on one node (`node-64-768g-14`). Native total time was **1330.206 s**, the whole startup **2.594 s**, and the RI/C V-Matrix plus Cholesky operation printed **0.1 s**.

The captured current input differs only in four-to-eight MPI ranks and the three declared downstream local-correlation cutoffs. The entire printed SHARK startup-settings block is byte-identical: 8 centers, 189 orbital functions/71 shells, 718 auxiliary-C functions/176 shells, and 3000 MiB MaxCore. Geometry, basis assets, SCF settings, ghost treatment and effective state are unchanged. Current placement is four ranks on each of two allocated nodes.

The known workload does **not** explain more than 113 minutes at this pre-SCF operation. The active CPU snapshot could represent a busy MPI wait; it does not prove useful progress or identify a particular transport defect. Eight versus four ranks also differs, so placement is not uniquely established as the cause.

A separately preserved **single-node, eight-MPI resource-only recovery is justified** if root closes the anomalous attempt. Preserve the exact tighter input, default grant-derived MaxCore, original admission and CP threshold. This is a prospective recommendation only: no cancellation, new job, live polling or method change occurred. The incomplete CP cycle remains unavailable.
