Communication scheme tests
==========================

This directory contains tests created after adding multiple communication schemes
to PSBLAS.

The goal is to exercise and compare communication patterns at different layers:

- direct halo exchange (`psi_swapdata`)
- overlap exchange for transpose/SpMV workflows (`psi_swaptran` + `psi_swapdata`)
- full Krylov solver runs (`CG`) using different comm schemes.

Communication schemes covered in this area:

- `psb_comm_isend_irecv_` (baseline point-to-point)
- `psb_comm_ineighbor_alltoallv_` (neighbor collective)
- `psb_comm_persistent_ineighbor_alltoallv_` (persistent neighbor collective)

See:

- `swapdata/` for a direct halo-exchange test.
- `spmv/` for an overlap SpMV test that uses different communication schemes.
- `cg/` for a conjugate-gradient solve-time comparison across the three schemes.
