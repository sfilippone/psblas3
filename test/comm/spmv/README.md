spmv overlap communication test
===============================

This test was added after introducing different communication schemes in PSBLAS.

It exercises the overlapped SpMV communication path inside `psb_spmm`.

Communication pattern:

- split exchange/computation flow (`start` + local compute + `wait`)
- halo/overlap update through internal swap routines used by SpMV kernels
- same matrix/vector workload repeated across schemes for timing comparison

Communication schemes compared:

- `psb_comm_isend_irecv_`
- `psb_comm_ineighbor_alltoallv_`
- `psb_comm_persistent_ineighbor_alltoallv_`

Unlike `swapdata/`, which checks direct halo exchange, this test covers the
overlapped SpMV workflow.

Run options
-----------

- Default PDE-generated matrix: `./runs/psb_spmv_kernel [--gpu=TRUE|FALSE] [--nooverlap]`
- External matrix: `./runs/psb_spmv_kernel [--gpu=TRUE|FALSE] --matrix=<path> [--fmt=MM|HB] [--nooverlap]`

The overlap path is enabled by default; pass `--nooverlap` to force the non-overlapped halo-update path.

When `--matrix` is provided, the benchmark reads and distributes that matrix instead of generating the 3D PDE test matrix.
