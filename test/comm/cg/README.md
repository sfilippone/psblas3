CG no-preconditioner communication test
=======================================

This test lives under `test/comm/cg` and builds a local executable:

- source: `psb_comm_cg_test.F90`
- executable: `runs/psb_comm_cg_test`

Behavior:

- generates a 3D PDE matrix using local `psb_d_gen_pde3d` (double precision)
- solves with `CG`
- uses preconditioner `NONE`
- runs CG three times, changing communication scheme of `x` each run

Communication pattern used in this test:

1. reset solution vector
2. set communication scheme on `x%v%comm_handle`
3. run full `psb_krylov('CG', ...)`
4. collect and compare solve time

Schemes compared:

- `psb_comm_isend_irecv_`
- `psb_comm_ineighbor_alltoallv_`
- `psb_comm_persistent_ineighbor_alltoallv_`

How to run
----------

From this directory:

- `make run` (defaults: `NP=4`, `IDIM=40`)
- `make run NP=8 IDIM=80`

Default PDE-generated matrix:

- `./runs/psb_comm_cg_test [idim] [nrep] [nwarm] [itmax] [--gpu=TRUE|FALSE]`

External matrix input:

- `./runs/psb_comm_cg_test [idim] [nrep] [nwarm] [itmax] --matrix=<path> [--fmt=MM|HB] [--gpu=TRUE|FALSE]`

When `--matrix` is provided, the test reads and distributes that matrix and ignores the PDE generator.
