swapdata communication test
============================

This test was added after introducing different communication schemes in PSBLAS.

It focuses on direct halo exchange through the `swapdata` path:

- index list type: halo (`psb_comm_halo_`)
- exchange API: `psi_swapdata`
- phases: `start`, `wait`, and `sync` (depending on test section)

Communication patterns exercised:

- baseline point-to-point (`isend/irecv`)
- neighbor collective (`ineighbor_alltoallv`)
- persistent neighbor collective (`persistent_ineighbor_alltoallv`)

This test validates the low-level communication behavior in isolation, without
the full SpMV overlap pipeline.

Run options
-----------

- Synthetic 3D stencil descriptor (default):
  - `./runs/psb_comm_test [--dim N] [--iters N] [--mode both|baseline|neighbor|persistent]`
- External matrix descriptor:
  - `./runs/psb_comm_test [--iters N] [--mode both|baseline|neighbor|persistent] --matrix <path> [--fmt MM|HB]`

With `--matrix`, the descriptor is built from the distributed input matrix pattern.
