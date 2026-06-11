# Nested (block-structured / MATNEST) matrices in PSBLAS

Author: Simone Staccone (Stack-1)

This directory contains the tests for the **nested matrix** support added to PSBLAS: a block-structured distributed operator

```
      [ A11  A12  ... ]
  M = [ A21  A22  ... ]
      [ ...       ... ]
```

whose blocks are kept as separate sparse matrices (one per field) but which presents itself to Krylov solvers and preconditioners as a **single ordinary distributed matrix**. It is the PSBLAS analogue of PETSc's `MATNEST`.

The motivating case is the **saddle-point** system

```
  M = [ A    B^T ]
      [ B     0  ]
```

(symmetric indefinite, with the (2,2) block absent), but the implementation supports any square multi-field block operator with possibly **rectangular**
sub-blocks.


## 1. Concepts

* **Field** — a contiguous index space (e.g. velocity `V` and pressure `Q` in a saddle-point problem). Each field has its own `psb_desc_type` distribution.
* **Block (i,j)** — the sub-matrix coupling field `i` (rows) with field `j` (columns). It may be rectangular (`|field i| /= |field j|`) and may be absent.
* **Global operator** — the blocks are concatenated into a single **square** operator `M` of size `sum(field_sizes)`, distributed over one **composed global descriptor** with a **union halo** (one halo exchange per matrix-vector product, covering all blocks of a given column field at once).
* **Rectangular blocks** — PSBLAS does not support rectangular *distributed* matrices, but it does support rectangular *local* CSR/COO matrices. The rectangular product therefore happens only in the **local** block `csmv`; the only object carrying a descriptor (and hence communication) is the global operator, which is always square.

The global operator (`a_glob`) and global descriptor (`desc_glob`) can be passed unchanged to `psb_spmm`, `psb_krylov`, and the standard preconditioners.


## 2. Recommended API: `psb_d_nest_matrix`

The easy way to build a nested matrix is the `psb_d_nest_matrix` type (module `psb_d_nest_builder_mod`, re-exported by the umbrella `psb_d_nest_mod`), which follows the usual PSBLAS `init` / `ins` / `asb` pattern and hides all the descriptor / halo / compose / setup boilerplate:

```fortran
use psb_d_nest_mod

type(psb_d_nest_matrix) :: nested_matrix
integer(psb_lpk_)       :: n1, n2

! 1) declare the field structure: 2 fields of global size n1, n2
call nested_matrix%init(ctxt, [n1, n2], info)

! 2) insert the block values, owned rows only (PSBLAS convention).
!    ins(block_row, block_col, n_entries, entry_rows, entry_cols, entry_vals, info)
!    rows are GLOBAL indices in field block_row, columns in field block_col.
call nested_matrix%ins(1, 1, nz_A,  iaA,  jaA,  valA,  info)   ! A   = block (1,1)
call nested_matrix%ins(1, 2, nz_Bt, iaBt, jaBt, valBt, info)   ! B^T = block (1,2)
call nested_matrix%ins(2, 1, nz_B,  iaB,  jaB,  valB,  info)   ! B   = block (2,1)
!    (the (2,2) block is simply not inserted)

! 3) assemble: builds nested_matrix%a_glob and nested_matrix%desc_glob
call nested_matrix%asb(info)

! 4) from here on it is an ordinary distributed matrix/descriptor
call psb_geall(x, nested_matrix%desc_glob, info)
...
call psb_krylov('CG', nested_matrix%a_glob, prec, b, x, eps, &
     & nested_matrix%desc_glob, info, itmax=..., iter=..., err=...)

! 5) release
call nested_matrix%free(info)
```

Notes:

* To know which rows it owns in a field, a process can query the per-field descriptor exposed as `nested_matrix%field_desc(i)` (e.g. `nested_matrix%field_desc(1)%get_local_rows()` and `%l2g(...)`), exactly as it would with a plain `psb_cdall` descriptor.
* Off-diagonal blocks may be rectangular: the cross-field column indices are registered into the union halo automatically by `ins`.
* The CG solver requires an SPD operator; a genuine saddle-point operator is indefinite and needs MINRES/GMRES (plus, eventually, a block preconditioner).
* **Do not copy/move** a `psb_d_nest_matrix` after `asb`: the wrapped operator holds internal pointers into the object.


## 3. Low-level path (advanced)

`psb_d_nest_matrix` is built on three lower-level pieces, available directly for advanced use (see `psb_d_nest_cg_test.F90` for an end-to-end example):

* `psb_cd_nest_compose(grid_desc, desc_glob, info)` — compose the per-field descriptors into the single global descriptor with the union halo.
* `psb_d_nest_base_setup(nest_op, block_storage, grid_desc, desc_glob, info)` — set up the `psb_d_nest_base_mat` operator (implements the local `csmv`).
* `psb_d_nest_rect_block(blk, nz, ia, ja, val, desc_row, desc_col, info)` — build a single (possibly rectangular) local block from global triplets, with rows localized against `desc_row` and columns against `desc_col`.

A field-split interface (`psb_d_nest_get_block`, `psb_d_nest_get_field_desc`,
`psb_d_nest_restrict_field`, `psb_d_nest_prolong_field`,
`psb_d_nest_apply_block`) is exposed on `psb_d_nest_base_mat` as the hook for a future block (field-split / Schur) preconditioner.


## 4. Tests

| Test                         | What it checks |
|------------------------------|----------------|
| `psb_d_nest_glob_test`       | Square 2×2 operator built with `psb_d_nest_matrix`; the nested `psb_spmm` is compared bit-for-bit against the same matrix assembled monolithically in CSR. |
| `psb_d_nest_rect_test`       | Same, with fields of different size (`|V| = 2|Q|`) and genuinely **rectangular** off-diagonal blocks. |
| `psb_d_nest_cg_test`         | Standard PSBLAS **CG** on an SPD, ill-conditioned operator (1D Laplacian reordered red-black), built on the **low-level path**; the solution is recovered to machine precision over hundreds of matvecs. |
| `psb_d_nest_builder_test`    | Same CG solve as above but built through the `psb_d_nest_matrix` utility (high-level path). |

All tests run both serially and in parallel, and the result is invariant with respect to the number of MPI processes.

### Build and run

The PSBLAS library must be built/installed first (from the repository root):

```sh
make            # or the CMake build
```

Then, from this directory:

```sh
make                                   # builds the executables into ./runs
./runs/psb_d_nest_glob_test            # serial
mpirun -np 4 ./runs/psb_d_nest_rect_test
mpirun -np 4 ./runs/psb_d_nest_cg_test
mpirun -np 4 ./runs/psb_d_nest_builder_test
```

Each test prints a single `[PASS]` / `[FAIL]` line (printed by rank 0).


## 5. Source files

Library (under `base/modules/`):

* `desc/psb_desc_nest_mod.f90`        — `psb_desc_nest_type` (grid of per-field descriptors)
* `serial/psb_d_nest_mat_mod.f90`     — `psb_d_nest_sparse_mat` (block storage)
* `serial/psb_d_nest_base_mat_mod.F90`— `psb_d_nest_base_mat` (the MATNEST operator + `csmv`)
* `tools/psb_cd_nest_tools_mod.F90`   — descriptor tools (`psb_cd_nest_compose`, ...)
* `tools/psb_d_nest_tools_mod.F90`    — block tools (`psb_d_nest_rect_block`, ...)
* `tools/psb_d_nest_builder_mod.F90`  — `psb_d_nest_matrix` frontend (init/ins/asb)
* `psb_d_nest_mod.f90`                — umbrella module (`use psb_d_nest_mod`)
