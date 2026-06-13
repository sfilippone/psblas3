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
* **Block (i,j)** — the sub-matrix coupling field `i` (rows) with field `j` (columns). It may be rectangular (different field sizes) and may be absent.
* **Global operator** — the blocks are concatenated into a single **square** operator `M` of size `sum(field_sizes)`, distributed over one **composed global descriptor** with a **union halo** (one halo exchange per matrix-vector product, covering all blocks of a given column field at once).
* **Rectangular blocks** — PSBLAS does not support rectangular *distributed* matrices, but it does support rectangular *local* CSR/COO matrices. The rectangular product therefore happens only in the **local** block `csmv`; the only object carrying a descriptor (and hence communication) is the global operator, which is always square.

The global operator (`a_glob`) and global descriptor (`desc_glob`) can be passed unchanged to `psb_spmm`, `psb_krylov`, and the standard preconditioners.


## 2. Quick start: `psb_d_nest_matrix`

The easy way to build a nested matrix is the `psb_d_nest_matrix` type (module `psb_d_nest_builder_mod`, re-exported by the umbrella `psb_d_nest_mod`), which follows the usual PSBLAS `init` / `ins` / `asb` pattern and hides all the descriptor / halo / compose / setup boilerplate:

```fortran
use psb_d_nest_mod

type(psb_d_nest_matrix) :: nested_matrix
integer(psb_lpk_)       :: n1, n2

! 1) declare the field structure: 2 fields of global size n1, n2
call nested_matrix%init(ctxt, [n1, n2], info)

! 2) insert the block values, owned rows only (PSBLAS convention).
call nested_matrix%ins(1, 1, nz_A,  iaA,  jaA,  valA,  info)   ! A   = block (1,1)
call nested_matrix%ins(1, 2, nz_Bt, iaBt, jaBt, valBt, info)   ! B^T = block (1,2)
call nested_matrix%ins(2, 1, nz_B,  iaB,  jaB,  valB,  info)   ! B   = block (2,1)
!    (the (2,2) block is simply not inserted)

! 3) assemble: builds nested_matrix%a_glob and nested_matrix%desc_glob
call nested_matrix%asb(info)

! 4) from here on it is an ordinary distributed matrix/descriptor
call psb_geall(x, nested_matrix%desc_glob, info)
...
call prec%init(ctxt, 'BJAC', info)
call prec%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
call psb_krylov('CG', nested_matrix%a_glob, prec, b, x, eps, &
     & nested_matrix%desc_glob, info, itmax=..., iter=..., err=...)

! 5) release
call nested_matrix%free(info)
```


## 3. User API reference

All of the public API is available through the umbrella module:

```fortran
use psb_d_nest_mod
```

### 3.1 `type(psb_d_nest_matrix)` — the nested matrix (recommended)

| Member | Meaning |
|--------|---------|
| `a_glob` | `type(psb_dspmat_type)` — the assembled global operator; pass it to `psb_spmm`, `psb_krylov`, `prec%build` |
| `desc_glob` | `type(psb_desc_type)` — the composed global descriptor; pass it wherever a descriptor is expected |
| `field_desc(i)` | `type(psb_desc_type)` — the descriptor of field `i` (advanced use; for the common queries see `get_owned_rows` below) |
| `n_fields` | number of fields |

To know which rows it must insert, a process asks the matrix directly — no
descriptor jargon needed:

```fortran
integer(psb_lpk_), allocatable :: my_rows(:)

my_rows = nested_matrix%get_owned_rows(1)     ! global rows of field 1 owned here
do k = 1, size(my_rows)
   global_row = my_rows(k)
   ...                                        ! build the entries of this row
end do
```

| Query | Result |
|-------|--------|
| `nested_matrix%get_owned_rows(i_field)` | `integer(psb_lpk_), allocatable (:)` — the GLOBAL indices (in the field index space, 1..field size) of the rows of field `i_field` owned by this process |
| `nested_matrix%get_owned_row_count(i_field)` | `integer(psb_ipk_)` — how many rows of field `i_field` this process owns |

Methods (collective over the communicator unless noted):

#### `call nested_matrix%init(ctxt, field_sizes, info)`

Create the field structure. One descriptor per field is created with a block
row distribution; the total size is independent of the number of processes.

| Argument | Type | Intent | Meaning |
|----------|------|--------|---------|
| `ctxt` | `type(psb_ctxt_type)` | in | parallel context from `psb_init` |
| `field_sizes(:)` | `integer(psb_lpk_)` | in | global size of each field, e.g. `[n1, n2]` |
| `info` | `integer(psb_ipk_)` | out | return code, `psb_success_` on success |

#### `call nested_matrix%ins(block_row, block_col, n_entries, entry_rows, entry_cols, entry_vals, info)`

Insert a batch of coefficients into block `(block_row, block_col)`. May be
called any number of times per block, in any order, before `asb`. Each process
inserts only the rows it owns (PSBLAS convention); cross-field columns are
registered into the union halo automatically.

| Argument | Type | Intent | Meaning |
|----------|------|--------|---------|
| `block_row` | `integer(psb_ipk_)` | in | row-field index of the block (1..n_fields) |
| `block_col` | `integer(psb_ipk_)` | in | column-field index of the block (1..n_fields) |
| `n_entries` | `integer(psb_ipk_)` | in | number of triplets in this batch |
| `entry_rows(:)` | `integer(psb_lpk_)` | in | GLOBAL row indices in field `block_row` (1..field size) |
| `entry_cols(:)` | `integer(psb_lpk_)` | in | GLOBAL column indices in field `block_col` (1..field size) |
| `entry_vals(:)` | `real(psb_dpk_)` | in | coefficient values |
| `info` | `integer(psb_ipk_)` | out | return code |

#### `call nested_matrix%asb(info [, type] [, mold])`

Assemble: builds the per-field halos, the (possibly rectangular) local blocks,
the composed global descriptor `desc_glob` and the global operator `a_glob`.
After `asb` no further `ins` is allowed, and the object must not be
copied/moved (the operator holds internal pointers into it).

The optional arguments select the **storage format of the blocks**:

| Argument | Type | Meaning |
|----------|------|---------|
| `type` | `character(len=*)` | a base format name: `'CSR'` (default), `'CSC'`, `'COO'` |
| `mold` | `class(psb_d_base_sparse_mat)` | any format class, e.g. `psb_d_ell_sparse_mat` / `psb_d_hll_sparse_mat` from `psb_ext` |

The nested operator is format-agnostic: every operation delegates to the
blocks' own methods, so each block runs its native kernels.

#### `call nested_matrix%free(info)`

Release every internal object (blocks, descriptors, global operator).

### 3.2 Solvers and preconditioners

`a_glob` / `desc_glob` work with the standard PSBLAS infrastructure:

* **Krylov methods** — `psb_krylov('CG' | 'BICGSTAB' | 'GMRES' | ..., nested_matrix%a_glob, prec, b, x, eps, nested_matrix%desc_glob, info, ...)`. Remember that CG requires an SPD operator; a genuine saddle-point operator is indefinite and needs MINRES/GMRES.
* **Preconditioners** — all the stock PSBLAS one-level preconditioners can be built directly on the nested operator:
  * `'NONE'` — identity;
  * `'DIAG'` / `'JACOBI'` — diagonal scaling (served by the nested `get_diag`, which concatenates the diagonals of the diagonal blocks; absent blocks contribute zeros);
  * `'BJAC'` — block Jacobi with ILU factorization of the local rows (served by the nested `csgetrow`, which extracts the local rows of the global operator across all blocks).

```fortran
call prec%init(ctxt, 'BJAC', info)
call prec%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
```

### 3.3 Implemented base-class contract

The nested operator (`psb_d_nest_base_mat`) implements the standard
`psb_d_base_sparse_mat` contract by delegation to the blocks, so it can be used
wherever an assembled PSBLAS matrix is expected:

* **Products** — `csmv` (also transposed, `trans='T'`), `csmm` (multi-RHS),
  `vect_mv` (encapsulated vectors: gathers/scatters through the vectors' own
  `gth`/`sct` and runs each block through its `vect_mv`, so device block
  formats execute their device kernels).
* **Access/conversions** — `get_diag`, `csgetrow` (and `csget`/`csgetblk`
  through the base generics), `cp_to_coo`/`mv_to_coo` (and `cscnv`, `csclip`,
  `tril`/`triu`, ... through the base generics built on the COO route).
* **Reductions** — `rowsum`/`arwsum`, `colsum`/`aclsum`, `maxval`,
  `spnmi` (infinity norm), `spnm1` (1-norm).
* **Mutation/bookkeeping** — `scal` (left/right) and `scals` (the operator is a
  view: scaling acts on the blocks), `clone` (shares the blocks, re-owns the
  private index maps), `mold`, `sizeof`, `free`, `get_nzeros`, `get_fmt`.

Intentionally **not** implemented (they fail with the standard "missing
override" error): `cp_from_coo`/`mv_from_coo` (a nested operator cannot be
built from a flat matrix without the field structure), `csput` (insertions go
to the blocks before assembly), `cssv`/`cssm` (a triangular solve is undefined
for a block operator).

### 3.4 Low-level API (advanced)

`psb_d_nest_matrix` is built on lower-level pieces, available directly:

* `psb_cd_nest_compose(grid_desc, desc_glob, info)` — compose the per-field descriptors into the single global descriptor with the union halo.
* `psb_d_nest_base_setup(nest_op, block_storage, grid_desc, desc_glob, info)` — set up the `psb_d_nest_base_mat` operator (implements the local `csmv`, `get_diag`, `csgetrow`).
* `psb_d_nest_rect_block(blk, nz, ia, ja, val, desc_row, desc_col, info)` — build a single (possibly rectangular) local block from global triplets, with rows localized against `desc_row` and columns against `desc_col`.

A field-split interface (`psb_d_nest_get_block`, `psb_d_nest_get_field_desc`, `psb_d_nest_restrict_field`, `psb_d_nest_prolong_field`, `psb_d_nest_apply_block`) is exposed on `psb_d_nest_base_mat` as the hook for a future block (field-split / Schur) preconditioner.


## 4. Tests

| Test                         | What it checks |
|------------------------------|----------------|
| `psb_d_nest_glob_test`       | Square 2×2 operator built with `psb_d_nest_matrix`; the nested `psb_spmm` is compared bit-for-bit against the same matrix assembled monolithically in CSR. |
| `psb_d_nest_rect_test`       | Same, with fields of different size (`nV = 2 nQ`) and genuinely **rectangular** off-diagonal blocks. |
| `psb_d_nest_cg_test`         | Standard PSBLAS **CG** on an SPD, ill-conditioned operator (1D Laplacian reordered red-black), solved under every stock preconditioner (`NONE`, `DIAG`, `BJAC`/ILU(0)); requires convergence to machine precision for all of them, and that `DIAG` reproduces the `NONE` iteration count exactly (a bit-precise check of the nested `get_diag`, since the diagonal is the constant `2I`). |

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
```

Each test prints a single `[PASS]` / `[FAIL]` line (printed by rank 0).


## 5. Source files

Library (under `base/modules/`):

* `desc/psb_desc_nest_mod.f90`        — `psb_desc_nest_type` (grid of per-field descriptors)
* `serial/psb_d_nest_mat_mod.f90`     — `psb_d_nest_sparse_mat` (block storage)
* `serial/psb_d_nest_base_mat_mod.F90`— `psb_d_nest_base_mat` (the MATNEST operator: `csmv`, `get_diag`, `csgetrow`)
* `tools/psb_cd_nest_tools_mod.F90`   — descriptor tools (`psb_cd_nest_compose`, ...)
* `tools/psb_d_nest_tools_mod.F90`    — block tools (`psb_d_nest_rect_block`, ...)
* `tools/psb_d_nest_builder_mod.F90`  — `psb_d_nest_matrix` frontend (init/ins/asb)
* `psb_d_nest_mod.f90`                — umbrella module (`use psb_d_nest_mod`)
