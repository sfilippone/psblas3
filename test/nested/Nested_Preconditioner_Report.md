# Nested Preconditioner Report

## Executive Summary

This report documents the double-precision nested preconditioner functionality added to PSBLAS. The implementation provides a PSBLAS-native field-split preconditioner for matrices represented through the nested matrix infrastructure. It is implemented primarily in:

- `prec/psb_d_nestedprec.f90`
- `prec/impl/psb_d_prec_type_impl.f90`
- `base/modules/serial/psb_d_nest_base_mat_mod.F90`
- nested test drivers under `test/nested`

The nested preconditioner is intended for block-structured systems, especially multiphysics and field-split problems, where the global matrix is naturally written as:

```text
      [ A_11  A_12  ... A_1m ]
  A = [ A_21  A_22  ... A_2m ]
      [  ...   ...  ...  ... ]
      [ A_m1  A_m2  ... A_mm ]
```

Each field owns a PSBLAS descriptor and each block is a normal sparse PSBLAS matrix. The preconditioner extracts the field diagonal blocks, builds ordinary PSBLAS sub-preconditioners on those blocks, and combines the field solves through additive, multiplicative, symmetric multiplicative, or two-field Schur-style compositions.

The implementation now also includes independent per-field inner Krylov solve contexts. These are configured through the existing `prec%set` option mechanism and are applied inside field solves when requested. The current inner methods are local nested-preconditioner kernels for `CG` and `BICGSTAB`, rather than calls into the top-level `psb_krylov` driver.

```
Outer Krylov solver (psb_krylov)
    |
    +-- Preconditioner / Block solve
            |
            +-- Field 1 solve
            |      |
            |      +-- Inner CG/BICGSTAB (optional)
            |
            +-- Field 2 solve
                   |
                   +-- Inner CG/BICGSTAB (optional)
```
## PSBLAS Background

PSBLAS separates sparse linear algebra into a few recurring concepts:

- a sparse matrix object, such as `psb_dspmat_type`
- a communication descriptor, `psb_desc_type`, that describes distributed row ownership and halo data
- dense vector storage, either raw arrays or PSBLAS vector objects
- a preconditioner wrapper, `psb_dprec_type`, that owns a polymorphic implementation derived from `psb_d_base_prec_type`
- Krylov solvers, such as `psb_krylov`, that operate on a matrix, a descriptor, and an optional preconditioner

The usual preconditioner lifecycle is:

```fortran
call prec%init('...', info)
call prec%set('...', value, info)
call prec%build(a, desc_a, info)
call psb_krylov(method, a, prec, b, x, eps, desc_a, info, ...)
call prec%free(info)
```

The nested preconditioner follows this same lifecycle. From the outside it is just another `psb_dprec_type` implementation selected by:

```fortran
call prec%init('NEST', info)
```

The outer Krylov method does not need to know that the preconditioner is block structured. It applies `prec` through the standard PSBLAS preconditioner virtual interface.

## Motivation

Many application matrices are block systems, not just scalar sparse matrices. Examples include:

- coupled PDE systems
- saddle point systems
- fluid-structure and multiphysics discretizations
- systems reordered by physical field or variable type
- matrices where different fields require different smoothers or approximate solvers

A scalar preconditioner treats the matrix as a single undifferentiated sparse operator. That is often too weak or too inflexible for block systems. A field-split preconditioner can use the mathematical structure of the problem:

- each field can use a different sub-preconditioner
- off-diagonal coupling can be handled additively, multiplicatively, or through Schur approximations
- expensive exact block solves can be replaced by approximate preconditioned field iterations
- the outer solver can remain unchanged

The goal of this implementation is to provide that behavior inside PSBLAS without introducing a separate solver framework or requiring users to leave the existing `prec%init`, `prec%set`, `prec%build`, and `psb_krylov` flow.

## Nested Matrix Infrastructure

The nested preconditioner depends on the nested matrix support already added to the repository. A nested matrix stores a global operator as a block matrix over fields. The important pieces are:

- `psb_d_nest_matrix`: user-facing builder/helper for nested matrices
- `psb_d_nest_base_mat`: sparse matrix backend that adapts the nested object to the standard PSBLAS sparse matrix interface
- `a_glob`: global sparse wrapper used by ordinary PSBLAS solver calls
- `desc_glob`: global descriptor for the composed vector layout
- per-field descriptors for field-local vectors and block operations
- helper routines to restrict/prolong between global and field vector layouts

The nested preconditioner uses the following helper operations from `psb_d_nest_base_mat_mod`:

```fortran
psb_d_nest_get_n_fields
psb_d_nest_get_field_owned
psb_d_nest_get_block
psb_d_nest_get_field_desc
psb_d_nest_restrict_field
psb_d_nest_restrict_field_local
psb_d_nest_prolong_field
psb_d_nest_apply_block
```

These routines are what keep the preconditioner PSBLAS-native. The preconditioner does not invent its own distributed layout. It asks the nested matrix infrastructure how fields map into the global vector and uses the existing descriptors for communication.

## Design Goals

The main design goals were:

- expose nested preconditioning through the standard PSBLAS preconditioner API
- reuse existing PSBLAS preconditioner implementations for diagonal field blocks
- preserve parallel correctness by using PSBLAS descriptors and nested field helpers
- support practical field-split compositions first
- allow per-field configuration of block preconditioners
- allow per-field inner Krylov solves without creating a dependency cycle in the library
- keep the first implementation scoped to double precision

The implementation is deliberately conservative. It does not attempt to become a separate nonlinear or block-solver framework. It provides preconditioner application paths that can be used by existing PSBLAS Krylov solvers.

## Mathematical Model

Assume a nested matrix with `m` fields:

```text
      [ A_11  A_12  ... A_1m ]
  A = [ A_21  A_22  ... A_2m ]
      [  ...   ...  ...  ... ]
      [ A_m1  A_m2  ... A_mm ]
```

For field `i`, let `B_i` denote the field preconditioner built from the diagonal block `A_ii`. A global vector `x` is restricted into field components:

```text
x -> (x_1, x_2, ..., x_m)
```

The nested preconditioner applies field solves and then prolongs field corrections back into the global vector layout:

```text
z_i ~= B_i r_i
z   = prolong(z_1, z_2, ..., z_m)
```

The different compositions change how field residuals are formed and how off-diagonal blocks are used.

## Additive And Diagonal Composition

The additive composition is the block-Jacobi field split:

```text
z_i = B_i r_i,  i = 1, ..., m
```

All fields are solved independently from the original restricted residual. Off-diagonal blocks are ignored during the preconditioner application. This is the simplest and most parallel composition.

For a three-field system:

```text
      [ A_11  A_12  A_13 ]
  A = [ A_21  A_22  A_23 ]
      [ A_31  A_32  A_33 ]
```

and residual:

```text
    [ r_1 ]
r = [ r_2 ]
    [ r_3 ]
```

the additive preconditioner applies:

```text
z_1 = B_1 r_1
z_2 = B_2 r_2
z_3 = B_3 r_3
```

Graphically:

```text
r_1 --> B_1 --> z_1

r_2 --> B_2 --> z_2

r_3 --> B_3 --> z_3
```

There is no dependency between field solves. Field 1 does not wait for field 2 or field 3; field 2 does not use corrections from field 1 or field 3; field 3 is also independent. In matrix terms, additive composition approximates the inverse of the block diagonal part:

```text
        [ A_11  0     0    ]^{-1}
P^{-1} ~[ 0     A_22  0    ]
        [ 0     0     A_33 ]
```

but with `B_i` in place of exact `A_ii^{-1}`:

```text
        [ B_1  0    0   ]
P^{-1} =[ 0    B_2  0   ]
        [ 0    0    B_3 ]
```

This is why additive composition is robust as a first preconditioner: it is cheap, parallel, and easy to reason about. Its weakness is that it does not directly account for coupling terms like `A_12`, `A_21`, `A_23`, or `A_32`. The outer Krylov method must handle those couplings.

`DIAGONAL` and `ADDITIVE` are treated as equivalent names in the implementation.

## Multiplicative Composition

The multiplicative composition is a field Gauss-Seidel sweep. Fields are visited in order. For field `i`, the residual is corrected using already computed field updates:

```text
r_i^hat = r_i - sum_{j < i} A_ij z_j
z_i     = B_i r_i^hat
```

The implementation uses nested block application for the off-diagonal products and prolongs intermediate field corrections into global work storage as needed. This lets later fields see the effect of earlier field corrections.

For a three-field forward sweep, the fields are processed as `1, 2, 3`.

Field 1 has nothing before it:

```text
rhs_1 = r_1
z_1   = B_1 rhs_1
```

Field 2 can use the field-1 correction:

```text
rhs_2 = r_2 - A_21 z_1
z_2   = B_2 rhs_2
```

Field 3 can use both earlier corrections:

```text
rhs_3 = r_3 - A_31 z_1 - A_32 z_2
z_3   = B_3 rhs_3
```

Graphically:

```text
r_1 ----------------------> B_1 ---> z_1
                                      |
                                      v
r_2 - A_21 z_1 ----------> B_2 ---> z_2
                                      |
                                      v
r_3 - A_31 z_1 - A_32 z_2 -> B_3 ---> z_3
```

This resembles a lower triangular block solve. In exact block form, replacing `B_i` by `A_ii^{-1}` would correspond to applying an approximate inverse of:

```text
[ A_11  0     0    ]
[ A_21  A_22  0    ]
[ A_31  A_32  A_33 ]
```

The actual preconditioner does not assemble this triangular matrix. It computes the needed off-diagonal products with the existing nested blocks and solves each field with `psb_d_nested_field_solve`.

## Symmetric Multiplicative Composition

The symmetric multiplicative composition performs a forward sweep followed by a backward sweep:

```text
forward:   i = 1, ..., m
backward:  i = m, ..., 1
```

The same sweep kernel is reused in both directions. This gives a symmetric-style block Gauss-Seidel preconditioner when the underlying block choices are compatible with the matrix properties.

For a three-field backward sweep, the fields are processed as `3, 2, 1`.

Field 3 has nothing above it in this ordering:

```text
rhs_3 = r_3
z_tilde_3 = B_3 rhs_3
```

Field 2 can now use the already computed field-3 correction:

```text
rhs_2 = r_2 - A_23 z_tilde_3
z_tilde_2 = B_2 rhs_2
```

Field 1 can use both later-field corrections:

```text
rhs_1 = r_1 - A_12 z_tilde_2 - A_13 z_tilde_3
z_tilde_1 = B_1 rhs_1
```

Graphically:

```text
r_3 ------------------------------> B_3 ---> z_tilde_3
                                              |
                                              v
r_2 - A_23 z_tilde_3 ------------> B_2 ---> z_tilde_2
                                              |
                                              v
r_1 - A_12 z_tilde_2 - A_13 z_tilde_3 -> B_1 ---> z_tilde_1
```

This resembles an upper triangular block solve:

```text
[ A_11  A_12  A_13 ]
[ 0     A_22  A_23 ]
[ 0     0     A_33 ]
```

The symmetric multiplicative composition combines the information flow from both directions. The forward sweep accounts for lower-block couplings such as `A_21`, `A_31`, and `A_32`; the backward sweep accounts for upper-block couplings such as `A_12`, `A_13`, and `A_23`. This is more expensive than additive composition because fields are no longer independent, but it can be much stronger when coupling terms dominate convergence.

## Schur-Style Composition

Schur-complement preconditioning starts from the two-field coupled system:

```text
    [ A_11  A_12 ]
A = [ A_21  A_22 ]
```

and:

```text
    [ x_1 ]        [ r_1 ]
x = [ x_2 ],  r = [ r_2 ]
```

The block system `A x = r` expands to:

```text
A_11 x_1 + A_12 x_2 = r_1
A_21 x_1 + A_22 x_2 = r_2
```

To eliminate field 1, solve the first equation for `x_1`:

```text
A_11 x_1 = r_1 - A_12 x_2
x_1      = A_11^{-1} (r_1 - A_12 x_2)
```

Insert this expression into the second equation:

```text
A_21 A_11^{-1} (r_1 - A_12 x_2) + A_22 x_2 = r_2
```

Rearranging gives:

```text
(A_22 - A_21 A_11^{-1} A_12) x_2 = r_2 - A_21 A_11^{-1} r_1
```

The matrix:

```text
S = A_22 - A_21 A_11^{-1} A_12
```

is the Schur complement.

This is important because the original coupled system has been reduced to three conceptual steps:

1. solve a field-1 problem
2. solve a field-2 Schur problem
3. recover or correct field 1

That pattern is the basis of many modern multiphysics and saddle-point preconditioners.

The exact block LU factorization is:

```text
    [ I              0 ] [ A_11  0 ] [ I  A_11^{-1} A_12 ]
A = [ A_21 A_11^{-1} I ] [ 0     S ] [ 0  I              ]
```

The diagonal blocks in this factorization are `A_11` and `S`; the remaining factors are triangular. This immediately suggests a preconditioner based on applying approximate triangular solves and approximate inverses for the diagonal factors.

In exact arithmetic, the inverse action would involve:

```text
P^{-1} = U^{-1} D^{-1} L^{-1}

D^{-1} = [ A_11^{-1}  0      ]
         [ 0          S^{-1} ]
```

In practice, nobody wants to explicitly compute the exact Schur complement for large sparse problems. The problem is `A_11^{-1}`. Even if `A_11`, `A_12`, and `A_21` are sparse, `A_11^{-1}` is generally dense, and therefore:

```text
A_21 A_11^{-1} A_12
```

can be dense or nearly dense. Assembling `S` would destroy sparsity and can make storage and computation infeasible for large systems.

The nested preconditioner therefore uses approximations:

```text
B_1 ~= A_11^{-1}
B_2 ~= A_22^{-1}  or  B_2 ~= S^{-1}
```

where `B_1` and `B_2` are field solves provided by `psb_d_nested_field_solve`. A field solve may be a direct block preconditioner application, or it may use the independent per-field inner Krylov context if `INNER_SOLVE` is enabled.

The implementation provides three Schur-style compositions.

### Lower Schur Variant

The lower Schur variant computes:

```text
z_1 = B_1 r_1
z_2 = B_2 (r_2 - A_21 z_1)
```

Substituting the first equation into the second gives:

```text
z_2 = B_2 (r_2 - A_21 B_1 r_1)
```

This is the approximate forward solve associated with the lower triangular block factor.

```text
r_1 --> B_1 --> z_1
                  |
                  v
              A_21 z_1
                  |
                  v
r_2 - A_21 z_1 --> B_2 --> z_2
```

### Upper Schur Variant

The upper Schur variant reverses the ordering:

```text
z_2 = B_2 r_2
z_1 = B_1 (r_1 - A_12 z_2)
```

This corresponds to applying the upper triangular block factor:

```text
r_2 --> B_2 --> z_2
                  |
                  v
              A_12 z_2
                  |
                  v
r_1 - A_12 z_2 --> B_1 --> z_1
```

### Full Schur Variant

The full Schur variant applies the lower-style solve and then corrects field 1 with the upper coupling:

```text
z_1 = B_1 r_1
z_2 = B_2 (r_2 - A_21 z_1)
z_1 = z_1 - B_1 A_12 z_2
```

This approximates the full `U^{-1} D^{-1} L^{-1}` block-LU inverse action. It is generally stronger than purely additive block-Jacobi or a single directional block sweep because it accounts for both off-diagonal couplings.

### Matrix-Free Schur Action

The most important implementation detail is that the code can apply a Schur action without assembling `S`.

Given a field-2 vector `x`, the matrix-free Schur action computes:

```text
y = S x
```

through the sequence:

```text
t = A_12 x
u = B_1 t
v = A_21 u
y = A_22 x - v
```

or equivalently:

```text
y = A_22 x - A_21 B_1 A_12 x
```

Graphically:

```text
x
|
v
A_12
|
v
B_1
|
v
A_21
|
v
subtract from A_22 x
|
v
y
```

No Schur matrix is ever stored. The code only needs block products and field solves.

### Matrix-Free Schur Solves

When `SCHUR_SOLVE = MATRIX_FREE`, the preconditioner needs an approximation to a Schur inverse even though it only has routines that apply Schur products. The implementation uses two matrix-free strategies, depending on the selected composition.

For the two-field Schur compositions (`SCHUR_LOWER`, `SCHUR_UPPER`, and `SCHUR_FULL`), the matrix-free Schur solve uses a small preconditioned Richardson iteration:

```text
z_{k+1} = z_k + M^{-1} (b - S z_k)
```

where `M^{-1}` is the field-2 solve, implemented through `psb_d_nested_field_solve` and usually acting like `B_2`.

Each Richardson step uses:

1. the matrix-free Schur product `S z_k`
2. a Schur residual `b - S z_k`
3. the field-2 preconditioner or field-2 inner solve to correct `z_k`

For the three-field PDE-control composition (`SCHUR_PDE_CONTROL`), the Schur field is field 3 and the KKT Schur complement is negative for the standard distributed-control block signs. The implementation therefore solves the positive matrix-free equation

```text
T x_3 = -r_3,       T = -S
```

with a small conjugate-gradient iteration. The product `T x_3` is applied without assembling `T`; each CG iteration calls the PDE-control Schur action and uses `psb_gedot` with the field-3 descriptor for global reductions. The field-3 block preconditioner is still used as an optional left preconditioner, but when `A_33 = 0` it is only an identity fallback.

This gives a stronger Schur inverse approximation for KKT optimal-control systems where the diagonal field-3 block is zero and a Richardson update based on `A_33` is ineffective.

The implementation does not assemble an exact sparse Schur complement. Instead, it applies Schur actions using:

- diagonal field solves through `psb_d_nested_field_solve`
- off-diagonal products through `psb_d_nest_apply_block`
- global/field work vectors managed by the nested matrix helpers

Supported Schur compositions are:

- `SCHUR_LOWER`
- `SCHUR_UPPER`
- `SCHUR_FULL`
- `SCHUR_PDE_CONTROL` (aliases: `PDE_CONTROL_SCHUR`, `SCHUR_2PLUS1`, `SCHUR3`, `SCHUR_3FIELD`, `SCHUR_THREE_FIELD`)

The Schur solve mode is controlled by `SCHUR_SOLVE`:

- `A22`: use the field-2 block preconditioner as the Schur approximation
- `MATRIX_FREE`: use a matrix-free Schur solve; two-field Schur compositions use Richardson, while `SCHUR_PDE_CONTROL` uses CG on the positive field-3 Schur action
- `SELFP` or `SELF`: build a PETSc-style lumped Schur approximation from `A_21`, the diagonal of `A_11`, and `A_12`, then build a normal PSBLAS `BJAC` Schur preconditioner on that approximation

For two-field Schur compositions, `SELFP` approximates

```text
S_p ~= A_22 - A_21 diag(A_11)^{-1} A_12
```

The assembled `S_p` is not the exact Schur complement. It is a sparse local approximation intended for use as a preconditioner, analogous in spirit to PETSc field-split `selfp` with a lumped inverse for `A_11`. After assembly, the nested preconditioner attaches the existing block-preconditioner machinery to `S_p`: currently this is a `BJAC` block preconditioner, and routed sub-options such as `SUB_SOLVE=ILU` and `SUB_FILLIN` are replayed onto it. During Schur application, `SCHUR_SOLVE=SELFP` applies this built Schur preconditioner directly instead of running Richardson on the matrix-free Schur action.

The `MATRIX_FREE` mode is controlled by `SCHUR_MAXIT` and `SCHUR_TOL`. The `SELFP` mode is controlled primarily by the sub-preconditioner options routed to the Schur approximation, especially `SUB_SOLVE` and `SUB_FILLIN`.

## Implementation Structure

The core implementation lives in `prec/psb_d_nestedprec.f90`, which defines `psb_d_nestedprec` and the concrete preconditioner type:

```fortran
type :: psb_d_nested_block_prec
  character(len=16) :: ptype
  class(psb_d_base_prec_type), allocatable :: pc
end type psb_d_nested_block_prec

type :: psb_d_nested_krylov_context
  logical :: enabled
  character(len=16) :: method
  integer(psb_ipk_) :: itmax
  integer(psb_ipk_) :: itrace
  integer(psb_ipk_) :: istop
  real(psb_dpk_) :: tol
end type psb_d_nested_krylov_context

type, extends(psb_d_base_prec_type) :: psb_d_nested_prec_type
  integer(psb_ipk_) :: composition
  character(len=32) :: composition_name
  character(len=16) :: default_block_ptype
  integer(psb_ipk_) :: schur_solve
  character(len=32) :: schur_solve_name
  integer(psb_ipk_) :: schur_maxit
  real(psb_dpk_) :: schur_tol
  type(psb_d_nested_krylov_context) :: default_krylov
  integer(psb_ipk_) :: nfields
  type(psb_d_nest_base_mat), pointer :: nest_op
  type(psb_d_nested_block_prec), allocatable :: blocks(:)
  character(len=16), allocatable :: field_block_ptype(:)
  type(psb_d_nested_krylov_context), allocatable :: field_krylov(:)
  type(psb_d_nested_iopt), allocatable :: field_iopts(:)
  type(psb_d_nested_ropt), allocatable :: field_ropts(:)
  type(psb_d_nested_copt), allocatable :: field_copts(:)
contains
  procedure, pass(prec) :: d_apply_v => psb_d_nested_apply_vect
  procedure, pass(prec) :: d_apply   => psb_d_nested_apply
  procedure, pass(prec) :: precbld   => psb_d_nested_precbld
  procedure, pass(prec) :: precinit  => psb_d_nested_precinit
  procedure, pass(prec) :: precseti  => psb_d_nested_precseti
  procedure, pass(prec) :: precsetr  => psb_d_nested_precsetr
  procedure, pass(prec) :: precsetc  => psb_d_nested_precsetc
  procedure, pass(prec) :: precdescr => psb_d_nested_precdescr
  procedure, pass(prec) :: dump      => psb_d_nested_dump
  procedure, pass(prec) :: clone     => psb_d_nested_clone
  procedure, pass(prec) :: free      => psb_d_nested_precfree
  procedure, pass(prec) :: sizeof    => psb_d_nested_sizeof
  procedure, pass(prec) :: get_nzeros => psb_d_nested_get_nzeros
end type psb_d_nested_prec_type
```

The defaults set by `psb_d_nested_precinit` are:

- composition: `ADDITIVE`
- default block solve: `DIAG`
- Schur solve: `A22`
- Schur maximum iterations: `4`
- Schur tolerance: `0.0`
- default inner Krylov disabled
- default inner Krylov method: `CG`
- default inner Krylov maximum iterations: `20`
- default inner Krylov trace: `-1`
- default inner Krylov stopping rule: `2`
- default inner Krylov tolerance: `1.0e-6`

The preconditioner stores:

- the selected field-split composition
- the default block preconditioner type
- optional per-field block preconditioner types
- Schur solve configuration
- default and per-field inner Krylov configuration
- pending per-field sub-preconditioner options
- one built block preconditioner per field
- a pointer to the nested matrix backend

The preconditioner does not own the nested matrix. It stores a pointer to the nested operator exposed by the global sparse matrix wrapper.

## Build Phase

The build phase is implemented by `psb_d_nested_precbld`.

At build time, the preconditioner:

1. releases previously built field preconditioners while keeping user configuration
2. verifies that the input sparse matrix is backed by a nested base matrix
3. stores a pointer to the nested operator
4. queries the number of fields
5. allocates one `psb_d_nested_block_prec` entry per field
6. obtains each diagonal block `A_ii`
7. obtains each field descriptor
8. allocates the selected PSBLAS block preconditioner type for that field
9. replays routed per-field options onto the block preconditioner
10. builds the field block preconditioner on `A_ii`

The build path deliberately uses existing PSBLAS preconditioner implementations. The nested preconditioner is a composition layer over ordinary field preconditioners, not a replacement for ILU, diagonal, approximate inverse, or block-Jacobi preconditioners.

## Field Block Preconditioners

The nested block solve type is configured through `BLOCK_SOLVE` or `NEST_BLOCK_SOLVE`.

Supported block preconditioner type names currently include:

- `DIAG`
- `BJAC`
- `NONE`

The implementation allocates the corresponding concrete PSBLAS preconditioner for each diagonal field block. `idx` can be used on `prec%set` calls to override a particular field:

```fortran
call prec%set('BLOCK_SOLVE', 'BJAC', info, idx=1)
call prec%set('BLOCK_SOLVE', 'DIAG', info, idx=2)
```

Field-level sub-preconditioner tuning is also routed through `psb_d_prec_type_impl.f90`. For nested preconditioners, options such as `SUB_SOLVE`, `SUB_FILLIN`, `SUB_ILUTHRS`, `INV_FILLIN`, and `INV_THRESH` are stored on the nested preconditioner and replayed onto the selected field block preconditioners during build.

This matters because the field block preconditioners do not exist until `precbld`. The routed option storage lets users configure fields before the nested object has been built.

## Independent Per-Field Inner Krylov Contexts

The implementation includes independent per-field inner Krylov solve contexts:

```fortran
type(psb_d_nested_krylov_context) :: default_krylov
type(psb_d_nested_krylov_context), allocatable :: field_krylov(:)
```

The default context applies to all fields when `idx` is not supplied. Field-specific contexts are allocated and updated when `idx` is supplied.

Supported inner methods are:

- `CG`
- `BICGSTAB`
- `NONE`, which disables the inner Krylov path for that field

The inner Krylov options are:

- `INNER_SOLVE`, also accepted as `KRYLOV_SOLVE` or `FIELD_SOLVE`
- `INNER_MAXIT`, also accepted as `INNER_ITMAX`, `KRYLOV_MAXIT`, or `KRYLOV_ITMAX`
- `INNER_TOL`, also accepted as `KRYLOV_TOL`
- `INNER_ITRACE`, also accepted as `KRYLOV_ITRACE`
- `INNER_ISTOP`, also accepted as `KRYLOV_ISTOP`

Example:

```fortran
call prec%set('INNER_SOLVE', 'CG', info, idx=1)
call prec%set('INNER_TOL', 1.0d-8, info, idx=1)
call prec%set('INNER_MAXIT', 50, info, idx=1)

call prec%set('INNER_SOLVE', 'BICGSTAB', info, idx=2)
call prec%set('INNER_TOL', 1.0d-6, info, idx=2)
call prec%set('INNER_MAXIT', 30, info, idx=2)
```

When a field has inner Krylov enabled, `psb_d_nested_field_solve` dispatches to:

- `psb_d_nested_inner_cg`
- `psb_d_nested_inner_bicgstab`

Those routines apply the field diagonal matrix through `psb_d_nested_field_matvec` and use the field block preconditioner as the inner preconditioner. When inner Krylov is disabled, `psb_d_nested_field_solve` directly applies the field block preconditioner.

The implementation does not call the top-level `psb_krylov` routine for inner solves. That choice avoids a build dependency cycle: the Krylov solvers live above the preconditioner layer, while `prec` must compile without depending on `linsolve`. The result is a PSBLAS-local inner Krylov implementation that provides independent per-field solve contexts without conflicting with the existing framework.

## Apply Phase

The apply phase is implemented by:

- `psb_d_nested_apply_vect` for PSBLAS vector objects
- `psb_d_nested_apply` for raw double-precision arrays

The vector-object path copies into raw array buffers, calls the array path, and copies back. The array path handles:

- optional transpose rejection
- alpha/beta scaling
- work-vector management
- dispatch to the selected composition

The selected composition is then applied by one of:

- `psb_d_nested_add_apply`
- `psb_d_nested_sweep`
- `psb_d_nested_apply_schur`

All paths eventually solve field systems through `psb_d_nested_field_solve`, so per-field block-preconditioner settings and per-field inner Krylov settings are shared by additive, multiplicative, symmetric, and Schur-style compositions.

## Additive Apply Path

`psb_d_nested_add_apply` performs the following for each field:

1. restrict the global residual to the field vector
2. solve the field problem through `psb_d_nested_field_solve`
3. prolong the field correction into the global output vector

Because the additive path does not use off-diagonal blocks, all field solves are independent at the mathematical level.

## Multiplicative Apply Path

`psb_d_nested_sweep` implements both forward and backward field sweeps. For each field:

1. restrict the original residual to the current field
2. subtract the effect of already computed field corrections through off-diagonal block products
3. solve the corrected field residual
4. prolong the correction into global work storage

The same routine is used for:

- `MULTIPLICATIVE`: one forward sweep
- `SYMMETRIC_MULTIPLICATIVE`: a forward sweep followed by a backward sweep

The sweep implementation uses field descriptors and global work vectors so that off-diagonal products see the correct local and halo data.

## Schur Apply Path

`psb_d_nested_apply_schur` is defined for exactly two fields. It implements lower, upper, and full Schur-style block factorizations using field solves and block products.

Supporting routines are:

- `psb_d_nested_schur_action`
- `psb_d_nested_schur_solve`

`psb_d_nested_schur_action` applies the matrix-free Schur operator:

```text
S x_2 = A_22 x_2 - A_21 B_1 A_12 x_2
```

where `B_1` is the field-1 block solve or inner field solve. `psb_d_nested_schur_solve` either applies the field-2 block solve directly (`A22`) or runs a short preconditioned Richardson iteration using the matrix-free Schur action (`MATRIX_FREE`).

## PDE-Control Schur Apply Path

`psb_d_nested_apply_pde_control_schur` is defined for exactly three fields. It targets KKT systems with the block structure

```text
[ A_11   0    A_13 ]
[  0    A_22  A_23 ]
[ A_31  A_32  A_33 ]
```

Fields 1 and 2 are grouped into the leading diagonal block `D = diag(A_11,A_22)`, and field 3 is the Schur field. The preconditioner applies the usual block factorization steps:

1. solve the leading fields with `B_1` and `B_2`
2. update the field-3 residual with `A_31 x_1 + A_32 x_2`
3. approximately solve the field-3 Schur equation
4. correct fields 1 and 2 with the upper coupling through `A_13` and `A_23`

The matrix-free positive PDE-control Schur action is

```text
T x_3 = A_31 B_1 A_13 x_3 + A_32 B_2 A_23 x_3 - A_33 x_3
```

For a standard distributed-control KKT matrix such as

```text
[ M      0      K  ]
[ 0    alpha M -M  ]
[ K     -M      0  ]
```

this operator is `T = -S`, where `S` is the field-3 KKT Schur complement. `psb_d_nested_pde_control_schur_solve` solves `T x_3 = -r_3` with a matrix-free CG iteration controlled by `SCHUR_MAXIT` and `SCHUR_TOL`. This is important when `A_33 = 0`: the old field-3 diagonal-block approximation cannot precondition the Schur complement, but the matrix-free CG path can still use the actual Schur action.

The `SCHUR_PDE_CONTROL_DIAG` composition applies the block-diagonal action

```text
z_1 = B_1 r_1,   z_2 = B_2 r_2,   z_3 = B_S r_3.
```

It reuses the PDE-control Schur solver for field 3 but reverses the KKT Schur
sign so that the third block approximates the positive operator `T^{-1}`.
Unlike `SCHUR_PDE_CONTROL`, it performs no lower or upper coupling corrections.
The composition is therefore structurally symmetric.  MINRES also requires
each selected block action to be fixed and SPD; in particular, the
residual-dependent inner CG used by `SCHUR_SOLVE=MATRIX_FREE` is not a fixed
linear preconditioner.  A fixed SPD field-3 solve such as the `A33`/identity
fallback is MINRES-compatible, while a future explicit SPD Schur
preconditioner can be attached through the same block-diagonal composition.

## Public Options

The nested preconditioner exposes integer option keys in `psb_d_nestedprec`:

```fortran
integer(psb_ipk_), parameter :: psb_d_nested_composition_ = 9101
integer(psb_ipk_), parameter :: psb_d_nested_block_solve_ = 9102
integer(psb_ipk_), parameter :: psb_d_nested_schur_solve_ = 9103
integer(psb_ipk_), parameter :: psb_d_nested_schur_maxit_ = 9104
integer(psb_ipk_), parameter :: psb_d_nested_schur_tol_ = 9105
integer(psb_ipk_), parameter :: psb_d_nested_inner_solve_ = 9106
integer(psb_ipk_), parameter :: psb_d_nested_inner_maxit_ = 9107
integer(psb_ipk_), parameter :: psb_d_nested_inner_tol_ = 9108
integer(psb_ipk_), parameter :: psb_d_nested_inner_itrace_ = 9109
integer(psb_ipk_), parameter :: psb_d_nested_inner_istop_ = 9110
```

The user-facing string options routed by `psb_d_prec_type_impl.f90` are:

| Option | Type | Values | Scope |
| --- | --- | --- | --- |
| `COMPOSITION`, `NEST_COMPOSITION` | character | `ADDITIVE`, `DIAGONAL`, `MULTIPLICATIVE`, `SYMMETRIC_MULTIPLICATIVE`, `SCHUR_LOWER`, `SCHUR_UPPER`, `SCHUR_FULL`, `SCHUR_PDE_CONTROL`, `SCHUR_PDE_CONTROL_DIAG` | global |
| `BLOCK_SOLVE`, `NEST_BLOCK_SOLVE` | character | `DIAG`, `BJAC`, `NONE` | global or `idx` field |
| `SCHUR_SOLVE`, `NEST_SCHUR_SOLVE` | character | `A22`, `A33`, `MATRIX_FREE`, `SELFP`, `SELF` | global |
| `SCHUR_MAXIT`, `NEST_SCHUR_MAXIT` | integer | nonnegative iteration count | global |
| `SCHUR_TOL`, `NEST_SCHUR_TOL` | real | nonnegative tolerance | global |
| `INNER_SOLVE`, `KRYLOV_SOLVE`, `FIELD_SOLVE` | character | `NONE`, `CG`, `BICGSTAB` | global or `idx` field |
| `INNER_MAXIT`, `INNER_ITMAX`, `KRYLOV_MAXIT`, `KRYLOV_ITMAX` | integer | positive iteration count | global or `idx` field |
| `INNER_TOL`, `KRYLOV_TOL` | real | nonnegative tolerance | global or `idx` field |
| `INNER_ITRACE`, `KRYLOV_ITRACE` | integer | trace level | global or `idx` field |
| `INNER_ISTOP`, `KRYLOV_ISTOP` | integer | stopping rule code | global or `idx` field |
| `SUB_SOLVE` | character | routed field sub-preconditioner selection | global or `idx` field |
| `SUB_FILLIN` | integer | routed ILU fill setting | global or `idx` field |
| `SUB_ILUTHRS` | real | routed ILUT threshold | global or `idx` field |
| `INV_FILLIN` | integer | routed inverse fill setting | global or `idx` field |
| `INV_THRESH` | real | routed inverse threshold | global or `idx` field |

## Example Usage

A typical nested preconditioner setup looks like:

```fortran
call prec%init('NEST', info)
call prec%set('COMPOSITION', 'SCHUR_FULL', info)
call prec%set('SCHUR_SOLVE', 'MATRIX_FREE', info)
call prec%set('SCHUR_MAXIT', 6, info)
call prec%set('SCHUR_TOL', 1.0d-8, info)

! Alternative two-field Stokes-style Schur approximation:
! call prec%set('SCHUR_SOLVE', 'SELFP', info)
! call prec%set('BLOCK_SOLVE', 'BJAC', info)
! call prec%set('SUB_SOLVE', 'ILU', info)
! call prec%set('SUB_FILLIN', 0, info)

call prec%set('BLOCK_SOLVE', 'BJAC', info, idx=1)
call prec%set('BLOCK_SOLVE', 'DIAG', info, idx=2)

call prec%set('INNER_SOLVE', 'CG', info, idx=1)
call prec%set('INNER_MAXIT', 40, info, idx=1)
call prec%set('INNER_TOL', 1.0d-8, info, idx=1)

call prec%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
call psb_krylov('BICGSTAB', nested_matrix%a_glob, prec, b, x, eps, &
     & nested_matrix%desc_glob, info)
```

This keeps the outer solve in the existing PSBLAS API while allowing field-specific behavior inside the nested preconditioner.

For a three-field optimal-control KKT system, the PDE-control Schur composition is configured as:

```fortran
call prec%init('NEST', info)
call prec%set('COMPOSITION', 'SCHUR_PDE_CONTROL', info)
call prec%set('SCHUR_SOLVE', 'MATRIX_FREE', info)
call prec%set('SCHUR_MAXIT', 200, info)
call prec%set('SCHUR_TOL', 0.0d0, info)

call prec%set('BLOCK_SOLVE', 'DIAG', info, idx=1)
call prec%set('BLOCK_SOLVE', 'DIAG', info, idx=2)
call prec%set('BLOCK_SOLVE', 'DIAG', info, idx=3)

call prec%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
call psb_krylov('BICGSTAB', nested_matrix%a_glob, prec, b, x, eps, &
     & nested_matrix%desc_glob, info)
```

The field-3 `DIAG` setting is harmless for zero `A_33` blocks because PSBLAS diagonal preconditioning treats zero diagonal entries as identity. The useful inverse approximation comes from the matrix-free CG solve on the PDE-control Schur action.

## Factory And API Integration

`psb_dprecinit` recognizes the nested preconditioner name:

```fortran
call prec%init('NEST', info)
```

and allocates:

```fortran
psb_d_nested_prec_type
```

The generic `prec%set` wrappers in `prec/impl/psb_d_prec_type_impl.f90` route nested-specific options only when the underlying concrete preconditioner is `psb_d_nested_prec_type`. For non-nested preconditioners, the existing behavior is preserved where appropriate.

This means nested-specific settings do not change the behavior of ordinary PSBLAS preconditioners.

## Parallel Layout Considerations

Parallel correctness depends on using PSBLAS descriptors and nested helper routines consistently.

Each field has:

- its own descriptor
- local owned rows
- possible halo entries needed by block products
- mappings between field-local vectors and global nested vectors

The preconditioner uses:

- `psb_d_nest_restrict_field` to extract field data from global vectors
- `psb_d_nest_restrict_field_local` where only local field data are needed
- `psb_d_nest_prolong_field` to insert field corrections into global vectors
- `psb_d_nest_apply_block` for diagonal and off-diagonal block products

Multiplicative and Schur compositions require more care than additive composition because they use off-diagonal products. Intermediate field corrections are prolonged into global work storage before being consumed by later block products. This ensures that distributed halo data are consistent with the nested matrix layout.

## Validation

The repository contains nested preconditioner tests under:

```text
test/nested
```

The main nested preconditioner drivers are:

- `test/nested/psb_d_nest_prec_cg_test.F90`
- `test/nested/psb_d_nest_mult_prec_cg_test.F90`
- `test/nested/psb_d_nest_schur_prec_cg_test.F90`

These tests exercise:

- nested matrix construction
- nested preconditioner initialization
- additive/diagonal composition
- multiplicative and symmetric multiplicative composition
- Schur lower, upper, full, and PDE-control composition
- matrix-free Schur solve mode, including CG for `SCHUR_PDE_CONTROL`
- per-field inner Krylov configuration and apply through the Schur preconditioner
- serial and MPI execution paths in the nested examples

Validation of the nested preconditioner examples showed convergence for the additive, multiplicative, and Schur-style configurations on the nested Laplacian test problems. The AMG4PSBLAS nested optimal-control MatrixMarket sample also converges with `SCHUR_PDE_CONTROL`, `SCHUR_SOLVE=MATRIX_FREE`, and `SCHUR_MAXIT=200` on serial and two-rank MPI runs.

The PETSc binary Stokes nested test was also used to compare the original matrix-free Schur solve against the new `SELFP` approximation on two MPI ranks. With `SCHUR_FULL`, `BICGSTAB`, `BLOCK_SOLVE=BJAC`, and `SUB_SOLVE=ILU`, the observed runs were:

| Schur solve | Sub fill | Iterations | Relative residual | Solve time |
| --- | ---: | ---: | ---: | ---: |
| `MATRIX_FREE` | 1 | 532 | `8.9510E-05` | about `4.95s` |
| `SELFP` | 0 | 412 | `9.5632E-05` | about `0.50s` |
| `SELFP` | 1 | 393 | `9.5120E-05` | about `0.51s` |
| `SELFP` | 2 | 570 | `9.5620E-05` | about `0.80s` |

For this case, `SELFP` reduced the iteration count compared with the matrix-free Schur solve and made each solve substantially cheaper. `SUB_FILLIN=1` gave the lowest iteration count, while `SUB_FILLIN=0` was the best practical default in the measured runs because it was close in iterations and cheapest to build/apply. Larger fill was not automatically better for the lumped Schur approximation.

The Schur nested preconditioner test includes an additional `SCHUR_FULL` / `A22+INNER` case. That case configures independent per-field inner Krylov contexts, builds the nested preconditioner, applies it directly to the right-hand side, and verifies a successful finite preconditioner application. The existing Schur cases continue to check outer BiCGSTAB convergence.

The serial and two-rank runs both passed after rebuilding the nested preconditioner object, the factory object, the preconditioner archive, and the test executable.

## Current Limitations

The current implementation has the following limitations:

- double precision only
- nested preconditioner name currently exposed as `NEST`
- transposed preconditioner application is rejected
- two-field Schur compositions are implemented by `SCHUR_LOWER`, `SCHUR_UPPER`, and `SCHUR_FULL`; the specialized three-field Schur composition is `SCHUR_PDE_CONTROL` only
- two-field matrix-free Schur solves use Richardson; PDE-control matrix-free Schur solves use CG on the positive field-3 Schur action
- `SELFP` is currently a two-field Schur approximation for `SCHUR_LOWER`, `SCHUR_UPPER`, and `SCHUR_FULL`; it assembles a local lumped Schur preconditioning matrix, not an exact distributed Schur complement
- inner Krylov methods currently include `CG` and `BICGSTAB`

## Future Work

Recommended next steps are:

- add single, complex, and complex double precision variants
- expose a dedicated Schur preconditioner for the PDE-control CG solve, instead of relying on the field-3 block preconditioner fallback

## Summary

The nested preconditioner adds a field-split preconditioning layer to PSBLAS while preserving the existing solver and preconditioner API. It reuses standard PSBLAS block preconditioners on diagonal field blocks and combines them through additive, multiplicative, symmetric multiplicative, two-field Schur, and three-field PDE-control Schur compositions.

The implementation integrates with:

- nested matrix storage and field mapping
- the `psb_dprec_type` factory
- the generic `prec%set` configuration path
- the ordinary PSBLAS preconditioner build/apply lifecycle
- the existing outer `psb_krylov` solver interface

The newer per-field inner Krylov support gives each field an independent local solve context, configured through the same PSBLAS preconditioner option API. It does not conflict with the existing framework because it is implemented inside the preconditioner layer and avoids adding a dependency from `prec` back to the higher-level Krylov solver layer.
