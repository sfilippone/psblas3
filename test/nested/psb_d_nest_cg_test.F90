!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific prior written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
! File: psb_d_nest_cg_test.F90
!
! Program: psb_d_nest_cg_test
! Author:  Simone Staccone (Stack-1)
!
! Solves a linear system with the nested operator using the standard PSBLAS CG
! (psb_krylov('CG', ...)) under every stock one-level preconditioner, to show
! that the nested operator plugs into the PSBLAS preconditioning infrastructure:
!   NONE  (operator only),
!   DIAG  (exercises the nested get_diag),
!   BJAC  (ILU(0), exercises the nested csgetrow through the ILU build).
!
! CG needs a SYMMETRIC POSITIVE DEFINITE operator and, to stress the test
! (hundreds of matvecs), an ILL-CONDITIONED one.  We use a real case: the 1D
! Laplacian tridiag(-1, 2, -1) on m = 2*field_size nodes, REORDERED red-black
! (odd nodes -> field 1, even nodes -> field 2).  Under this reordering the
! Laplacian becomes exactly
!
!   M = [ 2I    C  ]   C(r,r) = -1 ,  C(r,r-1) = -1   (the Laplacian edges)
!       [ C^T   2I ]   C^T = exact transpose
!
! (odd nodes are not adjacent to each other -> diagonal blocks = 2I; every -1
! edge of the Laplacian becomes the coupling C).  M is therefore the 1D
! Laplacian up to a permutation: SPD but with lambda_min ~ (pi/m)^2 => cond ~
! N^2 => CG performs O(N) iterations that GROW with N.
!
! The operator is built with the psb_d_nest_matrix utility.  The test passes if
! every solve converges to the exact solution and DIAG reproduces the NONE
! iteration count exactly (with the constant diagonal 2I, Jacobi is a pure
! rescaling, so any mismatch would expose a wrong nested get_diag).
!
! Run: ./psb_d_nest_cg_test ; mpirun -np 4 ./psb_d_nest_cg_test
!
program psb_d_nest_cg_test
  use psb_base_mod
  use psb_util_mod
  use psb_prec_mod
  use psb_linsolve_mod
  use psb_d_nest_mod          ! umbrella: includes psb_d_nest_matrix (builder)
  implicit none

  type(psb_ctxt_type)        :: context
  integer(psb_ipk_)          :: my_rank, num_procs, info, i_local_row, entry_idx
  integer(psb_ipk_)          :: field1_local_rows, field2_local_rows
  integer(psb_lpk_)          :: field1_global_row, field2_global_row, field_size

  type(psb_d_nest_matrix)    :: nested_matrix
  type(psb_dprec_type)       :: preconditioner
  type(psb_d_vect_type)      :: x_solution, rhs, x_exact

  integer(psb_lpk_), allocatable :: entry_rows(:), entry_cols(:)
  integer(psb_lpk_), allocatable :: field1_rows(:), field2_rows(:)
  real(psb_dpk_),    allocatable :: entry_vals(:)

  ! solver parameters
  real(psb_dpk_)             :: diag_value, stop_tol, final_residual, norm_x_exact, solution_error
  integer(psb_ipk_)          :: max_iter, trace_level, n_iter, stop_criterion
  real(psb_dpk_), parameter  :: solution_tol = 1.0e-6_psb_dpk_

  ! stock preconditioners to exercise on the nested operator
  integer(psb_ipk_), parameter :: n_precs = 3
  character(len=6),  parameter :: prec_names(n_precs) = ['NONE  ', 'DIAG  ', 'BJAC  ']
  integer(psb_ipk_)            :: i_prec, iter_none, iter_diag
  logical                      :: all_passed

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  field_size  = 512                         ! global rows per field (global N = 2*field_size)
  diag_value  = 2.0_psb_dpk_                ! Laplacian diagonal (diagonal blocks = diag*I)
  stop_tol    = 1.0e-9_psb_dpk_
  max_iter    = 4000
  trace_level = 0
  stop_criterion = 2                        ! stop on the relative residual

  !---------------------------------------------------------------
  ! 1) create the nested operator: 2 fields of global size field_size
  !---------------------------------------------------------------
  call nested_matrix%init(context, [field_size, field_size], info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: nested_matrix%init info=', info; goto 9999
  end if
  field1_rows = nested_matrix%get_owned_rows(1)
  field2_rows = nested_matrix%get_owned_rows(2)
  field1_local_rows = size(field1_rows)
  field2_local_rows = size(field2_rows)

  !---------------------------------------------------------------
  ! 2) insert the blocks (owned rows only)
  !---------------------------------------------------------------
  ! block (1,1) = diag*I
  allocate(entry_rows(field1_local_rows), entry_cols(field1_local_rows), entry_vals(field1_local_rows))
  do i_local_row = 1, field1_local_rows
    field1_global_row = field1_rows(i_local_row)
    entry_rows(i_local_row) = field1_global_row
    entry_cols(i_local_row) = field1_global_row
    entry_vals(i_local_row) = diag_value
  end do
  call nested_matrix%ins(1, 1, field1_local_rows, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! block (2,2) = diag*I
  allocate(entry_rows(field2_local_rows), entry_cols(field2_local_rows), entry_vals(field2_local_rows))
  do i_local_row = 1, field2_local_rows
    field2_global_row = field2_rows(i_local_row)
    entry_rows(i_local_row) = field2_global_row
    entry_cols(i_local_row) = field2_global_row
    entry_vals(i_local_row) = diag_value
  end do
  call nested_matrix%ins(2, 2, field2_local_rows, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! block (1,2) = C : rows field1, cols field2 ; C(r,r)=-1, C(r,r-1)=-1
  allocate(entry_rows(2*field1_local_rows), entry_cols(2*field1_local_rows), entry_vals(2*field1_local_rows))
  entry_idx = 0
  do i_local_row = 1, field1_local_rows
    field1_global_row = field1_rows(i_local_row)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = field1_global_row
    entry_cols(entry_idx) = field1_global_row
    entry_vals(entry_idx) = -1.0_psb_dpk_
    if (field1_global_row > 1) then
      entry_idx = entry_idx + 1
      entry_rows(entry_idx) = field1_global_row
      entry_cols(entry_idx) = field1_global_row - 1_psb_lpk_
      entry_vals(entry_idx) = -1.0_psb_dpk_
    end if
  end do
  call nested_matrix%ins(1, 2, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! block (2,1) = C^T : rows field2, cols field1 ; C^T(s,s)=-1, C^T(s,s+1)=-1
  allocate(entry_rows(2*field2_local_rows), entry_cols(2*field2_local_rows), entry_vals(2*field2_local_rows))
  entry_idx = 0
  do i_local_row = 1, field2_local_rows
    field2_global_row = field2_rows(i_local_row)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = field2_global_row
    entry_cols(entry_idx) = field2_global_row
    entry_vals(entry_idx) = -1.0_psb_dpk_
    if (field2_global_row < field_size) then
      entry_idx = entry_idx + 1
      entry_rows(entry_idx) = field2_global_row
      entry_cols(entry_idx) = field2_global_row + 1_psb_lpk_
      entry_vals(entry_idx) = -1.0_psb_dpk_
    end if
  end do
  call nested_matrix%ins(2, 1, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  !---------------------------------------------------------------
  ! 3) assemble: nested_matrix%a_glob / nested_matrix%desc_glob are ready for Krylov
  !---------------------------------------------------------------
  call nested_matrix%asb(info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: nested_matrix%asb info=', info; goto 9999
  end if

  !---------------------------------------------------------------
  ! 4) consistent RHS: x_exact = 1, rhs = M * x_exact (via the nested operator)
  !---------------------------------------------------------------
  call psb_geall(x_exact, nested_matrix%desc_glob, info)
  call psb_geasb(x_exact, nested_matrix%desc_glob, info)
  call x_exact%set(done)                                   ! x_exact = 1 everywhere

  call psb_geall(rhs, nested_matrix%desc_glob, info); call psb_geasb(rhs, nested_matrix%desc_glob, info)
  call psb_spmm(done, nested_matrix%a_glob, x_exact, dzero, rhs, nested_matrix%desc_glob, info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_spmm (RHS) info=', info
    goto 9999
  end if
  norm_x_exact = psb_genrm2(x_exact, nested_matrix%desc_glob, info)

  !---------------------------------------------------------------
  ! 5) solve with the standard PSBLAS CG under every stock preconditioner
  !---------------------------------------------------------------
  if (my_rank == 0) write(*,'(a,i0,a,i0)') ' np=', num_procs, '  N(global)=', 2*field_size
  all_passed = .true.
  iter_none  = 0
  iter_diag  = -1
  do i_prec = 1, n_precs
    call preconditioner%init(context, trim(prec_names(i_prec)), info)
    call preconditioner%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
    if (info /= psb_success_) then
      if (my_rank == 0) write(*,*) 'FAIL: prec%build (', trim(prec_names(i_prec)), ') info=', info
      all_passed = .false.; exit
    end if

    call psb_geall(x_solution, nested_matrix%desc_glob, info)
    call psb_geasb(x_solution, nested_matrix%desc_glob, info)
    call psb_krylov('CG', nested_matrix%a_glob, preconditioner, rhs, x_solution, stop_tol, &
         & nested_matrix%desc_glob, info, &
         & itmax=max_iter, iter=n_iter, err=final_residual, itrace=trace_level, istop=stop_criterion)
    if (info /= psb_success_) then
      if (my_rank == 0) write(*,*) 'FAIL: psb_krylov(CG,', trim(prec_names(i_prec)), ') info=', info
      all_passed = .false.; exit
    end if

    ! solution error: || x_solution - x_exact || / || x_exact ||
    call psb_geaxpby(-done, x_exact, done, x_solution, nested_matrix%desc_glob, info)
    solution_error = psb_genrm2(x_solution, nested_matrix%desc_glob, info) / norm_x_exact

    if (my_rank == 0) then
      write(*,'(a,a6,a,i6,a,es12.4,a,es12.4)') ' prec=', prec_names(i_prec), &
           & '  CG iterations=', n_iter, '  residual=', final_residual, &
           & '  ||x-x_ex||/||x_ex||=', solution_error
    end if
    if ((n_iter >= max_iter) .or. (solution_error > solution_tol)) all_passed = .false.
    if (trim(prec_names(i_prec)) == 'NONE') iter_none = n_iter
    if (trim(prec_names(i_prec)) == 'DIAG') iter_diag = n_iter

    call psb_gefree(x_solution, nested_matrix%desc_glob, info)
    call preconditioner%free(info)
  end do

  !---------------------------------------------------------------
  ! 6) verdict: every preconditioner converges to the right solution, and DIAG
  !     reproduces the NONE iteration count exactly (Jacobi on the constant
  !     diagonal 2I is a pure rescaling -> exactness check of the nested get_diag)
  !---------------------------------------------------------------
  if (my_rank == 0) then
    if (all_passed .and. (iter_diag == iter_none)) then
      write(*,*) '[PASS] CG converges on the nested operator with NONE/DIAG/BJAC'
    else
      write(*,*) '[FAIL] preconditioned CG on the nested operator (tol ', solution_tol, ')'
    end if
  end if

  call nested_matrix%free(info)

9999 continue
  call psb_exit(context)

end program psb_d_nest_cg_test
