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
! File: psb_d_nest_prec_cg_test.F90
!
! Program: psb_d_nest_prec_cg_test
!
! Same matrix/RHS/CG check as psb_d_nest_cg_test, but the preconditioner is the
! specialized nested preconditioner:
!
!   call preconditioner%init(context, 'NEST', info)
!   call preconditioner%set('COMPOSITION', 'ADDITIVE', info)
!   call preconditioner%set('BLOCK_SOLVE', 'DIAG', info)
!
! The matrix is the red-black reordered 1D Laplacian:
!
!   M = [ 2I    C  ]
!       [ C^T   2I ]
!
! With diagonal block solves this NEST/ADDITIVE preconditioner is a scalar
! rescaling of the identity on this problem, so CG should converge to the exact
! solution and reproduce the unpreconditioned iteration count.
!
program psb_d_nest_prec_cg_test
  use psb_base_mod
  use psb_util_mod
  use psb_prec_mod
  use psb_linsolve_mod
  use psb_d_nest_mod
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

  real(psb_dpk_)             :: diag_value, stop_tol, final_residual
  real(psb_dpk_)             :: norm_x_exact, solution_error
  integer(psb_ipk_)          :: max_iter, trace_level, n_iter, n_iter_none
  integer(psb_ipk_)          :: stop_criterion
  real(psb_dpk_), parameter  :: solution_tol = 1.0e-6_psb_dpk_
  logical                    :: all_passed

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  field_size     = 512
  diag_value     = 2.0_psb_dpk_
  stop_tol       = 1.0e-9_psb_dpk_
  max_iter       = 4000
  trace_level    = 0
  stop_criterion = 2
  all_passed     = .true.

  call nested_matrix%init(context, [field_size, field_size], info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: nested_matrix%init info=', info
    goto 9999
  end if
  field1_rows = nested_matrix%get_owned_rows(1)
  field2_rows = nested_matrix%get_owned_rows(2)
  field1_local_rows = size(field1_rows)
  field2_local_rows = size(field2_rows)

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

  ! block (1,2) = C
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

  ! block (2,1) = C^T
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

  call nested_matrix%asb(info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: nested_matrix%asb info=', info
    goto 9999
  end if

  call psb_geall(x_exact, nested_matrix%desc_glob, info)
  call psb_geasb(x_exact, nested_matrix%desc_glob, info)
  call x_exact%set(done)

  call psb_geall(rhs, nested_matrix%desc_glob, info)
  call psb_geasb(rhs, nested_matrix%desc_glob, info)
  call psb_spmm(done, nested_matrix%a_glob, x_exact, dzero, rhs, nested_matrix%desc_glob, info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_spmm (RHS) info=', info
    goto 9999
  end if
  norm_x_exact = psb_genrm2(x_exact, nested_matrix%desc_glob, info)

  if (my_rank == 0) write(*,'(a,i0,a,i0)') ' np=', num_procs, '  N(global)=', 2*field_size

  ! Baseline: no preconditioner.
  call preconditioner%init(context, 'NONE', info)
  call preconditioner%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
  call psb_geall(x_solution, nested_matrix%desc_glob, info)
  call psb_geasb(x_solution, nested_matrix%desc_glob, info)
  call psb_krylov('CG', nested_matrix%a_glob, preconditioner, rhs, x_solution, stop_tol, &
       & nested_matrix%desc_glob, info, &
       & itmax=max_iter, iter=n_iter_none, err=final_residual, itrace=trace_level, istop=stop_criterion)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_krylov(CG,NONE) info=', info
    all_passed = .false.
  end if
  if (my_rank == 0) write(*,'(a,i6,a,es12.4)') &
       & ' prec=NONE    CG iterations=', n_iter_none, '  residual=', final_residual
  call psb_gefree(x_solution, nested_matrix%desc_glob, info)
  call preconditioner%free(info)

  ! Specialized nested preconditioner: additive block-diagonal field split.
  call preconditioner%init(context, 'NEST', info)
  call preconditioner%set('COMPOSITION', 'ADDITIVE', info)
  call preconditioner%set('BLOCK_SOLVE', 'DIAG', info)
  call preconditioner%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: nested prec%build info=', info
    all_passed = .false.
    goto 9999
  end if

  call psb_geall(x_solution, nested_matrix%desc_glob, info)
  call psb_geasb(x_solution, nested_matrix%desc_glob, info)
  call psb_krylov('CG', nested_matrix%a_glob, preconditioner, rhs, x_solution, stop_tol, &
       & nested_matrix%desc_glob, info, &
       & itmax=max_iter, iter=n_iter, err=final_residual, itrace=trace_level, istop=stop_criterion)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_krylov(CG,NEST) info=', info
    all_passed = .false.
  end if

  call psb_geaxpby(-done, x_exact, done, x_solution, nested_matrix%desc_glob, info)
  solution_error = psb_genrm2(x_solution, nested_matrix%desc_glob, info) / norm_x_exact
  if (my_rank == 0) then
    write(*,'(a,i6,a,es12.4,a,es12.4)') &
         & ' prec=NEST    CG iterations=', n_iter, '  residual=', final_residual, &
         & '  ||x-x_ex||/||x_ex||=', solution_error
  end if

  if ((n_iter >= max_iter) .or. (solution_error > solution_tol)) all_passed = .false.
  if (n_iter /= n_iter_none) all_passed = .false.

  if (my_rank == 0) then
    if (all_passed) then
      write(*,*) '[PASS] CG converges with the specialized NEST preconditioner'
    else
      write(*,*) '[FAIL] CG with specialized NEST preconditioner'
    end if
  end if

  call psb_gefree(x_solution, nested_matrix%desc_glob, info)
  call preconditioner%free(info)
  call nested_matrix%free(info)

9999 continue
  call psb_exit(context)

end program psb_d_nest_prec_cg_test
