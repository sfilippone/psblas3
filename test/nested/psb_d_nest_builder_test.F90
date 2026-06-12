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
! File: psb_d_nest_builder_test.F90
!
! Program: psb_d_nest_builder_test
! Author:  Simone Staccone (Stack-1)
!
! Same operator as the low-level CG test (1D Laplacian reordered red-black, SPD
! and ill-conditioned) but built with the psb_d_nest_matrix utility: the user
! declares nested_matrix, gives the field sizes, inserts the block values and
! calls asb.  All the setup (per-field descriptors, union halo, compose, setup,
! wrap) is handled by the utility.  Solved with CG and checked against the
! exact solution.
!
!   M = [ 2I    C  ]   C(r,r) = -1 ,  C(r,r-1) = -1   (the Laplacian edges)
!       [ C^T   2I ]
!
! Run: ./psb_d_nest_builder_test ; mpirun -np 4 ./psb_d_nest_builder_test
!
program psb_d_nest_builder_test
  use psb_base_mod
  use psb_prec_mod
  use psb_linsolve_mod
  use psb_d_nest_mod          ! umbrella: includes psb_d_nest_matrix (builder)
  implicit none

  type(psb_ctxt_type)        :: context
  integer(psb_ipk_)          :: my_rank, num_procs, info, i_local_row, entry_idx
  integer(psb_ipk_)          :: field1_local_rows, field2_local_rows
  integer(psb_lpk_)          :: field1_global_row, field2_global_row, field_size

  type(psb_d_nest_matrix)    :: nested_matrix       ! the only object needed
  type(psb_dprec_type)       :: preconditioner
  type(psb_d_vect_type)      :: x_solution, rhs, x_exact
  real(psb_dpk_)             :: insert_value(1)

  integer(psb_lpk_), allocatable :: entry_rows(:), entry_cols(:)
  integer(psb_lpk_), allocatable :: field1_rows(:), field2_rows(:)
  real(psb_dpk_),    allocatable :: entry_vals(:)

  real(psb_dpk_)             :: stop_tol, final_residual, norm_x_exact, solution_error
  integer(psb_ipk_)          :: max_iter, n_iter, stop_criterion
  real(psb_dpk_), parameter  :: solution_tol = 1.0e-6_psb_dpk_

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  field_size     = 512                      ! global size of each field (N = 2*field_size)
  stop_tol       = 1.0e-9_psb_dpk_
  max_iter       = 4000
  stop_criterion = 2

  !---------------------------------------------------------------
  ! 1) create the nested operator: 2 fields of global size field_size
  !---------------------------------------------------------------
  call nested_matrix%init(context, [field_size, field_size], info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL init info=', info; goto 9999
  end if

  ! rows owned by this process in each field
  field1_rows = nested_matrix%get_owned_rows(1)
  field2_rows = nested_matrix%get_owned_rows(2)
  field1_local_rows = size(field1_rows)
  field2_local_rows = size(field2_rows)

  !---------------------------------------------------------------
  ! 2) insert the values, one block at a time (owned rows only)
  !---------------------------------------------------------------
  ! block (1,1) = 2I
  allocate(entry_rows(field1_local_rows), entry_cols(field1_local_rows), entry_vals(field1_local_rows))
  do i_local_row = 1, field1_local_rows
    field1_global_row = field1_rows(i_local_row)
    entry_rows(i_local_row)=field1_global_row; entry_cols(i_local_row)=field1_global_row
    entry_vals(i_local_row)=2.0_psb_dpk_
  end do
  call nested_matrix%ins(1, 1, field1_local_rows, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! block (2,2) = 2I
  allocate(entry_rows(field2_local_rows), entry_cols(field2_local_rows), entry_vals(field2_local_rows))
  do i_local_row = 1, field2_local_rows
    field2_global_row = field2_rows(i_local_row)
    entry_rows(i_local_row)=field2_global_row; entry_cols(i_local_row)=field2_global_row
    entry_vals(i_local_row)=2.0_psb_dpk_
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
  ! 3) assemble: from here nested_matrix%a_glob / nested_matrix%desc_glob are ready for Krylov
  !---------------------------------------------------------------
  call nested_matrix%asb(info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL asb info=', info; goto 9999
  end if

  !---------------------------------------------------------------
  ! 4) consistent RHS x_exact=1, rhs = M*x_exact, then solve with standard CG
  !---------------------------------------------------------------
  call psb_geall(x_exact, nested_matrix%desc_glob, info)
  do i_local_row = 1, nested_matrix%desc_glob%get_local_rows()
    call nested_matrix%desc_glob%l2g(i_local_row, field1_global_row, info)
    insert_value(1) = 1.0_psb_dpk_
    call psb_geins(1, [field1_global_row], insert_value, x_exact, nested_matrix%desc_glob, info)
  end do
  call psb_geasb(x_exact, nested_matrix%desc_glob, info)

  call psb_geall(rhs, nested_matrix%desc_glob, info); call psb_geasb(rhs, nested_matrix%desc_glob, info)
  call psb_spmm(done, nested_matrix%a_glob, x_exact, dzero, rhs, nested_matrix%desc_glob, info)
  norm_x_exact = psb_genrm2(x_exact, nested_matrix%desc_glob, info)

  call preconditioner%init(context, 'NONE', info)
  call preconditioner%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)

  call psb_geall(x_solution, nested_matrix%desc_glob, info); call psb_geasb(x_solution, nested_matrix%desc_glob, info)
  call psb_krylov('CG', nested_matrix%a_glob, preconditioner, rhs, x_solution, stop_tol, nested_matrix%desc_glob, info, &
       & itmax=max_iter, iter=n_iter, err=final_residual, istop=stop_criterion)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL krylov info=', info; goto 9999
  end if

  call psb_geaxpby(-done, x_exact, done, x_solution, nested_matrix%desc_glob, info)
  solution_error = psb_genrm2(x_solution, nested_matrix%desc_glob, info) / norm_x_exact

  if (my_rank == 0) then
    write(*,'(a,i0,a,i0)')  ' np=', num_procs, '  N(global)=', 2*field_size
    write(*,'(a,i0)')       ' CG iterations            = ', n_iter
    write(*,'(a,es12.4)')   ' CG relative residual     = ', final_residual
    write(*,'(a,es12.4)')   ' ||x - x_exact||/||x_ex|| = ', solution_error
    if ((n_iter < max_iter) .and. (solution_error <= solution_tol)) then
      write(*,*) '[PASS] nested matrix built with the utility, solved with CG'
    else
      write(*,*) '[FAIL] tol ', solution_tol
    end if
  end if

  call nested_matrix%free(info)

9999 continue
  call psb_exit(context)

end program psb_d_nest_builder_test
