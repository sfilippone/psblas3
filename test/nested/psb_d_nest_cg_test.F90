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
! Solves a linear system with the GLOBAL nested operator using the standard
! PSBLAS CG (psb_krylov('CG', ...)).  This test builds the operator on the
! LOW-LEVEL path (per-field descriptors, blocks, compose, setup, wrap) to
! directly validate the machinery the psb_d_nest_matrix utility relies on; the
! same solve through the utility is in psb_d_nest_builder_test.
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
! Run: ./psb_d_nest_cg_test ; mpirun -np 4 ./psb_d_nest_cg_test
!
program psb_d_nest_cg_test
  use psb_base_mod
  use psb_util_mod
  use psb_prec_mod
  use psb_linsolve_mod
  use psb_d_nest_mod      
  implicit none

  type(psb_ctxt_type)             :: context
  integer(psb_ipk_)               :: my_rank, num_procs, info, i_local_row, entry_idx, field_local_rows
  integer(psb_lpk_)               :: field1_global_row, field2_global_row, field_size

  ! per-field descriptors + blocks
  type(psb_desc_type)             :: field1_desc, field2_desc
  type(psb_dspmat_type)           :: diag_block1, coupling_12, coupling_21, diag_block2

  ! nested storage + grid descriptor + composed global path
  type(psb_d_nest_sparse_mat)     :: block_storage
  type(psb_desc_nest_type)        :: grid_desc
  type(psb_desc_type)             :: desc_global
  type(psb_d_nest_base_mat)       :: nest_operator
  type(psb_dspmat_type)           :: global_operator

  ! preconditioner + vectors
  type(psb_dprec_type)            :: preconditioner
  type(psb_d_vect_type)           :: x_solution, rhs, x_exact
  real(psb_dpk_)                  :: insert_value(1)

  ! global triplets for the coupling blocks
  integer(psb_lpk_), allocatable  :: entry_rows(:), entry_cols(:)
  real(psb_dpk_),    allocatable  :: entry_vals(:)

  ! solver parameters
  real(psb_dpk_)                  :: diag_value, stop_tol, final_residual, norm_x_exact, solution_error
  integer(psb_ipk_)               :: max_iter, trace_level, n_iter, stop_criterion
  real(psb_dpk_), parameter       :: solution_tol = 1.0e-6_psb_dpk_

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  field_size  = 512                         ! global rows per field (global N = 2*field_size)
  diag_value  = 2.0_psb_dpk_                ! Laplacian diagonal (diagonal blocks = diag*I)
  stop_tol    = 1.0e-9_psb_dpk_
  max_iter    = 4000
  trace_level = 0
  stop_criterion = 2                        ! stop on the relative residual

  !---------------------------------------------------------------
  ! 1) per-field descriptors: block distribution of field_size global rows
  !    over num_procs processes (total size independent of num_procs)
  !---------------------------------------------------------------
  field_local_rows = int(field_size / int(num_procs, psb_lpk_), psb_ipk_)
  if (int(my_rank, psb_lpk_) < mod(field_size, int(num_procs, psb_lpk_))) &
       & field_local_rows = field_local_rows + 1
  call psb_cdall(context, field1_desc, info, nl=field_local_rows)
  call psb_cdall(context, field2_desc, info, nl=field_local_rows)

  !---------------------------------------------------------------
  ! 2) diagonal blocks A = B = diag*I  (odd/even nodes of the red-black
  !    reordered Laplacian are not adjacent to each other)
  !---------------------------------------------------------------
  call psb_spall(diag_block1, field1_desc, info, nnz=field1_desc%get_local_rows())
  call psb_spall(diag_block2, field2_desc, info, nnz=field2_desc%get_local_rows())

  do i_local_row = 1, field1_desc%get_local_rows()
    call field1_desc%l2g(i_local_row, field1_global_row, info)
    insert_value(1) = diag_value
    call psb_spins(1,[field1_global_row],[field1_global_row],insert_value,diag_block1,field1_desc,info)
  end do
  do i_local_row = 1, field2_desc%get_local_rows()
    call field2_desc%l2g(i_local_row, field2_global_row, info)
    insert_value(1) = diag_value
    call psb_spins(1,[field2_global_row],[field2_global_row],insert_value,diag_block2,field2_desc,info)
  end do

  !---------------------------------------------------------------
  ! 3) register, in the union halo, the cross-field columns of the coupling blocks
  !    C   (row field1, col field2): columns {r, r-1} in field2  -> into field2_desc
  !    C^T (row field2, col field1): columns {s, s+1} in field1  -> into field1_desc
  !---------------------------------------------------------------
  do i_local_row = 1, field1_desc%get_local_rows()
    call field1_desc%l2g(i_local_row, field1_global_row, info)
    call psb_cdins(1, [field1_global_row], field2_desc, info)
    if (field1_global_row > 1) call psb_cdins(1, [field1_global_row-1_psb_lpk_], field2_desc, info)
  end do
  do i_local_row = 1, field2_desc%get_local_rows()
    call field2_desc%l2g(i_local_row, field2_global_row, info)
    call psb_cdins(1, [field2_global_row], field1_desc, info)
    if (field2_global_row < field_size) call psb_cdins(1, [field2_global_row+1_psb_lpk_], field1_desc, info)
  end do

  call psb_cdasb(field1_desc, info)
  call psb_cdasb(field2_desc, info)
  call psb_spasb(diag_block1, field1_desc, info, dupl=psb_dupl_add_)
  call psb_spasb(diag_block2, field2_desc, info, dupl=psb_dupl_add_)

  !---------------------------------------------------------------
  ! 4) coupling C (1,2): rows field1 (field1_desc), columns field2 (field2_desc)
  !    C(r,r) = -1 ,  C(r,r-1) = -1   (odd node 2r-1 -> even nodes 2r and 2r-2)
  !---------------------------------------------------------------
  allocate(entry_rows(2*field1_desc%get_local_rows()), entry_cols(2*field1_desc%get_local_rows()), &
       &   entry_vals(2*field1_desc%get_local_rows()))
  entry_idx = 0
  do i_local_row = 1, field1_desc%get_local_rows()
    call field1_desc%l2g(i_local_row, field1_global_row, info)
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
  call psb_d_nest_rect_block(coupling_12, entry_idx, entry_rows, entry_cols, entry_vals, field1_desc, field2_desc, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  !---------------------------------------------------------------
  ! 5) coupling C^T (2,1) = exact transpose of C:
  !    rows field2 (field2_desc), columns field1 (field1_desc)
  !    C^T(s,s) = -1 ,  C^T(s,s+1) = -1   (even node 2s -> odd nodes 2s-1 and 2s+1)
  !---------------------------------------------------------------
  allocate(entry_rows(2*field2_desc%get_local_rows()), entry_cols(2*field2_desc%get_local_rows()), &
       &   entry_vals(2*field2_desc%get_local_rows()))
  entry_idx = 0
  do i_local_row = 1, field2_desc%get_local_rows()
    call field2_desc%l2g(i_local_row, field2_global_row, info)
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
  call psb_d_nest_rect_block(coupling_21, entry_idx, entry_rows, entry_cols, entry_vals, field2_desc, field1_desc, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  !---------------------------------------------------------------
  ! 6) nested grid (all four blocks present)
  !---------------------------------------------------------------
  block_storage%nrblocks = 2
  block_storage%ncblocks = 2
  allocate(block_storage%mats(2,2))
  call psb_move_alloc(diag_block1, block_storage%mats(1,1), info)
  call psb_move_alloc(coupling_12, block_storage%mats(1,2), info)
  call psb_move_alloc(coupling_21, block_storage%mats(2,1), info)
  call psb_move_alloc(diag_block2, block_storage%mats(2,2), info)

  grid_desc%nrblocks = 2
  grid_desc%ncblocks = 2
  allocate(grid_desc%descs(2,2))
  call field1_desc%clone(grid_desc%descs(1,1), info)
  call field2_desc%clone(grid_desc%descs(1,2), info)
  call field1_desc%clone(grid_desc%descs(2,1), info)
  call field2_desc%clone(grid_desc%descs(2,2), info)

  !---------------------------------------------------------------
  ! 7) composed global operator (what CG will use as its matrix)
  !---------------------------------------------------------------
  call psb_cd_nest_compose(grid_desc, desc_global, info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_cd_nest_compose info=', info
    goto 9999
  end if
  call psb_d_nest_base_setup(nest_operator, block_storage, grid_desc, desc_global, info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_d_nest_base_setup info=', info
    goto 9999
  end if
  allocate(global_operator%a, source=nest_operator)
  call global_operator%set_nrows(desc_global%get_local_rows())
  call global_operator%set_ncols(desc_global%get_local_cols())
  call global_operator%set_asb()

  !---------------------------------------------------------------
  ! 8) consistent RHS: x_exact = 1, rhs = M * x_exact (via the nested operator)
  !---------------------------------------------------------------
  call psb_geall(x_exact, desc_global, info)
  do i_local_row = 1, desc_global%get_local_rows()
    call desc_global%l2g(i_local_row, field1_global_row, info)
    insert_value(1) = 1.0_psb_dpk_
    call psb_geins(1, [field1_global_row], insert_value, x_exact, desc_global, info)
  end do
  call psb_geasb(x_exact, desc_global, info)

  call psb_geall(rhs, desc_global, info); call psb_geasb(rhs, desc_global, info)
  call psb_spmm(done, global_operator, x_exact, dzero, rhs, desc_global, info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_spmm (RHS) info=', info
    goto 9999
  end if

  norm_x_exact = psb_genrm2(x_exact, desc_global, info)

  !---------------------------------------------------------------
  ! 9) identity preconditioner (NONE): CG exercises only the operator
  !---------------------------------------------------------------
  call preconditioner%init(context, 'NONE', info)
  call preconditioner%build(global_operator, desc_global, info)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: preconditioner%build info=', info
    goto 9999
  end if

  !---------------------------------------------------------------
  ! 10) solve with the standard PSBLAS CG
  !---------------------------------------------------------------
  call psb_geall(x_solution, desc_global, info); call psb_geasb(x_solution, desc_global, info)
  call psb_krylov('CG', global_operator, preconditioner, rhs, x_solution, stop_tol, desc_global, info, &
       & itmax=max_iter, iter=n_iter, err=final_residual, itrace=trace_level, istop=stop_criterion)
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_krylov(CG) info=', info
    goto 9999
  end if

  !---------------------------------------------------------------
  ! 11) solution error: || x_solution - x_exact || / || x_exact ||
  !---------------------------------------------------------------
  call psb_geaxpby(-done, x_exact, done, x_solution, desc_global, info)   ! x_solution <- x_solution - x_exact
  solution_error = psb_genrm2(x_solution, desc_global, info) / norm_x_exact

  if (my_rank == 0) then
    write(*,'(a,i0,a,i0)')  ' np=', num_procs, '  N(global)=', 2*field_size
    write(*,'(a,i0)')       ' CG iterations            = ', n_iter
    write(*,'(a,es12.4)')   ' CG relative residual     = ', final_residual
    write(*,'(a,es12.4)')   ' ||x - x_exact||/||x_ex|| = ', solution_error
    if ((n_iter < max_iter) .and. (solution_error <= solution_tol)) then
      write(*,*) '[PASS] CG converges on the global nested operator'
    else
      write(*,*) '[FAIL] CG does not converge / wrong solution (tol ', solution_tol, ')'
    end if
  end if

9999 continue
  call psb_exit(context)

end program psb_d_nest_cg_test
