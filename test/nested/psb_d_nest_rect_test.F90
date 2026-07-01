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
! File: psb_d_nest_rect_test.F90
!
! Program: psb_d_nest_rect_test
! Author:  Simone Staccone (Stack-1)
!
! Like psb_d_nest_glob_test but with fields of DIFFERENT size (|V| = 2|Q|) and
! GENUINELY RECTANGULAR off-diagonal blocks.  The operator is built with the
! psb_d_nest_matrix utility (init/ins/asb): the user only inserts the values,
! while the cross-field halo registration and the construction of the
! rectangular local blocks are handled internally.  Compared against the
! monolithic CSR oracle.
!
!   [ A   B^T ]   A   : V x V  tridiag(-1,2,-1)
!   [ B    0  ]   B^T : V x Q  rectangular  (row r -> col mod(r-1,nQ)+1, val 0.5)
!                 B   : Q x V  rectangular  (row q -> cols q and q+nQ, val 0.3)
!                 (2,2) absent
!
! Run: ./psb_d_nest_rect_test ; mpirun -np 4 ./psb_d_nest_rect_test
!
program psb_d_nest_rect_test
  use psb_base_mod
  use psb_util_mod
  use psb_d_nest_mod
  implicit none

  type(psb_ctxt_type)             :: context
  integer(psb_ipk_)               :: my_rank, num_procs, info, i_local_row
  integer(psb_ipk_)               :: entry_idx, v_local_rows, q_local_rows
  integer(psb_lpk_)               :: v_global_row, q_global_row, q_col, v_size, q_size

  type(psb_d_nest_matrix)         :: nested_matrix
  type(psb_dspmat_type)           :: monolithic_ref
  type(psb_d_vect_type)           :: x_vec, y_nested, y_monolithic
  real(psb_dpk_)                  :: insert_value(1)

  integer(psb_lpk_), allocatable  :: entry_rows(:), entry_cols(:)
  integer(psb_lpk_), allocatable  :: v_rows(:), q_rows(:)
  real(psb_dpk_),    allocatable  :: entry_vals(:)
  real(psb_dpk_)                  :: mismatch_norm
  real(psb_dpk_), parameter       :: tolerance = 1.0e-10_psb_dpk_

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  q_size = 16                               ! global size of field Q
  v_size = 2*q_size                         ! global size of field V (|V| = 2|Q|)

  !---------------------------------------------------------------
  ! 1) build the 2x2 nested operator (fields V, Q)
  !---------------------------------------------------------------
  call nested_matrix%init(context, [v_size, q_size], info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: nested_matrix%init info=', info; goto 9999
  end if
  v_rows = nested_matrix%get_owned_rows(1)
  q_rows = nested_matrix%get_owned_rows(2)
  v_local_rows = size(v_rows)
  q_local_rows = size(q_rows)

  !---------------------------------------------------------------
  ! 2) insert the blocks (owned rows only)
  !---------------------------------------------------------------
  ! A = tridiag(-1,2,-1)  -> block (1,1), V x V
  allocate(entry_rows(3*v_local_rows), entry_cols(3*v_local_rows), entry_vals(3*v_local_rows))
  entry_idx = 0
  do i_local_row = 1, v_local_rows
    v_global_row = v_rows(i_local_row)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = v_global_row
    entry_cols(entry_idx) = v_global_row
    entry_vals(entry_idx) = 2.0_psb_dpk_
    if (v_global_row > 1) then
      entry_idx = entry_idx + 1
      entry_rows(entry_idx) = v_global_row
      entry_cols(entry_idx) = v_global_row - 1_psb_lpk_
      entry_vals(entry_idx) = -1.0_psb_dpk_
    end if
    if (v_global_row < v_size) then
      entry_idx = entry_idx + 1
      entry_rows(entry_idx) = v_global_row
      entry_cols(entry_idx) = v_global_row + 1_psb_lpk_
      entry_vals(entry_idx) = -1.0_psb_dpk_
    end if
  end do
  call nested_matrix%ins(1, 1, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! B^T rectangular  -> block (1,2), V x Q : row r -> col mod(r-1,nQ)+1, val 0.5
  allocate(entry_rows(v_local_rows), entry_cols(v_local_rows), entry_vals(v_local_rows))
  entry_idx = 0
  do i_local_row = 1, v_local_rows
    v_global_row = v_rows(i_local_row)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = v_global_row
    entry_cols(entry_idx) = mod(v_global_row-1_psb_lpk_, q_size)+1
    entry_vals(entry_idx) = 0.5_psb_dpk_
  end do
  call nested_matrix%ins(1, 2, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! B rectangular  -> block (2,1), Q x V : row q -> cols q and q+nQ, val 0.3
  allocate(entry_rows(2*q_local_rows), entry_cols(2*q_local_rows), entry_vals(2*q_local_rows))
  entry_idx = 0
  do i_local_row = 1, q_local_rows
    q_global_row = q_rows(i_local_row)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = q_global_row
    entry_cols(entry_idx) = q_global_row
    entry_vals(entry_idx) = 0.3_psb_dpk_
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = q_global_row
    entry_cols(entry_idx) = q_global_row + q_size
    entry_vals(entry_idx) = 0.3_psb_dpk_
  end do
  call nested_matrix%ins(2, 1, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! assemble with the blocks stored in CSC instead of the CSR default:
  ! exercises the configurable block storage on a base format
  call nested_matrix%asb(info, type='CSC')
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: nested_matrix%asb info=', info; goto 9999
  end if

  !---------------------------------------------------------------
  ! 3) monolithic oracle on nested_matrix%desc_glob (global offsets: V -> g, Q -> v_size + g)
  !---------------------------------------------------------------
  call psb_spall(monolithic_ref, nested_matrix%desc_glob, info, &
       &         nnz=6*nested_matrix%desc_glob%get_local_rows())
  do i_local_row = 1, v_local_rows                     ! V rows
    v_global_row = v_rows(i_local_row)
    insert_value(1)=2.0_psb_dpk_
    call psb_spins(1,[v_global_row],[v_global_row],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
    if (v_global_row>1) then
      insert_value(1)=-1.0_psb_dpk_
      call psb_spins(1,[v_global_row],[v_global_row-1_psb_lpk_],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
    end if
    if (v_global_row<v_size) then
      insert_value(1)=-1.0_psb_dpk_
      call psb_spins(1,[v_global_row],[v_global_row+1_psb_lpk_],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
    end if
    q_col = v_size + (mod(v_global_row-1_psb_lpk_, q_size) + 1)
    insert_value(1)=0.5_psb_dpk_                                                    ! B^T
    call psb_spins(1,[v_global_row],[q_col],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
  end do
  do i_local_row = 1, q_local_rows                     ! Q rows
    q_global_row = q_rows(i_local_row)
    insert_value(1)=0.3_psb_dpk_
    call psb_spins(1,[v_size+q_global_row],[q_global_row],       insert_value,monolithic_ref,nested_matrix%desc_glob,info)  ! col q
    call psb_spins(1,[v_size+q_global_row],[q_global_row+q_size],insert_value,monolithic_ref,nested_matrix%desc_glob,info)  ! col q+nQ
  end do
  call psb_spasb(monolithic_ref, nested_matrix%desc_glob, info, dupl=psb_dupl_add_)

  !---------------------------------------------------------------
  ! 4) compare the two matrix-vector products on x[g] = g
  !---------------------------------------------------------------
  call psb_geall(x_vec, nested_matrix%desc_glob, info)
  do i_local_row = 1, nested_matrix%desc_glob%get_local_rows()
    call nested_matrix%desc_glob%l2g(i_local_row, v_global_row, info)
    insert_value(1) = real(v_global_row, psb_dpk_)
    call psb_geins(1, [v_global_row], insert_value, x_vec, nested_matrix%desc_glob, info)
  end do
  call psb_geasb(x_vec, nested_matrix%desc_glob, info)
  call psb_geall(y_nested,     nested_matrix%desc_glob, info); call psb_geasb(y_nested,     nested_matrix%desc_glob, info)
  call psb_geall(y_monolithic, nested_matrix%desc_glob, info); call psb_geasb(y_monolithic, nested_matrix%desc_glob, info)

  call psb_spmm(done, nested_matrix%a_glob, x_vec, dzero, y_nested, nested_matrix%desc_glob, info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL spmm nested info=', info; goto 9999
  end if
  call psb_spmm(done, monolithic_ref, x_vec, dzero, y_monolithic, nested_matrix%desc_glob, info)

  call psb_geaxpby(done, y_nested, -done, y_monolithic, nested_matrix%desc_glob, info)
  mismatch_norm = psb_genrm2(y_monolithic, nested_matrix%desc_glob, info)

  if (my_rank == 0) then
    write(*,'(a,i0,a,i0,a,i0)') ' np=', num_procs, '  |V|=', v_size, '  |Q|=', q_size
    write(*,'(a,es12.4)')       ' ||nested - monolithic||_2 = ', mismatch_norm
    if (mismatch_norm <= tolerance) then
      write(*,*) '[PASS] rectangular nested operator matches monolithic CSR'
    else
      write(*,*) '[FAIL] mismatch above tolerance ', tolerance
    end if
  end if

  call nested_matrix%free(info)

9999 continue
  call psb_exit(context)

end program psb_d_nest_rect_test
