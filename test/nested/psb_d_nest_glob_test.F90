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
! File: psb_d_nest_glob_test.F90
!
! Program: psb_d_nest_glob_test
! Author:  Simone Staccone (Stack-1)
!
! Validates the "global nested operator" path built through the
! psb_d_nest_matrix utility (init/ins/asb): the user only supplies the field
! sizes and the block values, and obtains nested_matrix%a_glob /
! nested_matrix%desc_glob ready for psb_spmm.  The result is compared against
! the SAME matrix assembled monolithically in CSR on the same global
! descriptor (oracle).
!
! 2x2 operator (fields of size field_size):
!   [ A   B^T ]   A   = tridiag(-1, 2, -1)   (block 1,1)
!   [ B    0  ]   B^T = 0.5 * I              (block 1,2)
!                 B   = 0.3 * I              (block 2,1)
!                 (2,2) absent
!
! Run: ./psb_d_nest_glob_test               (serial)
!      mpirun -np 4 ./psb_d_nest_glob_test
!
program psb_d_nest_glob_test
  use psb_base_mod
  use psb_util_mod
  use psb_d_nest_mod
  use psb_d_hll_mat_mod, only : psb_d_hll_sparse_mat   ! psb_ext format for the blocks
  implicit none

  type(psb_d_hll_sparse_mat) :: hll_mold

  type(psb_ctxt_type)             :: context
  integer(psb_ipk_)               :: my_rank, num_procs, info, i_local_row
  integer(psb_ipk_)               :: entry_idx, field1_local_rows, field2_local_rows
  integer(psb_lpk_)               :: global_row, global_col, field_size

  type(psb_d_nest_matrix)         :: nested_matrix   ! the nested operator (init/ins/asb)
  type(psb_dspmat_type)           :: monolithic_ref  ! monolithic CSR oracle
  type(psb_d_vect_type)           :: x_vec, y_nested, y_monolithic

  integer(psb_lpk_), allocatable  :: entry_rows(:), entry_cols(:)
  real(psb_dpk_),    allocatable  :: entry_vals(:)
  real(psb_dpk_)                  :: insert_value(1)
  real(psb_dpk_)                  :: mismatch_norm
  real(psb_dpk_), parameter       :: tolerance = 1.0e-10_psb_dpk_

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  field_size = 32                           ! global size of each field

  !---------------------------------------------------------------
  ! 1) build the 2x2 nested operator through the utility
  !---------------------------------------------------------------
  call nested_matrix%init(context, [field_size, field_size], info)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: nested_matrix%init info=', info; goto 9999
  end if
  field1_local_rows = nested_matrix%field_desc(1)%get_local_rows()
  field2_local_rows = nested_matrix%field_desc(2)%get_local_rows()

  !---------------------------------------------------------------
  ! 2) insert the block values (owned rows only)
  !---------------------------------------------------------------
  ! A = tridiag(-1,2,-1)  -> block (1,1)
  allocate(entry_rows(3*field1_local_rows), entry_cols(3*field1_local_rows), &
       &   entry_vals(3*field1_local_rows))
  entry_idx = 0
  do i_local_row = 1, field1_local_rows
    call nested_matrix%field_desc(1)%l2g(i_local_row, global_row, info)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = global_row
    entry_cols(entry_idx) = global_row
    entry_vals(entry_idx) = 2.0_psb_dpk_
    if (global_row > 1) then
      entry_idx = entry_idx + 1
      entry_rows(entry_idx) = global_row
      entry_cols(entry_idx) = global_row - 1_psb_lpk_
      entry_vals(entry_idx) = -1.0_psb_dpk_
    end if
    if (global_row < field_size) then
      entry_idx = entry_idx + 1
      entry_rows(entry_idx) = global_row
      entry_cols(entry_idx) = global_row + 1_psb_lpk_
      entry_vals(entry_idx) = -1.0_psb_dpk_
    end if
  end do
  call nested_matrix%ins(1, 1, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! B^T = 0.5 I  -> block (1,2): rows in field 1, columns in field 2
  allocate(entry_rows(field1_local_rows), entry_cols(field1_local_rows), entry_vals(field1_local_rows))
  entry_idx = 0
  do i_local_row = 1, field1_local_rows
    call nested_matrix%field_desc(1)%l2g(i_local_row, global_row, info)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = global_row
    entry_cols(entry_idx) = global_row
    entry_vals(entry_idx) = 0.5_psb_dpk_
  end do
  call nested_matrix%ins(1, 2, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! B = 0.3 I  -> block (2,1): rows in field 2, columns in field 1
  allocate(entry_rows(field2_local_rows), entry_cols(field2_local_rows), entry_vals(field2_local_rows))
  entry_idx = 0
  do i_local_row = 1, field2_local_rows
    call nested_matrix%field_desc(2)%l2g(i_local_row, global_row, info)
    entry_idx = entry_idx + 1
    entry_rows(entry_idx) = global_row
    entry_cols(entry_idx) = global_row
    entry_vals(entry_idx) = 0.3_psb_dpk_
  end do
  call nested_matrix%ins(2, 1, entry_idx, entry_rows, entry_cols, entry_vals, info)
  deallocate(entry_rows, entry_cols, entry_vals)

  ! assemble with the blocks stored in HLL (psb_ext format): exercises the
  ! configurable block storage and the format-agnostic nested matvec
  call nested_matrix%asb(info, mold=hll_mold)
  if (info /= psb_success_) then
    if (my_rank==0) write(*,*) 'FAIL: nested_matrix%asb info=', info; goto 9999
  end if

  !---------------------------------------------------------------
  ! 3) monolithic oracle on nested_matrix%desc_glob (global offsets:
  !    field 1 -> g ;  field 2 -> field_size + g)
  !---------------------------------------------------------------
  call psb_spall(monolithic_ref, nested_matrix%desc_glob, info, &
       &         nnz=5*nested_matrix%desc_glob%get_local_rows())
  do i_local_row = 1, field1_local_rows                ! field-1 rows
    call nested_matrix%field_desc(1)%l2g(i_local_row, global_row, info)
    insert_value(1) = 2.0_psb_dpk_
    call psb_spins(1,[global_row],[global_row],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
    if (global_row > 1) then
      insert_value(1)=-1.0_psb_dpk_
      call psb_spins(1,[global_row],[global_row-1_psb_lpk_],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
    end if
    if (global_row < field_size) then
      insert_value(1)=-1.0_psb_dpk_
      call psb_spins(1,[global_row],[global_row+1_psb_lpk_],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
    end if
    global_col = field_size + global_row
    insert_value(1) = 0.5_psb_dpk_                                                  ! B^T
    call psb_spins(1,[global_row],[global_col],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
  end do
  do i_local_row = 1, field2_local_rows                ! field-2 rows
    call nested_matrix%field_desc(2)%l2g(i_local_row, global_row, info)
    global_col = global_row
    insert_value(1) = 0.3_psb_dpk_                                                  ! B
    call psb_spins(1,[field_size+global_row],[global_col],insert_value,monolithic_ref,nested_matrix%desc_glob,info)
  end do
  call psb_spasb(monolithic_ref, nested_matrix%desc_glob, info, dupl=psb_dupl_add_)

  !---------------------------------------------------------------
  ! 4) compare the two matrix-vector products on a distinct-valued x (x[g] = g)
  !---------------------------------------------------------------
  call psb_geall(x_vec, nested_matrix%desc_glob, info)
  do i_local_row = 1, nested_matrix%desc_glob%get_local_rows()
    call nested_matrix%desc_glob%l2g(i_local_row, global_row, info)
    insert_value(1) = real(global_row, psb_dpk_)
    call psb_geins(1, [global_row], insert_value, x_vec, nested_matrix%desc_glob, info)
  end do
  call psb_geasb(x_vec, nested_matrix%desc_glob, info)

  call psb_geall(y_nested,     nested_matrix%desc_glob, info); call psb_geasb(y_nested,     nested_matrix%desc_glob, info)
  call psb_geall(y_monolithic, nested_matrix%desc_glob, info); call psb_geasb(y_monolithic, nested_matrix%desc_glob, info)

  call psb_spmm(done, nested_matrix%a_glob, x_vec, dzero, y_nested,     nested_matrix%desc_glob, info)  ! via nested csmv
  if (info /= psb_success_) then
    if (my_rank == 0) write(*,*) 'FAIL: psb_spmm (nested) info=', info
    goto 9999
  end if
  call psb_spmm(done, monolithic_ref,      x_vec, dzero, y_monolithic, nested_matrix%desc_glob, info)  ! CSR oracle

  call psb_geaxpby(done, y_nested, -done, y_monolithic, nested_matrix%desc_glob, info)
  mismatch_norm = psb_genrm2(y_monolithic, nested_matrix%desc_glob, info)

  if (my_rank == 0) then
    write(*,'(a,i0,a,i0)') ' np=', num_procs, '  N(field)=', field_size
    write(*,'(a,es12.4)')  ' ||nested - monolithic||_2 = ', mismatch_norm
    if (mismatch_norm <= tolerance) then
      write(*,*) '[PASS] nested global operator matches monolithic CSR'
    else
      write(*,*) '[FAIL] mismatch above tolerance ', tolerance
    end if
  end if

  call nested_matrix%free(info)

9999 continue
  call psb_exit(context)

end program psb_d_nest_glob_test
