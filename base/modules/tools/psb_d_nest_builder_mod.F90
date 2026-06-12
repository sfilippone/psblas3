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
! File: psb_d_nest_builder_mod.F90
!
! Module: psb_d_nest_builder_mod
! Author: Simone Staccone (Stack-1)
!
! User-friendly frontend to build a nested (MATNEST) operator without manually
! managing per-field descriptors, the union halo, composition and setup.
!
! All the boilerplate (identical for every nested operator) is hidden behind a
! single type, psb_d_nest_matrix, with the usual PSBLAS init/ins/asb pattern:
!
!   type(psb_d_nest_matrix) :: nested_matrix
!   call nested_matrix%init(ctxt, [n1, n2], info)          ! 2 fields of global size n1, n2
!   call nested_matrix%ins(1,1, n, rows, cols, vals, info) ! values of block (1,1) = A
!   call nested_matrix%ins(1,2, n, rows, cols, vals, info) ! values of block (1,2) = B^T
!   call nested_matrix%ins(2,1, n, rows, cols, vals, info) ! values of block (2,1) = B
!   ...                                                    ! (absent blocks = not inserted)
!   call nested_matrix%asb(info)                           ! assemble: builds a_glob, desc_glob
!
!   ! from here on nested_matrix%a_glob and nested_matrix%desc_glob are an
!   ! ordinary distributed matrix/descriptor:
!   call psb_geall(x, nested_matrix%desc_glob, info)
!   call psb_krylov('CG', nested_matrix%a_glob, prec, b, x, eps, nested_matrix%desc_glob, info, ...)
!
! Indices: in ins(block_row, block_col, ...) the rows live in the index space of
! field block_row, the columns in the index space of field block_col (GLOBAL
! field indices, 1..field_size).  Each process inserts only the rows it owns
! (PSBLAS convention).  Off-diagonal blocks may be rectangular.
!
! NOTE: after asb the object holds consistent internal pointers (a_glob%a points
! to block_storage / grid_desc): do not copy/move the object after assembly.
!
module psb_d_nest_builder_mod
  use psb_const_mod
  use psb_error_mod,           only : psb_errpush
  use psb_penv_mod,            only : psb_ctxt_type, psb_info
  use psb_desc_mod,            only : psb_desc_type
  use psb_d_mat_mod,           only : psb_dspmat_type
  use psb_d_base_mat_mod,      only : psb_d_base_sparse_mat
  use psb_cd_tools_mod,        only : psb_cdall, psb_cdins, psb_cdasb
  use psb_desc_nest_mod,       only : psb_desc_nest_type
  use psb_d_nest_mat_mod,      only : psb_d_nest_sparse_mat
  use psb_d_nest_base_mat_mod, only : psb_d_nest_base_mat, psb_d_nest_base_setup
  use psb_cd_nest_tools_mod,   only : psb_cd_nest_compose
  use psb_d_nest_tools_mod,    only : psb_d_nest_rect_block
  implicit none

  ! growing triplet buffer for a single block
  type :: psb_d_nest_block_buffer
    integer(psb_ipk_)              :: n_entries = 0
    integer(psb_lpk_), allocatable :: entry_rows(:), entry_cols(:)
    real(psb_dpk_),    allocatable :: entry_vals(:)
  end type psb_d_nest_block_buffer

  type :: psb_d_nest_matrix
    type(psb_ctxt_type)                          :: context
    integer(psb_ipk_)                            :: n_fields  = 0
    logical                                      :: assembled = .false.
    ! construction state
    type(psb_desc_type),           allocatable   :: field_desc(:)        ! one descriptor per field
    type(psb_d_nest_block_buffer), allocatable   :: block_buffer(:,:)    ! triplets per block (i,j)
    ! products (owned; the pointers in a_glob%a point in here)
    type(psb_d_nest_sparse_mat)                  :: block_storage
    type(psb_desc_nest_type)                     :: grid_desc
    type(psb_dspmat_type)                        :: a_glob               ! the matrix to hand to Krylov
    type(psb_desc_type)                          :: desc_glob            ! the global descriptor
  contains
    procedure, pass(op) :: init => psb_d_nest_op_init
    procedure, pass(op) :: ins  => psb_d_nest_op_ins
    procedure, pass(op) :: asb  => psb_d_nest_op_asb
    procedure, pass(op) :: free => psb_d_nest_op_free
  end type psb_d_nest_matrix

  private
  public :: psb_d_nest_matrix

contains

  ! init: create one descriptor per field (block distribution from the global sizes)
  subroutine psb_d_nest_op_init(op, context, field_sizes, info)
    class(psb_d_nest_matrix), intent(inout) :: op
    type(psb_ctxt_type),        intent(in)    :: context
    integer(psb_lpk_),          intent(in)    :: field_sizes(:)
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_)  :: my_rank, num_procs, n_fields, i_field, field_local_rows
    integer(psb_lpk_)  :: field_global_size
    character(len=24)  :: name

    info = psb_success_
    name = 'psb_d_nest_op_init'

    call psb_info(context, my_rank, num_procs)
    n_fields     = size(field_sizes)
    op%context   = context
    op%n_fields  = n_fields
    op%assembled = .false.

    allocate(op%field_desc(n_fields), op%block_buffer(n_fields,n_fields), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
    end if

    do i_field = 1, n_fields
      field_global_size = field_sizes(i_field)
      ! block distribution: field_global_size rows over num_procs processes (total size invariant)
      field_local_rows = int(field_global_size / int(num_procs, psb_lpk_), psb_ipk_)
      if (int(my_rank, psb_lpk_) < mod(field_global_size, int(num_procs, psb_lpk_))) &
           & field_local_rows = field_local_rows + 1
      call psb_cdall(context, op%field_desc(i_field), info, nl=field_local_rows)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdall'); return
      end if
    end do
  end subroutine psb_d_nest_op_init

  ! ins: accumulate the triplets into block (block_row,block_col) and register the
  !      columns (field block_col index space) into that descriptor's union halo
  subroutine psb_d_nest_op_ins(op, block_row, block_col, n_entries, entry_rows, entry_cols, entry_vals, info)
    class(psb_d_nest_matrix), intent(inout) :: op
    integer(psb_ipk_),          intent(in)    :: block_row, block_col, n_entries
    integer(psb_lpk_),          intent(in)    :: entry_rows(:), entry_cols(:)
    real(psb_dpk_),             intent(in)    :: entry_vals(:)
    integer(psb_ipk_),          intent(out)   :: info
    character(len=24) :: name

    info = psb_success_
    name = 'psb_d_nest_op_ins'

    if (op%assembled) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='operator already assembled'); return
    end if
    if (block_row < 1 .or. block_row > op%n_fields .or. &
      & block_col < 1 .or. block_col > op%n_fields) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='block index out of range'); return
    end if
    if (n_entries <= 0) return

    call block_buffer_append(op%block_buffer(block_row,block_col), n_entries, &
         &                   entry_rows, entry_cols, entry_vals, info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
    end if

    ! the columns of block (block_row,block_col) live in field block_col ->
    ! register their indices into that descriptor's union halo
    ! (this also applies when block_col == block_row)
    call psb_cdins(n_entries, entry_cols(1:n_entries), op%field_desc(block_col), info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdins'); return
    end if
  end subroutine psb_d_nest_op_ins

  ! asb: assemble the descriptors, build the blocks, compose the global
  !      descriptor, set up the operator and wrap it into a_glob.
  !      The optional type ('CSR'/'CSC'/'COO', default 'CSR') or mold (any
  !      class extending psb_d_base_sparse_mat, e.g. the psb_ext ELL/HLL or
  !      the psb_cuda device formats) selects the storage format of the blocks.
  subroutine psb_d_nest_op_asb(op, info, type, mold)
    class(psb_d_nest_matrix), intent(inout), target :: op
    integer(psb_ipk_),          intent(out)           :: info
    character(len=*),           intent(in), optional  :: type
    class(psb_d_base_sparse_mat), intent(in), optional :: mold

    type(psb_d_nest_base_mat) :: nest_operator
    integer(psb_ipk_)         :: n_fields, i_field, j_field
    character(len=24)         :: name

    info = psb_success_
    name = 'psb_d_nest_op_asb'
    n_fields = op%n_fields

    ! 1) assemble the per-field descriptors (with the union halo accumulated in ins)
    do i_field = 1, n_fields
      call psb_cdasb(op%field_desc(i_field), info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdasb'); return
      end if
    end do

    ! 2) build the local blocks (generally rectangular) from the triplets
    op%block_storage%nrblocks = n_fields
    op%block_storage%ncblocks = n_fields
    allocate(op%block_storage%mats(n_fields,n_fields), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
    end if
    do j_field = 1, n_fields
      do i_field = 1, n_fields
        if (op%block_buffer(i_field,j_field)%n_entries > 0) then
          call psb_d_nest_rect_block(op%block_storage%mats(i_field,j_field),         &
               & op%block_buffer(i_field,j_field)%n_entries,                         &
               & op%block_buffer(i_field,j_field)%entry_rows,                        &
               & op%block_buffer(i_field,j_field)%entry_cols,                        &
               & op%block_buffer(i_field,j_field)%entry_vals,                        &
               & op%field_desc(i_field), op%field_desc(j_field), info,               &
               & type=type, mold=mold)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_from_subroutine_, name, a_err='rect_block'); return
          end if
        end if
      end do
    end do

    ! 3) descriptor grid: descs(i,j) = descriptor of field j
    op%grid_desc%nrblocks = n_fields
    op%grid_desc%ncblocks = n_fields
    allocate(op%grid_desc%descs(n_fields,n_fields), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
    end if
    do j_field = 1, n_fields
      do i_field = 1, n_fields
        call op%field_desc(j_field)%clone(op%grid_desc%descs(i_field,j_field), info)
      end do
    end do

    ! 4) composed global descriptor + operator setup
    call psb_cd_nest_compose(op%grid_desc, op%desc_glob, info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_, name, a_err='cd_nest_compose'); return
    end if
    call psb_d_nest_base_setup(nest_operator, op%block_storage, op%grid_desc, op%desc_glob, info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_, name, a_err='nest_base_setup'); return
    end if

    ! 5) wrap into the standard matrix object (the pointers keep pointing at op%*)
    allocate(op%a_glob%a, source=nest_operator, stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
    end if
    call op%a_glob%set_nrows(op%desc_glob%get_local_rows())
    call op%a_glob%set_ncols(op%desc_glob%get_local_cols())
    call op%a_glob%set_asb()

    ! 6) the triplet buffers are no longer needed
    do j_field = 1, n_fields
      do i_field = 1, n_fields
        call block_buffer_free(op%block_buffer(i_field,j_field))
      end do
    end do
    op%assembled = .true.
  end subroutine psb_d_nest_op_asb

  ! free: release everything
  subroutine psb_d_nest_op_free(op, info)
    class(psb_d_nest_matrix), intent(inout) :: op
    integer(psb_ipk_),          intent(out)   :: info
    integer(psb_ipk_) :: i_field, j_field, local_info

    info = psb_success_
    if (allocated(op%block_buffer)) then
      do j_field = 1, size(op%block_buffer,2)
        do i_field = 1, size(op%block_buffer,1)
          call block_buffer_free(op%block_buffer(i_field,j_field))
        end do
      end do
      deallocate(op%block_buffer, stat=local_info)
    end if
    if (op%assembled) then
      call op%a_glob%free()
      call op%desc_glob%free(local_info)
      call op%grid_desc%free(local_info)
    end if
    if (allocated(op%field_desc)) then
      do i_field = 1, size(op%field_desc)
        call op%field_desc(i_field)%free(local_info)
      end do
      deallocate(op%field_desc, stat=local_info)
    end if
    op%n_fields  = 0
    op%assembled = .false.
  end subroutine psb_d_nest_op_free

  !-----------------------------------------------------------------
  ! private helpers: growing triplet buffer
  !-----------------------------------------------------------------
  subroutine block_buffer_append(buffer, n_entries, entry_rows, entry_cols, entry_vals, info)
    type(psb_d_nest_block_buffer), intent(inout) :: buffer
    integer(psb_ipk_),             intent(in)    :: n_entries
    integer(psb_lpk_),             intent(in)    :: entry_rows(:), entry_cols(:)
    real(psb_dpk_),                intent(in)    :: entry_vals(:)
    integer(psb_ipk_),             intent(out)   :: info
    integer(psb_ipk_) :: required_size

    info = psb_success_
    required_size = buffer%n_entries + n_entries
    call ensure_capacity_lpk(buffer%entry_rows, required_size, info); if (info /= 0) return
    call ensure_capacity_lpk(buffer%entry_cols, required_size, info); if (info /= 0) return
    call ensure_capacity_dpk(buffer%entry_vals, required_size, info); if (info /= 0) return
    buffer%entry_rows(buffer%n_entries+1:required_size) = entry_rows(1:n_entries)
    buffer%entry_cols(buffer%n_entries+1:required_size) = entry_cols(1:n_entries)
    buffer%entry_vals(buffer%n_entries+1:required_size) = entry_vals(1:n_entries)
    buffer%n_entries = required_size
  end subroutine block_buffer_append

  subroutine ensure_capacity_lpk(array, required_size, info)
    integer(psb_lpk_), allocatable, intent(inout) :: array(:)
    integer(psb_ipk_),              intent(in)    :: required_size
    integer(psb_ipk_),              intent(out)   :: info
    integer(psb_lpk_), allocatable :: grown(:)
    integer(psb_ipk_) :: capacity

    info = 0
    if (.not. allocated(array)) then
      allocate(array(max(required_size,16)), stat=info); return
    end if
    capacity = size(array)
    if (required_size <= capacity) return
    allocate(grown(max(2*capacity, required_size)), stat=info); if (info /= 0) return
    grown(1:capacity) = array(1:capacity)
    call move_alloc(grown, array)
  end subroutine ensure_capacity_lpk

  subroutine ensure_capacity_dpk(array, required_size, info)
    real(psb_dpk_), allocatable, intent(inout) :: array(:)
    integer(psb_ipk_),           intent(in)    :: required_size
    integer(psb_ipk_),           intent(out)   :: info
    real(psb_dpk_), allocatable :: grown(:)
    integer(psb_ipk_) :: capacity

    info = 0
    if (.not. allocated(array)) then
      allocate(array(max(required_size,16)), stat=info); return
    end if
    capacity = size(array)
    if (required_size <= capacity) return
    allocate(grown(max(2*capacity, required_size)), stat=info); if (info /= 0) return
    grown(1:capacity) = array(1:capacity)
    call move_alloc(grown, array)
  end subroutine ensure_capacity_dpk

  subroutine block_buffer_free(buffer)
    type(psb_d_nest_block_buffer), intent(inout) :: buffer
    if (allocated(buffer%entry_rows)) deallocate(buffer%entry_rows)
    if (allocated(buffer%entry_cols)) deallocate(buffer%entry_cols)
    if (allocated(buffer%entry_vals)) deallocate(buffer%entry_vals)
    buffer%n_entries = 0
  end subroutine block_buffer_free

end module psb_d_nest_builder_mod
