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
! File: psb_d_nest_base_mat_mod.F90
!
! Module: psb_d_nest_base_mat_mod
! Author: Simone Staccone (Stack-1)
!
! Adapter that makes a block-structured (nested) operator look like a standard
! local sparse matrix to PSBLAS: psb_d_nest_base_mat EXTENDS
! psb_d_base_sparse_mat and implements csmv (the local matrix-vector product).
! Wrapped in a psb_dspmat_type and paired with the composed global descriptor
! (see psb_cd_nest_compose), the nested operator can then be fed to psb_spmm,
! psb_krylov and the AMG4PSBLAS preconditioners unchanged (MATNEST-style).
!
! The local vector handed to csmv lives in the GLOBAL local layout produced by
! psb_cd_nest_compose: the owned entries of all fields are concatenated, followed
! by the global halo. For each field we precompute field_map(field)%global_local_pos,
! the positions in that global local vector of the field's own local vector
! (owned entries first, then the field's ghosts), so we can gather the field
! input sub-vector and scatter the field output sub-vector without further
! communication (the halo exchange is done once by psb_spmm on the global desc).
!
module psb_d_nest_base_mat_mod
  use psb_const_mod
  use psb_error_mod
  use psb_d_base_mat_mod, only : psb_d_base_sparse_mat
  use psb_desc_mod,       only : psb_desc_type
  use psb_desc_nest_mod,  only : psb_desc_nest_type
  use psb_d_nest_mat_mod, only : psb_d_nest_sparse_mat
  use psb_d_mat_mod,      only : psb_dspmat_type
  implicit none

  ! Per-field gather/scatter map into the global local vector.
  ! global_local_pos(1 : n_owned)        -> the field's owned entries
  ! global_local_pos(n_owned+1 : size)   -> the field's ghost (halo) entries
  type :: psb_d_nest_field_map
    integer(psb_ipk_)              :: n_owned = 0
    integer(psb_ipk_), allocatable :: global_local_pos(:)
  end type psb_d_nest_field_map

  type, extends(psb_d_base_sparse_mat) :: psb_d_nest_base_mat
    integer(psb_ipk_)                          :: n_fields = 0
    type(psb_d_nest_sparse_mat), pointer       :: block_storage => null()  ! blocks (not owned)
    type(psb_desc_nest_type),    pointer       :: grid_desc     => null()  ! per-field descriptors (not owned)
    type(psb_d_nest_field_map), allocatable    :: field_map(:)
  contains
    procedure, pass(a) :: csmv       => psb_d_nest_base_csmv
    procedure, pass(a) :: get_nzeros => psb_d_nest_base_get_nzeros
    procedure, nopass  :: get_fmt    => psb_d_nest_base_get_fmt
    procedure, pass(a) :: free       => psb_d_nest_base_free
  end type psb_d_nest_base_mat

  private
  public :: psb_d_nest_base_mat, psb_d_nest_base_setup, psb_d_nest_apply_block
  ! field-split interface (for the block preconditioner)
  public :: psb_d_nest_get_n_fields, psb_d_nest_get_field_owned, &
       &    psb_d_nest_get_block, psb_d_nest_get_field_desc,      &
       &    psb_d_nest_restrict_field, psb_d_nest_prolong_field

contains

  function psb_d_nest_base_get_fmt() result(format_name)
    character(len=5) :: format_name
    format_name = 'NEST'
  end function psb_d_nest_base_get_fmt

  ! free: the nested operator does NOT own block_storage / grid_desc (they are
  ! pointers into the caller), so we only detach them and release the field maps.
  subroutine psb_d_nest_base_free(a)
    class(psb_d_nest_base_mat), intent(inout) :: a
    a%block_storage => null()
    a%grid_desc     => null()
    if (allocated(a%field_map)) deallocate(a%field_map)
    a%n_fields = 0
    call a%set_null()
  end subroutine psb_d_nest_base_free

  function psb_d_nest_base_get_nzeros(a) result(total_nzeros)
    class(psb_d_nest_base_mat), intent(in) :: a
    integer(psb_ipk_) :: total_nzeros
    integer(psb_ipk_) :: i_block_row, j_block_col
    total_nzeros = 0
    if (associated(a%block_storage)) then
      do j_block_col = 1, a%block_storage%ncblocks
        do i_block_row = 1, a%block_storage%nrblocks
          if (a%block_storage%has_block(i_block_row, j_block_col)) &
               & total_nzeros = total_nzeros + &
               &   a%block_storage%mats(i_block_row, j_block_col)%get_nzeros()
        end do
      end do
    end if
  end function psb_d_nest_base_get_nzeros

  ! Build the per-field gather maps and set the local dimensions, from the nested
  ! grid descriptor (per-field distribution desc_grid%descs(1,field)) and the
  ! composed global descriptor desc_global (produced by psb_cd_nest_compose).
  subroutine psb_d_nest_base_setup(nest_op, block_storage, desc_grid, desc_global, info)
    type(psb_d_nest_base_mat),           intent(inout) :: nest_op
    type(psb_d_nest_sparse_mat), target, intent(in)    :: block_storage
    type(psb_desc_nest_type),    target, intent(in)    :: desc_grid
    type(psb_desc_type),                 intent(in)    :: desc_global
    integer(psb_ipk_),                   intent(out)   :: info

    integer(psb_ipk_)              :: n_fields, i_field, i_entry
    integer(psb_ipk_)              :: n_owned, n_local, n_ghost, owned_offset, local_pos
    integer(psb_lpk_)              :: global_idx
    integer(psb_lpk_), allocatable :: field_global_offset(:)
    character(len=24)              :: name

    info = psb_success_
    name = 'psb_d_nest_base_setup'

    if (desc_grid%nrblocks /= desc_grid%ncblocks) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='nested block structure must be square')
      return
    end if

    n_fields          = desc_grid%ncblocks
    nest_op%n_fields  = n_fields
    nest_op%grid_desc => desc_grid
    nest_op%block_storage => block_storage

    ! global field offsets (used to form ghost global indices)
    allocate(field_global_offset(n_fields+1), nest_op%field_map(n_fields), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
    end if
    field_global_offset(1) = 0
    do i_field = 1, n_fields
      field_global_offset(i_field+1) = field_global_offset(i_field) &
           & + desc_grid%descs(1,i_field)%get_global_rows()
    end do

    owned_offset = 0   ! running owned-local offset in the global local vector
    do i_field = 1, n_fields
      n_owned = desc_grid%descs(1,i_field)%get_local_rows()
      n_local = desc_grid%descs(1,i_field)%get_local_cols()
      n_ghost = n_local - n_owned
      nest_op%field_map(i_field)%n_owned = n_owned
      allocate(nest_op%field_map(i_field)%global_local_pos(n_local), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
      end if
      ! owned entries: contiguous in the global local vector
      do i_entry = 1, n_owned
        nest_op%field_map(i_field)%global_local_pos(i_entry) = owned_offset + i_entry
      end do
      ! ghost entries: locate the field's ghost global index in the global descriptor
      do i_entry = 1, n_ghost
        call desc_grid%descs(1,i_field)%l2g(n_owned + i_entry, global_idx, info)
        if (info /= 0) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='l2g'); return
        end if
        call desc_global%g2l(field_global_offset(i_field) + global_idx, local_pos, info)
        if (info /= 0) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='g2l'); return
        end if
        nest_op%field_map(i_field)%global_local_pos(n_owned + i_entry) = local_pos
      end do
      owned_offset = owned_offset + n_owned
    end do

    call nest_op%set_nrows(desc_global%get_local_rows())
    call nest_op%set_ncols(desc_global%get_local_cols())
    call nest_op%set_asb()

  end subroutine psb_d_nest_base_setup

  ! Local block matrix-vector product: y = alpha * A_nest * x + beta * y.
  ! x is in the global local layout (owned fields concatenated + global halo);
  ! y holds the owned entries (global local rows).
  subroutine psb_d_nest_base_csmv(alpha, a, x, beta, y, info, trans)
    real(psb_dpk_),             intent(in)    :: alpha, beta, x(:)
    class(psb_d_nest_base_mat), intent(in)    :: a
    real(psb_dpk_),             intent(inout) :: y(:)
    integer(psb_ipk_),          intent(out)   :: info
    character, optional,        intent(in)    :: trans

    real(psb_dpk_), allocatable :: x_field(:), y_field(:)
    integer(psb_ipk_)           :: i_block_row, j_block_col, i_entry
    integer(psb_ipk_)           :: n_local_col_field, n_owned_row_field
    character                   :: trans_op
    character(len=24)           :: name

    info = psb_success_
    name = 'psb_d_nest_base_csmv'
    trans_op = 'N'
    if (present(trans)) trans_op = trans
    if (trans_op /= 'N' .and. trans_op /= 'n') then
      ! Transposed nested product is not implemented (would swap block indices
      ! and need the transposed halo); reject explicitly. See P5.
      info = psb_err_transpose_not_n_unsupported_
      call psb_errpush(info, name)
      return
    end if
    if (.not. associated(a%block_storage)) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='nested operator not set up')
      return
    end if

    ! y <- beta * y
    if (beta == dzero) then
      y(:) = dzero
    else if (beta /= done) then
      y(:) = beta * y(:)
    end if

    do j_block_col = 1, a%n_fields
      n_local_col_field = size(a%field_map(j_block_col)%global_local_pos)
      allocate(x_field(n_local_col_field), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
      end if
      ! gather the column-field input sub-vector (owned + that field's ghosts)
      do i_entry = 1, n_local_col_field
        x_field(i_entry) = x(a%field_map(j_block_col)%global_local_pos(i_entry))
      end do

      do i_block_row = 1, a%n_fields
        if (a%block_storage%has_block(i_block_row, j_block_col)) then
          n_owned_row_field = a%field_map(i_block_row)%n_owned
          allocate(y_field(n_owned_row_field), stat=info)
          if (info /= 0) then
            info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
          end if
          ! current row-field output sub-vector (owned)
          do i_entry = 1, n_owned_row_field
            y_field(i_entry) = y(a%field_map(i_block_row)%global_local_pos(i_entry))
          end do
          ! y_field <- alpha * A(i_block_row, j_block_col) * x_field + y_field
          call a%block_storage%mats(i_block_row, j_block_col)%a%csmv( &
               & alpha, x_field, done, y_field, info, trans_op)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_from_subroutine_, name, a_err='block csmv')
            return
          end if
          ! scatter the row-field output sub-vector back into y
          do i_entry = 1, n_owned_row_field
            y(a%field_map(i_block_row)%global_local_pos(i_entry)) = y_field(i_entry)
          end do
          deallocate(y_field)
        end if
      end do
      deallocate(x_field)
    end do

  end subroutine psb_d_nest_base_csmv

  ! Selective (regime 2) application of a SINGLE block:
  !   y_field = alpha * A(i_block_row, j_block_col) * x_field + beta * y_field
  ! x_field is the column-field local vector (owned + ghosts) ALREADY halo-exchanged
  ! by the caller; y_field is the row-field owned local vector. The caller chooses
  ! the exchange regime (the union halo, or just this block's halo), so this
  ! routine is purely local. It is FORMAT-AGNOSTIC: it dispatches to the block's
  ! own polymorphic csmv, so the block may be CSR, COO, ... independently of the
  ! other blocks. (The full-operator matvec, regime 1, is psb_d_nest_base_csmv.)
  subroutine psb_d_nest_apply_block(nest_op, i_block_row, j_block_col, alpha, x_field, beta, y_field, info)
    type(psb_d_nest_base_mat), intent(in)    :: nest_op
    integer(psb_ipk_),         intent(in)    :: i_block_row, j_block_col
    real(psb_dpk_),            intent(in)    :: alpha, beta, x_field(:)
    real(psb_dpk_),            intent(inout) :: y_field(:)
    integer(psb_ipk_),         intent(out)   :: info
    character(len=24) :: name

    info = psb_success_
    name = 'psb_d_nest_apply_block'

    if (.not. associated(nest_op%block_storage)) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='nested operator not set up')
      return
    end if
    if (.not. nest_op%block_storage%has_block(i_block_row, j_block_col)) then
      ! absent block contributes zero: y_field <- beta * y_field
      if (beta == dzero) then
        y_field(:) = dzero
      else if (beta /= done) then
        y_field(:) = beta * y_field(:)
      end if
      return
    end if

    ! polymorphic dispatch: the block applies its own format (CSR/COO/...)
    call nest_op%block_storage%mats(i_block_row, j_block_col)%a%csmv( &
         & alpha, x_field, beta, y_field, info)
    if (info /= psb_success_) &
         & call psb_errpush(psb_err_from_subroutine_, name, a_err='block csmv')

  end subroutine psb_d_nest_apply_block

  ! ====================================================================
  !  Field-split interface (for the block preconditioner).
  !  Exposes the field structure so a fieldsplit/Schur preconditioner can:
  !   - know how many fields there are and their owned sizes;
  !   - get a block as a standard psb_dspmat_type (sub-preconditioner on A,
  !     Schur-complement matvecs with B / B^T);
  !   - get a field descriptor (run a field-level Krylov / halo exchange);
  !   - restrict the global vector to a field sub-vector and prolong it back.
  ! ====================================================================

  function psb_d_nest_get_n_fields(nest_op) result(n_fields)
    type(psb_d_nest_base_mat), intent(in) :: nest_op
    integer(psb_ipk_) :: n_fields
    n_fields = nest_op%n_fields
  end function psb_d_nest_get_n_fields

  function psb_d_nest_get_field_owned(nest_op, field) result(n_owned)
    type(psb_d_nest_base_mat), intent(in) :: nest_op
    integer(psb_ipk_),         intent(in) :: field
    integer(psb_ipk_) :: n_owned
    n_owned = 0
    if (allocated(nest_op%field_map) .and. field >= 1 .and. field <= nest_op%n_fields) &
         & n_owned = nest_op%field_map(field)%n_owned
  end function psb_d_nest_get_field_owned

  ! Pointer to block (i,j) as a standard psb_dspmat_type (null if absent).
  function psb_d_nest_get_block(nest_op, i_block_row, j_block_col) result(block_ptr)
    type(psb_d_nest_base_mat), target, intent(in) :: nest_op
    integer(psb_ipk_),                 intent(in) :: i_block_row, j_block_col
    type(psb_dspmat_type), pointer :: block_ptr
    block_ptr => null()
    if (associated(nest_op%block_storage)) then
      if (nest_op%block_storage%has_block(i_block_row, j_block_col)) &
           & block_ptr => nest_op%block_storage%mats(i_block_row, j_block_col)
    end if
  end function psb_d_nest_get_block

  ! Pointer to field k's descriptor (null if not set up).
  function psb_d_nest_get_field_desc(nest_op, field) result(desc_ptr)
    type(psb_d_nest_base_mat), target, intent(in) :: nest_op
    integer(psb_ipk_),                 intent(in) :: field
    type(psb_desc_type), pointer :: desc_ptr
    desc_ptr => null()
    if (associated(nest_op%grid_desc) .and. field >= 1 .and. field <= nest_op%n_fields) &
         & desc_ptr => nest_op%grid_desc%descs(1, field)
  end function psb_d_nest_get_field_desc

  ! Restrict: extract field k's OWNED sub-vector from the global local vector.
  subroutine psb_d_nest_restrict_field(nest_op, field, x_global, x_field, info)
    type(psb_d_nest_base_mat), intent(in)  :: nest_op
    integer(psb_ipk_),         intent(in)  :: field
    real(psb_dpk_),            intent(in)  :: x_global(:)
    real(psb_dpk_),            intent(out) :: x_field(:)
    integer(psb_ipk_),         intent(out) :: info
    integer(psb_ipk_) :: i_entry, n_owned
    info = psb_success_
    if (field < 1 .or. field > nest_op%n_fields) then
      info = psb_err_invalid_input_; return
    end if
    n_owned = nest_op%field_map(field)%n_owned
    do i_entry = 1, n_owned
      x_field(i_entry) = x_global(nest_op%field_map(field)%global_local_pos(i_entry))
    end do
  end subroutine psb_d_nest_restrict_field

  ! Prolong: insert field k's OWNED sub-vector into the global local vector.
  subroutine psb_d_nest_prolong_field(nest_op, field, x_field, x_global, info)
    type(psb_d_nest_base_mat), intent(in)    :: nest_op
    integer(psb_ipk_),         intent(in)    :: field
    real(psb_dpk_),            intent(in)    :: x_field(:)
    real(psb_dpk_),            intent(inout) :: x_global(:)
    integer(psb_ipk_),         intent(out)   :: info
    integer(psb_ipk_) :: i_entry, n_owned
    info = psb_success_
    if (field < 1 .or. field > nest_op%n_fields) then
      info = psb_err_invalid_input_; return
    end if
    n_owned = nest_op%field_map(field)%n_owned
    do i_entry = 1, n_owned
      x_global(nest_op%field_map(field)%global_local_pos(i_entry)) = x_field(i_entry)
    end do
  end subroutine psb_d_nest_prolong_field

end module psb_d_nest_base_mat_mod
