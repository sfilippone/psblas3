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
  use psb_realloc_mod,    only : psb_ensure_size
  use psb_d_base_mat_mod, only : psb_d_base_sparse_mat
  use psb_d_base_vect_mod, only : psb_d_base_vect_type
  use psb_i_vect_mod,     only : psb_i_vect_type
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
    ! same positions as an encapsulated index vector, for the device-capable
    ! gather/scatter (gth/sct) used by vect_mv; pointer so that its target can
    ! be synced even when the operator dummy argument is intent(in)
    type(psb_i_vect_type), pointer :: gather_pos => null()
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
    ! enable the stock PSBLAS preconditioners on the nested operator:
    ! get_diag is used by DIAG/JACOBI, csgetrow by BJAC (ILU factorizations
    ! go through the format-agnostic csget path)
    procedure, pass(a) :: get_diag   => psb_d_nest_base_get_diag
    procedure, pass(a) :: csgetrow   => psb_d_nest_base_csgetrow
    ! device-capable matvec on encapsulated vectors: gathers/scatters through
    ! the vectors' own gth/sct and runs each block through its vect_mv, so
    ! device block formats execute their device kernels
    procedure, pass(a) :: vect_mv    => psb_d_nest_base_vect_mv
    ! full base-class contract (delegating to the blocks):
    procedure, pass(a) :: csmm       => psb_d_nest_base_csmm
    procedure, pass(a) :: cp_to_coo  => psb_d_nest_base_cp_to_coo
    procedure, pass(a) :: mv_to_coo  => psb_d_nest_base_mv_to_coo
    procedure, pass(a) :: rowsum     => psb_d_nest_base_rowsum
    procedure, pass(a) :: arwsum     => psb_d_nest_base_arwsum
    procedure, pass(a) :: colsum     => psb_d_nest_base_colsum
    procedure, pass(a) :: aclsum     => psb_d_nest_base_aclsum
    procedure, pass(a) :: maxval     => psb_d_nest_base_maxval
    procedure, pass(a) :: spnmi      => psb_d_nest_base_csnmi
    procedure, pass(a) :: spnm1      => psb_d_nest_base_csnm1
    procedure, pass(a) :: scals      => psb_d_nest_base_scals
    procedure, pass(a) :: scalv      => psb_d_nest_base_scal
    procedure, pass(a) :: clone      => psb_d_nest_base_clone
    procedure, pass(a) :: mold       => psb_d_nest_base_mold
    procedure, pass(a) :: sizeof     => psb_d_nest_base_sizeof
    ! NOT implemented on purpose (base error 700 is the intended behaviour):
    ! cp_from_coo / mv_from_coo (a nested operator cannot be built from a flat
    ! matrix without the field structure), csput (insertions go to the blocks
    ! before assembly), cssv/cssm (triangular solve is undefined for a block
    ! operator)
  end type psb_d_nest_base_mat

  private
  public :: psb_d_nest_base_mat, psb_d_nest_base_setup, psb_d_nest_apply_block
  ! field-split interface (for the block preconditioner)
  public :: psb_d_nest_get_n_fields, psb_d_nest_get_field_owned, &
       &    psb_d_nest_get_block, psb_d_nest_get_field_desc,      &
       &    psb_d_nest_restrict_field, psb_d_nest_restrict_field_local, &
       &    psb_d_nest_prolong_field

contains

  function psb_d_nest_base_get_fmt() result(format_name)
    character(len=5) :: format_name
    format_name = 'NEST'
  end function psb_d_nest_base_get_fmt

  ! free: the nested operator does NOT own block_storage / grid_desc (they are
  ! pointers into the caller), so we only detach them and release the field maps.
  subroutine psb_d_nest_base_free(a)
    class(psb_d_nest_base_mat), intent(inout) :: a
    integer(psb_ipk_) :: i_field, local_info
    if (allocated(a%field_map)) then
      do i_field = 1, size(a%field_map)
        if (associated(a%field_map(i_field)%gather_pos)) then
          call a%field_map(i_field)%gather_pos%free(local_info)
          deallocate(a%field_map(i_field)%gather_pos)
          a%field_map(i_field)%gather_pos => null()
        end if
      end do
      deallocate(a%field_map)
    end if
    a%block_storage => null()
    a%grid_desc     => null()
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

  ! get_diag: diagonal of the global operator.  In the global-local layout the
  ! owned entries of field i occupy positions owned_offset+1..owned_offset+n_owned,
  ! and for owned indices the field-local column k maps to the same global-local
  ! position as row k, so the global diagonal is the concatenation of the
  ! diagonals of the diagonal blocks (i,i); absent blocks contribute zeros
  ! (e.g. the (2,2) block of a saddle-point operator).
  subroutine psb_d_nest_base_get_diag(a, d, info)
    class(psb_d_nest_base_mat), intent(in) :: a
    real(psb_dpk_),             intent(out) :: d(:)
    integer(psb_ipk_),          intent(out) :: info

    real(psb_dpk_), allocatable :: block_diag(:)
    integer(psb_ipk_)           :: i_field, n_owned, owned_offset
    character(len=24)           :: name

    info = psb_success_
    name = 'psb_d_nest_get_diag'

    if (.not. (associated(a%block_storage) .and. allocated(a%field_map))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, name, a_err='nested operator not set up')
      return
    end if
    if (size(d) < a%get_nrows()) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='d too small')
      return
    end if

    d(1:a%get_nrows()) = dzero
    owned_offset = 0
    do i_field = 1, a%n_fields
      n_owned = a%field_map(i_field)%n_owned
      if (a%block_storage%has_block(i_field, i_field)) then
        allocate(block_diag(n_owned), stat=info)
        if (info /= 0) then
          info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
        end if
        call a%block_storage%mats(i_field,i_field)%a%get_diag(block_diag, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='block get_diag')
          return
        end if
        d(owned_offset+1 : owned_offset+n_owned) = block_diag(1:n_owned)
        deallocate(block_diag)
      end if
      owned_offset = owned_offset + n_owned
    end do
  end subroutine psb_d_nest_base_get_diag

  ! csgetrow: extract local rows imin..imax of the global operator as COO
  ! triplets, with columns in the global-local layout (the operator's column
  ! space).  Each global-local row r belongs to one field i (row k within the
  ! field); its entries are the union over j of row k of block (i,j), with the
  ! block-local column c remapped through field_map(j)%global_local_pos(c).
  ! This is the format-agnostic access path used by the ILU factorizations of
  ! the BJAC preconditioner (via csget/csgetblk).
  subroutine psb_d_nest_base_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
    class(psb_d_nest_base_mat), intent(in) :: a
    integer(psb_ipk_), intent(in)                  :: imin,imax
    integer(psb_ipk_), intent(out)                 :: nz
    integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
    real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer(psb_ipk_),intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer(psb_ipk_), intent(in), optional        :: iren(:)
    integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
    logical, intent(in), optional        :: rscale,cscale,chksz

    integer(psb_ipk_), allocatable :: block_row_ia(:), block_row_ja(:)
    real(psb_dpk_),    allocatable :: block_row_val(:)
    integer(psb_ipk_) :: jmin_, jmax_, nzin_, out_pos
    integer(psb_ipk_) :: r_row, i_field, j_field, k_in_field, owned_offset
    integer(psb_ipk_) :: block_nz, t_entry, global_local_col
    logical           :: append_
    character(len=24) :: name

    info = psb_success_
    name = 'psb_d_nest_csgetrow'

    if (.not. (associated(a%block_storage) .and. allocated(a%field_map))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, name, a_err='nested operator not set up')
      return
    end if
    if (present(iren)) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='iren not supported'); return
    end if
    if (present(rscale)) then
      if (rscale) then
        info = psb_err_invalid_input_
        call psb_errpush(info, name, a_err='rscale not supported'); return
      end if
    end if
    if (present(cscale)) then
      if (cscale) then
        info = psb_err_invalid_input_
        call psb_errpush(info, name, a_err='cscale not supported'); return
      end if
    end if

    jmin_ = 1
    jmax_ = a%get_ncols()
    if (present(jmin)) jmin_ = jmin
    if (present(jmax)) jmax_ = jmax
    append_ = .false.
    if (present(append)) append_ = append
    nzin_ = 0
    if (append_ .and. present(nzin)) nzin_ = nzin

    nz      = 0
    out_pos = nzin_

    do r_row = max(imin, 1), min(imax, a%get_nrows())
      ! locate the field owning global-local row r_row
      owned_offset = 0
      i_field      = 0
      do while (i_field < a%n_fields)
        i_field = i_field + 1
        if (r_row <= owned_offset + a%field_map(i_field)%n_owned) exit
        owned_offset = owned_offset + a%field_map(i_field)%n_owned
      end do
      k_in_field = r_row - owned_offset

      do j_field = 1, a%n_fields
        if (.not. a%block_storage%has_block(i_field, j_field)) cycle
        call a%block_storage%mats(i_field,j_field)%a%csgetrow(k_in_field, k_in_field, &
             & block_nz, block_row_ia, block_row_ja, block_row_val, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='block csgetrow')
          return
        end if
        do t_entry = 1, block_nz
          global_local_col = a%field_map(j_field)%global_local_pos(block_row_ja(t_entry))
          if ((global_local_col < jmin_) .or. (global_local_col > jmax_)) cycle
          out_pos = out_pos + 1
          call psb_ensure_size(out_pos, ia,  info)
          if (info == psb_success_) call psb_ensure_size(out_pos, ja,  info)
          if (info == psb_success_) call psb_ensure_size(out_pos, val, info)
          if (info /= psb_success_) then
            info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
          end if
          ia(out_pos)  = r_row
          ja(out_pos)  = global_local_col
          val(out_pos) = block_row_val(t_entry)
          nz = nz + 1
        end do
      end do
    end do
  end subroutine psb_d_nest_base_csgetrow

  ! vect_mv: matvec on encapsulated vectors (the path taken by psb_spmm with
  ! psb_d_vect_type).  Instead of falling back to the host-array csmv, it
  ! (1) gathers each column-field sub-vector through the vector's own gth with
  !     an encapsulated index vector (a device kernel on device vectors),
  ! (2) runs each block through its vect_mv (device formats execute their own
  !     device kernels), with per-field work vectors allocated with mold=x so
  !     they share the dynamic type of the incoming vectors,
  ! (3) scatters each row-field result back through the vector's own sct.
  ! Host/device traffic is limited to the compact field buffers; on plain host
  ! vectors this is exactly equivalent to the array csmv.
  subroutine psb_d_nest_base_vect_mv(alpha, a, x, beta, y, info, trans)
    class(psb_d_nest_base_mat), intent(in)     :: a
    real(psb_dpk_),             intent(in)     :: alpha, beta
    class(psb_d_base_vect_type), intent(inout) :: x
    class(psb_d_base_vect_type), intent(inout) :: y
    integer(psb_ipk_),          intent(out)    :: info
    character, optional,        intent(in)     :: trans

    class(psb_d_base_vect_type), allocatable :: x_field_vec, y_field_vec
    real(psb_dpk_), allocatable :: x_field_buf(:), y_field_buf(:)
    real(psb_dpk_)    :: block_beta
    integer(psb_ipk_) :: i_field, j_field, n_owned, n_local, local_info
    logical           :: row_has_blocks
    character         :: trans_
    character(len=24) :: name

    info = psb_success_
    name = 'psb_d_nest_vect_mv'

    trans_ = 'N'
    if (present(trans)) trans_ = trans
    if (.not. (associated(a%block_storage) .and. allocated(a%field_map))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, name, a_err='nested operator not set up')
      return
    end if
    if (trans_ /= 'N' .and. trans_ /= 'n') then
      ! transposed product: fall back to host arrays (rare path)
      block
        real(psb_dpk_), allocatable :: x_host(:), y_host(:)
        x_host = x%get_vect()
        y_host = y%get_vect()
        call psb_d_nest_base_csmv_t(alpha, a, x_host, beta, y_host, info, trans_)
        call y%bld(y_host)
      end block
      return
    end if

    ! work vectors share the dynamic type of the incoming vectors
    allocate(x_field_vec, mold=x, stat=info)
    if (info == 0) allocate(y_field_vec, mold=y, stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
    end if

    do i_field = 1, a%n_fields
      n_owned = a%field_map(i_field)%n_owned
      call psb_ensure_size(n_owned, y_field_buf, info)
      if (info /= psb_success_) then
        info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
      end if

      row_has_blocks = .false.
      block_beta     = dzero
      do j_field = 1, a%n_fields
        if (.not. a%block_storage%has_block(i_field, j_field)) cycle

        ! gather the column-field sub-vector (owned + ghosts) from x
        n_local = size(a%field_map(j_field)%global_local_pos)
        call psb_ensure_size(n_local, x_field_buf, info)
        if (info /= psb_success_) then
          info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
        end if
        call x%gth(ione, int(n_local, psb_mpk_), &
             &     a%field_map(j_field)%gather_pos%v, x_field_buf)
        call x_field_vec%free(local_info)
        call x_field_vec%bld(x_field_buf(1:n_local))

        if (.not. row_has_blocks) then
          ! first block of this row field: (re)build the accumulator at the
          ! right size, zeroed
          y_field_buf(1:n_owned) = dzero
          call y_field_vec%free(local_info)
          call y_field_vec%bld(y_field_buf(1:n_owned))
          row_has_blocks = .true.
        end if

        ! y_field = alpha * A(i,j) * x_field + block_beta * y_field
        call a%block_storage%mats(i_field,j_field)%a%spmm(alpha, x_field_vec, &
             &                                            block_beta, y_field_vec, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='block vect_mv')
          return
        end if
        block_beta = done
      end do

      ! scatter the row-field result into y (beta applied on the owned rows);
      ! a row field with no blocks still rescales its rows by beta
      if (row_has_blocks) then
        y_field_buf(1:n_owned) = y_field_vec%get_vect()
      else
        y_field_buf(1:n_owned) = dzero
      end if
      call y%sct(ione, int(n_owned, psb_mpk_), &
           &     a%field_map(i_field)%gather_pos%v, y_field_buf, beta)
    end do

    call x_field_vec%free(local_info)
    call y_field_vec%free(local_info)
  end subroutine psb_d_nest_base_vect_mv

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
      ! encapsulated copy of the positions for the device-capable gth/sct
      allocate(nest_op%field_map(i_field)%gather_pos, stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
      end if
      call nest_op%field_map(i_field)%gather_pos%bld(nest_op%field_map(i_field)%global_local_pos)
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
      ! transposed product: the block structure of A^T is the transpose of the
      ! block grid, handled by the dedicated kernel below ('T' or 'C')
      call psb_d_nest_base_csmv_t(alpha, a, x, beta, y, info, trans_op)
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

  ! Transposed matvec kernel: y = alpha * A^T * x + beta * y.
  ! The block structure of A^T is the transpose of the block grid:
  !   y(cols of field j) += alpha * sum_i A(i,j)^T * x(owned rows of field i).
  ! x is read on the owned rows of each row field; the result lands on ALL the
  ! local columns of each column field (owned + ghosts); the distributed caller
  ! (psb_spmm with trans='T') then accumulates the ghost contributions to their
  ! owners through the transposed halo exchange.
  subroutine psb_d_nest_base_csmv_t(alpha, a, x, beta, y, info, trans)
    real(psb_dpk_),             intent(in)    :: alpha, beta, x(:)
    class(psb_d_nest_base_mat), intent(in)    :: a
    real(psb_dpk_),             intent(inout) :: y(:)
    integer(psb_ipk_),          intent(out)   :: info
    character,                  intent(in)    :: trans

    real(psb_dpk_), allocatable :: x_field(:), y_field(:)
    integer(psb_ipk_)           :: i_block_row, j_block_col, i_entry
    integer(psb_ipk_)           :: n_local_col_field, n_owned_row_field
    character(len=24)           :: name

    info = psb_success_
    name = 'psb_d_nest_base_csmv_t'

    if (.not. associated(a%block_storage)) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='nested operator not set up')
      return
    end if

    ! y <- beta * y  (on the whole column space)
    if (beta == dzero) then
      y(:) = dzero
    else if (beta /= done) then
      y(:) = beta * y(:)
    end if

    do j_block_col = 1, a%n_fields
      n_local_col_field = size(a%field_map(j_block_col)%global_local_pos)
      allocate(y_field(n_local_col_field), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
      end if
      ! current column-field output sub-vector (owned + ghosts)
      do i_entry = 1, n_local_col_field
        y_field(i_entry) = y(a%field_map(j_block_col)%global_local_pos(i_entry))
      end do

      do i_block_row = 1, a%n_fields
        if (a%block_storage%has_block(i_block_row, j_block_col)) then
          n_owned_row_field = a%field_map(i_block_row)%n_owned
          allocate(x_field(n_owned_row_field), stat=info)
          if (info /= 0) then
            info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
          end if
          ! gather the row-field input sub-vector (owned rows only)
          do i_entry = 1, n_owned_row_field
            x_field(i_entry) = x(a%field_map(i_block_row)%global_local_pos(i_entry))
          end do
          ! y_field <- alpha * A(i,j)^T (or ^H) * x_field + y_field
          call a%block_storage%mats(i_block_row, j_block_col)%a%csmv( &
               & alpha, x_field, done, y_field, info, trans)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_from_subroutine_, name, a_err='block csmv T')
            return
          end if
          deallocate(x_field)
        end if
      end do

      ! scatter the column-field output sub-vector back into y
      do i_entry = 1, n_local_col_field
        y(a%field_map(j_block_col)%global_local_pos(i_entry)) = y_field(i_entry)
      end do
      deallocate(y_field)
    end do
  end subroutine psb_d_nest_base_csmv_t

  ! csmm: multi-RHS product, the 2D analogue of csmv (same gather/scatter
  ! per field, the block product is the block's own csmm)
  subroutine psb_d_nest_base_csmm(alpha, a, x, beta, y, info, trans)
    class(psb_d_nest_base_mat), intent(in)    :: a
    real(psb_dpk_),             intent(in)    :: alpha, beta, x(:,:)
    real(psb_dpk_),             intent(inout) :: y(:,:)
    integer(psb_ipk_),          intent(out)   :: info
    character, optional,        intent(in)    :: trans

    real(psb_dpk_), allocatable :: x_field(:,:), y_field(:,:)
    integer(psb_ipk_)           :: i_block_row, j_block_col, i_entry
    integer(psb_ipk_)           :: n_local_col_field, n_owned_row_field, n_rhs
    character                   :: trans_op
    character(len=24)           :: name

    info = psb_success_
    name = 'psb_d_nest_base_csmm'
    trans_op = 'N'
    if (present(trans)) trans_op = trans
    if (trans_op /= 'N' .and. trans_op /= 'n') then
      info = psb_err_transpose_not_n_unsupported_
      call psb_errpush(info, name); return
    end if
    if (.not. associated(a%block_storage)) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='nested operator not set up')
      return
    end if
    n_rhs = min(size(x,2), size(y,2))

    if (beta == dzero) then
      y(:,:) = dzero
    else if (beta /= done) then
      y(:,:) = beta * y(:,:)
    end if

    do j_block_col = 1, a%n_fields
      n_local_col_field = size(a%field_map(j_block_col)%global_local_pos)
      allocate(x_field(n_local_col_field, n_rhs), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
      end if
      do i_entry = 1, n_local_col_field
        x_field(i_entry, 1:n_rhs) = x(a%field_map(j_block_col)%global_local_pos(i_entry), 1:n_rhs)
      end do

      do i_block_row = 1, a%n_fields
        if (a%block_storage%has_block(i_block_row, j_block_col)) then
          n_owned_row_field = a%field_map(i_block_row)%n_owned
          allocate(y_field(n_owned_row_field, n_rhs), stat=info)
          if (info /= 0) then
            info = psb_err_alloc_dealloc_; call psb_errpush(info, name); return
          end if
          do i_entry = 1, n_owned_row_field
            y_field(i_entry, 1:n_rhs) = y(a%field_map(i_block_row)%global_local_pos(i_entry), 1:n_rhs)
          end do
          call a%block_storage%mats(i_block_row, j_block_col)%a%csmm( &
               & alpha, x_field, done, y_field, info, trans_op)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_from_subroutine_, name, a_err='block csmm')
            return
          end if
          do i_entry = 1, n_owned_row_field
            y(a%field_map(i_block_row)%global_local_pos(i_entry), 1:n_rhs) = y_field(i_entry, 1:n_rhs)
          end do
          deallocate(y_field)
        end if
      end do
      deallocate(x_field)
    end do
  end subroutine psb_d_nest_base_csmm

  ! cp_to_coo: assemble all the blocks into a single local COO in the
  ! global-local layout (rows = concatenated owned rows, columns = the
  ! operator's column space).  This is the core conversion hook: the generic
  ! base-class machinery builds cscnv, csclip, tril/triu, ... on top of it.
  subroutine psb_d_nest_base_cp_to_coo(a, b, info)
    use psb_d_base_mat_mod, only : psb_d_coo_sparse_mat
    class(psb_d_nest_base_mat),  intent(in)    :: a
    class(psb_d_coo_sparse_mat), intent(inout) :: b
    integer(psb_ipk_),           intent(out)   :: info

    type(psb_d_coo_sparse_mat) :: block_coo
    integer(psb_ipk_) :: i_field, j_field, k_entry, n_entries, out_pos, owned_offset
    character(len=24) :: name

    info = psb_success_
    name = 'psb_d_nest_cp_to_coo'

    if (.not. (associated(a%block_storage) .and. allocated(a%field_map))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, name, a_err='nested operator not set up')
      return
    end if

    call b%allocate(a%get_nrows(), a%get_ncols(), a%get_nzeros())
    out_pos      = 0
    owned_offset = 0
    do i_field = 1, a%n_fields
      do j_field = 1, a%n_fields
        if (.not. a%block_storage%has_block(i_field, j_field)) cycle
        call a%block_storage%mats(i_field,j_field)%a%cp_to_coo(block_coo, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='block cp_to_coo')
          return
        end if
        n_entries = block_coo%get_nzeros()
        do k_entry = 1, n_entries
          b%ia(out_pos+k_entry)  = owned_offset + block_coo%ia(k_entry)
          b%ja(out_pos+k_entry)  = a%field_map(j_field)%global_local_pos(block_coo%ja(k_entry))
          b%val(out_pos+k_entry) = block_coo%val(k_entry)
        end do
        out_pos = out_pos + n_entries
        call block_coo%free()
      end do
      owned_offset = owned_offset + a%field_map(i_field)%n_owned
    end do
    call b%set_nzeros(out_pos)
    call b%set_dupl(psb_dupl_add_)
    call b%fix(info)
    if (info /= psb_success_) &
         & call psb_errpush(psb_err_from_subroutine_, name, a_err='coo fix')
  end subroutine psb_d_nest_base_cp_to_coo

  ! mv_to_coo: the adapter does not own the blocks, so "move" degenerates to
  ! copy + detach of the adapter
  subroutine psb_d_nest_base_mv_to_coo(a, b, info)
    use psb_d_base_mat_mod, only : psb_d_coo_sparse_mat
    class(psb_d_nest_base_mat),  intent(inout) :: a
    class(psb_d_coo_sparse_mat), intent(inout) :: b
    integer(psb_ipk_),           intent(out)   :: info

    call a%cp_to_coo(b, info)
    if (info == psb_success_) call a%free()
  end subroutine psb_d_nest_base_mv_to_coo

  ! rowsum: row sums (matrix-valued type), accumulated across the blocks of
  ! each row field; d is in the global-local row layout
  subroutine psb_d_nest_base_rowsum(d, a)
    class(psb_d_nest_base_mat), intent(in)  :: a
    real(psb_dpk_),             intent(out) :: d(:)

    real(psb_dpk_), allocatable :: block_sums(:)
    integer(psb_ipk_) :: i_field, j_field, k_entry, n_owned, owned_offset

    d(:) = dzero
    if (.not. associated(a%block_storage)) return
    owned_offset = 0
    do i_field = 1, a%n_fields
      n_owned = a%field_map(i_field)%n_owned
      allocate(block_sums(n_owned))
      do j_field = 1, a%n_fields
        if (.not. a%block_storage%has_block(i_field, j_field)) cycle
        call a%block_storage%mats(i_field,j_field)%a%rowsum(block_sums)
        do k_entry = 1, n_owned
          d(owned_offset+k_entry) = d(owned_offset+k_entry) + block_sums(k_entry)
        end do
      end do
      deallocate(block_sums)
      owned_offset = owned_offset + n_owned
    end do
  end subroutine psb_d_nest_base_rowsum

  ! arwsum: absolute row sums (always real-valued)
  subroutine psb_d_nest_base_arwsum(d, a)
    class(psb_d_nest_base_mat), intent(in)  :: a
    real(psb_dpk_),             intent(out) :: d(:)

    real(psb_dpk_), allocatable :: block_sums(:)
    integer(psb_ipk_) :: i_field, j_field, k_entry, n_owned, owned_offset

    d(:) = dzero
    if (.not. associated(a%block_storage)) return
    owned_offset = 0
    do i_field = 1, a%n_fields
      n_owned = a%field_map(i_field)%n_owned
      allocate(block_sums(n_owned))
      do j_field = 1, a%n_fields
        if (.not. a%block_storage%has_block(i_field, j_field)) cycle
        call a%block_storage%mats(i_field,j_field)%a%arwsum(block_sums)
        do k_entry = 1, n_owned
          d(owned_offset+k_entry) = d(owned_offset+k_entry) + block_sums(k_entry)
        end do
      end do
      deallocate(block_sums)
      owned_offset = owned_offset + n_owned
    end do
  end subroutine psb_d_nest_base_arwsum

  ! colsum: column sums (matrix-valued type) in the operator's column space,
  ! accumulated across the blocks of each column field
  subroutine psb_d_nest_base_colsum(d, a)
    class(psb_d_nest_base_mat), intent(in)  :: a
    real(psb_dpk_),             intent(out) :: d(:)

    real(psb_dpk_), allocatable :: field_sums(:), block_sums(:)
    integer(psb_ipk_) :: i_field, j_field, k_entry, n_local

    d(:) = dzero
    if (.not. associated(a%block_storage)) return
    do j_field = 1, a%n_fields
      n_local = size(a%field_map(j_field)%global_local_pos)
      allocate(field_sums(n_local), block_sums(n_local))
      field_sums(:) = dzero
      do i_field = 1, a%n_fields
        if (.not. a%block_storage%has_block(i_field, j_field)) cycle
        call a%block_storage%mats(i_field,j_field)%a%colsum(block_sums)
        field_sums(1:n_local) = field_sums(1:n_local) + block_sums(1:n_local)
      end do
      do k_entry = 1, n_local
        d(a%field_map(j_field)%global_local_pos(k_entry)) = field_sums(k_entry)
      end do
      deallocate(field_sums, block_sums)
    end do
  end subroutine psb_d_nest_base_colsum

  ! aclsum: absolute column sums (always real-valued)
  subroutine psb_d_nest_base_aclsum(d, a)
    class(psb_d_nest_base_mat), intent(in)  :: a
    real(psb_dpk_),             intent(out) :: d(:)

    real(psb_dpk_), allocatable :: field_sums(:), block_sums(:)
    integer(psb_ipk_) :: i_field, j_field, k_entry, n_local

    d(:) = dzero
    if (.not. associated(a%block_storage)) return
    do j_field = 1, a%n_fields
      n_local = size(a%field_map(j_field)%global_local_pos)
      allocate(field_sums(n_local), block_sums(n_local))
      field_sums(:) = dzero
      do i_field = 1, a%n_fields
        if (.not. a%block_storage%has_block(i_field, j_field)) cycle
        call a%block_storage%mats(i_field,j_field)%a%aclsum(block_sums)
        field_sums(1:n_local) = field_sums(1:n_local) + block_sums(1:n_local)
      end do
      do k_entry = 1, n_local
        d(a%field_map(j_field)%global_local_pos(k_entry)) = field_sums(k_entry)
      end do
      deallocate(field_sums, block_sums)
    end do
  end subroutine psb_d_nest_base_aclsum

  ! maxval / infinity norm / 1-norm, by delegation/accumulation over blocks
  function psb_d_nest_base_maxval(a) result(res)
    class(psb_d_nest_base_mat), intent(in) :: a
    real(psb_dpk_) :: res
    integer(psb_ipk_) :: i_field, j_field
    res = dzero
    if (.not. associated(a%block_storage)) return
    do j_field = 1, a%n_fields
      do i_field = 1, a%n_fields
        if (a%block_storage%has_block(i_field, j_field)) &
             & res = max(res, a%block_storage%mats(i_field,j_field)%a%maxval())
      end do
    end do
  end function psb_d_nest_base_maxval

  function psb_d_nest_base_csnmi(a) result(res)
    class(psb_d_nest_base_mat), intent(in) :: a
    real(psb_dpk_) :: res
    real(psb_dpk_), allocatable :: row_sums(:)
    res = dzero
    if (a%get_nrows() <= 0) return
    allocate(row_sums(a%get_nrows()))
    call psb_d_nest_base_arwsum(row_sums, a)
    res = maxval(row_sums)
  end function psb_d_nest_base_csnmi

  function psb_d_nest_base_csnm1(a) result(res)
    class(psb_d_nest_base_mat), intent(in) :: a
    real(psb_dpk_) :: res
    real(psb_dpk_), allocatable :: col_sums(:)
    res = dzero
    if (a%get_ncols() <= 0) return
    allocate(col_sums(a%get_ncols()))
    call psb_d_nest_base_aclsum(col_sums, a)
    res = maxval(col_sums)
  end function psb_d_nest_base_csnm1

  ! scals/scal: scaling acts on the underlying blocks (the operator is a view)
  subroutine psb_d_nest_base_scals(d, a, info)
    class(psb_d_nest_base_mat), intent(inout) :: a
    real(psb_dpk_),             intent(in)    :: d
    integer(psb_ipk_),          intent(out)   :: info
    integer(psb_ipk_) :: i_field, j_field
    character(len=24) :: name
    info = psb_success_
    name = 'psb_d_nest_scals'
    if (.not. associated(a%block_storage)) then
      info = psb_err_invalid_mat_state_; call psb_errpush(info, name); return
    end if
    do j_field = 1, a%n_fields
      do i_field = 1, a%n_fields
        if (.not. a%block_storage%has_block(i_field, j_field)) cycle
        call a%block_storage%mats(i_field,j_field)%a%scal(d, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='block scal'); return
        end if
      end do
    end do
  end subroutine psb_d_nest_base_scals

  subroutine psb_d_nest_base_scal(d, a, info, side)
    class(psb_d_nest_base_mat), intent(inout) :: a
    real(psb_dpk_),             intent(in)    :: d(:)
    integer(psb_ipk_),          intent(out)   :: info
    character, intent(in), optional           :: side

    real(psb_dpk_), allocatable :: d_field(:)
    integer(psb_ipk_) :: i_field, j_field, k_entry, n_owned, n_local, owned_offset
    character          :: side_
    character(len=24)  :: name

    info  = psb_success_
    name  = 'psb_d_nest_scal'
    side_ = 'L'
    if (present(side)) side_ = side
    if (.not. associated(a%block_storage)) then
      info = psb_err_invalid_mat_state_; call psb_errpush(info, name); return
    end if

    if (side_ == 'L' .or. side_ == 'l') then
      ! row scaling: each row field uses its owned slice of d
      owned_offset = 0
      do i_field = 1, a%n_fields
        n_owned = a%field_map(i_field)%n_owned
        do j_field = 1, a%n_fields
          if (.not. a%block_storage%has_block(i_field, j_field)) cycle
          call a%block_storage%mats(i_field,j_field)%a%scal( &
               & d(owned_offset+1:owned_offset+n_owned), info, side='L')
          if (info /= psb_success_) then
            call psb_errpush(psb_err_from_subroutine_, name, a_err='block scal L'); return
          end if
        end do
        owned_offset = owned_offset + n_owned
      end do
    else
      ! column scaling: each column field gathers its slice of d
      do j_field = 1, a%n_fields
        n_local = size(a%field_map(j_field)%global_local_pos)
        allocate(d_field(n_local))
        do k_entry = 1, n_local
          d_field(k_entry) = d(a%field_map(j_field)%global_local_pos(k_entry))
        end do
        do i_field = 1, a%n_fields
          if (.not. a%block_storage%has_block(i_field, j_field)) cycle
          call a%block_storage%mats(i_field,j_field)%a%scal(d_field, info, side='R')
          if (info /= psb_success_) then
            call psb_errpush(psb_err_from_subroutine_, name, a_err='block scal R'); return
          end if
        end do
        deallocate(d_field)
      end do
    end if
  end subroutine psb_d_nest_base_scal

  ! clone: the adapter is a view, so the clone shares the blocks and the grid
  ! descriptor (pointers) while re-owning its private gather index vectors
  subroutine psb_d_nest_base_clone(a, b, info)
    class(psb_d_nest_base_mat),                intent(inout) :: a
    class(psb_d_base_sparse_mat), allocatable, intent(inout) :: b
    integer(psb_ipk_),                         intent(out)   :: info
    integer(psb_ipk_) :: i_field

    info = psb_success_
    if (allocated(b)) deallocate(b)
    allocate(b, source=a, stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, 'psb_d_nest_clone'); return
    end if
    select type (b_nest => b)
    type is (psb_d_nest_base_mat)
      if (allocated(b_nest%field_map)) then
        do i_field = 1, size(b_nest%field_map)
          ! the sourced copy shares a's gather_pos targets: re-own fresh copies
          b_nest%field_map(i_field)%gather_pos => null()
          allocate(b_nest%field_map(i_field)%gather_pos, stat=info)
          if (info /= 0) then
            info = psb_err_alloc_dealloc_; call psb_errpush(info, 'psb_d_nest_clone'); return
          end if
          call b_nest%field_map(i_field)%gather_pos%bld( &
               & b_nest%field_map(i_field)%global_local_pos)
        end do
      end if
    end select
  end subroutine psb_d_nest_base_clone

  subroutine psb_d_nest_base_mold(a, b, info)
    class(psb_d_nest_base_mat),                intent(in)    :: a
    class(psb_d_base_sparse_mat), allocatable, intent(inout) :: b
    integer(psb_ipk_),                         intent(out)   :: info
    info = psb_success_
    if (allocated(b)) deallocate(b)
    allocate(b, mold=a, stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; call psb_errpush(info, 'psb_d_nest_mold')
    end if
  end subroutine psb_d_nest_base_mold

  ! sizeof: blocks + gather maps (the adapter does not own the descriptors)
  function psb_d_nest_base_sizeof(a) result(res)
    class(psb_d_nest_base_mat), intent(in) :: a
    integer(psb_epk_) :: res
    integer(psb_ipk_) :: i_field
    res = 8
    if (associated(a%block_storage)) res = res + a%block_storage%sizeof()
    if (allocated(a%field_map)) then
      do i_field = 1, size(a%field_map)
        if (allocated(a%field_map(i_field)%global_local_pos)) &
             & res = res + psb_sizeof_ip * size(a%field_map(i_field)%global_local_pos)
      end do
    end if
  end function psb_d_nest_base_sizeof

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

  ! Restrict: extract field k's full local sub-vector (owned + ghosts).
  subroutine psb_d_nest_restrict_field_local(nest_op, field, x_global, x_field, info)
    type(psb_d_nest_base_mat), intent(in)  :: nest_op
    integer(psb_ipk_),         intent(in)  :: field
    real(psb_dpk_),            intent(in)  :: x_global(:)
    real(psb_dpk_),            intent(out) :: x_field(:)
    integer(psb_ipk_),         intent(out) :: info
    integer(psb_ipk_) :: i_entry, n_local
    info = psb_success_
    if (field < 1 .or. field > nest_op%n_fields) then
      info = psb_err_invalid_input_; return
    end if
    n_local = size(nest_op%field_map(field)%global_local_pos)
    do i_entry = 1, n_local
      x_field(i_entry) = x_global(nest_op%field_map(field)%global_local_pos(i_entry))
    end do
  end subroutine psb_d_nest_restrict_field_local

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
