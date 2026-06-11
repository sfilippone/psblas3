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
!         software without specific written permission.
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
! Module: psb_d_nest_tools_mod
! Author: Simone Staccone (Stack-1)
!
!  Nested-specific assembly wrappers for PSBLAS3 — double precision matrix and vector routines
!

module psb_d_nest_tools_mod
  use psb_const_mod, only : psb_ipk_, psb_lpk_, psb_dpk_, psb_success_, psb_err_alloc_dealloc_, &
                             psb_err_invalid_input_, psb_err_from_subroutine_, &
                             psb_dupl_add_, psb_dupl_ovwrt_, psb_dupl_err_, psb_ctxt_type
  use psb_error_mod, only : psb_errpush
  use psb_d_tools_mod, only : psb_spall, psb_spins, psb_spasb, psb_spfree, psb_sprn, &
                               psb_geall, psb_geins, psb_geasb, psb_gefree
  use psb_desc_nest_mod, only : psb_desc_nest_type
  use psb_d_nest_mat_mod, only : psb_d_nest_sparse_mat
  use psb_d_mat_mod,      only : psb_dspmat_type
  use psb_d_base_mat_mod, only : psb_d_coo_sparse_mat
  use psb_desc_mod,       only : psb_desc_type
  implicit none

  private

  public :: psb_spall_nest, psb_spins_nest, psb_spasb_nest, psb_spfree_nest, psb_sprn_nest, &
            psb_d_nest_rect_block

contains

  ! Allocates all (nrblocks x ncblocks) sparse matrix blocks
  ! and marks all as present.  psb_spins_nest lazy-allocates individual
  ! blocks on first insertion; call psb_spall_nest instead when the
  ! full block structure is known up front.

  subroutine psb_spall_nest(a_nest, desc_nest, info, nnz)
    type(psb_d_nest_sparse_mat), intent(inout) :: a_nest
    type(psb_desc_nest_type),    intent(in)    :: desc_nest
    integer(psb_ipk_),           intent(out)   :: info
    integer(psb_ipk_),           intent(in), optional :: nnz

    integer(psb_ipk_) :: i_block_row, j_block_col, local_info
    character(len=20) :: name

    info = psb_success_
    name = 'psb_spall_nest'

    a_nest%nrblocks = desc_nest%nrblocks
    a_nest%ncblocks = desc_nest%ncblocks

    if (.not. allocated(a_nest%mats)) then
      allocate(a_nest%mats(a_nest%nrblocks, a_nest%ncblocks), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, name)
        return
      end if
    end if

    do i_block_row = 1, a_nest%nrblocks
      do j_block_col = 1, a_nest%ncblocks
        local_info = psb_success_
        if (present(nnz)) then
          call psb_spall(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), local_info, nnz=nnz)
        else
          call psb_spall(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), local_info)
        end if
        if (local_info /= psb_success_) then
          info = local_info
          call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spall')
          return
        end if
      end do
    end do

  end subroutine psb_spall_nest


  ! Inserts nz entries into block (blk_i, blk_j) of the nested matrix.
  ! The block is lazy-allocated on first insertion if psb_spall_nest
  ! was not called first.
  
  subroutine psb_spins_nest(block_row, block_col, n_entries, entry_rows, entry_cols, entry_vals, a_nest, desc_nest, info)
    integer(psb_ipk_),           intent(in)    :: block_row, block_col, n_entries
    integer(psb_lpk_),           intent(in)    :: entry_rows(:), entry_cols(:)
    real(psb_dpk_),              intent(in)    :: entry_vals(:)
    type(psb_d_nest_sparse_mat), intent(inout) :: a_nest
    type(psb_desc_nest_type),    intent(inout) :: desc_nest
    integer(psb_ipk_),           intent(out)   :: info

    integer(psb_ipk_) :: nnz_estimate
    character(len=20) :: name

    info = psb_success_
    name = 'psb_spins_nest'

    if (n_entries == 0) return

    if (block_row < 1 .or. block_row > a_nest%nrblocks .or. &
        block_col < 1 .or. block_col > a_nest%ncblocks) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='invalid block indices')
      return
    end if

    if (.not. allocated(a_nest%mats)) then
      allocate(a_nest%mats(a_nest%nrblocks, a_nest%ncblocks), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, name)
        return
      end if
    end if

    if (.not. allocated(a_nest%mats(block_row, block_col)%a)) then
      ! Estimate nnz: use n_entries + 50% buffer for future insertions
      nnz_estimate = max(n_entries, 10) + n_entries / 2
      call psb_spall(a_nest%mats(block_row, block_col), &
                     desc_nest%descs(block_row, block_col), info, nnz=nnz_estimate)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spall')
        return
      end if
    end if

    call psb_spins(n_entries, entry_rows, entry_cols, entry_vals, a_nest%mats(block_row, block_col), &
                   desc_nest%descs(block_row, block_col), info)
    if (info /= psb_success_) &
      call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spins')

  end subroutine psb_spins_nest

  ! Calls psb_spasb on all present block matrices.
  ! Must be called after psb_cdasb_nest.

  subroutine psb_spasb_nest(a_nest, desc_nest, info, dupl)
    type(psb_d_nest_sparse_mat), intent(inout) :: a_nest
    type(psb_desc_nest_type),    intent(inout) :: desc_nest
    integer(psb_ipk_),           intent(out)   :: info
    integer(psb_ipk_),           intent(in), optional :: dupl

    integer(psb_ipk_) :: i_block_row, j_block_col, dupl_mode, local_info
    character(len=20) :: name

    info      = psb_success_
    name      = 'psb_spasb_nest'
    dupl_mode = psb_dupl_add_
    if (present(dupl)) dupl_mode = dupl

    do i_block_row = 1, a_nest%nrblocks
      do j_block_col = 1, a_nest%ncblocks
        if (allocated(a_nest%mats(i_block_row, j_block_col)%a)) then
          local_info = psb_success_
          if (dupl_mode == psb_dupl_add_) then
            call psb_spasb(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), &
                           local_info, dupl=psb_dupl_add_)
          else if (dupl_mode == psb_dupl_ovwrt_) then
            call psb_spasb(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), &
                           local_info, dupl=psb_dupl_ovwrt_)
          else if (dupl_mode == psb_dupl_err_) then
            call psb_spasb(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), &
                           local_info, dupl=psb_dupl_err_)
          else
            call psb_spasb(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), &
                           local_info)
          end if
          if (local_info /= psb_success_) then
            info = local_info
            call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spasb')
            return
          end if
        end if
      end do
    end do

  end subroutine psb_spasb_nest

  ! Calls psb_spfree on every present block, then deallocates the
  ! mats array and resets nrblocks/ncblocks to 0.
  
  subroutine psb_spfree_nest(a_nest, desc_nest, info)
    type(psb_d_nest_sparse_mat), intent(inout) :: a_nest
    type(psb_desc_nest_type),    intent(in)    :: desc_nest
    integer(psb_ipk_),           intent(out)   :: info

    integer(psb_ipk_) :: i_block_row, j_block_col, local_info
    character(len=20) :: name

    info = psb_success_
    name = 'psb_spfree_nest'

    if (allocated(a_nest%mats)) then
      do i_block_row = 1, a_nest%nrblocks
        do j_block_col = 1, a_nest%ncblocks
          if (allocated(a_nest%mats(i_block_row, j_block_col)%a)) then
            local_info = psb_success_
            call psb_spfree(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), local_info)
            if (local_info /= psb_success_ .and. info == psb_success_) then
              info = local_info
              call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spfree')
            end if
          end if
        end do
      end do
      deallocate(a_nest%mats, stat=local_info)
      if (local_info /= 0 .and. info == psb_success_) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, name)
      end if
    end if

    a_nest%nrblocks = 0
    a_nest%ncblocks = 0

  end subroutine psb_spfree_nest

  ! Calls psb_sprn on every present block matrix, resetting it to
  ! the build state while preserving the sparsity pattern.
  
  subroutine psb_sprn_nest(a_nest, desc_nest, info, clear)
    type(psb_d_nest_sparse_mat), intent(inout) :: a_nest
    type(psb_desc_nest_type),    intent(in)    :: desc_nest
    integer(psb_ipk_),           intent(out)   :: info
    logical,                     intent(in), optional :: clear

    integer(psb_ipk_) :: i_block_row, j_block_col, local_info
    character(len=20) :: name

    info = psb_success_
    name = 'psb_sprn_nest'

    if (.not. allocated(a_nest%mats)) return

    do i_block_row = 1, a_nest%nrblocks
      do j_block_col = 1, a_nest%ncblocks
        if (allocated(a_nest%mats(i_block_row, j_block_col)%a)) then
          local_info = psb_success_
          if (present(clear)) then
            call psb_sprn(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), local_info, clear=clear)
          else
            call psb_sprn(a_nest%mats(i_block_row, j_block_col), desc_nest%descs(i_block_row, j_block_col), local_info)
          end if
          if (local_info /= psb_success_ .and. info == psb_success_) then
            info = local_info
            call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_sprn')
          end if
        end if
      end do
    end do

  end subroutine psb_sprn_nest



  ! psb_d_nest_rect_block
  !
  ! Build a local GENERAL (possibly rectangular) block A(i,j) of a nested
  ! operator, with rows in field i and columns in field j (field i /= field j,
  ! |field i| /= |field j| allowed). Rows are localized against the field-i
  ! (row) descriptor, columns against the field-j (column) descriptor — which
  ! must already carry the union halo of column j (cdall + cdins(all column-j
  ! blocks' columns) + cdasb). The result is a CSR block of shape
  !   (field-i owned rows) x (field-j local cols incl. halo)
  ! consumable directly by the nested csmv (psb_d_nest_base_mat).
  !
  ! A single-descriptor psb_spall/psb_spasb cannot express row-field /= col-field
  ! (it would force rows and columns into the same index space), hence the
  ! explicit COO build with separate row/column localization.
  !
  ! Arguments (this process's local contribution):
  !   blk            (out) the assembled block (CSR)
  !   nz             number of local entries
  !   ia_glob(:)     GLOBAL field-i row indices (owned by this process)
  !   ja_glob(:)     GLOBAL field-j column indices
  !   val(:)         values
  !   desc_row       field-i descriptor (rows)
  !   desc_col       field-j descriptor (columns, with union halo)
  !
  subroutine psb_d_nest_rect_block(blk, nz, ia_glob, ja_glob, val, desc_row, desc_col, info)
    type(psb_dspmat_type), intent(out) :: blk
    integer(psb_ipk_),     intent(in)  :: nz
    integer(psb_lpk_),     intent(in)  :: ia_glob(:), ja_glob(:)
    real(psb_dpk_),        intent(in)  :: val(:)
    type(psb_desc_type),   intent(in)  :: desc_row, desc_col
    integer(psb_ipk_),     intent(out) :: info

    type(psb_d_coo_sparse_mat) :: coo_block
    integer(psb_ipk_)          :: k_entry, n_loc_rows, n_loc_cols, loc_row, loc_col
    character(len=24)          :: name

    info = psb_success_
    name = 'psb_d_nest_rect_block'

    n_loc_rows = desc_row%get_local_rows()    ! owned rows of field i
    n_loc_cols = desc_col%get_local_cols()    ! field-j local cols (owned + halo)

    call coo_block%allocate(n_loc_rows, n_loc_cols, nz)
    do k_entry = 1, nz
      call desc_row%g2l(ia_glob(k_entry), loc_row, info)
      if (info /= 0 .or. loc_row < 1 .or. loc_row > n_loc_rows) then
        info = psb_err_invalid_input_
        call psb_errpush(info, name, a_err='row not owned / not localizable')
        return
      end if
      call desc_col%g2l(ja_glob(k_entry), loc_col, info)
      if (info /= 0 .or. loc_col < 1 .or. loc_col > n_loc_cols) then
        info = psb_err_invalid_input_
        call psb_errpush(info, name, a_err='column not in field-j descriptor (missing from union halo)')
        return
      end if
      coo_block%ia(k_entry)  = loc_row
      coo_block%ja(k_entry)  = loc_col
      coo_block%val(k_entry) = val(k_entry)
    end do
    call coo_block%set_nzeros(nz)
    call coo_block%set_dupl(psb_dupl_add_)
    call coo_block%fix(info)
    if (info /= 0) then
      call psb_errpush(psb_err_from_subroutine_, name, a_err='coo fix'); return
    end if
    call blk%mv_from(coo_block)
    call blk%cscnv(info, type='CSR')
    if (info /= 0) then
      call psb_errpush(psb_err_from_subroutine_, name, a_err='cscnv'); return
    end if
  end subroutine psb_d_nest_rect_block

end module psb_d_nest_tools_mod
