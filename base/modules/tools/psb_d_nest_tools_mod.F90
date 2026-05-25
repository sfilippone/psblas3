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
  use psb_d_nest_vect_mod, only : psb_d_nest_vect_type
  implicit none

  private

  public :: psb_spall_nest, psb_spins_nest, psb_spasb_nest, psb_spfree_nest, psb_sprn_nest, &
            psb_geall_nest, psb_geins_nest, psb_geasb_nest, psb_gefree_nest

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

    integer(psb_ipk_) :: i, j, linfo
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

    if (.not. allocated(a_nest%blk_present)) then
      allocate(a_nest%blk_present(a_nest%nrblocks, a_nest%ncblocks), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, name)
        return
      end if
      a_nest%blk_present = .false.
    end if

    do i = 1, a_nest%nrblocks
      do j = 1, a_nest%ncblocks
        linfo = psb_success_
        if (present(nnz)) then
          call psb_spall(a_nest%mats(i, j), desc_nest%descs(i, j), linfo, nnz=nnz)
        else
          call psb_spall(a_nest%mats(i, j), desc_nest%descs(i, j), linfo)
        end if
        if (linfo /= psb_success_) then
          info = linfo
          call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spall')
          return
        end if
        a_nest%blk_present(i, j) = .true.
      end do
    end do

  end subroutine psb_spall_nest


  ! Inserts nz entries into block (blk_i, blk_j) of the nested matrix.
  ! The block is lazy-allocated on first insertion if psb_spall_nest
  ! was not called first.
  
  subroutine psb_spins_nest(blk_i, blk_j, nz, ia, ja, val, a_nest, desc_nest, info)
    integer(psb_ipk_),           intent(in)    :: blk_i, blk_j, nz
    integer(psb_lpk_),           intent(in)    :: ia(:), ja(:)
    real(psb_dpk_),              intent(in)    :: val(:)
    type(psb_d_nest_sparse_mat), intent(inout) :: a_nest
    type(psb_desc_nest_type),    intent(inout) :: desc_nest
    integer(psb_ipk_),           intent(out)   :: info

    integer(psb_ipk_) :: nnz_est
    character(len=20) :: name

    info = psb_success_
    name = 'psb_spins_nest'

    if (nz == 0) return

    if (blk_i < 1 .or. blk_i > a_nest%nrblocks .or. &
        blk_j < 1 .or. blk_j > a_nest%ncblocks) then
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
      allocate(a_nest%blk_present(a_nest%nrblocks, a_nest%ncblocks), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, name)
        return
      end if
      a_nest%blk_present = .false.
    end if

    if (.not. a_nest%blk_present(blk_i, blk_j)) then
      ! Estimate nnz: use nz + 50% buffer for future insertions
      nnz_est = max(nz, 10) + nz / 2
      call psb_spall(a_nest%mats(blk_i, blk_j), &
                     desc_nest%descs(blk_i, blk_j), info, nnz=nnz_est)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spall')
        return
      end if
      a_nest%blk_present(blk_i, blk_j) = .true.
    end if

    call psb_spins(nz, ia, ja, val, a_nest%mats(blk_i, blk_j), &
                   desc_nest%descs(blk_i, blk_j), info)
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

    integer(psb_ipk_) :: i, j, dupl_, linfo
    character(len=20) :: name

    info  = psb_success_
    name  = 'psb_spasb_nest'
    dupl_ = psb_dupl_add_
    if (present(dupl)) dupl_ = dupl

    do i = 1, a_nest%nrblocks
      do j = 1, a_nest%ncblocks
        if (a_nest%blk_present(i, j)) then
          linfo = psb_success_
          if (dupl_ == psb_dupl_add_) then
            call psb_spasb(a_nest%mats(i, j), desc_nest%descs(i, j), linfo, &
                           dupl=psb_dupl_add_)
          else if (dupl_ == psb_dupl_ovwrt_) then
            call psb_spasb(a_nest%mats(i, j), desc_nest%descs(i, j), linfo, &
                           dupl=psb_dupl_ovwrt_)
          else if (dupl_ == psb_dupl_err_) then
            call psb_spasb(a_nest%mats(i, j), desc_nest%descs(i, j), linfo, &
                           dupl=psb_dupl_err_)
          else
            call psb_spasb(a_nest%mats(i, j), desc_nest%descs(i, j), linfo)
          end if
          if (linfo /= psb_success_) then
            info = linfo
            call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spasb')
            return
          end if
        end if
      end do
    end do

  end subroutine psb_spasb_nest

  ! Calls psb_spfree on every present block, then deallocates the
  ! mats and blk_present arrays and resets nrblocks/ncblocks to 0.
  
  subroutine psb_spfree_nest(a_nest, desc_nest, info)
    type(psb_d_nest_sparse_mat), intent(inout) :: a_nest
    type(psb_desc_nest_type),    intent(in)    :: desc_nest
    integer(psb_ipk_),           intent(out)   :: info

    integer(psb_ipk_) :: i, j, linfo
    character(len=20) :: name

    info = psb_success_
    name = 'psb_spfree_nest'

    if (allocated(a_nest%mats)) then
      do i = 1, a_nest%nrblocks
        do j = 1, a_nest%ncblocks
          if (allocated(a_nest%blk_present)) then
            if (a_nest%blk_present(i, j)) then
              linfo = psb_success_
              call psb_spfree(a_nest%mats(i, j), desc_nest%descs(i, j), linfo)
              if (linfo /= psb_success_ .and. info == psb_success_) then
                info = linfo
                call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_spfree')
              end if
            end if
          end if
        end do
      end do
      deallocate(a_nest%mats, stat=linfo)
      if (linfo /= 0 .and. info == psb_success_) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, name)
      end if
    end if

    if (allocated(a_nest%blk_present)) then
      deallocate(a_nest%blk_present, stat=linfo)
      if (linfo /= 0 .and. info == psb_success_) then
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

    integer(psb_ipk_) :: i, j, linfo
    character(len=20) :: name

    info = psb_success_
    name = 'psb_sprn_nest'

    if (.not. allocated(a_nest%mats) .or. .not. allocated(a_nest%blk_present)) return

    do i = 1, a_nest%nrblocks
      do j = 1, a_nest%ncblocks
        if (a_nest%blk_present(i, j)) then
          linfo = psb_success_
          if (present(clear)) then
            call psb_sprn(a_nest%mats(i, j), desc_nest%descs(i, j), linfo, clear=clear)
          else
            call psb_sprn(a_nest%mats(i, j), desc_nest%descs(i, j), linfo)
          end if
          if (linfo /= psb_success_ .and. info == psb_success_) then
            info = linfo
            call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_sprn')
          end if
        end if
      end do
    end do

  end subroutine psb_sprn_nest


  ! Allocates one sub-vector per block-row, using descs(i, 1) as
  ! the row descriptor for block i.  Must be called before psb_geins_nest.
  
  subroutine psb_geall_nest(x_nest, desc_nest, info)
    type(psb_d_nest_vect_type), intent(out) :: x_nest
    type(psb_desc_nest_type),   intent(in)  :: desc_nest
    integer(psb_ipk_),          intent(out) :: info

    integer(psb_ipk_) :: i, linfo
    character(len=20) :: name

    info = psb_success_
    name = 'psb_geall_nest'

    x_nest%nblocks = desc_nest%nrblocks
    allocate(x_nest%vects(x_nest%nblocks), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name)
      return
    end if

    do i = 1, x_nest%nblocks
      linfo = psb_success_
      call psb_geall(x_nest%vects(i), desc_nest%descs(i, 1), linfo)
      if (linfo /= psb_success_) then
        info = linfo
        call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_geall')
        return
      end if
    end do

  end subroutine psb_geall_nest

  ! Inserts m entries into block blk_i of the nested vector.

  subroutine psb_geins_nest(blk_i, m, irw, val, x_nest, desc_nest, info, local)
    integer(psb_ipk_),          intent(in)    :: blk_i, m
    integer(psb_lpk_),          intent(in)    :: irw(:)
    real(psb_dpk_),             intent(in)    :: val(:)
    type(psb_d_nest_vect_type), intent(inout) :: x_nest
    type(psb_desc_nest_type),   intent(in)    :: desc_nest
    integer(psb_ipk_),          intent(out)   :: info
    logical,                    intent(in), optional :: local

    character(len=20) :: name

    info = psb_success_
    name = 'psb_geins_nest'

    if (m == 0) return

    if (blk_i < 1 .or. blk_i > x_nest%nblocks) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='invalid block index')
      return
    end if

    if (present(local)) then
      call psb_geins(m, irw, val, x_nest%vects(blk_i), desc_nest%descs(blk_i, 1), info, &
                     local=local)
    else
      call psb_geins(m, irw, val, x_nest%vects(blk_i), desc_nest%descs(blk_i, 1), info)
    end if
    if (info /= psb_success_) &
      call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_geins')

  end subroutine psb_geins_nest

  ! Calls psb_geasb on every sub-vector.
  ! Must be called after psb_cdasb_nest and all psb_geins_nest calls.
  
  subroutine psb_geasb_nest(x_nest, desc_nest, info)
    type(psb_d_nest_vect_type), intent(inout) :: x_nest
    type(psb_desc_nest_type),   intent(in)    :: desc_nest
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i, linfo
    character(len=20) :: name

    info = psb_success_
    name = 'psb_geasb_nest'

    do i = 1, x_nest%nblocks
      linfo = psb_success_
      call psb_geasb(x_nest%vects(i), desc_nest%descs(i, 1), linfo)
      if (linfo /= psb_success_ .and. info == psb_success_) then
        info = linfo
        call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_geasb')
      end if
    end do

  end subroutine psb_geasb_nest

  ! Calls psb_gefree on every sub-vector, then deallocates the
  ! vects array and resets nblocks to 0.
  
  subroutine psb_gefree_nest(x_nest, desc_nest, info)
    type(psb_d_nest_vect_type), intent(inout) :: x_nest
    type(psb_desc_nest_type),   intent(in)    :: desc_nest
    integer(psb_ipk_),          intent(out)   :: info

    integer(psb_ipk_) :: i, linfo
    character(len=20) :: name

    info = psb_success_
    name = 'psb_gefree_nest'

    if (allocated(x_nest%vects)) then
      do i = 1, x_nest%nblocks
        linfo = psb_success_
        call psb_gefree(x_nest%vects(i), desc_nest%descs(i, 1), linfo)
        if (linfo /= psb_success_ .and. info == psb_success_) then
          info = linfo
          call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_gefree')
        end if
      end do
      deallocate(x_nest%vects, stat=linfo)
      if (linfo /= 0 .and. info == psb_success_) then
        info = psb_err_alloc_dealloc_
        call psb_errpush(info, name)
      end if
    end if

    x_nest%nblocks = 0

  end subroutine psb_gefree_nest

end module psb_d_nest_tools_mod
