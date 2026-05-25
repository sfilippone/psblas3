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
!  Nested-specific assembly wrappers for PSBLAS3 — descriptor routines
!

module psb_cd_nest_tools_mod
  use psb_const_mod, only : psb_ipk_, psb_lpk_, psb_success_, psb_err_alloc_dealloc_, &
                             psb_err_invalid_input_, psb_err_no_optional_arg_, psb_err_from_subroutine_, &
                             psb_ctxt_type
  use psb_error_mod, only : psb_errpush
  use psb_cd_tools_mod, only : psb_cdall, psb_cdasb, psb_cdins, psb_cdcpy, psb_cdprt
  use psb_desc_nest_mod, only : psb_desc_nest_type
  implicit none

  private

  public :: psb_cdall_nest, psb_cdins_nest, psb_cdins_nest_rc, &
            psb_cdasb_nest, psb_cdfree_nest, psb_cdcpy_nest, psb_cdprt_nest

  ! Column-only form: (blk_j, nz, ja, desc_nest, info [,mask, lidx])
  ! Row+column form: (blk_i, blk_j, nz, ia, ja, desc_nest, info)
  interface psb_cdins_nest
#if defined(PSB_IPK4) && defined(PSB_LPK8)
    module procedure psb_cdins_nest_c
    module procedure psb_cdins_nest_rc_sub
#endif
    module procedure psb_lcdins_nest_c
    module procedure psb_lcdins_nest_rc
  end interface

  ! Row+column form: (blk_i, blk_j, nz, ia, ja, desc_nest, info)
  interface psb_cdins_nest_rc
#if defined(PSB_IPK4) && defined(PSB_LPK8)
    module procedure psb_cdins_nest_rc_sub
#endif
    module procedure psb_lcdins_nest_rc
  end interface

contains

  ! Allocates the nested descriptor structure and creates block
  ! descriptors. The first block of each row uses psb_cdall with
  ! the given local row count; subsequent blocks in the same row
  ! are clones of the first block (same row distribution).
  !
  ! Arguments:
  !   ctxt        - PSBLAS context
  !   desc_nest   - nested descriptor (output)
  !   info        - error code (output)
  !   nrblocks    - number of block rows (optional, default 2)
  !   ncblocks    - number of block columns (optional, default 2)
  !   nl          - local row count per process (required for first blocks)

  subroutine psb_cdall_nest(ctxt, desc_nest, info, nrblocks, ncblocks, nl)
    type(psb_ctxt_type), intent(in) :: ctxt
    type(psb_desc_nest_type), intent(out) :: desc_nest
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in), optional :: nrblocks, ncblocks, nl

    integer(psb_ipk_) :: i, j, nr, nc, nl_
    character(len=20) :: name

    info = psb_success_
    name = 'psb_cdall_nest'

    ! Set default dimensions
    nr = 2
    nc = 2
    if (present(nrblocks)) nr = nrblocks
    if (present(ncblocks)) nc = ncblocks

    if (.not. present(nl)) then
      info = psb_err_no_optional_arg_
      call psb_errpush(info, name, a_err='nl (local row count)')
      return
    end if
    nl_ = nl

    ! Allocate nested descriptor structure
    desc_nest%nrblocks = nr
    desc_nest%ncblocks = nc
    allocate(desc_nest%descs(nr, nc), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name)
      return
    end if

    ! Build descriptors: each block gets its own independent psb_cdall.
    ! Cloning a build-state descriptor shares its base_desc pointer; when
    ! psb_cdasb_nest assembles both the original and the clone the shared
    ! base_desc is rebuilt twice, corrupting the global-to-local mapping of
    ! every block in that row.  Independent allocations avoid this entirely.
    do i = 1, nr
      do j = 1, nc
        call psb_cdall(ctxt, desc_nest%descs(i, j), info, nl=nl_)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_, name)
          return
        end if
      end do
    end do

  end subroutine psb_cdall_nest


#if defined(PSB_IPK4) && defined(PSB_LPK8)
  ! psb_cdins_nest_rc_sub: row+col form, ipk_ nz — only when ipk_ /= lpk_
  subroutine psb_cdins_nest_rc_sub(blk_i, blk_j, nz, ia, ja, desc_nest, info)
    integer(psb_ipk_), intent(in)    :: blk_i, blk_j, nz
    integer(psb_lpk_), intent(in)    :: ia(:), ja(:)
    type(psb_desc_nest_type), intent(inout) :: desc_nest
    integer(psb_ipk_), intent(out)   :: info

    character(len=20) :: name

    info = psb_success_
    name = 'psb_cdins_nest'

    if (nz == 0) return

    if (blk_i < 1 .or. blk_i > desc_nest%nrblocks .or. &
        blk_j < 1 .or. blk_j > desc_nest%ncblocks) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='invalid block indices')
      return
    end if

    call psb_cdins(nz, ia, ja, desc_nest%descs(blk_i, blk_j), info)
    if (info /= psb_success_) &
      call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdins')

  end subroutine psb_cdins_nest_rc_sub


  ! psb_cdins_nest_c: col-only form, ipk_ nz — only when ipk_ /= lpk_
  subroutine psb_cdins_nest_c(blk_j, nz, ja, desc_nest, info, mask, lidx)
    integer(psb_ipk_), intent(in)    :: blk_j, nz
    integer(psb_lpk_), intent(in)    :: ja(:)
    type(psb_desc_nest_type), intent(inout) :: desc_nest
    integer(psb_ipk_), intent(out)   :: info
    logical,           intent(in), optional, target :: mask(:)
    integer(psb_ipk_), intent(in), optional         :: lidx(:)

    integer(psb_ipk_) :: i, linfo
    character(len=20) :: name

    info = psb_success_
    name = 'psb_cdins_nest'

    if (nz == 0) return

    if (blk_j < 1 .or. blk_j > desc_nest%ncblocks) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='invalid block column index')
      return
    end if

    do i = 1, desc_nest%nrblocks
      linfo = psb_success_
      if (present(mask)) then
        if (present(lidx)) then
          call psb_cdins(nz, ja, desc_nest%descs(i, blk_j), linfo, mask=mask, lidx=lidx)
        else
          call psb_cdins(nz, ja, desc_nest%descs(i, blk_j), linfo, mask=mask)
        end if
      else
        if (present(lidx)) then
          call psb_cdins(nz, ja, desc_nest%descs(i, blk_j), linfo, lidx=lidx)
        else
          call psb_cdins(nz, ja, desc_nest%descs(i, blk_j), linfo)
        end if
      end if
      if (linfo /= psb_success_ .and. info == psb_success_) then
        info = linfo
        call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdins')
      end if
    end do

  end subroutine psb_cdins_nest_c
#endif


  ! psb_lcdins_nest_rc: row+col form, lpk_ nz
  !
  ! When entries in block (blk_i, blk_j) reference columns owned by other
  ! processes, use the col-only form afterwards to broadcast those column
  ! indices across all row-blocks in block-col blk_j.
  subroutine psb_lcdins_nest_rc(blk_i, blk_j, nz, ia, ja, desc_nest, info)
    integer(psb_ipk_), intent(in)    :: blk_i, blk_j
    integer(psb_lpk_), intent(in)    :: nz, ia(:), ja(:)
    type(psb_desc_nest_type), intent(inout) :: desc_nest
    integer(psb_ipk_), intent(out)   :: info

    character(len=20) :: name

    info = psb_success_
    name = 'psb_cdins_nest'

    if (nz == 0) return

    if (blk_i < 1 .or. blk_i > desc_nest%nrblocks .or. &
        blk_j < 1 .or. blk_j > desc_nest%ncblocks) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='invalid block indices')
      return
    end if

    call psb_cdins(nz, ia, ja, desc_nest%descs(blk_i, blk_j), info)
    if (info /= psb_success_) &
      call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdins')

  end subroutine psb_lcdins_nest_rc


  ! psb_lcdins_nest_c: col-only form, lpk_ nz
  !
  ! Registers nz global column indices ja into the descriptor for
  ! block column blk_j across all row-blocks (descs(i, blk_j) for
  ! i = 1..nrblocks).  mask and lidx are forwarded to psb_cdins.
  subroutine psb_lcdins_nest_c(blk_j, nz, ja, desc_nest, info, mask, lidx)
    integer(psb_ipk_), intent(in)    :: blk_j
    integer(psb_lpk_), intent(in)    :: nz, ja(:)
    type(psb_desc_nest_type), intent(inout) :: desc_nest
    integer(psb_ipk_), intent(out)   :: info
    logical,           intent(in), optional, target :: mask(:)
    integer(psb_ipk_), intent(in), optional         :: lidx(:)

    integer(psb_ipk_) :: i, linfo
    character(len=20) :: name

    info = psb_success_
    name = 'psb_cdins_nest'

    if (nz == 0) return

    if (blk_j < 1 .or. blk_j > desc_nest%ncblocks) then
      info = psb_err_invalid_input_
      call psb_errpush(info, name, a_err='invalid block column index')
      return
    end if

    do i = 1, desc_nest%nrblocks
      linfo = psb_success_
      if (present(mask)) then
        if (present(lidx)) then
          call psb_cdins(nz, ja, desc_nest%descs(i, blk_j), linfo, mask=mask, lidx=lidx)
        else
          call psb_cdins(nz, ja, desc_nest%descs(i, blk_j), linfo, mask=mask)
        end if
      else
        if (present(lidx)) then
          call psb_cdins(nz, ja, desc_nest%descs(i, blk_j), linfo, lidx=lidx)
        else
          call psb_cdins(nz, ja, desc_nest%descs(i, blk_j), linfo)
        end if
      end if
      if (linfo /= psb_success_ .and. info == psb_success_) then
        info = linfo
        call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdins')
      end if
    end do

  end subroutine psb_lcdins_nest_c


  ! psb_cdasb_nest: Finalize all nested descriptors
  ! 
  ! Calls psb_cdasb on all block descriptors in the nested structure.
  ! This must be called after all psb_cdins_nest calls and
  ! before psb_spasb_nest.
  !
  ! Arguments:
  !   desc_nest   - nested descriptor (input/output)
  !   info        - error code (output)

  subroutine psb_cdasb_nest(desc_nest, info)
    type(psb_desc_nest_type), intent(inout) :: desc_nest
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, j
    character(len=20) :: name

    info = psb_success_
    name = 'psb_cdasb_nest'

    do i = 1, desc_nest%nrblocks
      do j = 1, desc_nest%ncblocks
        call psb_cdasb(desc_nest%descs(i, j), info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdasb')
          return
        end if
      end do
    end do

  end subroutine psb_cdasb_nest


  ! psb_cdfree_nest: Free all nested descriptors
  ! 
  ! Calls psb_cdfree on every block descriptor in the nested
  ! structure, then deallocates the descriptor array and resets
  ! nrblocks/ncblocks to 0.  Mirrors what psb_cdfree does for a
  ! single psb_desc_type.
  !
  ! Arguments:
  !   desc_nest   - nested descriptor (input/output)
  !   info        - error code (output)
  ! 
  subroutine psb_cdfree_nest(desc_nest, info)
    type(psb_desc_nest_type), intent(inout) :: desc_nest
    integer(psb_ipk_), intent(out) :: info

    character(len=20) :: name

    info = psb_success_
    name = 'psb_cdfree_nest'

    call desc_nest%free(info)
    if (info /= psb_success_) &
      call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_desc_nest_free')

  end subroutine psb_cdfree_nest


  ! psb_cdcpy_nest: Deep copy (clone) a nested descriptor
  !
  ! Allocates desc_out and clones each block descriptor from desc_in
  ! using psb_cdcpy, preserving the full row/column block structure.
  !
  ! Arguments:
  !   desc_in   - source nested descriptor (inout — clone may need to read internal state)
  !   desc_out  - destination nested descriptor (output)
  !   info      - error code (output)

  subroutine psb_cdcpy_nest(desc_in, desc_out, info)
    type(psb_desc_nest_type), intent(inout) :: desc_in
    type(psb_desc_nest_type), intent(out)   :: desc_out
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, j
    character(len=20) :: name

    info = psb_success_
    name = 'psb_cdcpy_nest'

    desc_out%nrblocks = desc_in%nrblocks
    desc_out%ncblocks = desc_in%ncblocks
    allocate(desc_out%descs(desc_in%nrblocks, desc_in%ncblocks), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name)
      return
    end if

    do i = 1, desc_in%nrblocks
      do j = 1, desc_in%ncblocks
        call psb_cdcpy(desc_in%descs(i, j), desc_out%descs(i, j), info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_, name, a_err='psb_cdcpy')
          return
        end if
      end do
    end do

  end subroutine psb_cdcpy_nest


  ! psb_cdprt_nest: Print all block descriptors (debugging)
  !
  ! Loops over all (i,j) block descriptors in the nested structure
  ! and calls psb_cdprt on each, prefixing the output with the block
  ! index.  All optional arguments are forwarded unchanged.
  !
  ! Arguments:
  !   iout      - output unit
  !   desc_nest - nested descriptor (input)
  !   glob      - passed to psb_cdprt (optional)
  !   short     - passed to psb_cdprt (optional)
  !   verbosity - passed to psb_cdprt (optional)

  subroutine psb_cdprt_nest(iout, desc_nest, glob, short, verbosity)
    integer(psb_ipk_), intent(in)            :: iout
    type(psb_desc_nest_type), intent(in)     :: desc_nest
    logical,           intent(in), optional  :: glob, short
    integer(psb_ipk_), intent(in), optional  :: verbosity

    integer(psb_ipk_) :: i, j

    do i = 1, desc_nest%nrblocks
      do j = 1, desc_nest%ncblocks
        write(iout, '(a,i0,a,i0,a)') 'Block (', i, ',', j, '):'
        if (present(glob)) then
          if (present(short)) then
            if (present(verbosity)) then
              call psb_cdprt(iout, desc_nest%descs(i,j), glob=glob, short=short, verbosity=verbosity)
            else
              call psb_cdprt(iout, desc_nest%descs(i,j), glob=glob, short=short)
            end if
          else
            if (present(verbosity)) then
              call psb_cdprt(iout, desc_nest%descs(i,j), glob=glob, verbosity=verbosity)
            else
              call psb_cdprt(iout, desc_nest%descs(i,j), glob=glob)
            end if
          end if
        else
          if (present(short)) then
            if (present(verbosity)) then
              call psb_cdprt(iout, desc_nest%descs(i,j), short=short, verbosity=verbosity)
            else
              call psb_cdprt(iout, desc_nest%descs(i,j), short=short)
            end if
          else
            if (present(verbosity)) then
              call psb_cdprt(iout, desc_nest%descs(i,j), verbosity=verbosity)
            else
              call psb_cdprt(iout, desc_nest%descs(i,j))
            end if
          end if
        end if
      end do
    end do

  end subroutine psb_cdprt_nest

end module psb_cd_nest_tools_mod
