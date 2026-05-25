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
! module: psb_desc_nest_mod
!
! Defines psb_desc_nest_type: a 2-D array of psb_desc_type objects,
! one per block matrix entry in an nrblocks x ncblocks block system.
!
!
module psb_desc_nest_mod
  use psb_desc_mod
  implicit none

  type :: psb_desc_nest_type
    integer(psb_ipk_)                       :: nrblocks = 0
    integer(psb_ipk_)                       :: ncblocks = 0
    type(psb_desc_type), allocatable        :: descs(:,:)
  contains
    procedure :: get_nrblocks => psb_desc_nest_get_nrblocks
    procedure :: get_ncblocks => psb_desc_nest_get_ncblocks
    procedure :: get_desc     => psb_desc_nest_get_desc
    procedure :: is_valid     => psb_desc_nest_is_valid
    procedure :: sizeof       => psb_desc_nest_sizeof
    procedure :: free         => psb_desc_nest_free
  end type psb_desc_nest_type

contains


  !  get_nrblocks / get_ncblocks
  function psb_desc_nest_get_nrblocks(d) result(nb)
    class(psb_desc_nest_type), intent(in) :: d
    integer(psb_ipk_) :: nb
    nb = d%nrblocks
  end function psb_desc_nest_get_nrblocks

  function psb_desc_nest_get_ncblocks(d) result(nb)
    class(psb_desc_nest_type), intent(in) :: d
    integer(psb_ipk_) :: nb
    nb = d%ncblocks
  end function psb_desc_nest_get_ncblocks

  !  get_desc: copy descriptor (i,j) into the output argument
  subroutine psb_desc_nest_get_desc(d, i, j, desc, info)
    class(psb_desc_nest_type), intent(in) :: d
    integer(psb_ipk_),         intent(in) :: i, j
    type(psb_desc_type),       intent(out):: desc
    integer(psb_ipk_),         intent(out):: info

    info = 0
    if (i < 1 .or. i > d%nrblocks .or. j < 1 .or. j > d%ncblocks) then
      info = -1
      return
    end if
    desc = d%descs(i,j)
  end subroutine psb_desc_nest_get_desc

  !  is_valid: true if all diagonal sub-descriptors are valid
  function psb_desc_nest_is_valid(d) result(valid)
    class(psb_desc_nest_type), intent(in) :: d
    logical :: valid
    integer(psb_ipk_) :: i

    valid = (d%nrblocks >= 1) .and. (d%ncblocks >= 1) .and. allocated(d%descs)
    if (valid) then
      do i = 1, min(d%nrblocks, d%ncblocks)
        if (.not. d%descs(i,i)%is_valid()) then
          valid = .false.
          return
        end if
      end do
    end if
  end function psb_desc_nest_is_valid

  !  sizeof: total memory (bytes) of all sub-descriptors
  function psb_desc_nest_sizeof(d) result(s)
    class(psb_desc_nest_type), intent(in) :: d
    integer(psb_epk_) :: s
    integer(psb_ipk_) :: i, j

    s = 0_psb_epk_
    if (allocated(d%descs)) then
      do j = 1, d%ncblocks
        do i = 1, d%nrblocks
          s = s + d%descs(i,j)%sizeof()
        end do
      end do
    end if
  end function psb_desc_nest_sizeof

  !  free: release all sub-descriptors and reset
  subroutine psb_desc_nest_free(d, info)
    class(psb_desc_nest_type), intent(inout) :: d
    integer(psb_ipk_),         intent(out)   :: info

    integer(psb_ipk_) :: i, j, linfo

    info = 0
    if (allocated(d%descs)) then
      do j = 1, d%ncblocks
        do i = 1, d%nrblocks
          call d%descs(i,j)%free(linfo)
          if (linfo /= 0 .and. info == 0) info = linfo
        end do
      end do
      deallocate(d%descs, stat=linfo)
      if (linfo /= 0 .and. info == 0) info = linfo
    end if
    d%nrblocks = 0
    d%ncblocks = 0
  end subroutine psb_desc_nest_free

end module psb_desc_nest_mod
