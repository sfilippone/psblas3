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
! module: psb_d_nest_vect_mod
!
! Defines psb_d_nest_vect_type: a block-structured distributed dense
! vector for double precision real arithmetic.  Each sub-vector is a
! standard psb_d_vect_type assembled under its own descriptor.
!
! Parallel BLAS operations (nrm2, dot, axpby) are exposed as module
! subroutines/functions in psb_d_nest_psblas_mod so that they can
! exploit a single global reduction per call.
!
module psb_d_nest_vect_mod
  use psb_d_vect_mod
  use psb_desc_mod
  implicit none

  type :: psb_d_nest_vect_type
    integer(psb_ipk_)                       :: nblocks = 0
    type(psb_d_vect_type), allocatable      :: vects(:)
  contains
    procedure :: get_nblocks => psb_d_nest_vect_get_nblocks
    procedure :: zero        => psb_d_nest_vect_zero
    procedure :: sizeof      => psb_d_nest_vect_sizeof
    procedure :: free        => psb_d_nest_vect_free
  end type psb_d_nest_vect_type

contains

  !  get_nblocks
  function psb_d_nest_vect_get_nblocks(x) result(nb)
    class(psb_d_nest_vect_type), intent(in) :: x
    integer(psb_ipk_) :: nb
    nb = x%nblocks
  end function psb_d_nest_vect_get_nblocks

  !  zero: set all sub-vectors to zero (local, no halo zeroing needed)
  subroutine psb_d_nest_vect_zero(x)
    class(psb_d_nest_vect_type), intent(inout) :: x
    integer(psb_ipk_) :: i
    if (allocated(x%vects)) then
      do i = 1, x%nblocks
        call x%vects(i)%zero()
      end do
    end if
  end subroutine psb_d_nest_vect_zero

  !  sizeof: total bytes across all sub-vectors
  function psb_d_nest_vect_sizeof(x) result(s)
    class(psb_d_nest_vect_type), intent(in) :: x
    integer(psb_epk_) :: s
    integer(psb_ipk_) :: i
    s = 0_psb_epk_
    if (allocated(x%vects)) then
      do i = 1, x%nblocks
        s = s + x%vects(i)%sizeof()
      end do
    end if
  end function psb_d_nest_vect_sizeof

  !  free: release all sub-vectors
  subroutine psb_d_nest_vect_free(x, info)
    class(psb_d_nest_vect_type), intent(inout) :: x
    integer(psb_ipk_),           intent(out)   :: info

    integer(psb_ipk_) :: i, linfo

    info = 0
    if (allocated(x%vects)) then
      do i = 1, x%nblocks
        call x%vects(i)%free(linfo)
        if (linfo /= 0 .and. info == 0) info = linfo
      end do
      deallocate(x%vects, stat=linfo)
      if (linfo /= 0 .and. info == 0) info = linfo
    end if
    x%nblocks = 0
  end subroutine psb_d_nest_vect_free

end module psb_d_nest_vect_mod
