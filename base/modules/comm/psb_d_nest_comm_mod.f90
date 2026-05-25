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
! module: psb_d_nest_comm_mod
!
! Communication operations for nested (block-structured) double precision
! real vectors.
!
! psb_d_nest_halo
!   Halo exchange for all column blocks of a nested vector.
!   Calls psb_halo(x(j), descs(1,j)) for each column block j.
!   All descriptors descs(i,j) for fixed j are equivalent;
!   Called once before block SpMM to populate ghost entries of x.
!
! psb_d_nest_ovrl
!   Overlap update for all row blocks of a nested vector.
!   Calls psb_ovrl(x(i), descs(i,i)) for each row block i using the
!   diagonal descriptor.
!   Called after operations that contribute to overlapping rows
!   (e.g. FEM assembly).
!
module psb_d_nest_comm_mod
  use psb_desc_nest_mod
  use psb_d_nest_vect_mod
  use psb_d_comm_mod,  only : psb_halo, psb_ovrl
  use psb_const_mod,   only : psb_ipk_
  implicit none

  private
  public :: psb_d_nest_halo, psb_d_nest_ovrl

contains

  subroutine psb_d_nest_halo(xnest, descs, info, tran, mode, data)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    character, optional,        intent(in)    :: tran
    integer(psb_ipk_), optional, intent(in)   :: mode, data

    integer(psb_ipk_) :: j

    info = 0
    do j = 1, xnest%nblocks
      call psb_halo(xnest%vects(j), descs%descs(1,j), info, tran=tran, mode=mode, data=data)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_halo

  subroutine psb_d_nest_ovrl(xnest, descs, info, update, mode)
    type(psb_d_nest_vect_type), intent(inout) :: xnest
    type(psb_desc_nest_type),   intent(in)    :: descs
    integer(psb_ipk_),          intent(out)   :: info
    integer(psb_ipk_), optional, intent(in)   :: update, mode

    integer(psb_ipk_) :: i

    info = 0
    do i = 1, xnest%nblocks
      call psb_ovrl(xnest%vects(i), descs%descs(i,i), info, update=update, mode=mode)
      if (info /= 0) return
    end do
  end subroutine psb_d_nest_ovrl

end module psb_d_nest_comm_mod
