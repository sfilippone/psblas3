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
Module psb_i_tools_mod
  use psb_desc_mod, only : psb_desc_type, psb_ipk_, psb_lpk_, psb_success_
  use psb_i_vect_mod, only : psb_i_base_vect_type, psb_i_vect_type
  use psb_m_tools_a_mod
  use psb_e_tools_a_mod
  use psb_l_vect_mod, only : psb_l_vect_type
  use psb_i_multivect_mod, only : psb_i_base_multivect_type, psb_i_multivect_type
  use psi_mod, only : psb_snd, psb_rcv ! Needed only for psb_getelem

  interface  psb_geall
    subroutine psb_ialloc_vect(x, desc_a,info, dupl, bldmode)
      import
      implicit none
      type(psb_i_vect_type), intent(out)  :: x
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_),intent(out)             :: info
      integer(psb_ipk_), optional, intent(in) :: dupl, bldmode
    end subroutine psb_ialloc_vect
    subroutine psb_ialloc_vect_r2(x, desc_a,info,n,lb, dupl, bldmode)
      import
      implicit none
      type(psb_i_vect_type), allocatable, intent(out)  :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_),intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)   :: n, lb
      integer(psb_ipk_), optional, intent(in) :: dupl, bldmode
    end subroutine psb_ialloc_vect_r2
    subroutine psb_ialloc_multivect(x, desc_a,info,n, dupl, bldmode)
      import
      implicit none
      type(psb_i_multivect_type), intent(out)  :: x
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_),intent(out)             :: info
      integer(psb_ipk_), optional, intent(in)   :: n
      integer(psb_ipk_), optional, intent(in) :: dupl, bldmode
    end subroutine psb_ialloc_multivect
  end interface


  interface psb_geasb
    subroutine psb_iasb_vect(x, desc_a, info,mold, scratch)
      import
      implicit none
      type(psb_desc_type), intent(in)      ::  desc_a
      type(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)                 ::  info
      class(psb_i_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine psb_iasb_vect
    subroutine psb_iasb_vect_r2(x, desc_a, info,mold, scratch)
      import
      implicit none
      type(psb_desc_type), intent(in)      ::  desc_a
      type(psb_i_vect_type), intent(inout) :: x(:)
      integer(psb_ipk_), intent(out)                 ::  info
      class(psb_i_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine psb_iasb_vect_r2
    subroutine psb_iasb_multivect(x, desc_a, info,mold, scratch, n)
      import
      implicit none
      type(psb_desc_type), intent(in)      ::  desc_a
      type(psb_i_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)                 ::  info
      class(psb_i_base_multivect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
      integer(psb_ipk_), optional, intent(in)   :: n
    end subroutine psb_iasb_multivect
  end interface

  interface psb_gefree
    subroutine psb_ifree_vect(x, desc_a, info)
      import
      implicit none
      type(psb_desc_type), intent(in)  ::  desc_a
      type(psb_i_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             ::  info
    end subroutine psb_ifree_vect
    subroutine psb_ifree_vect_r2(x, desc_a, info)
      import
      implicit none
      type(psb_desc_type), intent(in)  ::  desc_a
      type(psb_i_vect_type), allocatable, intent(inout) :: x(:)
      integer(psb_ipk_), intent(out)             ::  info
    end subroutine psb_ifree_vect_r2
    subroutine psb_ifree_multivect(x, desc_a, info)
      import
      implicit none
      type(psb_desc_type), intent(in)  ::  desc_a
      type(psb_i_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             ::  info
    end subroutine psb_ifree_multivect
  end interface


  interface psb_geins
    subroutine psb_iins_vect(m,irw,val,x,desc_a,info,local)
      import
      implicit none
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_i_vect_type), intent(inout) :: x
      integer(psb_lpk_), intent(in)              :: irw(:)
      integer(psb_ipk_), intent(in)    :: val(:)
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: local
    end subroutine psb_iins_vect
    subroutine psb_iins_vect_v(m,irw,val,x,desc_a,info,local)
      import
      implicit none
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_i_vect_type), intent(inout) :: x
      type(psb_l_vect_type), intent(inout)       :: irw
      type(psb_i_vect_type), intent(inout)    :: val
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: local
    end subroutine psb_iins_vect_v
    subroutine psb_iins_vect_r2(m,irw,val,x,desc_a,info,local)
      import
      implicit none
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_i_vect_type), intent(inout) :: x(:)
      integer(psb_lpk_), intent(in)              :: irw(:)
      integer(psb_ipk_), intent(in)    :: val(:,:)
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: local
    end subroutine psb_iins_vect_r2
    subroutine psb_iins_multivect(m,irw,val,x,desc_a,info,local)
      import
      implicit none
      integer(psb_ipk_), intent(in)              :: m
      type(psb_desc_type), intent(in)  :: desc_a
      type(psb_i_multivect_type), intent(inout) :: x
      integer(psb_lpk_), intent(in)              :: irw(:)
      integer(psb_ipk_), intent(in)    :: val(:,:)
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: local
    end subroutine psb_iins_multivect
  end interface

end module psb_i_tools_mod
