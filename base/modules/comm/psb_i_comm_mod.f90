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
module psb_i_comm_mod
  use psb_desc_mod, only : psb_desc_type, psb_ipk_, psb_lpk_, psb_epk_, psb_mpk_
  
  use psb_i_vect_mod, only : psb_i_vect_type, psb_i_base_vect_type
  use psb_i_multivect_mod, only : psb_i_multivect_type, psb_i_base_multivect_type

  interface psb_ovrl
    subroutine psb_iovrl_vect(x,desc_a,info,work,update,mode)
      import
      implicit none
      type(psb_i_vect_type), intent(inout)    :: x
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), intent(inout), optional, target :: work(:)
      integer(psb_ipk_), intent(in), optional           :: update,mode
    end subroutine psb_iovrl_vect
    subroutine psb_iovrl_multivect(x,desc_a,info,work,update,mode)
      import
      implicit none
      type(psb_i_multivect_type), intent(inout)    :: x
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), intent(inout), optional, target :: work(:)
      integer(psb_ipk_), intent(in), optional           :: update,mode
    end subroutine psb_iovrl_multivect
  end interface psb_ovrl

  interface psb_halo
    subroutine psb_ihalo_vect(x,desc_a,info,work,tran,mode,data)
      import
      implicit none
      type(psb_i_vect_type), intent(inout)   :: x
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), target, optional, intent(inout) :: work(:)
      integer(psb_ipk_), intent(in), optional           :: mode,data
      character, intent(in), optional         :: tran
    end subroutine psb_ihalo_vect
    subroutine psb_ihalo_multivect(x,desc_a,info,work,tran,mode,data)
      import
      implicit none
      type(psb_i_multivect_type), intent(inout)   :: x
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)                    :: info
      integer(psb_ipk_), target, optional, intent(inout) :: work(:)
      integer(psb_ipk_), intent(in), optional           :: mode,data
      character, intent(in), optional         :: tran
    end subroutine psb_ihalo_multivect
  end interface psb_halo


  interface psb_scatter
    subroutine  psb_iscatter_vect(globx, locx, desc_a, info, root, mold)
      import
      implicit none
      type(psb_i_vect_type), intent(inout) :: locx
      integer(psb_ipk_), intent(in)  :: globx(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), intent(in), optional    :: root
      class(psb_i_base_vect_type), intent(in), optional :: mold  
    end subroutine psb_iscatter_vect
  end interface psb_scatter

  interface psb_gather
    subroutine psb_igather_vect(globx, locx, desc_a, info, root)
      import
      implicit none
      type(psb_i_vect_type), intent(inout) :: locx
      integer(psb_ipk_), intent(out), allocatable :: globx(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: root
    end subroutine psb_igather_vect
    subroutine psb_igather_multivect(globx, locx, desc_a, info, root)
      import
      implicit none
      type(psb_i_multivect_type), intent(inout) :: locx
      integer(psb_ipk_), intent(out), allocatable :: globx(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: root
    end subroutine psb_igather_multivect
  end interface psb_gather

end module psb_i_comm_mod
