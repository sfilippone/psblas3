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
module psi_i_comm_v_mod
  use psi_penv_mod, only : psb_ctxt_type
  use psb_desc_mod, only : psb_desc_type, psb_ipk_, psb_mpk_, &
       & psb_lpk_, psb_epk_, psb_i2pk_
  use psb_i_base_vect_mod, only : psb_i_base_vect_type 
  use psb_i_base_multivect_mod, only : psb_i_base_multivect_type 

  interface psi_swapdata
    subroutine psi_iswapdata_vect(flag,beta,y,desc_a,work,info,data)
      import 
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_vect_type) :: y
      integer(psb_ipk_)           :: beta 
      integer(psb_ipk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional :: data
    end subroutine psi_iswapdata_vect
    subroutine psi_iswapdata_multivect(flag,beta,y,desc_a,work,info,data)
      import 
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_multivect_type)    :: y
      integer(psb_ipk_)           :: beta 
      integer(psb_ipk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional :: data
    end subroutine psi_iswapdata_multivect
    subroutine psi_iswap_vidx_vect(ictxt,iicomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import 
      type(psb_ctxt_type), intent(in)         :: ictxt
      integer(psb_mpk_), intent(in)           :: iicomm
      integer(psb_ipk_), intent(in)           :: flag
      integer(psb_ipk_), intent(out)          :: info
      class(psb_i_base_vect_type)             :: y
      integer(psb_ipk_)                       :: beta
      integer(psb_ipk_), target               :: work(:)
      class(psb_i_base_vect_type), intent(inout) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_iswap_vidx_vect
    subroutine psi_iswap_vidx_multivect(ictxt,iicomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import 
      type(psb_ctxt_type), intent(in)       :: ictxt
      integer(psb_mpk_), intent(in)         :: iicomm
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_multivect_type)    :: y
      integer(psb_ipk_)                       :: beta
      integer(psb_ipk_), target               :: work(:)
      class(psb_i_base_vect_type), intent(inout) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_iswap_vidx_multivect
  end interface psi_swapdata


  interface psi_swaptran
    subroutine psi_iswaptran_vect(flag,beta,y,desc_a,work,info,data)
      import 
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_vect_type) :: y
      integer(psb_ipk_)           :: beta
      integer(psb_ipk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_iswaptran_vect
    subroutine psi_iswaptran_multivect(flag,beta,y,desc_a,work,info,data)
      import 
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_multivect_type) :: y
      integer(psb_ipk_)           :: beta
      integer(psb_ipk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_iswaptran_multivect
    subroutine psi_itran_vidx_vect(ictxt,iicomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import 
      type(psb_ctxt_type), intent(in)       :: ictxt
      integer(psb_mpk_), intent(in)         :: iicomm
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_vect_type)          :: y
      integer(psb_ipk_)                       :: beta
      integer(psb_ipk_), target               :: work(:)
      class(psb_i_base_vect_type), intent(inout) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_itran_vidx_vect
    subroutine psi_itran_vidx_multivect(ictxt,iicomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import 
      type(psb_ctxt_type), intent(in)       :: ictxt
      integer(psb_mpk_), intent(in)         :: iicomm
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_multivect_type)      :: y
      integer(psb_ipk_)                       :: beta
      integer(psb_ipk_), target               :: work(:)
      class(psb_i_base_vect_type), intent(inout) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_itran_vidx_multivect
  end interface psi_swaptran

  interface psi_ovrl_upd
    subroutine  psi_iovrl_upd_vect(x,desc_a,update,info)
      import 
      class(psb_i_base_vect_type)       :: x
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(in)               :: update
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psi_iovrl_upd_vect
    subroutine  psi_iovrl_upd_multivect(x,desc_a,update,info)
      import 
      class(psb_i_base_multivect_type)   :: x
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(in)        :: update
      integer(psb_ipk_), intent(out)       :: info
    end subroutine psi_iovrl_upd_multivect
  end interface psi_ovrl_upd

  interface psi_ovrl_save
    subroutine  psi_iovrl_save_vect(x,xs,desc_a,info)
      import 
      class(psb_i_base_vect_type)     :: x
      integer(psb_ipk_), allocatable  :: xs(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_iovrl_save_vect
    subroutine  psi_iovrl_save_multivect(x,xs,desc_a,info)
      import 
      class(psb_i_base_multivect_type)     :: x
      integer(psb_ipk_), allocatable  :: xs(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_iovrl_save_multivect
  end interface psi_ovrl_save

  interface psi_ovrl_restore
    subroutine  psi_iovrl_restr_vect(x,xs,desc_a,info)
      import 
      class(psb_i_base_vect_type)     :: x
      integer(psb_ipk_)               :: xs(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_iovrl_restr_vect
    subroutine  psi_iovrl_restr_multivect(x,xs,desc_a,info)
      import 
      class(psb_i_base_multivect_type)     :: x
      integer(psb_ipk_)               :: xs(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_iovrl_restr_multivect
  end interface psi_ovrl_restore

end module psi_i_comm_v_mod

