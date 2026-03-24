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
module psi_i2_comm_v_mod
  use psi_penv_mod, only : psb_ctxt_type
  use psb_desc_mod, only : psb_desc_type, psb_ipk_, psb_mpk_, &
       & psb_lpk_, psb_epk_, psb_i2pk_
  use psb_desc_mod, only : psb_desc_type, psb_ipk_, psb_mpk_, psb_lpk_, psb_epk_, psb_i_base_vect_type
  use psb_i2_base_vect_mod, only : psb_i2_base_vect_type 
  use psb_i2_base_multivect_mod, only : psb_i2_base_multivect_type 

  interface psi_swapdata
    module subroutine psi_i2swapdata_vect(flag,beta,y,desc_a,work,info,data)
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i2_base_vect_type) :: y
      integer(psb_i2pk_)           :: beta 
      integer(psb_i2pk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional :: data
    end subroutine psi_i2swapdata_vect
    module subroutine psi_i2swapdata_multivect(flag,beta,y,desc_a,work,info,data)
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i2_base_multivect_type)    :: y
      integer(psb_i2pk_)           :: beta 
      integer(psb_i2pk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional :: data
    end subroutine psi_i2swapdata_multivect
    module subroutine psi_i2swap_vidx_vect(ctxt,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      type(psb_ctxt_type), intent(in)         :: ctxt
      integer(psb_ipk_), intent(in)           :: flag
      integer(psb_ipk_), intent(out)          :: info
      class(psb_i2_base_vect_type)             :: y
      integer(psb_i2pk_)                       :: beta
      integer(psb_i2pk_), target               :: work(:)
      class(psb_i_base_vect_type), intent(inout) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_i2swap_vidx_vect
    module subroutine psi_i2swap_vidx_multivect(ctxt,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      type(psb_ctxt_type), intent(in)       :: ctxt
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i2_base_multivect_type)    :: y
      integer(psb_i2pk_)                       :: beta
      integer(psb_i2pk_), target               :: work(:)
      class(psb_i_base_vect_type), intent(inout) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_i2swap_vidx_multivect
  end interface psi_swapdata


  interface psi_swaptran
    module subroutine psi_i2swaptran_vect(flag,beta,y,desc_a,work,info,data)
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i2_base_vect_type) :: y
      integer(psb_i2pk_)           :: beta
      integer(psb_i2pk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_i2swaptran_vect
    module subroutine psi_i2swaptran_multivect(flag,beta,y,desc_a,work,info,data)
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i2_base_multivect_type) :: y
      integer(psb_i2pk_)           :: beta
      integer(psb_i2pk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_i2swaptran_multivect
    module subroutine psi_i2tran_vidx_vect(ctxt,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      type(psb_ctxt_type), intent(in)       :: ctxt
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i2_base_vect_type)          :: y
      integer(psb_i2pk_)                       :: beta
      integer(psb_i2pk_), target               :: work(:)
      class(psb_i_base_vect_type), intent(inout) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_i2tran_vidx_vect
    module subroutine psi_i2tran_vidx_multivect(ctxt,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      type(psb_ctxt_type), intent(in)       :: ctxt
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i2_base_multivect_type)      :: y
      integer(psb_i2pk_)                       :: beta
      integer(psb_i2pk_), target               :: work(:)
      class(psb_i_base_vect_type), intent(inout) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_i2tran_vidx_multivect
  end interface psi_swaptran

  interface psi_ovrl_upd
    module subroutine  psi_i2ovrl_upd_vect(x,desc_a,update,info)
      class(psb_i2_base_vect_type)       :: x
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(in)               :: update
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psi_i2ovrl_upd_vect
    module subroutine  psi_i2ovrl_upd_multivect(x,desc_a,update,info)
      class(psb_i2_base_multivect_type)   :: x
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(in)        :: update
      integer(psb_ipk_), intent(out)       :: info
    end subroutine psi_i2ovrl_upd_multivect
  end interface psi_ovrl_upd

  interface psi_ovrl_save
    module subroutine  psi_i2ovrl_save_vect(x,xs,desc_a,info)
      class(psb_i2_base_vect_type)     :: x
      integer(psb_i2pk_), allocatable  :: xs(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_i2ovrl_save_vect
    module subroutine  psi_i2ovrl_save_multivect(x,xs,desc_a,info)
      class(psb_i2_base_multivect_type)     :: x
      integer(psb_i2pk_), allocatable  :: xs(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_i2ovrl_save_multivect
  end interface psi_ovrl_save

  interface psi_ovrl_restore
    module subroutine  psi_i2ovrl_restr_vect(x,xs,desc_a,info)
      class(psb_i2_base_vect_type)     :: x
      integer(psb_i2pk_)               :: xs(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_i2ovrl_restr_vect
    module subroutine  psi_i2ovrl_restr_multivect(x,xs,desc_a,info)
      class(psb_i2_base_multivect_type)     :: x
      integer(psb_i2pk_)               :: xs(:,:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_i2ovrl_restr_multivect
  end interface psi_ovrl_restore

end module psi_i2_comm_v_mod

