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
module psi_i2_comm_a_mod
  use psi_penv_mod, only : psb_ctxt_type
  use psb_desc_mod, only : psb_desc_type, psb_mpk_, psb_ipk_, psb_epk_

  interface psi_swapdata
    module subroutine psi_i2swapdatam(flag,n,beta,y,desc_a,work,info,data)
      integer(psb_mpk_), intent(in)         :: n
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_i2pk_)           :: y(:,:), beta
      integer(psb_i2pk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_i2swapdatam
    module subroutine psi_i2swapdatav(flag,beta,y,desc_a,work,info,data)
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_i2pk_)           :: y(:), beta 
      integer(psb_i2pk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_i2swapdatav
    module subroutine psi_i2swapidxm(ctxt,flag,n,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      type(psb_ctxt_type), intent(in) :: ctxt
      integer(psb_mpk_), intent(in)   :: n
      integer(psb_ipk_), intent(in)   :: flag
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_i2pk_)        :: y(:,:), beta
      integer(psb_i2pk_),target :: work(:)
      integer(psb_ipk_), intent(in)      :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_i2swapidxm
    module subroutine psi_i2swapidxv(ctxt,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      type(psb_ctxt_type), intent(in) :: ctxt
      integer(psb_ipk_), intent(in)   :: flag
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_i2pk_)        :: y(:), beta
      integer(psb_i2pk_),target :: work(:)
      integer(psb_ipk_), intent(in)      :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_i2swapidxv
  end interface psi_swapdata


  interface psi_swaptran
    module subroutine psi_i2swaptranm(flag,n,beta,y,desc_a,work,info,data)
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_Mpk_), intent(in)         :: n
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_i2pk_)           :: y(:,:), beta
      integer(psb_i2pk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_i2swaptranm
    module subroutine psi_i2swaptranv(flag,beta,y,desc_a,work,info,data)
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_i2pk_)           :: y(:), beta
      integer(psb_i2pk_),target    :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_i2swaptranv
    module subroutine psi_i2tranidxm(ctxt,flag,n,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      type(psb_ctxt_type), intent(in) :: ctxt
      integer(psb_mpk_), intent(in)   :: n
      integer(psb_ipk_), intent(in)   :: flag
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_i2pk_)        :: y(:,:), beta
      integer(psb_i2pk_),target :: work(:)
      integer(psb_ipk_), intent(in)       :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_i2tranidxm
    module subroutine psi_i2tranidxv(ctxt,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      type(psb_ctxt_type), intent(in) :: ctxt
      integer(psb_ipk_), intent(in)   :: flag
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_i2pk_)        :: y(:), beta
      integer(psb_i2pk_),target :: work(:)
      integer(psb_ipk_), intent(in)      :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_i2tranidxv
  end interface psi_swaptran
 
  interface psi_ovrl_upd
    module subroutine  psi_i2ovrl_updr1(x,desc_a,update,info)
      integer(psb_i2pk_), intent(inout), target :: x(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(in)    :: update
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psi_i2ovrl_updr1
    module subroutine  psi_i2ovrl_updr2(x,desc_a,update,info)
      integer(psb_i2pk_), intent(inout), target :: x(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(in)      :: update
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_i2ovrl_updr2
  end interface psi_ovrl_upd

  interface psi_ovrl_save
    module subroutine  psi_i2ovrl_saver1(x,xs,desc_a,info)
      integer(psb_i2pk_), intent(inout) :: x(:)
      integer(psb_i2pk_), allocatable   :: xs(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psi_i2ovrl_saver1
    module subroutine  psi_i2ovrl_saver2(x,xs,desc_a,info)
      integer(psb_i2pk_), intent(inout) :: x(:,:)
      integer(psb_i2pk_), allocatable   :: xs(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psi_i2ovrl_saver2
  end interface psi_ovrl_save

  interface psi_ovrl_restore
    module subroutine  psi_i2ovrl_restrr1(x,xs,desc_a,info)
      integer(psb_i2pk_), intent(inout)  :: x(:)
      integer(psb_i2pk_)                 :: xs(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psi_i2ovrl_restrr1
    module subroutine  psi_i2ovrl_restrr2(x,xs,desc_a,info)
      integer(psb_i2pk_), intent(inout) :: x(:,:)
      integer(psb_i2pk_)                :: xs(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)   :: info
    end subroutine psi_i2ovrl_restrr2
  end interface psi_ovrl_restore

end module psi_i2_comm_a_mod

