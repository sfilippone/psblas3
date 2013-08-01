!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
module psi_d_mod
  use psb_desc_mod, only : psb_desc_type, psb_ipk_, psb_dpk_, psb_i_base_vect_type
  use psb_d_vect_mod, only : psb_d_base_vect_type 

  interface psi_swapdata
    subroutine psi_dswapdatam(flag,n,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)         :: flag, n
      integer(psb_ipk_), intent(out)        :: info
      real(psb_dpk_)              :: y(:,:), beta
      real(psb_dpk_),target       :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_dswapdatam
    subroutine psi_dswapdatav(flag,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      real(psb_dpk_)              :: y(:), beta 
      real(psb_dpk_),target       :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_dswapdatav
    subroutine psi_dswapdata_vect(flag,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_d_base_vect_type) :: y
      real(psb_dpk_)              :: beta 
      real(psb_dpk_),target       :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_dswapdata_vect
    subroutine psi_dswapidxm(ictxt,icomm,flag,n,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)   :: ictxt,icomm,flag, n
      integer(psb_ipk_), intent(out)  :: info
      real(psb_dpk_)        :: y(:,:), beta
      real(psb_dpk_),target :: work(:)
      integer(psb_ipk_), intent(in)   :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dswapidxm
    subroutine psi_dswapidxv(ictxt,icomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)   :: ictxt,icomm,flag
      integer(psb_ipk_), intent(out)  :: info
      real(psb_dpk_)        :: y(:), beta
      real(psb_dpk_),target :: work(:)
      integer(psb_ipk_), intent(in)   :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dswapidxv
    subroutine psi_dswapidx_vect(ictxt,icomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)   :: ictxt,icomm,flag
      integer(psb_ipk_), intent(out)  :: info
      class(psb_d_base_vect_type) :: y
      real(psb_dpk_)        :: beta
      real(psb_dpk_),target :: work(:)
      integer(psb_ipk_), intent(in)   :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dswapidx_vect
    subroutine psi_dswap_vidx_vect(iictxt,iicomm,flag,beta,y,idx,totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type, psb_i_base_vect_type
      integer(psb_ipk_), intent(in)           :: iictxt,iicomm,flag
      integer(psb_ipk_), intent(out)          :: info
      class(psb_d_base_vect_type)             :: y
      real(psb_dpk_)                          :: beta
      real(psb_dpk_), target                  :: work(:)
      class(psb_i_base_vect_type), intent(in) :: idx
      integer(psb_ipk_), intent(in)           :: totxch,totsnd, totrcv
    end subroutine psi_dswap_vidx_vect
  end interface


  interface psi_swaptran
    subroutine psi_dswaptranm(flag,n,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)         :: flag, n
      integer(psb_ipk_), intent(out)        :: info
      real(psb_dpk_)              :: y(:,:), beta
      real(psb_dpk_),target       :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_dswaptranm
    subroutine psi_dswaptranv(flag,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      real(psb_dpk_)              :: y(:), beta
      real(psb_dpk_),target       :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_dswaptranv
    subroutine psi_dswaptran_vect(flag,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_d_base_vect_type) :: y
      real(psb_dpk_)              :: beta
      real(psb_dpk_),target       :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_dswaptran_vect
    subroutine psi_dtranidxm(ictxt,icomm,flag,n,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)   :: ictxt,icomm,flag, n
      integer(psb_ipk_), intent(out)  :: info
      real(psb_dpk_)        :: y(:,:), beta
      real(psb_dpk_),target :: work(:)
      integer(psb_ipk_), intent(in)   :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dtranidxm
    subroutine psi_dtranidxv(ictxt,icomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)   :: ictxt,icomm,flag
      integer(psb_ipk_), intent(out)  :: info
      real(psb_dpk_)        :: y(:), beta
      real(psb_dpk_),target :: work(:)
      integer(psb_ipk_), intent(in)   :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dtranidxv
    subroutine psi_dtranidx_vect(ictxt,icomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      integer(psb_ipk_), intent(in)   :: ictxt,icomm,flag
      integer(psb_ipk_), intent(out)  :: info
      class(psb_d_base_vect_type) :: y
      real(psb_dpk_)        :: beta
      real(psb_dpk_),target :: work(:)
      integer(psb_ipk_), intent(in)   :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dtranidx_vect
  end interface

  interface psi_ovrl_upd
    subroutine  psi_dovrl_updr1(x,desc_a,update,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      real(psb_dpk_), intent(inout), target :: x(:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(in)                     :: update
      integer(psb_ipk_), intent(out)                    :: info
    end subroutine psi_dovrl_updr1
    subroutine  psi_dovrl_updr2(x,desc_a,update,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      real(psb_dpk_), intent(inout), target :: x(:,:)
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(in)                     :: update
      integer(psb_ipk_), intent(out)                    :: info
    end subroutine psi_dovrl_updr2
    subroutine  psi_dovrl_upd_vect(x,desc_a,update,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      class(psb_d_base_vect_type)       :: x
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(in)               :: update
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psi_dovrl_upd_vect
  end interface

  interface psi_ovrl_save
    subroutine  psi_dovrl_saver1(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      real(psb_dpk_), intent(inout)  :: x(:)
      real(psb_dpk_), allocatable    :: xs(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_dovrl_saver1
    subroutine  psi_dovrl_saver2(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      real(psb_dpk_), intent(inout)  :: x(:,:)
      real(psb_dpk_), allocatable    :: xs(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_dovrl_saver2
    subroutine  psi_dovrl_save_vect(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      class(psb_d_base_vect_type)     :: x
      real(psb_dpk_), allocatable     :: xs(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_dovrl_save_vect
  end interface

  interface psi_ovrl_restore
    subroutine  psi_dovrl_restrr1(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      real(psb_dpk_), intent(inout)  :: x(:)
      real(psb_dpk_)                 :: xs(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_dovrl_restrr1
    subroutine  psi_dovrl_restrr2(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      real(psb_dpk_), intent(inout)  :: x(:,:)
      real(psb_dpk_)                 :: xs(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_dovrl_restrr2
    subroutine  psi_dovrl_restr_vect(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_, psb_dpk_, psb_d_base_vect_type
      class(psb_d_base_vect_type)     :: x
      real(psb_dpk_)                  :: xs(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_dovrl_restr_vect
  end interface

end module psi_d_mod

