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
module psi_l_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_lpk_

  interface psb_gelp
    ! 2-D version
    subroutine psb_lgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_lpk_), intent(inout)     ::  x(:,:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_lgelp
    subroutine psb_lgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_lpk_), intent(inout)     ::  x(:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_lgelpv
  end interface psb_gelp
  
  interface psb_geaxpby 
    subroutine psi_laxpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m, n
      integer(psb_lpk_), intent (in)       ::  x(:,:)
      integer(psb_lpk_), intent (inout)    ::  y(:,:)
      integer(psb_lpk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_laxpby
    subroutine psi_laxpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m
      integer(psb_lpk_), intent (in)       ::  x(:)
      integer(psb_lpk_), intent (inout)    ::  y(:)
      integer(psb_lpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_laxpbyv
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_lgthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_lpk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_lgthmv
    subroutine psi_lgthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_lpk_) :: x(:), y(:),alpha,beta
    end subroutine psi_lgthv
    subroutine psi_lgthzmv(n,k,idx,x,y)
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_lpk_) :: x(:,:), y(:)
      
    end subroutine psi_lgthzmv
    subroutine psi_lgthzmm(n,k,idx,x,y)
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_lpk_) :: x(:,:), y(:,:)
      
    end subroutine psi_lgthzmm
    subroutine psi_lgthzv(n,idx,x,y)
      import :: psb_ipk_, psb_lpk_
      implicit none 
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_lpk_) :: x(:), y(:)
    end subroutine psi_lgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_lsctmm(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_lpk_) :: beta, x(:,:), y(:,:)
    end subroutine psi_lsctmm
    subroutine psi_lsctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_lpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_lpk_) :: beta, x(:), y(:,:)
    end subroutine psi_lsctmv
    subroutine psi_lsctv(n,idx,x,beta,y)
      import :: psb_ipk_, psb_lpk_
      implicit none
      
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_lpk_) :: beta, x(:), y(:)
    end subroutine psi_lsctv
  end interface psi_sct

end module psi_l_serial_mod
