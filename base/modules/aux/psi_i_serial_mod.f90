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
module psi_i_serial_mod
  use psb_const_mod, only : psb_ipk_

  interface psb_gelp
    ! 2-D version
    subroutine psb_igelp(trans,iperm,x,info)
      import :: psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(inout)     ::  x(:,:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_igelp
    subroutine psb_igelpv(trans,iperm,x,info)
      import :: psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(inout)     ::  x(:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_igelpv
  end interface psb_gelp
  
  interface psb_geaxpby 
    subroutine psi_iaxpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m, n
      integer(psb_ipk_), intent (in)       ::  x(:,:)
      integer(psb_ipk_), intent (inout)    ::  y(:,:)
      integer(psb_ipk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_iaxpby
    subroutine psi_iaxpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m
      integer(psb_ipk_), intent (in)       ::  x(:)
      integer(psb_ipk_), intent (inout)    ::  y(:)
      integer(psb_ipk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_iaxpbyv
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_igthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_ipk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_igthmv
    subroutine psi_igthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_ipk_) :: x(:), y(:),alpha,beta
    end subroutine psi_igthv
    subroutine psi_igthzmv(n,k,idx,x,y)
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_ipk_) :: x(:,:), y(:)
      
    end subroutine psi_igthzmv
    subroutine psi_igthzmm(n,k,idx,x,y)
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_ipk_) :: x(:,:), y(:,:)
      
    end subroutine psi_igthzmm
    subroutine psi_igthzv(n,idx,x,y)
      import :: psb_ipk_
      implicit none 
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_ipk_) :: x(:), y(:)
    end subroutine psi_igthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_isctmm(n,k,idx,x,beta,y)
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_ipk_) :: beta, x(:,:), y(:,:)
    end subroutine psi_isctmm
    subroutine psi_isctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_ipk_) :: beta, x(:), y(:,:)
    end subroutine psi_isctmv
    subroutine psi_isctv(n,idx,x,beta,y)
      import :: psb_ipk_
      implicit none
      
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_ipk_) :: beta, x(:), y(:)
    end subroutine psi_isctv
  end interface psi_sct

end module psi_i_serial_mod
