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
module psi_i2_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_

  interface psb_gelp
    ! 2-D version
    subroutine psb_i2gelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none 
      integer(psb_i2pk_), intent(inout)     ::  x(:,:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_i2gelp
    subroutine psb_i2gelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none 
      integer(psb_i2pk_), intent(inout)     ::  x(:)
      integer(psb_ipk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_i2gelpv
  end interface psb_gelp
  
  interface psb_geaxpby 
    subroutine psi_i2axpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m, n
      integer(psb_i2pk_), intent (in)       ::  x(:,:)
      integer(psb_i2pk_), intent (inout)    ::  y(:,:)
      integer(psb_i2pk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_i2axpby
    subroutine psi_i2axpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none 
      integer(psb_ipk_), intent(in)      :: m
      integer(psb_i2pk_), intent (in)       ::  x(:)
      integer(psb_i2pk_), intent (inout)    ::  y(:)
      integer(psb_i2pk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_i2axpbyv
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_i2gthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_i2pk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_i2gthmv
    subroutine psi_i2gthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_i2pk_) :: x(:), y(:),alpha,beta
    end subroutine psi_i2gthv
    subroutine psi_i2gthzmv(n,k,idx,x,y)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_i2pk_) :: x(:,:), y(:)
      
    end subroutine psi_i2gthzmv
    subroutine psi_i2gthzmm(n,k,idx,x,y)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_i2pk_) :: x(:,:), y(:,:)
      
    end subroutine psi_i2gthzmm
    subroutine psi_i2gthzv(n,idx,x,y)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none 
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_i2pk_) :: x(:), y(:)
    end subroutine psi_i2gthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_i2sctmm(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_i2pk_) :: beta, x(:,:), y(:,:)
    end subroutine psi_i2sctmm
    subroutine psi_i2sctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      integer(psb_i2pk_) :: beta, x(:), y(:,:)
    end subroutine psi_i2sctmv
    subroutine psi_i2sctv(n,idx,x,beta,y)
      import :: psb_ipk_, psb_lpk_,psb_mpk_, psb_epk_
      implicit none
      
      integer(psb_ipk_) :: n, idx(:)
      integer(psb_i2pk_) :: beta, x(:), y(:)
    end subroutine psi_i2sctv
  end interface psi_sct

end module psi_i2_serial_mod
