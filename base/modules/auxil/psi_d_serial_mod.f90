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
module psi_d_serial_mod
  use psb_const_mod, only :  psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_, psb_dpk_

  interface psb_gelp 
    ! 2-D version
    subroutine psb_m_dgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_mpk_, psb_dpk_
      implicit none
      real(psb_dpk_), intent(inout)     ::  x(:,:)
      integer(psb_mpk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_m_dgelp
    subroutine psb_m_dgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_mpk_,psb_dpk_
      implicit none
      real(psb_dpk_), intent(inout)     ::  x(:)
      integer(psb_mpk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_m_dgelpv
    subroutine psb_e_dgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_epk_, psb_dpk_
      implicit none
      real(psb_dpk_), intent(inout)     ::  x(:,:)
      integer(psb_epk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_e_dgelp
    subroutine psb_e_dgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_epk_, psb_dpk_
      implicit none
      real(psb_dpk_), intent(inout)     ::  x(:)
      integer(psb_epk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_e_dgelpv
  end interface psb_gelp

  interface psb_geaxpby
    subroutine psi_daxpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)      :: m, n
      real(psb_dpk_), intent (in)       ::  x(:,:)
      real(psb_dpk_), intent (inout)    ::  y(:,:)
      real(psb_dpk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_daxpby
    subroutine psi_daxpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)      :: m
      real(psb_dpk_), intent (in)       ::  x(:)
      real(psb_dpk_), intent (inout)    ::  y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_daxpbyv
    subroutine psi_daxpbyv2(m,alpha, x, beta, y, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)      :: m
      real(psb_dpk_), intent (in)       ::  x(:)
      real(psb_dpk_), intent (in)       ::  y(:)
      real(psb_dpk_), intent (inout)    ::  z(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_daxpbyv2
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_dgthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_dpk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_dgthmv
    subroutine psi_dgthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      real(psb_dpk_) :: x(:), y(:),alpha,beta
    end subroutine psi_dgthv
    subroutine psi_dgthzmv(n,k,idx,x,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_dpk_) :: x(:,:), y(:)

    end subroutine psi_dgthzmv
    subroutine psi_dgthzmm(n,k,idx,x,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_dpk_) :: x(:,:), y(:,:)

    end subroutine psi_dgthzmm
    subroutine psi_dgthzv(n,idx,x,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      real(psb_dpk_) :: x(:), y(:)
    end subroutine psi_dgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_dsctmm(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_dpk_) :: beta, x(:,:), y(:,:)
    end subroutine psi_dsctmm
    subroutine psi_dsctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      real(psb_dpk_) :: beta, x(:), y(:,:)
    end subroutine psi_dsctmv
    subroutine psi_dsctv(n,idx,x,beta,y)
      import :: psb_ipk_, psb_dpk_
      implicit none

      integer(psb_ipk_) :: n, idx(:)
      real(psb_dpk_) :: beta, x(:), y(:)
    end subroutine psi_dsctv
  end interface psi_sct

end module psi_d_serial_mod
