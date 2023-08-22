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
module psi_c_serial_mod
  use psb_const_mod, only :  psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_, psb_spk_

  interface psb_gelp 
    ! 2-D version
    subroutine psb_m_cgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_mpk_, psb_spk_
      implicit none
      complex(psb_spk_), intent(inout)     ::  x(:,:)
      integer(psb_mpk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_m_cgelp
    subroutine psb_m_cgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_mpk_,psb_spk_
      implicit none
      complex(psb_spk_), intent(inout)     ::  x(:)
      integer(psb_mpk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_m_cgelpv
    subroutine psb_e_cgelp(trans,iperm,x,info)
      import :: psb_ipk_, psb_epk_, psb_spk_
      implicit none
      complex(psb_spk_), intent(inout)     ::  x(:,:)
      integer(psb_epk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_e_cgelp
    subroutine psb_e_cgelpv(trans,iperm,x,info)
      import :: psb_ipk_, psb_epk_, psb_spk_
      implicit none
      complex(psb_spk_), intent(inout)     ::  x(:)
      integer(psb_epk_), intent(in)      ::  iperm(:)
      integer(psb_ipk_), intent(out)     ::  info
      character, intent(in)              :: trans
    end subroutine psb_e_cgelpv
  end interface psb_gelp

  interface psb_geaxpby
    subroutine psi_caxpby(m,n,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)      :: m, n
      complex(psb_spk_), intent (in)       ::  x(:,:)
      complex(psb_spk_), intent (inout)    ::  y(:,:)
      complex(psb_spk_), intent (in)       ::  alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_caxpby
    subroutine psi_caxpbyv(m,alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)      :: m
      complex(psb_spk_), intent (in)       ::  x(:)
      complex(psb_spk_), intent (inout)    ::  y(:)
      complex(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_caxpbyv
    subroutine psi_caxpbyv2(m,alpha, x, beta, y, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)      :: m
      complex(psb_spk_), intent (in)       ::  x(:)
      complex(psb_spk_), intent (in)       ::  y(:)
      complex(psb_spk_), intent (inout)    ::  z(:)
      complex(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psi_caxpbyv2
  end interface psb_geaxpby

  interface psi_gth
    subroutine psi_cgthmv(n,k,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_spk_) :: x(:,:), y(:),alpha,beta
    end subroutine psi_cgthmv
    subroutine psi_cgthv(n,idx,alpha,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      complex(psb_spk_) :: x(:), y(:),alpha,beta
    end subroutine psi_cgthv
    subroutine psi_cgthzmv(n,k,idx,x,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_spk_) :: x(:,:), y(:)

    end subroutine psi_cgthzmv
    subroutine psi_cgthzmm(n,k,idx,x,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_spk_) :: x(:,:), y(:,:)

    end subroutine psi_cgthzmm
    subroutine psi_cgthzv(n,idx,x,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, idx(:)
      complex(psb_spk_) :: x(:), y(:)
    end subroutine psi_cgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_csctmm(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_spk_) :: beta, x(:,:), y(:,:)
    end subroutine psi_csctmm
    subroutine psi_csctmv(n,k,idx,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_) :: n, k, idx(:)
      complex(psb_spk_) :: beta, x(:), y(:,:)
    end subroutine psi_csctmv
    subroutine psi_csctv(n,idx,x,beta,y)
      import :: psb_ipk_, psb_spk_
      implicit none

      integer(psb_ipk_) :: n, idx(:)
      complex(psb_spk_) :: beta, x(:), y(:)
    end subroutine psi_csctv
  end interface psi_sct

  interface psi_exscan
    subroutine psi_c_exscanv(n,x,info,shift)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)      :: n
      complex(psb_spk_), intent (inout)    :: x(:)
      integer(psb_ipk_), intent(out)     :: info
      complex(psb_spk_), intent(in), optional :: shift
    end subroutine psi_c_exscanv
  end interface psi_exscan
  
end module psi_c_serial_mod
