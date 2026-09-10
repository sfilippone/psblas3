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
module psi_s_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_, psb_spk_

  interface psb_gelp 
    ! 2-D version
    subroutine psb_m_sgelp(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_mpk_), intent(in)   :: iperm(:)
      real(psb_spk_), intent(inout)   :: x(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_m_sgelp

    subroutine psb_m_sgelpv(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_mpk_), intent(in)   :: iperm(:)
      real(psb_spk_), intent(inout)   :: x(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_m_sgelpv

    subroutine psb_e_sgelp(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_epk_), intent(in)   :: iperm(:)
      real(psb_spk_), intent(inout)   :: x(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_e_sgelp

    subroutine psb_e_sgelpv(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_epk_), intent(in)   :: iperm(:)
      real(psb_spk_), intent(inout)   :: x(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_e_sgelpv
  end interface psb_gelp

  interface psb_geaxpby
    subroutine psi_saxpby(m, n, alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: x(:, :)
      real(psb_spk_), intent(inout)   :: y(:, :)
      real(psb_spk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_saxpby

    subroutine psi_saxpby2(m, n, alpha, x, beta, y, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: x(:, :)
      real(psb_spk_), intent(in)      :: y(:, :)
      real(psb_spk_), intent(inout)   :: z(:, :)
      real(psb_spk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_saxpby2

    subroutine psi_saxpby3(m, n, alpha, x, beta, y, gamma, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: x(:, :)
      real(psb_spk_), intent(in)      :: y(:, :)
      real(psb_spk_), intent(inout)   :: z(:, :)
      real(psb_spk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_saxpby3

    subroutine psi_saxpbyv(m, alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(inout)   :: y(:)
      real(psb_spk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_saxpbyv

    subroutine psi_saxpbyv2(m, alpha, x, beta, y, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(in)      :: y(:)
      real(psb_spk_), intent(inout)   :: z(:)
      real(psb_spk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_saxpbyv2

    subroutine psi_saxpbyv3(m, alpha, x, beta, y, gamma, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(in)      :: y(:)
      real(psb_spk_), intent(inout)   :: z(:)
      real(psb_spk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_saxpbyv3

    subroutine psi_saxpbyv3_out(m, alpha, x, beta, y, gamma, z, w, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(in)      :: y(:)
      real(psb_spk_), intent(in)      :: z(:)
      real(psb_spk_), intent(inout)   :: w(:)
      real(psb_spk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_saxpbyv3_out
    
    subroutine psi_saxpbymvc(m, n, alpha, x, beta, y, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(inout)   :: y(:, :)
      real(psb_spk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_saxpbymvc
  end interface psb_geaxpby

  interface psb_gemlt
    subroutine psi_smlt(m, n, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:, :)
      real(psb_spk_), intent(inout)   :: y(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_smlt

    subroutine psi_smlt2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:, :)
      real(psb_spk_), intent(in)      :: y(:, :)
      real(psb_spk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_smlt2

    subroutine psi_smltv(m, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(inout)   :: y(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_smltv

    subroutine psi_smltv2(m, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(in)      :: y(:)
      real(psb_spk_), intent(inout)   :: z(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_smltv2

    subroutine psi_smltx(m, n, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(inout)   :: y(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_smltx

    subroutine psi_smltx2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(in)      :: y(:, :)
      real(psb_spk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_smltx2

    subroutine psi_smlte2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_spk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(in)      :: y(:)
      real(psb_spk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_smlte2
  end interface psb_gemlt

  interface psi_upd_xyz
    subroutine psi_s_upd_xyz(m, alpha, beta, gamma, delta, x, y, z, info)
      import
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(inout)   :: y(:)
      real(psb_spk_), intent(inout)   :: z(:)
      real(psb_spk_), intent(in)      :: alpha, beta, gamma, delta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_s_upd_xyz
  end interface psi_upd_xyz
  
  interface psi_xyzw
    subroutine psi_sxyzw(m, a, b, c, d, e, f, x, y, z, w, info)
      import
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(inout)   :: y(:)
      real(psb_spk_), intent(inout)   :: z(:)
      real(psb_spk_), intent(inout)   :: w(:)
      real(psb_spk_), intent(in)      :: a, b, c, d, e, f
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_sxyzw
  end interface psi_xyzw
  
  interface psi_gth
    subroutine psi_sgthmv(n, k, idx, alpha, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_)    :: alpha, x(:, :), beta, y(:)
    end subroutine psi_sgthmv

    subroutine psi_sgthv(n, idx, alpha, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_)    :: alpha, x(:), beta, y(:)
    end subroutine psi_sgthv

    subroutine psi_sgthzmv(n, k, idx, x, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_)    :: x(:, :), y(:)
    end subroutine psi_sgthzmv

    subroutine psi_sgthzmm(n, k, idx, x, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_)    :: x(:, :), y(:, :)
    end subroutine psi_sgthzmm

    subroutine psi_sgthzv(n, idx, x, y)
      import
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_)    :: x(:), y(:)
    end subroutine psi_sgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_ssctmm(n, k, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_)    :: x(:, :), beta, y(:, :)
    end subroutine psi_ssctmm

    subroutine psi_ssctmv(n, k, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_)    :: x(:), beta, y(:, :)
    end subroutine psi_ssctmv

    subroutine psi_ssctv(n, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_)    :: x(:), beta, y(:)
    end subroutine psi_ssctv
  end interface psi_sct

  interface psi_exscan
    subroutine psi_s_exscanv(n, x, info, shift)
      import
      implicit none
      integer(psb_ipk_), intent(in) :: n
      real(psb_spk_), intent(inout) :: x(:)
      integer(psb_ipk_), intent(out) :: info
      real(psb_spk_), intent(in), optional :: shift
    end subroutine psi_s_exscanv
  end interface psi_exscan

  !Dispach of axpy-like operations.
  interface get_axpbylike_code
    module procedure get_axpbylike_code1
    module procedure get_axpbylike_code2
    module procedure get_axpbylike_code3
  end interface get_axpbylike_code

contains
  function get_axpbylike_code1(var) result(code)
    real(psb_spk_), intent(in) :: var
    integer(psb_ipk_) :: code

    code = 0_psb_ipk_
    if(var ==  sone) code = 1_psb_ipk_
    if(var == szero) code = 2_psb_ipk_
    if(var == -sone) code = 3_psb_ipk_
  end function get_axpbylike_code1

  function get_axpbylike_code2(alpha, beta) result(code)
    real(psb_spk_), intent(in) :: alpha, beta
    integer(psb_ipk_) :: code

    code = get_axpbylike_code1(alpha) &
            + ishft(get_axpbylike_code1(beta), 2)
  end function get_axpbylike_code2

  function get_axpbylike_code3(alpha, beta, gamma) result(code)
    real(psb_spk_), intent(in) :: alpha, beta, gamma
    integer(psb_ipk_) :: code

    code = get_axpbylike_code1(alpha) &
            + ishft(get_axpbylike_code1(beta), 2) &
            + ishft(get_axpbylike_code1(gamma), 4)
  end function get_axpbylike_code3
end module psi_s_serial_mod