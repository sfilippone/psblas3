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
module psi_z_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_, psb_dpk_

  interface psb_gelp 
    ! 2-D version
    subroutine psb_m_zgelp(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_mpk_), intent(in)   :: iperm(:)
      complex(psb_dpk_), intent(inout)   :: x(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_m_zgelp

    subroutine psb_m_zgelpv(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_mpk_), intent(in)   :: iperm(:)
      complex(psb_dpk_), intent(inout)   :: x(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_m_zgelpv

    subroutine psb_e_zgelp(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_epk_), intent(in)   :: iperm(:)
      complex(psb_dpk_), intent(inout)   :: x(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_e_zgelp

    subroutine psb_e_zgelpv(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_epk_), intent(in)   :: iperm(:)
      complex(psb_dpk_), intent(inout)   :: x(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_e_zgelpv
  end interface psb_gelp

  interface psb_geaxpby
    subroutine psi_zaxpby(m, n, alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: x(:, :)
      complex(psb_dpk_), intent(inout)   :: y(:, :)
      complex(psb_dpk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zaxpby

    subroutine psi_zaxpby2(m, n, alpha, x, beta, y, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: x(:, :)
      complex(psb_dpk_), intent(in)      :: y(:, :)
      complex(psb_dpk_), intent(inout)   :: z(:, :)
      complex(psb_dpk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zaxpby2

    subroutine psi_zaxpby3(m, n, alpha, x, beta, y, gamma, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: x(:, :)
      complex(psb_dpk_), intent(in)      :: y(:, :)
      complex(psb_dpk_), intent(inout)   :: z(:, :)
      complex(psb_dpk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zaxpby3

    subroutine psi_zaxpbyv(m, alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(inout)   :: y(:)
      complex(psb_dpk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zaxpbyv

    subroutine psi_zaxpbyv2(m, alpha, x, beta, y, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(in)      :: y(:)
      complex(psb_dpk_), intent(inout)   :: z(:)
      complex(psb_dpk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zaxpbyv2

    subroutine psi_zaxpbyv3(m, alpha, x, beta, y, gamma, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(in)      :: y(:)
      complex(psb_dpk_), intent(inout)   :: z(:)
      complex(psb_dpk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zaxpbyv3

    subroutine psi_zaxpbyv3_out(m, alpha, x, beta, y, gamma, z, w, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(in)      :: y(:)
      complex(psb_dpk_), intent(in)      :: z(:)
      complex(psb_dpk_), intent(inout)   :: w(:)
      complex(psb_dpk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zaxpbyv3_out
    
    subroutine psi_zaxpbymvc(m, n, alpha, x, beta, y, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(inout)   :: y(:, :)
      complex(psb_dpk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zaxpbymvc
  end interface psb_geaxpby

  interface psb_gemlt
    subroutine psi_zmlt(m, n, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: alpha, beta
      complex(psb_dpk_), intent(in)      :: x(:, :)
      complex(psb_dpk_), intent(inout)   :: y(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zmlt

    subroutine psi_zmlt2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: alpha, beta
      complex(psb_dpk_), intent(in)      :: x(:, :)
      complex(psb_dpk_), intent(in)      :: y(:, :)
      complex(psb_dpk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zmlt2

    subroutine psi_zmltv(m, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      complex(psb_dpk_), intent(in)      :: alpha, beta
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(inout)   :: y(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zmltv

    subroutine psi_zmltv2(m, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      complex(psb_dpk_), intent(in)      :: alpha, beta
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(in)      :: y(:)
      complex(psb_dpk_), intent(inout)   :: z(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zmltv2

    subroutine psi_zmltx(m, n, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: alpha, beta
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(inout)   :: y(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zmltx

    subroutine psi_zmltx2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: alpha, beta
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(in)      :: y(:, :)
      complex(psb_dpk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zmltx2

    subroutine psi_zmlte2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_dpk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      complex(psb_dpk_), intent(in)      :: alpha, beta
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(in)      :: y(:)
      complex(psb_dpk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zmlte2
  end interface psb_gemlt

  interface psi_upd_xyz
    subroutine psi_z_upd_xyz(m, alpha, beta, gamma, delta, x, y, z, info)
      import
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(inout)   :: y(:)
      complex(psb_dpk_), intent(inout)   :: z(:)
      complex(psb_dpk_), intent(in)      :: alpha, beta, gamma, delta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_z_upd_xyz
  end interface psi_upd_xyz
  
  interface psi_xyzw
    subroutine psi_zxyzw(m, a, b, c, d, e, f, x, y, z, w, info)
      import
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      complex(psb_dpk_), intent(in)      :: x(:)
      complex(psb_dpk_), intent(inout)   :: y(:)
      complex(psb_dpk_), intent(inout)   :: z(:)
      complex(psb_dpk_), intent(inout)   :: w(:)
      complex(psb_dpk_), intent(in)      :: a, b, c, d, e, f
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_zxyzw
  end interface psi_xyzw
  
  interface psi_gth
    subroutine psi_zgthmv(n, k, idx, alpha, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      complex(psb_dpk_)    :: alpha, x(:, :), beta, y(:)
    end subroutine psi_zgthmv

    subroutine psi_zgthv(n, idx, alpha, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_dpk_)    :: alpha, x(:), beta, y(:)
    end subroutine psi_zgthv

    subroutine psi_zgthzmv(n, k, idx, x, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      complex(psb_dpk_)    :: x(:, :), y(:)
    end subroutine psi_zgthzmv

    subroutine psi_zgthzmm(n, k, idx, x, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      complex(psb_dpk_)    :: x(:, :), y(:, :)
    end subroutine psi_zgthzmm

    subroutine psi_zgthzv(n, idx, x, y)
      import
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_dpk_)    :: x(:), y(:)
    end subroutine psi_zgthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_zsctmm(n, k, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      complex(psb_dpk_)    :: x(:, :), beta, y(:, :)
    end subroutine psi_zsctmm

    subroutine psi_zsctmv(n, k, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      complex(psb_dpk_)    :: x(:), beta, y(:, :)
    end subroutine psi_zsctmv

    subroutine psi_zsctv(n, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_dpk_)    :: x(:), beta, y(:)
    end subroutine psi_zsctv
  end interface psi_sct

  interface psi_exscan
    subroutine psi_z_exscanv(n, x, info, shift)
      import
      implicit none
      integer(psb_ipk_), intent(in) :: n
      complex(psb_dpk_), intent(inout) :: x(:)
      integer(psb_ipk_), intent(out) :: info
      complex(psb_dpk_), intent(in), optional :: shift
    end subroutine psi_z_exscanv
  end interface psi_exscan

  !Dispach of axpy-like operations.
  interface get_axpbylike_code
    module procedure get_axpbylike_code1
    module procedure get_axpbylike_code2
    module procedure get_axpbylike_code3
  end interface get_axpbylike_code

contains
  function get_axpbylike_code1(var) result(code)
    complex(psb_dpk_), intent(in) :: var
    integer(psb_ipk_) :: code

    code = 0_psb_ipk_
    if(var ==  zone) code = 1_psb_ipk_
    if(var == zzero) code = 2_psb_ipk_
    if(var == -zone) code = 3_psb_ipk_
  end function get_axpbylike_code1

  function get_axpbylike_code2(alpha, beta) result(code)
    complex(psb_dpk_), intent(in) :: alpha, beta
    integer(psb_ipk_) :: code

    code = get_axpbylike_code1(alpha) &
            + ishft(get_axpbylike_code1(beta), 2)
  end function get_axpbylike_code2

  function get_axpbylike_code3(alpha, beta, gamma) result(code)
    complex(psb_dpk_), intent(in) :: alpha, beta, gamma
    integer(psb_ipk_) :: code

    code = get_axpbylike_code1(alpha) &
            + ishft(get_axpbylike_code1(beta), 2) &
            + ishft(get_axpbylike_code1(gamma), 4)
  end function get_axpbylike_code3
end module psi_z_serial_mod