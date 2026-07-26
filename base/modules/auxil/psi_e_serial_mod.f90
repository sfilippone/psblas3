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
module psi_e_serial_mod
  use psb_const_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_, psb_i2pk_

  interface psb_gelp 
    ! 2-D version
    subroutine psb_m_egelp(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_mpk_), intent(in)   :: iperm(:)
      integer(psb_epk_), intent(inout)   :: x(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_m_egelp

    subroutine psb_m_egelpv(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_mpk_), intent(in)   :: iperm(:)
      integer(psb_epk_), intent(inout)   :: x(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_m_egelpv

    subroutine psb_e_egelp(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_epk_), intent(in)   :: iperm(:)
      integer(psb_epk_), intent(inout)   :: x(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_e_egelp

    subroutine psb_e_egelpv(trans, iperm, x, info)
      import
      implicit none
      character, intent(in)           :: trans
      integer(psb_epk_), intent(in)   :: iperm(:)
      integer(psb_epk_), intent(inout)   :: x(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_e_egelpv
  end interface psb_gelp

  interface psb_geaxpby
    subroutine psi_eaxpby(m, n, alpha, x, beta, y, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      integer(psb_epk_), intent(in)      :: x(:, :)
      integer(psb_epk_), intent(inout)   :: y(:, :)
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_eaxpby

    subroutine psi_eaxpby2(m, n, alpha, x, beta, y, z, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      integer(psb_epk_), intent(in)      :: x(:, :)
      integer(psb_epk_), intent(in)      :: y(:, :)
      integer(psb_epk_), intent(inout)   :: z(:, :)
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_eaxpby2

    subroutine psi_eaxpby3(m, n, alpha, x, beta, y, gamma, z, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      integer(psb_epk_), intent(in)      :: x(:, :)
      integer(psb_epk_), intent(in)      :: y(:, :)
      integer(psb_epk_), intent(inout)   :: z(:, :)
      integer(psb_epk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_eaxpby3

    subroutine psi_eaxpbyv(m, alpha, x, beta, y, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(inout)   :: y(:)
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_eaxpbyv

    subroutine psi_eaxpbyv2(m, alpha, x, beta, y, z, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(in)      :: y(:)
      integer(psb_epk_), intent(inout)   :: z(:)
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_eaxpbyv2

    subroutine psi_eaxpbyv3(m, alpha, x, beta, y, gamma, z, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(in)      :: y(:)
      integer(psb_epk_), intent(inout)   :: z(:)
      integer(psb_epk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_eaxpbyv3

    subroutine psi_eaxpbyv3_out(m, alpha, x, beta, y, gamma, z, w, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(in)      :: y(:)
      integer(psb_epk_), intent(in)      :: z(:)
      integer(psb_epk_), intent(inout)   :: w(:)
      integer(psb_epk_), intent(in)      :: alpha, beta, gamma
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_eaxpbyv3_out
  end interface psb_geaxpby

  interface psb_gemlt
    subroutine psi_emlt(m, n, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_epk_), intent(in)      :: x(:, :)
      integer(psb_epk_), intent(inout)   :: y(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_emlt

    subroutine psi_emlt2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_epk_), intent(in)      :: x(:, :)
      integer(psb_epk_), intent(in)      :: y(:, :)
      integer(psb_epk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_emlt2

    subroutine psi_emltv(m, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(inout)   :: y(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_emltv

    subroutine psi_emltv2(m, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(in)      :: y(:)
      integer(psb_epk_), intent(inout)   :: z(:)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_emltv2

    subroutine psi_emltx(m, n, alpha, x, y, beta, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(inout)   :: y(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_emltx

    subroutine psi_emltx2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(in)      :: y(:, :)
      integer(psb_epk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_emltx2

    subroutine psi_emlte2(m, n, alpha, x, y, beta, z, info)
      import :: psb_ipk_, psb_epk_
      implicit none
      integer(psb_ipk_), intent(in)   :: m, n
      integer(psb_epk_), intent(in)      :: alpha, beta
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(in)      :: y(:)
      integer(psb_epk_), intent(inout)   :: z(:, :)
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_emlte2
  end interface psb_gemlt

  interface psi_upd_xyz
    subroutine psi_e_upd_xyz(m, alpha, beta, gamma, delta, x, y, z, info)
      import
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(inout)   :: y(:)
      integer(psb_epk_), intent(inout)   :: z(:)
      integer(psb_epk_), intent(in)      :: alpha, beta, gamma, delta
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_e_upd_xyz
  end interface psi_upd_xyz
  
  interface psi_xyzw
    subroutine psi_exyzw(m, a, b, c, d, e, f, x, y, z, w, info)
      import
      implicit none
      integer(psb_ipk_), intent(in)   :: m
      integer(psb_epk_), intent(in)      :: x(:)
      integer(psb_epk_), intent(inout)   :: y(:)
      integer(psb_epk_), intent(inout)   :: z(:)
      integer(psb_epk_), intent(inout)   :: w(:)
      integer(psb_epk_), intent(in)      :: a, b, c, d, e, f
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psi_exyzw
  end interface psi_xyzw
  
  interface psi_gth
    subroutine psi_egthmv(n, k, idx, alpha, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      integer(psb_epk_)    :: alpha, x(:, :), beta, y(:)
    end subroutine psi_egthmv

    subroutine psi_egthv(n, idx, alpha, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_epk_)    :: alpha, x(:), beta, y(:)
    end subroutine psi_egthv

    subroutine psi_egthzmv(n, k, idx, x, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      integer(psb_epk_)    :: x(:, :), y(:)
    end subroutine psi_egthzmv

    subroutine psi_egthzmm(n, k, idx, x, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      integer(psb_epk_)    :: x(:, :), y(:, :)
    end subroutine psi_egthzmm

    subroutine psi_egthzv(n, idx, x, y)
      import
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_epk_)    :: x(:), y(:)
    end subroutine psi_egthzv
  end interface psi_gth

  interface psi_sct
    subroutine psi_esctmm(n, k, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      integer(psb_epk_)    :: x(:, :), beta, y(:, :)
    end subroutine psi_esctmm

    subroutine psi_esctmv(n, k, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n, k
      integer(psb_ipk_) :: idx(:)
      integer(psb_epk_)    :: x(:), beta, y(:, :)
    end subroutine psi_esctmv

    subroutine psi_esctv(n, idx, x, beta, y)
      import
      implicit none
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_epk_)    :: x(:), beta, y(:)
    end subroutine psi_esctv
  end interface psi_sct

  interface psi_exscan
    subroutine psi_e_exscanv(n, x, info, shift)
      import
      implicit none
      integer(psb_ipk_), intent(in) :: n
      integer(psb_epk_), intent(inout) :: x(:)
      integer(psb_ipk_), intent(out) :: info
      integer(psb_epk_), intent(in), optional :: shift
    end subroutine psi_e_exscanv
  end interface psi_exscan

  !Dispach of axpy-like operations.
  interface get_axpbylike_code
    module procedure get_axpbylike_code1
    module procedure get_axpbylike_code2
    module procedure get_axpbylike_code3
  end interface get_axpbylike_code

contains
  function get_axpbylike_code1(var) result(code)
    integer(psb_epk_), intent(in) :: var
    integer(psb_ipk_) :: code

    code = 0_psb_ipk_
    if(var ==  eone) code = 1_psb_ipk_
    if(var == ezero) code = 2_psb_ipk_
    if(var == -eone) code = 3_psb_ipk_
  end function get_axpbylike_code1

  function get_axpbylike_code2(alpha, beta) result(code)
    integer(psb_epk_), intent(in) :: alpha, beta
    integer(psb_ipk_) :: code

    code = get_axpbylike_code1(alpha) &
            + ishft(get_axpbylike_code1(beta), 2)
  end function get_axpbylike_code2

  function get_axpbylike_code3(alpha, beta, gamma) result(code)
    integer(psb_epk_), intent(in) :: alpha, beta, gamma
    integer(psb_ipk_) :: code

    code = get_axpbylike_code1(alpha) &
            + ishft(get_axpbylike_code1(beta), 2) &
            + ishft(get_axpbylike_code1(gamma), 4)
  end function get_axpbylike_code3
end module psi_e_serial_mod