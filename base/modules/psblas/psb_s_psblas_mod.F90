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
module psb_s_psblas_mod
  use psb_desc_mod, only : psb_desc_type, psb_spk_, psb_ipk_, psb_lpk_
  use psb_s_vect_mod, only : psb_s_vect_type
  use psb_s_multivect_mod, only : psb_s_multivect_type
  use psb_s_mat_mod, only : psb_sspmat_type

  interface psb_gedot
    function psb_sdot_vect(x, y, desc_a,info, global) result(res)
      import :: psb_s_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_s_vect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_sdot_vect

    function psb_sdotv(x, y, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: x(:), y(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_sdotv
    end function psb_sdotv

    function psb_sdot(x, y, desc_a, info, jx, jy, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: x(:, :), y(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx, jy
      logical, intent(in), optional           :: global
      real(psb_spk_)  :: psb_sdot
    end function psb_sdot
  end interface

  interface psb_gedots
    subroutine psb_sdotvs(res, x, y, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res
      real(psb_spk_), intent(in)      :: x(:), y(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional   :: global
    end subroutine psb_sdotvs

    subroutine psb_smdots(res, x, y, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res(:)
      real(psb_spk_), intent(in)      :: x(:, :), y(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional   :: global
    end subroutine psb_smdots

    ! mvect dot products now available only as subroutines. 
    ! Maybe worth to implemented them also as allocatable-output functions
    subroutine psb_sdot_mvect(x, y, xty, desc_a, info, global)
      import :: psb_s_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_s_multivect_type), intent(inout) :: x, y
      real(psb_spk_), intent(out)               :: xty(:, :)
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional             :: global
    end subroutine psb_sdot_mvect

    subroutine psb_sdot_mvect_vect(x, y, xty, desc_a, info, global)
      import :: psb_s_multivect_type, psb_s_vect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      type(psb_s_multivect_type), intent(inout) :: x
      type(psb_s_vect_type), intent(inout)      :: y
      real(psb_spk_), intent(out)               :: xty(:)
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional             :: global
    end subroutine psb_sdot_mvect_vect
  end interface

  interface psb_geaxpby
    subroutine psb_saxpby_vect(alpha, x, beta, y, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_saxpby_vect

    subroutine psb_saxpby_extract_c(alpha, x, idx_x, beta, y, desc_a, info)
      import :: psb_s_multivect_type, psb_s_vect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_s_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: idx_x
      type(psb_s_vect_type), intent(inout)      :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_extract_c

    subroutine psb_saxpby_mv_v_full(alpha, x, beta, y, desc_a, info)
      import :: psb_s_multivect_type, psb_s_vect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_s_vect_type), intent(inout)      :: x
      type(psb_s_multivect_type), intent(inout) :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_v_full

    subroutine psb_saxpby_mv_v_idxs(alpha, x, beta, y, idx_y, desc_a, info)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_s_vect_type), intent(inout)      :: x
      type(psb_s_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(in)             :: idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_v_idxs

    subroutine psb_saxpby_mv_m_full(alpha, x, beta, y, desc_a, info)
      import :: psb_s_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_s_multivect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_m_full

    subroutine psb_saxpby_mv_m_full_out(alpha, x, beta, y, z, desc_a, info)
      import :: psb_s_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      type(psb_s_multivect_type), intent(inout) :: x, y, z
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_m_full_out

    subroutine psb_saxpby_mv_m_idxs(alpha, x, idx_x, beta, y, idx_y, desc_a, info)
      import :: psb_s_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_s_multivect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_m_idxs

    subroutine psb_saxpby_mv_vv(alpha, x, beta, y, gamma, z, idx_z, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
              & psb_s_vect_type, psb_s_multivect_type
      real(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_s_vect_type), intent(inout)      :: x, y
      type(psb_s_multivect_type), intent(inout) :: z
      integer(psb_ipk_), intent(in)             :: idx_z
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_vv

    subroutine psb_saxpby_mv_mv(alpha, x, beta, y, idx_y, gamma, z, idx_z, desc_a, info)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_s_vect_type), intent(inout)      :: x
      type(psb_s_multivect_type), intent(inout) :: y, z
      integer(psb_ipk_), intent(in)             :: idx_y, idx_z
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_mv

    subroutine psb_saxpby_mv_mm_idxs(alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
              & psb_s_multivect_type
      real(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_s_multivect_type), intent(inout) :: x, y, z
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y, idx_z
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_mm_idxs

    subroutine psb_saxpby_mv_mm_full(alpha, x, beta, y, gamma, z, desc_a, info)
      import :: psb_s_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_s_multivect_type), intent(inout) :: x, y, z
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_mm_full

    subroutine psb_saxpby_mv_mm_out(alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, w, idx_w, desc_a, info)
      import :: psb_s_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_s_multivect_type), intent(inout) :: x, y, z, w
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y, idx_z, idx_w
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_saxpby_mv_mm_out

    subroutine psb_saxpby_mv_cspan1D(x, coeff, y, desc_a, info, upd_flag)
      import :: psb_s_multivect_type, psb_s_vect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      type(psb_s_multivect_type), intent(inout) :: x
      real(psb_spk_), intent(in)                :: coeff(:)
      type(psb_s_vect_type), intent(inout)      :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional             :: upd_flag
    end subroutine psb_saxpby_mv_cspan1D

    subroutine psb_saxpby_mv_cspan2D(x, coeff, y, desc_a, info, upd_flag)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
              & psb_s_multivect_type
      type(psb_s_multivect_type), intent(inout) :: x, y
      real(psb_spk_), intent(in)                :: coeff(:, :)
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional             :: upd_flag
    end subroutine psb_saxpby_mv_cspan2D

    subroutine psb_saxpby_vect_out(alpha, x, beta, y, z, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_s_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_saxpby_vect_out
    
    subroutine psb_saxpbyv(alpha, x, beta, y, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:)
      real(psb_spk_), intent(inout)   :: y(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_saxpbyv
    
    subroutine psb_saxpbyvout(alpha, x, beta, y, z, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:), y(:)
      real(psb_spk_), intent(inout)   :: z(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_saxpbyvout
    
    subroutine psb_saxpby(alpha, x, beta, y, desc_a, info, n, jx, jy)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: alpha, beta
      real(psb_spk_), intent(in)      :: x(:, :)
      real(psb_spk_), intent(inout)   :: y(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: n, jx, jy
    end subroutine psb_saxpby
  end interface

  interface psb_upd_xyz
    subroutine psb_s_upd_xyz_vect(alpha, beta, gamma, delta, x, y, z, &
                                & desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)            :: alpha, beta, gamma, delta
      type(psb_s_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_s_upd_xyz_vect
  end interface psb_upd_xyz
  
  interface psb_geamax
    function psb_samax(x, desc_a, info, jx, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx
      logical, intent(in), optional           :: global
      real(psb_spk_)  :: psb_samax
    end function psb_samax

    function psb_samaxv(x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_samaxv
    end function psb_samaxv

    function psb_samax_vect(x, desc_a, info, global) result(res)
      import :: psb_desc_type, psb_s_vect_type, psb_spk_, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_samax_vect
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_genrmi
    procedure psb_samax, psb_samaxv, psb_samax_vect
  end interface

  interface psb_normi
    procedure psb_samax, psb_samaxv, psb_samax_vect
  end interface
#endif

  interface psb_geamaxs
    subroutine psb_samaxvs(res, x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res
      real(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
    end subroutine psb_samaxvs

    subroutine psb_smamaxs(res, x, desc_a, info, jx, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res(:)
      real(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx
      logical, intent(in), optional           :: global
    end subroutine psb_smamaxs
  end interface

  interface psb_gemin
    function psb_smin_vect(x, desc_a, info, global) result(res)
      import :: psb_s_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_smin_vect
  end interface

  interface psb_geasum
    function psb_sasum(x, desc_a, info, jx, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx
      logical, intent(in), optional           :: global
      real(psb_spk_)  :: psb_sasum
    end function psb_sasum

    function psb_sasumv(x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_sasumv
    end function psb_sasumv
  
    function psb_sasum_vect(x, desc_a, info, global) result(res)
      import :: psb_s_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_sasum_vect
  end interface

  interface psb_geasums
    subroutine psb_sasumvs(res, x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res
      real(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
    end subroutine psb_sasumvs

    subroutine psb_smasum(res, x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res(:)
      real(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
    end subroutine psb_smasum
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_genrm1
    procedure psb_sasum, psb_sasumv, psb_sasum_vect
  end interface
  
  interface psb_norm1
    procedure psb_sasum, psb_sasumv, psb_sasum_vect
  end interface
#endif

  interface psb_genrm2
    function psb_snrm2(x, desc_a, info, jx, global) result(res)
      import :: psb_spk_, psb_ipk_, psb_desc_type
      real(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx
      logical, intent(in), optional           :: global
      real(psb_spk_)  :: res
    end function psb_snrm2

    function psb_snrm2v(x, desc_a, info, global) result(res)
      import :: psb_spk_, psb_ipk_, psb_desc_type
      real(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_snrm2v

    function psb_snrm2_vect(x, desc_a, info, global) result(res)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_s_vect_type), intent(inout)  :: x
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_snrm2_vect

    function psb_snrm2_mvect_full(x, desc_a, info, global) result(res)
      import :: psb_s_multivect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_s_multivect_type), intent(inout) :: x
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional :: global
      real(psb_spk_), allocatable :: res(:)
    end function psb_snrm2_mvect_full

    function psb_snrm2_mvect_idxs(x, idx_x, desc_a, info, global) result(res)
      import :: psb_s_multivect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_s_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: idx_x
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_snrm2_mvect_idxs

    function psb_snrm2_weight_vect(x, w, desc_a, info, global, aux) result(res)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_s_vect_type), intent(inout)  :: x
      type(psb_s_vect_type), intent(inout)  :: w
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional                   :: global
      type(psb_s_vect_type), intent(inout), optional  :: aux
      real(psb_spk_)  :: res
    end function psb_snrm2_weight_vect

    function psb_snrm2_weightmask_vect(x, w, idv, desc_a, info, global, aux) result(res)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_s_vect_type), intent(inout)  :: x
      type(psb_s_vect_type), intent(inout)  :: w
      type(psb_s_vect_type), intent(inout)  :: idv
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional                   :: global
      type(psb_s_vect_type), intent(inout), optional  :: aux
      real(psb_spk_)  :: res
    end function psb_snrm2_weightmask_vect
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_norm2
    procedure psb_snrm2, psb_snrm2v, psb_snrm2_vect, psb_snrm2_mvect_full, psb_snrm2_mvect_idxs, & 
            & psb_snrm2_weight_vect, psb_snrm2_weightmask_vect
  end interface
#endif

  interface psb_genrm2s
    subroutine psb_snrm2vs(res, x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res
      real(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
    end subroutine psb_snrm2vs
  end interface

  interface psb_spnrmi
    function psb_snrmi(a, desc_a, info, global)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_snrmi
    end function psb_snrmi
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_normi
    procedure psb_snrmi
  end interface
#endif

  interface psb_spnrm1
    function psb_sspnrm1(a, desc_a, info, global)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_sspnrm1
    end function psb_sspnrm1
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_norm1
    procedure psb_sspnrm1
  end interface
#endif

  interface psb_spmm
    subroutine psb_sspmm(alpha, a, x, beta, y, desc_a, info, &
                        & trans, k, jx, jy, work, doswap)
      import :: psb_desc_type, psb_sspmat_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_sspmat_type), intent(in)     :: a
      real(psb_spk_), intent(inout), target :: x(:, :), y(:, :)
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans
      real(psb_spk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), optional, intent(in)         :: k, jx, jy
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_sspmm
    
    subroutine psb_sspmv(alpha, a, x, beta, y, desc_a, info, &
                        & trans, work, doswap)
      import :: psb_desc_type, psb_sspmat_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_sspmat_type), intent(in)     :: a
      real(psb_spk_), intent(inout), target :: x(:), y(:)
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans
      real(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_sspmv
    
    subroutine psb_sspmv_vect(alpha, a, x, beta, y, desc_a, info, &
                              & trans, work, doswap)
      import :: psb_desc_type, psb_sspmat_type, psb_s_vect_type, &
              & psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_sspmat_type), intent(in)     :: a
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans
      real(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_sspmv_vect
    
    subroutine psb_sspmv_mv(alpha, a, x, beta, y, idx_y, desc_a, info, &
                            & trans, work, doswap)
      import :: psb_desc_type, psb_sspmat_type, &
              & psb_s_vect_type, psb_s_multivect_type, &
              & psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_sspmat_type), intent(in)         :: a
      type(psb_s_vect_type), intent(inout)      :: x
      type(psb_s_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(in)             :: idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans
      real(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_sspmv_mv

    subroutine psb_sspmv_vm(alpha, a, x, idx_x, beta, y, desc_a, info, &
                            & trans, work, doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_sspmat_type, &
              & psb_s_vect_type, psb_s_multivect_type
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_sspmat_type), intent(in)         :: a
      type(psb_s_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: idx_x
      type(psb_s_vect_type), intent(inout)      :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans
      real(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_sspmv_vm

    subroutine psb_sspmv_mm_idxs(alpha, a, x, idx_x, beta, y, idx_y, desc_a, info, &
                                & trans, work, doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_sspmat_type, psb_s_multivect_type
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_sspmat_type), intent(in)         :: a
      type(psb_s_multivect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans
      real(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_sspmv_mm_idxs

    subroutine psb_sspmv_mm_full(alpha, a, x, beta, y, desc_a, info, &
                                & trans, work, doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_s_multivect_type, psb_sspmat_type
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_sspmat_type), intent(in)         :: a
      type(psb_s_multivect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans
      real(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_sspmv_mm_full
  end interface

  interface psb_spsm
    subroutine psb_sspsm(alpha, t, x, beta, y, desc_a, info, &
                        & trans, scale, choice, diag, n, jx, jy, work)
      import :: psb_desc_type, psb_sspmat_type, psb_spk_, psb_ipk_ 
      implicit none
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_sspmat_type), intent(in)     :: t
      real(psb_spk_), intent(in), target    :: x(:, :)
      real(psb_spk_), intent(inout), target :: y(:, :)
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      real(psb_spk_), optional, intent(in), target    :: diag(:)
      integer(psb_ipk_), optional, intent(in)         :: n, jx, jy
      real(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_sspsm

    subroutine psb_sspsv(alpha, t, x, beta, y, desc_a, info, &
                        & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_sspmat_type, psb_spk_, psb_ipk_
      implicit none
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_sspmat_type), intent(in)     :: t
      real(psb_spk_), intent(in), target    :: x(:)
      real(psb_spk_), intent(inout), target :: y(:)
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      real(psb_spk_), optional, intent(in), target    :: diag(:)
      real(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_sspsv

    subroutine psb_sspsv_vect(alpha, t, x, beta, y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_s_vect_type, psb_sspmat_type, psb_spk_, psb_ipk_
      implicit none
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_sspmat_type), intent(inout)  :: t
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_s_vect_type), intent(inout), optional  :: diag
      real(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_sspsv_vect

    subroutine psb_sspsv_mv(alpha, t, x, idx_x, beta, y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_sspmat_type, &
              & psb_s_multivect_type, psb_s_vect_type
      implicit none
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_sspmat_type), intent(inout)      :: t
      type(psb_s_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: idx_x
      type(psb_s_vect_type), intent(inout)      :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_s_vect_type), intent(inout), optional  :: diag
      real(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_sspsv_mv

    subroutine psb_sspsv_vm(alpha, t, x, beta, y, idx_y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_sspmat_type, &
              & psb_s_multivect_type, psb_s_vect_type
      implicit none
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_sspmat_type), intent(inout)      :: t
      type(psb_s_vect_type), intent(inout)      :: x
      type(psb_s_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(in)             :: idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_s_vect_type), intent(inout), optional  :: diag
      real(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_sspsv_vm

    subroutine psb_sspsv_mm_full(alpha, t, x, beta, y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_sspmat_type, &
              & psb_s_multivect_type, psb_s_vect_type
      implicit none
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_sspmat_type), intent(inout)      :: t
      type(psb_s_multivect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_s_vect_type), intent(inout), optional  :: diag
      real(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_sspsv_mm_full

    subroutine psb_sspsv_mm_idxs(alpha, t, x, idx_x, beta, y, idx_y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_sspmat_type, &
              & psb_s_multivect_type, psb_s_vect_type
      implicit none
      real(psb_spk_), intent(in)                :: alpha, beta
      type(psb_sspmat_type), intent(inout)      :: t
      type(psb_s_multivect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_s_vect_type), intent(inout), optional  :: diag
      real(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_sspsv_mm_idxs
  end interface

  interface psb_gemlt
    subroutine psb_smlt_vect(x, y, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_smlt_vect

    subroutine psb_smlt_vect2(alpha, x, y, beta, z, desc_a, info, &
                            & conjgx, conjgy)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_s_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_vect2
    
    subroutine psb_smlt_mvect_v_full(alpha, x, y, beta, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_vect_type), intent(inout)       :: x
      class(psb_s_multivect_type), intent(inout)  :: y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_v_full

    subroutine psb_smlt_mvect_v_idxs(alpha, x, y, idx_y, beta, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_vect_type), intent(inout)       :: x
      class(psb_s_multivect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(in)               :: idx_y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_v_idxs
    
    subroutine psb_smlt_mvect_m_full(alpha, x, y, beta, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_s_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_multivect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_m_full

    subroutine psb_smlt_mvect_m_idxs(alpha, x, idx_x, y, idx_y, beta, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_s_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_multivect_type), intent(inout)  :: x, y
      integer(psb_ipk_), intent(in)               :: idx_x, idx_y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_m_idxs
    
    subroutine psb_smlt_mvect_vv_full_out(alpha, x, y, beta, z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_vect_type), intent(inout)       :: x, y
      class(psb_s_multivect_type), intent(inout)  :: z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_vv_full_out

    subroutine psb_smlt_mvect_vv_idxs_out(alpha, x, y, beta, z, idx_z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_vect_type), intent(inout)       :: x, y
      class(psb_s_multivect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(in)               :: idx_z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_vv_idxs_out

    subroutine psb_smlt_mvect_vm_full_out(alpha, x, y, beta, z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_vect_type), intent(inout)       :: x
      class(psb_s_multivect_type), intent(inout)  :: y, z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_vm_full_out

    subroutine psb_smlt_mvect_vm_idxs_out(alpha, x, y, idx_y, beta, z, idx_z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_vect_type), intent(inout)       :: x
      class(psb_s_multivect_type), intent(inout)  :: y, z
      integer(psb_ipk_), intent(in)               :: idx_y, idx_z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_vm_idxs_out

    subroutine psb_smlt_mvect_mm_full_out(alpha, x, y, beta, z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_s_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_multivect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_mm_full_out

    subroutine psb_smlt_mvect_mm_idxs_out(alpha, x, idx_x, y, idx_y, beta, z, idx_z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_s_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_multivect_type), intent(inout)  :: x, y, z
      integer(psb_ipk_), intent(in)               :: idx_x, idx_y, idx_z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_mm_idxs_out

    subroutine psb_smlt_mvect_vm_ext(alpha, x, y, idx_y, beta, z, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_vect_type), intent(inout)       :: x, z
      class(psb_s_multivect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(in)               :: idx_y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_vm_ext

    subroutine psb_smlt_mvect_mm_ext(alpha, x, idx_x, y, idx_y, beta, z, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_s_vect_type, psb_s_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_s_multivect_type), intent(inout)  :: x, y
      integer(psb_ipk_), intent(in)               :: idx_x, idx_y
      class(psb_s_vect_type), intent(inout)       :: z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_smlt_mvect_mm_ext
  end interface

  interface psb_gediv
    subroutine psb_sdiv_vect(x, y, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_sdiv_vect

    subroutine psb_sdiv_vect2(x, y, z, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_sdiv_vect2

    subroutine psb_sdiv_vect_check(x, y, desc_a, info, flag)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in) :: flag
    end subroutine psb_sdiv_vect_check

    subroutine psb_sdiv_vect2_check(x, y, z, desc_a, info, flag)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in)                   :: flag
    end subroutine psb_sdiv_vect2_check
  end interface

  interface psb_geinv
    subroutine psb_sinv_vect(x, y, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_sinv_vect

    subroutine psb_sinv_vect_check(x, y, desc_a, info, flag)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in) :: flag
    end subroutine psb_sinv_vect_check
  end interface

  interface psb_geabs
    subroutine psb_sabs_vect(x, y, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_sabs_vect
  end interface

  interface psb_gecmp
    subroutine psb_scmp_vect(x, c, z, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_s_vect_type), intent(inout) :: x, z
      real(psb_spk_), intent(in)           :: c
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
    end subroutine psb_scmp_vect

    subroutine psb_scmp_spmatval(a, val, tol, desc_a, res, info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_, psb_lpk_
      type(psb_sspmat_type), intent(inout)  :: a
      real(psb_spk_), intent(in)            :: val
      real(psb_spk_), intent(in)            :: tol
      type(psb_desc_type), intent(in)       :: desc_a
      logical, intent(out)                  :: res
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_scmp_spmatval

    subroutine psb_scmp_spmat(a, b, tol, desc_a, res, info)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_lpk_, psb_ipk_
      type(psb_sspmat_type), intent(inout)  :: a, b
      real(psb_spk_), intent(in)            :: tol
      type(psb_desc_type), intent(in)       :: desc_a
      logical, intent(out)                  :: res
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_scmp_spmat
  end interface

  interface psb_geaddconst
    subroutine psb_saddconst_vect(x, b, z, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_s_vect_type), intent(inout) :: x, z
      real(psb_spk_), intent(in)           :: b
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
    end subroutine psb_saddconst_vect
  end interface

  interface psb_mask
    subroutine psb_smask_vect(c, x, m, t, desc_a, info)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_s_vect_type), intent(inout)  :: c, x, m
      logical, intent(out)                  :: t
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_smask_vect
  end interface
  
  interface psb_minquotient
    function psb_sminquotient_vect(x, y, desc_a, info, global) result(res)
      import :: psb_s_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_s_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function
  end interface

  interface psb_nnz
    function psb_sget_nnz(a, desc_a, info) result(res)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_, psb_lpk_
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      integer(psb_lpk_) :: res
    end function
  end interface

  interface psb_is_matupd
    function psb_s_is_matupd(a, desc_a, info) result(res)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical :: res
    end function
  end interface

  interface psb_is_matasb
    function psb_s_is_matasb(a, desc_a, info) result(res)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical :: res
    end function
  end interface

  interface psb_is_matbld
    function psb_s_is_matbld(a, desc_a, info) result(res)
      import :: psb_sspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_sspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical :: res
    end function
  end interface
end module psb_s_psblas_mod