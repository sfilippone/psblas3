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
module psb_c_psblas_mod
  use psb_desc_mod, only : psb_desc_type, psb_spk_, psb_ipk_, psb_lpk_
  use psb_c_vect_mod, only : psb_c_vect_type
  use psb_c_multivect_mod, only : psb_c_multivect_type
  use psb_c_mat_mod, only : psb_cspmat_type

  interface psb_gedot
    function psb_cdot_vect(x, y, desc_a,info, global) result(res)
      import :: psb_c_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_c_vect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional :: global
      complex(psb_spk_)  :: res
    end function psb_cdot_vect

    function psb_cdotv(x, y, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: x(:), y(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
      complex(psb_spk_)  :: psb_cdotv
    end function psb_cdotv

    function psb_cdot(x, y, desc_a, info, jx, jy, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: x(:, :), y(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx, jy
      logical, intent(in), optional           :: global
      complex(psb_spk_)  :: psb_cdot
    end function psb_cdot
  end interface

  interface psb_gedots
    subroutine psb_cdotvs(res, x, y, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(out)     :: res
      complex(psb_spk_), intent(in)      :: x(:), y(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional   :: global
    end subroutine psb_cdotvs

    subroutine psb_cmdots(res, x, y, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(out)     :: res(:)
      complex(psb_spk_), intent(in)      :: x(:, :), y(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional   :: global
    end subroutine psb_cmdots

    ! mvect dot products now available only as subroutines. 
    ! Maybe worth to implemented them also as allocatable-output functions
    subroutine psb_cdot_mvect(x, y, xty, desc_a, info, global)
      import :: psb_c_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_c_multivect_type), intent(inout) :: x, y
      complex(psb_spk_), intent(out)               :: xty(:, :)
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional             :: global
    end subroutine psb_cdot_mvect

    subroutine psb_cdot_mvect_vect(x, y, xty, desc_a, info, global)
      import :: psb_c_multivect_type, psb_c_vect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      type(psb_c_multivect_type), intent(inout) :: x
      type(psb_c_vect_type), intent(inout)      :: y
      complex(psb_spk_), intent(out)               :: xty(:)
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional             :: global
    end subroutine psb_cdot_mvect_vect
  end interface

  interface psb_geaxpby
    subroutine psb_caxpby_vect(alpha, x, beta, y, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_caxpby_vect

    subroutine psb_caxpby_extract_c(alpha, x, idx_x, beta, y, desc_a, info)
      import :: psb_c_multivect_type, psb_c_vect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_c_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: idx_x
      type(psb_c_vect_type), intent(inout)      :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_extract_c

    subroutine psb_caxpby_mv_v_full(alpha, x, beta, y, desc_a, info)
      import :: psb_c_multivect_type, psb_c_vect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_c_vect_type), intent(inout)      :: x
      type(psb_c_multivect_type), intent(inout) :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_v_full

    subroutine psb_caxpby_mv_v_idxs(alpha, x, beta, y, idx_y, desc_a, info)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_c_vect_type), intent(inout)      :: x
      type(psb_c_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(in)             :: idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_v_idxs

    subroutine psb_caxpby_mv_m_full(alpha, x, beta, y, desc_a, info)
      import :: psb_c_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_c_multivect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_m_full

    subroutine psb_caxpby_mv_m_full_out(alpha, x, beta, y, z, desc_a, info)
      import :: psb_c_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      type(psb_c_multivect_type), intent(inout) :: x, y, z
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_m_full_out

    subroutine psb_caxpby_mv_m_idxs(alpha, x, idx_x, beta, y, idx_y, desc_a, info)
      import :: psb_c_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_c_multivect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_m_idxs

    subroutine psb_caxpby_mv_vv(alpha, x, beta, y, gamma, z, idx_z, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
              & psb_c_vect_type, psb_c_multivect_type
      complex(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_c_vect_type), intent(inout)      :: x, y
      type(psb_c_multivect_type), intent(inout) :: z
      integer(psb_ipk_), intent(in)             :: idx_z
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_vv

    subroutine psb_caxpby_mv_mv(alpha, x, beta, y, idx_y, gamma, z, idx_z, desc_a, info)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_c_vect_type), intent(inout)      :: x
      type(psb_c_multivect_type), intent(inout) :: y, z
      integer(psb_ipk_), intent(in)             :: idx_y, idx_z
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_mv

    subroutine psb_caxpby_mv_mm_idxs(alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
              & psb_c_multivect_type
      complex(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_c_multivect_type), intent(inout) :: x, y, z
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y, idx_z
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_mm_idxs

    subroutine psb_caxpby_mv_mm_full(alpha, x, beta, y, gamma, z, desc_a, info)
      import :: psb_c_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_c_multivect_type), intent(inout) :: x, y, z
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_mm_full

    subroutine psb_caxpby_mv_mm_out(alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, w, idx_w, desc_a, info)
      import :: psb_c_multivect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta, gamma
      type(psb_c_multivect_type), intent(inout) :: x, y, z, w
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y, idx_z, idx_w
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psb_caxpby_mv_mm_out

    subroutine psb_caxpby_mv_cspan1D(x, coeff, y, desc_a, info, upd_flag)
      import :: psb_c_multivect_type, psb_c_vect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      type(psb_c_multivect_type), intent(inout) :: x
      complex(psb_spk_), intent(in)                :: coeff(:)
      type(psb_c_vect_type), intent(inout)      :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional             :: upd_flag
    end subroutine psb_caxpby_mv_cspan1D

    subroutine psb_caxpby_mv_cspan2D(x, coeff, y, desc_a, info, upd_flag)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
              & psb_c_multivect_type
      type(psb_c_multivect_type), intent(inout) :: x, y
      complex(psb_spk_), intent(in)                :: coeff(:, :)
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional             :: upd_flag
    end subroutine psb_caxpby_mv_cspan2D

    subroutine psb_caxpby_vect_out(alpha, x, beta, y, z, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, &
              & psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_c_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_caxpby_vect_out
    
    subroutine psb_caxpbyv(alpha, x, beta, y, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: alpha, beta
      complex(psb_spk_), intent(in)      :: x(:)
      complex(psb_spk_), intent(inout)   :: y(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_caxpbyv
    
    subroutine psb_caxpbyvout(alpha, x, beta, y, z, desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: alpha, beta
      complex(psb_spk_), intent(in)      :: x(:), y(:)
      complex(psb_spk_), intent(inout)   :: z(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
    end subroutine psb_caxpbyvout
    
    subroutine psb_caxpby(alpha, x, beta, y, desc_a, info, n, jx, jy)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: alpha, beta
      complex(psb_spk_), intent(in)      :: x(:, :)
      complex(psb_spk_), intent(inout)   :: y(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: n, jx, jy
    end subroutine psb_caxpby
  end interface

  interface psb_upd_xyz
    subroutine psb_c_upd_xyz_vect(alpha, beta, gamma, delta, x, y, z, &
                                & desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)            :: alpha, beta, gamma, delta
      type(psb_c_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_c_upd_xyz_vect
  end interface psb_upd_xyz
  
  interface psb_geamax
    function psb_camax(x, desc_a, info, jx, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx
      logical, intent(in), optional           :: global
      real(psb_spk_)  :: psb_camax
    end function psb_camax

    function psb_camaxv(x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_camaxv
    end function psb_camaxv

    function psb_camax_vect(x, desc_a, info, global) result(res)
      import :: psb_desc_type, psb_c_vect_type, psb_spk_, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_camax_vect
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_genrmi
    procedure psb_camax, psb_camaxv, psb_camax_vect
  end interface

  interface psb_normi
    procedure psb_camax, psb_camaxv, psb_camax_vect
  end interface
#endif

  interface psb_geamaxs
    subroutine psb_camaxvs(res, x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res
      complex(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
    end subroutine psb_camaxvs

    subroutine psb_cmamaxs(res, x, desc_a, info, jx, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res(:)
      complex(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx
      logical, intent(in), optional           :: global
    end subroutine psb_cmamaxs
  end interface


  interface psb_geasum
    function psb_casum(x, desc_a, info, jx, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx
      logical, intent(in), optional           :: global
      real(psb_spk_)  :: psb_casum
    end function psb_casum

    function psb_casumv(x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_casumv
    end function psb_casumv
  
    function psb_casum_vect(x, desc_a, info, global) result(res)
      import :: psb_c_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_casum_vect
  end interface

  interface psb_geasums
    subroutine psb_casumvs(res, x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res
      complex(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
    end subroutine psb_casumvs

    subroutine psb_cmasum(res, x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res(:)
      complex(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
    end subroutine psb_cmasum
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_genrm1
    procedure psb_casum, psb_casumv, psb_casum_vect
  end interface
  
  interface psb_norm1
    procedure psb_casum, psb_casumv, psb_casum_vect
  end interface
#endif

  interface psb_genrm2
    function psb_cnrm2(x, desc_a, info, jx, global) result(res)
      import :: psb_spk_, psb_ipk_, psb_desc_type
      complex(psb_spk_), intent(in)      :: x(:, :)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      integer(psb_ipk_), optional, intent(in) :: jx
      logical, intent(in), optional           :: global
      real(psb_spk_)  :: res
    end function psb_cnrm2

    function psb_cnrm2v(x, desc_a, info, global) result(res)
      import :: psb_spk_, psb_ipk_, psb_desc_type
      complex(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_cnrm2v

    function psb_cnrm2_vect(x, desc_a, info, global) result(res)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_c_vect_type), intent(inout)  :: x
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_cnrm2_vect

    function psb_cnrm2_mvect_full(x, desc_a, info, global) result(res)
      import :: psb_c_multivect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_c_multivect_type), intent(inout) :: x
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional :: global
      real(psb_spk_), allocatable :: res(:)
    end function psb_cnrm2_mvect_full

    function psb_cnrm2_mvect_idxs(x, idx_x, desc_a, info, global) result(res)
      import :: psb_c_multivect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_c_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: idx_x
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: res
    end function psb_cnrm2_mvect_idxs

    function psb_cnrm2_weight_vect(x, w, desc_a, info, global, aux) result(res)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_c_vect_type), intent(inout)  :: x
      type(psb_c_vect_type), intent(inout)  :: w
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional                   :: global
      type(psb_c_vect_type), intent(inout), optional  :: aux
      real(psb_spk_)  :: res
    end function psb_cnrm2_weight_vect

    function psb_cnrm2_weightmask_vect(x, w, idv, desc_a, info, global, aux) result(res)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_c_vect_type), intent(inout)  :: x
      type(psb_c_vect_type), intent(inout)  :: w
      type(psb_c_vect_type), intent(inout)  :: idv
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional                   :: global
      type(psb_c_vect_type), intent(inout), optional  :: aux
      real(psb_spk_)  :: res
    end function psb_cnrm2_weightmask_vect
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_norm2
    procedure psb_cnrm2, psb_cnrm2v, psb_cnrm2_vect, psb_cnrm2_mvect_full, psb_cnrm2_mvect_idxs, & 
            & psb_cnrm2_weight_vect, psb_cnrm2_weightmask_vect
  end interface
#endif

  interface psb_genrm2s
    subroutine psb_cnrm2vs(res, x, desc_a, info, global)
      import :: psb_desc_type, psb_spk_, psb_ipk_
      real(psb_spk_), intent(out)     :: res
      complex(psb_spk_), intent(in)      :: x(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)  :: info
      logical, intent(in), optional :: global
    end subroutine psb_cnrm2vs
  end interface

  interface psb_spnrmi
    function psb_cnrmi(a, desc_a, info, global)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_cspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_cnrmi
    end function psb_cnrmi
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_normi
    procedure psb_cnrmi
  end interface
#endif

  interface psb_spnrm1
    function psb_cspnrm1(a, desc_a, info, global)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_cspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical, intent(in), optional :: global
      real(psb_spk_)  :: psb_cspnrm1
    end function psb_cspnrm1
  end interface

#if !defined(PSB_HAVE_BUGGY_GENERICS)
  interface psb_norm1
    procedure psb_cspnrm1
  end interface
#endif

  interface psb_spmm
    subroutine psb_cspmm(alpha, a, x, beta, y, desc_a, info, &
                        & trans, k, jx, jy, work, doswap)
      import :: psb_desc_type, psb_cspmat_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_cspmat_type), intent(in)     :: a
      complex(psb_spk_), intent(inout), target :: x(:, :), y(:, :)
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), optional, intent(in)         :: k, jx, jy
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_cspmm
    
    subroutine psb_cspmv(alpha, a, x, beta, y, desc_a, info, &
                        & trans, work, doswap)
      import :: psb_desc_type, psb_cspmat_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_cspmat_type), intent(in)     :: a
      complex(psb_spk_), intent(inout), target :: x(:), y(:)
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_cspmv
    
    subroutine psb_cspmv_vect(alpha, a, x, beta, y, desc_a, info, &
                              & trans, work, doswap)
      import :: psb_desc_type, psb_cspmat_type, psb_c_vect_type, &
              & psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_cspmat_type), intent(in)     :: a
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_cspmv_vect
    
    subroutine psb_cspmv_mv(alpha, a, x, beta, y, idx_y, desc_a, info, &
                            & trans, work, doswap)
      import :: psb_desc_type, psb_cspmat_type, &
              & psb_c_vect_type, psb_c_multivect_type, &
              & psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_cspmat_type), intent(in)         :: a
      type(psb_c_vect_type), intent(inout)      :: x
      type(psb_c_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(in)             :: idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_cspmv_mv

    subroutine psb_cspmv_vm(alpha, a, x, idx_x, beta, y, desc_a, info, &
                            & trans, work, doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_cspmat_type, &
              & psb_c_vect_type, psb_c_multivect_type
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_cspmat_type), intent(in)         :: a
      type(psb_c_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: idx_x
      type(psb_c_vect_type), intent(inout)      :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_cspmv_vm

    subroutine psb_cspmv_mm_idxs(alpha, a, x, idx_x, beta, y, idx_y, desc_a, info, &
                                & trans, work, doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_cspmat_type, psb_c_multivect_type
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_cspmat_type), intent(in)         :: a
      type(psb_c_multivect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_cspmv_mm_idxs

    subroutine psb_cspmv_mm_full(alpha, a, x, beta, y, desc_a, info, &
                                & trans, work, doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_c_multivect_type, psb_cspmat_type
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_cspmat_type), intent(in)         :: a
      type(psb_c_multivect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      logical, optional, intent(in)                   :: doswap
    end subroutine psb_cspmv_mm_full
  end interface

  interface psb_spsm
    subroutine psb_cspsm(alpha, t, x, beta, y, desc_a, info, &
                        & trans, scale, choice, diag, n, jx, jy, work)
      import :: psb_desc_type, psb_cspmat_type, psb_spk_, psb_ipk_ 
      implicit none
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_cspmat_type), intent(in)     :: t
      complex(psb_spk_), intent(in), target    :: x(:, :)
      complex(psb_spk_), intent(inout), target :: y(:, :)
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      complex(psb_spk_), optional, intent(in), target    :: diag(:)
      integer(psb_ipk_), optional, intent(in)         :: n, jx, jy
      complex(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_cspsm

    subroutine psb_cspsv(alpha, t, x, beta, y, desc_a, info, &
                        & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_cspmat_type, psb_spk_, psb_ipk_
      implicit none
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_cspmat_type), intent(in)     :: t
      complex(psb_spk_), intent(in), target    :: x(:)
      complex(psb_spk_), intent(inout), target :: y(:)
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      complex(psb_spk_), optional, intent(in), target    :: diag(:)
      complex(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_cspsv

    subroutine psb_cspsv_vect(alpha, t, x, beta, y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_c_vect_type, psb_cspmat_type, psb_spk_, psb_ipk_
      implicit none
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_cspmat_type), intent(inout)  :: t
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_c_vect_type), intent(inout), optional  :: diag
      complex(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_cspsv_vect

    subroutine psb_cspsv_mv(alpha, t, x, idx_x, beta, y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_cspmat_type, &
              & psb_c_multivect_type, psb_c_vect_type
      implicit none
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_cspmat_type), intent(inout)      :: t
      type(psb_c_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: idx_x
      type(psb_c_vect_type), intent(inout)      :: y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_c_vect_type), intent(inout), optional  :: diag
      complex(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_cspsv_mv

    subroutine psb_cspsv_vm(alpha, t, x, beta, y, idx_y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_cspmat_type, &
              & psb_c_multivect_type, psb_c_vect_type
      implicit none
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_cspmat_type), intent(inout)      :: t
      type(psb_c_vect_type), intent(inout)      :: x
      type(psb_c_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(in)             :: idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_c_vect_type), intent(inout), optional  :: diag
      complex(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_cspsv_vm

    subroutine psb_cspsv_mm_full(alpha, t, x, beta, y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_cspmat_type, &
              & psb_c_multivect_type, psb_c_vect_type
      implicit none
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_cspmat_type), intent(inout)      :: t
      type(psb_c_multivect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_c_vect_type), intent(inout), optional  :: diag
      complex(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_cspsv_mm_full

    subroutine psb_cspsv_mm_idxs(alpha, t, x, idx_x, beta, y, idx_y, desc_a, info, &
                            & trans, scale, choice, diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, psb_cspmat_type, &
              & psb_c_multivect_type, psb_c_vect_type
      implicit none
      complex(psb_spk_), intent(in)                :: alpha, beta
      type(psb_cspmat_type), intent(inout)      :: t
      type(psb_c_multivect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)             :: idx_x, idx_y
      type(psb_desc_type), intent(in)           :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      character, optional, intent(in)                 :: trans, scale
      integer(psb_ipk_), optional, intent(in)         :: choice
      type(psb_c_vect_type), intent(inout), optional  :: diag
      complex(psb_spk_), optional, intent(inout), target :: work(:)
    end subroutine psb_cspsv_mm_idxs
  end interface

  interface psb_gemlt
    subroutine psb_cmlt_vect(x, y, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_cmlt_vect

    subroutine psb_cmlt_vect2(alpha, x, y, beta, z, desc_a, info, &
                            & conjgx, conjgy)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_c_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_vect2
    
    subroutine psb_cmlt_mvect_v_full(alpha, x, y, beta, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_vect_type), intent(inout)       :: x
      class(psb_c_multivect_type), intent(inout)  :: y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_v_full

    subroutine psb_cmlt_mvect_v_idxs(alpha, x, y, idx_y, beta, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_vect_type), intent(inout)       :: x
      class(psb_c_multivect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(in)               :: idx_y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_v_idxs
    
    subroutine psb_cmlt_mvect_m_full(alpha, x, y, beta, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_c_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_multivect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_m_full

    subroutine psb_cmlt_mvect_m_idxs(alpha, x, idx_x, y, idx_y, beta, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_c_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_multivect_type), intent(inout)  :: x, y
      integer(psb_ipk_), intent(in)               :: idx_x, idx_y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_m_idxs
    
    subroutine psb_cmlt_mvect_vv_full_out(alpha, x, y, beta, z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_vect_type), intent(inout)       :: x, y
      class(psb_c_multivect_type), intent(inout)  :: z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_vv_full_out

    subroutine psb_cmlt_mvect_vv_idxs_out(alpha, x, y, beta, z, idx_z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_vect_type), intent(inout)       :: x, y
      class(psb_c_multivect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(in)               :: idx_z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_vv_idxs_out

    subroutine psb_cmlt_mvect_vm_full_out(alpha, x, y, beta, z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_vect_type), intent(inout)       :: x
      class(psb_c_multivect_type), intent(inout)  :: y, z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_vm_full_out

    subroutine psb_cmlt_mvect_vm_idxs_out(alpha, x, y, idx_y, beta, z, idx_z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_vect_type), intent(inout)       :: x
      class(psb_c_multivect_type), intent(inout)  :: y, z
      integer(psb_ipk_), intent(in)               :: idx_y, idx_z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_vm_idxs_out

    subroutine psb_cmlt_mvect_mm_full_out(alpha, x, y, beta, z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_c_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_multivect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_mm_full_out

    subroutine psb_cmlt_mvect_mm_idxs_out(alpha, x, idx_x, y, idx_y, beta, z, idx_z, desc_a, info, &
                                        & conjgx, conjgy)
      import :: psb_c_multivect_type, psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_multivect_type), intent(inout)  :: x, y, z
      integer(psb_ipk_), intent(in)               :: idx_x, idx_y, idx_z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_mm_idxs_out

    subroutine psb_cmlt_mvect_vm_ext(alpha, x, y, idx_y, beta, z, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_vect_type), intent(inout)       :: x, z
      class(psb_c_multivect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(in)               :: idx_y
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_vm_ext

    subroutine psb_cmlt_mvect_mm_ext(alpha, x, idx_x, y, idx_y, beta, z, desc_a, info, &
                                    & conjgx, conjgy)
      import :: psb_c_vect_type, psb_c_multivect_type, &
              & psb_desc_type, psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)                  :: alpha, beta
      class(psb_c_multivect_type), intent(inout)  :: x, y
      integer(psb_ipk_), intent(in)               :: idx_x, idx_y
      class(psb_c_vect_type), intent(inout)       :: z
      type(psb_desc_type), intent(in)             :: desc_a
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional  :: conjgx, conjgy
    end subroutine psb_cmlt_mvect_mm_ext
  end interface

  interface psb_gediv
    subroutine psb_cdiv_vect(x, y, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_cdiv_vect

    subroutine psb_cdiv_vect2(x, y, z, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_cdiv_vect2

    subroutine psb_cdiv_vect_check(x, y, desc_a, info, flag)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in) :: flag
    end subroutine psb_cdiv_vect_check

    subroutine psb_cdiv_vect2_check(x, y, z, desc_a, info, flag)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x, y, z
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in)                   :: flag
    end subroutine psb_cdiv_vect2_check
  end interface

  interface psb_geinv
    subroutine psb_cinv_vect(x, y, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_cinv_vect

    subroutine psb_cinv_vect_check(x, y, desc_a, info, flag)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in) :: flag
    end subroutine psb_cinv_vect_check
  end interface

  interface psb_geabs
    subroutine psb_cabs_vect(x, y, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_
      type(psb_c_vect_type), intent(inout)  :: x, y
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_cabs_vect
  end interface

  interface psb_gecmp
    subroutine psb_ccmp_vect(x, c, z, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_c_vect_type), intent(inout) :: x, z
      real(psb_spk_), intent(in)           :: c
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
    end subroutine psb_ccmp_vect

    subroutine psb_ccmp_spmatval(a, val, tol, desc_a, res, info)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_, psb_lpk_
      type(psb_cspmat_type), intent(inout)  :: a
      complex(psb_spk_), intent(in)            :: val
      real(psb_spk_), intent(in)            :: tol
      type(psb_desc_type), intent(in)       :: desc_a
      logical, intent(out)                  :: res
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_ccmp_spmatval

    subroutine psb_ccmp_spmat(a, b, tol, desc_a, res, info)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_lpk_, psb_ipk_
      type(psb_cspmat_type), intent(inout)  :: a, b
      real(psb_spk_), intent(in)            :: tol
      type(psb_desc_type), intent(in)       :: desc_a
      logical, intent(out)                  :: res
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_ccmp_spmat
  end interface

  interface psb_geaddconst
    subroutine psb_caddconst_vect(x, b, z, desc_a, info)
      import :: psb_c_vect_type, psb_desc_type, psb_ipk_, psb_spk_
      type(psb_c_vect_type), intent(inout) :: x, z
      real(psb_spk_), intent(in)           :: b
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
    end subroutine psb_caddconst_vect
  end interface


  interface psb_nnz
    function psb_cget_nnz(a, desc_a, info) result(res)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_, psb_lpk_
      type(psb_cspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      integer(psb_lpk_) :: res
    end function
  end interface

  interface psb_is_matupd
    function psb_c_is_matupd(a, desc_a, info) result(res)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_cspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical :: res
    end function
  end interface

  interface psb_is_matasb
    function psb_c_is_matasb(a, desc_a, info) result(res)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_cspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical :: res
    end function
  end interface

  interface psb_is_matbld
    function psb_c_is_matbld(a, desc_a, info) result(res)
      import :: psb_cspmat_type, psb_desc_type, psb_spk_, psb_ipk_
      type(psb_cspmat_type), intent(in) :: a
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical :: res
    end function
  end interface
end module psb_c_psblas_mod