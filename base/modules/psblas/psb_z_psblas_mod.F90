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
module psb_z_psblas_mod
  use psb_desc_mod, only : psb_desc_type, psb_dpk_, psb_ipk_
  use psb_z_vect_mod, only : psb_z_vect_type
  use psb_z_mat_mod, only : psb_zspmat_type

  interface psb_gedot
    function psb_zdot_vect(x, y, desc_a,info,global) result(res)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      complex(psb_dpk_)                    :: res
      type(psb_z_vect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_zdot_vect
    function psb_zdotv(x, y, desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      complex(psb_dpk_)                :: psb_zdotv
      complex(psb_dpk_), intent(in)    :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_zdotv
    function psb_zdot(x, y, desc_a, info, jx, jy,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      complex(psb_dpk_)                :: psb_zdot
      complex(psb_dpk_), intent(in)    :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), optional, intent(in)      :: jx, jy
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_zdot
  end interface


  interface psb_gedots
    subroutine  psb_zdotvs(res,x, y, desc_a, info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      complex(psb_dpk_), intent(out)      :: res
      complex(psb_dpk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end subroutine psb_zdotvs
    subroutine  psb_zmdots(res,x, y, desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      complex(psb_dpk_), intent(out)      :: res(:)
      complex(psb_dpk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end subroutine psb_zmdots
  end interface

  interface psb_geaxpby
    subroutine psb_zaxpby_vect(alpha, x, beta, y,&
         & desc_a, info)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_z_vect_type), intent (inout) :: y
      complex(psb_dpk_), intent (in)        :: alpha, beta
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_zaxpby_vect
    subroutine psb_zaxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      complex(psb_dpk_), intent (in)       ::  x(:)
      complex(psb_dpk_), intent (inout)    ::  y(:)
      complex(psb_dpk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_zaxpbyv
    subroutine psb_zaxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      complex(psb_dpk_), intent (in)       ::  x(:,:)
      complex(psb_dpk_), intent (inout)    ::  y(:,:)
      complex(psb_dpk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent(in) :: n, jx, jy
      integer(psb_ipk_), intent(out)      :: info
    end subroutine psb_zaxpby
  end interface

  interface psb_geamax
    function psb_zamax(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)   psb_zamax
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_zamax
    function psb_zamaxv(x, desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_) psb_zamaxv
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_zamaxv
    function psb_zamax_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)                        :: res
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_zamax_vect
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_genrmi
    procedure psb_zamax, psb_zamaxv, psb_zamax_vect
  end interface
  interface psb_normi
    procedure psb_zamax, psb_zamaxv, psb_zamax_vect
  end interface
#endif

  interface psb_geamaxs
    subroutine  psb_zamaxvs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)    :: x(:)
      type(psb_desc_type), intent (in)  :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical, intent(in), optional     :: global
    end subroutine psb_zamaxvs
    subroutine  psb_zmamaxs(res,x,desc_a,info,jx,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_), intent (out)       :: res(:)
      complex(psb_dpk_), intent (in)     :: x(:,:)
      type(psb_desc_type), intent (in)   :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      integer(psb_ipk_), optional, intent(in)      :: jx
      logical, intent(in), optional      :: global
    end subroutine psb_zmamaxs
  end interface

  interface psb_geasum
    function psb_zasum_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)                        :: res
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_zasum_vect
    function psb_zasum(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)   psb_zasum
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_zasum
    function psb_zasumv(x, desc_a, info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_) psb_zasumv
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_zasumv
  end interface

  interface psb_geasums
    subroutine  psb_zasumvs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_zasumvs
    subroutine  psb_zmasum(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_), intent (out)      :: res(:)
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_zmasum
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_genrm1
    procedure psb_zasum, psb_zasumv, psb_zasum_vect
  end interface
  interface psb_norm1
    procedure psb_zasum, psb_zasumv, psb_zasum_vect
  end interface
#endif

  interface psb_genrm2
    function psb_znrm2(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)   psb_znrm2
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_znrm2
    function psb_znrm2v(x, desc_a, info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_) psb_znrm2v
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)                :: info
      logical, intent(in), optional        :: global
    end function psb_znrm2v
    function psb_znrm2_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)                      :: res
      type(psb_z_vect_type), intent (inout)   :: x
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_znrm2_vect
    function psb_znrm2_weight_vect(x,w, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)                      :: res
      type(psb_z_vect_type), intent (inout)   :: x
      type(psb_z_vect_type), intent (inout)   :: w
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_znrm2_weight_vect
    function psb_znrm2_weightmask_vect(x,w,idv, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)                      :: res
      type(psb_z_vect_type), intent (inout)   :: x
      type(psb_z_vect_type), intent (inout)   :: w
      type(psb_z_vect_type), intent (inout)   :: idv
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_znrm2_weightmask_vect
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_norm2
    procedure psb_znrm2, psb_znrm2v, psb_znrm2_vect, psb_znrm2_weight_vect, psb_znrm2_weightmask_vect
  end interface
#endif

  interface psb_genrm2s
    subroutine  psb_znrm2vs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_znrm2vs
  end interface


  interface psb_spnrmi
    function psb_znrmi(a, desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)                    :: psb_znrmi
      type(psb_zspmat_type), intent (in) :: a
      type(psb_desc_type), intent (in)   :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_znrmi
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_normi
    procedure psb_znrmi
  end interface
#endif

  interface psb_spnrm1
    function psb_zspnrm1(a, desc_a,info,global)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      real(psb_dpk_)                     :: psb_zspnrm1
      type(psb_zspmat_type), intent (in) :: a
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_zspnrm1
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_norm1
    procedure psb_zspnrm1
  end interface
#endif

  interface psb_spmm
    subroutine psb_zspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      type(psb_zspmat_type), intent(in)        :: a
      complex(psb_dpk_), intent(inout), target :: x(:,:)
      complex(psb_dpk_), intent(inout), target :: y(:,:)
      complex(psb_dpk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans
      complex(psb_dpk_), optional, intent(inout),target :: work(:)
      integer(psb_ipk_), optional, intent(in)        :: k, jx, jy
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_zspmm
    subroutine psb_zspmv(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      type(psb_zspmat_type), intent(in)        :: a
      complex(psb_dpk_), intent(inout), target :: x(:)
      complex(psb_dpk_), intent(inout), target :: y(:)
      complex(psb_dpk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans
      complex(psb_dpk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_zspmv
    subroutine psb_zspmv_vect(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      type(psb_zspmat_type), intent(in)    :: a
      type(psb_z_vect_type), intent(inout) :: x
      type(psb_z_vect_type), intent(inout) :: y
      complex(psb_dpk_), intent(in)        :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      complex(psb_dpk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_zspmv_vect
  end interface

  interface psb_spsm
    subroutine psb_zspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,&
         & diag, n, jx, jy, work)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      type(psb_zspmat_type), intent(in)        :: t
      complex(psb_dpk_), intent(in), target    :: x(:,:)
      complex(psb_dpk_), intent(inout), target :: y(:,:)
      complex(psb_dpk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans, scale
      integer(psb_ipk_), optional, intent(in)            :: n, jx, jy
      integer(psb_ipk_), optional, intent(in)            :: choice
      complex(psb_dpk_), optional, intent(in), target :: diag(:)
      complex(psb_dpk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_zspsm
    subroutine psb_zspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,&
         & diag, work)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      type(psb_zspmat_type), intent(in)        :: t
      complex(psb_dpk_), intent(in), target    :: x(:)
      complex(psb_dpk_), intent(inout), target :: y(:)
      complex(psb_dpk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans, scale
      integer(psb_ipk_), optional, intent(in)            :: choice
      complex(psb_dpk_), optional, intent(in), target :: diag(:)
      complex(psb_dpk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_zspsv
    subroutine psb_zspsv_vect(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,&
         & diag, work)
      import :: psb_desc_type, psb_dpk_, psb_ipk_, &
           & psb_z_vect_type, psb_zspmat_type
      type(psb_zspmat_type), intent(inout)   :: t
      type(psb_z_vect_type), intent(inout)   :: x
      type(psb_z_vect_type), intent(inout)   :: y
      complex(psb_dpk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer(psb_ipk_), optional, intent(in)          :: choice
      type(psb_z_vect_type), intent(inout), optional :: diag
      complex(psb_dpk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_zspsv_vect
  end interface

  interface psb_gemlt
    subroutine psb_zmlt_vect(x,y,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_z_vect_type
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_z_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_zmlt_vect
  end interface

  interface psb_gediv
    subroutine psb_zdiv_vect(x,y,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_z_vect_type
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_z_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_zdiv_vect
    subroutine psb_zdiv_vect_check(x,y,desc_a,info,flag)
      import :: psb_desc_type, psb_ipk_, &
           & psb_z_vect_type
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_z_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in)                   :: flag
    end subroutine psb_zdiv_vect_check
  end interface

  interface psb_geinv
    subroutine psb_zinv_vect(x,y,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_z_vect_type
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_z_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_zinv_vect
    subroutine psb_zinv_vect_check(x,y,desc_a,info,flag)
      import :: psb_desc_type, psb_ipk_, &
           & psb_z_vect_type
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_z_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in)                   :: flag
    end subroutine psb_zinv_vect_check
  end interface

  interface psb_geabs
    subroutine psb_zabs_vect(x,y,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_z_vect_type
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_z_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_zabs_vect
  end interface

  interface psb_gecmp
    subroutine psb_zcmp_vect(x,c,z,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_z_vect_type, psb_dpk_
      type(psb_z_vect_type), intent (inout) :: x
      type(psb_z_vect_type), intent (inout) :: z
      real(psb_dpk_), intent(in)             :: c
      type(psb_desc_type), intent (in)        :: desc_a
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_zcmp_vect
  end interface

end module psb_z_psblas_mod
