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
module psb_s_psblas_mod
  use psb_desc_mod, only : psb_desc_type, psb_spk_, psb_ipk_
  use psb_s_vect_mod, only : psb_s_vect_type
  use psb_s_mat_mod, only : psb_sspmat_type

  interface psb_gedot
    function psb_sdot_vect(x, y, desc_a,info,global) result(res)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)                    :: res
      type(psb_s_vect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_sdot_vect
    function psb_sdotv(x, y, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)                :: psb_sdotv
      real(psb_spk_), intent(in)    :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_sdotv
    function psb_sdot(x, y, desc_a, info, jx, jy,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)                :: psb_sdot
      real(psb_spk_), intent(in)    :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), optional, intent(in)      :: jx, jy
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_sdot
  end interface


  interface psb_gedots
    subroutine  psb_sdotvs(res,x, y, desc_a, info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent(out)      :: res
      real(psb_spk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end subroutine psb_sdotvs
    subroutine  psb_smdots(res,x, y, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent(out)      :: res(:)
      real(psb_spk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end subroutine psb_smdots
  end interface

  interface psb_geaxpby
    subroutine psb_saxpby_vect(alpha, x, beta, y,&
         & desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_s_vect_type), intent (inout) :: y
      real(psb_spk_), intent (in)        :: alpha, beta
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_saxpby_vect
    subroutine psb_saxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent (in)       ::  x(:)
      real(psb_spk_), intent (inout)    ::  y(:)
      real(psb_spk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_saxpbyv
    subroutine psb_saxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent (in)       ::  x(:,:)
      real(psb_spk_), intent (inout)    ::  y(:,:)
      real(psb_spk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent(in) :: n, jx, jy
      integer(psb_ipk_), intent(out)      :: info
    end subroutine psb_saxpby
  end interface

  interface psb_geamax
    function psb_samax(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)   psb_samax
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_samax
    function psb_samaxv(x, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_) psb_samaxv
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_samaxv
    function psb_samax_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)                        :: res
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_samax_vect
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_genrmi
    procedure psb_samax, psb_samaxv, psb_samax_vect
  end interface
  interface psb_normi
    procedure psb_samax, psb_samaxv, psb_samax_vect
  end interface
#endif

  interface psb_geamaxs
    subroutine  psb_samaxvs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)    :: x(:)
      type(psb_desc_type), intent (in)  :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical, intent(in), optional     :: global
    end subroutine psb_samaxvs
    subroutine  psb_smamaxs(res,x,desc_a,info,jx,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent (out)       :: res(:)
      real(psb_spk_), intent (in)     :: x(:,:)
      type(psb_desc_type), intent (in)   :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      integer(psb_ipk_), optional, intent(in)      :: jx
      logical, intent(in), optional      :: global
    end subroutine psb_smamaxs
  end interface

  interface psb_geasum
    function psb_sasum_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)                        :: res
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_sasum_vect
    function psb_sasum(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)   psb_sasum
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_sasum
    function psb_sasumv(x, desc_a, info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_) psb_sasumv
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_sasumv
  end interface

  interface psb_geasums
    subroutine  psb_sasumvs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_sasumvs
    subroutine  psb_smasum(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent (out)      :: res(:)
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_smasum
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_genrm1
    procedure psb_sasum, psb_sasumv, psb_sasum_vect
  end interface
  interface psb_norm1
    procedure psb_sasum, psb_sasumv, psb_sasum_vect
  end interface
#endif

  interface psb_genrm2
    function psb_snrm2(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)   psb_snrm2
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_snrm2
    function psb_snrm2v(x, desc_a, info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_) psb_snrm2v
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)                :: info
      logical, intent(in), optional        :: global
    end function psb_snrm2v
    function psb_snrm2_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)                      :: res
      type(psb_s_vect_type), intent (inout)   :: x
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_snrm2_vect
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_norm2
    procedure psb_snrm2, psb_snrm2v, psb_snrm2_vect
  end interface
#endif

  interface psb_genrm2s
    subroutine  psb_snrm2vs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_snrm2vs
  end interface


  interface psb_spnrmi
    function psb_snrmi(a, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)                    :: psb_snrmi
      type(psb_sspmat_type), intent (in) :: a
      type(psb_desc_type), intent (in)   :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_snrmi
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_normi
    procedure psb_snrmi
  end interface
#endif

  interface psb_spnrm1
    function psb_sspnrm1(a, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      real(psb_spk_)                     :: psb_sspnrm1
      type(psb_sspmat_type), intent (in) :: a
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_sspnrm1
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_norm1
    procedure psb_sspnrm1
  end interface
#endif

  interface psb_spmm
    subroutine psb_sspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      type(psb_sspmat_type), intent(in)        :: a
      real(psb_spk_), intent(inout), target :: x(:,:)
      real(psb_spk_), intent(inout), target :: y(:,:)
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans
      real(psb_spk_), optional, intent(inout),target :: work(:)
      integer(psb_ipk_), optional, intent(in)        :: k, jx, jy
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_sspmm
    subroutine psb_sspmv(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      type(psb_sspmat_type), intent(in)        :: a
      real(psb_spk_), intent(inout), target :: x(:)
      real(psb_spk_), intent(inout), target :: y(:)
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans
      real(psb_spk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_sspmv
    subroutine psb_sspmv_vect(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      type(psb_sspmat_type), intent(in)    :: a
      type(psb_s_vect_type), intent(inout) :: x
      type(psb_s_vect_type), intent(inout) :: y
      real(psb_spk_), intent(in)        :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      real(psb_spk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_sspmv_vect
  end interface

  interface psb_spsm
    subroutine psb_sspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,&
         & diag, n, jx, jy, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      type(psb_sspmat_type), intent(in)        :: t
      real(psb_spk_), intent(in), target    :: x(:,:)
      real(psb_spk_), intent(inout), target :: y(:,:)
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans, scale
      integer(psb_ipk_), optional, intent(in)            :: n, jx, jy
      integer(psb_ipk_), optional, intent(in)            :: choice
      real(psb_spk_), optional, intent(in), target :: diag(:)
      real(psb_spk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_sspsm
    subroutine psb_sspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,&
         & diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      type(psb_sspmat_type), intent(in)        :: t
      real(psb_spk_), intent(in), target    :: x(:)
      real(psb_spk_), intent(inout), target :: y(:)
      real(psb_spk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans, scale
      integer(psb_ipk_), optional, intent(in)            :: choice
      real(psb_spk_), optional, intent(in), target :: diag(:)
      real(psb_spk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_sspsv
    subroutine psb_sspsv_vect(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,&
         & diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_s_vect_type, psb_sspmat_type
      type(psb_sspmat_type), intent(inout)   :: t
      type(psb_s_vect_type), intent(inout)   :: x
      type(psb_s_vect_type), intent(inout)   :: y
      real(psb_spk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer(psb_ipk_), optional, intent(in)          :: choice
      type(psb_s_vect_type), intent(inout), optional :: diag
      real(psb_spk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_sspsv_vect
  end interface

  interface psb_gemlt
    subroutine psb_smlt_vect(x,y,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_s_vect_type
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_s_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_smlt_vect
  end interface

  interface psb_gediv
    subroutine psb_sdiv_vect(x,y,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_s_vect_type
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_s_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_sdiv_vect
    subroutine psb_sdiv_vect_check(x,y,desc_a,info,flag)
      import :: psb_desc_type, psb_ipk_, &
           & psb_s_vect_type
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_s_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in)                   :: flag
    end subroutine psb_sdiv_vect_check
  end interface

  interface psb_geinv
    subroutine psb_sinv_vect(x,y,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_s_vect_type
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_s_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_sinv_vect
    subroutine psb_sinv_vect_check(x,y,desc_a,info,flag)
      import :: psb_desc_type, psb_ipk_, &
           & psb_s_vect_type
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_s_vect_type), intent (inout) :: y
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in)                   :: flag
    end subroutine psb_sinv_vect_check
  end interface

  interface psb_gecmp
    subroutine psb_scmp_vect(x,c,z,desc_a,info)
      import :: psb_desc_type, psb_ipk_, &
           & psb_s_vect_type, psb_spk_
      type(psb_s_vect_type), intent (inout) :: x
      type(psb_s_vect_type), intent (inout) :: z
      real(psb_spk_), intent(in)             :: c
      type(psb_desc_type), intent (in)        :: desc_a
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_scmp_vect
  end interface

end module psb_s_psblas_mod
