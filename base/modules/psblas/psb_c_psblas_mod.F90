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
module psb_c_psblas_mod
  use psb_desc_mod, only : psb_desc_type, psb_spk_, psb_ipk_
  use psb_c_vect_mod, only : psb_c_vect_type
  use psb_c_mat_mod, only : psb_cspmat_type

  interface psb_gedot
    function psb_cdot_vect(x, y, desc_a,info,global) result(res)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      complex(psb_spk_)                    :: res
      type(psb_c_vect_type), intent(inout) :: x, y
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_cdot_vect
    function psb_cdotv(x, y, desc_a,info,global) 
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      complex(psb_spk_)                :: psb_cdotv
      complex(psb_spk_), intent(in)    :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_cdotv
    function psb_cdot(x, y, desc_a, info, jx, jy,global) 
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      complex(psb_spk_)                :: psb_cdot
      complex(psb_spk_), intent(in)    :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), optional, intent(in)      :: jx, jy
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_cdot
  end interface


  interface psb_gedots
    subroutine  psb_cdotvs(res,x, y, desc_a, info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type 
      complex(psb_spk_), intent(out)      :: res
      complex(psb_spk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end subroutine psb_cdotvs
    subroutine  psb_cmdots(res,x, y, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type 
      complex(psb_spk_), intent(out)      :: res(:)
      complex(psb_spk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end subroutine psb_cmdots
  end interface

  interface psb_geaxpby
    subroutine psb_caxpby_vect(alpha, x, beta, y,&
         & desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      type(psb_c_vect_type), intent (inout) :: x
      type(psb_c_vect_type), intent (inout) :: y
      complex(psb_spk_), intent (in)        :: alpha, beta
      type(psb_desc_type), intent (in)      :: desc_a
      integer(psb_ipk_), intent(out)        :: info
    end subroutine psb_caxpby_vect
    subroutine psb_caxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      complex(psb_spk_), intent (in)       ::  x(:)
      complex(psb_spk_), intent (inout)    ::  y(:)
      complex(psb_spk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_caxpbyv
    subroutine psb_caxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      complex(psb_spk_), intent (in)       ::  x(:,:)
      complex(psb_spk_), intent (inout)    ::  y(:,:)
      complex(psb_spk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent(in) :: n, jx, jy
      integer(psb_ipk_), intent(out)      :: info
    end subroutine psb_caxpby
  end interface

  interface psb_geamax
    function psb_camax(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_)   psb_camax
      complex(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_camax
    function psb_camaxv(x, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_) psb_camaxv
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_camaxv
    function psb_camax_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_)                        :: res
      type(psb_c_vect_type), intent (inout) :: x
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_camax_vect
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_genrmi
    procedure psb_camax, psb_camaxv, psb_camax_vect
  end interface 
  interface psb_normi
    procedure psb_camax, psb_camaxv, psb_camax_vect
  end interface 
#endif

  interface psb_geamaxs
    subroutine  psb_camaxvs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_), intent (out)      :: res
      complex(psb_spk_), intent (in)    :: x(:)
      type(psb_desc_type), intent (in)  :: desc_a
      integer(psb_ipk_), intent(out)    :: info
      logical, intent(in), optional     :: global
    end subroutine psb_camaxvs
    subroutine  psb_cmamaxs(res,x,desc_a,info,jx,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_), intent (out)       :: res(:)
      complex(psb_spk_), intent (in)     :: x(:,:)
      type(psb_desc_type), intent (in)   :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      integer(psb_ipk_), optional, intent(in)      :: jx
      logical, intent(in), optional      :: global
    end subroutine psb_cmamaxs
  end interface

  interface psb_geasum
    function psb_casum_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_)                        :: res
      type(psb_c_vect_type), intent (inout) :: x
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_casum_vect
    function psb_casum(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_)   psb_casum
      complex(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_casum
    function psb_casumv(x, desc_a, info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_) psb_casumv
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_casumv
  end interface

  interface psb_geasums
    subroutine  psb_casumvs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_), intent (out)      :: res
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_casumvs
    subroutine  psb_cmasum(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_), intent (out)      :: res(:)
      complex(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_cmasum
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_genrm1
    procedure psb_casum, psb_casumv, psb_casum_vect
  end interface 
  interface psb_norm1
    procedure psb_casum, psb_casumv, psb_casum_vect
  end interface 
#endif

  interface psb_genrm2
    function psb_cnrm2(x, desc_a, info, jx,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_)   psb_cnrm2
      complex(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), optional, intent (in)      :: jx
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_cnrm2
    function psb_cnrm2v(x, desc_a, info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_) psb_cnrm2v
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)                :: info
      logical, intent(in), optional        :: global
    end function psb_cnrm2v
    function psb_cnrm2_vect(x, desc_a, info,global) result(res)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_)                      :: res
      type(psb_c_vect_type), intent (inout)   :: x
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end function psb_cnrm2_vect
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_norm2
    procedure psb_cnrm2, psb_cnrm2v, psb_cnrm2_vect
  end interface 
#endif

  interface psb_genrm2s
    subroutine  psb_cnrm2vs(res,x,desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_), intent (out)      :: res
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: global
    end subroutine psb_cnrm2vs
  end interface


  interface psb_spnrmi
    function psb_cnrmi(a, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_)                    :: psb_cnrmi
      type(psb_cspmat_type), intent (in) :: a
      type(psb_desc_type), intent (in)   :: desc_a
      integer(psb_ipk_), intent(out)     :: info
      logical, intent(in), optional      :: global
    end function psb_cnrmi
  end interface

#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_normi
    procedure psb_cnrmi
  end interface 
#endif

  interface psb_spnrm1
    function psb_cspnrm1(a, desc_a,info,global)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      real(psb_spk_)                     :: psb_cspnrm1
      type(psb_cspmat_type), intent (in) :: a
      type(psb_desc_type), intent (in)     :: desc_a
      integer(psb_ipk_), intent(out)       :: info
      logical, intent(in), optional        :: global
    end function psb_cspnrm1
  end interface
  
#if ! defined(HAVE_BUGGY_GENERICS)
  interface psb_norm1
    procedure psb_cspnrm1
  end interface 
#endif

  interface psb_spmm
    subroutine psb_cspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      type(psb_cspmat_type), intent(in)        :: a
      complex(psb_spk_), intent(inout), target :: x(:,:)
      complex(psb_spk_), intent(inout), target :: y(:,:)
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans
      complex(psb_spk_), optional, intent(inout),target :: work(:)
      integer(psb_ipk_), optional, intent(in)        :: k, jx, jy
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_cspmm
    subroutine psb_cspmv(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      type(psb_cspmat_type), intent(in)        :: a
      complex(psb_spk_), intent(inout), target :: x(:)
      complex(psb_spk_), intent(inout), target :: y(:)
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans
      complex(psb_spk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_cspmv
    subroutine psb_cspmv_vect(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      type(psb_cspmat_type), intent(in)    :: a
      type(psb_c_vect_type), intent(inout) :: x
      type(psb_c_vect_type), intent(inout) :: y
      complex(psb_spk_), intent(in)        :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      complex(psb_spk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer(psb_ipk_), intent(out)                 :: info
    end subroutine psb_cspmv_vect
  end interface

  interface psb_spsm
    subroutine psb_cspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, n, jx, jy, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      type(psb_cspmat_type), intent(in)        :: t
      complex(psb_spk_), intent(in), target    :: x(:,:)
      complex(psb_spk_), intent(inout), target :: y(:,:)
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans, scale
      integer(psb_ipk_), optional, intent(in)            :: n, jx, jy
      integer(psb_ipk_), optional, intent(in)            :: choice
      complex(psb_spk_), optional, intent(in), target :: diag(:)
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_cspsm
    subroutine psb_cspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      type(psb_cspmat_type), intent(in)        :: t
      complex(psb_spk_), intent(in), target    :: x(:)
      complex(psb_spk_), intent(inout), target :: y(:)
      complex(psb_spk_), intent(in)            :: alpha, beta
      type(psb_desc_type), intent(in)          :: desc_a
      character, optional, intent(in)          :: trans, scale
      integer(psb_ipk_), optional, intent(in)            :: choice
      complex(psb_spk_), optional, intent(in), target :: diag(:)
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_cspsv
    subroutine psb_cspsv_vect(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, work)
      import :: psb_desc_type, psb_spk_, psb_ipk_, &
           & psb_c_vect_type, psb_cspmat_type
      type(psb_cspmat_type), intent(inout)   :: t
      type(psb_c_vect_type), intent(inout)   :: x
      type(psb_c_vect_type), intent(inout)   :: y
      complex(psb_spk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer(psb_ipk_), optional, intent(in)          :: choice
      type(psb_c_vect_type), intent(inout), optional :: diag
      complex(psb_spk_), optional, intent(inout), target :: work(:)
      integer(psb_ipk_), intent(out)                   :: info
    end subroutine psb_cspsv_vect
  end interface

end module psb_c_psblas_mod
