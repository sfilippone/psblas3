!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
module psb_s_psblas_mod

  interface psb_gedot
    function psb_sdotv(x, y, desc_a,info) 
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_)                   :: psb_sdotv
      real(psb_spk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end function psb_sdotv
    function psb_sdot(x, y, desc_a, info, jx, jy) 
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_)                   :: psb_sdot
      real(psb_spk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, optional, intent(in)      :: jx, jy
      integer, intent(out)               :: info
    end function psb_sdot
  end interface


  interface psb_gedots
    subroutine  psb_sdotvs(res,x, y, desc_a, info) 
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent(out)      :: res
      real(psb_spk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_sdotvs
    subroutine  psb_smdots(res,x, y, desc_a,info) 
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent(out)      :: res(:)
      real(psb_spk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_smdots
  end interface

  interface psb_geaxpby
    subroutine psb_saxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent (in)       ::  x(:)
      real(psb_spk_), intent (inout)    ::  y(:)
      real(psb_spk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_saxpbyv
    subroutine psb_saxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent (in)       ::  x(:,:)
      real(psb_spk_), intent (inout)    ::  y(:,:)
      real(psb_spk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)  :: desc_a
      integer, optional, intent(in)     :: n, jx, jy
      integer, intent(out)              :: info
    end subroutine psb_saxpby
  end interface

  interface psb_geamax
    function psb_samax(x, desc_a, info, jx)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_)   psb_samax
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_samax
    function psb_samaxv(x, desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_) psb_samaxv
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_samaxv
  end interface

  interface psb_geamaxs
    subroutine  psb_samaxvs(res,x,desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in) :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_samaxvs
    subroutine  psb_smamaxs(res,x,desc_a,info,jx)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent (out)      :: res(:)
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
      integer, optional, intent(in)       :: jx
    end subroutine psb_smamaxs
  end interface

  interface psb_geasum
    function psb_sasum(x, desc_a, info, jx)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_)   psb_sasum
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_sasum
    function psb_sasumv(x, desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_) psb_sasumv
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_sasumv
  end interface

  interface psb_geasums
    subroutine  psb_sasumvs(res,x,desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_sasumvs
    subroutine  psb_smasum(res,x,desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent (out)      :: res(:)
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_smasum
  end interface


  interface psb_genrm2
    function psb_snrm2(x, desc_a, info, jx)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_)   psb_snrm2
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_snrm2
    function psb_snrm2v(x, desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_) psb_snrm2v
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_snrm2v
  end interface

  interface psb_genrm2s
    subroutine  psb_snrm2vs(res,x,desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_snrm2vs
  end interface


  interface psb_spnrmi
    function psb_snrmi(a, desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_sspmat_type
      real(psb_spk_)                      :: psb_snrmi
      type(psb_sspmat_type), intent (in) :: a
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_snrmi
  end interface

  interface psb_spmm
    subroutine psb_sspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_sspmat_type
      type(psb_sspmat_type), intent(in)   :: a
      real(psb_spk_), intent(inout)      :: x(:,:)
      real(psb_spk_), intent(inout)      :: y(:,:)
      real(psb_spk_), intent(in)         :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      real(psb_spk_), optional, intent(inout),target :: work(:)
      integer, optional, intent(in)        :: k, jx, jy
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_sspmm
    subroutine psb_sspmv(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_sspmat_type
      type(psb_sspmat_type), intent(in)   :: a
      real(psb_spk_), intent(inout)      :: x(:)
      real(psb_spk_), intent(inout)      :: y(:)
      real(psb_spk_), intent(in)         :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      real(psb_spk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_sspmv
  end interface

  interface psb_spsm
    subroutine psb_sspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, n, jx, jy, work)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_sspmat_type
      type(psb_sspmat_type), intent(in)   :: t
      real(psb_spk_), intent(in)           :: x(:,:)
      real(psb_spk_), intent(inout)        :: y(:,:)
      real(psb_spk_), intent(in)           :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans, scale
      integer, optional, intent(in)        :: n, jx, jy
      integer, optional, intent(in)        :: choice
      real(psb_spk_), optional, intent(in),target :: diag(:)
      real(psb_spk_), optional, intent(inout),target :: work(:)
      integer, intent(out)               :: info
    end subroutine psb_sspsm
    subroutine psb_sspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, work)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_   
      use psb_mat_mod, only : psb_sspmat_type
      type(psb_sspmat_type), intent(in)   :: t
      real(psb_spk_), intent(in)           :: x(:)
      real(psb_spk_), intent(inout)        :: y(:)
      real(psb_spk_), intent(in)           :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans, scale
      integer, optional, intent(in)        :: choice
      real(psb_spk_), optional, intent(in), target :: diag(:)
      real(psb_spk_), optional, intent(inout), target :: work(:)
      integer, intent(out)                   :: info
    end subroutine psb_sspsv
  end interface

end module psb_s_psblas_mod
