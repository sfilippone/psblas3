!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
module psb_psblas_mod

  interface psb_gedot
    function psb_sdotv(x, y, desc_a,info) 
      use psb_descriptor_type
      real(psb_spk_)                   :: psb_sdotv
      real(psb_spk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end function psb_sdotv
    function psb_sdot(x, y, desc_a, info, jx, jy) 
      use psb_descriptor_type
      real(psb_spk_)                   :: psb_sdot
      real(psb_spk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, optional, intent(in)      :: jx, jy
      integer, intent(out)               :: info
    end function psb_sdot
    function psb_ddotv(x, y, desc_a,info) 
      use psb_descriptor_type
      real(psb_dpk_)                   :: psb_ddotv
      real(psb_dpk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end function psb_ddotv
    function psb_ddot(x, y, desc_a, info, jx, jy) 
      use psb_descriptor_type
      real(psb_dpk_)                   :: psb_ddot
      real(psb_dpk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, optional, intent(in)      :: jx, jy
      integer, intent(out)               :: info
    end function psb_ddot
    function psb_cdotv(x, y, desc_a,info) 
      use psb_descriptor_type
      complex(psb_spk_)                :: psb_cdotv
      complex(psb_spk_), intent(in)    :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end function psb_cdotv
    function psb_cdot(x, y, desc_a, info, jx, jy) 
      use psb_descriptor_type
      complex(psb_spk_)                :: psb_cdot
      complex(psb_spk_), intent(in)    :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, optional, intent(in)      :: jx, jy
      integer, intent(out)               :: info
    end function psb_cdot
    function psb_zdotv(x, y, desc_a,info) 
      use psb_descriptor_type
      complex(psb_dpk_)                :: psb_zdotv
      complex(psb_dpk_), intent(in)    :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end function psb_zdotv
    function psb_zdot(x, y, desc_a, info, jx, jy) 
      use psb_descriptor_type
      complex(psb_dpk_)                :: psb_zdot
      complex(psb_dpk_), intent(in)    :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, optional, intent(in)      :: jx, jy
      integer, intent(out)               :: info
    end function psb_zdot
  end interface


  interface psb_gedots
    subroutine  psb_sdotvs(res,x, y, desc_a, info) 
      use psb_descriptor_type
      real(psb_spk_), intent(out)      :: res
      real(psb_spk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_sdotvs
    subroutine  psb_smdots(res,x, y, desc_a,info) 
      use psb_descriptor_type
      real(psb_spk_), intent(out)      :: res(:)
      real(psb_spk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_smdots
    subroutine  psb_ddotvs(res,x, y, desc_a, info) 
      use psb_descriptor_type
      real(psb_dpk_), intent(out)      :: res
      real(psb_dpk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_ddotvs
    subroutine  psb_dmdots(res,x, y, desc_a,info) 
      use psb_descriptor_type
      real(psb_dpk_), intent(out)      :: res(:)
      real(psb_dpk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_dmdots
    subroutine  psb_cdotvs(res,x, y, desc_a, info) 
      use psb_descriptor_type
      complex(psb_spk_), intent(out)      :: res
      complex(psb_spk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_cdotvs
    subroutine  psb_cmdots(res,x, y, desc_a,info) 
      use psb_descriptor_type
      complex(psb_spk_), intent(out)      :: res(:)
      complex(psb_spk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_cmdots
    subroutine  psb_zdotvs(res,x, y, desc_a, info) 
      use psb_descriptor_type
      complex(psb_dpk_), intent(out)      :: res
      complex(psb_dpk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_zdotvs
    subroutine  psb_zmdots(res,x, y, desc_a,info) 
      use psb_descriptor_type
      complex(psb_dpk_), intent(out)      :: res(:)
      complex(psb_dpk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_zmdots
  end interface

  interface psb_geaxpby
    subroutine psb_saxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      use psb_descriptor_type
      real(psb_spk_), intent (in)       ::  x(:)
      real(psb_spk_), intent (inout)    ::  y(:)
      real(psb_spk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_saxpbyv
    subroutine psb_saxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      use psb_descriptor_type
      real(psb_spk_), intent (in)       ::  x(:,:)
      real(psb_spk_), intent (inout)    ::  y(:,:)
      real(psb_spk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional :: n, jx, jy
      integer, intent(out)                :: info
    end subroutine psb_saxpby
    subroutine psb_daxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      use psb_descriptor_type
      real(psb_dpk_), intent (in)       ::  x(:)
      real(psb_dpk_), intent (inout)    ::  y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_daxpbyv
    subroutine psb_daxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      use psb_descriptor_type
      real(psb_dpk_), intent (in)       ::  x(:,:)
      real(psb_dpk_), intent (inout)    ::  y(:,:)
      real(psb_dpk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional :: n, jx, jy
      integer, intent(out)                :: info
    end subroutine psb_daxpby
    subroutine psb_caxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      use psb_descriptor_type
      complex(psb_spk_), intent (in)       ::  x(:)
      complex(psb_spk_), intent (inout)    ::  y(:)
      complex(psb_spk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_caxpbyv
    subroutine psb_caxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      use psb_descriptor_type
      complex(psb_spk_), intent (in)       ::  x(:,:)
      complex(psb_spk_), intent (inout)    ::  y(:,:)
      complex(psb_spk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional :: n, jx, jy
      integer, intent(out)                :: info
    end subroutine psb_caxpby
    subroutine psb_zaxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      use psb_descriptor_type
      complex(psb_dpk_), intent (in)       ::  x(:)
      complex(psb_dpk_), intent (inout)    ::  y(:)
      complex(psb_dpk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_zaxpbyv
    subroutine psb_zaxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      use psb_descriptor_type
      complex(psb_dpk_), intent (in)       ::  x(:,:)
      complex(psb_dpk_), intent (inout)    ::  y(:,:)
      complex(psb_dpk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional :: n, jx, jy
      integer, intent(out)                :: info
    end subroutine psb_zaxpby
  end interface

  interface psb_geamax
    function psb_samax(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_spk_)   psb_samax
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_samax
    function psb_samaxv(x, desc_a,info)
      use psb_descriptor_type
      real(psb_spk_) psb_samaxv
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_samaxv
    function psb_damax(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_dpk_)   psb_damax
      real(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_damax
    function psb_damaxv(x, desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_) psb_damaxv
      real(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_damaxv
    function psb_camax(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_spk_)   psb_camax
      complex(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_camax
    function psb_camaxv(x, desc_a,info)
      use psb_descriptor_type
      real(psb_spk_) psb_camaxv
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_camaxv
    function psb_zamax(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_dpk_)   psb_zamax
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_zamax
    function psb_zamaxv(x, desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_) psb_zamaxv
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_zamaxv
  end interface

  interface psb_geamaxs
    subroutine  psb_samaxvs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in) :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_samaxvs
    subroutine  psb_smamaxs(res,x,desc_a,info,jx)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res(:)
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
      integer, optional                   :: jx
    end subroutine psb_smamaxs
    subroutine  psb_damaxvs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res
      real(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in) :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_damaxvs
    subroutine  psb_dmamaxs(res,x,desc_a,info,jx)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res(:)
      real(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
      integer, optional                   :: jx
    end subroutine psb_dmamaxs
    subroutine  psb_camaxvs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in) :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_camaxvs
    subroutine  psb_cmamaxs(res,x,desc_a,info,jx)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res(:)
      complex(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
      integer, optional                   :: jx
    end subroutine psb_cmamaxs
    subroutine  psb_zamaxvs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in) :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_zamaxvs
    subroutine  psb_zmamaxs(res,x,desc_a,info,jx)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res(:)
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
      integer, optional                   :: jx
    end subroutine psb_zmamaxs
  end interface

  interface psb_geasum
    function psb_sasum(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_spk_)   psb_sasum
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_sasum
    function psb_sasumv(x, desc_a, info)
      use psb_descriptor_type
      real(psb_spk_) psb_sasumv
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_sasumv
    function psb_dasum(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_dpk_)   psb_dasum
      real(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_dasum
    function psb_dasumv(x, desc_a, info)
      use psb_descriptor_type
      real(psb_dpk_) psb_dasumv
      real(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_dasumv
    function psb_casum(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_spk_)   psb_casum
      complex(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_casum
    function psb_casumv(x, desc_a, info)
      use psb_descriptor_type
      real(psb_spk_) psb_casumv
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_casumv
    function psb_zasum(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_dpk_)   psb_zasum
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_zasum
    function psb_zasumv(x, desc_a, info)
      use psb_descriptor_type
      real(psb_dpk_) psb_zasumv
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_zasumv
  end interface

  interface psb_geasums
    subroutine  psb_sasumvs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_sasumvs
    subroutine  psb_smasum(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res(:)
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_smasum
    subroutine  psb_dasumvs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res
      real(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_dasumvs
    subroutine  psb_dmasum(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res(:)
      real(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)  :: desc_a
      integer, intent(out)              :: info
    end subroutine psb_dmasum
    subroutine  psb_casumvs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res
      complex(psb_spk_), intent (in)    :: x(:)
      type(psb_desc_type), intent (in)  :: desc_a
      integer, intent(out)              :: info
    end subroutine psb_casumvs
    subroutine  psb_cmasum(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res(:)
      complex(psb_spk_), intent (in)    :: x(:,:)
      type(psb_desc_type), intent (in)  :: desc_a
      integer, intent(out)              :: info
    end subroutine psb_cmasum
    subroutine  psb_zasumvs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_zasumvs
    subroutine  psb_zmasum(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res(:)
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_zmasum
  end interface


  interface psb_genrm2
    function psb_snrm2(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_spk_)   psb_snrm2
      real(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_snrm2
    function psb_snrm2v(x, desc_a, info)
      use psb_descriptor_type
      real(psb_spk_) psb_snrm2v
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_snrm2v
    function psb_dnrm2(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_dpk_)   psb_dnrm2
      real(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_dnrm2
    function psb_dnrm2v(x, desc_a, info)
      use psb_descriptor_type
      real(psb_dpk_) psb_dnrm2v
      real(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_dnrm2v
    function psb_cnrm2(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_spk_)   psb_snrm2
      complex(psb_spk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_cnrm2
    function psb_cnrm2v(x, desc_a, info)
      use psb_descriptor_type
      real(psb_spk_) psb_cnrm2v
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_cnrm2v
    function psb_znrm2(x, desc_a, info, jx)
      use psb_descriptor_type
      real(psb_dpk_)   psb_znrm2
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_znrm2
    function psb_znrm2v(x, desc_a, info)
      use psb_descriptor_type
      real(psb_dpk_) psb_znrm2v
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_znrm2v
  end interface

  interface psb_genrm2s
    subroutine  psb_snrm2vs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res
      real(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_snrm2vs
    subroutine  psb_dnrm2vs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res
      real(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_dnrm2vs
    subroutine  psb_cnrm2vs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_spk_), intent (out)      :: res
      complex(psb_spk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_cnrm2vs
    subroutine  psb_znrm2vs(res,x,desc_a,info)
      use psb_descriptor_type
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_znrm2vs
  end interface


  interface psb_spnrmi
    function psb_snrmi(a, desc_a,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      real(psb_spk_)                      :: psb_snrmi
      type(psb_s_sparse_mat), intent (in) :: a
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_snrmi
    function psb_dnrmi(a, desc_a,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      real(psb_dpk_)                      :: psb_dnrmi
      type(psb_d_sparse_mat), intent (in) :: a
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_dnrmi
    function psb_cnrmi(a, desc_a,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      real(psb_spk_)                    :: psb_cnrmi
      type(psb_c_sparse_mat), intent (in) :: a
      type(psb_desc_type), intent (in)   :: desc_a
      integer, intent(out)                :: info
    end function psb_cnrmi
    function psb_znrmi(a, desc_a,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      real(psb_dpk_)                    :: psb_znrmi
      type(psb_z_sparse_mat), intent (in) :: a
      type(psb_desc_type), intent (in)   :: desc_a
      integer, intent(out)                :: info
    end function psb_znrmi
  end interface

  interface psb_spmm
    subroutine psb_sspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_s_sparse_mat), intent(in)   :: a
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
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_s_sparse_mat), intent(in)   :: a
      real(psb_spk_), intent(inout)      :: x(:)
      real(psb_spk_), intent(inout)      :: y(:)
      real(psb_spk_), intent(in)         :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      real(psb_spk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_sspmv
    subroutine psb_dspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_d_sparse_mat), intent(in)    :: a
      real(psb_dpk_), intent(inout)        :: x(:,:)
      real(psb_dpk_), intent(inout)        :: y(:,:)
      real(psb_dpk_), intent(in)           :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      real(psb_dpk_), optional, intent(inout),target :: work(:)
      integer, optional, intent(in)        :: k, jx, jy
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_dspmm
    subroutine psb_dspmv(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_d_sparse_mat), intent(in)   :: a
      real(psb_dpk_), intent(inout)       :: x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      real(psb_dpk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)     :: desc_a
      character, optional, intent(in)     :: trans
      real(psb_dpk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_dspmv
    subroutine psb_cspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_c_sparse_mat), intent(in)    :: a
      complex(psb_spk_), intent(inout)     :: x(:,:)
      complex(psb_spk_), intent(inout)     :: y(:,:)
      complex(psb_spk_), intent(in)        :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      complex(psb_spk_), optional, intent(inout),target :: work(:)
      integer, optional, intent(in)        :: k, jx, jy
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_cspmm
    subroutine psb_cspmv(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_c_sparse_mat), intent(in)    :: a
      complex(psb_spk_), intent(inout)     :: x(:)
      complex(psb_spk_), intent(inout)     :: y(:)
      complex(psb_spk_), intent(in)        :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      complex(psb_spk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_cspmv
    subroutine psb_zspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_z_sparse_mat), intent(in)    :: a
      complex(psb_dpk_), intent(inout)     :: x(:,:)
      complex(psb_dpk_), intent(inout)     :: y(:,:)
      complex(psb_dpk_), intent(in)        :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      complex(psb_dpk_), optional, intent(inout),target :: work(:)
      integer, optional, intent(in)        :: k, jx, jy
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_zspmm
    subroutine psb_zspmv(alpha, a, x, beta, y,&
         & desc_a, info, trans, work,doswap)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_z_sparse_mat), intent(in)    :: a
      complex(psb_dpk_), intent(inout)     :: x(:)
      complex(psb_dpk_), intent(inout)     :: y(:)
      complex(psb_dpk_), intent(in)        :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans
      complex(psb_dpk_), optional, intent(inout),target :: work(:)
      logical, optional, intent(in)        :: doswap
      integer, intent(out)                 :: info
    end subroutine psb_zspmv
  end interface

  interface psb_spsm
    subroutine psb_sspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, n, jx, jy, work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_s_sparse_mat), intent(in)   :: t
      real(psb_spk_), intent(in)           :: x(:,:)
      real(psb_spk_), intent(inout)        :: y(:,:)
      real(psb_spk_), intent(in)           :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans, scale
      integer, optional, intent(in)        :: n, jx, jy
      integer, optional, intent(in)        :: choice
      real(psb_spk_), optional, intent(in),target :: work(:), diag(:)
      integer, intent(out)               :: info
    end subroutine psb_sspsm
    subroutine psb_sspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_s_sparse_mat), intent(in)   :: t
      real(psb_spk_), intent(in)           :: x(:)
      real(psb_spk_), intent(inout)        :: y(:)
      real(psb_spk_), intent(in)           :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans, scale
      integer, optional, intent(in)        :: choice
      real(psb_spk_), optional, intent(in),target :: work(:), diag(:)
      integer, intent(out)                   :: info
    end subroutine psb_sspsv
    subroutine psb_dspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, n, jx, jy, work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_d_sparse_mat), intent(in)    :: t
      real(psb_dpk_), intent(in)           :: x(:,:)
      real(psb_dpk_), intent(inout)        :: y(:,:)
      real(psb_dpk_), intent(in)           :: alpha, beta
      type(psb_desc_type), intent(in)      :: desc_a
      character, optional, intent(in)      :: trans, scale
      integer, optional, intent(in)        :: n, jx, jy
      integer, optional, intent(in)        :: choice
      real(psb_dpk_), optional, intent(in),target :: work(:), diag(:)
      integer, intent(out)               :: info
    end subroutine psb_dspsm
    subroutine psb_dspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_d_sparse_mat), intent(in)    :: t
      real(psb_dpk_), intent(in)             :: x(:)
      real(psb_dpk_), intent(inout)          :: y(:)
      real(psb_dpk_), intent(in)             :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer, optional, intent(in)          :: choice
      real(psb_dpk_), optional, intent(in),target :: work(:), diag(:)
      integer, intent(out)                   :: info
    end subroutine psb_dspsv
    subroutine psb_cspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, n, jx, jy, work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_c_sparse_mat), intent(in)      :: t
      complex(psb_spk_), intent(in)          :: x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      complex(psb_spk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer, optional, intent(in)          :: n, jx, jy
      integer, optional, intent(in)          :: choice
      complex(psb_spk_), optional, intent(in),target :: work(:), diag(:)
      integer, intent(out)               :: info
    end subroutine psb_cspsm
    subroutine psb_cspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_c_sparse_mat), intent(in)      :: t
      complex(psb_spk_), intent(in)          :: x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      complex(psb_spk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer, optional, intent(in)          :: choice
      complex(psb_spk_), optional, intent(in),target :: work(:), diag(:)
      integer, intent(out)                   :: info
    end subroutine psb_cspsv
    subroutine psb_zspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, n, jx, jy, work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_z_sparse_mat), intent(in)      :: t
      complex(psb_dpk_), intent(in)          :: x(:,:)
      complex(psb_dpk_), intent(inout)       :: y(:,:)
      complex(psb_dpk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer, optional, intent(in)          :: n, jx, jy
      integer, optional, intent(in)          :: choice
      complex(psb_dpk_), optional, intent(in),target :: work(:), diag(:)
      integer, intent(out)               :: info
    end subroutine psb_zspsm
    subroutine psb_zspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, work)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_mat_mod
      type(psb_z_sparse_mat), intent(in)      :: t
      complex(psb_dpk_), intent(in)          :: x(:)
      complex(psb_dpk_), intent(inout)       :: y(:)
      complex(psb_dpk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer, optional, intent(in)          :: choice
      complex(psb_dpk_), optional, intent(in),target :: work(:), diag(:)
      integer, intent(out)                   :: info
    end subroutine psb_zspsv
  end interface

end module psb_psblas_mod
