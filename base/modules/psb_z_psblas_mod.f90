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
module psb_z_psblas_mod

  interface psb_gedot
    function psb_zdotv(x, y, desc_a,info) 
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      complex(psb_dpk_)                :: psb_zdotv
      complex(psb_dpk_), intent(in)    :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end function psb_zdotv
    function psb_zdot(x, y, desc_a, info, jx, jy) 
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      complex(psb_dpk_)                :: psb_zdot
      complex(psb_dpk_), intent(in)    :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, optional, intent(in)      :: jx, jy
      integer, intent(out)               :: info
    end function psb_zdot
  end interface


  interface psb_gedots
    subroutine  psb_zdotvs(res,x, y, desc_a, info) 
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      complex(psb_dpk_), intent(out)      :: res
      complex(psb_dpk_), intent(in)       :: x(:), y(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_zdotvs
    subroutine  psb_zmdots(res,x, y, desc_a,info) 
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      complex(psb_dpk_), intent(out)      :: res(:)
      complex(psb_dpk_), intent(in)       :: x(:,:), y(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer, intent(out)               :: info
    end subroutine psb_zmdots
  end interface

  interface psb_geaxpby
    subroutine psb_zaxpbyv(alpha, x, beta, y,&
         & desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      complex(psb_dpk_), intent (in)       ::  x(:)
      complex(psb_dpk_), intent (inout)    ::  y(:)
      complex(psb_dpk_), intent (in)       :: alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_zaxpbyv
    subroutine psb_zaxpby(alpha, x, beta, y,&
         & desc_a, info, n, jx, jy)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      complex(psb_dpk_), intent (in)       ::  x(:,:)
      complex(psb_dpk_), intent (inout)    ::  y(:,:)
      complex(psb_dpk_), intent (in)       ::  alpha, beta
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent(in) :: n, jx, jy
      integer, intent(out)                :: info
    end subroutine psb_zaxpby
  end interface

  interface psb_geamax
    function psb_zamax(x, desc_a, info, jx)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_)   psb_zamax
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_zamax
    function psb_zamaxv(x, desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_) psb_zamaxv
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_zamaxv
  end interface

  interface psb_geamaxs
    subroutine  psb_zamaxvs(res,x,desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)    :: x(:)
      type(psb_desc_type), intent (in)  :: desc_a
      integer, intent(out)              :: info
    end subroutine psb_zamaxvs
    subroutine  psb_zmamaxs(res,x,desc_a,info,jx)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_), intent (out)       :: res(:)
      complex(psb_dpk_), intent (in)     :: x(:,:)
      type(psb_desc_type), intent (in)   :: desc_a
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: jx
    end subroutine psb_zmamaxs
  end interface

  interface psb_geasum
    function psb_zasum(x, desc_a, info, jx)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_)   psb_zasum
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_zasum
    function psb_zasumv(x, desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_) psb_zasumv
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_zasumv
  end interface

  interface psb_geasums
    subroutine  psb_zasumvs(res,x,desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_zasumvs
    subroutine  psb_zmasum(res,x,desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_), intent (out)      :: res(:)
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_zmasum
  end interface


  interface psb_genrm2
    function psb_znrm2(x, desc_a, info, jx)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_)   psb_znrm2
      complex(psb_dpk_), intent (in)       :: x(:,:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, optional, intent (in)      :: jx
      integer, intent(out)                :: info
    end function psb_znrm2
    function psb_znrm2v(x, desc_a, info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_) psb_znrm2v
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end function psb_znrm2v
  end interface

  interface psb_genrm2s
    subroutine  psb_znrm2vs(res,x,desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      real(psb_dpk_), intent (out)      :: res
      complex(psb_dpk_), intent (in)       :: x(:)
      type(psb_desc_type), intent (in)    :: desc_a
      integer, intent(out)                :: info
    end subroutine psb_znrm2vs
  end interface


  interface psb_spnrmi
    function psb_znrmi(a, desc_a,info)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_z_sparse_mat
      real(psb_dpk_)                    :: psb_znrmi
      type(psb_z_sparse_mat), intent (in) :: a
      type(psb_desc_type), intent (in)   :: desc_a
      integer, intent(out)                :: info
    end function psb_znrmi
  end interface

  interface psb_spmm
    subroutine psb_zspmm(alpha, a, x, beta, y, desc_a, info,&
         &trans, k, jx, jy,work,doswap)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_z_sparse_mat
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
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_z_sparse_mat
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
    subroutine psb_zspsm(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, n, jx, jy, work)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_
      use psb_mat_mod, only : psb_z_sparse_mat
      type(psb_z_sparse_mat), intent(in)      :: t
      complex(psb_dpk_), intent(in)          :: x(:,:)
      complex(psb_dpk_), intent(inout)       :: y(:,:)
      complex(psb_dpk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer, optional, intent(in)          :: n, jx, jy
      integer, optional, intent(in)          :: choice
      complex(psb_dpk_), optional, intent(in), target :: diag(:)
      complex(psb_dpk_), optional, intent(inout), target :: work(:)
      integer, intent(out)               :: info
    end subroutine psb_zspsm
    subroutine psb_zspsv(alpha, t, x, beta, y,&
         & desc_a, info, trans, scale, choice,& 
         & diag, work)
      use psb_descriptor_type, only : psb_desc_type, psb_spk_, psb_dpk_ 
      use psb_mat_mod, only : psb_z_sparse_mat
      type(psb_z_sparse_mat), intent(in)      :: t
      complex(psb_dpk_), intent(in)          :: x(:)
      complex(psb_dpk_), intent(inout)       :: y(:)
      complex(psb_dpk_), intent(in)          :: alpha, beta
      type(psb_desc_type), intent(in)        :: desc_a
      character, optional, intent(in)        :: trans, scale
      integer, optional, intent(in)          :: choice
      complex(psb_dpk_), optional, intent(in), target :: diag(:)
      complex(psb_dpk_), optional, intent(inout), target :: work(:)
      integer, intent(out)                   :: info
    end subroutine psb_zspsv
  end interface

end module psb_z_psblas_mod
