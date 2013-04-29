!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
module psb_d_serial_mod
  use psb_const_mod
  use psb_error_mod

  interface psb_amax
    function psb_damax_s(n, x) result(val)
      import :: psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      real(psb_dpk_), intent(in)  :: x(:)
      real(psb_dpk_)                 :: val
    end function psb_damax_s
  end interface psb_amax

  interface psb_asum
    function psb_dasum_s(n, x) result(val)
      import :: psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      real(psb_dpk_), intent(in)  :: x(:)
      real(psb_dpk_)                 :: val
    end function psb_dasum_s
  end interface psb_asum

  interface psb_symbmm
    subroutine psb_dsymbmm(a,b,c,info)
      use psb_d_mat_mod, only : psb_dspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_dspmat_type), intent(in)  :: a,b
      type(psb_dspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_dsymbmm
    subroutine psb_dbase_symbmm(a,b,c,info)
      use psb_d_mat_mod, only : psb_d_base_sparse_mat, psb_d_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_dbase_symbmm
  end interface psb_symbmm

  interface psb_numbmm
    subroutine psb_dnumbmm(a,b,c)
      use psb_d_mat_mod, only : psb_dspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_dspmat_type), intent(in) :: a,b
      type(psb_dspmat_type), intent(inout)  :: c
    end subroutine psb_dnumbmm
    subroutine psb_dbase_numbmm(a,b,c)
      use psb_d_mat_mod, only : psb_d_base_sparse_mat, psb_d_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_dbase_numbmm
  end interface psb_numbmm

  interface psb_rwextd
    subroutine psb_drwextd(nr,a,info,b,rowscale)
      use psb_d_mat_mod, only : psb_dspmat_type
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                          :: nr
      type(psb_dspmat_type), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                          :: info
      type(psb_dspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_drwextd
    subroutine psb_dbase_rwextd(nr,a,info,b,rowscale)
      use psb_d_mat_mod, only : psb_d_base_sparse_mat
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                                :: nr
      class(psb_d_base_sparse_mat), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                                :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_dbase_rwextd
  end interface psb_rwextd


  interface psb_geprt
    subroutine psb_dgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      real(psb_dpk_), intent(in)    :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_dgeprtn2
    subroutine psb_dgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      real(psb_dpk_), intent(in)    :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_dgeprtn1
    subroutine psb_dgeprt2(iout,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      real(psb_dpk_), intent(in)     :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_dgeprt2
    subroutine psb_dgeprt1(iout,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      real(psb_dpk_), intent(in)     :: a(:)
      character(len=*), optional     :: head
    end subroutine psb_dgeprt1
  end interface psb_geprt

  interface psb_csprt
    module procedure psb_dcsprt, psb_dcsprtn
  end interface psb_csprt

  interface psb_spdot_srtd
    function psb_d_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1,nv2
      integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
      real(psb_dpk_), intent(in) :: v1(*),v2(*)
      real(psb_dpk_)      :: dot
    end function psb_d_spdot_srtd
  end interface psb_spdot_srtd


  interface psb_spge_dot
    function psb_d_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1
      integer(psb_ipk_), intent(in) :: iv1(*)
      real(psb_dpk_), intent(in) :: v1(*),v2(*)
      real(psb_dpk_)      :: dot
    end function psb_d_spge_dot
  end interface psb_spge_dot


  interface psb_aspxpby
    subroutine psb_d_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod, only : psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in)               :: nx
      integer(psb_ipk_), intent(in)               :: ix(:)
      real(psb_dpk_), intent (in)       :: x(:)
      real(psb_dpk_), intent (inout)    :: y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_d_aspxpby
  end interface psb_aspxpby

contains

  subroutine psb_dcsprt(iout,a,iv,head,ivr,ivc)
    use psb_d_mat_mod, only : psb_dspmat_type
    integer(psb_ipk_), intent(in)       :: iout
    type(psb_dspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_dcsprt

  subroutine psb_dcsprtn(fname,a,iv,head,ivr,ivc)
    use psb_d_mat_mod, only : psb_dspmat_type
    character(len=*), intent(in)  :: fname   
    type(psb_dspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_dcsprtn

end module psb_d_serial_mod

