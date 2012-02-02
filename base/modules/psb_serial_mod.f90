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
module psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_string_mod
  use psb_sort_mod

  use psi_serial_mod, &
       & psb_gth => psi_gth,&
       & psb_sct => psi_sct

  interface psb_nrm1
    module procedure psb_snrm1, psb_dnrm1, psb_cnrm1, psb_znrm1
  end interface psb_nrm1

  interface psb_amax
    function psb_samax_s(n, x) result(val)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      real(psb_spk_), intent(in)  :: x(:)
      real(psb_spk_)                 :: val
    end function psb_samax_s
    function psb_damax_s(n, x) result(val)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      real(psb_dpk_), intent(in)  :: x(:)
      real(psb_dpk_)                 :: val
    end function psb_damax_s
    function psb_camax_s(n, x) result(val)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      complex(psb_spk_), intent(in)  :: x(:)
      real(psb_spk_)                 :: val
    end function psb_camax_s
    function psb_zamax_s(n, x) result(val)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      complex(psb_dpk_), intent(in)  :: x(:)
      real(psb_dpk_)                 :: val
    end function psb_zamax_s
  end interface psb_amax

  interface psb_asum
    function psb_sasum_s(n, x) result(val)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      real(psb_spk_), intent(in)  :: x(:)
      real(psb_spk_)                 :: val
    end function psb_sasum_s
    function psb_dasum_s(n, x) result(val)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      real(psb_dpk_), intent(in)  :: x(:)
      real(psb_dpk_)                 :: val
    end function psb_dasum_s
    function psb_casum_s(n, x) result(val)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      complex(psb_spk_), intent(in)  :: x(:)
      real(psb_spk_)                 :: val
    end function psb_casum_s
    function psb_zasum_s(n, x) result(val)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      complex(psb_dpk_), intent(in)  :: x(:)
      real(psb_dpk_)                 :: val
    end function psb_zasum_s
  end interface psb_asum

  interface psb_symbmm
    subroutine psb_ssymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_sspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_sspmat_type), intent(in)  :: a,b
      type(psb_sspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_ssymbmm
    subroutine psb_sbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_s_base_sparse_mat, psb_s_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_s_base_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_sbase_symbmm
    subroutine psb_dsymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_dspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_dspmat_type), intent(in)  :: a,b
      type(psb_dspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_dsymbmm
    subroutine psb_dbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_d_base_sparse_mat, psb_d_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_dbase_symbmm
    subroutine psb_csymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_cspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_cspmat_type), intent(in)  :: a,b
      type(psb_cspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_csymbmm
    subroutine psb_cbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_c_base_sparse_mat, psb_c_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_c_base_sparse_mat), intent(in) :: a,b
      type(psb_c_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_cbase_symbmm
    subroutine psb_zsymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_zspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_zspmat_type), intent(in)  :: a,b
      type(psb_zspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_zsymbmm
    subroutine psb_zbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_z_base_sparse_mat, psb_z_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_z_base_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_zbase_symbmm
  end interface psb_symbmm

  interface psb_numbmm
    subroutine psb_snumbmm(a,b,c)
      use psb_mat_mod, only : psb_sspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_sspmat_type), intent(in) :: a,b
      type(psb_sspmat_type), intent(inout)  :: c
    end subroutine psb_snumbmm
    subroutine psb_sbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_s_base_sparse_mat, psb_s_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_s_base_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_sbase_numbmm
    subroutine psb_dnumbmm(a,b,c)
      use psb_mat_mod, only : psb_dspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_dspmat_type), intent(in) :: a,b
      type(psb_dspmat_type), intent(inout)  :: c
    end subroutine psb_dnumbmm
    subroutine psb_dbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_d_base_sparse_mat, psb_d_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_dbase_numbmm
    subroutine psb_cnumbmm(a,b,c)
      use psb_mat_mod, only : psb_cspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_cspmat_type), intent(in) :: a,b
      type(psb_cspmat_type), intent(inout)  :: c
    end subroutine psb_cnumbmm
    subroutine psb_cbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_c_base_sparse_mat, psb_c_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_c_base_sparse_mat), intent(in) :: a,b
      type(psb_c_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_cbase_numbmm
    subroutine psb_znumbmm(a,b,c)
      use psb_mat_mod, only : psb_zspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_zspmat_type), intent(in) :: a,b
      type(psb_zspmat_type), intent(inout)  :: c
    end subroutine psb_znumbmm
    subroutine psb_zbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_z_base_sparse_mat, psb_z_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_z_base_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_zbase_numbmm
  end interface psb_numbmm

  interface psb_rwextd
    subroutine psb_srwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_sspmat_type
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                          :: nr
      type(psb_sspmat_type), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                          :: info
      type(psb_sspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_srwextd
    subroutine psb_sbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_s_base_sparse_mat
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                                :: nr
      class(psb_s_base_sparse_mat), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                                :: info
      class(psb_s_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_sbase_rwextd
    subroutine psb_drwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_dspmat_type
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                          :: nr
      type(psb_dspmat_type), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                          :: info
      type(psb_dspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_drwextd
    subroutine psb_dbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_d_base_sparse_mat
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                                :: nr
      class(psb_d_base_sparse_mat), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                                :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_dbase_rwextd
    subroutine psb_crwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_cspmat_type
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                          :: nr
      type(psb_cspmat_type), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                          :: info
      type(psb_cspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_crwextd
    subroutine psb_cbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_c_base_sparse_mat
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                                :: nr
      class(psb_c_base_sparse_mat), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                                :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_cbase_rwextd
    subroutine psb_zrwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_zspmat_type
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                          :: nr
      type(psb_zspmat_type), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                          :: info
      type(psb_zspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_zrwextd
    subroutine psb_zbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_z_base_sparse_mat
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                                :: nr
      class(psb_z_base_sparse_mat), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                                :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_zbase_rwextd
  end interface psb_rwextd


  interface psb_geprt
    subroutine psb_sgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_spk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      real(psb_spk_), intent(in)    :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_sgeprtn2
    subroutine psb_sgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_spk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      real(psb_spk_), intent(in)    :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_sgeprtn1
    subroutine psb_sgeprt2(iout,a,head)
      use psb_const_mod, only : psb_spk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      real(psb_spk_), intent(in)     :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_sgeprt2
    subroutine psb_sgeprt1(iout,a,head)
      use psb_const_mod, only : psb_spk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      real(psb_spk_), intent(in)     :: a(:)
      character(len=*), optional     :: head
    end subroutine psb_sgeprt1
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
    subroutine psb_cgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_spk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      complex(psb_spk_), intent(in) :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_cgeprtn2
    subroutine psb_cgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_spk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      complex(psb_spk_), intent(in) :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_cgeprtn1
    subroutine psb_cgeprt2(iout,a,head)
      use psb_const_mod, only : psb_spk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      complex(psb_spk_), intent(in)  :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_cgeprt2
    subroutine psb_cgeprt1(iout,a,head)
      use psb_const_mod, only : psb_spk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      complex(psb_spk_), intent(in)  :: a(:)
      character(len=*), optional     :: head
    end subroutine psb_cgeprt1
    subroutine psb_zgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      complex(psb_dpk_), intent(in) :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_zgeprtn2
    subroutine psb_zgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      complex(psb_dpk_), intent(in) :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_zgeprtn1
    subroutine psb_zgeprt2(iout,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      complex(psb_dpk_), intent(in)  :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_zgeprt2
    subroutine psb_zgeprt1(iout,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      complex(psb_dpk_), intent(in)  :: a(:)
      character(len=*), optional     :: head
    end subroutine psb_zgeprt1
  end interface psb_geprt

  interface psb_csprt
    module procedure psb_scsprt, psb_scsprtn, psb_dcsprt, psb_dcsprtn, &
         & psb_ccsprt, psb_ccsprtn, psb_zcsprt, psb_zcsprtn
  end interface psb_csprt

  interface psb_spdot_srtd
    function psb_s_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1,nv2
      integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
      real(psb_spk_), intent(in) :: v1(*),v2(*)
      real(psb_spk_)      :: dot
    end function psb_s_spdot_srtd

    function psb_d_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1,nv2
      integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
      real(psb_dpk_), intent(in) :: v1(*), v2(*)
      real(psb_dpk_)      :: dot
    end function psb_d_spdot_srtd

    function psb_c_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1,nv2
      integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
      complex(psb_spk_), intent(in) :: v1(*),v2(*)
      complex(psb_spk_)      :: dot
    end function psb_c_spdot_srtd

    function psb_z_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1,nv2
      integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
      complex(psb_dpk_), intent(in) :: v1(*),v2(*)
      complex(psb_dpk_)      :: dot
    end function psb_z_spdot_srtd
  end interface psb_spdot_srtd


  interface psb_spge_dot
    function psb_s_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1
      integer(psb_ipk_), intent(in) :: iv1(*)
      real(psb_spk_), intent(in) :: v1(*),v2(*)
      real(psb_spk_)      :: dot
    end function psb_s_spge_dot

    function psb_d_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1
      integer(psb_ipk_), intent(in) :: iv1(*)
      real(psb_dpk_), intent(in) :: v1(*),v2(*)
      real(psb_dpk_)      :: dot
    end function psb_d_spge_dot

    function psb_c_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1
      integer(psb_ipk_), intent(in) :: iv1(*)
      complex(psb_spk_), intent(in) :: v1(*),v2(*)
      complex(psb_spk_)      :: dot
    end function psb_c_spge_dot

    function psb_z_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1
      integer(psb_ipk_), intent(in) :: iv1(*)
      complex(psb_dpk_), intent(in) :: v1(*),v2(*)
      complex(psb_dpk_)      :: dot
    end function psb_z_spge_dot
  end interface psb_spge_dot


  interface psb_nspaxpby
    subroutine psb_d_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(out)              :: nz
      integer(psb_ipk_), intent(out)              :: iz(:)
      real(psb_dpk_), intent (out)      :: z(:)
      integer(psb_ipk_), intent(in)               :: nx, ny
      integer(psb_ipk_), intent(in)               :: ix(:), iy(:)
      real(psb_dpk_), intent (in)       :: x(:), y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_d_nspaxpby
  end interface psb_nspaxpby

  interface psb_aspxpby
    subroutine psb_s_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)               :: nx
      integer(psb_ipk_), intent(in)               :: ix(:)
      real(psb_spk_), intent (in)       :: x(:)
      real(psb_spk_), intent (inout)    :: y(:)
      real(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_s_aspxpby
    subroutine psb_d_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)               :: nx
      integer(psb_ipk_), intent(in)               :: ix(:)
      real(psb_dpk_), intent (in)       :: x(:)
      real(psb_dpk_), intent (inout)    :: y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_d_aspxpby
    subroutine psb_c_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)               :: nx
      integer(psb_ipk_), intent(in)               :: ix(:)
      complex(psb_spk_), intent (in)    :: x(:)
      complex(psb_spk_), intent (inout) :: y(:)
      complex(psb_spk_), intent (in)    :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_c_aspxpby
    subroutine psb_z_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(in)               :: nx
      integer(psb_ipk_), intent(in)               :: ix(:)
      complex(psb_dpk_), intent (in)    :: x(:)
      complex(psb_dpk_), intent (inout) :: y(:)
      complex(psb_dpk_), intent (in)    :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_z_aspxpby
  end interface psb_aspxpby

contains

  elemental function psb_snrm1(x) result(res)
    real(psb_spk_), intent(in)  :: x
    real(psb_spk_)              :: val
    val = abs( x )
  end function psb_snrm1

  elemental function psb_dnrm1(x) result(res)
    real(psb_dpk_), intent(in) :: x
    real(psb_dpk_)             :: val
    val = abs(  x )
  end function psb_dnrm1

  elemental function psb_cnrm1(x) result(res)
    complex(psb_spk_), intent(in)  :: x
    real(psb_spk_)                 :: val
    val = abs( real( x ) ) + abs( aimag( x ) )  
  end function psb_cnrm1

  elemental function psb_znrm1(x) result(res)
    complex(psb_dpk_), intent(in)  :: x
    real(psb_dpk_)                 :: val
    val = abs( real( x ) ) + abs( aimag( x ) )  
  end function psb_znrm1


  subroutine crot( n, cx, incx, cy, incy, c, s )
    !
    !  -- lapack auxiliary routine (version 3.0) --
    !     univ. of tennessee, univ. of california berkeley, nag ltd.,
    !     courant institute, argonne national lab, and rice university
    !     october 31, 1992
    !
    !     .. scalar arguments ..
    integer(psb_mpik_) :: incx, incy, n
    real(psb_spk_)    c
    complex(psb_spk_)   s
    !     ..
    !     .. array arguments ..
    complex(psb_spk_) cx( * ), cy( * )
    !     ..
    !
    !  purpose
    !  == = ====
    !
    !  zrot   applies a plane rotation, where the cos (c) is real and the
    !  sin (s) is complex, and the vectors cx and cy are complex.
    !
    !  arguments
    !  == = ======
    !
    !  n       (input) integer
    !          the number of elements in the vectors cx and cy.
    !
    !  cx      (input/output) complex*16 array, dimension (n)
    !          on input, the vector x.
    !          on output, cx is overwritten with c*x + s*y.
    !
    !  incx    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  cy      (input/output) complex*16 array, dimension (n)
    !          on input, the vector y.
    !          on output, cy is overwritten with -conjg(s)*x + c*y.
    !
    !  incy    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  c       (input) double precision
    !  s       (input) complex*16
    !          c and s define a rotation
    !             [  c          s  ]
    !             [ -conjg(s)   c  ]
    !          where c*c + s*conjg(s) = 1.0.
    !
    ! == = ==================================================================
    !
    !     .. local scalars ..
    integer(psb_mpik_) :: i, ix, iy
    complex(psb_spk_)         stemp
    !     ..
    !     .. intrinsic functions ..
    !     ..
    !     .. executable statements ..
    !
    if( n <= 0 ) return
    if( incx == 1 .and. incy == 1 ) then 
      !
      !     code for both increments equal to 1
      !
      do  i = 1, n
        stemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - conjg(s)*cx(i)
        cx(i) = stemp
      end do
    else
      !
      !     code for unequal increments or equal increments not equal to 1
      !
      ix = 1
      iy = 1
      if( incx < 0 )ix = ( -n+1 )*incx + 1
      if( incy < 0 )iy = ( -n+1 )*incy + 1
      do  i = 1, n
        stemp  = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - conjg(s)*cx(ix)
        cx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    return
  end subroutine crot
  !
  !
  subroutine crotg(ca,cb,c,s)
    complex(psb_spk_) ca,cb,s
    real(psb_spk_) c
    real(psb_spk_) norm,scale
    complex(psb_spk_) alpha
    !
    if (cabs(ca) == 0.0) then 
      !
      c = 0.0d0
      s = (1.0,0.0)
      ca = cb
      return
    end if
    !

    scale = cabs(ca) + cabs(cb)
    norm = scale*sqrt((cabs(ca/cmplx(scale,0.0)))**2 +&
         &   (cabs(cb/cmplx(scale,0.0)))**2)
    alpha = ca /cabs(ca)
    c = cabs(ca) / norm
    s = alpha * conjg(cb) / norm
    ca = alpha * norm
    !

    return
  end subroutine crotg



  subroutine zrot( n, cx, incx, cy, incy, c, s )
    !
    !  -- lapack auxiliary routine (version 3.0) --
    !     univ. of tennessee, univ. of california berkeley, nag ltd.,
    !     courant institute, argonne national lab, and rice university
    !     october 31, 1992
    !
    !     .. scalar arguments ..
    integer(psb_mpik_) :: incx, incy, n
    real(psb_dpk_)    c
    complex(psb_dpk_)   s
    !     ..
    !     .. array arguments ..
    complex(psb_dpk_) cx( * ), cy( * )
    !     ..
    !
    !  purpose
    !  == = ====
    !
    !  zrot   applies a plane rotation, where the cos (c) is real and the
    !  sin (s) is complex, and the vectors cx and cy are complex.
    !
    !  arguments
    !  == = ======
    !
    !  n       (input) integer
    !          the number of elements in the vectors cx and cy.
    !
    !  cx      (input/output) complex*16 array, dimension (n)
    !          on input, the vector x.
    !          on output, cx is overwritten with c*x + s*y.
    !
    !  incx    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  cy      (input/output) complex*16 array, dimension (n)
    !          on input, the vector y.
    !          on output, cy is overwritten with -conjg(s)*x + c*y.
    !
    !  incy    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  c       (input) double precision
    !  s       (input) complex*16
    !          c and s define a rotation
    !             [  c          s  ]
    !             [ -conjg(s)   c  ]
    !          where c*c + s*conjg(s) = 1.0.
    !
    ! == = ==================================================================
    !
    !     .. local scalars ..
    integer(psb_mpik_) :: i, ix, iy
    complex(psb_dpk_)         stemp
    !     ..
    !     .. intrinsic functions ..
    intrinsic          dconjg
    !     ..
    !     .. executable statements ..
    !
    if( n <= 0 ) return
    if( incx == 1 .and. incy == 1 ) then 
      !
      !     code for both increments equal to 1
      !
      do  i = 1, n
        stemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - dconjg(s)*cx(i)
        cx(i) = stemp
      end do
    else
      !
      !     code for unequal increments or equal increments not equal to 1
      !
      ix = 1
      iy = 1
      if( incx < 0 )ix = ( -n+1 )*incx + 1
      if( incy < 0 )iy = ( -n+1 )*incy + 1
      do  i = 1, n
        stemp  = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - dconjg(s)*cx(ix)
        cx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    return
    return
  end subroutine zrot
  !
  !
  subroutine zrotg(ca,cb,c,s)
    complex(psb_dpk_) ca,cb,s
    real(psb_dpk_) c
    real(psb_dpk_) norm,scale
    complex(psb_dpk_) alpha
    !
    if (cdabs(ca) == 0.0d0) then 
      !
      c = 0.0d0
      s = (1.0d0,0.0d0)
      ca = cb
      return
    end if
    !

    scale = cdabs(ca) + cdabs(cb)
    norm = scale*dsqrt((cdabs(ca/cmplx(scale,0.0d0,kind=psb_dpk_)))**2 +&
         &   (cdabs(cb/cmplx(scale,0.0d0,kind=psb_dpk_)))**2)
    alpha = ca /cdabs(ca)
    c = cdabs(ca) / norm
    s = alpha * conjg(cb) / norm
    ca = alpha * norm
    !

    return
  end subroutine zrotg


  subroutine psb_scsprt(iout,a,iv,head,ivr,ivc)
    use psb_mat_mod
    integer(psb_ipk_), intent(in)       :: iout
    type(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_scsprt

  subroutine psb_scsprtn(fname,a,iv,head,ivr,ivc)
    use psb_mat_mod
    character(len=*), intent(in)  :: fname   
    type(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_scsprtn

  subroutine psb_dcsprt(iout,a,iv,head,ivr,ivc)
    use psb_mat_mod
    integer(psb_ipk_), intent(in)       :: iout
    type(psb_dspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_dcsprt

  subroutine psb_dcsprtn(fname,a,iv,head,ivr,ivc)
    use psb_mat_mod
    character(len=*), intent(in)  :: fname   
    type(psb_dspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_dcsprtn

  subroutine psb_ccsprt(iout,a,iv,head,ivr,ivc)
    use psb_mat_mod
    integer(psb_ipk_), intent(in)       :: iout
    type(psb_cspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_ccsprt

  subroutine psb_ccsprtn(fname,a,iv,head,ivr,ivc)
    use psb_mat_mod
    character(len=*), intent(in)  :: fname   
    type(psb_cspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_ccsprtn

  subroutine psb_zcsprt(iout,a,iv,head,ivr,ivc)
    use psb_mat_mod
    integer(psb_ipk_), intent(in)       :: iout
    type(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_zcsprt

  subroutine psb_zcsprtn(fname,a,iv,head,ivr,ivc)
    use psb_mat_mod
    character(len=*), intent(in)  :: fname   
    type(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_zcsprtn

end module psb_serial_mod

