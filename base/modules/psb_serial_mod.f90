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


  interface psb_symbmm
    subroutine psb_ssymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_sspmat_type
      implicit none 
      type(psb_sspmat_type), intent(in)  :: a,b
      type(psb_sspmat_type), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_ssymbmm
    subroutine psb_sbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_s_base_sparse_mat, psb_s_csr_sparse_mat
      implicit none 
      class(psb_s_base_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_sbase_symbmm
    subroutine psb_dsymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_dspmat_type
      implicit none 
      type(psb_dspmat_type), intent(in)  :: a,b
      type(psb_dspmat_type), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_dsymbmm
    subroutine psb_dbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_d_base_sparse_mat, psb_d_csr_sparse_mat
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_dbase_symbmm
    subroutine psb_csymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_cspmat_type
      implicit none 
      type(psb_cspmat_type), intent(in)  :: a,b
      type(psb_cspmat_type), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_csymbmm
    subroutine psb_cbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_c_base_sparse_mat, psb_c_csr_sparse_mat
      implicit none 
      class(psb_c_base_sparse_mat), intent(in) :: a,b
      type(psb_c_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_cbase_symbmm
    subroutine psb_zsymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_zspmat_type
      implicit none 
      type(psb_zspmat_type), intent(in)  :: a,b
      type(psb_zspmat_type), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_zsymbmm
    subroutine psb_zbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_z_base_sparse_mat, psb_z_csr_sparse_mat
      implicit none 
      class(psb_z_base_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_zbase_symbmm
  end interface psb_symbmm

  interface psb_numbmm
    subroutine psb_snumbmm(a,b,c)
      use psb_mat_mod, only : psb_sspmat_type
      implicit none 
      type(psb_sspmat_type), intent(in) :: a,b
      type(psb_sspmat_type), intent(inout)  :: c
    end subroutine psb_snumbmm
    subroutine psb_sbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_s_base_sparse_mat, psb_s_csr_sparse_mat
      implicit none 
      class(psb_s_base_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_sbase_numbmm
    subroutine psb_dnumbmm(a,b,c)
      use psb_mat_mod, only : psb_dspmat_type
      implicit none 
      type(psb_dspmat_type), intent(in) :: a,b
      type(psb_dspmat_type), intent(inout)  :: c
    end subroutine psb_dnumbmm
    subroutine psb_dbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_d_base_sparse_mat, psb_d_csr_sparse_mat
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_dbase_numbmm
    subroutine psb_cnumbmm(a,b,c)
      use psb_mat_mod, only : psb_cspmat_type
      implicit none 
      type(psb_cspmat_type), intent(in) :: a,b
      type(psb_cspmat_type), intent(inout)  :: c
    end subroutine psb_cnumbmm
    subroutine psb_cbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_c_base_sparse_mat, psb_c_csr_sparse_mat
      implicit none 
      class(psb_c_base_sparse_mat), intent(in) :: a,b
      type(psb_c_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_cbase_numbmm
    subroutine psb_znumbmm(a,b,c)
      use psb_mat_mod, only : psb_zspmat_type
      implicit none 
      type(psb_zspmat_type), intent(in) :: a,b
      type(psb_zspmat_type), intent(inout)  :: c
    end subroutine psb_znumbmm
    subroutine psb_zbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_z_base_sparse_mat, psb_z_csr_sparse_mat
      implicit none 
      class(psb_z_base_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_zbase_numbmm
  end interface psb_numbmm

  interface psb_rwextd
    subroutine psb_srwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_sspmat_type
      implicit none
      integer, intent(in)                          :: nr
      type(psb_sspmat_type), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_sspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_srwextd
    subroutine psb_sbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_s_base_sparse_mat
      implicit none
      integer, intent(in)                                :: nr
      class(psb_s_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_s_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_sbase_rwextd
    subroutine psb_drwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_dspmat_type
      implicit none
      integer, intent(in)                          :: nr
      type(psb_dspmat_type), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_dspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_drwextd
    subroutine psb_dbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_d_base_sparse_mat
      implicit none
      integer, intent(in)                                :: nr
      class(psb_d_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_dbase_rwextd
    subroutine psb_crwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_cspmat_type
      implicit none
      integer, intent(in)                          :: nr
      type(psb_cspmat_type), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_cspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_crwextd
    subroutine psb_cbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_c_base_sparse_mat
      implicit none
      integer, intent(in)                                :: nr
      class(psb_c_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_cbase_rwextd
    subroutine psb_zrwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_zspmat_type
      implicit none
      integer, intent(in)                          :: nr
      type(psb_zspmat_type), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_zspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_zrwextd
    subroutine psb_zbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_z_base_sparse_mat
      implicit none
      integer, intent(in)                                :: nr
      class(psb_z_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_zbase_rwextd
  end interface psb_rwextd


  interface psb_geprt
    subroutine psb_sgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_spk_
      character(len=*), intent(in)  :: fname   
      real(psb_spk_), intent(in)    :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_sgeprtn2
    subroutine psb_sgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_spk_
      character(len=*), intent(in)  :: fname   
      real(psb_spk_), intent(in)    :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_sgeprtn1
    subroutine psb_sgeprt2(iout,a,head)
      use psb_const_mod, only : psb_spk_
      integer, intent(in)            :: iout
      real(psb_spk_), intent(in)     :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_sgeprt2
    subroutine psb_sgeprt1(iout,a,head)
      use psb_const_mod, only : psb_spk_
      integer, intent(in)            :: iout
      real(psb_spk_), intent(in)     :: a(:)
      character(len=*), optional     :: head
    end subroutine psb_sgeprt1
    subroutine psb_dgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_dpk_
      character(len=*), intent(in)  :: fname   
      real(psb_dpk_), intent(in)    :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_dgeprtn2
    subroutine psb_dgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_dpk_
      character(len=*), intent(in)  :: fname   
      real(psb_dpk_), intent(in)    :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_dgeprtn1
    subroutine psb_dgeprt2(iout,a,head)
      use psb_const_mod, only : psb_dpk_
      integer, intent(in)            :: iout
      real(psb_dpk_), intent(in)     :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_dgeprt2
    subroutine psb_dgeprt1(iout,a,head)
      use psb_const_mod, only : psb_dpk_
      integer, intent(in)            :: iout
      real(psb_dpk_), intent(in)     :: a(:)
      character(len=*), optional     :: head
    end subroutine psb_dgeprt1
    subroutine psb_cgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_spk_
      character(len=*), intent(in)  :: fname   
      complex(psb_spk_), intent(in) :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_cgeprtn2
    subroutine psb_cgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_spk_
      character(len=*), intent(in)  :: fname   
      complex(psb_spk_), intent(in) :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_cgeprtn1
    subroutine psb_cgeprt2(iout,a,head)
      use psb_const_mod, only : psb_spk_
      integer, intent(in)            :: iout
      complex(psb_spk_), intent(in)  :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_cgeprt2
    subroutine psb_cgeprt1(iout,a,head)
      use psb_const_mod, only : psb_spk_
      integer, intent(in)            :: iout
      complex(psb_spk_), intent(in)  :: a(:)
      character(len=*), optional     :: head
    end subroutine psb_cgeprt1
    subroutine psb_zgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_dpk_
      character(len=*), intent(in)  :: fname   
      complex(psb_dpk_), intent(in) :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_zgeprtn2
    subroutine psb_zgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_dpk_
      character(len=*), intent(in)  :: fname   
      complex(psb_dpk_), intent(in) :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_zgeprtn1
    subroutine psb_zgeprt2(iout,a,head)
      use psb_const_mod, only : psb_dpk_
      integer, intent(in)            :: iout
      complex(psb_dpk_), intent(in)  :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_zgeprt2
    subroutine psb_zgeprt1(iout,a,head)
      use psb_const_mod, only : psb_dpk_
      integer, intent(in)            :: iout
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
      use psb_const_mod
      integer, intent(in) :: nv1,nv2
      integer, intent(in) :: iv1(*), iv2(*)
      real(psb_spk_), intent(in) :: v1(*),v2(*)
      real(psb_spk_)      :: dot
    end function psb_s_spdot_srtd

    function psb_d_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod
      integer, intent(in) :: nv1,nv2
      integer, intent(in) :: iv1(*), iv2(*)
      real(psb_dpk_), intent(in) :: v1(*), v2(*)
      real(psb_dpk_)      :: dot
    end function psb_d_spdot_srtd

    function psb_c_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod
      integer, intent(in) :: nv1,nv2
      integer, intent(in) :: iv1(*), iv2(*)
      complex(psb_spk_), intent(in) :: v1(*),v2(*)
      complex(psb_spk_)      :: dot
    end function psb_c_spdot_srtd

    function psb_z_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod
      integer, intent(in) :: nv1,nv2
      integer, intent(in) :: iv1(*), iv2(*)
      complex(psb_dpk_), intent(in) :: v1(*),v2(*)
      complex(psb_dpk_)      :: dot
    end function psb_z_spdot_srtd
  end interface psb_spdot_srtd


  interface psb_spge_dot
    function psb_s_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod
      integer, intent(in) :: nv1
      integer, intent(in) :: iv1(*)
      real(psb_spk_), intent(in) :: v1(*),v2(*)
      real(psb_spk_)      :: dot
    end function psb_s_spge_dot

    function psb_d_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod
      integer, intent(in) :: nv1
      integer, intent(in) :: iv1(*)
      real(psb_dpk_), intent(in) :: v1(*),v2(*)
      real(psb_dpk_)      :: dot
    end function psb_d_spge_dot

    function psb_c_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod
      integer, intent(in) :: nv1
      integer, intent(in) :: iv1(*)
      complex(psb_spk_), intent(in) :: v1(*),v2(*)
      complex(psb_spk_)      :: dot
    end function psb_c_spge_dot

    function psb_z_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod
      integer, intent(in) :: nv1
      integer, intent(in) :: iv1(*)
      complex(psb_dpk_), intent(in) :: v1(*),v2(*)
      complex(psb_dpk_)      :: dot
    end function psb_z_spge_dot
  end interface psb_spge_dot


  interface psb_nspaxpby
    subroutine psb_d_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
      use psb_const_mod
      integer, intent(out)              :: nz
      integer, intent(out)              :: iz(:)
      real(psb_dpk_), intent (out)      :: z(:)
      integer, intent(in)               :: nx, ny
      integer, intent(in)               :: ix(:), iy(:)
      real(psb_dpk_), intent (in)       :: x(:), y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer, intent(out)              :: info
    end subroutine psb_d_nspaxpby
  end interface psb_nspaxpby

  interface psb_aspxpby
    subroutine psb_s_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod
      integer, intent(in)               :: nx
      integer, intent(in)               :: ix(:)
      real(psb_spk_), intent (in)       :: x(:)
      real(psb_spk_), intent (inout)    :: y(:)
      real(psb_spk_), intent (in)       :: alpha, beta
      integer, intent(out)              :: info
    end subroutine psb_s_aspxpby
    subroutine psb_d_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod
      integer, intent(in)               :: nx
      integer, intent(in)               :: ix(:)
      real(psb_dpk_), intent (in)       :: x(:)
      real(psb_dpk_), intent (inout)    :: y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer, intent(out)              :: info
    end subroutine psb_d_aspxpby
    subroutine psb_c_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod
      integer, intent(in)               :: nx
      integer, intent(in)               :: ix(:)
      complex(psb_spk_), intent (in)    :: x(:)
      complex(psb_spk_), intent (inout) :: y(:)
      complex(psb_spk_), intent (in)    :: alpha, beta
      integer, intent(out)              :: info
    end subroutine psb_c_aspxpby
    subroutine psb_z_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod
      integer, intent(in)               :: nx
      integer, intent(in)               :: ix(:)
      complex(psb_dpk_), intent (in)    :: x(:)
      complex(psb_dpk_), intent (inout) :: y(:)
      complex(psb_dpk_), intent (in)    :: alpha, beta
      integer, intent(out)              :: info
    end subroutine psb_z_aspxpby
  end interface psb_aspxpby

contains

  subroutine psb_scsprt(iout,a,iv,head,ivr,ivc)
    use psb_mat_mod
    integer, intent(in)       :: iout
    type(psb_sspmat_type), intent(in) :: a
    integer, intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer, intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_scsprt

  subroutine psb_scsprtn(fname,a,iv,head,ivr,ivc)
    use psb_mat_mod
    character(len=*), intent(in)  :: fname   
    type(psb_sspmat_type), intent(in) :: a
    integer, intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer, intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_scsprtn

  subroutine psb_dcsprt(iout,a,iv,head,ivr,ivc)
    use psb_mat_mod
    integer, intent(in)       :: iout
    type(psb_dspmat_type), intent(in) :: a
    integer, intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer, intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_dcsprt

  subroutine psb_dcsprtn(fname,a,iv,head,ivr,ivc)
    use psb_mat_mod
    character(len=*), intent(in)  :: fname   
    type(psb_dspmat_type), intent(in) :: a
    integer, intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer, intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_dcsprtn

  subroutine psb_ccsprt(iout,a,iv,head,ivr,ivc)
    use psb_mat_mod
    integer, intent(in)       :: iout
    type(psb_cspmat_type), intent(in) :: a
    integer, intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer, intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_ccsprt

  subroutine psb_ccsprtn(fname,a,iv,head,ivr,ivc)
    use psb_mat_mod
    character(len=*), intent(in)  :: fname   
    type(psb_cspmat_type), intent(in) :: a
    integer, intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer, intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_ccsprtn

  subroutine psb_zcsprt(iout,a,iv,head,ivr,ivc)
    use psb_mat_mod
    integer, intent(in)       :: iout
    type(psb_zspmat_type), intent(in) :: a
    integer, intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer, intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_zcsprt

  subroutine psb_zcsprtn(fname,a,iv,head,ivr,ivc)
    use psb_mat_mod
    character(len=*), intent(in)  :: fname   
    type(psb_zspmat_type), intent(in) :: a
    integer, intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer, intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_zcsprtn

end module psb_serial_mod

