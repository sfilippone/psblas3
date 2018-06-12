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
module psb_z_serial_mod
  use psb_const_mod
  use psb_error_mod

  interface psb_amax
    function psb_zamax_s(n, x) result(val)
      import :: psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      complex(psb_dpk_), intent(in)  :: x(:)
      real(psb_dpk_)                 :: val
    end function psb_zamax_s
  end interface psb_amax

  interface psb_asum
    function psb_zasum_s(n, x) result(val)
      import :: psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in)  :: n
      complex(psb_dpk_), intent(in)  :: x(:)
      real(psb_dpk_)                 :: val
    end function psb_zasum_s
  end interface psb_asum

  interface psb_spspmm
    subroutine psb_zspspmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_zspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_zspmat_type), intent(in)  :: a,b
      type(psb_zspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_zspspmm
    subroutine psb_zcsrspspmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_z_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_z_csr_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(out) :: c
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_zcsrspspmm
    subroutine psb_zcscspspmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_z_csc_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_z_csc_sparse_mat), intent(in) :: a,b
      type(psb_z_csc_sparse_mat), intent(out) :: c
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_zcscspspmm
  end interface

  interface psb_symbmm
    subroutine psb_zsymbmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_zspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_zspmat_type), intent(in)  :: a,b
      type(psb_zspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_zsymbmm
    subroutine psb_zbase_symbmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_z_base_sparse_mat, psb_z_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_z_base_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_zbase_symbmm
  end interface psb_symbmm

  interface psb_numbmm
    subroutine psb_znumbmm(a,b,c)
      use psb_z_mat_mod, only : psb_zspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_zspmat_type), intent(in) :: a,b
      type(psb_zspmat_type), intent(inout)  :: c
    end subroutine psb_znumbmm
    subroutine psb_zbase_numbmm(a,b,c)
      use psb_z_mat_mod, only : psb_z_base_sparse_mat, psb_z_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_z_base_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_zbase_numbmm
  end interface psb_numbmm

  interface psb_rwextd
    subroutine psb_zrwextd(nr,a,info,b,rowscale)
      use psb_z_mat_mod, only : psb_zspmat_type
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_), intent(in)                          :: nr
      type(psb_zspmat_type), intent(inout)        :: a
      integer(psb_ipk_),intent(out)                          :: info
      type(psb_zspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_zrwextd
    subroutine psb_zbase_rwextd(nr,a,info,b,rowscale)
      use psb_z_mat_mod, only : psb_z_base_sparse_mat
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
    subroutine psb_zgeprtn2(fname,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      complex(psb_dpk_), intent(in)    :: a(:,:)
      character(len=*), optional    :: head
    end subroutine psb_zgeprtn2
    subroutine psb_zgeprtn1(fname,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      character(len=*), intent(in)  :: fname   
      complex(psb_dpk_), intent(in)    :: a(:)
      character(len=*), optional    :: head
    end subroutine psb_zgeprtn1
    subroutine psb_zgeprt2(iout,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      complex(psb_dpk_), intent(in)     :: a(:,:)
      character(len=*), optional     :: head
    end subroutine psb_zgeprt2
    subroutine psb_zgeprt1(iout,a,head)
      use psb_const_mod, only : psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)            :: iout
      complex(psb_dpk_), intent(in)     :: a(:)
      character(len=*), optional     :: head
    end subroutine psb_zgeprt1
  end interface psb_geprt

  interface psb_csprt
    module procedure psb_zcsprt, psb_zcsprtn
  end interface psb_csprt

  interface psb_spdot_srtd
    function psb_z_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1,nv2
      integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
      complex(psb_dpk_), intent(in) :: v1(*),v2(*)
      complex(psb_dpk_)      :: dot
    end function psb_z_spdot_srtd
  end interface psb_spdot_srtd


  interface psb_spge_dot
    function psb_z_spge_dot(nv1,iv1,v1,v2) result(dot) 
      use psb_const_mod, only : psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in) :: nv1
      integer(psb_ipk_), intent(in) :: iv1(*)
      complex(psb_dpk_), intent(in) :: v1(*),v2(*)
      complex(psb_dpk_)      :: dot
    end function psb_z_spge_dot
  end interface psb_spge_dot


  interface psb_aspxpby
    subroutine psb_z_aspxpby(alpha, nx, ix, x, beta, y, info)
      use psb_const_mod, only : psb_ipk_, psb_dpk_
      integer(psb_ipk_), intent(in)               :: nx
      integer(psb_ipk_), intent(in)               :: ix(:)
      complex(psb_dpk_), intent (in)       :: x(:)
      complex(psb_dpk_), intent (inout)    :: y(:)
      complex(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_z_aspxpby
  end interface psb_aspxpby

  interface psb_spspmm
    subroutine psb_lzspspmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_lzspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_lzspmat_type), intent(in)  :: a,b
      type(psb_lzspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_lzspspmm
    subroutine psb_lzcsrspspmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_lz_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_lz_csr_sparse_mat), intent(in) :: a,b
      type(psb_lz_csr_sparse_mat), intent(out) :: c
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psb_lzcsrspspmm
!!$    subroutine psb_zcscspspmm(a,b,c,info)
!!$      use psb_z_mat_mod, only : psb_z_csc_sparse_mat
!!$      import :: psb_ipk_
!!$      implicit none 
!!$      class(psb_z_csc_sparse_mat), intent(in) :: a,b
!!$      type(psb_z_csc_sparse_mat), intent(out) :: c
!!$      integer(psb_ipk_), intent(out)          :: info
!!$    end subroutine psb_zcscspspmm
  end interface psb_spspmm

  interface psb_symbmm
    subroutine psb_lzsymbmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_lzspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_lzspmat_type), intent(in)  :: a,b
      type(psb_lzspmat_type), intent(out) :: c
      integer(psb_ipk_), intent(out)                :: info
    end subroutine psb_lzsymbmm
    subroutine psb_lzbase_symbmm(a,b,c,info)
      use psb_z_mat_mod, only : psb_lz_base_sparse_mat, psb_lz_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_lz_base_sparse_mat), intent(in) :: a,b
      type(psb_lz_csr_sparse_mat), intent(out)  :: c
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_lzbase_symbmm
  end interface psb_symbmm

  interface psb_numbmm
    subroutine psb_lznumbmm(a,b,c)
      use psb_z_mat_mod, only : psb_lzspmat_type
      import :: psb_ipk_
      implicit none 
      type(psb_lzspmat_type), intent(in) :: a,b
      type(psb_lzspmat_type), intent(inout)  :: c
    end subroutine psb_lznumbmm
    subroutine psb_lzbase_numbmm(a,b,c)
      use psb_z_mat_mod, only : psb_lz_base_sparse_mat, psb_lz_csr_sparse_mat
      import :: psb_ipk_
      implicit none 
      class(psb_lz_base_sparse_mat), intent(in) :: a,b
      type(psb_lz_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_lzbase_numbmm
  end interface psb_numbmm
  
contains

  subroutine psb_zcsprt(iout,a,iv,head,ivr,ivc)
    use psb_z_mat_mod, only : psb_zspmat_type
    integer(psb_ipk_), intent(in)       :: iout
    type(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(iout,iv,head,ivr,ivc)

  end subroutine psb_zcsprt

  subroutine psb_zcsprtn(fname,a,iv,head,ivr,ivc)
    use psb_z_mat_mod, only : psb_zspmat_type
    character(len=*), intent(in)  :: fname   
    type(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_), intent(in), optional :: iv(:)
    character(len=*), optional    :: head
    integer(psb_ipk_), intent(in), optional :: ivr(:),ivc(:)

    call a%print(fname,iv,head,ivr,ivc)

  end subroutine psb_zcsprtn

end module psb_z_serial_mod

