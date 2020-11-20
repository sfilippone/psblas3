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
!    Moved here from AMG4PSBLAS, original copyright below.
!
!
!                       AMG-AINV: Approximate Inverse plugin for
!                             AMG4PSBLAS version 1.0
!
!    (C) Copyright 2020
!
!                        Salvatore Filippone  University of Rome Tor Vergata
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
! File: psb_z_invk_fact_mod.f90
!
! Module: psb_z_invk_fact_mod
!
!  This module defines some interfaces used internally by the implementation of
!  psb_z_invk_solver, but not visible to the end user.
!
!
module psb_z_invk_fact_mod

  use psb_base_mod
  use psb_prec_const_mod

  interface psb_invk_fact
    subroutine psb_z_invk_bld(a,fill1, fill2,lmat,d,umat,desc,info,blck)
      ! import
      import psb_zspmat_type, psb_ipk_, psb_dpk_, psb_desc_type
      ! Arguments
      type(psb_zspmat_type), intent(in), target   :: a
      integer(psb_ipk_), intent(in)               :: fill1, fill2
      type(psb_zspmat_type), intent(inout)        :: lmat, umat
      complex(psb_dpk_), allocatable                 :: d(:)
      Type(psb_desc_type), Intent(inout)          :: desc
      integer(psb_ipk_), intent(out)              :: info
      type(psb_zspmat_type), intent(in), optional :: blck
    end subroutine
  end interface

  ! Inner workings
  interface psb_sparse_invk
    subroutine psb_zsparse_invk(n,a,z,fill_in,info,inlevs)
      ! Import
      import psb_ipk_, psb_zspmat_type
      ! Arguments
      integer(psb_ipk_), intent(in)           :: n
      type(psb_zspmat_type), intent(in)       :: a
      type(psb_zspmat_type), intent(inout)    :: z
      integer(psb_ipk_), intent(in)           :: fill_in
      integer(psb_ipk_), intent(out)          :: info
      integer(psb_ipk_), intent(in), optional :: inlevs(:)
    end subroutine psb_zsparse_invk
  end interface

  interface psb_invk_inv
    subroutine psb_zinvk_inv(fill_in,i,row,rowlevs,heap,uia1,uia2,uaspk,uplevs,&
         & nidx,idxs,info)

      use psb_base_mod, only : psb_zspmat_type, psb_dpk_, psb_i_heap, psb_ipk_
      implicit none

      ! Arguments
      type(psb_i_heap), intent(inout)               :: heap
      integer(psb_ipk_), intent(in)                 :: i, fill_in
      integer(psb_ipk_), intent(inout)              :: nidx,info
      integer(psb_ipk_), intent(inout)              :: rowlevs(:)
      integer(psb_ipk_), allocatable, intent(inout) :: idxs(:)
      integer(psb_ipk_), intent(in)                 :: uia1(:),uia2(:),uplevs(:)
      complex(psb_dpk_), intent(in)                    :: uaspk(:)
      complex(psb_dpk_), intent(inout)                 :: row(:)


    end subroutine psb_zinvk_inv
  end interface

  interface  psb_invk_copyin
    subroutine psb_z_invk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,&
         & ktrw,trw,info,sign,inlevs)
      ! Import
      use psb_base_mod, only : psb_z_csr_sparse_mat, psb_z_coo_sparse_mat,&
           & psb_dpk_, psb_i_heap, psb_ipk_
      implicit none
      ! Arguments
      type(psb_z_csr_sparse_mat), intent(in)    :: a
      type(psb_z_coo_sparse_mat), intent(inout) :: trw
      integer(psb_ipk_), intent(in)             :: i,m,jmin,jmax
      integer(psb_ipk_), intent(inout)          :: ktrw,info
      integer(psb_ipk_), intent(inout)          :: rowlevs(:)
      complex(psb_dpk_), intent(inout)             :: row(:)
      type(psb_i_heap), intent(inout)           :: heap
      real(psb_dpk_), optional, intent(in)      :: sign
      integer(psb_ipk_), intent(in), optional   :: inlevs(:)
    end subroutine psb_z_invk_copyin
  end interface

  interface psb_invk_copyout
    subroutine psb_z_invk_copyout(fill_in,i,m,row,rowlevs,nidx,idxs,&
         &  l2,uia1,uia2,uaspk,info)
      use psb_base_mod, only : psb_zspmat_type, psb_dpk_, psb_i_heap, psb_ipk_
      implicit none
      ! Arguments
      integer(psb_ipk_), intent(in)                 :: fill_in, i, m, nidx
      integer(psb_ipk_), intent(inout)              :: l2, info
      integer(psb_ipk_), intent(inout)              :: rowlevs(:), idxs(:)
      integer(psb_ipk_), allocatable, intent(inout) :: uia1(:), uia2(:)
      complex(psb_dpk_), allocatable, intent(inout)    :: uaspk(:)
      complex(psb_dpk_), intent(inout)                 :: row(:)
    end subroutine psb_z_invk_copyout
  end interface

end module
