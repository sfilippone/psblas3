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
! File: psb_c_invt_fact_mod.f90
!
! Module: psb_c_invt_fact_mod
!
!  This module defines some interfaces used internally by the implementation of
!  psb_c_invt_solver, but not visible to the end user.
!
!
module psb_c_invt_fact_mod

  use psb_base_mod
  use psb_prec_const_mod

  interface psb_invt_fact
      subroutine psb_c_invt_bld(a,fillin,invfill,thresh,invthresh,&
         & lmat,d,umat,desc,info,blck)
         ! Import
         import psb_cspmat_type, psb_spk_, psb_ipk_, psb_desc_type
         ! Arguments
         type(psb_cspmat_type), intent(in), target   :: a
         integer(psb_ipk_), intent(in)               :: fillin,invfill
         real(psb_spk_), intent(in)                :: thresh
         real(psb_spk_), intent(in)                  :: invthresh
         type(psb_cspmat_type), intent(inout)        :: lmat, umat
         complex(psb_spk_), allocatable                 :: d(:)
         Type(psb_desc_type), Intent(inout)          :: desc
         integer(psb_ipk_), intent(out)              :: info
         type(psb_cspmat_type), intent(in), optional :: blck
      end subroutine psb_c_invt_bld
  end interface

  ! Interfaces for the inner workings

  interface
      subroutine psb_csparse_invt(n,a,z,nzrmax,sp_thresh,info)
         ! Import
         import psb_cspmat_type, psb_spk_, psb_ipk_
         ! Arguments
         integer(psb_ipk_), intent(in)        :: n
         type(psb_cspmat_type), intent(in)    :: a
         type(psb_cspmat_type), intent(inout) :: z
         integer(psb_ipk_), intent(in)        :: nzrmax
         real(psb_spk_), intent(in)           :: sp_thresh
         integer(psb_ipk_), intent(out)       :: info
      end subroutine psb_csparse_invt
  end interface

  interface
     subroutine psb_c_invt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
          & irwt,ktrw,trw,info,sign)
       ! Import
       import psb_c_csr_sparse_mat, psb_c_coo_sparse_mat, psb_spk_, &
         & psb_ipk_, psb_i_heap
       ! Arguments
       type(psb_c_csr_sparse_mat), intent(in)    :: a
       type(psb_c_coo_sparse_mat), intent(inout) :: trw
       integer(psb_ipk_), intent(in)             :: i, m,jmin,jmax,jd
       integer(psb_ipk_), intent(inout)          :: ktrw,nlw,nup,jmaxup,info
       integer(psb_ipk_), intent(inout)          :: irwt(:)
       real(psb_spk_), intent(inout)          :: nrmi
       complex(psb_spk_), intent(inout)            :: row(:)
       type(psb_i_heap), intent(inout)           :: heap
       real(psb_spk_), intent(in), optional    :: sign
    end subroutine psb_c_invt_copyin
  end interface

  interface
     subroutine psb_c_invt_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
          & nidx,idxs,l2,ja,irp,val,info)
      ! Import
      import psb_spk_, psb_ipk_
      ! Arguments
      integer(psb_ipk_), intent(in)                 :: fill_in,i,m,nidx,nlw,nup,jmaxup
      integer(psb_ipk_), intent(in)                 :: idxs(:)
      integer(psb_ipk_), intent(inout)              :: l2, info
      integer(psb_ipk_), allocatable, intent(inout) :: ja(:),irp(:)
      real(psb_spk_), intent(in)                    :: thres,nrmi
      complex(psb_spk_),allocatable, intent(inout)     :: val(:)
      complex(psb_spk_), intent(inout)                 :: row(:)
     end subroutine psb_c_invt_copyout
  end interface

  interface psb_invt_inv
    subroutine psb_c_invt_inv(thres,i,nrmi,row,heap,irwt,ja,irp,val, &
      & nidx,idxs,info)
      ! import
      import psb_spk_, psb_i_heap, psb_ipk_
      ! Arguments
      type(psb_i_heap), intent(inout)     :: heap
      integer(psb_ipk_), intent(in)       :: i
      integer(psb_ipk_), intent(inout)    :: nidx,info
      integer(psb_ipk_), intent(inout)    :: irwt(:)
      real(psb_spk_), intent(in)          :: thres,nrmi
      integer(psb_ipk_), allocatable, intent(inout) :: idxs(:)
      integer(psb_ipk_), intent(in)       :: ja(:),irp(:)
      complex(psb_spk_), intent(in)          :: val(:)
      complex(psb_spk_), intent(inout)       :: row(:)
    end subroutine
  end interface

contains

end module psb_c_invt_fact_mod
