!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006, 2010, 2015, 2017
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
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
module psb_renum_mod
  use psb_base_mod

  integer(psb_ipk_), parameter :: psb_mat_renum_identity_ = 0
  integer(psb_ipk_), parameter :: psb_mat_renum_gps_ = 456
  integer(psb_ipk_), parameter :: psb_mat_renum_amd_ = psb_mat_renum_gps_ + 1 
  
  
  interface psb_mat_renum
    subroutine psb_d_mat_renums(alg,mat,info,perm)
      import :: psb_ipk_, psb_dspmat_type
      character(len=*), intent(in) :: alg
      type(psb_dspmat_type), intent(inout) :: mat
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_d_mat_renums
    subroutine psb_d_mat_renum(alg,mat,info,perm)
      import :: psb_ipk_, psb_dspmat_type
      integer(psb_ipk_), intent(in) :: alg
      type(psb_dspmat_type), intent(inout) :: mat
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_d_mat_renum
    subroutine psb_s_mat_renums(alg,mat,info,perm)
      import :: psb_ipk_, psb_sspmat_type
      character(len=*), intent(in) :: alg
      type(psb_sspmat_type), intent(inout) :: mat
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_s_mat_renums
    subroutine psb_s_mat_renum(alg,mat,info,perm)
      import :: psb_ipk_, psb_sspmat_type
      integer(psb_ipk_), intent(in) :: alg
      type(psb_sspmat_type), intent(inout) :: mat
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_s_mat_renum
    subroutine psb_z_mat_renums(alg,mat,info,perm)
      import :: psb_ipk_, psb_zspmat_type
      character(len=*), intent(in) :: alg
      type(psb_zspmat_type), intent(inout) :: mat
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_z_mat_renums
    subroutine psb_z_mat_renum(alg,mat,info,perm)
      import :: psb_ipk_, psb_zspmat_type
      integer(psb_ipk_), intent(in) :: alg
      type(psb_zspmat_type), intent(inout) :: mat
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_z_mat_renum
    subroutine psb_c_mat_renums(alg,mat,info,perm)
      import :: psb_ipk_, psb_cspmat_type
      character(len=*), intent(in) :: alg
      type(psb_cspmat_type), intent(inout) :: mat
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_c_mat_renums
    subroutine psb_c_mat_renum(alg,mat,info,perm)
      import :: psb_ipk_, psb_cspmat_type
      integer(psb_ipk_), intent(in) :: alg
      type(psb_cspmat_type), intent(inout) :: mat
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_c_mat_renum
  end interface psb_mat_renum

  
  interface psb_cmp_bwpf
    subroutine psb_s_cmp_bwpf(mat,bwl,bwu,prf,info)
      import :: psb_ipk_, psb_sspmat_type
      type(psb_sspmat_type), intent(in) :: mat
      integer(psb_ipk_), intent(out) :: bwl, bwu
      integer(psb_ipk_), intent(out) :: prf
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_cmp_bwpf
    subroutine psb_d_cmp_bwpf(mat,bwl,bwu,prf,info)
      import :: psb_ipk_, psb_dspmat_type
      type(psb_dspmat_type), intent(in) :: mat
      integer(psb_ipk_), intent(out) :: bwl, bwu
      integer(psb_ipk_), intent(out) :: prf
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_d_cmp_bwpf
    subroutine psb_c_cmp_bwpf(mat,bwl,bwu,prf,info)
      import :: psb_ipk_, psb_cspmat_type
      type(psb_cspmat_type), intent(in) :: mat
      integer(psb_ipk_), intent(out) :: bwl, bwu
      integer(psb_ipk_), intent(out) :: prf
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_c_cmp_bwpf
    subroutine psb_z_cmp_bwpf(mat,bwl,bwu,prf,info)
      import :: psb_ipk_, psb_zspmat_type
      type(psb_zspmat_type), intent(in) :: mat
      integer(psb_ipk_), intent(out) :: bwl, bwu
      integer(psb_ipk_), intent(out) :: prf
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_z_cmp_bwpf
  end interface psb_cmp_bwpf


end module psb_renum_mod
