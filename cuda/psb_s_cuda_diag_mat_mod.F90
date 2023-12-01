!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
  

module psb_s_cuda_diag_mat_mod

  use iso_c_binding
  use psb_base_mod
  use psb_s_dia_mat_mod

  type, extends(psb_s_dia_sparse_mat) :: psb_s_cuda_diag_sparse_mat
    !
    ! ITPACK/HLL format, extended.
    ! We are adding here the routines to create a copy of the data
    ! into the GPU. 
    ! If HAVE_SPGPU is undefined this is just
    ! a copy of HLL, indistinguishable.
    ! 
    type(c_ptr) :: deviceMat = c_null_ptr

  contains
    procedure, nopass  :: get_fmt       => s_cuda_diag_get_fmt
    procedure, pass(a) :: sizeof        => s_cuda_diag_sizeof
    procedure, pass(a) :: vect_mv       => psb_s_cuda_diag_vect_mv
!    procedure, pass(a) :: csmm          => psb_s_cuda_diag_csmm
    procedure, pass(a) :: csmv          => psb_s_cuda_diag_csmv
!    procedure, pass(a) :: in_vect_sv    => psb_s_cuda_diag_inner_vect_sv
!    procedure, pass(a) :: scals         => psb_s_cuda_diag_scals
!    procedure, pass(a) :: scalv         => psb_s_cuda_diag_scal
!    procedure, pass(a) :: reallocate_nz => psb_s_cuda_diag_reallocate_nz
!    procedure, pass(a) :: allocate_mnnz => psb_s_cuda_diag_allocate_mnnz
    ! Note: we do *not* need the TO methods, because the parent type
    ! methods will work. 
    procedure, pass(a) :: cp_from_coo   => psb_s_cuda_cp_diag_from_coo
!    procedure, pass(a) :: cp_from_fmt   => psb_s_cuda_cp_diag_from_fmt
    procedure, pass(a) :: mv_from_coo   => psb_s_cuda_mv_diag_from_coo
!    procedure, pass(a) :: mv_from_fmt   => psb_s_cuda_mv_diag_from_fmt
    procedure, pass(a) :: free          => s_cuda_diag_free
    procedure, pass(a) :: mold          => psb_s_cuda_diag_mold
    procedure, pass(a) :: to_gpu        => psb_s_cuda_diag_to_gpu
    final              :: s_cuda_diag_finalize
  end type psb_s_cuda_diag_sparse_mat

  private :: s_cuda_diag_get_nzeros, s_cuda_diag_free,  s_cuda_diag_get_fmt, &
       & s_cuda_diag_get_size, s_cuda_diag_sizeof, s_cuda_diag_get_nz_row


  interface 
    subroutine psb_s_cuda_diag_vect_mv(alpha,a,x,beta,y,info,trans) 
      import :: psb_s_cuda_diag_sparse_mat, psb_spk_, psb_s_base_vect_type, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(in)    :: a
      real(psb_spk_), intent(in)                 :: alpha, beta
      class(psb_s_base_vect_type), intent(inout) :: x
      class(psb_s_base_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)             :: info
      character, optional, intent(in)            :: trans
    end subroutine psb_s_cuda_diag_vect_mv
  end interface

  interface 
    subroutine psb_s_cuda_diag_inner_vect_sv(alpha,a,x,beta,y,info,trans) 
      import :: psb_ipk_, psb_s_cuda_diag_sparse_mat, psb_spk_,  psb_s_base_vect_type
      class(psb_s_cuda_diag_sparse_mat), intent(in)    :: a
      real(psb_spk_), intent(in)                 :: alpha, beta
      class(psb_s_base_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(out)             :: info
      character, optional, intent(in)            :: trans
    end subroutine psb_s_cuda_diag_inner_vect_sv
  end interface

  interface
    subroutine  psb_s_cuda_diag_reallocate_nz(nz,a) 
      import :: psb_s_cuda_diag_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in)              :: nz
      class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a
    end subroutine psb_s_cuda_diag_reallocate_nz
  end interface

  interface
    subroutine  psb_s_cuda_diag_allocate_mnnz(m,n,a,nz) 
      import :: psb_s_cuda_diag_sparse_mat, psb_ipk_
      integer(psb_ipk_), intent(in)              :: m,n
      class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in), optional    :: nz
    end subroutine psb_s_cuda_diag_allocate_mnnz
  end interface

  interface 
    subroutine psb_s_cuda_diag_mold(a,b,info) 
      import :: psb_s_cuda_diag_sparse_mat, psb_s_base_sparse_mat, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(in)                  :: a
      class(psb_s_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_s_cuda_diag_mold
  end interface

  interface 
    subroutine psb_s_cuda_diag_to_gpu(a,info, nzrm) 
      import :: psb_s_cuda_diag_sparse_mat, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)             :: info
      integer(psb_ipk_), intent(in), optional    :: nzrm
    end subroutine psb_s_cuda_diag_to_gpu
  end interface

  interface 
    subroutine psb_s_cuda_cp_diag_from_coo(a,b,info) 
      import :: psb_s_cuda_diag_sparse_mat, psb_s_coo_sparse_mat, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_s_cuda_cp_diag_from_coo
  end interface
  
  interface 
    subroutine psb_s_cuda_cp_diag_from_fmt(a,b,info) 
      import :: psb_s_cuda_diag_sparse_mat, psb_s_base_sparse_mat, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_s_cuda_cp_diag_from_fmt
  end interface
  
  interface 
    subroutine psb_s_cuda_mv_diag_from_coo(a,b,info) 
      import :: psb_s_cuda_diag_sparse_mat, psb_s_coo_sparse_mat, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_s_cuda_mv_diag_from_coo
  end interface
  

  interface 
    subroutine psb_s_cuda_mv_diag_from_fmt(a,b,info) 
      import :: psb_s_cuda_diag_sparse_mat, psb_s_base_sparse_mat, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(inout)  :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_s_cuda_mv_diag_from_fmt
  end interface
  
  interface 
    subroutine psb_s_cuda_diag_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_s_cuda_diag_sparse_mat, psb_spk_, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)              :: alpha, beta, x(:)
      real(psb_spk_), intent(inout)           :: y(:)
      integer(psb_ipk_), intent(out)          :: info
      character, optional, intent(in)         :: trans
    end subroutine psb_s_cuda_diag_csmv
  end interface
  interface 
    subroutine psb_s_cuda_diag_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_s_cuda_diag_sparse_mat, psb_spk_, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)              :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout)           :: y(:,:)
      integer(psb_ipk_), intent(out)          :: info
      character, optional, intent(in)         :: trans
    end subroutine psb_s_cuda_diag_csmm
  end interface
  
  interface 
    subroutine psb_s_cuda_diag_scal(d,a,info, side) 
      import :: psb_s_cuda_diag_sparse_mat, psb_spk_, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)                 :: d(:)
      integer(psb_ipk_), intent(out)             :: info
      character, intent(in), optional            :: side
    end subroutine psb_s_cuda_diag_scal
  end interface
  
  interface
    subroutine psb_s_cuda_diag_scals(d,a,info) 
      import :: psb_s_cuda_diag_sparse_mat, psb_spk_, psb_ipk_
      class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)                 :: d
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_s_cuda_diag_scals
  end interface
  

contains 

  ! == ===================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  ! == ===================================

  
  function s_cuda_diag_sizeof(a) result(res)
    implicit none 
    class(psb_s_cuda_diag_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res

    res = 8 
    res = res + psb_sizeof_sp  * size(a%data)
    res = res + psb_sizeof_ip * size(a%offset)

    ! Should we account for the shadow data structure
    ! on the GPU device side? 
    ! res = 2*res
      
  end function s_cuda_diag_sizeof

  function s_cuda_diag_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'DIAG'
  end function s_cuda_diag_get_fmt
  


  ! == ===================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  ! == ===================================  

  subroutine  s_cuda_diag_free(a) 
    use diagdev_mod
    implicit none 
    integer(psb_ipk_) :: info
    class(psb_s_cuda_diag_sparse_mat), intent(inout) :: a

    if (c_associated(a%deviceMat)) &
         & call freeDiagDevice(a%deviceMat)
    a%deviceMat = c_null_ptr
    call a%psb_s_dia_sparse_mat%free()
    
    return

  end subroutine s_cuda_diag_free

  subroutine  s_cuda_diag_finalize(a) 
    use diagdev_mod
    implicit none 
    type(psb_s_cuda_diag_sparse_mat), intent(inout) :: a

    if (c_associated(a%deviceMat)) &
         &  call freeDiagDevice(a%deviceMat)
    a%deviceMat = c_null_ptr
    
    return
  end subroutine s_cuda_diag_finalize

end module psb_s_cuda_diag_mat_mod
