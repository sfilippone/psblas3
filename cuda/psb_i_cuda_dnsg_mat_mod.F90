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
  

module psb_i_cuda_dnsg_mat_mod

  use iso_c_binding
  use psb_i_mat_mod 
  use psb_i_dns_mat_mod
  use psb_cuda_env_mod  
  use dnsdev_mod
  
  type, extends(psb_i_dns_sparse_mat) :: psb_i_cuda_dnsg_sparse_mat
    !
    ! ITPACK/DNS format, extended.
    ! We are adding here the routines to create a copy of the data
    ! into the GPU. 
    ! If HAVE_SPGPU is undefined this is just
    ! a copy of DNS, indistinguishable.
    ! 
    type(c_ptr) :: deviceMat = c_null_ptr

  contains
    procedure, nopass  :: get_fmt       => i_cuda_dnsg_get_fmt
    ! procedure, pass(a) :: sizeof        => i_cuda_dnsg_sizeof
    procedure, pass(a) :: vect_mv       => psb_i_cuda_dnsg_vect_mv
!!$    procedure, pass(a) :: csmm          => psb_i_cuda_dnsg_csmm
!!$    procedure, pass(a) :: csmv          => psb_i_cuda_dnsg_csmv
!!$    procedure, pass(a) :: in_vect_sv    => psb_i_cuda_dnsg_inner_vect_sv
!!$    procedure, pass(a) :: scals         => psb_i_cuda_dnsg_scals
!!$    procedure, pass(a) :: scalv         => psb_i_cuda_dnsg_scal
!!$    procedure, pass(a) :: reallocate_nz => psb_i_cuda_dnsg_reallocate_nz
!!$    procedure, pass(a) :: allocate_mnnz => psb_i_cuda_dnsg_allocate_mnnz
    ! Note: we *do* need the TO methods, because of the need to invoke SYNC
    ! 
    procedure, pass(a) :: cp_from_coo   => psb_i_cuda_cp_dnsg_from_coo
    procedure, pass(a) :: cp_from_fmt   => psb_i_cuda_cp_dnsg_from_fmt
    procedure, pass(a) :: mv_from_coo   => psb_i_cuda_mv_dnsg_from_coo
    procedure, pass(a) :: mv_from_fmt   => psb_i_cuda_mv_dnsg_from_fmt
    procedure, pass(a) :: free          => i_cuda_dnsg_free
    procedure, pass(a) :: mold          => psb_i_cuda_dnsg_mold
    procedure, pass(a) :: to_gpu        => psb_i_cuda_dnsg_to_gpu
    final              :: i_cuda_dnsg_finalize
  end type psb_i_cuda_dnsg_sparse_mat

  private :: i_cuda_dnsg_get_nzeros, i_cuda_dnsg_free,  i_cuda_dnsg_get_fmt, &
       & i_cuda_dnsg_get_size, i_cuda_dnsg_get_nz_row


  interface 
    subroutine psb_i_cuda_dnsg_vect_mv(alpha,a,x,beta,y,info,trans) 
      import :: psb_i_cuda_dnsg_sparse_mat, psb_ipk_, psb_i_base_vect_type, psb_ipk_
      class(psb_i_cuda_dnsg_sparse_mat), intent(in)    :: a
      integer(psb_ipk_), intent(in)                 :: alpha, beta
      class(psb_i_base_vect_type), intent(inout) :: x
      class(psb_i_base_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)             :: info
      character, optional, intent(in)            :: trans
    end subroutine psb_i_cuda_dnsg_vect_mv
  end interface
!!$
!!$  interface 
!!$    subroutine psb_i_cuda_dnsg_inner_vect_sv(alpha,a,x,beta,y,info,trans) 
!!$      import :: psb_ipk_, psb_i_cuda_dnsg_sparse_mat, psb_ipk_,  psb_i_base_vect_type
!!$      class(psb_i_cuda_dnsg_sparse_mat), intent(in)    :: a
!!$      integer(psb_ipk_), intent(in)                 :: alpha, beta
!!$      class(psb_i_base_vect_type), intent(inout) :: x, y
!!$      integer(psb_ipk_), intent(out)             :: info
!!$      character, optional, intent(in)            :: trans
!!$    end subroutine psb_i_cuda_dnsg_inner_vect_sv
!!$  end interface

!!$  interface
!!$    subroutine  psb_i_cuda_dnsg_reallocate_nz(nz,a) 
!!$      import :: psb_i_cuda_dnsg_sparse_mat, psb_ipk_
!!$      integer(psb_ipk_), intent(in)              :: nz
!!$      class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a
!!$    end subroutine psb_i_cuda_dnsg_reallocate_nz
!!$  end interface
!!$
!!$  interface
!!$    subroutine  psb_i_cuda_dnsg_allocate_mnnz(m,n,a,nz) 
!!$      import :: psb_i_cuda_dnsg_sparse_mat, psb_ipk_
!!$      integer(psb_ipk_), intent(in)              :: m,n
!!$      class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a
!!$      integer(psb_ipk_), intent(in), optional    :: nz
!!$    end subroutine psb_i_cuda_dnsg_allocate_mnnz
!!$  end interface

  interface 
    subroutine psb_i_cuda_dnsg_mold(a,b,info) 
      import :: psb_i_cuda_dnsg_sparse_mat, psb_i_base_sparse_mat, psb_ipk_
      class(psb_i_cuda_dnsg_sparse_mat), intent(in)                  :: a
      class(psb_i_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out)                           :: info
    end subroutine psb_i_cuda_dnsg_mold
  end interface

  interface 
    subroutine psb_i_cuda_dnsg_to_gpu(a,info) 
      import :: psb_i_cuda_dnsg_sparse_mat, psb_ipk_
      class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_dnsg_to_gpu
  end interface

  interface 
    subroutine psb_i_cuda_cp_dnsg_from_coo(a,b,info) 
      import :: psb_i_cuda_dnsg_sparse_mat, psb_i_coo_sparse_mat, psb_ipk_
      class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a
      class(psb_i_coo_sparse_mat), intent(in)    :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_cp_dnsg_from_coo
  end interface
  
  interface 
    subroutine psb_i_cuda_cp_dnsg_from_fmt(a,b,info) 
      import :: psb_i_cuda_dnsg_sparse_mat, psb_i_base_sparse_mat, psb_ipk_
      class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a
      class(psb_i_base_sparse_mat), intent(in)   :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_cp_dnsg_from_fmt
  end interface
  
  interface 
    subroutine psb_i_cuda_mv_dnsg_from_coo(a,b,info) 
      import :: psb_i_cuda_dnsg_sparse_mat, psb_i_coo_sparse_mat, psb_ipk_
      class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a
      class(psb_i_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_i_cuda_mv_dnsg_from_coo
  end interface
  

  interface 
    subroutine psb_i_cuda_mv_dnsg_from_fmt(a,b,info) 
      import :: psb_i_cuda_dnsg_sparse_mat, psb_i_base_sparse_mat, psb_ipk_
      class(psb_i_cuda_dnsg_sparse_mat), intent(inout)  :: a
      class(psb_i_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_i_cuda_mv_dnsg_from_fmt
  end interface
  
!!$  interface 
!!$    subroutine psb_i_cuda_dnsg_csmv(alpha,a,x,beta,y,info,trans) 
!!$      import :: psb_i_cuda_dnsg_sparse_mat, psb_ipk_, psb_ipk_
!!$      class(psb_i_cuda_dnsg_sparse_mat), intent(in) :: a
!!$      integer(psb_ipk_), intent(in)              :: alpha, beta, x(:)
!!$      integer(psb_ipk_), intent(inout)           :: y(:)
!!$      integer(psb_ipk_), intent(out)          :: info
!!$      character, optional, intent(in)         :: trans
!!$    end subroutine psb_i_cuda_dnsg_csmv
!!$  end interface
!!$  interface 
!!$    subroutine psb_i_cuda_dnsg_csmm(alpha,a,x,beta,y,info,trans) 
!!$      import :: psb_i_cuda_dnsg_sparse_mat, psb_ipk_, psb_ipk_
!!$      class(psb_i_cuda_dnsg_sparse_mat), intent(in) :: a
!!$      integer(psb_ipk_), intent(in)              :: alpha, beta, x(:,:)
!!$      integer(psb_ipk_), intent(inout)           :: y(:,:)
!!$      integer(psb_ipk_), intent(out)          :: info
!!$      character, optional, intent(in)         :: trans
!!$    end subroutine psb_i_cuda_dnsg_csmm
!!$  end interface
!!$  
!!$  interface 
!!$    subroutine psb_i_cuda_dnsg_scal(d,a,info, side) 
!!$      import :: psb_i_cuda_dnsg_sparse_mat, psb_ipk_, psb_ipk_
!!$      class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a
!!$      integer(psb_ipk_), intent(in)                 :: d(:)
!!$      integer(psb_ipk_), intent(out)             :: info
!!$      character, intent(in), optional            :: side
!!$    end subroutine psb_i_cuda_dnsg_scal
!!$  end interface
!!$  
!!$  interface
!!$    subroutine psb_i_cuda_dnsg_scals(d,a,info) 
!!$      import :: psb_i_cuda_dnsg_sparse_mat, psb_ipk_, psb_ipk_
!!$      class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a
!!$      integer(psb_ipk_), intent(in)                 :: d
!!$      integer(psb_ipk_), intent(out)             :: info
!!$    end subroutine psb_i_cuda_dnsg_scals
!!$  end interface
!!$  

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

  

  function i_cuda_dnsg_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'DNSG'
  end function i_cuda_dnsg_get_fmt
  


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

  subroutine  i_cuda_dnsg_free(a) 
    use dnsdev_mod
    implicit none 
    integer(psb_ipk_) :: info
    class(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a

    if (c_associated(a%deviceMat)) then
      call trackCudaFree(' i_dnsg ',a%sizeof())         
      call freeDnsDevice(a%deviceMat)
    end if
    a%deviceMat = c_null_ptr
    call a%psb_i_dns_sparse_mat%free()
    
    return

  end subroutine i_cuda_dnsg_free

  subroutine  i_cuda_dnsg_finalize(a) 
    use dnsdev_mod
    implicit none 
    type(psb_i_cuda_dnsg_sparse_mat), intent(inout) :: a

    if (c_associated(a%deviceMat)) &
         & call freeDnsDevice(a%deviceMat)
    a%deviceMat = c_null_ptr
    
    return
  end subroutine i_cuda_dnsg_finalize

end module psb_i_cuda_dnsg_mat_mod
