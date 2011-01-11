module psb_d_elg_mat_mod

  use iso_c_binding
  use psb_d_base_mat_mod
  use psb_d_ell_mat_mod

  type, extends(psb_d_ell_sparse_mat) :: psb_d_elg_sparse_mat
    !
    ! ITPACK/ELL format, extended.
    ! We are adding here the routines to create a copy of the data
    ! into the GPU. 
    ! If HAVE_ELL_GPU is undefined this is just
    ! a copy of ELL, indistinguishable.
    ! 
#ifdef HAVE_ELL_GPU
    type(c_ptr) :: deviceMat = c_null_ptr

  contains
    procedure, pass(a) :: get_fmt      => d_elg_get_fmt
    procedure, pass(a) :: sizeof       => d_elg_sizeof
    procedure, pass(a) :: d_csmm       => psb_d_elg_csmm
    procedure, pass(a) :: d_csmv       => psb_d_elg_csmv
    procedure, pass(a) :: d_scals      => psb_d_elg_scals
    procedure, pass(a) :: d_scal       => psb_d_elg_scal
    procedure, pass(a) :: reallocate_nz => psb_d_elg_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_d_elg_allocate_mnnz
    procedure, pass(a) :: cp_from_coo  => psb_d_cp_elg_from_coo
    procedure, pass(a) :: cp_from_fmt  => psb_d_cp_elg_from_fmt
    procedure, pass(a) :: mv_from_coo  => psb_d_mv_elg_from_coo
    procedure, pass(a) :: mv_from_fmt  => psb_d_mv_elg_from_fmt
    procedure, pass(a) :: free         => d_elg_free
    procedure, pass(a) :: mold         => psb_d_elg_mold
    procedure, pass(a) :: psb_d_elg_cp_from
    generic, public    :: cp_from => psb_d_elg_cp_from
    procedure, pass(a) :: psb_d_elg_mv_from
    generic, public    :: mv_from => psb_d_elg_mv_from
#endif
  end type psb_d_elg_sparse_mat

#ifdef HAVE_ELL_GPU
  private :: d_elg_get_nzeros, d_elg_free,  d_elg_get_fmt, &
       & d_elg_get_size, d_elg_sizeof, d_elg_get_nz_row

  interface
    subroutine  psb_d_elg_reallocate_nz(nz,a) 
      import :: psb_d_elg_sparse_mat
      integer, intent(in) :: nz
      class(psb_d_elg_sparse_mat), intent(inout) :: a
    end subroutine psb_d_elg_reallocate_nz
  end interface

  interface
    subroutine  psb_d_elg_allocate_mnnz(m,n,a,nz) 
      import :: psb_d_elg_sparse_mat
      integer, intent(in) :: m,n
      class(psb_d_elg_sparse_mat), intent(inout) :: a
      integer, intent(in), optional :: nz
    end subroutine psb_d_elg_allocate_mnnz
  end interface

  interface 
    subroutine psb_d_elg_mold(a,b,info) 
      import :: psb_d_elg_sparse_mat, psb_d_base_sparse_mat, psb_long_int_k_
      class(psb_d_elg_sparse_mat), intent(in)               :: a
      class(psb_d_base_sparse_mat), intent(out), allocatable :: b
      integer, intent(out)                                 :: info
    end subroutine psb_d_elg_mold
  end interface

  interface 
    subroutine psb_d_cp_elg_from_coo(a,b,info) 
      import :: psb_d_elg_sparse_mat, psb_d_coo_sparse_mat
      class(psb_d_elg_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_cp_elg_from_coo
  end interface
  
  interface 
    subroutine psb_d_cp_elg_from_fmt(a,b,info) 
      import :: psb_d_elg_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_elg_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(in)   :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_cp_elg_from_fmt
  end interface
  
  interface 
    subroutine psb_d_mv_elg_from_coo(a,b,info) 
      import :: psb_d_elg_sparse_mat, psb_d_coo_sparse_mat
      class(psb_d_elg_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_mv_elg_from_coo
  end interface
  

  interface 
    subroutine psb_d_mv_elg_from_fmt(a,b,info) 
      import :: psb_d_elg_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_elg_sparse_mat), intent(inout)  :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine psb_d_mv_elg_from_fmt
  end interface
  
  interface 
    subroutine psb_d_elg_cp_from(a,b)
      import :: psb_d_elg_sparse_mat, psb_dpk_
      class(psb_d_elg_sparse_mat), intent(inout) :: a
      type(psb_d_elg_sparse_mat), intent(in)   :: b
    end subroutine psb_d_elg_cp_from
  end interface
  
  interface 
    subroutine psb_d_elg_mv_from(a,b)
      import :: psb_d_elg_sparse_mat, psb_dpk_
      class(psb_d_elg_sparse_mat), intent(inout)  :: a
      type(psb_d_elg_sparse_mat), intent(inout) :: b
    end subroutine psb_d_elg_mv_from
  end interface
  
  interface 
    subroutine psb_d_elg_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_d_elg_sparse_mat, psb_dpk_
      class(psb_d_elg_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_elg_csmv
  end interface
  interface 
    subroutine psb_d_elg_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_d_elg_sparse_mat, psb_dpk_
      class(psb_d_elg_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_elg_csmm
  end interface
  
  interface 
    subroutine psb_d_elg_scal(d,a,info) 
      import :: psb_d_elg_sparse_mat, psb_dpk_
      class(psb_d_elg_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_d_elg_scal
  end interface
  
  interface
    subroutine psb_d_elg_scals(d,a,info) 
      import :: psb_d_elg_sparse_mat, psb_dpk_
      class(psb_d_elg_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d
      integer, intent(out)            :: info
    end subroutine psb_d_elg_scals
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

  
  function d_elg_sizeof(a) result(res)
    implicit none 
    class(psb_d_elg_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8 
    res = res + psb_sizeof_dp  * size(a%val)
    res = res + psb_sizeof_int * size(a%irn)
    res = res + psb_sizeof_int * size(a%idiag)
    res = res + psb_sizeof_int * size(a%ja)
    ! Should we account for the shadow data structure
    ! on the GPU device side? 
    ! res = 2*res
      
  end function d_elg_sizeof

  function d_elg_get_fmt(a) result(res)
    implicit none 
    class(psb_d_elg_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'ELG'
  end function d_elg_get_fmt
  
!!$  function d_elg_get_nzeros(a) result(res)
!!$    implicit none 
!!$    class(psb_d_elg_sparse_mat), intent(in) :: a
!!$    integer :: res
!!$    res = sum(a%irn(1:a%get_nrows()))
!!$  end function d_elg_get_nzeros
!!$
!!$  function d_elg_get_size(a) result(res)
!!$    implicit none 
!!$    class(psb_d_elg_sparse_mat), intent(in) :: a
!!$    integer :: res
!!$
!!$    res = -1
!!$    
!!$    if (allocated(a%ja)) then 
!!$      if (res >= 0) then 
!!$        res = min(res,size(a%ja))
!!$      else 
!!$        res = size(a%ja)
!!$      end if
!!$    end if
!!$    if (allocated(a%val)) then 
!!$      if (res >= 0) then 
!!$        res = min(res,size(a%val))
!!$      else 
!!$        res = size(a%val)
!!$      end if
!!$    end if
!!$
!!$  end function d_elg_get_size
!!$
!!$
!!$
!!$  function  d_elg_get_nz_row(idx,a) result(res)
!!$
!!$    implicit none
!!$    
!!$    class(psb_d_elg_sparse_mat), intent(in) :: a
!!$    integer, intent(in)                  :: idx
!!$    integer                              :: res
!!$    
!!$    res = 0 
!!$ 
!!$    if ((1<=idx).and.(idx<=a%get_nrows())) then 
!!$      res = a%irn(idx)
!!$    end if
!!$    
!!$  end function d_elg_get_nz_row
!!$


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

  subroutine  d_elg_free(a) 
#ifdef HAVE_ELL_GPU
    use elldev_mod
#endif
    implicit none 
    integer :: info

    class(psb_d_elg_sparse_mat), intent(inout) :: a

    if (allocated(a%idiag)) deallocate(a%idiag)
    if (allocated(a%irn)) deallocate(a%irn)
    if (allocated(a%ja))  deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
#ifdef HAVE_ELL_GPU 
    call freeEllDevice(a%deviceMat)
#endif
    
    return

  end subroutine d_elg_free

#endif

end module psb_d_elg_mat_mod
