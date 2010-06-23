module psb_d_csc_mat_mod

  use psb_d_base_mat_mod

  type, extends(psb_d_base_sparse_mat) :: psb_d_csc_sparse_mat

    integer, allocatable :: icp(:), ia(:)
    real(psb_dpk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: get_size     => d_csc_get_size
    procedure, pass(a) :: get_nzeros   => d_csc_get_nzeros
    procedure, pass(a) :: get_fmt      => d_csc_get_fmt
    procedure, pass(a) :: sizeof       => d_csc_sizeof
    procedure, pass(a) :: d_csmm       => psb_d_csc_csmm
    procedure, pass(a) :: d_csmv       => psb_d_csc_csmv
    procedure, pass(a) :: d_inner_cssm => psb_d_csc_cssm
    procedure, pass(a) :: d_inner_cssv => psb_d_csc_cssv
    procedure, pass(a) :: d_scals      => psb_d_csc_scals
    procedure, pass(a) :: d_scal       => psb_d_csc_scal
    procedure, pass(a) :: csnmi        => psb_d_csc_csnmi
    procedure, pass(a) :: csnm1        => psb_d_csc_csnm1
    procedure, pass(a) :: rowsum       => psb_d_csc_rowsum
    procedure, pass(a) :: arwsum       => psb_d_csc_arwsum
    procedure, pass(a) :: colsum       => psb_d_csc_colsum
    procedure, pass(a) :: aclsum       => psb_d_csc_aclsum
    procedure, pass(a) :: reallocate_nz => psb_d_csc_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_d_csc_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_d_cp_csc_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_d_cp_csc_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_d_cp_csc_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_d_cp_csc_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_d_mv_csc_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_d_mv_csc_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_d_mv_csc_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_d_mv_csc_from_fmt
    procedure, pass(a) :: csput        => psb_d_csc_csput
    procedure, pass(a) :: get_diag     => psb_d_csc_get_diag
    procedure, pass(a) :: csgetptn     => psb_d_csc_csgetptn
    procedure, pass(a) :: d_csgetrow   => psb_d_csc_csgetrow
    procedure, pass(a) :: get_nz_col   => d_csc_get_nz_col
    procedure, pass(a) :: reinit       => psb_d_csc_reinit
    procedure, pass(a) :: trim         => psb_d_csc_trim
    procedure, pass(a) :: print        => psb_d_csc_print
    procedure, pass(a) :: free         => d_csc_free
    procedure, pass(a) :: psb_d_csc_cp_from
    generic, public    :: cp_from => psb_d_csc_cp_from
    procedure, pass(a) :: psb_d_csc_mv_from
    generic, public    :: mv_from => psb_d_csc_mv_from

  end type psb_d_csc_sparse_mat

 private :: d_csc_get_nzeros, d_csc_free,  d_csc_get_fmt, &
       & d_csc_get_size, d_csc_sizeof, d_csc_get_nz_col

  interface
    subroutine  psb_d_csc_reallocate_nz(nz,a) 
      import psb_d_csc_sparse_mat
      integer, intent(in) :: nz
      class(psb_d_csc_sparse_mat), intent(inout) :: a
    end subroutine psb_d_csc_reallocate_nz
  end interface
  
  interface 
    subroutine psb_d_csc_reinit(a,clear)
      import psb_d_csc_sparse_mat
      class(psb_d_csc_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_d_csc_reinit
  end interface
  
  interface
    subroutine  psb_d_csc_trim(a)
      import psb_d_csc_sparse_mat
      class(psb_d_csc_sparse_mat), intent(inout) :: a
    end subroutine psb_d_csc_trim
  end interface
  
  interface
    subroutine  psb_d_csc_allocate_mnnz(m,n,a,nz) 
      import psb_d_csc_sparse_mat
      integer, intent(in) :: m,n
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      integer, intent(in), optional :: nz
    end subroutine psb_d_csc_allocate_mnnz
  end interface
  
  interface
    subroutine psb_d_csc_print(iout,a,iv,eirs,eics,head,ivr,ivc)
      import psb_d_csc_sparse_mat
      integer, intent(in)               :: iout
      class(psb_d_csc_sparse_mat), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_d_csc_print
  end interface
  
  interface 
    subroutine psb_d_cp_csc_to_coo(a,b,info) 
      import psb_d_coo_sparse_mat, psb_d_csc_sparse_mat
      class(psb_d_csc_sparse_mat), intent(in) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_d_cp_csc_to_coo
  end interface
  
  interface 
    subroutine psb_d_cp_csc_from_coo(a,b,info) 
      import psb_d_csc_sparse_mat, psb_d_coo_sparse_mat
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_cp_csc_from_coo
  end interface
  
  interface 
    subroutine psb_d_cp_csc_to_fmt(a,b,info) 
      import psb_d_csc_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_csc_sparse_mat), intent(in)   :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                       :: info
    end subroutine psb_d_cp_csc_to_fmt
  end interface
  
  interface 
    subroutine psb_d_cp_csc_from_fmt(a,b,info) 
      import psb_d_csc_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(in)   :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_cp_csc_from_fmt
  end interface
  
  interface 
    subroutine psb_d_mv_csc_to_coo(a,b,info) 
      import psb_d_csc_sparse_mat, psb_d_coo_sparse_mat
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout)   :: b
      integer, intent(out)            :: info
    end subroutine psb_d_mv_csc_to_coo
  end interface
  
  interface 
    subroutine psb_d_mv_csc_from_coo(a,b,info) 
      import psb_d_csc_sparse_mat, psb_d_coo_sparse_mat
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_mv_csc_from_coo
  end interface
  
  interface 
    subroutine psb_d_mv_csc_to_fmt(a,b,info) 
      import psb_d_csc_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(inout)  :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_mv_csc_to_fmt
  end interface
  
  interface 
    subroutine psb_d_mv_csc_from_fmt(a,b,info) 
      import psb_d_csc_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_csc_sparse_mat), intent(inout)  :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine psb_d_mv_csc_from_fmt
  end interface
  
  interface 
    subroutine psb_d_csc_cp_from(a,b)
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      type(psb_d_csc_sparse_mat), intent(in)   :: b
    end subroutine psb_d_csc_cp_from
  end interface
  
  interface 
    subroutine psb_d_csc_mv_from(a,b)
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(inout)  :: a
      type(psb_d_csc_sparse_mat), intent(inout) :: b
    end subroutine psb_d_csc_mv_from
  end interface
  
  
  interface 
    subroutine psb_d_csc_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psb_d_csc_csput
  end interface
  
  interface 
    subroutine psb_d_csc_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_csc_csgetptn
  end interface
  
  interface 
    subroutine psb_d_csc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_csc_csgetrow
  end interface

  interface 
    subroutine psb_d_csc_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import psb_d_csc_sparse_mat, psb_dpk_, psb_d_coo_sparse_mat
      class(psb_d_csc_sparse_mat), intent(in) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(in)                  :: imin,imax
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_csc_csgetblk
  end interface
    
  interface 
    subroutine psb_d_csc_cssv(alpha,a,x,beta,y,info,trans) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_csc_cssv
    subroutine psb_d_csc_cssm(alpha,a,x,beta,y,info,trans) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_csc_cssm
  end interface
  
  interface 
    subroutine psb_d_csc_csmv(alpha,a,x,beta,y,info,trans) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_csc_csmv
    subroutine psb_d_csc_csmm(alpha,a,x,beta,y,info,trans) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_csc_csmm
  end interface
  
  
  interface 
    function psb_d_csc_csnmi(a) result(res)
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_d_csc_csnmi
  end interface
  
  interface 
    function psb_d_csc_csnm1(a) result(res)
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_d_csc_csnm1
  end interface

  interface 
    subroutine psb_d_csc_rowsum(d,a) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_csc_rowsum
  end interface

  interface 
    subroutine psb_d_csc_arwsum(d,a) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_csc_arwsum
  end interface
  
  interface 
    subroutine psb_d_csc_colsum(d,a) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_csc_colsum
  end interface

  interface 
    subroutine psb_d_csc_aclsum(d,a) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_csc_aclsum
  end interface
    
  interface 
    subroutine psb_d_csc_get_diag(a,d,info) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)     :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_d_csc_get_diag
  end interface
  
  interface 
    subroutine psb_d_csc_scal(d,a,info) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_d_csc_scal
  end interface
  
  interface
    subroutine psb_d_csc_scals(d,a,info) 
      import psb_d_csc_sparse_mat, psb_dpk_
      class(psb_d_csc_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d
      integer, intent(out)            :: info
    end subroutine psb_d_csc_scals
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

  
  function d_csc_sizeof(a) result(res)
    implicit none 
    class(psb_d_csc_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8 
    res = res + psb_sizeof_dp  * size(a%val)
    res = res + psb_sizeof_int * size(a%icp)
    res = res + psb_sizeof_int * size(a%ia)
      
  end function d_csc_sizeof

  function d_csc_get_fmt(a) result(res)
    implicit none 
    class(psb_d_csc_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'CSC'
  end function d_csc_get_fmt
  
  function d_csc_get_nzeros(a) result(res)
    implicit none 
    class(psb_d_csc_sparse_mat), intent(in) :: a
    integer :: res
    res = a%icp(a%get_ncols()+1)-1
  end function d_csc_get_nzeros

  function d_csc_get_size(a) result(res)
    implicit none 
    class(psb_d_csc_sparse_mat), intent(in) :: a
    integer :: res

    res = -1
    
    if (allocated(a%ia)) then 
      if (res >= 0) then 
        res = min(res,size(a%ia))
      else 
        res = size(a%ia)
      end if
    end if
    if (allocated(a%val)) then 
      if (res >= 0) then 
        res = min(res,size(a%val))
      else 
        res = size(a%val)
      end if
    end if

  end function d_csc_get_size



  function  d_csc_get_nz_col(idx,a) result(res)
    use psb_const_mod
    implicit none
    
    class(psb_d_csc_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: idx
    integer                              :: res
    
    res = 0 
 
    if ((1<=idx).and.(idx<=a%get_ncols())) then 
      res = a%icp(idx+1)-a%icp(idx)
    end if
    
  end function d_csc_get_nz_col



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


  subroutine  d_csc_free(a) 
    implicit none 

    class(psb_d_csc_sparse_mat), intent(inout) :: a

    if (allocated(a%icp)) deallocate(a%icp)
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
    
    return

  end subroutine d_csc_free

end module psb_d_csc_mat_mod
