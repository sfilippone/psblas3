module psb_c_csr_mat_mod

  use psb_c_base_mat_mod

  type, extends(psb_c_base_sparse_mat) :: psb_c_csr_sparse_mat

    integer, allocatable :: irp(:), ja(:)
    complex(psb_spk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: get_size     => c_csr_get_size
    procedure, pass(a) :: get_nzeros   => c_csr_get_nzeros
    procedure, pass(a) :: get_fmt      => c_csr_get_fmt
    procedure, pass(a) :: sizeof       => c_csr_sizeof
    procedure, pass(a) :: c_csmm       => psb_c_csr_csmm
    procedure, pass(a) :: c_csmv       => psb_c_csr_csmv
    procedure, pass(a) :: c_inner_cssm => psb_c_csr_cssm
    procedure, pass(a) :: c_inner_cssv => psb_c_csr_cssv
    procedure, pass(a) :: c_scals      => psb_c_csr_scals
    procedure, pass(a) :: c_scal       => psb_c_csr_scal
    procedure, pass(a) :: csnmi        => psb_c_csr_csnmi
    procedure, pass(a) :: reallocate_nz => psb_c_csr_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_c_csr_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_c_cp_csr_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_c_cp_csr_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_c_cp_csr_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_c_cp_csr_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_c_mv_csr_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_c_mv_csr_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_c_mv_csr_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_c_mv_csr_from_fmt
    procedure, pass(a) :: csput        => psb_c_csr_csput
    procedure, pass(a) :: get_diag     => psb_c_csr_get_diag
    procedure, pass(a) :: csgetptn     => psb_c_csr_csgetptn
    procedure, pass(a) :: c_csgetrow   => psb_c_csr_csgetrow
    procedure, pass(a) :: get_nc_row   => c_csr_get_nc_row
    procedure, pass(a) :: reinit       => psb_c_csr_reinit
    procedure, pass(a) :: trim         => psb_c_csr_trim
    procedure, pass(a) :: print        => psb_c_csr_print
    procedure, pass(a) :: free         => c_csr_free
    procedure, pass(a) :: psb_c_csr_cp_from
    generic, public    :: cp_from => psb_c_csr_cp_from
    procedure, pass(a) :: psb_c_csr_mv_from
    generic, public    :: mv_from => psb_c_csr_mv_from

  end type psb_c_csr_sparse_mat

  private :: c_csr_get_nzeros, c_csr_free,  c_csr_get_fmt, &
       & c_csr_get_size, c_csr_sizeof, c_csr_get_nc_row

  interface
    subroutine  psb_c_csr_reallocate_nz(nz,a) 
      import psb_c_csr_sparse_mat
      integer, intent(in) :: nz
      class(psb_c_csr_sparse_mat), intent(inout) :: a
    end subroutine psb_c_csr_reallocate_nz
  end interface
  
  interface 
    subroutine psb_c_csr_reinit(a,clear)
      import psb_c_csr_sparse_mat
      class(psb_c_csr_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_c_csr_reinit
  end interface
  
  interface
    subroutine  psb_c_csr_trim(a)
      import psb_c_csr_sparse_mat
      class(psb_c_csr_sparse_mat), intent(inout) :: a
    end subroutine psb_c_csr_trim
  end interface
  
  interface
    subroutine  psb_c_csr_allocate_mnnz(m,n,a,nz) 
      import psb_c_csr_sparse_mat
      integer, intent(in) :: m,n
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      integer, intent(in), optional :: nz
    end subroutine psb_c_csr_allocate_mnnz
  end interface
  
  interface
    subroutine psb_c_csr_print(iout,a,iv,eirs,eics,head,ivr,ivc)
      import psb_c_csr_sparse_mat
      integer, intent(in)               :: iout
      class(psb_c_csr_sparse_mat), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_c_csr_print
  end interface
  
  interface 
    subroutine psb_c_cp_csr_to_coo(a,b,info) 
      import psb_c_coo_sparse_mat, psb_c_csr_sparse_mat
      class(psb_c_csr_sparse_mat), intent(in) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_c_cp_csr_to_coo
  end interface
  
  interface 
    subroutine psb_c_cp_csr_from_coo(a,b,info) 
      import psb_c_csr_sparse_mat, psb_c_coo_sparse_mat
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine psb_c_cp_csr_from_coo
  end interface
  
  interface 
    subroutine psb_c_cp_csr_to_fmt(a,b,info) 
      import psb_c_csr_sparse_mat, psb_c_base_sparse_mat
      class(psb_c_csr_sparse_mat), intent(in)   :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                       :: info
    end subroutine psb_c_cp_csr_to_fmt
  end interface
  
  interface 
    subroutine psb_c_cp_csr_from_fmt(a,b,info) 
      import psb_c_csr_sparse_mat, psb_c_base_sparse_mat
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(in)   :: b
      integer, intent(out)                        :: info
    end subroutine psb_c_cp_csr_from_fmt
  end interface
  
  interface 
    subroutine psb_c_mv_csr_to_coo(a,b,info) 
      import psb_c_csr_sparse_mat, psb_c_coo_sparse_mat
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout)   :: b
      integer, intent(out)            :: info
    end subroutine psb_c_mv_csr_to_coo
  end interface
  
  interface 
    subroutine psb_c_mv_csr_from_coo(a,b,info) 
      import psb_c_csr_sparse_mat, psb_c_coo_sparse_mat
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine psb_c_mv_csr_from_coo
  end interface
  
  interface 
    subroutine psb_c_mv_csr_to_fmt(a,b,info) 
      import psb_c_csr_sparse_mat, psb_c_base_sparse_mat
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(inout)  :: b
      integer, intent(out)                        :: info
    end subroutine psb_c_mv_csr_to_fmt
  end interface
  
  interface 
    subroutine psb_c_mv_csr_from_fmt(a,b,info) 
      import psb_c_csr_sparse_mat, psb_c_base_sparse_mat
      class(psb_c_csr_sparse_mat), intent(inout)  :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine psb_c_mv_csr_from_fmt
  end interface
  
  interface 
    subroutine psb_c_csr_cp_from(a,b)
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      type(psb_c_csr_sparse_mat), intent(in)   :: b
    end subroutine psb_c_csr_cp_from
  end interface
  
  interface 
    subroutine psb_c_csr_mv_from(a,b)
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(inout)  :: a
      type(psb_c_csr_sparse_mat), intent(inout) :: b
    end subroutine psb_c_csr_mv_from
  end interface
  
  
  interface 
    subroutine psb_c_csr_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psb_c_csr_csput
  end interface
  
  interface 
    subroutine psb_c_csr_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_csr_csgetptn
  end interface
  
  interface 
    subroutine psb_c_csr_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_csr_csgetrow
  end interface

  interface 
    subroutine psb_c_csr_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import psb_c_csr_sparse_mat, psb_spk_, psb_c_coo_sparse_mat
      class(psb_c_csr_sparse_mat), intent(in) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer, intent(in)                  :: imin,imax
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_csr_csgetblk
  end interface
    
  interface 
    subroutine psb_c_csr_cssv(alpha,a,x,beta,y,info,trans) 
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_c_csr_cssv
    subroutine psb_c_csr_cssm(alpha,a,x,beta,y,info,trans) 
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_c_csr_cssm
  end interface
  
  interface 
    subroutine psb_c_csr_csmv(alpha,a,x,beta,y,info,trans) 
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_c_csr_csmv
    subroutine psb_c_csr_csmm(alpha,a,x,beta,y,info,trans) 
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_c_csr_csmm
  end interface
  
  
  interface 
    function psb_c_csr_csnmi(a) result(res)
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_c_csr_csnmi
  end interface
  
  interface 
    subroutine psb_c_csr_get_diag(a,d,info) 
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)     :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_c_csr_get_diag
  end interface
  
  interface 
    subroutine psb_c_csr_scal(d,a,info) 
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_c_csr_scal
  end interface
  
  interface
    subroutine psb_c_csr_scals(d,a,info) 
      import psb_c_csr_sparse_mat, psb_spk_
      class(psb_c_csr_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: d
      integer, intent(out)            :: info
    end subroutine psb_c_csr_scals
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

  
  function c_csr_sizeof(a) result(res)
    implicit none 
    class(psb_c_csr_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8 
    res = res + 2 * psb_sizeof_sp  * size(a%val)
    res = res + psb_sizeof_int * size(a%irp)
    res = res + psb_sizeof_int * size(a%ja)
      
  end function c_csr_sizeof

  function c_csr_get_fmt(a) result(res)
    implicit none 
    class(psb_c_csr_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'CSR'
  end function c_csr_get_fmt
  
  function c_csr_get_nzeros(a) result(res)
    implicit none 
    class(psb_c_csr_sparse_mat), intent(in) :: a
    integer :: res
    res = a%irp(a%get_nrows()+1)-1
  end function c_csr_get_nzeros

  function c_csr_get_size(a) result(res)
    implicit none 
    class(psb_c_csr_sparse_mat), intent(in) :: a
    integer :: res

    res = -1
    
    if (allocated(a%ja)) then 
      if (res >= 0) then 
        res = min(res,size(a%ja))
      else 
        res = size(a%ja)
      end if
    end if
    if (allocated(a%val)) then 
      if (res >= 0) then 
        res = min(res,size(a%val))
      else 
        res = size(a%val)
      end if
    end if

  end function c_csr_get_size



  function  c_csr_get_nc_row(idx,a) result(res)

    implicit none
    
    class(psb_c_csr_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: idx
    integer                              :: res
    
    res = 0 
 
    if ((1<=idx).and.(idx<=a%get_nrows())) then 
      res = a%irp(idx+1)-a%irp(idx)
    end if
    
  end function c_csr_get_nc_row



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

  subroutine  c_csr_free(a) 
    implicit none 

    class(psb_c_csr_sparse_mat), intent(inout) :: a

    if (allocated(a%irp)) deallocate(a%irp)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
    
    return

  end subroutine c_csr_free


end module psb_c_csr_mat_mod
