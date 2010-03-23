module psb_s_base_mat_mod
  
  use psb_base_mat_mod

  type, extends(psb_base_sparse_mat) :: psb_s_base_sparse_mat
  contains
    procedure, pass(a) :: s_csmv       => psb_s_base_csmv
    procedure, pass(a) :: s_csmm       => psb_s_base_csmm
    generic, public    :: csmm         => s_csmm, s_csmv
    procedure, pass(a) :: s_inner_cssv => psb_s_base_inner_cssv    
    procedure, pass(a) :: s_inner_cssm => psb_s_base_inner_cssm
    generic, public    :: inner_cssm   => s_inner_cssm, s_inner_cssv
    procedure, pass(a) :: s_cssv       => psb_s_base_cssv
    procedure, pass(a) :: s_cssm       => psb_s_base_cssm
    generic, public    :: cssm         => s_cssm, s_cssv
    procedure, pass(a) :: s_scals      => psb_s_base_scals
    procedure, pass(a) :: s_scal       => psb_s_base_scal
    generic, public    :: scal         => s_scals, s_scal 
    procedure, pass(a) :: csnmi        => psb_s_base_csnmi
    procedure, pass(a) :: get_diag     => psb_s_base_get_diag
    
    procedure, pass(a) :: csput       => psb_s_base_csput  
    procedure, pass(a) :: s_csgetrow  => psb_s_base_csgetrow
    procedure, pass(a) :: s_csgetblk  => psb_s_base_csgetblk
    generic, public    :: csget       => s_csgetrow, s_csgetblk 
    procedure, pass(a) :: csclip      => psb_s_base_csclip 
    procedure, pass(a) :: cp_to_coo   => psb_s_base_cp_to_coo   
    procedure, pass(a) :: cp_from_coo => psb_s_base_cp_from_coo 
    procedure, pass(a) :: cp_to_fmt   => psb_s_base_cp_to_fmt   
    procedure, pass(a) :: cp_from_fmt => psb_s_base_cp_from_fmt 
    procedure, pass(a) :: mv_to_coo   => psb_s_base_mv_to_coo   
    procedure, pass(a) :: mv_from_coo => psb_s_base_mv_from_coo 
    procedure, pass(a) :: mv_to_fmt   => psb_s_base_mv_to_fmt   
    procedure, pass(a) :: mv_from_fmt => psb_s_base_mv_from_fmt 
    procedure, pass(a) :: s_base_cp_from
    generic, public    :: cp_from => s_base_cp_from
    procedure, pass(a) :: s_base_mv_from
    generic, public    :: mv_from => s_base_mv_from
    
    procedure, pass(a) :: transp_1mat => psb_s_base_transp_1mat
    procedure, pass(a) :: transp_2mat => psb_s_base_transp_2mat
    procedure, pass(a) :: transc_1mat => psb_s_base_transc_1mat
    procedure, pass(a) :: transc_2mat => psb_s_base_transc_2mat
    
  end type psb_s_base_sparse_mat
  
  private :: s_base_cssv, s_base_cssm, s_base_cp_from, s_base_mv_from
  
  
  type, extends(psb_s_base_sparse_mat) :: psb_s_coo_sparse_mat
    
    integer              :: nnz
    integer, allocatable :: ia(:), ja(:)
    real(psb_spk_), allocatable :: val(:)
    
  contains
    
    procedure, pass(a) :: get_size     => s_coo_get_size
    procedure, pass(a) :: get_nzeros   => s_coo_get_nzeros
    procedure, pass(a) :: set_nzeros   => s_coo_set_nzeros
    procedure, pass(a) :: get_fmt      => s_coo_get_fmt
    procedure, pass(a) :: sizeof       => s_coo_sizeof
    procedure, pass(a) :: s_csmm       => psb_s_coo_csmm
    procedure, pass(a) :: s_csmv       => psb_s_coo_csmv
    procedure, pass(a) :: s_inner_cssm => psb_s_coo_cssm
    procedure, pass(a) :: s_inner_cssv => psb_s_coo_cssv
    procedure, pass(a) :: s_scals      => psb_s_coo_scals
    procedure, pass(a) :: s_scal       => psb_s_coo_scal
    procedure, pass(a) :: csnmi        => psb_s_coo_csnmi
    procedure, pass(a) :: reallocate_nz => psb_s_coo_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_s_coo_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_s_cp_coo_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_s_cp_coo_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_s_cp_coo_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_s_cp_coo_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_s_mv_coo_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_s_mv_coo_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_s_mv_coo_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_s_mv_coo_from_fmt
    procedure, pass(a) :: csput        => psb_s_coo_csput
    procedure, pass(a) :: get_diag     => psb_s_coo_get_diag
    procedure, pass(a) :: s_csgetrow   => psb_s_coo_csgetrow
    procedure, pass(a) :: csgetptn     => psb_s_coo_csgetptn
    procedure, pass(a) :: get_nz_row   => psb_s_coo_get_nz_row
    procedure, pass(a) :: reinit       => psb_s_coo_reinit
    procedure, pass(a) :: fix          => psb_s_fix_coo
    procedure, pass(a) :: trim         => psb_s_coo_trim
    procedure, pass(a) :: print        => psb_s_coo_print
    procedure, pass(a) :: free         => s_coo_free
    procedure, pass(a) :: psb_s_coo_cp_from
    generic, public    :: cp_from => psb_s_coo_cp_from
    procedure, pass(a) :: psb_s_coo_mv_from
    generic, public    :: mv_from => psb_s_coo_mv_from
    procedure, pass(a) :: transp_1mat => s_coo_transp_1mat
    procedure, pass(a) :: transc_1mat => s_coo_transc_1mat
    
  end type psb_s_coo_sparse_mat
  
  private :: s_coo_get_nzeros, s_coo_set_nzeros, &
       & s_coo_get_fmt,  s_coo_free, s_coo_sizeof, &
       & s_coo_transp_1mat, s_coo_transc_1mat
  
  
  
  !===================
  !
  ! BASE interfaces
  !
  !===================
  
  
  interface 
    subroutine psb_s_base_csmm(alpha,a,x,beta,y,info,trans)
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_base_csmm
  end interface
  
  interface 
    subroutine psb_s_base_csmv(alpha,a,x,beta,y,info,trans) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:)
      real(psb_spk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_base_csmv
  end interface
  
  interface 
    subroutine psb_s_base_inner_cssm(alpha,a,x,beta,y,info,trans) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_base_inner_cssm
  end interface
  
  interface 
    subroutine psb_s_base_inner_cssv(alpha,a,x,beta,y,info,trans) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:)
      real(psb_spk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_base_inner_cssv
  end interface
  
  interface 
    subroutine psb_s_base_cssm(alpha,a,x,beta,y,info,trans,scale,d)
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      real(psb_spk_), intent(in), optional :: d(:)
    end subroutine psb_s_base_cssm
  end interface
  
  interface 
    subroutine psb_s_base_cssv(alpha,a,x,beta,y,info,trans,scale,d)
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)    :: alpha, beta, x(:)
      real(psb_spk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      real(psb_spk_), intent(in), optional :: d(:)
    end subroutine psb_s_base_cssv
  end interface
  
  interface 
    subroutine psb_s_base_scals(d,a,info) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: d
      integer, intent(out)            :: info
    end subroutine psb_s_base_scals
  end interface
  
  interface 
    subroutine psb_s_base_scal(d,a,info) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_s_base_scal
  end interface
  
  interface 
    function psb_s_base_csnmi(a) result(res)
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_s_base_csnmi
  end interface
  
  interface 
    subroutine psb_s_base_get_diag(a,d,info) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)     :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_s_base_get_diag
  end interface
  
  interface 
    subroutine psb_s_base_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psb_s_base_csput
  end interface
  
  interface 
    subroutine psb_s_base_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_base_csgetrow
  end interface
  
  interface 
    subroutine psb_s_base_csgetblk(imin,imax,a,b,info,&
         & jmin,jmax,iren,append,rscale,cscale)
      import psb_s_base_sparse_mat, psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer, intent(in)                  :: imin,imax
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_base_csgetblk
  end interface
  
  
  interface 
    subroutine psb_s_base_csclip(a,b,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
      import psb_s_base_sparse_mat, psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      class(psb_s_coo_sparse_mat), intent(out) :: b
      integer,intent(out)                  :: info
      integer, intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_base_csclip
  end interface
  
  
  interface 
    subroutine psb_s_base_cp_to_coo(a,b,info) 
      import psb_s_base_sparse_mat, psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_base_cp_to_coo
  end interface
  
  interface 
    subroutine psb_s_base_cp_from_coo(a,b,info) 
      import psb_s_base_sparse_mat, psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(in) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_base_cp_from_coo
  end interface
  
  interface 
    subroutine psb_s_base_cp_to_fmt(a,b,info) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(in) :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_base_cp_to_fmt
  end interface
  
  interface 
    subroutine psb_s_base_cp_from_fmt(a,b,info) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(in) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_base_cp_from_fmt
  end interface
  
  interface 
    subroutine psb_s_base_mv_to_coo(a,b,info) 
      import psb_s_base_sparse_mat, psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_base_mv_to_coo
  end interface
  
  interface 
    subroutine psb_s_base_mv_from_coo(a,b,info) 
      import psb_s_base_sparse_mat, psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_base_mv_from_coo
  end interface
  
  interface 
    subroutine psb_s_base_mv_to_fmt(a,b,info) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_base_mv_to_fmt
  end interface
  
  interface 
    subroutine psb_s_base_mv_from_fmt(a,b,info) 
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_base_mv_from_fmt
  end interface
  
  interface 
    subroutine psb_s_base_transp_2mat(a,b)
      import psb_s_base_sparse_mat, psb_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(out) :: a
      class(psb_base_sparse_mat), intent(in)   :: b
    end subroutine psb_s_base_transp_2mat
  end interface
  
  interface  
    subroutine psb_s_base_transc_2mat(a,b)
      import psb_s_base_sparse_mat, psb_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(out) :: a
      class(psb_base_sparse_mat), intent(in)   :: b
    end subroutine psb_s_base_transc_2mat
  end interface
  
  interface 
    subroutine psb_s_base_transp_1mat(a)
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
    end subroutine psb_s_base_transp_1mat
  end interface
  
  interface 
    subroutine psb_s_base_transc_1mat(a)
      import psb_s_base_sparse_mat, psb_spk_
      class(psb_s_base_sparse_mat), intent(inout) :: a
    end subroutine psb_s_base_transc_1mat
  end interface
  
  
  
  
  !=================
  !
  ! COO interfaces
  !
  !=================
  
  interface
    subroutine  psb_s_coo_reallocate_nz(nz,a) 
      import psb_s_coo_sparse_mat
      integer, intent(in) :: nz
      class(psb_s_coo_sparse_mat), intent(inout) :: a
    end subroutine psb_s_coo_reallocate_nz
  end interface
  
  interface 
    subroutine psb_s_coo_reinit(a,clear)
      import psb_s_coo_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_s_coo_reinit
  end interface
  
  interface
    subroutine  psb_s_coo_trim(a)
      import psb_s_coo_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout) :: a
    end subroutine psb_s_coo_trim
  end interface
  
  interface
    subroutine  psb_s_coo_allocate_mnnz(m,n,a,nz) 
      import psb_s_coo_sparse_mat
      integer, intent(in) :: m,n
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      integer, intent(in), optional :: nz
    end subroutine psb_s_coo_allocate_mnnz
  end interface
  
  interface
    subroutine psb_s_coo_print(iout,a,iv,eirs,eics,head,ivr,ivc)
      import psb_s_coo_sparse_mat
      integer, intent(in)               :: iout
      class(psb_s_coo_sparse_mat), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_s_coo_print
  end interface
  
  
  interface 
    function  psb_s_coo_get_nz_row(idx,a) result(res)
      import psb_s_coo_sparse_mat
      class(psb_s_coo_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: idx
      integer                              :: res
    end function psb_s_coo_get_nz_row
  end interface
  
  
  interface 
    subroutine psb_s_fix_coo_inner(nzin,dupl,ia,ja,val,nzout,info,idir) 
      import psb_spk_
      integer, intent(in)           :: nzin,dupl
      integer, intent(inout)        :: ia(:), ja(:)
      real(psb_spk_), intent(inout) :: val(:)
      integer, intent(out)          :: nzout, info
      integer, intent(in), optional :: idir
    end subroutine psb_s_fix_coo_inner
  end interface

  interface 
    subroutine psb_s_fix_coo(a,info,idir) 
      import psb_s_coo_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      integer, intent(out)                :: info
      integer, intent(in), optional :: idir
    end subroutine psb_s_fix_coo
  end interface
  
  interface 
    subroutine psb_s_cp_coo_to_coo(a,b,info) 
      import psb_s_coo_sparse_mat
      class(psb_s_coo_sparse_mat), intent(in) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_s_cp_coo_to_coo
  end interface
  
  interface 
    subroutine psb_s_cp_coo_from_coo(a,b,info) 
      import psb_s_coo_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine psb_s_cp_coo_from_coo
  end interface
  
  interface 
    subroutine psb_s_cp_coo_to_fmt(a,b,info) 
      import psb_s_coo_sparse_mat, psb_s_base_sparse_mat
      class(psb_s_coo_sparse_mat), intent(in)   :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                       :: info
    end subroutine psb_s_cp_coo_to_fmt
  end interface
  
  interface 
    subroutine psb_s_cp_coo_from_fmt(a,b,info) 
      import psb_s_coo_sparse_mat, psb_s_base_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(in)   :: b
      integer, intent(out)                        :: info
    end subroutine psb_s_cp_coo_from_fmt
  end interface
  
  interface 
    subroutine psb_s_mv_coo_to_coo(a,b,info) 
      import psb_s_coo_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout)   :: b
      integer, intent(out)            :: info
    end subroutine psb_s_mv_coo_to_coo
  end interface
  
  interface 
    subroutine psb_s_mv_coo_from_coo(a,b,info) 
      import psb_s_coo_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine psb_s_mv_coo_from_coo
  end interface
  
  interface 
    subroutine psb_s_mv_coo_to_fmt(a,b,info) 
      import psb_s_coo_sparse_mat, psb_s_base_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(inout)  :: b
      integer, intent(out)                        :: info
    end subroutine psb_s_mv_coo_to_fmt
  end interface
  
  interface 
    subroutine psb_s_mv_coo_from_fmt(a,b,info) 
      import psb_s_coo_sparse_mat, psb_s_base_sparse_mat
      class(psb_s_coo_sparse_mat), intent(inout)  :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine psb_s_mv_coo_from_fmt
  end interface
  
  interface 
    subroutine psb_s_coo_cp_from(a,b)
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      type(psb_s_coo_sparse_mat), intent(in)   :: b
    end subroutine psb_s_coo_cp_from
  end interface
  
  interface 
    subroutine psb_s_coo_mv_from(a,b)
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(inout)  :: a
      type(psb_s_coo_sparse_mat), intent(inout) :: b
    end subroutine psb_s_coo_mv_from
  end interface
  
  
  interface 
    subroutine psb_s_coo_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psb_s_coo_csput
  end interface
  
  interface 
    subroutine psb_s_coo_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_coo_csgetptn
  end interface
  
  interface 
    subroutine psb_s_coo_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_s_coo_csgetrow
  end interface
  
  interface 
    subroutine psb_s_coo_cssv(alpha,a,x,beta,y,info,trans) 
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:)
      real(psb_spk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_s_coo_cssv
    subroutine psb_s_coo_cssm(alpha,a,x,beta,y,info,trans) 
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_s_coo_cssm
  end interface
  
  interface 
    subroutine psb_s_coo_csmv(alpha,a,x,beta,y,info,trans) 
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:)
      real(psb_spk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_s_coo_csmv
    subroutine psb_s_coo_csmm(alpha,a,x,beta,y,info,trans) 
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_s_coo_csmm
  end interface
  
  
  interface 
    function psb_s_coo_csnmi(a) result(res)
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_s_coo_csnmi
  end interface
  
  interface 
    subroutine psb_s_coo_get_diag(a,d,info) 
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(out)     :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_s_coo_get_diag
  end interface
  
  interface 
    subroutine psb_s_coo_scal(d,a,info) 
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_s_coo_scal
  end interface
  
  interface
    subroutine psb_s_coo_scals(d,a,info) 
      import psb_s_coo_sparse_mat, psb_spk_
      class(psb_s_coo_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: d
      integer, intent(out)            :: info
    end subroutine psb_s_coo_scals
  end interface
  
  
contains 
  
  
  subroutine s_base_mv_from(a,b)
    
    implicit none 
    
    class(psb_s_base_sparse_mat), intent(out)   :: a
    type(psb_s_base_sparse_mat), intent(inout) :: b
    
    
    ! No new things here, very easy
    call a%psb_base_sparse_mat%mv_from(b%psb_base_sparse_mat)
    
    return
    
  end subroutine s_base_mv_from
  
  subroutine s_base_cp_from(a,b)
    implicit none 
    
    class(psb_s_base_sparse_mat), intent(out) :: a
    type(psb_s_base_sparse_mat), intent(in)  :: b
    
    ! No new things here, very easy
    call a%psb_base_sparse_mat%cp_from(b%psb_base_sparse_mat)
    
    return
    
  end subroutine s_base_cp_from
  
 
  
  !====================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  !====================================
  
  
  
  function s_coo_sizeof(a) result(res)
    implicit none 
    class(psb_s_coo_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8 + 1
    res = res + psb_sizeof_sp  * size(a%val)
    res = res + psb_sizeof_int * size(a%ia)
    res = res + psb_sizeof_int * size(a%ja)
    
  end function s_coo_sizeof
  
  
  function s_coo_get_fmt(a) result(res)
    implicit none 
    class(psb_s_coo_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'COO'
  end function s_coo_get_fmt
  
  
  function s_coo_get_size(a) result(res)
    implicit none 
    class(psb_s_coo_sparse_mat), intent(in) :: a
    integer :: res
    res = -1
    
    if (allocated(a%ia)) res = size(a%ia)
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
  end function s_coo_get_size
  
  
  function s_coo_get_nzeros(a) result(res)
    implicit none 
    class(psb_s_coo_sparse_mat), intent(in) :: a
    integer :: res
    res  = a%nnz
  end function s_coo_get_nzeros
  
  
  !====================================
  !
  !
  !
  ! Setters 
  !
  !
  !
  !
  !
  !
  !====================================
  
  subroutine  s_coo_set_nzeros(nz,a)
    implicit none 
    integer, intent(in) :: nz
    class(psb_s_coo_sparse_mat), intent(inout) :: a
    
    a%nnz = nz
    
  end subroutine s_coo_set_nzeros
  
  !====================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  !====================================
  
  
  
  subroutine  s_coo_free(a) 
    implicit none 
    
    class(psb_s_coo_sparse_mat), intent(inout) :: a
    
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
    call a%set_nzeros(0)
    
    return
    
  end subroutine s_coo_free
  
  
  
  !====================================
  !
  !
  !
  ! Computational routines
  !
  !
  !
  !
  !
  !
  !====================================
  subroutine s_coo_transp_1mat(a)
    implicit none 
    
    class(psb_s_coo_sparse_mat), intent(inout) :: a
    
    integer, allocatable :: itemp(:) 
    integer :: info
    
    call a%psb_s_base_sparse_mat%psb_base_sparse_mat%transp()
    call move_alloc(a%ia,itemp)
    call move_alloc(a%ja,a%ia)
    call move_alloc(itemp,a%ja)
    
    call a%fix(info)
    
    return
    
  end subroutine s_coo_transp_1mat
  
  subroutine s_coo_transc_1mat(a)
    implicit none 
    
    class(psb_s_coo_sparse_mat), intent(inout) :: a
    
    call a%transp() 
  end subroutine s_coo_transc_1mat



end module psb_s_base_mat_mod



