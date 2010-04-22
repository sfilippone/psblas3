module psb_d_mat_mod

  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, only : psb_d_csr_sparse_mat
  use psb_d_csc_mat_mod, only : psb_d_csc_sparse_mat

  type :: psb_d_sparse_mat

    class(psb_d_base_sparse_mat), allocatable  :: a 

  contains
    ! Getters
    procedure, pass(a) :: get_nrows
    procedure, pass(a) :: get_ncols
    procedure, pass(a) :: get_nzeros
    procedure, pass(a) :: get_nz_row
    procedure, pass(a) :: get_size
    procedure, pass(a) :: get_state
    procedure, pass(a) :: get_dupl
    procedure, pass(a) :: is_null
    procedure, pass(a) :: is_bld
    procedure, pass(a) :: is_upd
    procedure, pass(a) :: is_asb
    procedure, pass(a) :: is_sorted
    procedure, pass(a) :: is_upper
    procedure, pass(a) :: is_lower
    procedure, pass(a) :: is_triangle
    procedure, pass(a) :: is_unit
    procedure, pass(a) :: get_fmt => psb_d_get_fmt
    procedure, pass(a) :: sizeof => psb_d_sizeof

    ! Setters
    procedure, pass(a) :: set_nrows    => psb_d_set_nrows
    procedure, pass(a) :: set_ncols    => psb_d_set_ncols
    procedure, pass(a) :: set_dupl     => psb_d_set_dupl
    procedure, pass(a) :: set_state    => psb_d_set_state
    procedure, pass(a) :: set_null     => psb_d_set_null
    procedure, pass(a) :: set_bld      => psb_d_set_bld
    procedure, pass(a) :: set_upd      => psb_d_set_upd
    procedure, pass(a) :: set_asb      => psb_d_set_asb
    procedure, pass(a) :: set_sorted   => psb_d_set_sorted
    procedure, pass(a) :: set_upper    => psb_d_set_upper
    procedure, pass(a) :: set_lower    => psb_d_set_lower
    procedure, pass(a) :: set_triangle => psb_d_set_triangle
    procedure, pass(a) :: set_unit     => psb_d_set_unit

    ! Memory/data management 
    procedure, pass(a) :: csall => psb_d_csall
    procedure, pass(a) :: free => psb_d_free
    procedure, pass(a) :: trim => psb_d_trim
    procedure, pass(a) :: csput  => psb_d_csput 
    procedure, pass(a) :: d_csgetptn => psb_d_csgetptn
    procedure, pass(a) :: d_csgetrow => psb_d_csgetrow
    procedure, pass(a) :: d_csgetblk => psb_d_csgetblk
    generic, public    :: csget => d_csgetptn, d_csgetrow, d_csgetblk 
    procedure, pass(a) :: d_csclip   => psb_d_csclip
    procedure, pass(a) :: d_b_csclip => psb_d_b_csclip
    generic, public    :: csclip => d_b_csclip, d_csclip
    procedure, pass(a) :: d_clip_d_ip => psb_d_clip_d_ip
    procedure, pass(a) :: d_clip_d    => psb_d_clip_d
    generic, public    :: clip_diag => d_clip_d_ip, d_clip_d
    procedure, pass(a) :: reall => psb_d_reallocate_nz
    procedure, pass(a) :: get_neigh    => psb_d_get_neigh
    procedure, pass(a) :: d_cscnv      => psb_d_cscnv
    procedure, pass(a) :: d_cscnv_ip   => psb_d_cscnv_ip
    procedure, pass(a) :: d_cscnv_base => psb_d_cscnv_base
    generic, public    :: cscnv => d_cscnv, d_cscnv_ip, d_cscnv_base
    procedure, pass(a) :: reinit => psb_d_reinit
    procedure, pass(a) :: print  => psb_d_sparse_print
    procedure, pass(a) :: d_mv_from => psb_d_mv_from
    generic, public    :: mv_from => d_mv_from
    procedure, pass(a) :: d_mv_to => psb_d_mv_to
    generic, public    :: mv_to => d_mv_to
    procedure, pass(a) :: d_cp_from => psb_d_cp_from
    generic, public    :: cp_from => d_cp_from
    procedure, pass(a) :: d_cp_to => psb_d_cp_to
    generic, public    :: cp_to => d_cp_to
    procedure, pass(a) :: d_transp_1mat => psb_d_transp_1mat
    procedure, pass(a) :: d_transp_2mat => psb_d_transp_2mat
    generic, public    :: transp => d_transp_1mat, d_transp_2mat
    procedure, pass(a) :: d_transc_1mat => psb_d_transc_1mat
    procedure, pass(a) :: d_transc_2mat => psb_d_transc_2mat
    generic, public    :: transc => d_transc_1mat, d_transc_2mat

    
    
    ! Computational routines 
    procedure, pass(a) :: get_diag => psb_d_get_diag
    procedure, pass(a) :: csnmi    => psb_d_csnmi
    procedure, pass(a) :: d_csmv   => psb_d_csmv
    procedure, pass(a) :: d_csmm   => psb_d_csmm
    generic, public    :: csmm => d_csmm, d_csmv
    procedure, pass(a) :: d_scals  => psb_d_scals
    procedure, pass(a) :: d_scal   => psb_d_scal
    generic, public    :: scal  => d_scals, d_scal 
    procedure, pass(a) :: d_cssv   => psb_d_cssv
    procedure, pass(a) :: d_cssm   => psb_d_cssm
    generic, public    :: cssm => d_cssm, d_cssv

  end type psb_d_sparse_mat

  private :: get_nrows, get_ncols, get_nzeros, get_size, &
       & get_state, get_dupl, is_null, is_bld, is_upd, &
       & is_asb, is_sorted, is_upper, is_lower, is_triangle

  interface psb_sizeof
    module procedure psb_d_sizeof
  end interface


  ! == ===================================
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
  ! == ===================================


  interface 
    subroutine  psb_d_set_nrows(m,a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(in) :: m
    end subroutine psb_d_set_nrows
  end interface
  
  interface 
    subroutine psb_d_set_ncols(n,a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_d_set_ncols
  end interface
  
  interface 
    subroutine  psb_d_set_state(n,a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_d_set_state
  end interface
  
  interface 
    subroutine  psb_d_set_dupl(n,a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_d_set_dupl
  end interface
  
  interface 
    subroutine psb_d_set_null(a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_set_null
  end interface
  
  interface 
    subroutine psb_d_set_bld(a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_set_bld
  end interface
  
  interface 
    subroutine psb_d_set_upd(a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_set_upd
  end interface
  
  interface 
    subroutine psb_d_set_asb(a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_set_asb
  end interface
  
  interface 
    subroutine psb_d_set_sorted(a,val) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_d_set_sorted
  end interface
  
  interface 
    subroutine psb_d_set_triangle(a,val) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_d_set_triangle
  end interface
  
  interface 
    subroutine psb_d_set_unit(a,val) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_d_set_unit
  end interface
  
  interface 
    subroutine psb_d_set_lower(a,val) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_d_set_lower
  end interface
  
  interface 
    subroutine psb_d_set_upper(a,val) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_d_set_upper
  end interface
  
  
  interface 
    subroutine psb_d_sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
      import psb_d_sparse_mat
      integer, intent(in)               :: iout
      class(psb_d_sparse_mat), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_d_sparse_print
  end interface
  
  interface 
    subroutine psb_d_get_neigh(a,idx,neigh,n,info,lev)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(in) :: a   
      integer, intent(in)                :: idx 
      integer, intent(out)               :: n   
      integer, allocatable, intent(out)  :: neigh(:)
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: lev 
    end subroutine psb_d_get_neigh
  end interface
  
  interface 
    subroutine psb_d_csall(nr,nc,a,info,nz) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(out) :: a
      integer, intent(in)             :: nr,nc
      integer, intent(out)            :: info
      integer, intent(in), optional   :: nz
    end subroutine psb_d_csall
  end interface
  
  interface 
    subroutine psb_d_reallocate_nz(nz,a) 
      import psb_d_sparse_mat
      integer, intent(in) :: nz
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_reallocate_nz
  end interface
  
  interface 
    subroutine psb_d_free(a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_free
  end interface
  
  interface 
    subroutine psb_d_trim(a) 
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_trim
  end interface
  
  interface 
    subroutine psb_d_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psb_d_csput
  end interface
  
  interface 
    subroutine psb_d_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_csgetptn
  end interface
  
  interface 
    subroutine psb_d_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_csgetrow
  end interface
  
  interface 
    subroutine psb_d_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      class(psb_d_sparse_mat), intent(out) :: b
      integer, intent(in)                  :: imin,imax
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_csgetblk
  end interface
  
  interface 
    subroutine psb_d_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      class(psb_d_sparse_mat), intent(out) :: b
      integer,intent(out)                  :: info
      integer, intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_csclip
  end interface
  
  interface 
    subroutine psb_d_b_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import psb_d_sparse_mat, psb_dpk_, psb_d_coo_sparse_mat
      class(psb_d_sparse_mat), intent(in) :: a
      type(psb_d_coo_sparse_mat), intent(out) :: b
      integer,intent(out)                  :: info
      integer, intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_b_csclip
  end interface
  
  interface 
    subroutine psb_d_cscnv(a,b,info,type,mold,upd,dupl)
      import psb_d_sparse_mat, psb_dpk_, psb_d_base_sparse_mat
      class(psb_d_sparse_mat), intent(in)    :: a
      class(psb_d_sparse_mat), intent(out)   :: b
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: type
      class(psb_d_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_d_cscnv
  end interface
  

  interface 
    subroutine psb_d_cscnv_ip(a,iinfo,type,mold,dupl)
      import psb_d_sparse_mat, psb_dpk_, psb_d_base_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      integer, intent(out)                   :: iinfo
      integer,optional, intent(in)           :: dupl
      character(len=*), optional, intent(in) :: type
      class(psb_d_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_d_cscnv_ip
  end interface
  

  interface 
    subroutine psb_d_cscnv_base(a,b,info,dupl)
      import psb_d_sparse_mat, psb_dpk_, psb_d_base_sparse_mat
      class(psb_d_sparse_mat), intent(in)       :: a
      class(psb_d_base_sparse_mat), intent(out) :: b
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl
    end subroutine psb_d_cscnv_base
  end interface
  
  interface 
    subroutine psb_d_clip_d(a,b,info)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(in) :: a
      class(psb_d_sparse_mat), intent(out) :: b
      integer,intent(out)                  :: info
    end subroutine psb_d_clip_d
  end interface
  
  interface 
    subroutine psb_d_clip_d_ip(a,info)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      integer,intent(out)                  :: info
    end subroutine psb_d_clip_d_ip
  end interface
  
  interface 
    subroutine psb_d_mv_from(a,b)
      import psb_d_sparse_mat, psb_dpk_, psb_d_base_sparse_mat
      class(psb_d_sparse_mat), intent(out) :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
    end subroutine psb_d_mv_from
  end interface
  
  interface 
    subroutine psb_d_cp_from(a,b)
      import psb_d_sparse_mat, psb_dpk_, psb_d_base_sparse_mat
      class(psb_d_sparse_mat), intent(out) :: a
      class(psb_d_base_sparse_mat), intent(inout), allocatable :: b
    end subroutine psb_d_cp_from
  end interface
  
  interface 
    subroutine psb_d_mv_to(a,b)
      import psb_d_sparse_mat, psb_dpk_, psb_d_base_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(out) :: b
    end subroutine psb_d_mv_to
  end interface
  
  interface 
    subroutine psb_d_cp_to(a,b)
      import psb_d_sparse_mat, psb_dpk_, psb_d_base_sparse_mat    
      class(psb_d_sparse_mat), intent(in) :: a
      class(psb_d_base_sparse_mat), intent(out) :: b
    end subroutine psb_d_cp_to
  end interface
  
  interface psb_move_alloc 
    subroutine psb_d_sparse_mat_move(a,b,info)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
      class(psb_d_sparse_mat), intent(out)   :: b
      integer, intent(out)                   :: info
    end subroutine psb_d_sparse_mat_move
  end interface
  

  interface psb_clone
    subroutine psb_d_sparse_mat_clone(a,b,info)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(in)  :: a
      class(psb_d_sparse_mat), intent(out) :: b
      integer, intent(out)                 :: info
    end subroutine psb_d_sparse_mat_clone
  end interface
  
  interface 
    subroutine psb_d_transp_1mat(a)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_transp_1mat
  end interface
  
  interface 
    subroutine psb_d_transp_2mat(a,b)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(out) :: a
      class(psb_d_sparse_mat), intent(in)  :: b
    end subroutine psb_d_transp_2mat
  end interface
  
  interface 
    subroutine psb_d_transc_1mat(a)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a
    end subroutine psb_d_transc_1mat
  end interface
  
  interface 
    subroutine psb_d_transc_2mat(a,b)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(out) :: a
      class(psb_d_sparse_mat), intent(in)  :: b
    end subroutine psb_d_transc_2mat
  end interface
  
  interface 
    subroutine psb_d_reinit(a,clear)
      import psb_d_sparse_mat
      class(psb_d_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_d_reinit
    
  end interface
  


  ! == ===================================
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
  ! == ===================================

  interface psb_csmm
    subroutine psb_d_csmm(alpha,a,x,beta,y,info,trans) 
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_d_csmm
    subroutine psb_d_csmv(alpha,a,x,beta,y,info,trans) 
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_d_csmv
  end interface
  
  interface psb_cssm
    subroutine psb_d_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      real(psb_dpk_), intent(in), optional :: d(:)
    end subroutine psb_d_cssm
    subroutine psb_d_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      real(psb_dpk_), intent(in), optional :: d(:)
    end subroutine psb_d_cssv
  end interface
  
  interface 
    function psb_d_csnmi(a) result(res)
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_d_csnmi
  end interface
  
  interface 
    subroutine psb_d_get_diag(a,d,info)
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)          :: d(:)
      integer, intent(out)                 :: info
    end subroutine psb_d_get_diag
  end interface
  
  interface psb_scal
    subroutine psb_d_scal(d,a,info)
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)              :: d(:)
      integer, intent(out)                    :: info
    end subroutine psb_d_scal
    subroutine psb_d_scals(d,a,info)
      import psb_d_sparse_mat, psb_dpk_
      class(psb_d_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)              :: d
      integer, intent(out)                    :: info
    end subroutine psb_d_scals
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

  
  function psb_d_sizeof(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    
    res = 0
    if (allocated(a%a)) then 
      res = a%a%sizeof()
    end if
    
  end function psb_d_sizeof



  function psb_d_get_fmt(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    character(len=5) :: res

    if (allocated(a%a)) then 
      res = a%a%get_fmt()
    else
      res = 'NULL'
    end if

  end function psb_d_get_fmt



  function get_dupl(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psb_invalid_
    end if
  end function get_dupl


  function get_state(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_state()
    else
      res = psb_spmat_null_
    end if
  end function get_state

  function get_nrows(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function get_nrows

  function get_ncols(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_ncols()
    else
      res = 0
    end if

  end function get_ncols

  function is_triangle(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_triangle()
    else
      res = .false.
    end if

  end function is_triangle

  function is_unit(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_unit()
    else
      res = .false.
    end if

  end function is_unit

  function is_upper(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upper()
    else
      res = .false.
    end if

  end function is_upper

  function is_lower(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = .not. a%a%is_upper()
    else
      res = .false.
    end if

  end function is_lower

  function is_null(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_null() 
    else
      res = .true.
    end if

  end function is_null

  function is_bld(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_bld()
    else
      res = .false.
    end if

  end function is_bld

  function is_upd(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upd()
    else
      res = .false.
    end if

  end function is_upd

  function is_asb(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_asb()
    else
      res = .false.
    end if

  end function is_asb

  function is_sorted(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_sorted()
    else
      res = .false.
    end if

  end function is_sorted



  function get_nzeros(a) result(res)
    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    integer :: res

    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_nzeros()
    end if

  end function get_nzeros

  function get_size(a) result(res)

    implicit none 
    class(psb_d_sparse_mat), intent(in) :: a
    integer :: res


    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_size()
    end if

  end function get_size


  function get_nz_row(idx,a) result(res)
    implicit none 
    integer, intent(in)               :: idx
    class(psb_d_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act

    res = 0
    
    if (allocated(a%a)) res = a%a%get_nz_row(idx)

  end function get_nz_row


end module psb_d_mat_mod
