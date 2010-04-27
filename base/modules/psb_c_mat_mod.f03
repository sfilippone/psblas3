module psb_c_mat_mod

  use psb_c_base_mat_mod
  use psb_c_csr_mat_mod, only : psb_c_csr_sparse_mat
  use psb_c_csc_mat_mod, only : psb_c_csc_sparse_mat

  type :: psb_c_sparse_mat

    class(psb_c_base_sparse_mat), allocatable  :: a 

  contains
    ! Getters
    procedure, pass(a) :: get_nrows => psb_c_get_nrows
    procedure, pass(a) :: get_ncols => psb_c_get_ncols
    procedure, pass(a) :: get_nzeros => psb_c_get_nzeros
    procedure, pass(a) :: get_nz_row => psb_c_get_nz_row
    procedure, pass(a) :: get_size => psb_c_get_size
    procedure, pass(a) :: get_state => psb_c_get_state
    procedure, pass(a) :: get_dupl => psb_c_get_dupl
    procedure, pass(a) :: is_null => psb_c_is_null
    procedure, pass(a) :: is_bld => psb_c_is_bld
    procedure, pass(a) :: is_upd => psb_c_is_upd
    procedure, pass(a) :: is_asb => psb_c_is_asb
    procedure, pass(a) :: is_sorted => psb_c_is_sorted
    procedure, pass(a) :: is_upper => psb_c_is_upper
    procedure, pass(a) :: is_lower => psb_c_is_lower
    procedure, pass(a) :: is_triangle => psb_c_is_triangle
    procedure, pass(a) :: is_unit => psb_c_is_unit
    procedure, pass(a) :: get_fmt => psb_c_get_fmt
    procedure, pass(a) :: sizeof => psb_c_sizeof

    ! Setters
    procedure, pass(a) :: set_nrows    => psb_c_set_nrows
    procedure, pass(a) :: set_ncols    => psb_c_set_ncols
    procedure, pass(a) :: set_dupl     => psb_c_set_dupl
    procedure, pass(a) :: set_state    => psb_c_set_state
    procedure, pass(a) :: set_null     => psb_c_set_null
    procedure, pass(a) :: set_bld      => psb_c_set_bld
    procedure, pass(a) :: set_upd      => psb_c_set_upd
    procedure, pass(a) :: set_asb      => psb_c_set_asb
    procedure, pass(a) :: set_sorted   => psb_c_set_sorted
    procedure, pass(a) :: set_upper    => psb_c_set_upper
    procedure, pass(a) :: set_lower    => psb_c_set_lower
    procedure, pass(a) :: set_triangle => psb_c_set_triangle
    procedure, pass(a) :: set_unit     => psb_c_set_unit

    ! Memory/data management 
    procedure, pass(a) :: csall => psb_c_csall
    procedure, pass(a) :: free => psb_c_free
    procedure, pass(a) :: trim => psb_c_trim
    procedure, pass(a) :: csput  => psb_c_csput 
    procedure, pass(a) :: c_csgetptn => psb_c_csgetptn
    procedure, pass(a) :: c_csgetrow => psb_c_csgetrow
    procedure, pass(a) :: c_csgetblk => psb_c_csgetblk
    generic, public    :: csget => c_csgetptn, c_csgetrow, c_csgetblk 
    procedure, pass(a) :: c_csclip   => psb_c_csclip
    procedure, pass(a) :: c_b_csclip => psb_c_b_csclip
    generic, public    :: csclip => c_b_csclip, c_csclip
    procedure, pass(a) :: c_clip_d_ip => psb_c_clip_d_ip
    procedure, pass(a) :: c_clip_d    => psb_c_clip_d
    generic, public    :: clip_diag => c_clip_d_ip, c_clip_d
    procedure, pass(a) :: reall => psb_c_reallocate_nz
    procedure, pass(a) :: get_neigh    => psb_c_get_neigh
    procedure, pass(a) :: c_cscnv      => psb_c_cscnv
    procedure, pass(a) :: c_cscnv_ip   => psb_c_cscnv_ip
    procedure, pass(a) :: c_cscnv_base => psb_c_cscnv_base
    generic, public    :: cscnv => c_cscnv, c_cscnv_ip, c_cscnv_base
    procedure, pass(a) :: reinit => psb_c_reinit
    procedure, pass(a) :: print  => psb_c_sparse_print
    procedure, pass(a) :: c_mv_from => psb_c_mv_from
    generic, public    :: mv_from => c_mv_from
    procedure, pass(a) :: c_mv_to => psb_c_mv_to
    generic, public    :: mv_to => c_mv_to
    procedure, pass(a) :: c_cp_from => psb_c_cp_from
    generic, public    :: cp_from => c_cp_from
    procedure, pass(a) :: c_cp_to => psb_c_cp_to
    generic, public    :: cp_to => c_cp_to
    procedure, pass(a) :: c_transp_1mat => psb_c_transp_1mat
    procedure, pass(a) :: c_transp_2mat => psb_c_transp_2mat
    generic, public    :: transp => c_transp_1mat, c_transp_2mat
    procedure, pass(a) :: c_transc_1mat => psb_c_transc_1mat
    procedure, pass(a) :: c_transc_2mat => psb_c_transc_2mat
    generic, public    :: transc => c_transc_1mat, c_transc_2mat

    
    
    ! Computational routines 
    procedure, pass(a) :: get_diag => psb_c_get_diag
    procedure, pass(a) :: csnmi    => psb_c_csnmi
    procedure, pass(a) :: c_csmv   => psb_c_csmv
    procedure, pass(a) :: c_csmm   => psb_c_csmm
    generic, public    :: csmm => c_csmm, c_csmv
    procedure, pass(a) :: c_scals  => psb_c_scals
    procedure, pass(a) :: c_scal   => psb_c_scal
    generic, public    :: scal  => c_scals, c_scal 
    procedure, pass(a) :: c_cssv   => psb_c_cssv
    procedure, pass(a) :: c_cssm   => psb_c_cssm
    generic, public    :: cssm => c_cssm, c_cssv

  end type psb_c_sparse_mat

  private :: psb_c_get_nrows, psb_c_get_ncols, psb_c_get_nzeros, psb_c_get_size, &
       & psb_c_get_state, psb_c_get_dupl, psb_c_is_null, psb_c_is_bld, psb_c_is_upd, &
       & psb_c_is_asb, psb_c_is_sorted, psb_c_is_upper, psb_c_is_lower, psb_c_is_triangle,&
       & psb_c_get_nz_row

  interface psb_sizeof
    module procedure psb_c_sizeof
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
    subroutine  psb_c_set_nrows(m,a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      integer, intent(in) :: m
    end subroutine psb_c_set_nrows
  end interface
  
  interface 
    subroutine psb_c_set_ncols(n,a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_c_set_ncols
  end interface
  
  interface 
    subroutine  psb_c_set_state(n,a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_c_set_state
  end interface
  
  interface 
    subroutine  psb_c_set_dupl(n,a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      integer, intent(in) :: n
    end subroutine psb_c_set_dupl
  end interface
  
  interface 
    subroutine psb_c_set_null(a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_set_null
  end interface
  
  interface 
    subroutine psb_c_set_bld(a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_set_bld
  end interface
  
  interface 
    subroutine psb_c_set_upd(a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_set_upd
  end interface
  
  interface 
    subroutine psb_c_set_asb(a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_set_asb
  end interface
  
  interface 
    subroutine psb_c_set_sorted(a,val) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_c_set_sorted
  end interface
  
  interface 
    subroutine psb_c_set_triangle(a,val) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_c_set_triangle
  end interface
  
  interface 
    subroutine psb_c_set_unit(a,val) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_c_set_unit
  end interface
  
  interface 
    subroutine psb_c_set_lower(a,val) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_c_set_lower
  end interface
  
  interface 
    subroutine psb_c_set_upper(a,val) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      logical, intent(in), optional :: val
    end subroutine psb_c_set_upper
  end interface
  
  
  interface 
    subroutine psb_c_sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
      import psb_c_sparse_mat
      integer, intent(in)               :: iout
      class(psb_c_sparse_mat), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_c_sparse_print
  end interface
  
  interface 
    subroutine psb_c_get_neigh(a,idx,neigh,n,info,lev)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(in) :: a   
      integer, intent(in)                :: idx 
      integer, intent(out)               :: n   
      integer, allocatable, intent(out)  :: neigh(:)
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: lev 
    end subroutine psb_c_get_neigh
  end interface
  
  interface 
    subroutine psb_c_csall(nr,nc,a,info,nz) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(out) :: a
      integer, intent(in)             :: nr,nc
      integer, intent(out)            :: info
      integer, intent(in), optional   :: nz
    end subroutine psb_c_csall
  end interface
  
  interface 
    subroutine psb_c_reallocate_nz(nz,a) 
      import psb_c_sparse_mat
      integer, intent(in) :: nz
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_reallocate_nz
  end interface
  
  interface 
    subroutine psb_c_free(a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_free
  end interface
  
  interface 
    subroutine psb_c_trim(a) 
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_trim
  end interface
  
  interface 
    subroutine psb_c_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psb_c_csput
  end interface
  
  interface 
    subroutine psb_c_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_csgetptn
  end interface
  
  interface 
    subroutine psb_c_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_csgetrow
  end interface
  
  interface 
    subroutine psb_c_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      class(psb_c_sparse_mat), intent(out) :: b
      integer, intent(in)                  :: imin,imax
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_csgetblk
  end interface
  
  interface 
    subroutine psb_c_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      class(psb_c_sparse_mat), intent(out) :: b
      integer,intent(out)                  :: info
      integer, intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_csclip
  end interface
  
  interface 
    subroutine psb_c_b_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
      import psb_c_sparse_mat, psb_spk_, psb_c_coo_sparse_mat
      class(psb_c_sparse_mat), intent(in) :: a
      type(psb_c_coo_sparse_mat), intent(out) :: b
      integer,intent(out)                  :: info
      integer, intent(in), optional        :: imin,imax,jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_c_b_csclip
  end interface
  
  interface 
    subroutine psb_c_cscnv(a,b,info,type,mold,upd,dupl)
      import psb_c_sparse_mat, psb_spk_, psb_c_base_sparse_mat
      class(psb_c_sparse_mat), intent(in)    :: a
      class(psb_c_sparse_mat), intent(out)   :: b
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: type
      class(psb_c_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_c_cscnv
  end interface
  

  interface 
    subroutine psb_c_cscnv_ip(a,iinfo,type,mold,dupl)
      import psb_c_sparse_mat, psb_spk_, psb_c_base_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      integer, intent(out)                   :: iinfo
      integer,optional, intent(in)           :: dupl
      character(len=*), optional, intent(in) :: type
      class(psb_c_base_sparse_mat), intent(in), optional :: mold
    end subroutine psb_c_cscnv_ip
  end interface
  

  interface 
    subroutine psb_c_cscnv_base(a,b,info,dupl)
      import psb_c_sparse_mat, psb_spk_, psb_c_base_sparse_mat
      class(psb_c_sparse_mat), intent(in)       :: a
      class(psb_c_base_sparse_mat), intent(out) :: b
      integer, intent(out)                   :: info
      integer,optional, intent(in)           :: dupl
    end subroutine psb_c_cscnv_base
  end interface
  
  interface 
    subroutine psb_c_clip_d(a,b,info)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(in) :: a
      class(psb_c_sparse_mat), intent(out) :: b
      integer,intent(out)                  :: info
    end subroutine psb_c_clip_d
  end interface
  
  interface 
    subroutine psb_c_clip_d_ip(a,info)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      integer,intent(out)                  :: info
    end subroutine psb_c_clip_d_ip
  end interface
  
  interface 
    subroutine psb_c_mv_from(a,b)
      import psb_c_sparse_mat, psb_spk_, psb_c_base_sparse_mat
      class(psb_c_sparse_mat), intent(out) :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
    end subroutine psb_c_mv_from
  end interface
  
  interface 
    subroutine psb_c_cp_from(a,b)
      import psb_c_sparse_mat, psb_spk_, psb_c_base_sparse_mat
      class(psb_c_sparse_mat), intent(out) :: a
      class(psb_c_base_sparse_mat), intent(inout), allocatable :: b
    end subroutine psb_c_cp_from
  end interface
  
  interface 
    subroutine psb_c_mv_to(a,b)
      import psb_c_sparse_mat, psb_spk_, psb_c_base_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(out) :: b
    end subroutine psb_c_mv_to
  end interface
  
  interface 
    subroutine psb_c_cp_to(a,b)
      import psb_c_sparse_mat, psb_spk_, psb_c_base_sparse_mat    
      class(psb_c_sparse_mat), intent(in) :: a
      class(psb_c_base_sparse_mat), intent(out) :: b
    end subroutine psb_c_cp_to
  end interface
  
  interface psb_move_alloc 
    subroutine psb_c_sparse_mat_move(a,b,info)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
      class(psb_c_sparse_mat), intent(out)   :: b
      integer, intent(out)                   :: info
    end subroutine psb_c_sparse_mat_move
  end interface
  

  interface psb_clone
    subroutine psb_c_sparse_mat_clone(a,b,info)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(in)  :: a
      class(psb_c_sparse_mat), intent(out) :: b
      integer, intent(out)                 :: info
    end subroutine psb_c_sparse_mat_clone
  end interface
  
  interface 
    subroutine psb_c_transp_1mat(a)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_transp_1mat
  end interface
  
  interface 
    subroutine psb_c_transp_2mat(a,b)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(out) :: a
      class(psb_c_sparse_mat), intent(in)  :: b
    end subroutine psb_c_transp_2mat
  end interface
  
  interface 
    subroutine psb_c_transc_1mat(a)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a
    end subroutine psb_c_transc_1mat
  end interface
  
  interface 
    subroutine psb_c_transc_2mat(a,b)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(out) :: a
      class(psb_c_sparse_mat), intent(in)  :: b
    end subroutine psb_c_transc_2mat
  end interface
  
  interface 
    subroutine psb_c_reinit(a,clear)
      import psb_c_sparse_mat
      class(psb_c_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_c_reinit
    
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
    subroutine psb_c_csmm(alpha,a,x,beta,y,info,trans) 
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_c_csmm
    subroutine psb_c_csmv(alpha,a,x,beta,y,info,trans) 
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans
    end subroutine psb_c_csmv
  end interface
  
  interface psb_cssm
    subroutine psb_c_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout) :: y(:,:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      complex(psb_spk_), intent(in), optional :: d(:)
    end subroutine psb_c_cssm
    subroutine psb_c_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)    :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout) :: y(:)
      integer, intent(out)            :: info
      character, optional, intent(in) :: trans, scale
      complex(psb_spk_), intent(in), optional :: d(:)
    end subroutine psb_c_cssv
  end interface
  
  interface 
    function psb_c_csnmi(a) result(res)
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function psb_c_csnmi
  end interface
  
  interface 
    subroutine psb_c_get_diag(a,d,info)
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(out)          :: d(:)
      integer, intent(out)                 :: info
    end subroutine psb_c_get_diag
  end interface
  
  interface psb_scal
    subroutine psb_c_scal(d,a,info)
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)              :: d(:)
      integer, intent(out)                    :: info
    end subroutine psb_c_scal
    subroutine psb_c_scals(d,a,info)
      import psb_c_sparse_mat, psb_spk_
      class(psb_c_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)              :: d
      integer, intent(out)                    :: info
    end subroutine psb_c_scals
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

  
  function psb_c_sizeof(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    
    res = 0
    if (allocated(a%a)) then 
      res = a%a%sizeof()
    end if
    
  end function psb_c_sizeof



  function psb_c_get_fmt(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    character(len=5) :: res

    if (allocated(a%a)) then 
      res = a%a%get_fmt()
    else
      res = 'NULL'
    end if

  end function psb_c_get_fmt


  function psb_c_get_dupl(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psb_invalid_
    end if
  end function psb_c_get_dupl


  function psb_c_get_state(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_state()
    else
      res = psb_spmat_null_
    end if
  end function psb_c_get_state

  function psb_c_get_nrows(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function psb_c_get_nrows

  function psb_c_get_ncols(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_ncols()
    else
      res = 0
    end if

  end function psb_c_get_ncols

  function psb_c_is_triangle(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_triangle()
    else
      res = .false.
    end if

  end function psb_c_is_triangle

  function psb_c_is_unit(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_unit()
    else
      res = .false.
    end if

  end function psb_c_is_unit

  function psb_c_is_upper(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upper()
    else
      res = .false.
    end if

  end function psb_c_is_upper

  function psb_c_is_lower(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = .not. a%a%is_upper()
    else
      res = .false.
    end if

  end function psb_c_is_lower

  function psb_c_is_null(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_null() 
    else
      res = .true.
    end if

  end function psb_c_is_null

  function psb_c_is_bld(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_bld()
    else
      res = .false.
    end if

  end function psb_c_is_bld

  function psb_c_is_upd(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upd()
    else
      res = .false.
    end if

  end function psb_c_is_upd

  function psb_c_is_asb(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_asb()
    else
      res = .false.
    end if

  end function psb_c_is_asb

  function psb_c_is_sorted(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_sorted()
    else
      res = .false.
    end if

  end function psb_c_is_sorted



  function psb_c_get_nzeros(a) result(res)
    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    integer :: res

    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_nzeros()
    end if

  end function psb_c_get_nzeros

  function psb_c_get_size(a) result(res)

    implicit none 
    class(psb_c_sparse_mat), intent(in) :: a
    integer :: res


    res = 0
    if (allocated(a%a)) then 
      res = a%a%get_size()
    end if

  end function psb_c_get_size


  function psb_c_get_nz_row(idx,a) result(res)
    implicit none 
    integer, intent(in)               :: idx
    class(psb_c_sparse_mat), intent(in) :: a
    integer :: res
    
    res = 0
    
    if (allocated(a%a)) res = a%a%get_nz_row(idx)

  end function psb_c_get_nz_row


end module psb_c_mat_mod
