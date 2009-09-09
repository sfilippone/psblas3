module psbn_d_base_mat_mod

  use psbn_base_mat_mod

  type, extends(psbn_base_sparse_mat) :: psbn_d_base_sparse_mat
  contains
    procedure, pass(a) :: d_base_csmv
    procedure, pass(a) :: d_base_csmm
    generic, public    :: psbn_csmm => d_base_csmm, d_base_csmv
    procedure, pass(a) :: d_base_cssv
    procedure, pass(a) :: d_base_cssm
    generic, public    :: psbn_cssm => d_base_cssm, d_base_cssv
    procedure, pass(a) :: csins
    procedure, pass(a) :: cp_to_coo
    procedure, pass(a) :: cp_from_coo
    procedure, pass(a) :: cp_to_fmt
    procedure, pass(a) :: cp_from_fmt
    procedure, pass(a) :: mv_to_coo
    procedure, pass(a) :: mv_from_coo
    procedure, pass(a) :: mv_to_fmt
    procedure, pass(a) :: mv_from_fmt
  end type psbn_d_base_sparse_mat
  private :: d_base_csmv, d_base_csmm, d_base_cssv, d_base_cssm,&
       & csins, cp_to_coo, cp_from_coo, cp_to_fmt, cp_from_fmt, &
       & mv_to_coo, mv_from_coo, mv_to_fmt, mv_from_fmt
  
  
  type, extends(psbn_d_base_sparse_mat) :: psbn_d_coo_sparse_mat
    
    integer              :: nnz
    integer, allocatable :: ia(:), ja(:)
    real(psb_dpk_), allocatable :: val(:)
    
  contains
    
    procedure, pass(a)  :: get_size => d_coo_get_size
    procedure, pass(a)  :: get_nzeros => d_coo_get_nzeros
    procedure, pass(a)  :: set_nzeros => d_coo_set_nzeros
    procedure, pass(a)  :: d_base_csmm => d_coo_csmm
    procedure, pass(a)  :: d_base_csmv => d_coo_csmv
    procedure, pass(a)  :: d_base_cssm => d_coo_cssm
    procedure, pass(a)  :: d_base_cssv => d_coo_cssv
    procedure, pass(a)  :: csins => d_coo_csins
    procedure, pass(a)  :: reallocate_nz => d_coo_reallocate_nz
    procedure, pass(a)  :: allocate_mnnz => d_coo_allocate_mnnz
    procedure, pass(a)  :: allocate_mn => d_coo_allocate_mn
    procedure, pass(a)  :: cp_to_coo   => d_cp_coo_to_coo
    procedure, pass(a)  :: cp_from_coo => d_cp_coo_from_coo
    procedure, pass(a)  :: cp_to_fmt   => d_cp_coo_to_fmt
    procedure, pass(a)  :: cp_from_fmt => d_cp_coo_from_fmt
    procedure, pass(a)  :: fix      => d_fix_coo
    procedure, pass(a)  :: free     => d_coo_free
    procedure, pass(a)  :: print    => d_coo_print
    
  end type psbn_d_coo_sparse_mat
  private :: d_coo_get_nzeros, d_coo_set_nzeros, &
       & d_coo_csmm, d_coo_csmv, d_coo_cssm, d_coo_cssv, &
       & d_coo_csins, d_coo_reallocate_nz, d_coo_allocate_mnnz, &
       & d_coo_allocate_mn, d_fix_coo, d_coo_free, &
       & d_cp_coo_to_coo, d_cp_coo_from_coo, &
       & d_cp_coo_to_fmt, d_cp_coo_from_fmt
  
  
  interface 
    subroutine d_fix_coo_inner(nzin,dupl,ia,ja,val,nzout,info,idir) 
      use psb_const_mod
      integer, intent(in)           :: nzin,dupl
      integer, intent(inout)        :: ia(:), ja(:)
      real(psb_dpk_), intent(inout) :: val(:)
      integer, intent(out)          :: nzout, info
      integer, intent(in), optional :: idir
    end subroutine d_fix_coo_inner
  end interface

  interface 
    subroutine d_fix_coo_impl(a,info,idir) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(inout) :: a
      integer, intent(out)                :: info
      integer, intent(in), optional :: idir
    end subroutine d_fix_coo_impl
  end interface

  interface 
    subroutine d_cp_coo_to_coo_impl(a,b,info) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(in) :: a
      class(psbn_d_coo_sparse_mat), intent(out) :: b
      integer, intent(out)            :: info
    end subroutine d_cp_coo_to_coo_impl
  end interface
  
  interface 
    subroutine d_cp_coo_from_coo_impl(a,b,info) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(inout) :: a
      class(psbn_d_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine d_cp_coo_from_coo_impl
  end interface

  interface 
    subroutine d_cp_coo_to_fmt_impl(a,b,info) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat, psbn_d_base_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(in)   :: a
      class(psbn_d_base_sparse_mat), intent(out) :: b
      integer, intent(out)                       :: info
    end subroutine d_cp_coo_to_fmt_impl
  end interface
  
  interface 
    subroutine d_cp_coo_from_fmt_impl(a,b,info) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat, psbn_d_base_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(inout) :: a
      class(psbn_d_base_sparse_mat), intent(in)   :: b
      integer, intent(out)                        :: info
    end subroutine d_cp_coo_from_fmt_impl
  end interface

  interface 
    subroutine d_mv_coo_to_coo_impl(a,b,info) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(inout) :: a
      class(psbn_d_coo_sparse_mat), intent(out)   :: b
      integer, intent(out)            :: info
    end subroutine d_mv_coo_to_coo_impl
  end interface
  
  interface 
    subroutine d_mv_coo_from_coo_impl(a,b,info) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(inout) :: a
      class(psbn_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine d_mv_coo_from_coo_impl
  end interface

  interface 
    subroutine d_mv_coo_to_fmt_impl(a,b,info) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat, psbn_d_base_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(inout) :: a
      class(psbn_d_base_sparse_mat), intent(out)  :: b
      integer, intent(out)                        :: info
    end subroutine d_mv_coo_to_fmt_impl
  end interface
  
  interface 
    subroutine d_mv_coo_from_fmt_impl(a,b,info) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat, psbn_d_base_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(inout)  :: a
      class(psbn_d_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine d_mv_coo_from_fmt_impl
  end interface

  
  interface 
    subroutine d_coo_csins_impl(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine d_coo_csins_impl
  end interface
  
  interface d_coo_cssm_impl
    subroutine d_coo_cssv_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine d_coo_cssv_impl
    subroutine d_coo_cssm_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine d_coo_cssm_impl
  end interface
  
  interface d_coo_csmm_impl
    subroutine d_coo_csmv_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine d_coo_csmv_impl
    subroutine d_coo_csmm_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psbn_d_coo_sparse_mat
      class(psbn_d_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine d_coo_csmm_impl
  end interface
  
contains 
  
  
  subroutine cp_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(in) :: a
    class(psbn_d_coo_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine cp_to_coo
  
  subroutine cp_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine cp_from_coo
  
  
  subroutine cp_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(in) :: a
    class(psbn_d_base_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='to_fmt'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine cp_to_fmt
  
  subroutine cp_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(inout) :: a
    class(psbn_d_base_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='from_fmt'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine cp_from_fmt
  
  
  subroutine mv_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine mv_to_coo
  
  subroutine mv_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine mv_from_coo
  
  
  subroutine mv_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(inout) :: a
    class(psbn_d_base_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='to_fmt'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine mv_to_fmt
  
  subroutine mv_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(inout) :: a
    class(psbn_d_base_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='from_fmt'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine mv_from_fmt
  
  
  
  subroutine d_fix_coo(a,info,idir) 
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    integer, intent(out)                :: info
    integer, intent(in), optional :: idir
    Integer :: err_act
    character(len=20)  :: name='fix_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_fix_coo_impl(a,info,idir)
    
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
    
  end subroutine d_fix_coo
  
  
  subroutine csins(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(inout) :: a
    real(psb_dpk_), intent(in)      :: val(:)
    integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
    integer, intent(out)            :: info
    integer, intent(in), optional   :: gtl(:)
    
    Integer :: err_act
    character(len=20)  :: name='csins'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine csins
  
  subroutine d_base_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
    real(psb_dpk_), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    
    Integer :: err_act
    character(len=20)  :: name='d_base_csmm'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_base_csmm
  
  subroutine d_base_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    
    Integer :: err_act
    character(len=20)  :: name='d_base_csmv'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
    
  end subroutine d_base_csmv
  
  subroutine d_base_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
    real(psb_dpk_), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    
    Integer :: err_act
    character(len=20)  :: name='d_base_cssm'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_base_cssm
  
  subroutine d_base_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psbn_d_base_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
    real(psb_dpk_), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    
    Integer :: err_act
    character(len=20)  :: name='d_base_cssv'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    info = 700
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
    
  end subroutine d_base_cssv
  
  
  subroutine d_cp_coo_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    class(psbn_d_coo_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_cp_coo_to_coo_impl(a,b,info)
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_cp_coo_to_coo
  
  subroutine d_cp_coo_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_cp_coo_from_coo_impl(a,b,info)
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_cp_coo_from_coo
  
  subroutine d_cp_coo_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    class(psbn_d_base_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_cp_coo_to_fmt_impl(a,b,info)
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_cp_coo_to_fmt
  
  subroutine d_cp_coo_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    class(psbn_d_base_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_cp_coo_from_fmt_impl(a,b,info)
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_cp_coo_from_fmt


  
  subroutine d_mv_coo_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_mv_coo_to_coo_impl(a,b,info)
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_mv_coo_to_coo
  
  subroutine d_mv_coo_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_mv_coo_from_coo_impl(a,b,info)
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_mv_coo_from_coo
  
  subroutine d_mv_coo_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    class(psbn_d_base_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_mv_coo_to_fmt_impl(a,b,info)
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_mv_coo_to_fmt
  
  subroutine d_mv_coo_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    class(psbn_d_base_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info
    
    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    call d_mv_coo_from_fmt_impl(a,b,info)
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    call psb_errpush(info,name)
    
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
    
  end subroutine d_mv_coo_from_fmt


  
  subroutine  d_coo_reallocate_nz(nz,a) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: nz
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='d_coo_reallocate_nz'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    
    call psb_realloc(nz,a%ia,a%ja,a%val,info)
    
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
    end if
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_coo_reallocate_nz
  
  
  function d_coo_get_size(a) result(res)
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(in) :: a
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
  end function d_coo_get_size
  
  
  function d_coo_get_nzeros(a) result(res)
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    integer :: res
    res  = a%nnz
  end function d_coo_get_nzeros
  
  
  subroutine  d_coo_set_nzeros(nz,a)
    implicit none 
    integer, intent(in) :: nz
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    
    a%nnz = nz
    
  end subroutine d_coo_set_nzeros
  
  
  subroutine d_coo_csins(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    real(psb_dpk_), intent(in)      :: val(:)
    integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
    integer, intent(out)            :: info
    integer, intent(in), optional   :: gtl(:)
    
    
    Integer            :: err_act
    character(len=20)  :: name='d_coo_csins'
    logical, parameter :: debug=.false.
    integer            :: nza, i,j,k, nzl, isza, int_err(5)
    
    call psb_erractionsave(err_act)
    info = 0
    
    if (nz <= 0) then 
      info = 10
      int_err(1)=1
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    if (size(ia) < nz) then 
      info = 35
      int_err(1)=2
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    
    if (size(ja) < nz) then 
      info = 35
      int_err(1)=3
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    if (size(val) < nz) then 
      info = 35
      int_err(1)=4
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    
    if (nz == 0) return
    nza = a%get_nzeros()
!!$    write(0,*) 'On entry to csins: ',nza
    call d_coo_csins_impl(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_coo_csins
  
  
  subroutine d_coo_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
    real(psb_dpk_), intent(inout)       :: y(:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc
    real(psb_dpk_) :: acc
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='d_coo_csmv'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    
    if (.not.a%is_asb()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif
    
    
    call d_coo_csmm_impl(alpha,a,x,beta,y,info,trans) 
    
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_coo_csmv
  
  subroutine d_coo_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
    real(psb_dpk_), intent(inout)       :: y(:,:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc, nc
    real(psb_dpk_), allocatable  :: acc(:)
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='d_coo_csmm'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    
    
    
    call d_coo_csmm_impl(alpha,a,x,beta,y,info,trans) 
    
    if (info /= 0) goto 9999
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_coo_csmm
  
  
  subroutine d_coo_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
    real(psb_dpk_), intent(inout)       :: y(:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc
    real(psb_dpk_) :: acc
    real(psb_dpk_), allocatable :: tmp(:)
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='d_coo_cssv'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    
    if (.not.a%is_asb()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif
    
    
    if (.not. (a%is_triangle())) then 
      write(0,*) 'Called SM on a non-triangular mat!'
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call d_coo_cssm_impl(alpha,a,x,beta,y,info,trans) 
    
    call psb_erractionrestore(err_act)
    return
    
    
9999 continue
    call psb_erractionrestore(err_act)
    
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
    
  end subroutine d_coo_cssv
  
  
  
  subroutine d_coo_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
    real(psb_dpk_), intent(inout)       :: y(:,:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc, nc
    real(psb_dpk_) :: acc
    real(psb_dpk_), allocatable :: tmp(:,:)
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='d_coo_csmm'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    
    if (.not.a%is_asb()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif
    
    
    if (.not. (a%is_triangle())) then 
      write(0,*) 'Called SM on a non-triangular mat!'
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call d_coo_cssm_impl(alpha,a,x,beta,y,info,trans) 
    call psb_erractionrestore(err_act)
    return
    
    
9999 continue
    call psb_erractionrestore(err_act)
    
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_coo_cssm
  
  
  subroutine  d_coo_free(a) 
    implicit none 
    
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
    
    return
    
  end subroutine d_coo_free
  
  subroutine  d_coo_allocate_mnnz(m,n,nz,a) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: m,n,nz
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='allocate_mnz'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    if (m < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/1,0,0,0,0/))
      goto 9999
    endif
    if (n < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/2,0,0,0,0/))
      goto 9999
    endif
    if (nz < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/3,0,0,0,0/))
      goto 9999
    endif
    
    if (info == 0) call psb_realloc(nz,a%ia,info)
    if (info == 0) call psb_realloc(nz,a%ja,info)
    if (info == 0) call psb_realloc(nz,a%val,info)
    if (info == 0) then 
      call a%set_nrows(m)
      call a%set_ncols(n)
      call a%set_nzeros(0)
      call a%set_bld()
      call a%set_triangle(.false.)
      call a%set_dupl(psbn_dupl_def_)
    end if
    if (info /= 0) goto 9999 
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_coo_allocate_mnnz
  
  
  subroutine  d_coo_allocate_mn(m,n,a) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: m,n
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    Integer :: err_act, info, nz
    character(len=20)  :: name='allocate_mn'
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0
    if (m < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/1,0,0,0,0/))
      goto 9999
    endif
    if (n < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/2,0,0,0,0/))
      goto 9999
    endif
    
    nz = max(7*m,7*n,1)
    call a%allocate(m,n,nz)
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_coo_allocate_mn


  subroutine d_coo_print(iout,a,iv,eirs,eics,head,ivr,ivc)
    use psb_spmat_type
    use psb_string_mod
    implicit none 

    integer, intent(in)               :: iout
    class(psbn_d_coo_sparse_mat), intent(in) :: a   
    integer, intent(in), optional     :: iv(:)
    integer, intent(in), optional     :: eirs,eics
    character(len=*), optional        :: head
    integer, intent(in), optional     :: ivr(:), ivc(:)

    Integer :: err_act
    character(len=20)  :: name='d_coo_print'
    logical, parameter :: debug=.false.

    character(len=80)                 :: frmtv 
    integer  :: irs,ics,i,j, nmx, ni, nr, nc, nz
    
    if (present(eirs)) then 
      irs = eirs
    else
      irs = 0
    endif
    if (present(eics)) then 
      ics = eics
    else
      ics = 0
    endif
    
    if (present(head)) then 
      write(iout,'(a)') '%%MatrixMarket matrix coordinate real general'
      write(iout,'(a,a)') '% ',head 
      write(iout,'(a)') '%'    
      write(iout,'(a,a)') '% COO'
    endif
    
    nr = a%get_nrows()
    nc = a%get_ncols()
    nz = a%get_nzeros()
    nmx = max(nr,nc,1)
    ni  = floor(log10(1.0*nmx)) + 1
    
    write(frmtv,'(a,i3.3,a,i3.3,a)') '(2(i',ni,',1x),es26.18,1x,2(i',ni,',1x))'
    write(iout,*) nr, nc, nz 
    if(present(iv)) then 
      do j=1,a%get_nzeros()
        write(iout,frmtv) iv(a%ia(j)),iv(a%ja(j)),a%val(j)
      enddo
    else      
      if (present(ivr).and..not.present(ivc)) then 
            do j=1,a%get_nzeros()
          write(iout,frmtv) ivr(a%ia(j)),a%ja(j),a%val(j)
        enddo
      else if (present(ivr).and.present(ivc)) then 
            do j=1,a%get_nzeros()
          write(iout,frmtv) ivr(a%ia(j)),ivc(a%ja(j)),a%val(j)
        enddo
      else if (.not.present(ivr).and.present(ivc)) then 
            do j=1,a%get_nzeros()
          write(iout,frmtv) a%ia(j),ivc(a%ja(j)),a%val(j)
        enddo
      else if (.not.present(ivr).and..not.present(ivc)) then 
            do j=1,a%get_nzeros()
          write(iout,frmtv) a%ia(j),a%ja(j),a%val(j)
        enddo
      endif
    endif

  end subroutine d_coo_print

  
  
end module psbn_d_base_mat_mod



