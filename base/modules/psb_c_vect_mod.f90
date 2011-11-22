module psb_c_vect_mod

  use psb_c_base_vect_mod

  type psb_c_vect_type
    class(psb_c_base_vect_type), allocatable :: v 
  contains
    procedure, pass(x) :: get_nrows => c_vect_get_nrows
    procedure, pass(x) :: sizeof   => c_vect_sizeof
    procedure, pass(x) :: dot_v    => c_vect_dot_v
    procedure, pass(x) :: dot_a    => c_vect_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => c_vect_axpby_v
    procedure, pass(y) :: axpby_a  => c_vect_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    procedure, pass(y) :: mlt_v    => c_vect_mlt_v
    procedure, pass(y) :: mlt_a    => c_vect_mlt_a
    procedure, pass(z) :: mlt_a_2  => c_vect_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => c_vect_mlt_v_2
    procedure, pass(z) :: mlt_va   => c_vect_mlt_va
    procedure, pass(z) :: mlt_av   => c_vect_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
         & mlt_v_2, mlt_av, mlt_va
    procedure, pass(x) :: scal     => c_vect_scal
    procedure, pass(x) :: nrm2     => c_vect_nrm2
    procedure, pass(x) :: amax     => c_vect_amax
    procedure, pass(x) :: asum     => c_vect_asum
    procedure, pass(x) :: all      => c_vect_all
    procedure, pass(x) :: zero     => c_vect_zero
    procedure, pass(x) :: asb      => c_vect_asb
    procedure, pass(x) :: sync     => c_vect_sync
    procedure, pass(x) :: gthab    => c_vect_gthab
    procedure, pass(x) :: gthzv    => c_vect_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => c_vect_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => c_vect_free
    procedure, pass(x) :: ins      => c_vect_ins
    procedure, pass(x) :: bld_x    => c_vect_bld_x
    procedure, pass(x) :: bld_n    => c_vect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: getCopy  => c_vect_getCopy
    procedure, pass(x) :: cpy_vect => c_vect_cpy_vect
    generic, public    :: assignment(=) => cpy_vect
    procedure, pass(x) :: cnv      => c_vect_cnv
    procedure, pass(x) :: set_scal => c_vect_set_scal
    procedure, pass(x) :: set_vect => c_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
  end type psb_c_vect_type

  public  :: psb_c_vect
  private :: constructor, size_const
  interface psb_c_vect
    module procedure constructor, size_const
  end interface psb_c_vect

contains

  subroutine c_vect_bld_x(x,invect,mold)
    complex(psb_spk_), intent(in)          :: invect(:)
    class(psb_c_vect_type), intent(out) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    integer :: info

    if (present(mold)) then 
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_c_base_vect_type :: x%v,stat=info)
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine c_vect_bld_x


  subroutine c_vect_bld_n(x,n,mold)
    integer, intent(in) :: n
    class(psb_c_vect_type), intent(out) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    integer :: info

    if (present(mold)) then 
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_c_base_vect_type :: x%v,stat=info)
    endif
    if (info == psb_success_) call x%v%bld(n)

  end subroutine c_vect_bld_n

  function  c_vect_getCopy(x) result(res)
    class(psb_c_vect_type), intent(in)  :: x
    complex(psb_spk_), allocatable              :: res(:)
    integer :: info

    if (allocated(x%v)) res = x%v%getCopy()

  end function c_vect_getCopy

  subroutine c_vect_cpy_vect(res,x)
    complex(psb_spk_), allocatable, intent(out) :: res(:)
    class(psb_c_vect_type), intent(in)  :: x
    integer :: info

    if (allocated(x%v)) res = x%v

  end subroutine c_vect_cpy_vect

  subroutine c_vect_set_scal(x,val)
    class(psb_c_vect_type), intent(inout)  :: x
    complex(psb_spk_), intent(in) :: val
        
    integer :: info
    if (allocated(x%v)) call x%v%set(val)
    
  end subroutine c_vect_set_scal

  subroutine c_vect_set_vect(x,val)
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)         :: val(:)
        
    integer :: info
    if (allocated(x%v)) call x%v%set(val)
    
  end subroutine c_vect_set_vect


  function constructor(x) result(this)
    complex(psb_spk_)   :: x(:)
    type(psb_c_vect_type) :: this
    integer :: info

    allocate(psb_c_base_vect_type :: this%v, stat=info)

    if (info == 0) call this%v%bld(x)

    call this%asb(size(x),info)

  end function constructor


  function size_const(n) result(this)
    integer, intent(in) :: n
    type(psb_c_vect_type) :: this
    integer :: info

    allocate(psb_c_base_vect_type :: this%v, stat=info)
    call this%asb(n,info)

  end function size_const

  function c_vect_get_nrows(x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(in) :: x
    integer :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function c_vect_get_nrows

  function c_vect_sizeof(x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(in) :: x
    integer(psb_long_int_k_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function c_vect_sizeof

  function c_vect_dot_v(n,x,y) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x, y
    integer, intent(in)           :: n
    complex(psb_spk_)                :: res

    res = czero
    if (allocated(x%v).and.allocated(y%v)) &
         & res = x%v%dot(n,y%v)

  end function c_vect_dot_v

  function c_vect_dot_a(n,x,y) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    complex(psb_spk_), intent(in)    :: y(:)
    integer, intent(in)           :: n
    complex(psb_spk_)                :: res
    
    res = czero
    if (allocated(x%v)) &
         & res = x%v%dot(n,y)
    
  end function c_vect_dot_a
    
  subroutine c_vect_axpby_v(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer, intent(in)               :: m
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    complex(psb_spk_), intent (in)       :: alpha, beta
    integer, intent(out)              :: info
    
    if (allocated(x%v).and.allocated(y%v)) then 
      call y%v%axpby(m,alpha,x%v,beta,info)
    else
      info = psb_err_invalid_vect_state_
    end if

  end subroutine c_vect_axpby_v

  subroutine c_vect_axpby_a(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer, intent(in)               :: m
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: y
    complex(psb_spk_), intent (in)       :: alpha, beta
    integer, intent(out)              :: info
    
    if (allocated(y%v)) &
         & call y%v%axpby(m,alpha,x,beta,info)
    
  end subroutine c_vect_axpby_a

    
  subroutine c_vect_mlt_v(x, y, info)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%mlt(x%v,info)

  end subroutine c_vect_mlt_v

  subroutine c_vect_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: y
    integer, intent(out)              :: info
    integer :: i, n


    info = 0
    if (allocated(y%v)) &
         & call y%v%mlt(x,info)
    
  end subroutine c_vect_mlt_a


  subroutine c_vect_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)        :: alpha,beta
    complex(psb_spk_), intent(in)        :: y(:)
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info
    integer :: i, n

    info = 0    
    if (allocated(z%v)) &
         & call z%v%mlt(alpha,x,y,beta,info)
    
  end subroutine c_vect_mlt_a_2

  subroutine c_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)           :: alpha,beta
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    integer, intent(out)                   :: info    
    character(len=1), intent(in), optional :: conjgx, conjgy
    integer :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v).and.&
         & allocated(z%v)) &
         & call z%v%mlt(alpha,x%v,y%v,beta,info,conjgx,conjgy)

  end subroutine c_vect_mlt_v_2

  subroutine c_vect_mlt_av(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)        :: alpha,beta
    complex(psb_spk_), intent(in)        :: x(:)
    class(psb_c_vect_type), intent(inout)  :: y
    class(psb_c_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    if (allocated(z%v).and.allocated(y%v)) &
         & call z%v%mlt(alpha,x,y%v,beta,info)

  end subroutine c_vect_mlt_av

  subroutine c_vect_mlt_va(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    complex(psb_spk_), intent(in)        :: alpha,beta
    complex(psb_spk_), intent(in)        :: y(:)
    class(psb_c_vect_type), intent(inout)  :: x
    class(psb_c_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    
    if (allocated(z%v).and.allocated(x%v)) &
         & call z%v%mlt(alpha,x%v,y,beta,info)

  end subroutine c_vect_mlt_va

  subroutine c_vect_scal(alpha, x)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    complex(psb_spk_), intent (in)       :: alpha
    
    if (allocated(x%v)) call x%v%scal(alpha)

  end subroutine c_vect_scal


  function c_vect_nrm2(n,x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_spk_)                :: res
    
    if (allocated(x%v)) then 
      res = x%v%nrm2(n)
    else
      res = szero
    end if

  end function c_vect_nrm2
  
  function c_vect_amax(n,x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then 
      res = x%v%amax(n)
    else
      res = szero
    end if

  end function c_vect_amax

  function c_vect_asum(n,x) result(res)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_spk_)                :: res

    if (allocated(x%v)) then 
      res = x%v%asum(n)
    else
      res = szero
    end if

  end function c_vect_asum
  
  subroutine c_vect_all(n, x, info, mold)

    implicit none 
    integer, intent(in)                 :: n
    class(psb_c_vect_type), intent(out) :: x
    class(psb_c_base_vect_type), intent(in), optional :: mold
    integer, intent(out)                :: info
    
    if (present(mold)) then 
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_c_base_vect_type :: x%v,stat=info)
    endif
    if (info == 0) then 
      call x%v%all(n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine c_vect_all

  subroutine c_vect_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine c_vect_zero

  subroutine c_vect_asb(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in)              :: n
    class(psb_c_vect_type), intent(inout) :: x
    integer, intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(n,info)
    
  end subroutine c_vect_asb

  subroutine c_vect_sync(x)
    implicit none 
    class(psb_c_vect_type), intent(inout) :: x
    
    if (allocated(x%v)) &
         & call x%v%sync()
    
  end subroutine c_vect_sync

  subroutine c_vect_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer :: n, idx(:)
    complex(psb_spk_) :: alpha, beta, y(:)
    class(psb_c_vect_type) :: x
    
    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)
    
  end subroutine c_vect_gthab

  subroutine c_vect_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer :: n, idx(:)
    complex(psb_spk_) ::  y(:)
    class(psb_c_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)
    
  end subroutine c_vect_gthzv

  subroutine c_vect_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer :: n, idx(:)
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_vect_type) :: y
    
    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine c_vect_sctb

  subroutine c_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    integer, intent(out)              :: info
    
    info = 0
    if (allocated(x%v)) then 
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if
        
  end subroutine c_vect_free

  subroutine c_vect_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_type), intent(inout)  :: x
    integer, intent(in)               :: n, dupl
    integer, intent(in)               :: irl(:)
    complex(psb_spk_), intent(in)        :: val(:)
    integer, intent(out)              :: info

    integer :: i

    info = 0
    if (.not.allocated(x%v)) then 
      info = psb_err_invalid_vect_state_
      return
    end if
    
    call  x%v%ins(n,irl,val,dupl,info)
    
  end subroutine c_vect_ins


  subroutine c_vect_cnv(x,mold)
    class(psb_c_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(in) :: mold
    class(psb_c_base_vect_type), allocatable :: tmp
    complex(psb_spk_), allocatable          :: invect(:)
    integer :: info

    allocate(tmp,stat=info,mold=mold)
    call x%v%sync()
    if (info == psb_success_) call tmp%bld(x%v%v)
    call x%v%free(info)
    call move_alloc(tmp,x%v)

  end subroutine c_vect_cnv

end module psb_c_vect_mod
