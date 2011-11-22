module psb_d_vect_mod

  use psb_d_base_vect_mod

  type psb_d_vect_type
    class(psb_d_base_vect_type), allocatable :: v 
  contains
    procedure, pass(x) :: get_nrows => d_vect_get_nrows
    procedure, pass(x) :: sizeof   => d_vect_sizeof
    procedure, pass(x) :: dot_v    => d_vect_dot_v
    procedure, pass(x) :: dot_a    => d_vect_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => d_vect_axpby_v
    procedure, pass(y) :: axpby_a  => d_vect_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    procedure, pass(y) :: mlt_v    => d_vect_mlt_v
    procedure, pass(y) :: mlt_a    => d_vect_mlt_a
    procedure, pass(z) :: mlt_a_2  => d_vect_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => d_vect_mlt_v_2
    procedure, pass(z) :: mlt_va   => d_vect_mlt_va
    procedure, pass(z) :: mlt_av   => d_vect_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
         & mlt_v_2, mlt_av, mlt_va
    procedure, pass(x) :: scal     => d_vect_scal
    procedure, pass(x) :: nrm2     => d_vect_nrm2
    procedure, pass(x) :: amax     => d_vect_amax
    procedure, pass(x) :: asum     => d_vect_asum
    procedure, pass(x) :: all      => d_vect_all
    procedure, pass(x) :: zero     => d_vect_zero
    procedure, pass(x) :: asb      => d_vect_asb
    procedure, pass(x) :: sync     => d_vect_sync
    procedure, pass(x) :: gthab    => d_vect_gthab
    procedure, pass(x) :: gthzv    => d_vect_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => d_vect_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => d_vect_free
    procedure, pass(x) :: ins      => d_vect_ins
    procedure, pass(x) :: bld_x    => d_vect_bld_x
    procedure, pass(x) :: bld_n    => d_vect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: getCopy  => d_vect_getCopy
    procedure, pass(x) :: cpy_vect => d_vect_cpy_vect
    generic, public    :: assignment(=) => cpy_vect
    procedure, pass(x) :: cnv      => d_vect_cnv
    procedure, pass(x) :: set_scal => d_vect_set_scal
    procedure, pass(x) :: set_vect => d_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
  end type psb_d_vect_type

  public  :: psb_d_vect
  private :: constructor, size_const
  interface psb_d_vect
    module procedure constructor, size_const
  end interface psb_d_vect

contains

  subroutine d_vect_bld_x(x,invect,mold)
    real(psb_dpk_), intent(in)          :: invect(:)
    class(psb_d_vect_type), intent(out) :: x
    class(psb_d_base_vect_type), intent(in), optional :: mold
    integer :: info

    if (present(mold)) then 
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_d_base_vect_type :: x%v,stat=info)
    endif

    if (info == psb_success_) call x%v%bld(invect)

  end subroutine d_vect_bld_x


  subroutine d_vect_bld_n(x,n,mold)
    integer, intent(in) :: n
    class(psb_d_vect_type), intent(out) :: x
    class(psb_d_base_vect_type), intent(in), optional :: mold
    integer :: info

    if (present(mold)) then 
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_d_base_vect_type :: x%v,stat=info)
    endif
    if (info == psb_success_) call x%v%bld(n)

  end subroutine d_vect_bld_n

  function  d_vect_getCopy(x) result(res)
    class(psb_d_vect_type), intent(in)  :: x
    real(psb_dpk_), allocatable         :: res(:)
    integer :: info

    if (allocated(x%v)) res = x%v%getCopy()

  end function d_vect_getCopy

  subroutine d_vect_cpy_vect(res,x)
    real(psb_dpk_), allocatable, intent(out) :: res(:)
    class(psb_d_vect_type), intent(in)  :: x
    integer :: info

    if (allocated(x%v)) res = x%v

  end subroutine d_vect_cpy_vect

  subroutine d_vect_set_scal(x,val)
    class(psb_d_vect_type), intent(inout)  :: x
    real(psb_dpk_), intent(in) :: val
        
    integer :: info
    if (allocated(x%v)) call x%v%set(val)
    
  end subroutine d_vect_set_scal

  subroutine d_vect_set_vect(x,val)
    class(psb_d_vect_type), intent(inout) :: x
    real(psb_dpk_), intent(in)            :: val(:)
        
    integer :: info
    if (allocated(x%v)) call x%v%set(val)
    
  end subroutine d_vect_set_vect


  function constructor(x) result(this)
    real(psb_dpk_)   :: x(:)
    type(psb_d_vect_type) :: this
    integer :: info

    allocate(psb_d_base_vect_type :: this%v, stat=info)

    if (info == 0) call this%v%bld(x)

    call this%asb(size(x),info)

  end function constructor


  function size_const(n) result(this)
    integer, intent(in) :: n
    type(psb_d_vect_type) :: this
    integer :: info

    allocate(psb_d_base_vect_type :: this%v, stat=info)
    call this%asb(n,info)

  end function size_const

  function d_vect_get_nrows(x) result(res)
    implicit none 
    class(psb_d_vect_type), intent(in) :: x
    integer :: res
    res = 0
    if (allocated(x%v)) res = x%v%get_nrows()
  end function d_vect_get_nrows

  function d_vect_sizeof(x) result(res)
    implicit none 
    class(psb_d_vect_type), intent(in) :: x
    integer(psb_long_int_k_) :: res
    res = 0
    if (allocated(x%v)) res = x%v%sizeof()
  end function d_vect_sizeof

  function d_vect_dot_v(n,x,y) result(res)
    implicit none 
    class(psb_d_vect_type), intent(inout) :: x, y
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res
    real(psb_dpk_), external      :: ddot

    res = dzero
    if (allocated(x%v).and.allocated(y%v)) &
         & res = x%v%dot(n,y%v)

  end function d_vect_dot_v

  function d_vect_dot_a(n,x,y) result(res)
    implicit none 
    class(psb_d_vect_type), intent(inout) :: x
    real(psb_dpk_), intent(in)    :: y(:)
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res
    real(psb_dpk_), external      :: ddot
    
    res = dzero
    if (allocated(x%v)) &
         & res = x%v%dot(n,y)
    
  end function d_vect_dot_a
    
  subroutine d_vect_axpby_v(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer, intent(in)               :: m
    class(psb_d_vect_type), intent(inout)  :: x
    class(psb_d_vect_type), intent(inout)  :: y
    real(psb_dpk_), intent (in)       :: alpha, beta
    integer, intent(out)              :: info
    
    if (allocated(x%v).and.allocated(y%v)) then 
      call y%v%axpby(m,alpha,x%v,beta,info)
    else
      info = psb_err_invalid_vect_state_
    end if

  end subroutine d_vect_axpby_v

  subroutine d_vect_axpby_a(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer, intent(in)               :: m
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_vect_type), intent(inout)  :: y
    real(psb_dpk_), intent (in)       :: alpha, beta
    integer, intent(out)              :: info
    
    if (allocated(y%v)) &
         & call y%v%axpby(m,alpha,x,beta,info)
    
  end subroutine d_vect_axpby_a

    
  subroutine d_vect_mlt_v(x, y, info)
    use psi_serial_mod
    implicit none 
    class(psb_d_vect_type), intent(inout)  :: x
    class(psb_d_vect_type), intent(inout)  :: y
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v)) &
         & call y%v%mlt(x%v,info)

  end subroutine d_vect_mlt_v

  subroutine d_vect_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_vect_type), intent(inout)  :: y
    integer, intent(out)              :: info
    integer :: i, n


    info = 0
    if (allocated(y%v)) &
         & call y%v%mlt(x,info)
    
  end subroutine d_vect_mlt_a


  subroutine d_vect_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: y(:)
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info
    integer :: i, n

    info = 0    
    if (allocated(z%v)) &
         & call z%v%mlt(alpha,x,y,beta,info)
    
  end subroutine d_vect_mlt_a_2

  subroutine d_vect_mlt_v_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: alpha,beta
    class(psb_d_vect_type), intent(inout)  :: x
    class(psb_d_vect_type), intent(inout)  :: y
    class(psb_d_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    if (allocated(x%v).and.allocated(y%v).and.&
         & allocated(z%v)) &
         & call z%v%mlt(alpha,x%v,y%v,beta,info)

  end subroutine d_vect_mlt_v_2

  subroutine d_vect_mlt_av(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_vect_type), intent(inout)  :: y
    class(psb_d_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    if (allocated(z%v).and.allocated(y%v)) &
         & call z%v%mlt(alpha,x,y%v,beta,info)

  end subroutine d_vect_mlt_av

  subroutine d_vect_mlt_va(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: y(:)
    class(psb_d_vect_type), intent(inout)  :: x
    class(psb_d_vect_type), intent(inout)  :: z
    integer, intent(out)              :: info    
    integer :: i, n

    info = 0
    
    if (allocated(z%v).and.allocated(x%v)) &
         & call z%v%mlt(alpha,x%v,y,beta,info)

  end subroutine d_vect_mlt_va

  subroutine d_vect_scal(alpha, x)
    use psi_serial_mod
    implicit none 
    class(psb_d_vect_type), intent(inout)  :: x
    real(psb_dpk_), intent (in)       :: alpha
    
    if (allocated(x%v)) call x%v%scal(alpha)

  end subroutine d_vect_scal


  function d_vect_nrm2(n,x) result(res)
    implicit none 
    class(psb_d_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res
    
    if (allocated(x%v)) then 
      res = x%v%nrm2(n)
    else
      res = dzero
    end if

  end function d_vect_nrm2
  
  function d_vect_amax(n,x) result(res)
    implicit none 
    class(psb_d_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res

    if (allocated(x%v)) then 
      res = x%v%amax(n)
    else
      res = dzero
    end if

  end function d_vect_amax

  function d_vect_asum(n,x) result(res)
    implicit none 
    class(psb_d_vect_type), intent(inout) :: x
    integer, intent(in)           :: n
    real(psb_dpk_)                :: res

    if (allocated(x%v)) then 
      res = x%v%asum(n)
    else
      res = dzero
    end if

  end function d_vect_asum
  
  subroutine d_vect_all(n, x, info, mold)

    implicit none 
    integer, intent(in)                 :: n
    class(psb_d_vect_type), intent(out) :: x
    class(psb_d_base_vect_type), intent(in), optional :: mold
    integer, intent(out)                :: info
    
    if (present(mold)) then 
      allocate(x%v,stat=info,mold=mold)
    else
      allocate(psb_d_base_vect_type :: x%v,stat=info)
    endif
    if (info == 0) then 
      call x%v%all(n,info)
    else
      info = psb_err_alloc_dealloc_
    end if

  end subroutine d_vect_all

  subroutine d_vect_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_d_vect_type), intent(inout)    :: x

    if (allocated(x%v)) call x%v%zero()

  end subroutine d_vect_zero

  subroutine d_vect_asb(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in)              :: n
    class(psb_d_vect_type), intent(inout) :: x
    integer, intent(out)             :: info

    if (allocated(x%v)) &
         & call x%v%asb(n,info)
    
  end subroutine d_vect_asb

  subroutine d_vect_sync(x)
    implicit none 
    class(psb_d_vect_type), intent(inout) :: x
    
    if (allocated(x%v)) &
         & call x%v%sync()
    
  end subroutine d_vect_sync

  subroutine d_vect_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer :: n, idx(:)
    real(psb_dpk_) :: alpha, beta, y(:)
    class(psb_d_vect_type) :: x
    
    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,alpha,beta,y)
    
  end subroutine d_vect_gthab

  subroutine d_vect_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer :: n, idx(:)
    real(psb_dpk_) ::  y(:)
    class(psb_d_vect_type) :: x

    if (allocated(x%v)) &
         &  call x%v%gth(n,idx,y)
    
  end subroutine d_vect_gthzv

  subroutine d_vect_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer :: n, idx(:)
    real(psb_dpk_) :: beta, x(:)
    class(psb_d_vect_type) :: y
    
    if (allocated(y%v)) &
         &  call y%v%sct(n,idx,x,beta)

  end subroutine d_vect_sctb

  subroutine d_vect_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_d_vect_type), intent(inout)  :: x
    integer, intent(out)              :: info
    
    info = 0
    if (allocated(x%v)) then 
      call x%v%free(info)
      if (info == 0) deallocate(x%v,stat=info)
    end if
        
  end subroutine d_vect_free

  subroutine d_vect_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_d_vect_type), intent(inout)  :: x
    integer, intent(in)               :: n, dupl
    integer, intent(in)               :: irl(:)
    real(psb_dpk_), intent(in)        :: val(:)
    integer, intent(out)              :: info

    integer :: i

    info = 0
    if (.not.allocated(x%v)) then 
      info = psb_err_invalid_vect_state_
      return
    end if
    
    call  x%v%ins(n,irl,val,dupl,info)
    
  end subroutine d_vect_ins


  subroutine d_vect_cnv(x,mold)
    class(psb_d_vect_type), intent(inout)    :: x
    class(psb_d_base_vect_type), intent(in)  :: mold
    class(psb_d_base_vect_type), allocatable :: tmp
    real(psb_dpk_), allocatable              :: invect(:)
    integer :: info

    allocate(tmp,stat=info,mold=mold)
    call x%v%sync()
    if (info == psb_success_) call tmp%bld(x%v%v)
    call x%v%free(info)
    call move_alloc(tmp,x%v)

  end subroutine d_vect_cnv

end module psb_d_vect_mod
