module psb_z_oacc_vect_mod
  use iso_c_binding
  use psb_const_mod
  use psb_error_mod
  use psb_z_vect_mod
  use psb_i_vect_mod
  use psb_i_oacc_vect_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_z_base_vect_type) :: psb_z_vect_oacc
    integer     :: state = is_host

  contains
    procedure, pass(x) :: get_nrows   => z_oacc_get_nrows
    procedure, nopass :: get_fmt      => z_oacc_get_fmt

    procedure, pass(x) :: all         => z_oacc_vect_all
    procedure, pass(x) :: zero        => z_oacc_zero
    procedure, pass(x) :: asb_m       => z_oacc_asb_m
    procedure, pass(x) :: sync        => z_oacc_sync
    procedure, pass(x) :: sync_space  => z_oacc_sync_space
    procedure, pass(x) :: bld_x       => z_oacc_bld_x
    procedure, pass(x) :: bld_mn      => z_oacc_bld_mn
    procedure, pass(x) :: free        => z_oacc_vect_free
    procedure, pass(x) :: ins_a       => z_oacc_ins_a
    procedure, pass(x) :: ins_v       => z_oacc_ins_v
    procedure, pass(x) :: is_host     => z_oacc_is_host
    procedure, pass(x) :: is_dev      => z_oacc_is_dev
    procedure, pass(x) :: is_sync     => z_oacc_is_sync
    procedure, pass(x) :: set_host    => z_oacc_set_host
    procedure, pass(x) :: set_dev     => z_oacc_set_dev
    procedure, pass(x) :: set_sync    => z_oacc_set_sync
    procedure, pass(x) :: set_scal    => z_oacc_set_scal

    procedure, pass(x) :: gthzv_x     => z_oacc_gthzv_x
    procedure, pass(x) :: gthzbuf_x   => z_oacc_gthzbuf
    procedure, pass(y) :: sctb        => z_oacc_sctb
    procedure, pass(y) :: sctb_x      => z_oacc_sctb_x
    procedure, pass(y) :: sctb_buf    => z_oacc_sctb_buf

    procedure, pass(x) :: get_size    => z_oacc_get_size
    procedure, pass(x) :: dot_v       => z_oacc_vect_dot
    procedure, pass(x) :: dot_a       => z_oacc_dot_a
    procedure, pass(y) :: axpby_v     => z_oacc_axpby_v
    procedure, pass(y) :: axpby_a     => z_oacc_axpby_a
    procedure, pass(z) :: abgdxyz     => z_oacc_abgdxyz
    procedure, pass(y) :: mlt_a       => z_oacc_mlt_a
    procedure, pass(z) :: mlt_a_2     => z_oacc_mlt_a_2
    procedure, pass(y) :: mlt_v       => z_oacc_mlt_v
    procedure, pass(z) :: mlt_v_2     => z_oacc_mlt_v_2
    procedure, pass(x) :: scal        => z_oacc_scal 
    procedure, pass(x) :: nrm2        => z_oacc_nrm2
    procedure, pass(x) :: amax        => z_oacc_amax
    procedure, pass(x) :: asum        => z_oacc_asum
    procedure, pass(x) :: absval1     => z_oacc_absval1
    procedure, pass(x) :: absval2     => z_oacc_absval2

  end type psb_z_vect_oacc

  interface
    subroutine z_oacc_mlt_v(x, y, info)
      import
      implicit none 
      class(psb_z_base_vect_type), intent(inout) :: x
      class(psb_z_vect_oacc), intent(inout)       :: y
      integer(psb_ipk_), intent(out)             :: info
    end subroutine z_oacc_mlt_v
  end interface
  

  interface
    subroutine z_oacc_mlt_v_2(alpha, x, y, beta, z, info, conjgx, conjgy)
      import
      implicit none 
      complex(psb_dpk_), intent(in)                 :: alpha, beta
      class(psb_z_base_vect_type), intent(inout) :: x
      class(psb_z_base_vect_type), intent(inout) :: y
      class(psb_z_vect_oacc), intent(inout)      :: z
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), intent(in), optional     :: conjgx, conjgy
    end subroutine z_oacc_mlt_v_2
  end interface
  
contains

  subroutine z_oacc_absval1(x)
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: n, i

    if (x%is_host()) call x%sync_space()
    n = size(x%v)
    !$acc parallel loop
    do i = 1, n
      x%v(i) = abs(x%v(i))
    end do
    call x%set_dev()
  end subroutine z_oacc_absval1

  subroutine z_oacc_absval2(x, y)
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    class(psb_z_base_vect_type), intent(inout) :: y
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: i

    n = min(size(x%v), size(y%v))
    select type (yy => y)
    class is (psb_z_vect_oacc)
        if (x%is_host()) call x%sync()
        if (yy%is_host()) call yy%sync()
        !$acc parallel loop
        do i = 1, n
            yy%v(i) = abs(x%v(i))
        end do
    class default
        if (x%is_dev()) call x%sync()
        if (y%is_dev()) call y%sync()
        call x%psb_z_base_vect_type%absval(y)
    end select
  end subroutine z_oacc_absval2

  subroutine z_oacc_scal(alpha, x)
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    complex(psb_dpk_), intent(in)          :: alpha
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: i

    if (x%is_host()) call x%sync_space()
    !$acc parallel loop
    do i = 1, size(x%v)
        x%v(i) = alpha * x%v(i)
    end do
    call x%set_dev()
  end subroutine z_oacc_scal

  function z_oacc_nrm2(n, x) result(res)
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_dpk_)                       :: res
    integer(psb_ipk_) :: info
    real(psb_dpk_) :: sum
    integer(psb_ipk_) :: i

    if (x%is_host()) call x%sync_space()
    sum = 0.0
    !$acc parallel loop reduction(+:sum)
    do i = 1, n
        sum = sum + abs(x%v(i))**2
    end do
    res = sqrt(sum)
  end function z_oacc_nrm2

  function z_oacc_amax(n, x) result(res)
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_dpk_)                       :: res
    integer(psb_ipk_) :: info
    real(psb_dpk_) :: max_val
    integer(psb_ipk_) :: i

    if (x%is_host()) call x%sync_space()
    max_val = -huge(0.0)
    !$acc parallel loop reduction(max:max_val)
    do i = 1, n
        if (abs(x%v(i)) > max_val) max_val = abs(x%v(i))
    end do
    res = max_val
  end function z_oacc_amax

  function z_oacc_asum(n, x) result(res)
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_dpk_)                       :: res
    integer(psb_ipk_) :: info
    complex(psb_dpk_) :: sum
    integer(psb_ipk_) :: i

    if (x%is_host()) call x%sync_space()
    sum = 0.0
    !$acc parallel loop reduction(+:sum)
    do i = 1, n
        sum = sum + abs(x%v(i))
    end do
    res = sum
  end function z_oacc_asum


  subroutine z_oacc_mlt_a(x, y, info)
    implicit none 
    complex(psb_dpk_), intent(in)           :: x(:)
    class(psb_z_vect_oacc), intent(inout) :: y
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: i, n
    
    info = 0    
    if (y%is_dev()) call y%sync_space()
    !$acc parallel loop
    do i = 1, size(x)
        y%v(i) = y%v(i) * x(i)
    end do
    call y%set_host()
  end subroutine z_oacc_mlt_a

  subroutine z_oacc_mlt_a_2(alpha, x, y, beta, z, info)
    implicit none 
    complex(psb_dpk_), intent(in)           :: alpha, beta
    complex(psb_dpk_), intent(in)           :: x(:)
    complex(psb_dpk_), intent(in)           :: y(:)
    class(psb_z_vect_oacc), intent(inout) :: z
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: i, n
    
    info = 0    
    if (z%is_dev()) call z%sync_space()
    !$acc parallel loop
    do i = 1, size(x)
        z%v(i) = alpha * x(i) * y(i) + beta * z%v(i)
    end do
    call z%set_host()
  end subroutine z_oacc_mlt_a_2


!!$  subroutine z_oacc_mlt_v(x, y, info)
!!$    implicit none 
!!$    class(psb_z_base_vect_type), intent(inout) :: x
!!$    class(psb_z_vect_oacc), intent(inout)       :: y
!!$    integer(psb_ipk_), intent(out)             :: info
!!$
!!$    integer(psb_ipk_) :: i, n
!!$    
!!$    info = 0    
!!$    n = min(x%get_nrows(), y%get_nrows())
!!$    select type(xx => x)
!!$    type is (psb_z_base_vect_type)
!!$        if (y%is_dev()) call y%sync()
!!$        !$acc parallel loop
!!$        do i = 1, n
!!$            y%v(i) = y%v(i) * xx%v(i)
!!$        end do
!!$        call y%set_host()
!!$    class default
!!$        if (xx%is_dev()) call xx%sync()
!!$        if (y%is_dev()) call y%sync()
!!$        !$acc parallel loop
!!$        do i = 1, n
!!$            y%v(i) = y%v(i) * xx%v(i)
!!$        end do
!!$        call y%set_host()
!!$    end select
!!$  end subroutine z_oacc_mlt_v
!!$
!!$  subroutine z_oacc_mlt_v_2(alpha, x, y, beta, z, info, conjgx, conjgy)
!!$    use psi_serial_mod
!!$    use psb_string_mod
!!$    implicit none 
!!$    complex(psb_dpk_), intent(in)                 :: alpha, beta
!!$    class(psb_z_base_vect_type), intent(inout) :: x
!!$    class(psb_z_base_vect_type), intent(inout) :: y
!!$    class(psb_z_vect_oacc), intent(inout)      :: z
!!$    integer(psb_ipk_), intent(out)             :: info
!!$    character(len=1), intent(in), optional     :: conjgx, conjgy
!!$    integer(psb_ipk_) :: i, n
!!$    logical :: conjgx_, conjgy_
!!$
!!$    conjgx_ = .false.
!!$    conjgy_ = .false.
!!$    if (present(conjgx)) conjgx_ = (psb_toupper(conjgx) == 'C')
!!$    if (present(conjgy)) conjgy_ = (psb_toupper(conjgy) == 'C')
!!$    
!!$    n = min(x%get_nrows(), y%get_nrows(), z%get_nrows())
!!$    
!!$    info = 0    
!!$    select type(xx => x)
!!$    class is (psb_z_vect_oacc)
!!$        select type (yy => y)
!!$        class is (psb_z_vect_oacc)
!!$            if (xx%is_host()) call xx%sync_space()
!!$            if (yy%is_host()) call yy%sync_space()
!!$            if ((beta /= zzero) .and. (z%is_host())) call z%sync_space()
!!$            !$acc parallel loop
!!$            do i = 1, n
!!$                z%v(i) = alpha * xx%v(i) * yy%v(i) + beta * z%v(i)
!!$            end do
!!$            call z%set_dev()
!!$        class default
!!$            if (xx%is_dev()) call xx%sync_space()
!!$            if (yy%is_dev()) call yy%sync()
!!$            if ((beta /= zzero) .and. (z%is_dev())) call z%sync_space()
!!$            !$acc parallel loop
!!$            do i = 1, n
!!$                z%v(i) = alpha * xx%v(i) * yy%v(i) + beta * z%v(i)
!!$            end do
!!$            call z%set_host()
!!$        end select
!!$    class default
!!$        if (x%is_dev()) call x%sync()
!!$        if (y%is_dev()) call y%sync()
!!$        if ((beta /= zzero) .and. (z%is_dev())) call z%sync_space()
!!$        !$acc parallel loop
!!$        do i = 1, n
!!$            z%v(i) = alpha * x%v(i) * y%v(i) + beta * z%v(i)
!!$        end do
!!$        call z%set_host()
!!$    end select
!!$  end subroutine z_oacc_mlt_v_2


  subroutine z_oacc_axpby_v(m, alpha, x, beta, y, info)
    !use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    class(psb_z_base_vect_type), intent(inout) :: x
    class(psb_z_vect_oacc), intent(inout) :: y
    complex(psb_dpk_), intent(in) :: alpha, beta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: nx, ny, i

    info = psb_success_

    select type(xx => x)
    type is (psb_z_vect_oacc)
        if ((beta /= zzero) .and. y%is_host()) call y%sync_space()
        if (xx%is_host()) call xx%sync_space()
        nx = size(xx%v)
        ny = size(y%v)
        if ((nx < m) .or. (ny < m)) then
            info = psb_err_internal_error_
        else
            !$acc parallel loop
            do i = 1, m
                y%v(i) = alpha * xx%v(i) + beta * y%v(i)
            end do
        end if
        call y%set_dev()
    class default
        if ((alpha /= zzero) .and. (x%is_dev())) call x%sync()
        call y%axpby(m, alpha, x%v, beta, info)
    end select
  end subroutine z_oacc_axpby_v

  subroutine z_oacc_axpby_a(m, alpha, x, beta, y, info)
    !use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    complex(psb_dpk_), intent(in) :: x(:)
    class(psb_z_vect_oacc), intent(inout) :: y
    complex(psb_dpk_), intent(in) :: alpha, beta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: i

    if ((beta /= zzero) .and. (y%is_dev())) call y%sync_space()
    !$acc parallel loop
    do i = 1, m
        y%v(i) = alpha * x(i) + beta * y%v(i)
    end do
    call y%set_host()
  end subroutine z_oacc_axpby_a

  subroutine z_oacc_abgdxyz(m, alpha, beta, gamma, delta, x, y, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    class(psb_z_base_vect_type), intent(inout) :: x
    class(psb_z_base_vect_type), intent(inout) :: y
    class(psb_z_vect_oacc), intent(inout) :: z
    complex(psb_dpk_), intent(in) :: alpha, beta, gamma, delta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: nx, ny, nz, i
    logical :: gpu_done

    info = psb_success_
    gpu_done = .false.

    select type(xx => x)
    class is (psb_z_vect_oacc)
        select type(yy => y)
        class is (psb_z_vect_oacc)
            select type(zz => z)
            class is (psb_z_vect_oacc)
                if ((beta /= zzero) .and. yy%is_host()) call yy%sync_space()
                if ((delta /= zzero) .and. zz%is_host()) call zz%sync_space()
                if (xx%is_host()) call xx%sync_space()
                nx = size(xx%v)
                ny = size(yy%v)
                nz = size(zz%v)
                if ((nx < m) .or. (ny < m) .or. (nz < m)) then
                    info = psb_err_internal_error_
                else
                    !$acc parallel loop
                    do i = 1, m
                        yy%v(i) = alpha * xx%v(i) + beta * yy%v(i)
                        zz%v(i) = gamma * yy%v(i) + delta * zz%v(i)
                    end do
                end if
                call yy%set_dev()
                call zz%set_dev()
                gpu_done = .true.
            end select
        end select
    end select

    if (.not. gpu_done) then
        if (x%is_host()) call x%sync()
        if (y%is_host()) call y%sync()
        if (z%is_host()) call z%sync()
        call y%axpby(m, alpha, x, beta, info)
        call z%axpby(m, gamma, y, delta, info)
    end if
  end subroutine z_oacc_abgdxyz

  subroutine z_oacc_sctb_buf(i, n, idx, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    complex(psb_dpk_) :: beta
    class(psb_z_vect_oacc) :: y
    integer(psb_ipk_) :: info

    if (.not.allocated(y%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_, 'sctb_buf')
      return
    end if

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync_space(info)
      if (y%is_host()) call y%sync_space()

      !$acc parallel loop
      do i = 1, n
        y%v(ii%v(i)) = beta * y%v(ii%v(i)) + y%combuf(i)
      end do

    class default
      !$acc parallel loop
      do i = 1, n
        y%v(idx%v(i)) = beta * y%v(idx%v(i)) + y%combuf(i)
      end do
    end select
  end subroutine z_oacc_sctb_buf

  subroutine z_oacc_sctb_x(i, n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_):: i, n
    class(psb_i_base_vect_type) :: idx
    complex(psb_dpk_) :: beta, x(:)
    class(psb_z_vect_oacc) :: y
    integer(psb_ipk_) :: info, ni

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync_space(info)
    class default
      call psb_errpush(info, 'z_oacc_sctb_x')
      return
    end select

    if (y%is_host()) call y%sync_space()

    !$acc parallel loop
    do i = 1, n
      y%v(idx%v(i)) = beta * y%v(idx%v(i)) + x(i)
    end do

    call y%set_dev()
  end subroutine z_oacc_sctb_x



  subroutine z_oacc_sctb(n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_dpk_) :: beta, x(:)
    class(psb_z_vect_oacc) :: y
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: i

    if (n == 0) return
    if (y%is_dev()) call y%sync_space()

    !$acc parallel loop
    do i = 1, n
      y%v(idx(i)) = beta * y%v(idx(i)) + x(i)
    end do

    call y%set_host()
  end subroutine z_oacc_sctb


  subroutine z_oacc_gthzbuf(i, n, idx, x)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    class(psb_z_vect_oacc) :: x
    integer(psb_ipk_) :: info

    info = 0
    if (.not.allocated(x%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_, 'gthzbuf')
      return
    end if

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync_space(info)
    class default
      call psb_errpush(info, 'z_oacc_gthzbuf')
      return
    end select

    if (x%is_host()) call x%sync_space()

    !$acc parallel loop
    do i = 1, n
      x%combuf(i) = x%v(idx%v(i))
    end do
  end subroutine z_oacc_gthzbuf

  subroutine z_oacc_gthzv_x(i, n, idx, x, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type):: idx
    complex(psb_dpk_) :: y(:)
    class(psb_z_vect_oacc):: x
    integer(psb_ipk_) :: info

    info = 0

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync_space(info)
    class default
      call psb_errpush(info, 'z_oacc_gthzv_x')
      return
    end select

    if (x%is_host()) call x%sync_space()

    !$acc parallel loop
    do i = 1, n
      y(i) = x%v(idx%v(i))
    end do
  end subroutine z_oacc_gthzv_x

  subroutine z_oacc_ins_v(n, irl, val, dupl, x, info)
    use psi_serial_mod
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n, dupl
    class(psb_i_base_vect_type), intent(inout) :: irl
    class(psb_z_base_vect_type), intent(inout) :: val
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, isz
    logical :: done_oacc

    info = 0
    if (psb_errstatus_fatal()) return

    done_oacc = .false.
    select type(virl => irl)
    type is (psb_i_vect_oacc)
      select type(vval => val)
      type is (psb_z_vect_oacc)
        if (vval%is_host()) call vval%sync_space()
        if (virl%is_host()) call virl%sync_space(info)
        if (x%is_host()) call x%sync_space()
        !$acc parallel loop
        do i = 1, n
          x%v(virl%v(i)) = vval%v(i)
        end do
        call x%set_dev()
        done_oacc = .true.
      end select
    end select

    if (.not.done_oacc) then
      select type(virl => irl)
      type is (psb_i_vect_oacc)
        if (virl%is_dev()) call virl%sync_space(info)
      end select
      select type(vval => val)
      type is (psb_z_vect_oacc)
        if (vval%is_dev()) call vval%sync_space()
      end select
      call x%ins(n, irl%v, val%v, dupl, info)
    end if

    if (info /= 0) then
      call psb_errpush(info, 'oacc_vect_ins')
      return
    end if

  end subroutine z_oacc_ins_v



  subroutine z_oacc_ins_a(n, irl, val, dupl, x, info)
    use psi_serial_mod
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n, dupl
    integer(psb_ipk_), intent(in) :: irl(:)
    complex(psb_dpk_), intent(in) :: val(:)
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i

    info = 0
    if (x%is_dev()) call x%sync_space()
    call x%psb_z_base_vect_type%ins(n, irl, val, dupl, info)
    call x%set_host()
    !$acc update device(x%v)

  end subroutine z_oacc_ins_a



  subroutine z_oacc_bld_mn(x, n)
    use psb_base_mod
    implicit none
    integer(psb_mpk_), intent(in) :: n
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%all(n, info)
    if (info /= 0) then
      call psb_errpush(info, 'z_oacc_bld_mn', i_err=(/n, n, n, n, n/))
    end if
    call x%set_host()
    !$acc update device(x%v)

  end subroutine z_oacc_bld_mn


  subroutine z_oacc_bld_x(x, this)
    use psb_base_mod
    implicit none
    complex(psb_dpk_), intent(in) :: this(:)
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this), x%v, info)
    if (info /= 0) then
      info = psb_err_alloc_request_
      call psb_errpush(info, 'z_oacc_bld_x', &
           i_err=(/size(this), izero, izero, izero, izero/))
      return
    end if

    x%v(:) = this(:)
    call x%set_host()
    !$acc update device(x%v)

  end subroutine z_oacc_bld_x


  subroutine z_oacc_asb_m(n, x, info)
    use psb_base_mod
    implicit none 
    integer(psb_mpk_), intent(in)        :: n
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_mpk_) :: nd

    info = psb_success_

    if (x%is_dev()) then
      nd = size(x%v)
      if (nd < n) then
        call x%sync()
        call x%psb_z_base_vect_type%asb(n, info)
        if (info == psb_success_) call x%sync_space()
        call x%set_host()
      end if
    else
      if (size(x%v) < n) then
        call x%psb_z_base_vect_type%asb(n, info)
        if (info == psb_success_) call x%sync_space()
        call x%set_host()
      end if
    end if
  end subroutine z_oacc_asb_m



  subroutine z_oacc_set_scal(x, val, first, last)
    class(psb_z_vect_oacc), intent(inout) :: x
    complex(psb_dpk_), intent(in)           :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: first_, last_
    first_ = 1
    last_  = x%get_nrows()
    if (present(first)) first_ = max(1, first)
    if (present(last))  last_  = min(last, last_)

    !$acc parallel loop
    do i = first_, last_
      x%v(i) = val
    end do
    !$acc end parallel loop

    call x%set_dev()
  end subroutine z_oacc_set_scal



  subroutine z_oacc_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x
    call x%set_dev()
    call x%set_scal(zzero)
  end subroutine z_oacc_zero

  function z_oacc_get_nrows(x) result(res)
    implicit none 
    class(psb_z_vect_oacc), intent(in) :: x
    integer(psb_ipk_) :: res

    if (allocated(x%v)) res = size(x%v)
  end function z_oacc_get_nrows

  function z_oacc_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = "zOACC"

  end function z_oacc_get_fmt

  function z_oacc_vect_dot(n, x, y) result(res)
    implicit none
    class(psb_z_vect_oacc), intent(inout) :: x
    class(psb_z_base_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(in) :: n
    complex(psb_dpk_) :: res
    complex(psb_dpk_), external :: ddot
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: i

    res = zzero

    select type(yy => y)
    type is (psb_z_base_vect_type)
        if (x%is_dev()) call x%sync()
        res = ddot(n, x%v, 1, yy%v, 1)
    type is (psb_z_vect_oacc)
        if (x%is_host()) call x%sync()
        if (yy%is_host()) call yy%sync()

        !$acc parallel loop reduction(+:res) present(x%v, yy%v)
        do i = 1, n
            res = res + x%v(i) * yy%v(i)
        end do
        !$acc end parallel loop

    class default
        call x%sync()
        res = y%dot(n, x%v)
    end select

  end function z_oacc_vect_dot




  function z_oacc_dot_a(n, x, y) result(res)
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x
    complex(psb_dpk_), intent(in) :: y(:)
    integer(psb_ipk_), intent(in) :: n
    complex(psb_dpk_)  :: res
    complex(psb_dpk_), external :: ddot

    if (x%is_dev()) call x%sync()
    res = ddot(n, y, 1, x%v, 1)

  end function z_oacc_dot_a

  ! subroutine z_oacc_set_vect(x,y)
  !   implicit none 
  !   class(psb_z_vect_oacc), intent(inout) :: x
  !   complex(psb_dpk_), intent(in) :: y(:)
  !   integer(psb_ipk_) :: info
  
  !   if (size(x%v) /= size(y)) then 
  !     call x%free(info)
  !     call x%all(size(y),info)
  !   end if
  !   x%v(:) = y(:)
  !   call x%set_host()
  ! end subroutine z_oacc_set_vect

  subroutine z_oacc_to_dev(v)
    implicit none
    complex(psb_dpk_) :: v(:)
    !$acc update device(v)          
  end subroutine z_oacc_to_dev

  subroutine z_oacc_to_host(v)
    implicit none
    complex(psb_dpk_) :: v(:)
    !$acc update self(v)          
  end subroutine z_oacc_to_host

  subroutine z_oacc_sync_space(x)
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x
    if (allocated(x%v)) then
      call z_oacc_create_dev(x%v)
    end if
  contains
    subroutine z_oacc_create_dev(v)
      implicit none
      complex(psb_dpk_) :: v(:)
      !$acc enter data copyin(v)          
    end subroutine z_oacc_create_dev
  end subroutine z_oacc_sync_space

  subroutine z_oacc_sync(x)
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x
    if (x%is_dev()) then
      call z_oacc_to_host(x%v)
    end if
    if (x%is_host()) then
      call z_oacc_to_dev(x%v)
    end if
    call x%set_sync()
  end subroutine z_oacc_sync

  subroutine z_oacc_set_host(x)
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x

    x%state = is_host
  end subroutine z_oacc_set_host

  subroutine z_oacc_set_dev(x)
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x

    x%state = is_dev
  end subroutine z_oacc_set_dev

  subroutine z_oacc_set_sync(x)
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x

    x%state = is_sync
  end subroutine z_oacc_set_sync

  function z_oacc_is_dev(x) result(res)
    implicit none 
    class(psb_z_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_dev)
  end function z_oacc_is_dev

  function z_oacc_is_host(x) result(res)
    implicit none 
    class(psb_z_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_host)
  end function z_oacc_is_host

  function z_oacc_is_sync(x) result(res)
    implicit none 
    class(psb_z_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_sync)
  end function z_oacc_is_sync

  subroutine z_oacc_vect_all(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)      :: n
    class(psb_z_vect_oacc), intent(out) :: x
    integer(psb_ipk_), intent(out)     :: info

    call psb_realloc(n, x%v, info)
    if (info == 0) then
      call x%set_host()
      !$acc enter data create(x%v)
      call x%sync_space()
    end if
    if (info /= 0) then 
      info = psb_err_alloc_request_
      call psb_errpush(info, 'z_oacc_all', &
           i_err=(/n, n, n, n, n/))
    end if
  end subroutine z_oacc_vect_all


  subroutine z_oacc_vect_free(x, info)
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)     :: info
    info = 0
    if (allocated(x%v)) then
      !$acc exit data delete(x%v) finalize
      deallocate(x%v, stat=info)
    end if

  end subroutine z_oacc_vect_free

  function z_oacc_get_size(x) result(res)
    implicit none 
    class(psb_z_vect_oacc), intent(inout) :: x
    integer(psb_ipk_)   :: res

    if (x%is_dev()) call x%sync()
    res = size(x%v)
  end function z_oacc_get_size

end module psb_z_oacc_vect_mod
