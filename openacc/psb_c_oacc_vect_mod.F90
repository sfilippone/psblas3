module psb_c_oacc_vect_mod
  use iso_c_binding
  use openacc
  use psb_const_mod
  use psb_error_mod
  use psb_c_vect_mod
  use psb_i_vect_mod
  use psb_i_oacc_vect_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_c_base_vect_type) :: psb_c_vect_oacc
    integer     :: state = is_host

  contains
    procedure, pass(x) :: get_nrows   => c_oacc_get_nrows
    procedure, nopass :: get_fmt      => c_oacc_get_fmt

    procedure, pass(x) :: all         => c_oacc_vect_all
    procedure, pass(x) :: zero        => c_oacc_zero
    procedure, pass(x) :: asb_m       => c_oacc_asb_m
    procedure, pass(x) :: sync        => c_oacc_sync
    procedure, pass(x) :: sync_space  => c_oacc_sync_space
    procedure, pass(x) :: bld_x       => c_oacc_bld_x
    procedure, pass(x) :: bld_mn      => c_oacc_bld_mn
    procedure, pass(x) :: free        => c_oacc_vect_free
    procedure, pass(x) :: ins_a       => c_oacc_ins_a
    procedure, pass(x) :: ins_v       => c_oacc_ins_v
    procedure, pass(x) :: is_host     => c_oacc_is_host
    procedure, pass(x) :: is_dev      => c_oacc_is_dev
    procedure, pass(x) :: is_sync     => c_oacc_is_sync
    procedure, pass(x) :: set_host    => c_oacc_set_host
    procedure, pass(x) :: set_dev     => c_oacc_set_dev
    procedure, pass(x) :: set_sync    => c_oacc_set_sync
    procedure, pass(x) :: set_scal    => c_oacc_set_scal

    procedure, pass(x) :: gthzv_x     => c_oacc_gthzv_x
    procedure, pass(x) :: gthzbuf_x   => c_oacc_gthzbuf
    procedure, pass(y) :: sctb        => c_oacc_sctb
    procedure, pass(y) :: sctb_x      => c_oacc_sctb_x
    procedure, pass(y) :: sctb_buf    => c_oacc_sctb_buf

    procedure, pass(x) :: get_size    => c_oacc_get_size

    procedure, pass(x) :: dot_v       => c_oacc_vect_dot
    procedure, pass(x) :: dot_a       => c_oacc_dot_a
    procedure, pass(y) :: axpby_v     => c_oacc_axpby_v
    procedure, pass(y) :: axpby_a     => c_oacc_axpby_a
    procedure, pass(z) :: upd_xyz     => c_oacc_upd_xyz
    procedure, pass(y) :: mlt_a       => c_oacc_mlt_a
    procedure, pass(z) :: mlt_a_2     => c_oacc_mlt_a_2
    procedure, pass(y) :: mlt_v       => psb_c_oacc_mlt_v
    procedure, pass(z) :: mlt_v_2     => psb_c_oacc_mlt_v_2
    procedure, pass(x) :: scal        => c_oacc_scal 
    procedure, pass(x) :: nrm2        => c_oacc_nrm2
    procedure, pass(x) :: amax        => c_oacc_amax
    procedure, pass(x) :: asum        => c_oacc_asum
    procedure, pass(x) :: absval1     => c_oacc_absval1
    procedure, pass(x) :: absval2     => c_oacc_absval2

  end type psb_c_vect_oacc

  interface
    subroutine psb_c_oacc_mlt_v(x, y, info)
      import
      implicit none 
      class(psb_c_base_vect_type), intent(inout) :: x
      class(psb_c_vect_oacc), intent(inout)       :: y
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_c_oacc_mlt_v
  end interface
  
  interface
    subroutine psb_c_oacc_mlt_v_2(alpha, x, y, beta, z, info, conjgx, conjgy)
      import
      implicit none 
      complex(psb_spk_), intent(in)                 :: alpha, beta
      class(psb_c_base_vect_type), intent(inout) :: x
      class(psb_c_base_vect_type), intent(inout) :: y
      class(psb_c_vect_oacc), intent(inout)      :: z
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), intent(in), optional     :: conjgx, conjgy
    end subroutine psb_c_oacc_mlt_v_2
  end interface

contains

  subroutine c_oacc_absval1(x)
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: n

    if (x%is_host()) call x%sync()
    n = size(x%v)
    call c_inner_oacc_absval1(n,x%v)
    call x%set_dev()
  contains 
    subroutine c_inner_oacc_absval1(n,x)
      implicit none
      complex(psb_spk_), intent(inout) :: x(:)
      integer(psb_ipk_) :: n
      integer(psb_ipk_) :: i
      !$acc parallel loop
      do i = 1, n
        x(i) = abs(x(i))
      end do
    end subroutine c_inner_oacc_absval1
  end subroutine c_oacc_absval1

  subroutine c_oacc_absval2(x, y)
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    class(psb_c_base_vect_type), intent(inout) :: y
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: i

    n = min(size(x%v), size(y%v))
    select type (yy => y)
    class is (psb_c_vect_oacc)
        if (x%is_host()) call x%sync()
        if (yy%is_host()) call yy%sync()
        call c_inner_oacc_absval2(n,x%v,yy%v)
    class default
        if (x%is_dev()) call x%sync()
        if (y%is_dev()) call y%sync()
        call x%psb_c_base_vect_type%absval(y)
    end select
  contains 
    subroutine c_inner_oacc_absval2(n,x,y)
      implicit none
      complex(psb_spk_), intent(inout) :: x(:),y(:)
      integer(psb_ipk_) :: n
      integer(psb_ipk_) :: i
      !$acc parallel loop
      do i = 1, n
        y(i) = abs(x(i))
      end do
    end subroutine c_inner_oacc_absval2
  end subroutine c_oacc_absval2

  subroutine c_oacc_scal(alpha, x)
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    complex(psb_spk_), intent(in)          :: alpha
    integer(psb_ipk_) :: info
    if (x%is_host()) call x%sync()
    call c_inner_oacc_scal(alpha, x%v)
    call x%set_dev()
  contains
    subroutine c_inner_oacc_scal(alpha, x)
      complex(psb_spk_), intent(in) :: alpha
      complex(psb_spk_), intent(inout) :: x(:)
      integer(psb_ipk_) :: i
      !$acc parallel loop
      do i = 1, size(x)
        x(i) = alpha * x(i)
      end do
    end subroutine c_inner_oacc_scal
  end subroutine c_oacc_scal

  function c_oacc_nrm2(n, x) result(res)
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                     :: res
    real(psb_spk_)                     :: mx
    integer(psb_ipk_) :: info

    if (x%is_host()) call x%sync()
    mx  = c_oacc_amax(n,x)
    res = c_inner_oacc_nrm2(n, mx, x%v)
  contains
    function c_inner_oacc_nrm2(n, mx,x) result(res)
      integer(psb_ipk_) :: n
      complex(psb_spk_) :: x(:)
      real(psb_spk_) :: mx, res
      real(psb_spk_) :: sum
      integer(psb_ipk_) :: i
      sum = 0.0
      !$acc parallel loop reduction(+:sum)
      do i = 1, n
        sum = sum + abs(x(i)/mx)**2
      end do
      res = mx*sqrt(sum)
    end function c_inner_oacc_nrm2
  end function c_oacc_nrm2

  function c_oacc_amax(n, x) result(res)
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                       :: res
    integer(psb_ipk_) :: info

    if (x%is_host()) call x%sync()
    res = c_inner_oacc_amax(n, x%v)
  contains
    function c_inner_oacc_amax(n, x) result(res)
      integer(psb_ipk_) :: n
      complex(psb_spk_) :: x(:)
      real(psb_spk_) :: res
      real(psb_spk_) :: max_val
      integer(psb_ipk_) :: i
      max_val = -huge(0.0)
      !$acc parallel loop reduction(max:max_val)
      do i = 1, n
        if (abs(x(i)) > max_val) max_val = abs(x(i))
      end do
      res = max_val
    end function c_inner_oacc_amax
  end function c_oacc_amax
  
  function c_oacc_asum(n, x) result(res)
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                       :: res
    integer(psb_ipk_) :: info
    complex(psb_spk_) :: sum
    integer(psb_ipk_) :: i
    if (x%is_host()) call x%sync()
    res =  c_inner_oacc_asum(n, x%v)
  contains
    function c_inner_oacc_asum(n, x) result(res)
      integer(psb_ipk_) :: n
      complex(psb_spk_) :: x(:)
      real(psb_spk_) :: res
      integer(psb_ipk_) :: i
      res = 0.0
      !$acc parallel loop reduction(+:res)
      do i = 1, n
        res = res + abs(x(i))
      end do
    end function c_inner_oacc_asum
  end function c_oacc_asum


  subroutine c_oacc_mlt_a(x, y, info)
    implicit none 
    complex(psb_spk_), intent(in)           :: x(:)
    class(psb_c_vect_oacc), intent(inout) :: y
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: i, n
    
    info = 0    
    if (y%is_dev()) call y%sync()
    !$acc parallel loop
    do i = 1, size(x)
        y%v(i) = y%v(i) * x(i)
    end do
    call y%set_host()
  end subroutine c_oacc_mlt_a

  subroutine c_oacc_mlt_a_2(alpha, x, y, beta, z, info)
    implicit none 
    complex(psb_spk_), intent(in)           :: alpha, beta
    complex(psb_spk_), intent(in)           :: x(:)
    complex(psb_spk_), intent(in)           :: y(:)
    class(psb_c_vect_oacc), intent(inout) :: z
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: i, n
    
    info = 0    
    if (z%is_dev()) call z%sync()
    !$acc parallel loop
    do i = 1, size(x)
        z%v(i) = alpha * x(i) * y(i) + beta * z%v(i)
    end do
    call z%set_host()
  end subroutine c_oacc_mlt_a_2


!!$  subroutine c_oacc_mlt_v(x, y, info)
!!$    implicit none 
!!$    class(psb_c_base_vect_type), intent(inout) :: x
!!$    class(psb_c_vect_oacc), intent(inout)       :: y
!!$    integer(psb_ipk_), intent(out)             :: info
!!$
!!$    integer(psb_ipk_) :: i, n
!!$    
!!$    info = 0    
!!$    n = min(x%get_nrows(), y%get_nrows())
!!$    select type(xx => x)
!!$    type is (psb_c_base_vect_type)
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
!!$  end subroutine c_oacc_mlt_v
!!$
!!$  subroutine c_oacc_mlt_v_2(alpha, x, y, beta, z, info, conjgx, conjgy)
!!$    use psi_serial_mod
!!$    use psb_string_mod
!!$    implicit none 
!!$    complex(psb_spk_), intent(in)                 :: alpha, beta
!!$    class(psb_c_base_vect_type), intent(inout) :: x
!!$    class(psb_c_base_vect_type), intent(inout) :: y
!!$    class(psb_c_vect_oacc), intent(inout)      :: z
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
!!$    class is (psb_c_vect_oacc)
!!$        select type (yy => y)
!!$        class is (psb_c_vect_oacc)
!!$            if (xx%is_host()) call xx%sync()
!!$            if (yy%is_host()) call yy%sync()
!!$            if ((beta /= czero) .and. (z%is_host())) call z%sync()
!!$            !$acc parallel loop
!!$            do i = 1, n
!!$                z%v(i) = alpha * xx%v(i) * yy%v(i) + beta * z%v(i)
!!$            end do
!!$            call z%set_dev()
!!$        class default
!!$            if (xx%is_dev()) call xx%sync()
!!$            if (yy%is_dev()) call yy%sync()
!!$            if ((beta /= czero) .and. (z%is_dev())) call z%sync()
!!$            !$acc parallel loop
!!$            do i = 1, n
!!$                z%v(i) = alpha * xx%v(i) * yy%v(i) + beta * z%v(i)
!!$            end do
!!$            call z%set_host()
!!$        end select
!!$    class default
!!$        if (x%is_dev()) call x%sync()
!!$        if (y%is_dev()) call y%sync()
!!$        if ((beta /= czero) .and. (z%is_dev())) call z%sync()
!!$        !$acc parallel loop
!!$        do i = 1, n
!!$            z%v(i) = alpha * x%v(i) * y%v(i) + beta * z%v(i)
!!$        end do
!!$        call z%set_host()
!!$    end select
!!$  end subroutine c_oacc_mlt_v_2


  subroutine c_oacc_axpby_v(m, alpha, x, beta, y, info)
    !use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    class(psb_c_base_vect_type), intent(inout) :: x
    class(psb_c_vect_oacc), intent(inout) :: y
    complex(psb_spk_), intent(in) :: alpha, beta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: nx, ny, i

    info = psb_success_

    select type(xx => x)
    type is (psb_c_vect_oacc)
        if ((beta /= czero) .and. y%is_host()) call y%sync()
        if (xx%is_host()) call xx%sync()
        nx = size(xx%v)
        ny = size(y%v)
        if ((nx < m) .or. (ny < m)) then
            info = psb_err_internal_error_
        else
          call c_inner_oacc_axpby(m, alpha, x%v, beta, y%v, info)
        end if
        call y%set_dev()
    class default
        if ((alpha /= czero) .and. (x%is_dev())) call x%sync()
        call y%axpby(m, alpha, x%v, beta, info)
      end select
    contains
      subroutine c_inner_oacc_axpby(m, alpha, x, beta, y, info)
      !use psi_serial_mod
      implicit none
      integer(psb_ipk_), intent(in) :: m
      complex(psb_spk_), intent(inout) :: x(:)
      complex(psb_spk_), intent(inout) :: y(:)
      complex(psb_spk_), intent(in) :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
      !$acc parallel 
      !$acc loop
      do i = 1, m
        y(i) = alpha * x(i) + beta * y(i)
      end do
      !$acc end parallel 
    end subroutine c_inner_oacc_axpby
  end subroutine c_oacc_axpby_v

  subroutine c_oacc_axpby_a(m, alpha, x, beta, y, info)
    !use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    complex(psb_spk_), intent(in) :: x(:)
    class(psb_c_vect_oacc), intent(inout) :: y
    complex(psb_spk_), intent(in) :: alpha, beta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: i

    if ((beta /= czero) .and. (y%is_dev())) call y%sync()
    !$acc parallel loop
    do i = 1, m
        y%v(i) = alpha * x(i) + beta * y%v(i)
    end do
    call y%set_host()
  end subroutine c_oacc_axpby_a

  subroutine c_oacc_upd_xyz(m, alpha, beta, gamma, delta, x, y, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    class(psb_c_base_vect_type), intent(inout) :: x
    class(psb_c_base_vect_type), intent(inout) :: y
    class(psb_c_vect_oacc), intent(inout) :: z
    complex(psb_spk_), intent(in) :: alpha, beta, gamma, delta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: nx, ny, nz, i
    logical :: gpu_done
    write(0,*)'upd_xyz'    
    info = psb_success_
    gpu_done = .false.

    select type(xx => x)
    class is (psb_c_vect_oacc)
        select type(yy => y)
        class is (psb_c_vect_oacc)
            select type(zz => z)
            class is (psb_c_vect_oacc)
                if ((beta /= czero) .and. yy%is_host()) call yy%sync()
                if ((delta /= czero) .and. zz%is_host()) call zz%sync()
                if (xx%is_host()) call xx%sync()
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
  end subroutine c_oacc_upd_xyz

  subroutine c_oacc_sctb_buf(i, n, idx, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    complex(psb_spk_) :: beta
    class(psb_c_vect_oacc) :: y
    integer(psb_ipk_) :: info

    if (.not.allocated(y%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_, 'sctb_buf')
      return
    end if

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
      if (y%is_host()) call y%sync()

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
  end subroutine c_oacc_sctb_buf

  subroutine c_oacc_sctb_x(i, n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_):: i, n
    class(psb_i_base_vect_type) :: idx
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_vect_oacc) :: y
    integer(psb_ipk_) :: info, ni

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
    class default
      call psb_errpush(info, 'c_oacc_sctb_x')
      return
    end select

    if (y%is_host()) call y%sync()

    !$acc parallel loop
    do i = 1, n
      y%v(idx%v(i)) = beta * y%v(idx%v(i)) + x(i)
    end do

    call y%set_dev()
  end subroutine c_oacc_sctb_x

  subroutine c_oacc_sctb(n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: idx(:)
    complex(psb_spk_) :: beta, x(:)
    class(psb_c_vect_oacc) :: y
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: i

    if (n == 0) return
    if (y%is_dev()) call y%sync()

    !$acc parallel loop
    do i = 1, n
      y%v(idx(i)) = beta * y%v(idx(i)) + x(i)
    end do

    call y%set_host()
  end subroutine c_oacc_sctb

  subroutine c_oacc_gthzbuf(i, n, idx, x)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    class(psb_c_vect_oacc) :: x
    integer(psb_ipk_) :: info

    info = 0
    if (.not.allocated(x%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_, 'gthzbuf')
      return
    end if

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
    class default
      call psb_errpush(info, 'c_oacc_gthzbuf')
      return
    end select

    if (x%is_host()) call x%sync()

    !$acc parallel loop
    do i = 1, n
      x%combuf(i) = x%v(idx%v(i))
    end do
  end subroutine c_oacc_gthzbuf

  subroutine c_oacc_gthzv_x(i, n, idx, x, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type):: idx
    complex(psb_spk_) :: y(:)
    class(psb_c_vect_oacc):: x
    integer(psb_ipk_) :: info

    info = 0

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
    class default
      call psb_errpush(info, 'c_oacc_gthzv_x')
      return
    end select

    if (x%is_host()) call x%sync()

    !$acc parallel loop
    do i = 1, n
      y(i) = x%v(idx%v(i))
    end do
  end subroutine c_oacc_gthzv_x

  subroutine c_oacc_ins_v(n, irl, val, dupl, x, info)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n, dupl
    class(psb_i_base_vect_type), intent(inout) :: irl
    class(psb_c_base_vect_type), intent(inout) :: val
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, isz
    logical :: done_oacc

    info = 0
    if (psb_errstatus_fatal()) return

    done_oacc = .false.
    select type(virl => irl)
    type is (psb_i_vect_oacc)
      select type(vval => val)
      type is (psb_c_vect_oacc)
        if (vval%is_host()) call vval%sync()
        if (virl%is_host()) call virl%sync()
        if (x%is_host()) call x%sync()
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
        if (virl%is_dev()) call virl%sync()
      end select
      select type(vval => val)
      type is (psb_c_vect_oacc)
        if (vval%is_dev()) call vval%sync()
      end select
      call x%ins(n, irl%v, val%v, dupl, info)
    end if

    if (info /= 0) then
      call psb_errpush(info, 'oacc_vect_ins')
      return
    end if

  end subroutine c_oacc_ins_v

  subroutine c_oacc_ins_a(n, irl, val, dupl, x, info)
    use psi_serial_mod
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n, dupl
    integer(psb_ipk_), intent(in) :: irl(:)
    complex(psb_spk_), intent(in) :: val(:)
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i

    info = 0
    if (x%is_dev()) call x%sync()
    call x%psb_c_base_vect_type%ins(n, irl, val, dupl, info)
    call x%set_host()
    !$acc update device(x%v)

  end subroutine c_oacc_ins_a

  subroutine c_oacc_bld_mn(x, n)
    use psb_base_mod
    implicit none
    integer(psb_mpk_), intent(in) :: n
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%all(n, info)
    if (info /= 0) then
      call psb_errpush(info, 'c_oacc_bld_mn', i_err=(/n, n, n, n, n/))
    end if
    call x%set_host()
    if (acc_is_present(x%v))  then
      !$acc exit data delete(x%v) finalize
    end if
    !$acc enter data copyin(x%v)

  end subroutine c_oacc_bld_mn


  subroutine c_oacc_bld_x(x, this)
    use psb_base_mod
    implicit none
    complex(psb_spk_), intent(in) :: this(:)
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this), x%v, info)
    if (info /= 0) then
      info = psb_err_alloc_request_
      call psb_errpush(info, 'c_oacc_bld_x', &
           i_err=(/size(this), izero, izero, izero, izero/))
      return
    end if

    x%v(:) = this(:)
    call x%set_host()
    if (acc_is_present(x%v))  then
      !$acc exit data delete(x%v) finalize
    end if
    !$acc enter data copyin(x%v)

  end subroutine c_oacc_bld_x

  subroutine c_oacc_asb_m(n, x, info)
    use psb_base_mod
    implicit none 
    integer(psb_mpk_), intent(in)        :: n
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_mpk_) :: nd

    info = psb_success_

    if (x%is_dev()) then
      nd = size(x%v)
      if (nd < n) then
        call x%sync()
        call x%psb_c_base_vect_type%asb(n, info)
        if (info == psb_success_) call x%sync()
        call x%set_host()
      end if
    else
      if (size(x%v) < n) then
        call x%psb_c_base_vect_type%asb(n, info)
        if (info == psb_success_) call x%sync()
        call x%set_host()
      end if
    end if
  end subroutine c_oacc_asb_m

  subroutine c_oacc_set_scal(x, val, first, last)
    class(psb_c_vect_oacc), intent(inout) :: x
    complex(psb_spk_), intent(in)           :: val
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
  end subroutine c_oacc_set_scal

  subroutine c_oacc_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x
    call x%set_dev()
    call x%set_scal(czero)
  end subroutine c_oacc_zero

  function c_oacc_get_nrows(x) result(res)
    implicit none 
    class(psb_c_vect_oacc), intent(in) :: x
    integer(psb_ipk_) :: res

    if (allocated(x%v)) res = size(x%v)
  end function c_oacc_get_nrows

  function c_oacc_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = "cOACC"

  end function c_oacc_get_fmt


  function c_oacc_vect_dot(n, x, y) result(res)
    implicit none
    class(psb_c_vect_oacc), intent(inout) :: x
    class(psb_c_base_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(in) :: n
    complex(psb_spk_) :: res
    complex(psb_spk_), external :: ddot
    integer(psb_ipk_) :: info

    res = czero
    !write(0,*) 'dot_v'
    select type(yy => y)
    type is (psb_c_base_vect_type)
        if (x%is_dev()) call x%sync()
        res = ddot(n, x%v, 1, yy%v, 1)
    type is (psb_c_vect_oacc)
        if (x%is_host()) call x%sync()
        if (yy%is_host()) call yy%sync()
        res = c_inner_oacc_dot(n, x%v, yy%v)
    class default
        call x%sync()
        res = y%dot(n, x%v)
    end select
  contains
    function c_inner_oacc_dot(n, x, y) result(res)
      implicit none 
      complex(psb_spk_), intent(in) :: x(:)
      complex(psb_spk_), intent(in) :: y(:)
      integer(psb_ipk_), intent(in) :: n
      complex(psb_spk_)  :: res
      integer(psb_ipk_) :: i
      
      !$acc parallel loop reduction(+:res) present(x, y)
      do i = 1, n
        res = res + x(i) * y(i)
      end do
      !$acc end parallel loop
    end function c_inner_oacc_dot
  end function c_oacc_vect_dot

  function c_oacc_dot_a(n, x, y) result(res)
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x
    complex(psb_spk_), intent(in) :: y(:)
    integer(psb_ipk_), intent(in) :: n
    complex(psb_spk_)  :: res
    complex(psb_spk_), external :: ddot

    if (x%is_dev()) call x%sync()
    res = ddot(n, y, 1, x%v, 1)

  end function c_oacc_dot_a

  ! subroutine c_oacc_set_vect(x,y)
  !   implicit none 
  !   class(psb_c_vect_oacc), intent(inout) :: x
  !   complex(psb_spk_), intent(in) :: y(:)
  !   integer(psb_ipk_) :: info
  
  !   if (size(x%v) /= size(y)) then 
  !     call x%free(info)
  !     call x%all(size(y),info)
  !   end if
  !   x%v(:) = y(:)
  !   call x%set_host()
  ! end subroutine c_oacc_set_vect

  subroutine c_oacc_to_dev(v)
    implicit none
    complex(psb_spk_) :: v(:)
    !$acc update device(v)          
  end subroutine c_oacc_to_dev

  subroutine c_oacc_to_host(v)
    implicit none
    complex(psb_spk_) :: v(:)
    !$acc update self(v)          
  end subroutine c_oacc_to_host

  subroutine c_oacc_sync_space(x)
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x
    if (allocated(x%v)) then
      if (.not.acc_is_present(x%v)) call c_oacc_create_dev(x%v)
    end if
  contains
    subroutine c_oacc_create_dev(v)
      implicit none
      complex(psb_spk_) :: v(:)
      !$acc enter data copyin(v)          
    end subroutine c_oacc_create_dev
  end subroutine c_oacc_sync_space

  subroutine c_oacc_sync(x)
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x
    if (x%is_dev()) then
      call c_oacc_to_host(x%v)
    end if
    if (x%is_host()) then
      call c_oacc_to_dev(x%v)
    end if
    call x%set_sync()
  end subroutine c_oacc_sync

  subroutine c_oacc_set_host(x)
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x

    x%state = is_host
  end subroutine c_oacc_set_host

  subroutine c_oacc_set_dev(x)
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x

    x%state = is_dev
  end subroutine c_oacc_set_dev

  subroutine c_oacc_set_sync(x)
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x

    x%state = is_sync
  end subroutine c_oacc_set_sync

  function c_oacc_is_dev(x) result(res)
    implicit none 
    class(psb_c_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_dev)
  end function c_oacc_is_dev

  function c_oacc_is_host(x) result(res)
    implicit none 
    class(psb_c_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_host)
  end function c_oacc_is_host

  function c_oacc_is_sync(x) result(res)
    implicit none 
    class(psb_c_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_sync)
  end function c_oacc_is_sync

  subroutine c_oacc_vect_all(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)      :: n
    class(psb_c_vect_oacc), intent(out) :: x
    integer(psb_ipk_), intent(out)     :: info

    call psb_realloc(n, x%v, info)
    if (info == 0) then
      call x%set_host()
      if (acc_is_present(x%v))  then
        !$acc exit data delete(x%v) finalize
      end if
      !$acc enter data create(x%v)
      call x%sync_space()
    end if
    if (info /= 0) then 
      info = psb_err_alloc_request_
      call psb_errpush(info, 'c_oacc_all', &
           i_err=(/n, n, n, n, n/))
    end if
  end subroutine c_oacc_vect_all

  subroutine c_oacc_vect_free(x, info)
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)     :: info
    info = 0
    if (allocated(x%v)) then
      if (acc_is_present(x%v)) then 
        !$acc exit data delete(x%v) finalize
      end if
      deallocate(x%v, stat=info)
    end if

  end subroutine c_oacc_vect_free

  function c_oacc_get_size(x) result(res)
    implicit none 
    class(psb_c_vect_oacc), intent(inout) :: x
    integer(psb_ipk_)   :: res

    if (x%is_dev()) call x%sync()
    res = size(x%v)
  end function c_oacc_get_size

end module psb_c_oacc_vect_mod
