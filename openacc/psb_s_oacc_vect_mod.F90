module psb_s_oacc_vect_mod
  use iso_c_binding
  use openacc
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_oacc_env_mod
  use psb_s_vect_mod
  use psb_i_vect_mod
  use psb_i_oacc_vect_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_s_base_vect_type) :: psb_s_vect_oacc
    integer     :: state = is_host

  contains
    procedure, pass(x) :: get_nrows    => s_oacc_get_nrows
    procedure, nopass :: get_fmt       => s_oacc_get_fmt

    procedure, pass(x) :: all          => s_oacc_vect_all
    procedure, pass(x) :: zero         => s_oacc_zero
    procedure, pass(x) :: asb_m        => s_oacc_asb_m
    procedure, pass(x) :: sync         => s_oacc_sync
    procedure, pass(x) :: sync_dev_space   => s_oacc_sync_dev_space
    procedure, pass(x) :: bld_x        => s_oacc_bld_x
    procedure, pass(x) :: bld_mn       => s_oacc_bld_mn
    procedure, pass(x) :: free         => s_oacc_vect_free
    procedure, pass(x) :: free_buffer  => s_oacc_vect_free_buffer
    procedure, pass(x) :: maybe_free_buffer => s_oacc_vect_maybe_free_buffer
    procedure, pass(x) :: ins_a        => s_oacc_ins_a
    procedure, pass(x) :: ins_v        => s_oacc_ins_v
    procedure, pass(x) :: is_host      => s_oacc_is_host
    procedure, pass(x) :: is_dev       => s_oacc_is_dev
    procedure, pass(x) :: is_sync      => s_oacc_is_sync
    procedure, pass(x) :: set_host     => s_oacc_set_host
    procedure, pass(x) :: set_dev      => s_oacc_set_dev
    procedure, pass(x) :: set_sync     => s_oacc_set_sync
    procedure, pass(x) :: set_scal     => s_oacc_set_scal

    procedure, pass(x) :: new_buffer   => s_oacc_new_buffer
    procedure, pass(x) :: gthzv_x      => s_oacc_gthzv_x
    procedure, pass(x) :: gthzbuf      => s_oacc_gthzbuf
    procedure, pass(y) :: sctb         => s_oacc_sctb
    procedure, pass(y) :: sctb_x       => s_oacc_sctb_x
    procedure, pass(y) :: sctb_buf     => s_oacc_sctb_buf
    procedure, nopass  :: device_wait  => s_oacc_device_wait

    procedure, pass(x) :: get_size     => s_oacc_get_size

    procedure, pass(x) :: dot_v        => s_oacc_vect_dot
    procedure, pass(x) :: dot_a        => s_oacc_dot_a
    procedure, pass(y) :: axpby_v      => s_oacc_axpby_v
    procedure, pass(y) :: axpby_a      => s_oacc_axpby_a
    procedure, pass(z) :: upd_xyz      => s_oacc_upd_xyz
    procedure, pass(y) :: mlt_a        => s_oacc_mlt_a
    procedure, pass(z) :: mlt_a_2      => s_oacc_mlt_a_2
    procedure, pass(y) :: mlt_v        => psb_s_oacc_mlt_v
    procedure, pass(z) :: mlt_v_2      => psb_s_oacc_mlt_v_2
    procedure, pass(x) :: scal         => s_oacc_scal 
    procedure, pass(x) :: nrm2         => s_oacc_nrm2
    procedure, pass(x) :: amax         => s_oacc_amax
    procedure, pass(x) :: asum         => s_oacc_asum
    procedure, pass(x) :: absval1      => s_oacc_absval1
    procedure, pass(x) :: absval2      => s_oacc_absval2
    final ::  s_oacc_final_vect_free
  end type psb_s_vect_oacc

  interface
    subroutine psb_s_oacc_mlt_v(x, y, info)
      import
      implicit none 
      class(psb_s_base_vect_type), intent(inout) :: x
      class(psb_s_vect_oacc), intent(inout)       :: y
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psb_s_oacc_mlt_v
  end interface
  
  interface
    subroutine psb_s_oacc_mlt_v_2(alpha, x, y, beta, z, info, conjgx, conjgy)
      import
      implicit none 
      real(psb_spk_), intent(in)                 :: alpha, beta
      class(psb_s_base_vect_type), intent(inout) :: x
      class(psb_s_base_vect_type), intent(inout) :: y
      class(psb_s_vect_oacc), intent(inout)      :: z
      integer(psb_ipk_), intent(out)             :: info
      character(len=1), intent(in), optional     :: conjgx, conjgy
    end subroutine psb_s_oacc_mlt_v_2
  end interface

contains

   subroutine s_oacc_device_wait()
     implicit none
     call acc_wait_all()
   end subroutine s_oacc_device_wait

  subroutine s_oacc_absval1(x)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: n

    if (x%is_host()) call x%sync()
    n = size(x%v)
    call s_inner_oacc_absval1(n,x%v)
    call x%set_dev()
  contains 
    subroutine s_inner_oacc_absval1(n,x)
      implicit none
      real(psb_spk_), intent(inout) :: x(:)
      integer(psb_ipk_) :: n
      integer(psb_ipk_) :: i
      !$acc parallel loop present(x)
      do i = 1, n
        x(i) = abs(x(i))
      end do
    end subroutine s_inner_oacc_absval1
  end subroutine s_oacc_absval1

  subroutine s_oacc_absval2(x, y)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    class(psb_s_base_vect_type), intent(inout) :: y
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: i

    n = min(size(x%v), size(y%v))
    select type (yy  => y)
    class is (psb_s_vect_oacc)
        if (x%is_host()) call x%sync()
        if (yy%is_host()) call yy%sync()
        call s_inner_oacc_absval2(n,x%v,yy%v)
    class default
        if (x%is_dev()) call x%sync()
        if (y%is_dev()) call y%sync()
        call x%psb_s_base_vect_type%absval(y)
    end select
  contains 
    subroutine s_inner_oacc_absval2(n,x,y)
      implicit none
      real(psb_spk_), intent(inout) :: x(:),y(:)
      integer(psb_ipk_) :: n
      integer(psb_ipk_) :: i
      !$acc parallel loop present(x,y)
      do i = 1, n
        y(i) = abs(x(i))
      end do
    end subroutine s_inner_oacc_absval2
  end subroutine s_oacc_absval2

  subroutine s_oacc_scal(alpha, x)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    real(psb_spk_), intent(in)          :: alpha
    integer(psb_ipk_) :: info
    if (x%is_host()) call x%sync()
    call s_inner_oacc_scal(alpha, x%v)
    call x%set_dev()
  contains
    subroutine s_inner_oacc_scal(alpha, x)
      real(psb_spk_), intent(in) :: alpha
      real(psb_spk_), intent(inout) :: x(:)
      integer(psb_ipk_) :: i
      !$acc parallel loop present(x)
      do i = 1, size(x)
        x(i) = alpha * x(i)
      end do
    end subroutine s_inner_oacc_scal
  end subroutine s_oacc_scal

  function s_oacc_nrm2(n, x) result(res)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                       :: res
    integer(psb_ipk_) :: info

    if (x%is_host()) call x%sync()
!!$    write(0,*)'oacc_nrm2'
    res = s_inner_oacc_nrm2(n, x%v)
  contains
    function s_inner_oacc_nrm2(n, x) result(res)
      integer(psb_ipk_) :: n
      real(psb_spk_) :: x(:)
      real(psb_spk_) :: res
      real(psb_spk_) :: sum, mx
      integer(psb_ipk_) :: i
      mx = szero
      !$acc parallel loop reduction(max:mx) present(x)
      do i = 1, n
        if (abs(x(i)) > mx) mx = abs(x(i))
      end do
      if (mx == szero) then
        res = mx
      else
        sum = szero
        !$acc parallel loop reduction(+:sum)  present(x)
        do i = 1, n
          sum = sum + abs(x(i)/mx)**2
        end do
        res = mx*sqrt(sum)
      end if
    end function s_inner_oacc_nrm2
  end function s_oacc_nrm2

  function s_oacc_amax(n, x) result(res)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                       :: res
    integer(psb_ipk_) :: info

    if (x%is_host()) call x%sync()
    res = s_inner_oacc_amax(n, x%v)
  contains
    function s_inner_oacc_amax(n, x) result(res)
      integer(psb_ipk_) :: n
      real(psb_spk_) :: x(:)
      real(psb_spk_) :: res
      real(psb_spk_) :: max_val
      integer(psb_ipk_) :: i
      max_val = szero
      !$acc parallel loop reduction(max:max_val) present(x) 
      do i = 1, n
        if (abs(x(i)) > max_val) max_val = abs(x(i))
      end do
      res = max_val
    end function s_inner_oacc_amax
  end function s_oacc_amax
  
  function s_oacc_asum(n, x) result(res)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                       :: res
    integer(psb_ipk_) :: info
    real(psb_spk_) :: sum
    integer(psb_ipk_) :: i
    if (x%is_host()) call x%sync()
    res =  s_inner_oacc_asum(n, x%v)
  contains
    function s_inner_oacc_asum(n, x) result(res)
      integer(psb_ipk_) :: n
      real(psb_spk_) :: x(:)
      real(psb_spk_) :: res
      integer(psb_ipk_) :: i
      res = szero
      !$acc parallel loop reduction(+:res) present(x)
      do i = 1, n
        res = res + abs(x(i))
      end do
    end function s_inner_oacc_asum
  end function s_oacc_asum


  subroutine s_oacc_mlt_a(x, y, info)
    implicit none 
    real(psb_spk_), intent(in)           :: x(:)
    class(psb_s_vect_oacc), intent(inout) :: y
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: i, n
    
    info = 0    
    if (y%is_dev()) call y%sync()
    !$acc parallel loop present(x,y)
    do i = 1, size(x)
        y%v(i) = y%v(i) * x(i)
    end do
    call y%set_host()
  end subroutine s_oacc_mlt_a

  subroutine s_oacc_mlt_a_2(alpha, x, y, beta, z, info)
    implicit none 
    real(psb_spk_), intent(in)           :: alpha, beta
    real(psb_spk_), intent(in)           :: x(:)
    real(psb_spk_), intent(in)           :: y(:)
    class(psb_s_vect_oacc), intent(inout) :: z
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: i, n
    
    info = 0    
    if (z%is_dev()) call z%sync()
    !$acc parallel loop  present(x,y,z%v)
    do i = 1, size(x)
        z%v(i) = alpha * x(i) * y(i) + beta * z%v(i)
    end do
    call z%set_host()
  end subroutine s_oacc_mlt_a_2

  subroutine s_oacc_axpby_v(m, alpha, x, beta, y, info)
    !use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    class(psb_s_base_vect_type), intent(inout) :: x
    class(psb_s_vect_oacc), intent(inout) :: y
    real(psb_spk_), intent(in) :: alpha, beta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: nx, ny, i

    info = psb_success_

    select type(xx  => x)
    type is (psb_s_vect_oacc)
        if ((beta /= szero) .and. y%is_host()) call y%sync()
        if (xx%is_host()) call xx%sync()
        nx = size(xx%v)
        ny = size(y%v)
        if ((nx < m) .or. (ny < m)) then
            info = psb_err_internal_error_
        else
          call s_inner_oacc_axpby(m, alpha, x%v, beta, y%v, info)
        end if
        call y%set_dev()
    class default
        if ((alpha /= szero) .and. (x%is_dev())) call x%sync()
        call y%axpby(m, alpha, x%v, beta, info)
      end select
    contains
      subroutine s_inner_oacc_axpby(m, alpha, x, beta, y, info)
      !use psi_serial_mod
      implicit none
      integer(psb_ipk_), intent(in) :: m
      real(psb_spk_), intent(inout) :: x(:)
      real(psb_spk_), intent(inout) :: y(:)
      real(psb_spk_), intent(in) :: alpha, beta
      integer(psb_ipk_), intent(out)  :: info
      !$acc parallel present(x,y)
      !$acc loop
      do i = 1, m
        y(i) = alpha * x(i) + beta * y(i)
      end do
      !$acc end parallel 
    end subroutine s_inner_oacc_axpby
  end subroutine s_oacc_axpby_v

  subroutine s_oacc_axpby_a(m, alpha, x, beta, y, info)
    !use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    real(psb_spk_), intent(in) :: x(:)
    class(psb_s_vect_oacc), intent(inout) :: y
    real(psb_spk_), intent(in) :: alpha, beta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: i

    if ((beta /= szero) .and. (y%is_dev())) call y%sync()

    do i = 1, m
        y%v(i) = alpha * x(i) + beta * y%v(i)
    end do
    call y%set_host()
  end subroutine s_oacc_axpby_a

  subroutine s_oacc_upd_xyz(m, alpha, beta, gamma, delta, x, y, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in) :: m
    class(psb_s_base_vect_type), intent(inout) :: x
    class(psb_s_base_vect_type), intent(inout) :: y
    class(psb_s_vect_oacc), intent(inout) :: z
    real(psb_spk_), intent(in) :: alpha, beta, gamma, delta
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: nx, ny, nz, i
    logical :: gpu_done

    info = psb_success_
    gpu_done = .false.

    select type(xx  => x)
    class is (psb_s_vect_oacc)
      select type(yy  => y)
      class is (psb_s_vect_oacc)
        select type(zz  => z)
        class is (psb_s_vect_oacc)
          if ((beta /= szero) .and. yy%is_host()) call yy%sync()
          if ((delta /= szero) .and. zz%is_host()) call zz%sync()
          if (xx%is_host()) call xx%sync()
          nx = size(xx%v)
          ny = size(yy%v)
          nz = size(zz%v)
          if ((nx < m) .or. (ny < m) .or. (nz < m)) then
            info = psb_err_internal_error_
          else
            !$acc parallel loop  present(xx%v,yy%v,zz%v)
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
  end subroutine s_oacc_upd_xyz

  subroutine s_oacc_sctb_buf(i, n, idx, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    real(psb_spk_) :: beta
    class(psb_s_vect_oacc) :: y
    integer(psb_ipk_) :: info, k
    logical           :: acc_done
    if (.not.allocated(y%combuf)) then
      write(0,*) 'allocation error for y%combuf '
      call psb_errpush(psb_err_alloc_dealloc_, 'sctb_buf')
      return
    end if

    acc_done = .false. 
    select type(ii  => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
      if (y%is_host()) call y%sync()
      call inner_sctb(n,y%combuf(i:i+n-1),beta,y%v,ii%v(i:i+n-1))
      call y%set_dev()
      acc_done = .true.
    end select

    if (.not.acc_done) then 
      if (idx%is_dev()) call idx%sync()
      if (y%is_dev()) call y%sync()
      do k = 1, n
        y%v(idx%v(k+i-1)) = beta * y%v(idx%v(k+i-1)) + y%combuf(k)
      end do
    end if

  contains
    subroutine inner_sctb(n,x,beta,y,idx)
      integer(psb_ipk_) :: n, idx(:)
      real(psb_spk_)    :: beta,x(:), y(:)
      integer(psb_ipk_) :: k
      !$acc update device(x(1:n)) 
      !$acc parallel loop present(x,y)
      do k = 1, n
        y(idx(k)) = x(k) + beta *y(idx(k))
      end do
      !$acc end parallel loop 
    end subroutine inner_sctb
    
  end subroutine s_oacc_sctb_buf

  subroutine s_oacc_sctb_x(i, n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_):: i, n
    class(psb_i_base_vect_type) :: idx
    real(psb_spk_) :: beta, x(:)
    class(psb_s_vect_oacc) :: y
    integer(psb_ipk_) :: info, ni, k 
    logical :: acc_done 

    acc_done = .false. 
    select type(ii  => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
      if (y%is_host()) call y%sync()
      if (acc_is_present(x)) then
        call inner_sctb(n,x(i:i+n-1),beta,y%v,idx%v(i:i+n-1))
        acc_done = .true.
        call y%set_dev()
      end if
    end select
    if (.not.acc_done) then
      if (idx%is_dev()) call idx%sync()
      if (y%is_dev()) call y%sync()
      do k = 1, n
        y%v(idx%v(k+i-1)) = beta * y%v(idx%v(k+i-1)) + x(k+i-1)
      end do
      call y%set_host()
    end if

  contains
    subroutine inner_sctb(n,x,beta,y,idx)
      integer(psb_ipk_) :: n, idx(:)
      real(psb_spk_) :: beta, x(:), y(:)
      integer(psb_ipk_) :: k
      !$acc update device(x(1:n)) 
      !$acc parallel loop  present(x,y)
      do k = 1, n
        y(idx(k)) = x(k) + beta *y(idx(k))
      end do
      !$acc end parallel loop 
    end subroutine inner_sctb
    
  end subroutine s_oacc_sctb_x

  subroutine s_oacc_sctb(n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: idx(:)
    real(psb_spk_) :: beta, x(:)
    class(psb_s_vect_oacc) :: y
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: i

    if (n == 0) return
    if (y%is_dev()) call y%sync()

    do i = 1, n
      y%v(idx(i)) = beta * y%v(idx(i)) + x(i)
    end do

    call y%set_host()
  end subroutine s_oacc_sctb

  subroutine s_oacc_gthzbuf(i, n, idx, x)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    class(psb_s_vect_oacc) :: x
    integer(psb_ipk_) :: info,k
    logical :: acc_done

    info = 0
    acc_done = .false.

    if (.not.allocated(x%combuf)) then
      write(0,*) 'oacc allocation error combuf gthzbuf '    
      call psb_errpush(psb_err_alloc_dealloc_, 'gthzbuf')
      return
    end if

    select type (ii  => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
      if (x%is_host()) call x%sync()
      call inner_gth(n,x%v,x%combuf(i:i+n-1),ii%v(i:i+n-1))
      acc_done = .true.      
    end select

    if (.not.acc_done) then
      if (idx%is_dev()) call idx%sync()
      if (x%is_dev()) call x%sync()
      do k = 1, n
        x%combuf(k+i-1) = x%v(idx%v(k+i-1))
      end do
    end if

  contains
    subroutine inner_gth(n,x,y,idx)
      integer(psb_ipk_) :: n, idx(:)
      real(psb_spk_)   :: x(:), y(:)
      integer(psb_ipk_) :: k
      !
      !$acc parallel loop present(x,y)
      do k = 1, n
        y(k) = x(idx(k))
      end do
      !$acc end parallel loop
      !$acc update self(y(1:n)) 
    end subroutine inner_gth
  end subroutine s_oacc_gthzbuf
  
  subroutine s_oacc_gthzv_x(i, n, idx, x, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type):: idx
    real(psb_spk_) :: y(:)
    class(psb_s_vect_oacc):: x
    integer(psb_ipk_) :: info, k 
    logical :: acc_done 

    info = 0
    acc_done = .false. 
    select type (ii  => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
      if (x%is_host()) call x%sync()
      if (acc_is_present(y)) then 
        call inner_gth(n,x%v,y(i:),ii%v(i:))
        acc_done=.true.
      end if
    end select
    if (.not.acc_done) then 
      if (x%is_dev()) call x%sync()
      if (idx%is_dev()) call idx%sync()
      do k = 1, n
        y(k+i-1) = x%v(idx%v(k+i-1))
        !write(0,*) 'oa gthzv ',k+i-1,idx%v(k+i-1),k,y(k)
      end do
    end if
  contains
    subroutine inner_gth(n,x,y,idx)
      integer(psb_ipk_) :: n, idx(:)
      real(psb_spk_) :: x(:), y(:)
      integer(psb_ipk_) :: k
      !
      !$acc parallel loop present(x,y)
      do k = 1, n
        y(k) = x(idx(k))
      end do
      !$acc end parallel loop
      !$acc update self(y(1:n))      
    end subroutine inner_gth
  end subroutine s_oacc_gthzv_x

  subroutine s_oacc_ins_v(n, irl, val, dupl, x, info)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n, dupl
    class(psb_i_base_vect_type), intent(inout) :: irl
    class(psb_s_base_vect_type), intent(inout) :: val
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, isz
    logical :: done_oacc

    info = 0
    if (psb_errstatus_fatal()) return

    done_oacc = .false.
    select type(virl  => irl)
    type is (psb_i_vect_oacc)
      select type(vval  => val)
      type is (psb_s_vect_oacc)
        if (vval%is_host()) call vval%sync()
        if (virl%is_host()) call virl%sync()
        if (x%is_host()) call x%sync()
        !$acc parallel loop  present(x%v,virl%v,vval%v)
        do i = 1, n
          x%v(virl%v(i)) = vval%v(i)
        end do
        call x%set_dev()
        done_oacc = .true.
      end select
    end select

    if (.not.done_oacc) then
      select type(virl  => irl)
      type is (psb_i_vect_oacc)
        if (virl%is_dev()) call virl%sync()
      end select
      select type(vval  => val)
      type is (psb_s_vect_oacc)
        if (vval%is_dev()) call vval%sync()
      end select
      call x%ins(n, irl%v, val%v, dupl, info)
    end if

    if (info /= 0) then
      call psb_errpush(info, 'oacc_vect_ins')
      return
    end if

  end subroutine s_oacc_ins_v

  subroutine s_oacc_ins_a(n, irl, val, dupl, x, info)
    use psi_serial_mod
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n, dupl
    integer(psb_ipk_), intent(in) :: irl(:)
    real(psb_spk_), intent(in) :: val(:)
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i

    info = 0
    if (x%is_dev()) call x%sync()
    call x%psb_s_base_vect_type%ins(n, irl, val, dupl, info)
    call x%set_host()


  end subroutine s_oacc_ins_a

  subroutine s_oacc_bld_mn(x, n)
    use psb_base_mod
    implicit none
    integer(psb_mpk_), intent(in) :: n
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%free(info)
    call x%all(n, info)
    if (info /= 0) then
      call psb_errpush(info, 's_oacc_bld_mn', i_err=(/n, n, n, n, n/))
    end if
    call x%set_host()
    call x%sync_dev_space()
    
  end subroutine s_oacc_bld_mn


  subroutine s_oacc_bld_x(x, this)
    use psb_base_mod
    implicit none
    real(psb_spk_), intent(in) :: this(:)
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%free(info)
    call psb_realloc(size(this), x%v, info)
    if (info /= 0) then
      info = psb_err_alloc_request_
      call psb_errpush(info, 's_oacc_bld_x', &
           i_err=(/size(this), izero, izero, izero, izero/))
      return
    end if
    x%v(:) = this(:)
    call x%set_host()
    call x%sync_dev_space()

  end subroutine s_oacc_bld_x

  subroutine s_oacc_asb_m(n, x, info)
    use psb_base_mod
    implicit none 
    integer(psb_mpk_), intent(in)        :: n
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_mpk_) :: nd

    info = psb_success_

    if (x%is_dev()) then
      nd = size(x%v)
      if (nd < n) then
        call x%sync()
        call x%psb_s_base_vect_type%asb(n, info)
        if (info == psb_success_) call x%sync()
        call x%set_host()
      end if
    else
      if (size(x%v) < n) then
        call x%psb_s_base_vect_type%asb(n, info)
        if (info == psb_success_) call x%sync()
        call x%set_host()
      end if
    end if
  end subroutine s_oacc_asb_m

  subroutine s_oacc_set_scal(x, val, first, last)
    class(psb_s_vect_oacc), intent(inout) :: x
    real(psb_spk_), intent(in)           :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: first_, last_
    first_ = 1
    last_  = x%get_nrows()
    if (present(first)) first_ = max(1, first)
    if (present(last))  last_  = min(last, last_)

    !$acc parallel loop present(x%v)
    do i = first_, last_
      x%v(i) = val
    end do
    !$acc end parallel loop

    call x%set_dev()
  end subroutine s_oacc_set_scal

  subroutine s_oacc_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x
    call x%set_dev()
    call x%set_scal(szero)
  end subroutine s_oacc_zero

  function s_oacc_get_nrows(x) result(res)
    implicit none 
    class(psb_s_vect_oacc), intent(in) :: x
    integer(psb_ipk_) :: res

    if (allocated(x%v)) res = size(x%v)
  end function s_oacc_get_nrows

  function s_oacc_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = "sOACC"

  end function s_oacc_get_fmt


  function s_oacc_vect_dot(n, x, y) result(res)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    class(psb_s_base_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(in) :: n
    real(psb_spk_) :: res
    real(psb_spk_), external :: ddot
    integer(psb_ipk_) :: info

    res = szero
!!$    write(0,*) 'oacc_dot_v'
    select type(yy  => y)
    type is (psb_s_base_vect_type)
        if (x%is_dev()) call x%sync()
        res = ddot(n, x%v, 1, yy%v, 1)
    type is (psb_s_vect_oacc)
        if (x%is_host()) call x%sync()
        if (yy%is_host()) call yy%sync()
        res = s_inner_oacc_dot(n, x%v, yy%v)
    class default
        call x%sync()
        res = y%dot(n, x%v)
    end select
  contains
    function s_inner_oacc_dot(n, x, y) result(res)
      implicit none 
      real(psb_spk_), intent(in) :: x(:)
      real(psb_spk_), intent(in) :: y(:)
      integer(psb_ipk_), intent(in) :: n
      real(psb_spk_)  :: res
      integer(psb_ipk_) :: i
      
      !$acc parallel loop reduction(+:res) present(x, y)
      do i = 1, n
        res = res + x(i) * y(i)
      end do
      !$acc end parallel loop
    end function s_inner_oacc_dot
  end function s_oacc_vect_dot

  function s_oacc_dot_a(n, x, y) result(res)
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x
    real(psb_spk_), intent(in) :: y(:)
    integer(psb_ipk_), intent(in) :: n
    real(psb_spk_)  :: res
    real(psb_spk_), external :: ddot

    if (x%is_dev()) call x%sync()
    res = ddot(n, y, 1, x%v, 1)

  end function s_oacc_dot_a


  subroutine s_oacc_new_buffer(n,x,info)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    !write(0,*) 'oacc new_buffer',n,psb_size(x%combuf)    
    if (n > psb_size(x%combuf)) then
      !write(0,*) 'oacc new_buffer: reallocating '
      if (allocated(x%combuf)) then
        !if (acc_is_present(x%combuf)) call acc_delete_finalize(x%combuf)
        !$acc exit data delete(x%combuf) 
      end if
      call x%psb_s_base_vect_type%new_buffer(n,info)
      !$acc enter data copyin(x%combuf)
      ! call acc_copyin(x%combuf)
    end if
  end subroutine s_oacc_new_buffer

  subroutine s_oacc_sync_dev_space(x)
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x
!!$    write(0,*) 'oacc sync_dev_space'    
    if (psb_size(x%v)>0) call acc_copyin(x%v)
  end subroutine s_oacc_sync_dev_space

  subroutine s_oacc_sync(x)
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x
    if (x%is_dev()) then
      if (psb_size(x%v)>0) call acc_update_self(x%v)
    end if
    if (x%is_host()) then
      if (.not.acc_is_present(x%v)) call s_oacc_sync_dev_space(x)      
      if (psb_size(x%v)>0) call acc_update_device(x%v)
    end if
    call x%set_sync()
  end subroutine s_oacc_sync

  subroutine s_oacc_set_host(x)
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x

    x%state = is_host
  end subroutine s_oacc_set_host

  subroutine s_oacc_set_dev(x)
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x

    x%state = is_dev
  end subroutine s_oacc_set_dev

  subroutine s_oacc_set_sync(x)
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x

    x%state = is_sync
  end subroutine s_oacc_set_sync

  function s_oacc_is_dev(x) result(res)
    implicit none 
    class(psb_s_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_dev)
  end function s_oacc_is_dev

  function s_oacc_is_host(x) result(res)
    implicit none 
    class(psb_s_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_host)
  end function s_oacc_is_host

  function s_oacc_is_sync(x) result(res)
    implicit none 
    class(psb_s_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_sync)
  end function s_oacc_is_sync

  subroutine s_oacc_vect_all(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)      :: n
    class(psb_s_vect_oacc), intent(out) :: x
    integer(psb_ipk_), intent(out)     :: info

    call psb_realloc(n, x%v, info)
    if (info /= 0) then 
      info = psb_err_alloc_request_
      call psb_errpush(info, 's_oacc_all', &
           i_err=(/n, n, n, n, n/))
    end if
    call x%set_host()
    call x%sync_dev_space()
  end subroutine s_oacc_vect_all

  subroutine s_oacc_final_vect_free(x)
    implicit none 
    type(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_)     :: info
    info = 0
!!$    write(0,*) 'oacc final_vect_free'
    call x%free_buffer(info)
    if (allocated(x%v)) then
      if (acc_is_present(x%v)) call acc_delete_finalize(x%v) 
      deallocate(x%v, stat=info)
    end if

  end subroutine s_oacc_final_vect_free

  subroutine s_oacc_vect_free(x, info)
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)     :: info
    info = 0
!!$    write(0,*) 'oacc vect_free'
    call x%free_buffer(info)
    if (acc_is_present(x%v)) call acc_delete_finalize(x%v)
    call x%psb_s_base_vect_type%free(info)
  end subroutine s_oacc_vect_free
  
  subroutine s_oacc_vect_maybe_free_buffer(x,info)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    info = 0
    if (psb_oacc_get_maybe_free_buffer()) then
      !write(0,*) 'psb_oacc_get_maybe_free_buffer() ',psb_oacc_get_maybe_free_buffer()
      call x%free_buffer(info)
    end if

  end subroutine s_oacc_vect_maybe_free_buffer
  
  subroutine s_oacc_vect_free_buffer(x,info)
    implicit none
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info
!    write(0,*) 'oacc free_buffer'    
    info = 0
    if (acc_is_present(x%combuf))  call acc_delete_finalize(x%combuf) 
    call x%psb_s_base_vect_type%free_buffer(info)

  end subroutine s_oacc_vect_free_buffer
  
  function s_oacc_get_size(x) result(res)
    implicit none 
    class(psb_s_vect_oacc), intent(inout) :: x
    integer(psb_ipk_)   :: res

    res = size(x%v)
  end function s_oacc_get_size

end module psb_s_oacc_vect_mod
