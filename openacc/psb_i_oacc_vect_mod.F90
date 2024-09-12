module psb_i_oacc_vect_mod
  use iso_c_binding
  use openacc
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_oacc_env_mod
  use psb_i_vect_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_i_base_vect_type) :: psb_i_vect_oacc
    integer     :: state = is_host

  contains
    procedure, pass(x) :: get_nrows    => i_oacc_get_nrows
    procedure, nopass :: get_fmt       => i_oacc_get_fmt

    procedure, pass(x) :: all          => i_oacc_vect_all
    procedure, pass(x) :: zero         => i_oacc_zero
    procedure, pass(x) :: asb_m        => i_oacc_asb_m
    procedure, pass(x) :: sync         => i_oacc_sync
    procedure, pass(x) :: sync_dev_space   => i_oacc_sync_dev_space
    procedure, pass(x) :: bld_x        => i_oacc_bld_x
    procedure, pass(x) :: bld_mn       => i_oacc_bld_mn
    procedure, pass(x) :: free         => i_oacc_vect_free
    procedure, pass(x) :: free_buffer  => i_oacc_vect_free_buffer
    procedure, pass(x) :: maybe_free_buffer => i_oacc_vect_maybe_free_buffer
    procedure, pass(x) :: ins_a        => i_oacc_ins_a
    procedure, pass(x) :: ins_v        => i_oacc_ins_v
    procedure, pass(x) :: is_host      => i_oacc_is_host
    procedure, pass(x) :: is_dev       => i_oacc_is_dev
    procedure, pass(x) :: is_sync      => i_oacc_is_sync
    procedure, pass(x) :: set_host     => i_oacc_set_host
    procedure, pass(x) :: set_dev      => i_oacc_set_dev
    procedure, pass(x) :: set_sync     => i_oacc_set_sync
    procedure, pass(x) :: set_scal     => i_oacc_set_scal

    procedure, pass(x) :: new_buffer   => i_oacc_new_buffer
    procedure, pass(x) :: gthzv_x      => i_oacc_gthzv_x
    procedure, pass(x) :: gthzbuf      => i_oacc_gthzbuf
    procedure, pass(y) :: sctb         => i_oacc_sctb
    procedure, pass(y) :: sctb_x       => i_oacc_sctb_x
    procedure, pass(y) :: sctb_buf     => i_oacc_sctb_buf
    procedure, nopass  :: device_wait  => i_oacc_device_wait

    procedure, pass(x) :: get_size     => i_oacc_get_size

    final ::  i_oacc_final_vect_free
  end type psb_i_vect_oacc


contains

   subroutine i_oacc_device_wait()
     implicit none
     call acc_wait_all()
   end subroutine i_oacc_device_wait


  subroutine i_oacc_sctb_buf(i, n, idx, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    integer(psb_ipk_) :: beta
    class(psb_i_vect_oacc) :: y
    integer(psb_ipk_) :: info, k
    logical           :: acc_done
    if (.not.allocated(y%combuf)) then
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
      integer(psb_ipk_)    :: beta,x(:), y(:)
      integer(psb_ipk_) :: k
      !$acc update device(x(1:n)) async 
      !$acc parallel loop 
      do k = 1, n
        y(idx(k)) = x(k) + beta *y(idx(k))
      end do
      !$acc end parallel loop 
    end subroutine inner_sctb
    
  end subroutine i_oacc_sctb_buf

  subroutine i_oacc_sctb_x(i, n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_):: i, n
    class(psb_i_base_vect_type) :: idx
    integer(psb_ipk_) :: beta, x(:)
    class(psb_i_vect_oacc) :: y
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
      integer(psb_ipk_) :: beta, x(:), y(:)
      integer(psb_ipk_) :: k
      !$acc update device(x(1:n)) async
      !$acc parallel loop 
      do k = 1, n
        y(idx(k)) = x(k) + beta *y(idx(k))
      end do
      !$acc end parallel loop 
    end subroutine inner_sctb
    
  end subroutine i_oacc_sctb_x

  subroutine i_oacc_sctb(n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_ipk_) :: beta, x(:)
    class(psb_i_vect_oacc) :: y
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: i

    if (n == 0) return
    if (y%is_dev()) call y%sync()

    do i = 1, n
      y%v(idx(i)) = beta * y%v(idx(i)) + x(i)
    end do

    call y%set_host()
  end subroutine i_oacc_sctb

  subroutine i_oacc_gthzbuf(i, n, idx, x)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    class(psb_i_vect_oacc) :: x
    integer(psb_ipk_) :: info,k
    logical :: acc_done

    info = 0
    acc_done = .false.

    if (.not.allocated(x%combuf)) then
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
      integer(psb_ipk_)   :: x(:), y(:)
      integer(psb_ipk_) :: k
      
      !$acc parallel loop present(y)
      do k = 1, n
        y(k) = x(idx(k))
      end do
      !$acc end parallel loop
      !$acc update self(y(1:n)) async
    end subroutine inner_gth
  end subroutine i_oacc_gthzbuf
  
  subroutine i_oacc_gthzv_x(i, n, idx, x, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type):: idx
    integer(psb_ipk_) :: y(:)
    class(psb_i_vect_oacc):: x
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
      integer(psb_ipk_) :: x(:), y(:)
      integer(psb_ipk_) :: k
      
      !$acc parallel loop present(y)
      do k = 1, n
        y(k) = x(idx(k))
      end do
      !$acc end parallel loop
      !$acc update self(y(1:n)) async     
    end subroutine inner_gth
  end subroutine i_oacc_gthzv_x

  subroutine i_oacc_ins_v(n, irl, val, dupl, x, info)
    use psi_serial_mod
    implicit none
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n, dupl
    class(psb_i_base_vect_type), intent(inout) :: irl
    class(psb_i_base_vect_type), intent(inout) :: val
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, isz
    logical :: done_oacc

    info = 0
    if (psb_errstatus_fatal()) return

    done_oacc = .false.
    select type(virl  => irl)
    type is (psb_i_vect_oacc)
      select type(vval  => val)
      type is (psb_i_vect_oacc)
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
      select type(virl  => irl)
      type is (psb_i_vect_oacc)
        if (virl%is_dev()) call virl%sync()
      end select
      select type(vval  => val)
      type is (psb_i_vect_oacc)
        if (vval%is_dev()) call vval%sync()
      end select
      call x%ins(n, irl%v, val%v, dupl, info)
    end if

    if (info /= 0) then
      call psb_errpush(info, 'oacc_vect_ins')
      return
    end if

  end subroutine i_oacc_ins_v

  subroutine i_oacc_ins_a(n, irl, val, dupl, x, info)
    use psi_serial_mod
    implicit none
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n, dupl
    integer(psb_ipk_), intent(in) :: irl(:)
    integer(psb_ipk_), intent(in) :: val(:)
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i

    info = 0
    if (x%is_dev()) call x%sync()
    call x%psb_i_base_vect_type%ins(n, irl, val, dupl, info)
    call x%set_host()


  end subroutine i_oacc_ins_a

  subroutine i_oacc_bld_mn(x, n)
    use psb_base_mod
    implicit none
    integer(psb_mpk_), intent(in) :: n
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%free(info)
    call x%all(n, info)
    if (info /= 0) then
      call psb_errpush(info, 'i_oacc_bld_mn', i_err=(/n, n, n, n, n/))
    end if
    call x%set_host()
    call x%sync_dev_space()
    
  end subroutine i_oacc_bld_mn


  subroutine i_oacc_bld_x(x, this)
    use psb_base_mod
    implicit none
    integer(psb_ipk_), intent(in) :: this(:)
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%free(info)
    call psb_realloc(size(this), x%v, info)
    if (info /= 0) then
      info = psb_err_alloc_request_
      call psb_errpush(info, 'i_oacc_bld_x', &
           i_err=(/size(this), izero, izero, izero, izero/))
      return
    end if
    x%v(:) = this(:)
    call x%set_host()
    call x%sync_dev_space()

  end subroutine i_oacc_bld_x

  subroutine i_oacc_asb_m(n, x, info)
    use psb_base_mod
    implicit none 
    integer(psb_mpk_), intent(in)        :: n
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_mpk_) :: nd

    info = psb_success_

    if (x%is_dev()) then
      nd = size(x%v)
      if (nd < n) then
        call x%sync()
        call x%psb_i_base_vect_type%asb(n, info)
        if (info == psb_success_) call x%sync()
        call x%set_host()
      end if
    else
      if (size(x%v) < n) then
        call x%psb_i_base_vect_type%asb(n, info)
        if (info == psb_success_) call x%sync()
        call x%set_host()
      end if
    end if
  end subroutine i_oacc_asb_m

  subroutine i_oacc_set_scal(x, val, first, last)
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: val
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
  end subroutine i_oacc_set_scal

  subroutine i_oacc_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    call x%set_dev()
    call x%set_scal(izero)
  end subroutine i_oacc_zero

  function i_oacc_get_nrows(x) result(res)
    implicit none 
    class(psb_i_vect_oacc), intent(in) :: x
    integer(psb_ipk_) :: res

    if (allocated(x%v)) res = size(x%v)
  end function i_oacc_get_nrows

  function i_oacc_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = "iOACC"

  end function i_oacc_get_fmt


  subroutine i_oacc_new_buffer(n,x,info)
    implicit none
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info
    if (n /= psb_size(x%combuf)) then 
      call x%psb_i_base_vect_type%new_buffer(n,info)
      !$acc enter data copyin(x%combuf)
    end if
  end subroutine i_oacc_new_buffer

  subroutine i_oacc_sync_dev_space(x)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    if (allocated(x%v)) call acc_create(x%v)
  end subroutine i_oacc_sync_dev_space

  subroutine i_oacc_sync(x)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    if (x%is_dev()) then
      call acc_update_self(x%v)
    end if
    if (x%is_host()) then
      call acc_update_device(x%v)
    end if
    call x%set_sync()
  end subroutine i_oacc_sync

  subroutine i_oacc_set_host(x)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x

    x%state = is_host
  end subroutine i_oacc_set_host

  subroutine i_oacc_set_dev(x)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x

    x%state = is_dev
  end subroutine i_oacc_set_dev

  subroutine i_oacc_set_sync(x)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x

    x%state = is_sync
  end subroutine i_oacc_set_sync

  function i_oacc_is_dev(x) result(res)
    implicit none 
    class(psb_i_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_dev)
  end function i_oacc_is_dev

  function i_oacc_is_host(x) result(res)
    implicit none 
    class(psb_i_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_host)
  end function i_oacc_is_host

  function i_oacc_is_sync(x) result(res)
    implicit none 
    class(psb_i_vect_oacc), intent(in) :: x
    logical  :: res

    res = (x%state == is_sync)
  end function i_oacc_is_sync

  subroutine i_oacc_vect_all(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)      :: n
    class(psb_i_vect_oacc), intent(out) :: x
    integer(psb_ipk_), intent(out)     :: info

    call psb_realloc(n, x%v, info)
    if (info /= 0) then 
      info = psb_err_alloc_request_
      call psb_errpush(info, 'i_oacc_all', &
           i_err=(/n, n, n, n, n/))
    end if
    call x%set_host()
    call x%sync_dev_space()
  end subroutine i_oacc_vect_all

  subroutine i_oacc_final_vect_free(x)
    implicit none 
    type(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_)     :: info
    info = 0
    if (allocated(x%v)) then
      if (acc_is_present(x%v)) call acc_delete_finalize(x%v) 
      deallocate(x%v, stat=info)
    end if

  end subroutine i_oacc_final_vect_free

  subroutine i_oacc_vect_free(x, info)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)     :: info
    info = 0
    if (acc_is_present(x%v)) call acc_delete_finalize(x%v)
    if (acc_is_present(x%combuf))  call acc_delete_finalize(x%combuf)
    call x%psb_i_base_vect_type%free(info)
  end subroutine i_oacc_vect_free
  
  subroutine i_oacc_vect_maybe_free_buffer(x,info)
    implicit none
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    info = 0
    if (psb_oacc_get_maybe_free_buffer())&
         &  call x%free_buffer(info)

  end subroutine i_oacc_vect_maybe_free_buffer
  
  subroutine i_oacc_vect_free_buffer(x,info)
    implicit none
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    info = 0
    if (acc_is_present(x%combuf))  call acc_delete_finalize(x%combuf) 
    call x%psb_i_base_vect_type%free_buffer(info)

  end subroutine i_oacc_vect_free_buffer
  
  function i_oacc_get_size(x) result(res)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_)   :: res

    if (x%is_dev()) call x%sync()
    res = size(x%v)
  end function i_oacc_get_size

end module psb_i_oacc_vect_mod
