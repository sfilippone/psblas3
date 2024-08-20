module psb_i_oacc_vect_mod
  use iso_c_binding
  use openacc
  use psb_const_mod
  use psb_error_mod
  use psb_i_vect_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_i_base_vect_type) :: psb_i_vect_oacc
    integer     :: state = is_host

  contains
    procedure, pass(x) :: get_nrows   => i_oacc_get_nrows
    procedure, nopass :: get_fmt      => i_oacc_get_fmt

    procedure, pass(x) :: all         => i_oacc_vect_all
    procedure, pass(x) :: zero        => i_oacc_zero
    procedure, pass(x) :: asb_m       => i_oacc_asb_m
    procedure, pass(x) :: sync        => i_oacc_sync
    procedure, pass(x) :: sync_space  => i_oacc_sync_space
    procedure, pass(x) :: bld_x       => i_oacc_bld_x
    procedure, pass(x) :: bld_mn      => i_oacc_bld_mn
    procedure, pass(x) :: free        => i_oacc_vect_free
    procedure, pass(x) :: ins_a       => i_oacc_ins_a
    procedure, pass(x) :: ins_v       => i_oacc_ins_v
    procedure, pass(x) :: is_host     => i_oacc_is_host
    procedure, pass(x) :: is_dev      => i_oacc_is_dev
    procedure, pass(x) :: is_sync     => i_oacc_is_sync
    procedure, pass(x) :: set_host    => i_oacc_set_host
    procedure, pass(x) :: set_dev     => i_oacc_set_dev
    procedure, pass(x) :: set_sync    => i_oacc_set_sync
    procedure, pass(x) :: set_scal    => i_oacc_set_scal

    procedure, pass(x) :: gthzv_x     => i_oacc_gthzv_x
    procedure, pass(x) :: gthzbuf_x   => i_oacc_gthzbuf
    procedure, pass(y) :: sctb        => i_oacc_sctb
    procedure, pass(y) :: sctb_x      => i_oacc_sctb_x
    procedure, pass(y) :: sctb_buf    => i_oacc_sctb_buf

    procedure, pass(x) :: get_size    => i_oacc_get_size


  end type psb_i_vect_oacc


contains


  subroutine i_oacc_sctb_buf(i, n, idx, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    integer(psb_ipk_) :: beta
    class(psb_i_vect_oacc) :: y
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
  end subroutine i_oacc_sctb_buf

  subroutine i_oacc_sctb_x(i, n, idx, x, beta, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_):: i, n
    class(psb_i_base_vect_type) :: idx
    integer(psb_ipk_) :: beta, x(:)
    class(psb_i_vect_oacc) :: y
    integer(psb_ipk_) :: info, ni

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
    class default
      call psb_errpush(info, 'i_oacc_sctb_x')
      return
    end select

    if (y%is_host()) call y%sync()

    !$acc parallel loop
    do i = 1, n
      y%v(idx%v(i)) = beta * y%v(idx%v(i)) + x(i)
    end do

    call y%set_dev()
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

    !$acc parallel loop
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
      call psb_errpush(info, 'i_oacc_gthzbuf')
      return
    end select

    if (x%is_host()) call x%sync()

    !$acc parallel loop
    do i = 1, n
      x%combuf(i) = x%v(idx%v(i))
    end do
  end subroutine i_oacc_gthzbuf

  subroutine i_oacc_gthzv_x(i, n, idx, x, y)
    use psb_base_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type):: idx
    integer(psb_ipk_) :: y(:)
    class(psb_i_vect_oacc):: x
    integer(psb_ipk_) :: info

    info = 0

    select type(ii => idx)
    class is (psb_i_vect_oacc)
      if (ii%is_host()) call ii%sync()
    class default
      call psb_errpush(info, 'i_oacc_gthzv_x')
      return
    end select

    if (x%is_host()) call x%sync()

    !$acc parallel loop
    do i = 1, n
      y(i) = x%v(idx%v(i))
    end do
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
    select type(virl => irl)
    type is (psb_i_vect_oacc)
      select type(vval => val)
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
      select type(virl => irl)
      type is (psb_i_vect_oacc)
        if (virl%is_dev()) call virl%sync()
      end select
      select type(vval => val)
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
    !$acc update device(x%v)

  end subroutine i_oacc_ins_a

  subroutine i_oacc_bld_mn(x, n)
    use psb_base_mod
    implicit none
    integer(psb_mpk_), intent(in) :: n
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%all(n, info)
    if (info /= 0) then
      call psb_errpush(info, 'i_oacc_bld_mn', i_err=(/n, n, n, n, n/))
    end if
    call x%set_host()
    if (acc_is_present(x%v))  then
      !$acc exit data delete(x%v) finalize
    end if
    !$acc enter data copyin(x%v)

  end subroutine i_oacc_bld_mn


  subroutine i_oacc_bld_x(x, this)
    use psb_base_mod
    implicit none
    integer(psb_ipk_), intent(in) :: this(:)
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this), x%v, info)
    if (info /= 0) then
      info = psb_err_alloc_request_
      call psb_errpush(info, 'i_oacc_bld_x', &
           i_err=(/size(this), izero, izero, izero, izero/))
      return
    end if

    x%v(:) = this(:)
    call x%set_host()
    if (acc_is_present(x%v))  then
      !$acc exit data delete(x%v) finalize
    end if
    !$acc enter data copyin(x%v)

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

  ! subroutine i_oacc_set_vect(x,y)
  !   implicit none 
  !   class(psb_i_vect_oacc), intent(inout) :: x
  !   integer(psb_ipk_), intent(in) :: y(:)
  !   integer(psb_ipk_) :: info
  
  !   if (size(x%v) /= size(y)) then 
  !     call x%free(info)
  !     call x%all(size(y),info)
  !   end if
  !   x%v(:) = y(:)
  !   call x%set_host()
  ! end subroutine i_oacc_set_vect

  subroutine i_oacc_to_dev(v)
    implicit none
    integer(psb_ipk_) :: v(:)
    !$acc update device(v)          
  end subroutine i_oacc_to_dev

  subroutine i_oacc_to_host(v)
    implicit none
    integer(psb_ipk_) :: v(:)
    !$acc update self(v)          
  end subroutine i_oacc_to_host

  subroutine i_oacc_sync_space(x)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    if (allocated(x%v)) then
      if (.not.acc_is_present(x%v)) call i_oacc_create_dev(x%v)
    end if
  contains
    subroutine i_oacc_create_dev(v)
      implicit none
      integer(psb_ipk_) :: v(:)
      !$acc enter data copyin(v)          
    end subroutine i_oacc_create_dev
  end subroutine i_oacc_sync_space

  subroutine i_oacc_sync(x)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    if (x%is_dev()) then
      call i_oacc_to_host(x%v)
    end if
    if (x%is_host()) then
      call i_oacc_to_dev(x%v)
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
      call psb_errpush(info, 'i_oacc_all', &
           i_err=(/n, n, n, n, n/))
    end if
  end subroutine i_oacc_vect_all

  subroutine i_oacc_vect_free(x, info)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_), intent(out)     :: info
    info = 0
    if (allocated(x%v)) then
      if (acc_is_present(x%v)) then 
        !$acc exit data delete(x%v) finalize
      end if
      deallocate(x%v, stat=info)
    end if

  end subroutine i_oacc_vect_free

  function i_oacc_get_size(x) result(res)
    implicit none 
    class(psb_i_vect_oacc), intent(inout) :: x
    integer(psb_ipk_)   :: res

    if (x%is_dev()) call x%sync()
    res = size(x%v)
  end function i_oacc_get_size

end module psb_i_oacc_vect_mod
