module psb_i_oacc_vect_mod
    use iso_c_binding
    use psb_const_mod
    use psb_error_mod
    use psb_i_vect_mod

    integer(psb_ipk_), parameter, private :: is_host = -1
    integer(psb_ipk_), parameter, private :: is_sync = 0 
    integer(psb_ipk_), parameter, private :: is_dev  = 1

    type, extends(psb_i_base_vect_type) :: psb_i_vect_oacc
        integer :: state = is_host
    contains
        procedure, pass(x) :: get_nrows => i_oacc_get_nrows
        procedure, nopass  :: get_fmt   => i_oacc_get_fmt

        procedure, pass(x) :: all      => i_oacc_all
        procedure, pass(x) :: zero     => i_oacc_zero
        procedure, pass(x) :: asb_m    => i_oacc_asb_m
        procedure, pass(x) :: sync     => i_oacc_sync
        procedure, pass(x) :: sync_space => i_oacc_sync_space
        procedure, pass(x) :: bld_x    => i_oacc_bld_x
        procedure, pass(x) :: bld_mn   => i_oacc_bld_mn
        procedure, pass(x) :: free     => i_oacc_free
        procedure, pass(x) :: ins_a    => i_oacc_ins_a
        procedure, pass(x) :: ins_v    => i_oacc_ins_v
        procedure, pass(x) :: is_host  => i_oacc_is_host
        procedure, pass(x) :: is_dev   => i_oacc_is_dev
        procedure, pass(x) :: is_sync  => i_oacc_is_sync
        procedure, pass(x) :: set_host => i_oacc_set_host
        procedure, pass(x) :: set_dev  => i_oacc_set_dev
        procedure, pass(x) :: set_sync => i_oacc_set_sync
        procedure, pass(x) :: set_scal => i_oacc_set_scal
        procedure, pass(x) :: gthzv_x  => i_oacc_gthzv_x
        procedure, pass(y) :: sctb     => i_oacc_sctb
        procedure, pass(y) :: sctb_x   => i_oacc_sctb_x
        procedure, pass(x) :: gthzbuf  => i_oacc_gthzbuf
        procedure, pass(y) :: sctb_buf => i_oacc_sctb_buf

        final              :: i_oacc_vect_finalize
    end type psb_i_vect_oacc

    public  :: psb_i_vect_oacc_
    private :: constructor
    interface psb_i_vect_oacc_
        module procedure constructor
    end interface psb_i_vect_oacc_

contains

    function constructor(x) result(this)
        integer(psb_ipk_)       :: x(:)
        type(psb_i_vect_oacc) :: this
        integer(psb_ipk_) :: info

        this%v = x
        call this%asb(size(x), info)
    end function constructor


    subroutine i_oacc_gthzv_x(i, n, idx, x, y)
        implicit none
        integer(psb_ipk_) :: i, n
        class(psb_i_base_vect_type) :: idx
        integer(psb_ipk_) :: y(:)
        class(psb_i_vect_oacc)  :: x
        integer :: info

        !$acc parallel loop
        do i = 1, n
            y(i) = x%v(idx%v(i))
        end do
    end subroutine i_oacc_gthzv_x

    subroutine i_oacc_gthzbuf(i, n, idx, x)
        implicit none
        integer(psb_ipk_) :: i, n
        class(psb_i_base_vect_type) :: idx
        class(psb_i_vect_oacc) :: x
        integer :: info

        if (.not.allocated(x%combuf)) then
            call psb_errpush(psb_err_alloc_dealloc_, 'gthzbuf')
            return
        end if

        !$acc parallel loop
        do i = 1, n
            x%combuf(i) = x%v(idx%v(i))
        end do
    end subroutine i_oacc_gthzbuf

    subroutine i_oacc_sctb(n, idx, x, beta, y)
        implicit none
        integer(psb_ipk_) :: n, idx(:)
        integer(psb_ipk_) :: beta, x(:)
        class(psb_i_vect_oacc) :: y
        integer(psb_ipk_) :: info
        integer :: i

        if (n == 0) return

        !$acc parallel loop
        do i = 1, n
            y%v(idx(i)) = beta * y%v(idx(i)) + x(i)
        end do
    end subroutine i_oacc_sctb

    subroutine i_oacc_sctb_x(i, n, idx, x, beta, y)
        implicit none
        integer(psb_ipk_) :: i, n
        class(psb_i_base_vect_type) :: idx
        integer(psb_ipk_) :: beta, x(:)
        class(psb_i_vect_oacc) :: y
        integer :: info
    
        select type(ii => idx)
        class is (psb_i_vect_oacc)
            if (ii%is_host()) call ii%sync_space(info)
            if (y%is_host()) call y%sync_space(info)
    
            !$acc parallel loop
            do i = 1, n
                y%v(ii%v(i)) = beta * y%v(ii%v(i)) + x(i)
            end do
    
        class default
            !$acc parallel loop
            do i = 1, n
                y%v(idx%v(i)) = beta * y%v(idx%v(i)) + x(i)
            end do
        end select
    end subroutine i_oacc_sctb_x
    
    subroutine i_oacc_sctb_buf(i, n, idx, beta, y)
        implicit none
        integer(psb_ipk_) :: i, n
        class(psb_i_base_vect_type) :: idx
        integer(psb_ipk_) :: beta
        class(psb_i_vect_oacc) :: y
        integer(psb_ipk_) :: info
    
        if (.not.allocated(y%v)) then
            call psb_errpush(psb_err_alloc_dealloc_, 'sctb_buf')
            return
        end if
    
        !$acc parallel loop
        do i = 1, n
            y%v(idx%v(i)) = beta * y%v(idx%v(i)) + y%v(i)
        end do
    end subroutine i_oacc_sctb_buf

    subroutine i_oacc_set_host(x)
        class(psb_i_vect_oacc), intent(inout) :: x
        x%state = is_host
    end subroutine i_oacc_set_host

    subroutine i_oacc_set_sync(x)
        class(psb_i_vect_oacc), intent(inout) :: x
        x%state = is_sync
    end subroutine i_oacc_set_sync

    subroutine i_oacc_set_dev(x)
        class(psb_i_vect_oacc), intent(inout) :: x
        x%state = is_dev
    end subroutine i_oacc_set_dev
    
    subroutine i_oacc_set_scal(x, val, first, last)
        class(psb_i_vect_oacc), intent(inout) :: x
        integer(psb_ipk_), intent(in) :: val
        integer(psb_ipk_), optional :: first, last
        
        integer(psb_ipk_) :: first_, last_
    
        first_ = 1
        last_ = size(x%v)
        if (present(first)) first_ = max(1, first)
        if (present(last)) last_ = min(size(x%v), last)
        
        !$acc parallel loop
        do i = first_, last_
            x%v(i) = val
        end do
        call x%set_dev()
    end subroutine i_oacc_set_scal

    function i_oacc_is_host(x) result(res)
        implicit none
        class(psb_i_vect_oacc), intent(in) :: x
        logical :: res
    
        res = (x%state == is_host)
    end function i_oacc_is_host
    
    function i_oacc_is_dev(x) result(res)
        implicit none
        class(psb_i_vect_oacc), intent(in) :: x
        logical :: res
    
        res = (x%state == is_dev)
    end function i_oacc_is_dev
    
    function i_oacc_is_sync(x) result(res)
        implicit none
        class(psb_i_vect_oacc), intent(in) :: x
        logical  :: res
    
        res = (x%state == is_sync)
    end function i_oacc_is_sync

    subroutine i_oacc_free(x, info)
        use psb_error_mod
        implicit none 
        class(psb_i_vect_oacc), intent(inout)  :: x
        integer(psb_ipk_), intent(out)        :: info
        
        info = 0  
        if (allocated(x%v)) deallocate(x%v, stat=info)
        if (info /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info, 'i_oacc_free')
        end if
        call x%set_sync()
    end subroutine i_oacc_free
    
    subroutine i_oacc_ins_a(n, irl, val, dupl, x, info)
        implicit none 
        class(psb_i_vect_oacc), intent(inout) :: x
        integer(psb_ipk_), intent(in)        :: n, dupl
        integer(psb_ipk_), intent(in)        :: irl(:)
        integer(psb_ipk_), intent(in)        :: val(:)
        integer(psb_ipk_), intent(out)       :: info
    
        integer(psb_ipk_) :: i
    
        info = 0
        if (x%is_dev()) call x%sync()
        call x%psb_i_base_vect_type%ins(n, irl, val, dupl, info)
        call x%set_host()
    end subroutine i_oacc_ins_a
    
    subroutine i_oacc_ins_v(n, irl, val, dupl, x, info)
        implicit none 
        class(psb_i_vect_oacc), intent(inout)        :: x
        integer(psb_ipk_), intent(in)               :: n, dupl
        class(psb_i_base_vect_type), intent(inout)  :: irl
        class(psb_i_base_vect_type), intent(inout)  :: val
        integer(psb_ipk_), intent(out)              :: info
    
        integer(psb_ipk_) :: i, isz
        logical :: done_oacc
    
        info = 0
        if (psb_errstatus_fatal()) return 
    
        done_oacc = .false. 
        select type(virl => irl)
        class is (psb_i_vect_oacc) 
          select type(vval => val)
          class is (psb_i_vect_oacc) 
            if (vval%is_host()) call vval%sync()
            if (virl%is_host()) call virl%sync()
            if (x%is_host())    call x%sync()
            ! Add the OpenACC kernel call here if needed
            call x%set_dev()
            done_oacc = .true.
          end select
        end select
    
        if (.not.done_oacc) then 
          if (irl%is_dev()) call irl%sync()
          if (val%is_dev()) call val%sync()
          call x%ins(n, irl%v, val%v, dupl, info)
        end if
    
        if (info /= 0) then 
          call psb_errpush(info,'i_oacc_ins_v')
          return
        end if
    end subroutine i_oacc_ins_v

    subroutine i_oacc_bld_x(x, this)
        use psb_error_mod
        implicit none 
        integer(psb_ipk_), intent(in)         :: this(:)
        class(psb_i_vect_oacc), intent(inout) :: x
        integer(psb_ipk_) :: info
    
        call psb_realloc(size(this), x%v, info)
        if (info /= 0) then 
            info = psb_err_alloc_request_
            call psb_errpush(info, 'i_oacc_bld_x', i_err = (/size(this), izero, izero, izero, izero/))
        end if
        x%v(:) = this(:)
        call x%set_host()
        call x%sync()
    end subroutine i_oacc_bld_x
    
    subroutine i_oacc_bld_mn(x, n)
        use psb_error_mod
        implicit none 
        integer(psb_mpk_), intent(in)         :: n
        class(psb_i_vect_oacc), intent(inout) :: x
        integer(psb_ipk_) :: info
    
        call x%all(n, info)
        if (info /= 0) then 
            call psb_errpush(info, 'i_oacc_bld_mn', i_err = (/n, n, n, n, n/))
        end if
    end subroutine i_oacc_bld_mn

    subroutine i_oacc_sync(x)
        use psb_error_mod
        implicit none 
        class(psb_i_vect_oacc), intent(inout) :: x
        integer(psb_ipk_) :: n, info
    
        info = 0
        if (x%is_host()) then 
            n = size(x%v)
            if (.not.allocated(x%v)) then 
                write(*, *) 'Incoherent situation : x%v not allocated'
                call psb_realloc(n, x%v, info)
            end if
            if ((n > size(x%v)) .or. (n > x%get_nrows())) then 
                write(*, *) 'Incoherent situation : sizes', n, size(x%v), x%get_nrows()
                call psb_realloc(n, x%v, info)
            end if
            !$acc update device(x%v)
        else if (x%is_dev()) then 
            n = size(x%v)
            if (.not.allocated(x%v)) then 
                write(*, *) 'Incoherent situation : x%v not allocated'
                call psb_realloc(n, x%v, info)
            end if
            !$acc update self(x%v)
        end if
        if (info == 0) call x%set_sync()
        if (info /= 0) then
            info = psb_err_internal_error_
            call psb_errpush(info, 'i_oacc_sync')
        end if
    end subroutine i_oacc_sync
    
    subroutine i_oacc_sync_space(x, info)
        use psb_error_mod
        implicit none 
        class(psb_i_vect_oacc), intent(inout) :: x
        integer(psb_ipk_), intent(out)        :: info 
        integer(psb_ipk_) :: nh, nd
    
        info = 0
        if (x%is_dev()) then 
            nh = size(x%v)
            nd = nh
            if (nh < nd) then 
                call psb_realloc(nd, x%v, info)
            end if
        else  
            nh = size(x%v)
            nd = nh
            if (nh < nd) then 
                call psb_realloc(nd, x%v, info)
            end if
        end if
    end subroutine i_oacc_sync_space

    function i_oacc_get_nrows(x) result(res)
        implicit none
        class(psb_i_vect_oacc), intent(in) :: x
        integer(psb_ipk_) :: res
    
        res = 0
        if (allocated(x%v)) res = size(x%v)
    end function i_oacc_get_nrows
    
    function i_oacc_get_fmt() result(res)
        implicit none
        character(len=5) :: res
        res = 'iOACC'
    end function i_oacc_get_fmt
    
    subroutine i_oacc_all(n, x, info)
        use psb_error_mod
        implicit none
        class(psb_i_vect_oacc), intent(out) :: x
        integer(psb_ipk_), intent(in) :: n
        integer(psb_ipk_), intent(out) :: info
    
        call psb_realloc(n, x%v, info)
        if (info == 0) call x%set_host()
        if (info == 0) call x%sync_space(info)
        if (info /= 0) then 
            info = psb_err_alloc_request_
            call psb_errpush(info, 'i_oacc_all', i_err=(/n, n, n, n, n/))
        end if
    end subroutine i_oacc_all
    
    subroutine i_oacc_zero(x)
        use psb_error_mod
        implicit none
        class(psb_i_vect_oacc), intent(inout) :: x
        ! Ensure zeroing on the GPU side
        call x%set_dev()
        x%v = 0
        !$acc update device(x%v)
    end subroutine i_oacc_zero
    
    subroutine i_oacc_asb_m(n, x, info)
        use psb_error_mod
        use psb_realloc_mod
        implicit none
        integer(psb_ipk_), intent(in) :: n
        class(psb_i_vect_oacc), intent(inout) :: x
        integer(psb_ipk_), intent(out) :: info
        integer(psb_ipk_) :: nh, nd
    
        info = 0
        if (x%is_dev()) then
            nd = size(x%v)
            if (nd < n) then
                call x%sync()
                call x%psb_i_base_vect_type%asb(n, info)
                if (info == psb_success_) call x%sync_space(info)
                call x%set_host()
            end if
        else
            nh = size(x%v)
            if (nh < n) then
                call x%psb_i_base_vect_type%asb(n, info)
                if (info == psb_success_) call x%sync_space(info)
                call x%set_host()
            end if
        end if
    end subroutine i_oacc_asb_m

    subroutine i_oacc_vect_finalize(x)
        use psi_serial_mod
        use psb_realloc_mod
        implicit none
        type(psb_i_vect_oacc), intent(inout) :: x
        integer(psb_ipk_) :: info

        info = 0
        call x%free(info)
    end subroutine i_oacc_vect_finalize

end module psb_i_oacc_vect_mod

    
    
    
    
    