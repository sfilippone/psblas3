module psbn_base_mat_mod
  
  use psb_const_mod 

!!$  integer, parameter :: psbn_invalid_ = -1 
!!$  integer, parameter :: psbn_spmat_null_=0, psbn_spmat_bld_=1
!!$  integer, parameter :: psbn_spmat_asb_=2, psbn_spmat_upd_=4
!!$
!!$  integer, parameter :: psbn_ireg_flgs_=10, psbn_ip2_=0
!!$  integer, parameter :: psbn_iflag_=2, psbn_ichk_=3
!!$  integer, parameter :: psbn_nnzt_=4, psbn_zero_=5,psbn_ipc_=6
!!$  ! Duplicate coefficients handling
!!$  ! These are usually set while calling spcnv as one of its
!!$  ! optional arugments.
!!$  integer, parameter :: psbn_dupl_ovwrt_ = 0
!!$  integer, parameter :: psbn_dupl_add_   = 1
!!$  integer, parameter :: psbn_dupl_err_   = 2
!!$  integer, parameter :: psbn_dupl_def_   = psbn_dupl_ovwrt_
!!$  ! Matrix update mode
!!$  integer, parameter :: psbn_upd_srch_   = 98764
!!$  integer, parameter :: psbn_upd_perm_   = 98765
!!$  integer, parameter :: psbn_upd_dflt_   = psbn_upd_srch_
!!$  integer, parameter :: psbn_maxjdrows_=8, psbn_minjdrows_=4
!!$  integer, parameter :: psbn_dbleint_=2


  type  :: psbn_base_sparse_mat
    integer              :: m, n
    integer, private     :: state, duplicate 
    logical, private     :: triangle, unitd, upper, sorted
  contains 

    ! ====================================
    !
    ! Getters 
    !
    !
    ! ====================================
    procedure, pass(a) :: get_nrows
    procedure, pass(a) :: get_ncols
    procedure, pass(a) :: get_nzeros
    procedure, pass(a) :: get_size
    procedure, pass(a) :: get_state
    procedure, pass(a) :: get_dupl
    procedure, pass(a) :: get_fmt
    procedure, pass(a) :: is_null
    procedure, pass(a) :: is_bld
    procedure, pass(a) :: is_upd
    procedure, pass(a) :: is_asb
    procedure, pass(a) :: is_sorted
    procedure, pass(a) :: is_upper
    procedure, pass(a) :: is_lower
    procedure, pass(a) :: is_triangle
    procedure, pass(a) :: is_unit
    
    ! ====================================
    !
    ! Setters 
    !
    ! ====================================
    procedure, pass(a) :: set_nrows
    procedure, pass(a) :: set_ncols
    procedure, pass(a) :: set_dupl
    procedure, pass(a) :: set_state
    procedure, pass(a) :: set_null
    procedure, pass(a) :: set_bld
    procedure, pass(a) :: set_upd
    procedure, pass(a) :: set_asb
    procedure, pass(a) :: set_sorted
    procedure, pass(a) :: set_upper
    procedure, pass(a) :: set_lower
    procedure, pass(a) :: set_triangle
    procedure, pass(a) :: set_unit



    ! ====================================
    !
    ! Data management
    !
    ! ====================================  
    procedure, pass(a) :: get_neigh
    procedure, pass(a) :: allocate_mnnz
    procedure, pass(a) :: reallocate_nz
    procedure, pass(a) :: free
    procedure, pass(a) :: trim
    generic,   public  :: allocate => allocate_mnnz
    generic,   public  :: reallocate => reallocate_nz
    procedure, pass(a) :: print => sparse_print
    procedure, pass(a) :: sizeof

  end type psbn_base_sparse_mat

  private :: set_nrows, set_ncols, set_dupl, set_state, &
       & set_null, set_bld, set_upd, set_asb, set_sorted, set_upper, &
       & set_lower, set_triangle, set_unit, get_nrows, get_ncols, &
       & get_nzeros, get_size, get_state, get_dupl, is_null, is_bld, &
       & is_upd, is_asb, is_sorted, is_upper, is_lower, is_triangle, &
       & is_unit, get_neigh, allocate_mn, allocate_mnnz, reallocate_nz, &
       & free, sparse_print, get_fmt, trim, sizeof
  

contains

  
  function sizeof(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8
  end function sizeof
 
 

  function get_fmt(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'NULL'
  end function get_fmt
  
  function get_dupl(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%duplicate
  end function get_dupl
 
 
  function get_state(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%state
  end function get_state
 
  function get_nrows(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%m
  end function get_nrows

  function get_ncols(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%n
  end function get_ncols

 
  subroutine  set_nrows(m,a) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: m
    a%m = m
  end subroutine set_nrows

  subroutine  set_ncols(n,a) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    a%n = n
  end subroutine set_ncols


  subroutine  set_state(n,a) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    a%state = n
  end subroutine set_state


  subroutine  set_dupl(n,a) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    a%duplicate = n
  end subroutine set_dupl

  subroutine  set_null(a) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a

    a%state = psbn_spmat_null_
  end subroutine set_null

  subroutine  set_bld(a) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a

    a%state = psbn_spmat_bld_
  end subroutine set_bld

  subroutine  set_upd(a) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a

    a%state = psbn_spmat_upd_
  end subroutine set_upd

  subroutine  set_asb(a) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a

    a%state = psbn_spmat_asb_
  end subroutine set_asb

  subroutine set_sorted(a,val) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%sorted = val
    else
      a%sorted = .true.
    end if
  end subroutine set_sorted

  subroutine set_triangle(a,val) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%triangle = val
    else
      a%triangle = .true.
    end if
  end subroutine set_triangle

  subroutine set_unit(a,val) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%unitd = val
    else
      a%unitd = .true.
    end if
  end subroutine set_unit

  subroutine set_lower(a,val) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%upper = .not.val
    else
      a%upper = .false.
    end if
  end subroutine set_lower

  subroutine set_upper(a,val) 
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%upper = val
    else
      a%upper = .true.
    end if
  end subroutine set_upper

  function is_triangle(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%triangle
  end function is_triangle

  function is_unit(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%unitd
  end function is_unit

  function is_upper(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%upper
  end function is_upper

  function is_lower(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = .not.a%upper
  end function is_lower

  function is_null(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_null_)
  end function is_null

  function is_bld(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_bld_)
  end function is_bld

  function is_upd(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_upd_)
  end function is_upd

  function is_asb(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_asb_)
  end function is_asb

  function is_sorted(a) result(res)
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%sorted
  end function is_sorted


  function get_nzeros(a) result(res)
    use psb_error_mod
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act
    character(len=20)  :: name='base_get_nzeros'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    res = -1
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name,a_err=a%get_fmt())
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end function get_nzeros

  function get_size(a) result(res)
    use psb_error_mod
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act
    character(len=20)  :: name='get_size'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    res = -1
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name,a_err=a%get_fmt())
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end function get_size

  subroutine sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
    use psb_error_mod
    implicit none 

    integer, intent(in)               :: iout
    class(psbn_base_sparse_mat), intent(in) :: a   
    integer, intent(in), optional     :: iv(:)
    integer, intent(in), optional     :: eirs,eics
    character(len=*), optional        :: head
    integer, intent(in), optional     :: ivr(:), ivc(:)

    Integer :: err_act, info
    character(len=20)  :: name='sparse_print'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    info = 700
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name,a_err=a%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine sparse_print


  subroutine get_neigh(a,idx,neigh,n,info,lev)
    use psb_error_mod
    implicit none 
    class(psbn_base_sparse_mat), intent(in) :: a   
    integer, intent(in)                :: idx 
    integer, intent(out)               :: n   
    integer, allocatable, intent(out)  :: neigh(:)
    integer, intent(out)               :: info
    integer, optional, intent(in)      :: lev 
    
    Integer :: err_act
    character(len=20)  :: name='get_neigh'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    info = 700
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name,a_err=a%get_fmt())

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine get_neigh

  subroutine  allocate_mnnz(m,n,a,nz) 
    use psb_error_mod
    implicit none 
    integer, intent(in) :: m,n
    class(psbn_base_sparse_mat), intent(inout) :: a
    integer, intent(in), optional  :: nz
    Integer :: err_act
    character(len=20)  :: name='allocate_mnz'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name,a_err=a%get_fmt())
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine allocate_mnnz

  subroutine  reallocate_nz(nz,a) 
    use psb_error_mod
    implicit none 
    integer, intent(in) :: nz
    class(psbn_base_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name,a_err=a%get_fmt())
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine reallocate_nz

  subroutine  free(a) 
    use psb_error_mod
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='free'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name,a_err=a%get_fmt())
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine free

  subroutine  trim(a) 
    use psb_error_mod
    implicit none 
    class(psbn_base_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='trim'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name,a_err=a%get_fmt())
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine trim

end module psbn_base_mat_mod

