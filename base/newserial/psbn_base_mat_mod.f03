module psbn_base_mat_mod
  
  use psb_const_mod 

  integer, parameter :: psbn_spmat_null_=0, psbn_spmat_bld_=1
  integer, parameter :: psbn_spmat_asb_=2, psbn_spmat_upd_=4

  integer, parameter :: psbn_ireg_flgs_=10, psbn_ip2_=0
  integer, parameter :: psbn_iflag_=2, psbn_ichk_=3
  integer, parameter :: psbn_nnzt_=4, psbn_zero_=5,psbn_ipc_=6
  ! Duplicate coefficients handling
  ! These are usually set while calling spcnv as one of its
  ! optional arugments.
  integer, parameter :: psbn_dupl_ovwrt_ = 0
  integer, parameter :: psbn_dupl_add_   = 1
  integer, parameter :: psbn_dupl_err_   = 2
  integer, parameter :: psbn_dupl_def_   = psbn_dupl_ovwrt_
  ! Matrix update mode
  integer, parameter :: psbn_upd_srch_   = 98764
  integer, parameter :: psbn_upd_perm_   = 98765
  integer, parameter :: psbn_upd_dflt_   = psbn_upd_srch_
  integer, parameter :: psbn_maxjdrows_=8, psbn_minjdrows_=4
  integer, parameter :: psbn_dbleint_=2
  character(len=5)   :: psbn_fidef_='CSR'


  type  :: psbn_base_sparse_mat
    integer              :: m, n
    integer, private     :: state, duplicate 
    logical, private     :: triangle, unitd, upper, sorted
  contains 
    procedure, pass(a) :: get_nrows
    procedure, pass(a) :: get_ncols
    procedure, pass(a) :: get_nzeros
    procedure, pass(a) :: get_size
    procedure, pass(a) :: get_state
    procedure, pass(a) :: get_dupl
    procedure, pass(a) :: is_null
    procedure, pass(a) :: is_bld
    procedure, pass(a) :: is_upd
    procedure, pass(a) :: is_asb
    procedure, pass(a) :: is_sorted
    procedure, pass(a) :: is_upper
    procedure, pass(a) :: is_lower
    procedure, pass(a) :: is_triangle
    procedure, pass(a) :: is_unit
    procedure, pass(a) :: get_neigh
    procedure, pass(a) :: allocate_mn
    procedure, pass(a) :: allocate_mnnz
    procedure, pass(a) :: reallocate_nz
    procedure, pass(a) :: free
    generic,   public  :: allocate => allocate_mn, allocate_mnnz
    generic,   public  :: reallocate => reallocate_nz
    
  end type psbn_base_sparse_mat
  
contains
 
  function get_dupl(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%duplicate
  end function get_dupl
 
 
  function get_state(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%state
  end function get_state
 
  function get_nrows(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%m
  end function get_nrows

  function get_ncols(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%n
  end function get_ncols

  function is_triangle(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%triangle
  end function is_triangle

  function is_unit(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%unitd
  end function is_unit

  function is_upper(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%upper
  end function is_upper

  function is_lower(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = .not.a%upper
  end function is_lower

  function is_null(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_null_)
  end function is_null

  function is_bld(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_bld_)
  end function is_bld

  function is_upd(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_upd_)
  end function is_upd

  function is_asb(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_asb_)
  end function is_asb

  function is_sorted(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%sorted
  end function is_sorted


  function get_nzeros(a) result(res)
    use psb_error_mod
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act
    character(len=20)  :: name='base_get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    res = -1
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end function get_nzeros

  function get_size(a) result(res)
    use psb_error_mod
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act
    character(len=20)  :: name='get_size'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    res = -1
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end function get_size


  subroutine get_neigh(a,idx,neigh,n,info,lev)
    use psb_error_mod
    class(psbn_base_sparse_mat), intent(in) :: a   
    integer, intent(in)                :: idx 
    integer, intent(out)               :: n   
    integer, allocatable, intent(out)  :: neigh(:)
    integer, intent(out)               :: info
    integer, optional, intent(in)      :: lev 
    
    Integer :: err_act
    character(len=20)  :: name='get_neigh'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 700
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name)

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine get_neigh

  subroutine  allocate_mn(m,n,a) 
    use psb_error_mod
    integer, intent(in) :: m,n
    class(psbn_base_sparse_mat), intent(inout) :: a

    Integer :: err_act
    character(len=20)  :: name='allocate_mn'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine allocate_mn

  subroutine  allocate_mnnz(m,n,nz,a) 
    use psb_error_mod
    integer, intent(in) :: m,n,nz
    class(psbn_base_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='allocate_mnz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine allocate_mnnz

  subroutine  reallocate_nz(nz,a) 
    use psb_error_mod
    integer, intent(in) :: nz
    class(psbn_base_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine reallocate_nz

  subroutine  free(a) 
    use psb_error_mod
    class(psbn_base_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='free'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    ! This is the base version. If we get here
    ! it means the derived class is incomplete,
    ! so we throw an error.
    call psb_errpush(700,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine free

end module psbn_base_mat_mod

