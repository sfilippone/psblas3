module psbn_base_mat_mod

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
  ! Mark a COO matrix with sorted entries.
  integer, parameter :: psbn_isrtdcoo_   = 98761
  integer, parameter :: psbn_maxjdrows_=8, psbn_minjdrows_=4
  integer, parameter :: psbn_dbleint_=2
  character(len=5)   :: psbn_fidef_='CSR'


  type  :: psbn_base_sparse_mat
    integer              :: m, n
    integer, private     :: state 
    logical, private     :: triangle, unitd, upper
  contains 
    procedure, pass(a) :: base_get_nrows
    procedure, pass(a) :: base_get_ncols
    procedure, pass(a) :: base_get_nzeros
    procedure, pass(a) :: base_get_size
    procedure, pass(a) :: base_get_state
    procedure, pass(a) :: base_is_bld
    procedure, pass(a) :: base_is_upd
    procedure, pass(a) :: base_is_asb
    procedure, pass(a) :: base_is_upper
    procedure, pass(a) :: base_is_lower
    procedure, pass(a) :: base_is_triangle
    procedure, pass(a) :: base_is_unit
    procedure, pass(a) :: base_get_neigh
    procedure, pass(a) :: base_allocate_mn
    procedure, pass(a) :: base_allocate_mnnz
    procedure, pass(a) :: base_free
    generic, public    :: allocate => base_allocate_mn, base_allocate_mnnz
    generic, public    :: get_nrows => base_get_nrows
    generic, public    :: get_ncols => base_get_ncols
    generic, public    :: get_nzeros => base_get_nzeros
    generic, public    :: get_size  => base_get_size
    generic, public    :: get_state => base_get_state
    generic, public    :: is_triangle => base_is_triangle
    generic, public    :: is_unit => base_is_unit
    generic, public    :: is_upper => base_is_upper
    generic, public    :: is_lower => base_is_lower
    generic, public    :: is_bld => base_is_bld
    generic, public    :: is_upd => base_is_upd
    generic, public    :: is_asb => base_is_asb
    generic, public    :: get_neigh => base_get_neigh
    generic, public    :: free => base_free
    
  end type psbn_base_sparse_mat
  
contains
 
  function base_get_state(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%state
  end function base_get_state
 
  function base_get_nrows(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%m
  end function base_get_nrows

  function base_get_ncols(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%n
  end function base_get_ncols

  function base_is_triangle(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%triangle
  end function base_is_triangle

  function base_is_unit(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%unitd
  end function base_is_unit

  function base_is_upper(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%upper
  end function base_is_upper

  function base_is_lower(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = .not.a%upper
  end function base_is_lower

  function base_is_bld(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_bld_)
  end function base_is_bld

  function base_is_upd(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_upd_)
  end function base_is_upd

  function base_is_asb(a) result(res)
    class(psbn_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psbn_spmat_asb_)
  end function base_is_asb


  function base_get_nzeros(a) result(res)
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

  end function base_get_nzeros

  function base_get_size(a) result(res)
    use psb_error_mod
    class(psbn_base_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act
    character(len=20)  :: name='base_get_size'
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

  end function base_get_size


  subroutine base_get_neigh(a,idx,neigh,n,info,lev)
    use psb_error_mod
    class(psbn_base_sparse_mat), intent(in) :: a   
    integer, intent(in)                :: idx 
    integer, intent(out)               :: n   
    integer, allocatable, intent(out)  :: neigh(:)
    integer, intent(out)               :: info
    integer, optional, intent(in)      :: lev 
    
    Integer :: err_act
    character(len=20)  :: name='base_get_neigh'
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

  end subroutine base_get_neigh

  subroutine  base_allocate_mn(m,n,a) 
    use psb_error_mod
    integer, intent(in) :: m,n
    class(psbn_base_sparse_mat), intent(inout) :: a

    Integer :: err_act
    character(len=20)  :: name='base_allocate_mn'
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

  end subroutine base_allocate_mn

  subroutine  base_allocate_mnnz(m,n,nz,a) 
    use psb_error_mod
    integer, intent(in) :: m,n,nz
    class(psbn_base_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='base_allocate_mnz'
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

  end subroutine base_allocate_mnnz

  subroutine  base_free(a) 
    use psb_error_mod
    class(psbn_base_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='base_free'
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

  end subroutine base_free

end module psbn_base_mat_mod

