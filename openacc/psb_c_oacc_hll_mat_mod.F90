module psb_c_oacc_hll_mat_mod
  use iso_c_binding
  use openacc
  use psb_c_mat_mod
  use psb_c_hll_mat_mod
  use psb_c_oacc_vect_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_c_hll_sparse_mat) :: psb_c_oacc_hll_sparse_mat
    integer(psb_ipk_) :: devstate = is_host
  contains
    procedure, nopass  :: get_fmt         => c_oacc_hll_get_fmt
    procedure, pass(a) :: sizeof          => c_oacc_hll_sizeof
    procedure, pass(a) :: is_host         => c_oacc_hll_is_host
    procedure, pass(a) :: is_sync         => c_oacc_hll_is_sync
    procedure, pass(a) :: is_dev          => c_oacc_hll_is_dev
    procedure, pass(a) :: set_host        => c_oacc_hll_set_host
    procedure, pass(a) :: set_sync        => c_oacc_hll_set_sync
    procedure, pass(a) :: set_dev         => c_oacc_hll_set_dev
    procedure, pass(a) :: sync_dev_space  => c_oacc_hll_sync_dev_space
    procedure, pass(a) :: sync            => c_oacc_hll_sync
    procedure, pass(a) :: free_dev_space  => c_oacc_hll_free_dev_space
    procedure, pass(a) :: free            => c_oacc_hll_free
    procedure, pass(a) :: vect_mv         => psb_c_oacc_hll_vect_mv
    procedure, pass(a) :: in_vect_sv      => psb_c_oacc_hll_inner_vect_sv
    procedure, pass(a) :: scals           => psb_c_oacc_hll_scals
    procedure, pass(a) :: scalv           => psb_c_oacc_hll_scal
    procedure, pass(a) :: reallocate_nz   => psb_c_oacc_hll_reallocate_nz
    procedure, pass(a) :: allocate_mnnz   => psb_c_oacc_hll_allocate_mnnz
    procedure, pass(a) :: cp_from_coo     => psb_c_oacc_hll_cp_from_coo
    procedure, pass(a) :: cp_from_fmt     => psb_c_oacc_hll_cp_from_fmt
    procedure, pass(a) :: mv_from_coo     => psb_c_oacc_hll_mv_from_coo
    procedure, pass(a) :: mv_from_fmt     => psb_c_oacc_hll_mv_from_fmt
    procedure, pass(a) :: mold            => psb_c_oacc_hll_mold

  end type psb_c_oacc_hll_sparse_mat

  interface 
    module subroutine psb_c_oacc_hll_mold(a,b,info)
      class(psb_c_oacc_hll_sparse_mat), intent(in) :: a
      class(psb_c_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_c_oacc_hll_mold
  end interface

  interface
    module subroutine psb_c_oacc_hll_cp_from_fmt(a,b,info)
      class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_c_oacc_hll_cp_from_fmt
  end interface

  interface 
    module subroutine psb_c_oacc_hll_mv_from_coo(a,b,info)
      class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_c_oacc_hll_mv_from_coo
  end interface

  interface
    module subroutine psb_c_oacc_hll_mv_from_fmt(a,b,info)
      class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_c_oacc_hll_mv_from_fmt
  end interface

  interface 
    module subroutine psb_c_oacc_hll_vect_mv(alpha, a, x, beta, y, info, trans)
      class(psb_c_oacc_hll_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in) :: alpha, beta
      class(psb_c_base_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: trans
    end subroutine psb_c_oacc_hll_vect_mv
  end interface

  interface
    module subroutine psb_c_oacc_hll_inner_vect_sv(alpha, a, x, beta, y, info, trans)
      class(psb_c_oacc_hll_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in) :: alpha, beta
      class(psb_c_base_vect_type), intent(inout) :: x,y
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: trans
    end subroutine psb_c_oacc_hll_inner_vect_sv
  end interface

  interface
    module subroutine psb_c_oacc_hll_scals(d, a, info)
      class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in) :: d
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_c_oacc_hll_scals
  end interface

  interface 
    module subroutine psb_c_oacc_hll_scal(d,a,info,side)
      class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in) :: d(:)
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: side
    end subroutine psb_c_oacc_hll_scal
  end interface

  interface
    module subroutine psb_c_oacc_hll_reallocate_nz(nz,a)
      class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: nz
    end subroutine psb_c_oacc_hll_reallocate_nz
  end interface

  interface
    module subroutine psb_c_oacc_hll_allocate_mnnz(m,n,a,nz)
      class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: m,n
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_c_oacc_hll_allocate_mnnz
  end interface

  interface 
    module subroutine psb_c_oacc_hll_cp_from_coo(a,b,info)
      class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_c_oacc_hll_cp_from_coo
  end interface

contains 

  subroutine c_oacc_hll_free_dev_space(a)
    use psb_base_mod
    implicit none 
    class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
    integer(psb_ipk_) :: info

    if (allocated(a%val))    call acc_delete_finalize(a%val)
    if (allocated(a%ja))     call acc_delete_finalize(a%ja)
    if (allocated(a%irn))    call acc_delete_finalize(a%irn)
    if (allocated(a%idiag))  call acc_delete_finalize(a%idiag)
    if (allocated(a%hkoffs)) call acc_delete_finalize(a%hkoffs)

    return
  end subroutine c_oacc_hll_free_dev_space

  subroutine c_oacc_hll_free(a)
    use psb_base_mod
    implicit none 
    class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a
    integer(psb_ipk_) :: info

    call a%free_dev_space()
    call a%psb_c_hll_sparse_mat%free()

    return
  end subroutine c_oacc_hll_free

  function c_oacc_hll_sizeof(a) result(res)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res

    if (a%is_dev()) call a%sync()

    res = 8
    res = res + psb_sizeof_dp * size(a%val)
    res = res + psb_sizeof_ip * size(a%ja)
    res = res + psb_sizeof_ip * size(a%irn)
    res = res + psb_sizeof_ip * size(a%idiag)
    res = res + psb_sizeof_ip * size(a%hkoffs)
  end function c_oacc_hll_sizeof



  function c_oacc_hll_is_host(a) result(res)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_host)
  end function c_oacc_hll_is_host

  function c_oacc_hll_is_sync(a) result(res)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_sync)
  end function c_oacc_hll_is_sync

  function c_oacc_hll_is_dev(a) result(res)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_dev)
  end function c_oacc_hll_is_dev

  subroutine c_oacc_hll_set_host(a)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a

    a%devstate = is_host
  end subroutine c_oacc_hll_set_host

  subroutine c_oacc_hll_set_sync(a)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a

    a%devstate = is_sync
  end subroutine c_oacc_hll_set_sync

  subroutine c_oacc_hll_set_dev(a)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a

    a%devstate = is_dev
  end subroutine c_oacc_hll_set_dev

  function c_oacc_hll_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'HLLOA'
  end function c_oacc_hll_get_fmt

  subroutine c_oacc_hll_sync_dev_space(a)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), intent(inout) :: a

    if (allocated(a%val))    call acc_create(a%val)
    if (allocated(a%ja))     call acc_create(a%ja)
    if (allocated(a%irn))    call acc_create(a%irn)
    if (allocated(a%idiag))  call acc_create(a%idiag)
    if (allocated(a%hkoffs)) call acc_create(a%hkoffs)
  end subroutine c_oacc_hll_sync_dev_space


  subroutine c_oacc_hll_sync(a)
    implicit none
    class(psb_c_oacc_hll_sparse_mat), target, intent(in) :: a
    class(psb_c_oacc_hll_sparse_mat), pointer :: tmpa
    integer(psb_ipk_) :: info

    tmpa  => a
    if (a%is_dev()) then
      call acc_update_self(a%val)
      call acc_update_self(a%ja)
      call acc_update_self(a%irn)
      call acc_update_self(a%idiag)
      call acc_update_self(a%hkoffs)
    else if (a%is_host()) then
      call acc_update_device(a%val)
      call acc_update_device(a%ja)
      call acc_update_device(a%irn)
      call acc_update_device(a%idiag)
      call acc_update_device(a%hkoffs)
    end if
    call tmpa%set_sync()
  end subroutine c_oacc_hll_sync

end module psb_c_oacc_hll_mat_mod
