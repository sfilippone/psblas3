module psb_s_oacc_hll_mat_mod
  use iso_c_binding
  use psb_s_mat_mod
  use psb_s_hll_mat_mod
  use psb_s_oacc_vect_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_s_hll_sparse_mat) :: psb_s_oacc_hll_sparse_mat
    integer(psb_ipk_) :: devstate = is_host
  contains
    procedure, nopass  :: get_fmt        => s_oacc_hll_get_fmt
    procedure, pass(a) :: sizeof         => s_oacc_hll_sizeof
    procedure, pass(a) :: is_host        => s_oacc_hll_is_host
    procedure, pass(a) :: is_sync        => s_oacc_hll_is_sync
    procedure, pass(a) :: is_dev         => s_oacc_hll_is_dev
    procedure, pass(a) :: set_host       => s_oacc_hll_set_host
    procedure, pass(a) :: set_sync       => s_oacc_hll_set_sync
    procedure, pass(a) :: set_dev        => s_oacc_hll_set_dev
    procedure, pass(a) :: sync_space     => s_oacc_hll_sync_space
    procedure, pass(a) :: sync           => s_oacc_hll_sync
    procedure, pass(a) :: free           => s_oacc_hll_free
    procedure, pass(a) :: vect_mv        => psb_s_oacc_hll_vect_mv
    procedure, pass(a) :: in_vect_sv     => psb_s_oacc_hll_inner_vect_sv
    procedure, pass(a) :: csmm           => psb_s_oacc_hll_csmm
    procedure, pass(a) :: csmv           => psb_s_oacc_hll_csmv
    procedure, pass(a) :: scals          => psb_s_oacc_hll_scals
    procedure, pass(a) :: scalv          => psb_s_oacc_hll_scal
    procedure, pass(a) :: reallocate_nz  => psb_s_oacc_hll_reallocate_nz
    procedure, pass(a) :: allocate_mnnz  => psb_s_oacc_hll_allocate_mnnz
    procedure, pass(a) :: cp_from_coo    => psb_s_oacc_hll_cp_from_coo
    procedure, pass(a) :: cp_from_fmt    => psb_s_oacc_hll_cp_from_fmt
    procedure, pass(a) :: mv_from_coo    => psb_s_oacc_hll_mv_from_coo
    procedure, pass(a) :: mv_from_fmt    => psb_s_oacc_hll_mv_from_fmt
    procedure, pass(a) :: mold           => psb_s_oacc_hll_mold

  end type psb_s_oacc_hll_sparse_mat

  interface 
    module subroutine psb_s_oacc_hll_mold(a,b,info)
      class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
      class(psb_s_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_hll_mold
  end interface

  interface
    module subroutine psb_s_oacc_hll_cp_from_fmt(a,b,info)
      class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_hll_cp_from_fmt
  end interface

  interface 
    module subroutine psb_s_oacc_hll_mv_from_coo(a,b,info)
      class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_hll_mv_from_coo
  end interface

  interface
    module subroutine psb_s_oacc_hll_mv_from_fmt(a,b,info)
      class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_hll_mv_from_fmt
  end interface

  interface 
    module subroutine psb_s_oacc_hll_vect_mv(alpha, a, x, beta, y, info, trans)
      class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in) :: alpha, beta
      class(psb_s_base_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_oacc_hll_vect_mv
  end interface

  interface
    module subroutine psb_s_oacc_hll_inner_vect_sv(alpha, a, x, beta, y, info, trans)
      class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in) :: alpha, beta
      class(psb_s_base_vect_type), intent(inout) :: x,y
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_oacc_hll_inner_vect_sv
  end interface

  interface
    module subroutine psb_s_oacc_hll_csmm(alpha, a, x, beta, y, info, trans)
      class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in) :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout) :: y(:,:)
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_oacc_hll_csmm
  end interface

  interface
    module subroutine psb_s_oacc_hll_csmv(alpha, a, x, beta, y, info, trans)
      class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in) :: alpha, beta, x(:)
      real(psb_spk_), intent(inout) :: y(:)
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_oacc_hll_csmv
  end interface

  interface
    module subroutine psb_s_oacc_hll_scals(d, a, info)
      class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in) :: d
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_hll_scals
  end interface

  interface 
    module subroutine psb_s_oacc_hll_scal(d,a,info,side)
      class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in) :: d(:)
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: side
    end subroutine psb_s_oacc_hll_scal
  end interface

  interface
    module subroutine psb_s_oacc_hll_reallocate_nz(nz,a)
      class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: nz
    end subroutine psb_s_oacc_hll_reallocate_nz
  end interface

  interface
    module subroutine psb_s_oacc_hll_allocate_mnnz(m,n,a,nz)
      class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: m,n
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_s_oacc_hll_allocate_mnnz
  end interface

  interface 
    module subroutine psb_s_oacc_hll_cp_from_coo(a,b,info)
      class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_hll_cp_from_coo
  end interface

contains 

  subroutine s_oacc_hll_free(a)
    use psb_base_mod
    implicit none 
    class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a
    integer(psb_ipk_) :: info

    if (allocated(a%val)) then
      !$acc exit data delete(a%val)
    end if
    if (allocated(a%ja)) then
      !$acc exit data delete(a%ja)
    end if
    if (allocated(a%irn)) then
      !$acc exit data delete(a%irn)
    end if
    if (allocated(a%idiag)) then
      !$acc exit data delete(a%idiag)
    end if
    if (allocated(a%hkoffs)) then
      !$acc exit data delete(a%hkoffs)
    end if

    call a%psb_s_hll_sparse_mat%free()

    return
  end subroutine s_oacc_hll_free

  function s_oacc_hll_sizeof(a) result(res)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res

    if (a%is_dev()) call a%sync()

    res = 8
    res = res + psb_sizeof_dp * size(a%val)
    res = res + psb_sizeof_ip * size(a%ja)
    res = res + psb_sizeof_ip * size(a%irn)
    res = res + psb_sizeof_ip * size(a%idiag)
    res = res + psb_sizeof_ip * size(a%hkoffs)
  end function s_oacc_hll_sizeof



  function s_oacc_hll_is_host(a) result(res)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_host)
  end function s_oacc_hll_is_host

  function s_oacc_hll_is_sync(a) result(res)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_sync)
  end function s_oacc_hll_is_sync

  function s_oacc_hll_is_dev(a) result(res)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_dev)
  end function s_oacc_hll_is_dev

  subroutine s_oacc_hll_set_host(a)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a

    a%devstate = is_host
  end subroutine s_oacc_hll_set_host

  subroutine s_oacc_hll_set_sync(a)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a

    a%devstate = is_sync
  end subroutine s_oacc_hll_set_sync

  subroutine s_oacc_hll_set_dev(a)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a

    a%devstate = is_dev
  end subroutine s_oacc_hll_set_dev

  function s_oacc_hll_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'HLL_oacc'
  end function s_oacc_hll_get_fmt

  subroutine s_oacc_hll_sync_space(a)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), intent(inout) :: a

    if (allocated(a%val)) then
      call s_oacc_create_dev(a%val)
    end if
    if (allocated(a%ja)) then
      call i_oacc_create_dev(a%ja)
    end if
    if (allocated(a%irn)) then
      call i_oacc_create_dev_scalar(a%irn)
    end if
    if (allocated(a%idiag)) then
      call i_oacc_create_dev_scalar(a%idiag)
    end if
    if (allocated(a%hkoffs)) then
      call i_oacc_create_dev_scalar(a%hkoffs)
    end if

  contains
    subroutine s_oacc_create_dev(v)
      implicit none
      real(psb_spk_), intent(in) :: v(:)
      !$acc enter data copyin(v)
    end subroutine s_oacc_create_dev

    subroutine i_oacc_create_dev(v)
      implicit none
      integer(psb_ipk_), intent(in) :: v(:)
      !$acc enter data copyin(v)
    end subroutine i_oacc_create_dev

    subroutine i_oacc_create_dev_scalar(v)
      implicit none
      integer(psb_ipk_), intent(in) :: v(:)
      !$acc enter data copyin(v)
    end subroutine i_oacc_create_dev_scalar

  end subroutine s_oacc_hll_sync_space


  subroutine s_oacc_hll_sync(a)
    implicit none
    class(psb_s_oacc_hll_sparse_mat), target, intent(in) :: a
    class(psb_s_oacc_hll_sparse_mat), pointer :: tmpa
    integer(psb_ipk_) :: info

    tmpa => a
    if (a%is_dev()) then
      call s_oacc_hll_to_host(a%val)
      call i_oacc_hll_to_host(a%ja)
      call i_oacc_hll_to_host_scalar(a%irn)
      call i_oacc_hll_to_host_scalar(a%idiag)
      call i_oacc_hll_to_host_scalar(a%hkoffs)
    else if (a%is_host()) then
      call s_oacc_hll_to_dev(a%val)
      call i_oacc_hll_to_dev(a%ja)
      call i_oacc_hll_to_dev_scalar(a%irn)
      call i_oacc_hll_to_dev_scalar(a%idiag)
      call i_oacc_hll_to_dev_scalar(a%hkoffs)
    end if
    call tmpa%set_sync()
  end subroutine s_oacc_hll_sync

  subroutine s_oacc_hll_to_host(v)
    implicit none
    real(psb_spk_), intent(in) :: v(:)
    !$acc update self(v)
  end subroutine s_oacc_hll_to_host

  subroutine s_oacc_hll_to_dev(v)
    implicit none
    real(psb_spk_), intent(in) :: v(:)
   !$acc update device(v)
  end subroutine s_oacc_hll_to_dev

  subroutine i_oacc_hll_to_host(v)
    implicit none
    integer(psb_ipk_), intent(in) :: v(:)
    !$acc update self(v)
  end subroutine i_oacc_hll_to_host

  subroutine i_oacc_hll_to_dev(v)
    implicit none
    integer(psb_ipk_), intent(in) :: v(:)
    !$acc update device(v)
  end subroutine i_oacc_hll_to_dev

  subroutine i_oacc_hll_to_host_scalar(v)
    implicit none
    integer(psb_ipk_), intent(in) :: v(:)
    !$acc update self(v)
  end subroutine i_oacc_hll_to_host_scalar

  subroutine i_oacc_hll_to_dev_scalar(v)
    implicit none
    integer(psb_ipk_), intent(in) :: v(:)
    !$acc update device(v)
  end subroutine i_oacc_hll_to_dev_scalar


end module psb_s_oacc_hll_mat_mod
