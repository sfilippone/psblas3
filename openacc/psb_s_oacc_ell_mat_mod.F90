module psb_s_oacc_ell_mat_mod
  use iso_c_binding
  use openacc
  use psb_s_mat_mod
  use psb_s_ell_mat_mod
  use psb_s_oacc_vect_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 

  type, extends(psb_s_ell_sparse_mat) :: psb_s_oacc_ell_sparse_mat
    integer(psb_ipk_) :: devstate = is_host
  contains
    procedure, nopass  :: get_fmt         => s_oacc_ell_get_fmt
    procedure, pass(a) :: sizeof          => s_oacc_ell_sizeof
    procedure, pass(a) :: is_host         => s_oacc_ell_is_host
    procedure, pass(a) :: is_sync         => s_oacc_ell_is_sync
    procedure, pass(a) :: is_dev          => s_oacc_ell_is_dev
    procedure, pass(a) :: set_host        => s_oacc_ell_set_host
    procedure, pass(a) :: set_sync        => s_oacc_ell_set_sync
    procedure, pass(a) :: set_dev         => s_oacc_ell_set_dev
    procedure, pass(a) :: sync_dev_space  => s_oacc_ell_sync_dev_space
    procedure, pass(a) :: sync            => s_oacc_ell_sync
    procedure, pass(a) :: free_dev_space  => s_oacc_ell_free_dev_space
    procedure, pass(a) :: free            => s_oacc_ell_free
    procedure, pass(a) :: vect_mv         => psb_s_oacc_ell_vect_mv
    procedure, pass(a) :: in_vect_sv      => psb_s_oacc_ell_inner_vect_sv
    procedure, pass(a) :: scals           => psb_s_oacc_ell_scals
    procedure, pass(a) :: scalv           => psb_s_oacc_ell_scal
    procedure, pass(a) :: reallocate_nz   => psb_s_oacc_ell_reallocate_nz
    procedure, pass(a) :: allocate_mnnz   => psb_s_oacc_ell_allocate_mnnz
    procedure, pass(a) :: cp_from_coo     => psb_s_oacc_ell_cp_from_coo
    procedure, pass(a) :: cp_from_fmt     => psb_s_oacc_ell_cp_from_fmt
    procedure, pass(a) :: mv_from_coo     => psb_s_oacc_ell_mv_from_coo
    procedure, pass(a) :: mv_from_fmt     => psb_s_oacc_ell_mv_from_fmt
    procedure, pass(a) :: mold            => psb_s_oacc_ell_mold

  end type psb_s_oacc_ell_sparse_mat

  interface 
    module subroutine psb_s_oacc_ell_mold(a,b,info)
      class(psb_s_oacc_ell_sparse_mat), intent(in) :: a
      class(psb_s_base_sparse_mat), intent(inout), allocatable :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_ell_mold
  end interface

  interface
    module subroutine psb_s_oacc_ell_cp_from_fmt(a,b,info)
      class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_ell_cp_from_fmt
  end interface

  interface 
    module subroutine psb_s_oacc_ell_mv_from_coo(a,b,info)
      class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_ell_mv_from_coo
  end interface

  interface
    module subroutine psb_s_oacc_ell_mv_from_fmt(a,b,info)
      class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_ell_mv_from_fmt
  end interface

  interface 
    module subroutine psb_s_oacc_ell_vect_mv(alpha, a, x, beta, y, info, trans)
      class(psb_s_oacc_ell_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in) :: alpha, beta
      class(psb_s_base_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_oacc_ell_vect_mv
  end interface

  interface
    module subroutine psb_s_oacc_ell_inner_vect_sv(alpha, a, x, beta, y, info, trans)
      class(psb_s_oacc_ell_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in) :: alpha, beta
      class(psb_s_base_vect_type), intent(inout) :: x,y
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: trans
    end subroutine psb_s_oacc_ell_inner_vect_sv
  end interface

  interface
    module subroutine psb_s_oacc_ell_scals(d, a, info)
      class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in) :: d
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_ell_scals
  end interface

  interface 
    module subroutine psb_s_oacc_ell_scal(d,a,info,side)
      class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in) :: d(:)
      integer(psb_ipk_), intent(out) :: info
      character, optional, intent(in) :: side
    end subroutine psb_s_oacc_ell_scal
  end interface

  interface
    module subroutine psb_s_oacc_ell_reallocate_nz(nz,a)
      class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: nz
    end subroutine psb_s_oacc_ell_reallocate_nz
  end interface

  interface
    module subroutine psb_s_oacc_ell_allocate_mnnz(m,n,a,nz)
      class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
      integer(psb_ipk_), intent(in) :: m,n
      integer(psb_ipk_), intent(in), optional :: nz
    end subroutine psb_s_oacc_ell_allocate_mnnz
  end interface

  interface 
    module subroutine psb_s_oacc_ell_cp_from_coo(a,b,info)
      class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(in) :: b
      integer(psb_ipk_), intent(out) :: info
    end subroutine psb_s_oacc_ell_cp_from_coo
  end interface

contains

  subroutine s_oacc_ell_free_dev_space(a)
    use psb_base_mod
    implicit none 
    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
    integer(psb_ipk_) :: info

    !
    ! Note: at least on GNU, if an array is allocated
    !       but with size 0, then CREATE,UPDATE and DELETE
    !       will fail
    !
    if (psb_size(a%val)>0)   call acc_delete_finalize(a%val)
    if (psb_size(a%ja)>0)    call acc_delete_finalize(a%ja)
    if (psb_size(a%irn)>0)   call acc_delete_finalize(a%irn)
    if (psb_size(a%idiag)>0) call acc_delete_finalize(a%idiag)
    return
  end subroutine s_oacc_ell_free_dev_space

  subroutine s_oacc_ell_free(a)
    use psb_base_mod
    implicit none 
    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a
    integer(psb_ipk_) :: info

    call a%free_dev_space()
    call a%psb_s_ell_sparse_mat%free()

    return
  end subroutine s_oacc_ell_free

  function s_oacc_ell_sizeof(a) result(res)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), intent(in) :: a
    integer(psb_epk_) :: res

    if (a%is_dev()) call a%sync()

    res = 8
    res = res + psb_sizeof_dp * size(a%val)
    res = res + psb_sizeof_ip * size(a%ja)
    res = res + psb_sizeof_ip * size(a%irn)
    res = res + psb_sizeof_ip * size(a%idiag)

  end function s_oacc_ell_sizeof

  subroutine s_oacc_ell_sync_dev_space(a)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a

    !
    ! Note: at least on GNU, if an array is allocated
    !       but with size 0, then CREATE,UPDATE and DELETE
    !       will fail
    !
    if (psb_size(a%val)>0)   call acc_create(a%val)
    if (psb_size(a%ja)>0)    call acc_create(a%ja)
    if (psb_size(a%irn)>0)   call acc_create(a%irn)
    if (psb_size(a%idiag)>0) call acc_create(a%idiag)
  end subroutine s_oacc_ell_sync_dev_space

  function s_oacc_ell_is_host(a) result(res)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_host)
  end function s_oacc_ell_is_host

  function s_oacc_ell_is_sync(a) result(res)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_sync)
  end function s_oacc_ell_is_sync

  function s_oacc_ell_is_dev(a) result(res)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), intent(in) :: a
    logical :: res

    res = (a%devstate == is_dev)
  end function s_oacc_ell_is_dev

  subroutine s_oacc_ell_set_host(a)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a

    a%devstate = is_host
  end subroutine s_oacc_ell_set_host

  subroutine s_oacc_ell_set_sync(a)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a

    a%devstate = is_sync
  end subroutine s_oacc_ell_set_sync

  subroutine s_oacc_ell_set_dev(a)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), intent(inout) :: a

    a%devstate = is_dev
  end subroutine s_oacc_ell_set_dev

  function s_oacc_ell_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'ELLOA'
  end function s_oacc_ell_get_fmt

  subroutine s_oacc_ell_sync(a)
    implicit none
    class(psb_s_oacc_ell_sparse_mat), target, intent(in) :: a
    class(psb_s_oacc_ell_sparse_mat), pointer :: tmpa
    integer(psb_ipk_) :: info

    tmpa  => a
    !
    ! Note: at least on GNU, if an array is allocated
    !       but with size 0, then CREATE,UPDATE and DELETE
    !       will fail
    !
    if (a%is_dev()) then
      if (psb_size(a%val)>0)   call acc_update_self(a%val)
      if (psb_size(a%ja)>0)    call acc_update_self(a%ja)
      if (psb_size(a%irn)>0)   call acc_update_self(a%irn)
      if (psb_size(a%idiag)>0) call acc_update_self(a%idiag)
    else if (a%is_host()) then
      if (psb_size(a%val)>0)         call acc_update_device(a%val)
      if (psb_size(a%ja)>0)          call acc_update_device(a%ja)
      if (psb_size(a%irn)>0)         call acc_update_device(a%irn)
      if (psb_size(a%idiag)>0)       call acc_update_device(a%idiag)
    end if
    call tmpa%set_sync()
  end subroutine s_oacc_ell_sync

end module psb_s_oacc_ell_mat_mod
