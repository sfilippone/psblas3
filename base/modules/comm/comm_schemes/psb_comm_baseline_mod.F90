module psb_comm_baseline_mod
  use psb_comm_schemes_mod, only: psb_comm_handle_type, psb_comm_isend_irecv_
  use psb_const_mod
  implicit none

  type, extends(psb_comm_handle_type) :: psb_comm_baseline_handle
    ! MPI request IDs for Isend/Irecv (dimension: num_neighbors x 2)
    ! First column: send requests, second column: recv requests
    integer(psb_ipk_), allocatable :: comid(:,:)
  contains
    procedure, pass :: init => psb_comm_baseline_init
    procedure, pass :: free => psb_comm_baseline_free
    procedure, pass :: set_swap_status => psb_comm_baseline_set_swap_status
    procedure, pass :: get_swap_status => psb_comm_baseline_get_swap_status
  end type psb_comm_baseline_handle

contains

  subroutine psb_comm_baseline_init(this, info)
    implicit none
    class(psb_comm_baseline_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info
    info = 0
    this%comm_type = psb_comm_isend_irecv_
    this%id = 0
    this%swap_status = 0
  end subroutine psb_comm_baseline_init

  subroutine psb_comm_baseline_free(this, info)
    class(psb_comm_baseline_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info
    info = 0
    ! Free MPI resources (comid)
    if (allocated(this%comid)) deallocate(this%comid, stat=info)
  end subroutine psb_comm_baseline_free

  subroutine psb_comm_baseline_set_swap_status(this, flag, info)
    class(psb_comm_baseline_handle), intent(inout) :: this
    integer(psb_ipk_), intent(in) :: flag
    integer(psb_ipk_), intent(out) :: info
    info = 0
    this%swap_status = flag
  end subroutine psb_comm_baseline_set_swap_status

  subroutine psb_comm_baseline_get_swap_status(this, flag, info)
    class(psb_comm_baseline_handle), intent(in) :: this
    integer(psb_ipk_), intent(out) :: flag
    integer(psb_ipk_), intent(out) :: info
    info = 0
    flag = this%swap_status
  end subroutine psb_comm_baseline_get_swap_status

  ! Allocate comid array for num_neighbors
  subroutine psb_comm_baseline_alloc_comid(this, n, info)
    implicit none
    class(psb_comm_baseline_handle), intent(inout) :: this
    integer(psb_ipk_), intent(in) :: n
    integer(psb_ipk_), intent(out) :: info
    allocate(this%comid(n, 2_psb_ipk_), stat=info)
  end subroutine psb_comm_baseline_alloc_comid

end module psb_comm_baseline_mod
