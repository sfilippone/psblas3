!
! psb_comm_mod - communication handle module
!
module psb_comm_schemes_mod
  use psb_const_mod
  implicit none

  ! Communication type enumeration (keeps compatibility with integer selectors)
  enum, bind(c)
    enumerator psb_comm_unknown_
    enumerator psb_comm_isend_irecv_
    enumerator psb_comm_ineighbor_alltoallv_
    enumerator psb_comm_persistent_ineighbor_alltoallv_
    enumerator psb_comm_rma_pull_
    enumerator psb_comm_rma_push_
  end enum

  enum, bind(c)
    enumerator psb_comm_status_unknown_
    enumerator psb_comm_status_start_
    enumerator psb_comm_status_wait_
    enumerator psb_comm_status_sync_ ! Used in order to exchange data in a synchronous way (Start and recv in the same call)
  end enum


  ! --- comm handle type  ---
  type, abstract :: psb_comm_handle_type
    integer(psb_ipk_) :: id = -1
    integer(psb_ipk_) :: comm_type = psb_comm_unknown_
    integer(psb_ipk_) :: swap_status = psb_comm_status_unknown_
  contains
    procedure(psb_comm_set), deferred             :: init
    procedure(psb_comm_free), deferred            :: free
    procedure(psb_comm_set_swap_status), deferred :: set_swap_status
    procedure(psb_comm_get_swap_status), deferred :: get_swap_status
  end type psb_comm_handle_type

  ! --- abstract interfaces ---
  abstract interface
    subroutine psb_comm_set(this, info)
      import :: psb_ipk_, psb_comm_handle_type
      class(psb_comm_handle_type), intent(inout)  :: this
      integer(psb_ipk_), intent(out)              :: info
    end subroutine

    subroutine psb_comm_free(this, info)
      import :: psb_ipk_, psb_comm_handle_type
      class(psb_comm_handle_type), intent(inout)  :: this
      integer(psb_ipk_), intent(out)              :: info
    end subroutine

    subroutine psb_comm_set_swap_status(this, flag, info)
      import :: psb_ipk_, psb_comm_handle_type
      class(psb_comm_handle_type), intent(inout)  :: this
      integer(psb_ipk_), intent(in)               :: flag
      integer(psb_ipk_), intent(out)              :: info
    end subroutine

    subroutine psb_comm_get_swap_status(this, flag, info)
      import :: psb_ipk_, psb_comm_handle_type
      class(psb_comm_handle_type), intent(in)     :: this
      integer(psb_ipk_), intent(out)              :: flag
      integer(psb_ipk_), intent(out)              :: info
    end subroutine
  end interface

contains

end module psb_comm_schemes_mod
