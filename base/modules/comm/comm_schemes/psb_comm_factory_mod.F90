module psb_comm_factory_mod
  use psb_const_mod
  use psb_comm_schemes_mod, only: psb_comm_handle_type, psb_comm_isend_irecv_, &
    & psb_comm_ineighbor_alltoallv_, psb_comm_persistent_ineighbor_alltoallv_, &
    & psb_comm_unknown_, psb_comm_status_start_, psb_comm_status_wait_, &
    & psb_comm_status_sync_, psb_comm_status_unknown_
  use psb_comm_baseline_mod, only: psb_comm_baseline_handle
  use psb_comm_neighbor_impl_mod, only: psb_comm_neighbor_handle
  implicit none

contains

  ! Allocatable-based factory routines (preferred names)
  subroutine psb_comm_init(comm_type, handle, info)
    implicit none
    integer(psb_ipk_), intent(in) :: comm_type
    class(psb_comm_handle_type), allocatable, intent(inout) :: handle
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: old_id, old_swap_status

    info = 0
    old_id = 0
    old_swap_status = psb_comm_status_unknown_

    if (allocated(handle)) then
      old_id = handle%id
      old_swap_status = handle%swap_status

      if (handle%comm_type == comm_type) then
        call handle%free(info)
        if (info /= 0) return
        call handle%init(info)
        if (info /= 0) return
        handle%id = old_id
        handle%swap_status = old_swap_status
        select type(h => handle)
        type is(psb_comm_neighbor_handle)
          h%comm_type = comm_type
          h%use_persistent_buffers = (comm_type == psb_comm_persistent_ineighbor_alltoallv_)
        class default
          ! nothing else to configure
        end select
        return
      else
        call psb_comm_free(handle, info)
        if (info /= 0) return
      end if
    end if

    select case(comm_type)
    case(psb_comm_ineighbor_alltoallv_, psb_comm_persistent_ineighbor_alltoallv_)
      allocate(psb_comm_neighbor_handle :: handle, stat=info)
      if (info /= 0) return
      call handle%init(info)
      if (info /= 0) return
      handle%id = old_id
      handle%swap_status = old_swap_status
      select type(h => handle)
      type is(psb_comm_neighbor_handle)
        h%comm_type = comm_type
        h%use_persistent_buffers = (comm_type == psb_comm_persistent_ineighbor_alltoallv_)
      end select
    case default
      allocate(psb_comm_baseline_handle :: handle, stat=info)
      if (info /= 0) return
      call handle%init(info)
      if (info /= 0) return
      handle%id = old_id
      handle%swap_status = old_swap_status
    end select
  end subroutine psb_comm_init

  subroutine psb_comm_free(handle, info)
    implicit none
    class(psb_comm_handle_type), allocatable, intent(inout) :: handle
    integer(psb_ipk_), intent(out) :: info

    info = 0
    if (.not. allocated(handle)) return
    call handle%free(info)
    if (allocated(handle)) then
      deallocate(handle)
    end if
  end subroutine psb_comm_free


  ! Allocatable-based factory routines 
  subroutine psb_comm_create(comm_type, handle, info)
    implicit none
    integer(psb_ipk_), intent(in) :: comm_type
    class(psb_comm_handle_type), allocatable, intent(inout) :: handle
    integer(psb_ipk_), intent(out) :: info

    call psb_comm_init(comm_type, handle, info)
  end subroutine psb_comm_create

  subroutine psb_comm_destroy(handle, info)
    implicit none
    class(psb_comm_handle_type), allocatable, intent(inout) :: handle
    integer(psb_ipk_), intent(out) :: info

    call psb_comm_free(handle, info)
  end subroutine psb_comm_destroy

  subroutine psb_comm_set_swap_status(handle, flag, info)
    implicit none
    class(psb_comm_handle_type), allocatable, intent(inout) :: handle
    integer(psb_ipk_), intent(in) :: flag
    integer(psb_ipk_), intent(out) :: info
    info = 0
    if (.not. allocated(handle)) then
      info = -1
      return
    end if
    call handle%set_swap_status(flag, info)
  end subroutine psb_comm_set_swap_status

  subroutine psb_comm_get_swap_status(handle, flag, info)
    implicit none
    class(psb_comm_handle_type), allocatable, intent(in) :: handle
    integer(psb_ipk_), intent(out) :: flag
    integer(psb_ipk_), intent(out) :: info
    info = 0
    if (.not. allocated(handle)) then
      flag = 0
      info = -1
      return
    end if
    call handle%get_swap_status(flag, info)
  end subroutine psb_comm_get_swap_status

end module psb_comm_factory_mod
