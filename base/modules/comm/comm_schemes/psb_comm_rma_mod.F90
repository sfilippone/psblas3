module psb_comm_rma_mod
  use psb_const_mod
  use psb_desc_const_mod, only: psb_proc_id_, psb_n_elem_recv_, psb_elem_recv_, &
       & psb_n_elem_send_, psb_elem_send_
  use psb_error_mod
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_f_pointer
#ifdef PSB_MPI_MOD
  use mpi
#endif
    use psb_comm_schemes_mod, only: psb_comm_handle_type, psb_comm_rma_pull_, psb_comm_rma_push_, &
      & psb_comm_unknown_
  implicit none
  integer(psb_mpk_), parameter :: psb_rma_meta_tag = 913

  !
  ! Window pool.
  !
  ! MPI_Win_create and MPI_Win_create_dynamic are both collective over the whole
  ! communicator, and profiling at 448 ranks put either at ~0.9 s per call: with
  ! the window owned by the communication handle, every handle that came and
  ! went paid one. Measured, that was ~1400-1700 s aggregated, against ~3.5 s for
  ! the MPI_Get that actually moves the data.
  !
  ! So the window is owned by the module and keyed by communicator: created once
  ! on first use and kept for the run. Handles attach and detach their own
  ! buffers to it, which are local operations. This is the reuse policy PETSc
  ! applies to its window flavors; the flavor alone buys nothing without it.
  !
  integer(psb_ipk_), parameter :: psb_rma_max_wins = 16
  integer(psb_mpk_), private, save :: rma_pool_comm(psb_rma_max_wins) = mpi_comm_null
  integer(psb_mpk_), private, save :: rma_pool_win(psb_rma_max_wins)  = mpi_win_null
  integer(psb_ipk_), private, save :: rma_pool_n = 0
  public :: psb_comm_rma_get_window, psb_comm_rma_free_windows
#ifdef PSB_MPI_H
  include 'mpif.h'
#endif

  type, extends(psb_comm_handle_type) :: psb_comm_rma_handle
    integer(psb_mpk_) :: win = mpi_win_null
    logical :: window_ready = .false.
    logical :: window_open = .false.
    type(c_ptr) :: win_base = c_null_ptr
    integer(psb_ipk_) :: win_nelem = 0
    integer(psb_mpk_) :: nbr_grp = mpi_group_null   ! neighbor MPI group, built once with the layout
    logical :: layout_ready = .false.
    integer(psb_ipk_) :: layout_nnbr = -1
    integer(psb_ipk_) :: layout_send = -1
    integer(psb_ipk_) :: layout_recv = -1
    integer(psb_ipk_), allocatable :: peer_proc(:)
    integer(psb_ipk_), allocatable :: peer_send_counts(:)
    integer(psb_ipk_), allocatable :: peer_recv_counts(:)
    integer(psb_ipk_), allocatable :: peer_send_displs(:)
    integer(psb_ipk_), allocatable :: peer_recv_displs(:)
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    integer(psb_ipk_), allocatable :: peer_remote_send_displs(:)
    integer(psb_ipk_), allocatable :: peer_remote_recv_displs(:)
    integer(psb_ipk_), allocatable :: peer_send_indexes(:)
    integer(psb_ipk_), allocatable :: peer_recv_indexes(:)
    integer(psb_mpk_), allocatable :: notify_buf(:)
    integer(psb_mpk_), allocatable :: notify_recv_reqs(:)
    integer(psb_mpk_), allocatable :: notify_send_reqs(:)
    !
    ! Dynamic-window support.
    !
    ! MPI_Win_create is collective over the whole communicator and has to be
    ! repeated whenever the exposed buffer changes. Profiling on 448 ranks put
    ! it at ~0.9 s per call, against milliseconds for the MPI_Get that actually
    ! moves the data. With a dynamic window that collective is paid once and
    ! buffers are attached and detached locally.
    !
    ! The price is that on a dynamic window the target displacement is an
    ! absolute address in the target's address space rather than an offset in
    ! window units, so every rank must publish the address of its own buffer to
    ! its neighbours each time it (re)attaches.
    !
    logical :: buf_attached = .false.
    integer(kind=MPI_ADDRESS_KIND) :: my_buf_addr = 0
    integer(kind=MPI_ADDRESS_KIND), allocatable :: peer_buf_addr(:)
  contains
    procedure, pass :: init => psb_comm_rma_init
    procedure, pass :: free => psb_comm_rma_free
    procedure, pass :: clear_memory_buffer_layout => psb_comm_rma_clear_memory_buffer_layout
    procedure, pass :: init_memory_buffer_layout => psb_comm_rma_ini_memory_buffer_layout
    procedure, pass :: init_memory_buffer_layout_tran => psb_comm_rma_ini_memory_buffer_layout_tran
    procedure, pass :: set_swap_status => psb_comm_rma_set_swap_status
    procedure, pass :: get_swap_status => psb_comm_rma_get_swap_status
    procedure, pass :: publish_buf_addr => psb_comm_rma_publish_buf_addr
  end type psb_comm_rma_handle

contains

  !
  ! Hand back the dynamic window for this communicator, creating it if this is
  ! the first request. The collective is paid here, once per communicator per
  ! run, instead of once per handle.
  !
  subroutine psb_comm_rma_get_window(icomm, win, info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: icomm
    integer(psb_mpk_), intent(out) :: win
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: k
    integer(psb_mpk_) :: iret

    info = psb_success_
    win  = mpi_win_null

    do k = 1, rma_pool_n
      if (rma_pool_comm(k) == icomm) then
        win = rma_pool_win(k)
        return
      end if
    end do

    if (rma_pool_n >= psb_rma_max_wins) then
      ! More distinct communicators than the pool can hold. Raising the bound is
      ! the fix; silently creating an unpooled window would reintroduce the very
      ! per-handle collective this exists to remove.
      info = psb_err_internal_error_
      return
    end if

    call mpi_win_create_dynamic(mpi_info_null, icomm, win, iret)
    if (iret /= mpi_success) then
      info = psb_err_mpi_error_
      win  = mpi_win_null
      return
    end if

    rma_pool_n = rma_pool_n + 1
    rma_pool_comm(rma_pool_n) = icomm
    rma_pool_win(rma_pool_n)  = win
  end subroutine psb_comm_rma_get_window

  !
  ! Release every pooled window. Collective on each communicator involved, so it
  ! belongs at teardown, before MPI_Finalize.
  !
  subroutine psb_comm_rma_free_windows(info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: k
    integer(psb_mpk_) :: iret

    info = psb_success_
    do k = 1, rma_pool_n
      if (rma_pool_win(k) /= mpi_win_null) then
        call mpi_win_free(rma_pool_win(k), iret)
        if (iret /= mpi_success) info = psb_err_mpi_error_
      end if
      rma_pool_win(k)  = mpi_win_null
      rma_pool_comm(k) = mpi_comm_null
    end do
    rma_pool_n = 0
  end subroutine psb_comm_rma_free_windows

  subroutine psb_comm_rma_init(this, info)
    class(psb_comm_rma_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info
    info = 0
    this%comm_type = psb_comm_unknown_
    this%id = 0
    this%swap_status = 0
    this%win = mpi_win_null
    this%window_ready = .false.
    this%window_open = .false.
    this%win_base = c_null_ptr
    this%win_nelem = 0
    this%buf_attached = .false.
    this%my_buf_addr = 0
    if (allocated(this%peer_buf_addr)) deallocate(this%peer_buf_addr)
    this%layout_ready = .false.
    this%layout_nnbr = -1
    this%layout_send = -1
    this%layout_recv = -1
    call this%clear_memory_buffer_layout(info)
  end subroutine psb_comm_rma_init

  subroutine psb_comm_rma_clear_memory_buffer_layout(this, info)
    class(psb_comm_rma_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info
    integer(psb_mpk_) :: iret

    info = psb_success_
    ! Neighbors change when the layout is rebuilt: free the cached group.
    if (this%nbr_grp /= mpi_group_null) then
      call mpi_group_free(this%nbr_grp, iret)
      this%nbr_grp = mpi_group_null
    end if
    if (allocated(this%peer_proc)) deallocate(this%peer_proc)
    if (allocated(this%peer_send_counts)) deallocate(this%peer_send_counts)
    if (allocated(this%peer_recv_counts)) deallocate(this%peer_recv_counts)
    if (allocated(this%peer_send_displs)) deallocate(this%peer_send_displs)
    if (allocated(this%peer_recv_displs)) deallocate(this%peer_recv_displs)
    if (allocated(this%peer_mpi_rank)) deallocate(this%peer_mpi_rank)
    if (allocated(this%peer_remote_send_displs)) deallocate(this%peer_remote_send_displs)
    if (allocated(this%peer_remote_recv_displs)) deallocate(this%peer_remote_recv_displs)
    if (allocated(this%peer_send_indexes)) deallocate(this%peer_send_indexes)
    if (allocated(this%peer_recv_indexes)) deallocate(this%peer_recv_indexes)
    if (allocated(this%notify_buf)) deallocate(this%notify_buf)
    if (allocated(this%notify_recv_reqs)) deallocate(this%notify_recv_reqs)
    if (allocated(this%notify_send_reqs)) deallocate(this%notify_send_reqs)

    this%layout_ready = .false.
    this%layout_nnbr = -1
    this%layout_send = -1
    this%layout_recv = -1
  end subroutine psb_comm_rma_clear_memory_buffer_layout

    subroutine psb_comm_rma_ini_memory_buffer_layout(this, info, comm_list, peer_mpi_rank, &
      & num_neighbors, total_send, total_recv, my_rank, icomm)
    class(psb_comm_rma_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in) :: comm_list(:)
    integer(psb_ipk_), intent(in) :: peer_mpi_rank(:)
    integer(psb_ipk_), intent(in) :: num_neighbors, total_send, total_recv
    integer(psb_ipk_), intent(in) :: my_rank
    integer(psb_mpk_), intent(in) :: icomm

    integer(psb_mpk_) :: iret, p2pstat(mpi_status_size), prc_rank
    integer(psb_ipk_) :: n_neighbors, send_total, recv_total
    integer(psb_ipk_) :: neighbor_idx, item_idx
    integer(psb_ipk_) :: list_pos, send_offset, recv_offset
    integer(psb_ipk_) :: proc_to_comm, recv_count, send_count
    integer(psb_ipk_) :: local_meta(4), remote_meta(4)
    integer(psb_mpk_) :: comm_grp, n_off, grp_idx
    integer(psb_mpk_), allocatable :: grp_ranks(:)

    call this%clear_memory_buffer_layout(info)
    if (info /= psb_success_) return

    n_neighbors = num_neighbors
    send_total = total_send
    recv_total = total_recv

    if (n_neighbors > 0) then
       allocate(this%peer_proc(n_neighbors), this%peer_send_counts(n_neighbors), &
         & this%peer_recv_counts(n_neighbors), this%peer_send_displs(n_neighbors), &
         & this%peer_recv_displs(n_neighbors), this%peer_mpi_rank(n_neighbors), &
         & this%peer_remote_send_displs(n_neighbors), this%peer_remote_recv_displs(n_neighbors), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
      allocate(this%notify_buf(n_neighbors), this%notify_recv_reqs(n_neighbors), &
        & this%notify_send_reqs(n_neighbors), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
      this%notify_buf = 0_psb_mpk_
      this%notify_recv_reqs = MPI_REQUEST_NULL
      this%notify_send_reqs = MPI_REQUEST_NULL
    end if

    if (send_total > 0) then
      allocate(this%peer_send_indexes(send_total), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
    end if

    if (recv_total > 0) then
      allocate(this%peer_recv_indexes(recv_total), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
    end if

    list_pos = 1
    send_offset = 0
    recv_offset = 0
    do neighbor_idx = 1, n_neighbors
      proc_to_comm = comm_list(list_pos + psb_proc_id_)
      recv_count = comm_list(list_pos + psb_n_elem_recv_)
      send_count = comm_list(list_pos + recv_count + psb_n_elem_send_)

      this%peer_proc(neighbor_idx) = proc_to_comm
      this%peer_recv_counts(neighbor_idx) = recv_count
      this%peer_send_counts(neighbor_idx) = send_count
      this%peer_recv_displs(neighbor_idx) = recv_offset
      this%peer_send_displs(neighbor_idx) = send_offset
      this%peer_mpi_rank(neighbor_idx) = peer_mpi_rank(neighbor_idx)

      if (recv_count > 0) then
        do item_idx = 1, recv_count
          this%peer_recv_indexes(recv_offset + item_idx) = comm_list(list_pos + psb_elem_recv_ + item_idx - 1)
        end do
      end if

      if (send_count > 0) then
        do item_idx = 1, send_count
          this%peer_send_indexes(send_offset + item_idx) = comm_list(list_pos + recv_count + psb_elem_send_ + item_idx - 1)
        end do
      end if

      recv_offset = recv_offset + recv_count
      send_offset = send_offset + send_count
      list_pos = list_pos + recv_count + send_count + 3
    end do

    do neighbor_idx = 1, n_neighbors
      proc_to_comm = this%peer_proc(neighbor_idx)
      if (proc_to_comm /= my_rank) then
        prc_rank = this%peer_mpi_rank(neighbor_idx)
        local_meta = (/ this%peer_send_displs(neighbor_idx)+1, this%peer_send_counts(neighbor_idx), &
             & this%peer_recv_displs(neighbor_idx)+1+send_total, this%peer_recv_counts(neighbor_idx) /)
           call mpi_sendrecv(local_meta, 4, psb_mpi_mpk_, prc_rank, psb_rma_meta_tag, &
             & remote_meta, 4, psb_mpi_mpk_, prc_rank, psb_rma_meta_tag, icomm, p2pstat, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          return
        end if
        this%peer_remote_send_displs(neighbor_idx) = remote_meta(1)
        this%peer_remote_recv_displs(neighbor_idx) = remote_meta(3)
      else
        this%peer_remote_send_displs(neighbor_idx) = this%peer_send_displs(neighbor_idx)+1
        this%peer_remote_recv_displs(neighbor_idx) = this%peer_recv_displs(neighbor_idx)+1+send_total
      end if
    end do

    if ((send_offset /= send_total) .or. (recv_offset /= recv_total)) then
      info = psb_err_internal_error_
      return
    end if

    ! Build the neighbor MPI group ONCE (off-diagonal peers); reused by every
    ! RMA exchange (PSCW post/start) instead of rebuilding it each call.
    n_off = 0
    do neighbor_idx = 1, n_neighbors
      if (this%peer_proc(neighbor_idx) /= my_rank) n_off = n_off + 1
    end do
    if (n_off > 0) then
      allocate(grp_ranks(n_off), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
      grp_idx = 0
      do neighbor_idx = 1, n_neighbors
        if (this%peer_proc(neighbor_idx) /= my_rank) then
          grp_idx = grp_idx + 1
          grp_ranks(grp_idx) = this%peer_mpi_rank(neighbor_idx)
        end if
      end do
      call mpi_comm_group(icomm, comm_grp, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        return
      end if
      call mpi_group_incl(comm_grp, n_off, grp_ranks, this%nbr_grp, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        return
      end if
      call mpi_group_free(comm_grp, iret)
      deallocate(grp_ranks, stat=iret)
    end if

    this%layout_nnbr = n_neighbors
    this%layout_send = send_total
    this%layout_recv = recv_total
    this%layout_ready = .true.
  end subroutine psb_comm_rma_ini_memory_buffer_layout

  subroutine psb_comm_rma_free(this, info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    class(psb_comm_rma_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info
    integer(psb_mpk_) :: iret
    real(psb_dpk_), pointer :: detach_p(:)

    info = 0
    if (this%window_open) then
      call mpi_win_unlock_all(this%win, iret)
      this%window_open = .false.
    end if
    !
    ! The window belongs to the module pool, not to this handle: it is not freed
    ! here. What must go is the attachment, because the buffer behind it is about
    ! to be deallocated with the vector, and leaving a dangling region attached
    ! would also make the next attach of the same address overlap.
    !
    ! MPI_Win_detach wants the buffer, not its address, and this module is
    ! type-generic; c_f_pointer on the stored c_ptr gives an array at the right
    ! address, which is all detach looks at.
    !
    if (this%buf_attached .and. c_associated(this%win_base) &
         & .and. (this%win /= mpi_win_null)) then
      call c_f_pointer(this%win_base, detach_p, [1])
      call mpi_win_detach(this%win, detach_p, iret)
      if (iret /= mpi_success) info = psb_err_mpi_error_
    end if
    this%buf_attached = .false.
    this%win = mpi_win_null
    this%win_base = c_null_ptr
    this%win_nelem = 0
    ! Freeing the window drops whatever is attached to it, so the attachment
    ! bookkeeping goes with it.
    this%buf_attached = .false.
    this%my_buf_addr = 0
    if (allocated(this%peer_buf_addr)) deallocate(this%peer_buf_addr)
    this%window_ready = .false.
    this%layout_ready = .false.
    this%layout_nnbr = -1
    this%layout_send = -1
    this%layout_recv = -1
    call this%clear_memory_buffer_layout(info)
  end subroutine psb_comm_rma_free

  !
  ! Publish the address of the locally attached buffer to the neighbours, and
  ! collect theirs. On a dynamic window MPI_Get/MPI_Put take an absolute address
  ! in the target's address space, so this has to be redone whenever the buffer
  ! is re-attached. Point-to-point over the neighbour list, reusing the metadata
  ! tag already used for the displacement exchange: no collective.
  !
  subroutine psb_comm_rma_publish_buf_addr(this, icomm, info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif
    class(psb_comm_rma_handle), intent(inout) :: this
    integer(psb_mpk_), intent(in)  :: icomm
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: k, nnbr
    integer(psb_mpk_) :: prc_rank, iret
    integer(psb_mpk_) :: p2pstat(mpi_status_size)

    info = psb_success_
    if (.not.allocated(this%peer_mpi_rank)) return
    nnbr = size(this%peer_mpi_rank)

    if (allocated(this%peer_buf_addr)) then
      if (size(this%peer_buf_addr) /= nnbr) deallocate(this%peer_buf_addr)
    end if
    if (.not.allocated(this%peer_buf_addr)) then
      allocate(this%peer_buf_addr(nnbr), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
    end if

    do k = 1, nnbr
      prc_rank = int(this%peer_mpi_rank(k), psb_mpk_)
      call mpi_sendrecv(this%my_buf_addr,     1, mpi_aint, prc_rank, psb_rma_meta_tag, &
           &            this%peer_buf_addr(k), 1, mpi_aint, prc_rank, psb_rma_meta_tag, &
           &            icomm, p2pstat, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        return
      end if
    end do
  end subroutine psb_comm_rma_publish_buf_addr

  ! Transpose variant: peer_send_* is filled from comm_list RECV area,
  ! peer_recv_* from comm_list SEND area. Metadata exchange tells peers our
  ! recv displacement so they can GET/PUT to the correct location in swaptran.
  ! layout_send = total_recv (effective send), layout_recv = total_send.
  subroutine psb_comm_rma_ini_memory_buffer_layout_tran(this, info, comm_list, peer_mpi_rank, &
      & num_neighbors, total_send, total_recv, my_rank, icomm)
    class(psb_comm_rma_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in) :: comm_list(:)
    integer(psb_ipk_), intent(in) :: peer_mpi_rank(:)
    integer(psb_ipk_), intent(in) :: num_neighbors, total_send, total_recv
    integer(psb_ipk_), intent(in) :: my_rank
    integer(psb_mpk_), intent(in) :: icomm

    integer(psb_mpk_) :: iret, p2pstat(mpi_status_size), prc_rank
    integer(psb_ipk_) :: n_neighbors, eff_send, eff_recv
    integer(psb_ipk_) :: neighbor_idx, item_idx
    integer(psb_ipk_) :: list_pos, send_offset, recv_offset
    integer(psb_ipk_) :: proc_to_comm, actual_recv_count, actual_send_count
    integer(psb_ipk_) :: local_meta(4), remote_meta(4)

    call this%clear_memory_buffer_layout(info)
    if (info /= psb_success_) return

    n_neighbors = num_neighbors
    ! For transpose: effective send = actual recv, effective recv = actual send
    eff_send = total_recv
    eff_recv = total_send

    if (n_neighbors > 0) then
      allocate(this%peer_proc(n_neighbors), this%peer_send_counts(n_neighbors), &
        & this%peer_recv_counts(n_neighbors), this%peer_send_displs(n_neighbors), &
        & this%peer_recv_displs(n_neighbors), this%peer_mpi_rank(n_neighbors), &
        & this%peer_remote_send_displs(n_neighbors), this%peer_remote_recv_displs(n_neighbors), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
      allocate(this%notify_buf(n_neighbors), this%notify_recv_reqs(n_neighbors), &
        & this%notify_send_reqs(n_neighbors), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
      this%notify_buf = 0_psb_mpk_
      this%notify_recv_reqs = MPI_REQUEST_NULL
      this%notify_send_reqs = MPI_REQUEST_NULL
    end if

    if (eff_send > 0) then
      allocate(this%peer_send_indexes(eff_send), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
    end if

    if (eff_recv > 0) then
      allocate(this%peer_recv_indexes(eff_recv), stat=iret)
      if (iret /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
    end if

    list_pos = 1
    send_offset = 0
    recv_offset = 0
    do neighbor_idx = 1, n_neighbors
      proc_to_comm = comm_list(list_pos + psb_proc_id_)
      actual_recv_count = comm_list(list_pos + psb_n_elem_recv_)
      actual_send_count = comm_list(list_pos + actual_recv_count + psb_n_elem_send_)

      this%peer_proc(neighbor_idx) = proc_to_comm
      ! Swap: effective send = actual recv, effective recv = actual send
      this%peer_send_counts(neighbor_idx) = actual_recv_count
      this%peer_recv_counts(neighbor_idx) = actual_send_count
      this%peer_send_displs(neighbor_idx) = send_offset
      this%peer_recv_displs(neighbor_idx) = recv_offset
      this%peer_mpi_rank(neighbor_idx) = peer_mpi_rank(neighbor_idx)

      ! peer_send_indexes from RECV area of comm_list (these are indices we gather for "sending")
      if (actual_recv_count > 0) then
        do item_idx = 1, actual_recv_count
          this%peer_send_indexes(send_offset + item_idx) = comm_list(list_pos + psb_elem_recv_ + item_idx - 1)
        end do
      end if

      ! peer_recv_indexes from SEND area of comm_list (these are indices we scatter into)
      if (actual_send_count > 0) then
        do item_idx = 1, actual_send_count
          this%peer_recv_indexes(recv_offset + item_idx) = comm_list(list_pos + actual_recv_count + psb_elem_send_ + item_idx - 1)
        end do
      end if

      send_offset = send_offset + actual_recv_count
      recv_offset = recv_offset + actual_send_count
      list_pos = list_pos + actual_recv_count + actual_send_count + 3
    end do

    do neighbor_idx = 1, n_neighbors
      proc_to_comm = this%peer_proc(neighbor_idx)
      if (proc_to_comm /= my_rank) then
        prc_rank = this%peer_mpi_rank(neighbor_idx)
        ! Tell peer: our effective-send (=actual recv) displacement and count,
        !            our effective-recv (=actual send) displacement (offset past eff_send area).
        local_meta = (/ this%peer_send_displs(neighbor_idx)+1, this%peer_send_counts(neighbor_idx), &
             & this%peer_recv_displs(neighbor_idx)+1+eff_send, this%peer_recv_counts(neighbor_idx) /)
        call mpi_sendrecv(local_meta, 4, psb_mpi_mpk_, prc_rank, psb_rma_meta_tag, &
          & remote_meta, 4, psb_mpi_mpk_, prc_rank, psb_rma_meta_tag, icomm, p2pstat, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          return
        end if
        this%peer_remote_send_displs(neighbor_idx) = remote_meta(1)
        this%peer_remote_recv_displs(neighbor_idx) = remote_meta(3)
      else
        this%peer_remote_send_displs(neighbor_idx) = this%peer_send_displs(neighbor_idx)+1
        this%peer_remote_recv_displs(neighbor_idx) = this%peer_recv_displs(neighbor_idx)+1+eff_send
      end if
    end do

    if ((send_offset /= eff_send) .or. (recv_offset /= eff_recv)) then
      info = psb_err_internal_error_
      return
    end if

    this%layout_nnbr = n_neighbors
    this%layout_send = eff_send
    this%layout_recv = eff_recv
    this%layout_ready = .true.
  end subroutine psb_comm_rma_ini_memory_buffer_layout_tran

  subroutine psb_comm_rma_set_swap_status(this, flag, info)
    class(psb_comm_rma_handle), intent(inout) :: this
    integer(psb_ipk_), intent(in) :: flag
    integer(psb_ipk_), intent(out) :: info
    info = 0
    this%swap_status = flag
  end subroutine psb_comm_rma_set_swap_status

  subroutine psb_comm_rma_get_swap_status(this, flag, info)
    class(psb_comm_rma_handle), intent(in) :: this
    integer(psb_ipk_), intent(out) :: flag
    integer(psb_ipk_), intent(out) :: info
    info = 0
    flag = this%swap_status
  end subroutine psb_comm_rma_get_swap_status

end module psb_comm_rma_mod