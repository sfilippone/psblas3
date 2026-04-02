! Merged neighbor-topology module
!
module psb_comm_neighbor_impl_mod
  use psb_const_mod
  use psb_desc_const_mod, only: psb_proc_id_, psb_n_elem_recv_, psb_elem_recv_, &
       & psb_n_elem_send_, psb_elem_send_
  use psb_error_mod
  use psb_comm_schemes_mod, only: psb_comm_handle_type, psb_comm_ineighbor_alltoallv_
#ifdef PSB_MPI_MOD
  use mpi
#endif
  implicit none
#ifdef PSB_MPI_H
  include 'mpif.h'
#endif

  type, extends(psb_comm_handle_type) :: psb_comm_neighbor_handle
    integer(psb_mpk_)  :: graph_comm  = mpi_comm_null
    integer(psb_ipk_)  :: num_neighbors = 0
    integer(psb_mpk_), allocatable :: send_counts(:),  recv_counts(:)
    integer(psb_mpk_), allocatable :: send_displs(:),  recv_displs(:)
    integer(psb_ipk_), allocatable :: send_indexes(:)
    integer(psb_ipk_), allocatable :: recv_indexes(:)
    integer(psb_ipk_)  :: total_send = 0
    integer(psb_ipk_)  :: total_recv = 0
    logical :: is_initialized = .false.
    logical :: use_persistent_buffers = .false.
    integer(psb_mpk_) :: comm_request = mpi_request_null
    integer(psb_mpk_) :: persistent_request = mpi_request_null
    logical :: persistent_request_ready = .false.
    integer(psb_ipk_) :: persistent_buffer_size = 0
  contains
    procedure, pass :: init => psb_comm_neighbor_init
    procedure, pass :: free => neighbor_topology_free
    procedure, pass :: set_swap_status => psb_comm_neighbor_set_swap_status
    procedure, pass :: get_swap_status => psb_comm_neighbor_get_swap_status
    procedure, pass :: topology_init => neighbor_topology_init
    procedure, pass :: sizeof  => neighbor_topology_sizeof
  end type psb_comm_neighbor_handle


contains

  ! ---------------------------------------------------------------
  !  neighbor_topology_init
  !
  !  Parse the halo index list (obtained via desc_a%get_list_p)
  !  and build:
  !    - MPI dist-graph communicator with only the true neighbors
  !    - per-neighbor send/recv counts and displacements
  !    - contiguous gather/scatter index arrays
  !
  !  The topology is stored inside the vector and lazily built
  !  on the first psi_swapdata call that uses the neighbor-alltoallv
  !  communication mode.
  !
  !  Arguments:
  !    topology         - the persistent state (output, intent inout)
  !    halo_index       - halo_index array  (from get_list_p, intent in)
  !    num_neighbors    - number of exchanges (from get_list_p)
  !    total_send_elems - total send count   (from get_list_p)
  !    total_recv_elems - total recv count   (from get_list_p)
  !    ctxt             - PSBLAS context
  !    icomm            - MPI communicator
  !    info             - error code (output)
  ! ---------------------------------------------------------------
  subroutine neighbor_topology_init(topology, halo_index, num_neighbors, &
       & total_send_elems, total_recv_elems, ctxt, icomm, info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif
    
    class(psb_comm_neighbor_handle), intent(inout)    :: topology
    integer(psb_ipk_), intent(in)                     :: halo_index(:)
    integer(psb_ipk_), intent(in)                     :: num_neighbors, total_send_elems, total_recv_elems
    type(psb_ctxt_type), intent(in)                   :: ctxt
    integer(psb_mpk_), intent(in)                     :: icomm
    integer(psb_ipk_), intent(out)                    :: info

    ! locals
    integer(psb_mpk_)                                 :: iret
    integer(psb_ipk_)                                 :: i, k, idx_ptr, num_elem_recv, num_elem_send, partner_proc
    integer(psb_ipk_)                                 :: neighbor_count, send_offset, recv_offset
    integer(psb_mpk_), allocatable                    :: source_ranks(:), dest_ranks(:)
    integer(psb_mpk_), allocatable                    :: source_weights(:), dest_weights(:)
    integer(psb_mpk_)                                 :: in_degree, out_degree
    character(len=40)                                 :: name
    integer(psb_ipk_)                                 :: proc_id
    integer(psb_ipk_)                                 :: position 
    integer(psb_ipk_)                                 :: err_act

    info = psb_success_
    name = 'neighbor_topology_init'
    call psb_erractionsave(err_act)

    ! Clean up any previous state
    call topology%free(info)

    ! ----------------------------------------------------------
    !  First pass: count neighbors (excluding self) and totals
    ! ----------------------------------------------------------
    topology%num_neighbors     = 0
    topology%total_send        = 0
    topology%total_recv        = 0

    if(size(halo_index) < 1) then
      call psb_errpush(psb_err_topology_invalid_args_,name)
      goto 9999
    end if

    allocate(source_ranks(num_neighbors), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Source ranks allocation failed')
      goto 9999
    end if

    allocate(dest_ranks(num_neighbors), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Destination ranks allocation failed')
      goto 9999
    end if

    allocate(source_weights(num_neighbors), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Source weights allocation failed')
      goto 9999
    end if

    allocate(dest_weights(num_neighbors), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Destination weights allocation failed')
      goto 9999
    end if

    allocate(topology%send_counts(num_neighbors), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Send counts allocation failed')
      goto 9999
    end if

    allocate(topology%recv_counts(num_neighbors), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Receive counts allocation failed')
      goto 9999
    end if

    allocate(topology%send_displs(num_neighbors), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Send displacements allocation failed')
      goto 9999
    end if

    allocate(topology%recv_displs(num_neighbors), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Receive displacements allocation failed')
      goto 9999
    end if


    ! -----------------------------------------------------------
    ! Allocate the gather/scatter index arrays
    ! -----------------------------------------------------------
    allocate(topology%send_indexes(total_send_elems), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Send indexes allocation failed')
      goto 9999
    end if

    allocate(topology%recv_indexes(total_recv_elems), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name, a_err='Recv indexes allocation failed')
      goto 9999
    end if

    ! -----------------------------------------------------------
    ! Fill neighbor ranks, weights, counts, displacements,
    ! and gather/scatter index arrays.
    !
    ! The halo_index layout per neighbor (starting at position):
    !   position + 0              : process id
    !   position + 1              : nerv (num recv elements)
    !   position + 2 .. +1+nerv   : recv element indexes
    !   position + 2+nerv         : nesd (num send elements)
    !   position + 3+nerv .. +2+nerv+nesd : send element indexes
    ! Total stride per neighbor: nerv + nesd + 3
    ! -----------------------------------------------------------
    send_offset = 0
    recv_offset = 0
    position    = 1

    do i = 1, num_neighbors
      proc_id       = halo_index(position)
      num_elem_recv = halo_index(position + 1)
      num_elem_send = halo_index(position + num_elem_recv + 2)

      ! Fill source/destination ranks and weights (weights are all 1 for now)
      source_ranks(i)   = int(proc_id, psb_mpk_)
      dest_ranks(i)     = int(proc_id, psb_mpk_)
      source_weights(i) = 1
      dest_weights(i)   = 1

      ! Counts and displacements (displs set BEFORE accumulating offset)
      topology%send_counts(i) = int(num_elem_send, psb_mpk_)
      topology%recv_counts(i) = int(num_elem_recv, psb_mpk_)
      topology%send_displs(i) = int(send_offset, psb_mpk_)
      topology%recv_displs(i) = int(recv_offset, psb_mpk_)

      ! Fill recv_indexes from halo_index(position+2 .. position+1+nerv)
      do k = 1, num_elem_recv
        topology%recv_indexes(recv_offset + k) = halo_index(position + psb_elem_recv_ + k - 1)
      end do

      ! Fill send_indexes from halo_index(position+3+nerv .. position+2+nerv+nesd)
      do k = 1, num_elem_send
        topology%send_indexes(send_offset + k) = halo_index(position + num_elem_recv + psb_elem_send_ + k - 1)
      end do

      send_offset = send_offset + num_elem_send
      recv_offset = recv_offset + num_elem_recv

      topology%num_neighbors  = topology%num_neighbors + 1
      topology%total_send     = topology%total_send + num_elem_send
      topology%total_recv     = topology%total_recv + num_elem_recv

      position = position + num_elem_recv + num_elem_send + 3
    end do

    ! ----------------------------------------------------------
    !  Sanity check: the totals computed from the neighbor list
    !  should match the totals returned by get_list_p.  
    ! ----------------------------------------------------------
    if (topology%total_send /= total_send_elems) then
      info = psb_err_topology_args_mismatch_
      call psb_errpush(info, name, a_err='Send elements mismatch')
      goto 9999
    end if

    if (topology%total_recv /= total_recv_elems) then
      info = psb_err_topology_args_mismatch_
      call psb_errpush(info, name, a_err='Receive elements mismatch')
      goto 9999
    end if

    if(topology%num_neighbors /= num_neighbors) then
      info = psb_err_topology_args_mismatch_
      call psb_errpush(info, name, a_err='Number of neighbors mismatch')
      goto 9999
    end if


    ! ----------------------------------------------------------
    !  Build the dist-graph communicator
    ! ----------------------------------------------------------
    in_degree  = topology%num_neighbors !! Just for clarity
    out_degree = topology%num_neighbors !! Just for clarity

    call mpi_dist_graph_create_adjacent(icomm, &
         & in_degree,  source_ranks,  source_weights, &
         & out_degree, dest_ranks,    dest_weights,   &
         & mpi_info_null, .false.,                    & ! Check this line for optimizations
         & topology%graph_comm, info)
    if (info /= mpi_success) then
      info = psb_err_topology_error_
      call psb_errpush(info, name)
      goto 9999
    end if

    topology%is_initialized = .true.

    ! TODO: Is it safe to deallocate these temporary arrays here, or do we need them for the gather/scatter indexes?
    ! deallocate(source_ranks, dest_ranks, source_weights, dest_weights)
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine neighbor_topology_init


  ! ---------------------------------------------------------------
  !  neighbor_topology_free
  !  Release all resources held by the persistent state.
  ! ---------------------------------------------------------------
  subroutine neighbor_topology_free(this, info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif
    class(psb_comm_neighbor_handle), intent(inout)   :: this
    integer(psb_ipk_), intent(out) :: info
    integer(psb_mpk_) :: iret

    info = psb_success_

    if (this%persistent_request_ready) then
      if (this%persistent_request /= mpi_request_null) then
        call mpi_request_free(this%persistent_request, iret)
      end if
      this%persistent_request = mpi_request_null
      this%persistent_request_ready = .false.
      this%persistent_buffer_size = 0
    end if

    if (this%graph_comm /= mpi_comm_null) then
      call mpi_comm_free(this%graph_comm, iret)
      this%graph_comm = mpi_comm_null
    end if

    if (allocated(this%send_counts))      deallocate(this%send_counts)
    if (allocated(this%recv_counts))      deallocate(this%recv_counts)
    if (allocated(this%send_displs))      deallocate(this%send_displs)
    if (allocated(this%recv_displs))      deallocate(this%recv_displs)
    if (allocated(this%send_indexes))     deallocate(this%send_indexes)
    if (allocated(this%recv_indexes))     deallocate(this%recv_indexes)

    this%num_neighbors     = 0
    this%total_send        = 0
    this%total_recv        = 0
    this%is_initialized    = .false.
    this%comm_request      = mpi_request_null

  end subroutine neighbor_topology_free


  ! ---------------------------------------------------------------
  !  neighbor_topology_sizeof
  !  Return approximate memory footprint in bytes.
  ! ---------------------------------------------------------------
  function neighbor_topology_sizeof(this) result(val)
    implicit none
    class(psb_comm_neighbor_handle), intent(in) :: this
    integer(psb_epk_)                             :: val

    val = 0
    val = val + psb_sizeof_ip * 6   ! scalar integers + logicals
    if (allocated(this%send_counts))      val = val + psb_sizeof_ip * size(this%send_counts)
    if (allocated(this%recv_counts))      val = val + psb_sizeof_ip * size(this%recv_counts)
    if (allocated(this%send_displs))      val = val + psb_sizeof_ip * size(this%send_displs)
    if (allocated(this%recv_displs))      val = val + psb_sizeof_ip * size(this%recv_displs)
    if (allocated(this%send_indexes))     val = val + psb_sizeof_ip * size(this%send_indexes)
    if (allocated(this%recv_indexes))     val = val + psb_sizeof_ip * size(this%recv_indexes)


  end function neighbor_topology_sizeof


  subroutine psb_comm_neighbor_create(this, comm_type, info)
    class(psb_comm_neighbor_handle), intent(inout) :: this
    integer(psb_ipk_), intent(in) :: comm_type
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    this%comm_type = comm_type
    this%id = 0
    this%swap_status = 0
    this%comm_request = mpi_request_null
    this%persistent_request = mpi_request_null
    this%persistent_request_ready = .false.
    this%persistent_buffer_size = 0

    call this%free(info)
  end subroutine psb_comm_neighbor_create


  subroutine psb_comm_neighbor_destroy(this, info)
    class(psb_comm_neighbor_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    this%comm_request = mpi_request_null
    this%persistent_request = mpi_request_null
    this%persistent_request_ready = .false.
    this%persistent_buffer_size = 0
    call this%free(info)
  end subroutine psb_comm_neighbor_destroy


  subroutine psb_comm_neighbor_set_swap_status(this, flag, info)
    class(psb_comm_neighbor_handle), intent(inout) :: this
    integer(psb_ipk_), intent(in) :: flag
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    this%swap_status = flag
  end subroutine psb_comm_neighbor_set_swap_status


  subroutine psb_comm_neighbor_get_swap_status(this, flag, info)
    class(psb_comm_neighbor_handle), intent(in) :: this
    integer(psb_ipk_), intent(out) :: flag
    integer(psb_ipk_), intent(out) :: info

    info = psb_success_
    flag = this%swap_status
  end subroutine psb_comm_neighbor_get_swap_status

  subroutine psb_comm_neighbor_init(this, info)
    class(psb_comm_neighbor_handle), intent(inout) :: this
    integer(psb_ipk_), intent(out) :: info
    info = 0
    this%comm_type = psb_comm_ineighbor_alltoallv_
    this%id = 0
    this%swap_status = 0
    this%is_initialized = .false.
    this%use_persistent_buffers = .false.
    this%comm_request = mpi_request_null
    this%persistent_request = mpi_request_null
    this%persistent_request_ready = .false.
    this%persistent_buffer_size = 0
  end subroutine psb_comm_neighbor_init


end module psb_comm_neighbor_impl_mod
