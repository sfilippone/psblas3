!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!
! Module: psb_neighbor_topology_mod
!   Provides a type to hold pre-built MPI neighborhood topology
!   information for persistent/repeated halo exchanges via
!   MPI_Neighbor_alltoallv  (MPI >= 3.0).
!
!   The topology is stored inside the vector type (psb_d_base_vect_type)
!   and lazily created on the first psi_swapdata call with the
!   neighbor-alltoallv communication mode.  Once built it is reused
!   for every subsequent halo exchange, avoiding the per-call overhead
!   of re-scanning the index list and allocating temporary arrays.
!
!   The graph communicator and per-neighbor counts/displacements
!   are built once and reused.
!
!   The gather/scatter index arrays (send_indexes, recv_indexes) record
!   which local vector positions must be packed / unpacked.
!
module psb_neighbor_topology_mod
  use psb_const_mod
  use psb_desc_const_mod
  use psb_error_mod
  !
  ! Only import mpi_comm_null here (needed for type default initializer).
  ! Full MPI access is done inside each contained subroutine so that
  ! MPI symbols do NOT leak into modules that use psb_neighbor_topology_mod.
  !
#ifdef PSB_MPI_MOD
  use mpi, only: mpi_comm_null
#endif
  implicit none
#ifdef PSB_MPI_H
  include 'mpif.h'
#endif

  type :: psb_neighbor_topology_type
    !
    ! MPI dist-graph communicator (only communicating neighbors).
    !
    integer(psb_mpk_)  :: graph_comm  = mpi_comm_null
    !
    ! Number of neighbors (processes I exchange with, excluding self).
    !
    integer(psb_ipk_)  :: num_neighbors = 0
    !
    ! Per-neighbor send/recv counts and displacements (units of
    ! single elements; for n-column multivectors multiply by n).
    !   send_counts(i)   = number of elements sent to i-th neighbor
    !   recv_counts(i)   = number of elements received from i-th neighbor
    !   send_displs(i)   = displacement into contiguous send buffer
    !   recv_displs(i)   = displacement into contiguous recv buffer
    !
    integer(psb_mpk_), allocatable :: send_counts(:),  recv_counts(:)
    integer(psb_mpk_), allocatable :: send_displs(:),  recv_displs(:)
    !
    ! Gather indexes: the k-th element of the send buffer is
    !   y%v( send_indexes(k) )
    ! Scatter indexes: the k-th element of the recv buffer goes to
    !   y%v( recv_indexes(k) )
    !
    integer(psb_ipk_), allocatable :: send_indexes(:)
    integer(psb_ipk_), allocatable :: recv_indexes(:)
    !
    ! Total number of elements to send / receive (per single column),
    ! excluding self-exchange.
    !
    integer(psb_ipk_)  :: total_send = 0
    integer(psb_ipk_)  :: total_recv = 0
    !

    ! Initialization flag.
    !
    logical             :: is_initialized = .false.
  contains
    procedure, pass(topology) :: init    => neighbor_topology_init
    procedure, pass(topology) :: free    => neighbor_topology_free
    procedure, pass(topology) :: sizeof  => neighbor_topology_sizeof
  end type psb_neighbor_topology_type

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
    
    class(psb_neighbor_topology_type), intent(inout)  :: topology
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
  subroutine neighbor_topology_free(topology, info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif
    class(psb_neighbor_topology_type), intent(inout) :: topology
    integer(psb_ipk_), intent(out) :: info
    integer(psb_mpk_) :: iret

    info = psb_success_

    if (topology%graph_comm /= mpi_comm_null) then
      call mpi_comm_free(topology%graph_comm, iret)
      topology%graph_comm = mpi_comm_null
    end if

    if (allocated(topology%send_counts))      deallocate(topology%send_counts)
    if (allocated(topology%recv_counts))      deallocate(topology%recv_counts)
    if (allocated(topology%send_displs))      deallocate(topology%send_displs)
    if (allocated(topology%recv_displs))      deallocate(topology%recv_displs)
    if (allocated(topology%send_indexes))     deallocate(topology%send_indexes)
    if (allocated(topology%recv_indexes))     deallocate(topology%recv_indexes)

    topology%num_neighbors     = 0
    topology%total_send        = 0
    topology%total_recv        = 0
    topology%is_initialized    = .false.

  end subroutine neighbor_topology_free


  ! ---------------------------------------------------------------
  !  neighbor_topology_sizeof
  !  Return approximate memory footprint in bytes.
  ! ---------------------------------------------------------------
  function neighbor_topology_sizeof(topology) result(val)
    implicit none
    class(psb_neighbor_topology_type), intent(in) :: topology
    integer(psb_epk_)                             :: val

    val = 0
    val = val + psb_sizeof_ip * 6   ! scalar integers + logicals
    if (allocated(topology%send_counts))      val = val + psb_sizeof_ip * size(topology%send_counts)
    if (allocated(topology%recv_counts))      val = val + psb_sizeof_ip * size(topology%recv_counts)
    if (allocated(topology%send_displs))      val = val + psb_sizeof_ip * size(topology%send_displs)
    if (allocated(topology%recv_displs))      val = val + psb_sizeof_ip * size(topology%recv_displs)
    if (allocated(topology%send_indexes))     val = val + psb_sizeof_ip * size(topology%send_indexes)
    if (allocated(topology%recv_indexes))     val = val + psb_sizeof_ip * size(topology%recv_indexes)


  end function neighbor_topology_sizeof

end module psb_neighbor_topology_mod