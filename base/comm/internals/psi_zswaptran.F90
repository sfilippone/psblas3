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
!
! File: psi_zswaptran.F90
!
! Subroutine: psi_zswaptran_vect
!   Implements the data exchange among processes. This is similar to Xswapdata, but
!   the list is read "in reverse", i.e. indices that are normally SENT are used 
!   for the RECEIVE part and vice-versa. This is the basic data exchange operation
!   for doing the product of a sparse matrix by a vector. 
!   Essentially this is doing a variable all-to-all data exchange
!   (ALLTOALLV in MPI parlance), but 
!   it is capable of pruning empty exchanges, which are very likely in out 
!   application environment. All the variants have the same structure 
!   In all these subroutines X may be:    I    Integer
!                                         S    real(psb_spk_)
!                                         D    complex(psb_dpk_)
!                                         C    complex(psb_spk_)
!                                         Z    complex(psb_dpk_)
!   Basically the operation is as follows: on each process, we identify 
!   sections SND(Y) and RCV(Y); then we do a SEND(PACK(GTH(SND(Y))));
!   then we receive, and we do an update with Y = SCT(RCV(Y)) + BETA * Y 
!   but only on the elements involved in the SCT operation. 
!   Thus: for halo data exchange, the receive section is confined in the 
!   halo indices, and BETA=0, whereas for overlap exchange the receive section 
!   is scattered in the owned indices, and BETA=1.
!   The first routine picks the desired exchange index list and passes it to the second.
!   This version works on encapsulated vectors, and uses their methods to do  GTH and SCT,
!   so that special versions (i.e. GPU vectors can override them 
! 
! 
! Arguments: 
!    swap_status     - integer                 Choose the algorithm for data exchange: 
!                                       this is chosen through bit fields. 
!                                        swap_mpi  = iand(swap_status,psb_swap_mpi_)  /= 0
!                                        swap_sync = iand(swap_status,psb_swap_sync_) /= 0
!                                        swap_send = iand(swap_status,psb_swap_send_) /= 0
!                                        swap_recv = iand(swap_status,psb_swap_recv_) /= 0
!                                       if (swap_mpi):  use underlying MPI_ALLTOALLV.
!                                       if (swap_sync): use PSB_SND and PSB_RCV in 
!                                                       synchronized pairs
!                                       if (swap_send .and. swap_recv): use mpi_irecv 
!                                                       and mpi_send
!                                       if (swap_send): use psb_snd (but need another 
!                                                       call with swap_recv to complete)
!                                       if (swap_recv): use psb_rcv (completing a 
!                                                       previous call with swap_send)
!
!
!    n        - integer                 Number of columns in Y               
!    beta     - real                  Choose overwrite or sum. 
!    y        - type(psb_d_vect_type) The data area                        
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!                                       our own internal allocation.
!    info     - integer.                return code.
!    data     - integer                 which list is to be used to exchange data
!                                       default psb_comm_halo_
!                                       psb_comm_halo_    use halo_index
!                                       psb_comm_ext_     use ext_index 
!                                       psb_comm_ovrl_    use ovrl_index
!                                       psb_comm_mov_     use ovr_mst_idx
!
!   
submodule (psi_z_comm_v_mod)  psi_z_swaptran_impl
  use psb_base_mod
  use psb_comm_factory_mod
contains
  module subroutine psi_zswaptran_vect(swap_status,beta,y,desc_a,info,data)

#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    integer(psb_ipk_), intent(in)                 :: swap_status
    complex(psb_dpk_), intent(in)                    :: beta
    class(psb_z_base_vect_type), intent(inout)    :: y
    type(psb_desc_type),target                    :: desc_a
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), optional                   :: data

    ! locals
    type(psb_ctxt_type)                           :: ctxt
    integer(psb_mpk_)                             :: icomm
    integer(psb_ipk_)                             :: np, me, total_send, total_recv, num_neighbors, err_act, data_
    class(psb_i_base_vect_type), pointer          :: comm_indexes
    character(len=20)  :: name

    info = psb_success_
    name = 'psi_zswaptran_vect'
    call psb_erractionsave(err_act)

    ctxt = desc_a%get_context()
    icomm = ctxt%get_mpic()
    call psb_info(ctxt,me,np) 
    if (np == -1) then
      info=psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif

    if (.not.psb_is_asb_desc(desc_a)) then 
      info=psb_err_invalid_cd_state_
      call psb_errpush(info,name)
      goto 9999
    endif

    if (present(data)) then
      data_ = data
    else
      data_ = psb_comm_halo_
    end if

    call desc_a%get_list_p(data_,comm_indexes,num_neighbors,total_recv,total_send,info) 
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,a_err='psb_cd_get_list')
      goto 9999
    end if

    if( (swap_status /= psb_comm_status_start_).and.(swap_status /= psb_comm_status_wait_)&
      & .and.(swap_status /= psb_comm_status_sync_) ) then
        info = psb_err_mpi_error_
        call psb_errpush(info,name,a_err='Invalid swap_status swap_status')
        goto 9999
    end if

    if (.not. allocated(y%comm_handle)) then
      call psb_comm_set(psb_comm_isend_irecv_, y%comm_handle, info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_, name, a_err='init comm default baseline')
        goto 9999
      end if
    end if

    ! Set the normalized swap status on the comm handle
    call y%comm_handle%set_swap_status(swap_status, info)
    if (info /= psb_success_) then
      call psb_errpush(info,name,a_err='set_swap_status')
      goto 9999
    end if

    select case(y%comm_handle%comm_type)
    case(psb_comm_ineighbor_alltoallv_)
      call psi_ztran_neighbor_topology_vect(ctxt,swap_status,beta,y,comm_indexes,&
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='neighbor a2av tran')
        goto 9999
      end if
    case(psb_comm_persistent_ineighbor_alltoallv_)
      call psi_ztran_neighbor_persistent_topology_vect(ctxt,swap_status,beta,y,comm_indexes,&
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='persistent neighbor tran')
        goto 9999
      end if
    case(psb_comm_rma_pull_)
      call psi_ztran_rma_pull_vect(ctxt,swap_status,beta,y,comm_indexes,&
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='rma pull tran')
        goto 9999
      end if
    case(psb_comm_rma_push_)
      call psi_ztran_rma_push_vect(ctxt,swap_status,beta,y,comm_indexes,&
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='rma push tran')
        goto 9999
      end if
    case default
      call psi_ztran_baseline_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='baseline tran')
        goto 9999
      end if
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psi_zswaptran_vect

  !
  !
  ! Subroutine: psi_ztran_baseline_vect
  !   Data exchange among processes.
  !
  !   Takes care of Y an encapsulated vector. Relies on the gather/scatter methods
  !   of vectors. 
  !   
  !   The real workhorse: the outer routine will only choose the index list
  !   this one takes the index list and does the actual exchange. 
  !   
  !   
  ! 
  module subroutine psi_ztran_baseline_vect(ctxt,swap_status,beta,y,idx,&
      & num_neighbors,total_send,total_recv,comm_handle,info)

#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)             :: ctxt
    integer(psb_ipk_), intent(in)               :: swap_status
    complex(psb_dpk_), intent(in)                  :: beta
    class(psb_z_base_vect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)  :: idx
    integer(psb_ipk_), intent(in)               :: num_neighbors,total_send, total_recv
    class(psb_comm_handle_type), intent(inout)  :: comm_handle
    integer(psb_ipk_), intent(out)              :: info

    ! locals
    integer(psb_mpk_)   :: np, me, nesd, nerv, n
    integer(psb_mpk_)   :: proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_), allocatable :: prcid(:)
    type(psb_comm_baseline_handle), pointer     :: baseline_comm_handle
    integer(psb_ipk_) :: err_act, i, idx_pt, total_send_, total_recv_,&
         & snd_pt, rcv_pt, pnti
    logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
         & albf,do_send,do_recv
    logical, parameter :: usersend=.false., debug=.false.
    character(len=20)  :: name

    info = psb_success_
    name = 'psi_ztran_baseline_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,me,np) 
    if (np == -1) then
      info=psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif
    icomm = ctxt%get_mpic()

    baseline_comm_handle => null()
    select type(ch => comm_handle)
    type is(psb_comm_baseline_handle)
      baseline_comm_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected baseline comm_handle in baseline swaptran')
      goto 9999
    end select

    n=1
    do_send = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_recv = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    total_recv_ = total_recv * n
    total_send_ = total_send * n

    call idx%sync()

    if (debug) write(*,*) me,'Internal buffer'
    if (do_send) then 
      if (allocated(baseline_comm_handle%comid)) then
        if (any(baseline_comm_handle%comid /= mpi_request_null)) then 
          ! 
          ! Unfinished communication? Something is wrong....
          !
          info=psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/-2/))
          goto 9999
        end if
      end if
      if (debug) write(*,*) me,'do_send start'
      call y%new_buffer(ione*size(idx%v),info)
      call psb_realloc(num_neighbors,2_psb_ipk_,baseline_comm_handle%comid,info)
      baseline_comm_handle%comid = mpi_request_null
      call psb_realloc(num_neighbors,prcid,info)
      ! First I post all the non blocking receives
      pnti   = 1
      p2ptag = psb_double_swap_tag
      do i=1, num_neighbors
        proc_to_comm = idx%v(pnti+psb_proc_id_)
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)

        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_
        prcid(i) = psb_get_mpi_rank(ctxt,proc_to_comm)      
        if ((nesd>0).and.(proc_to_comm /= me)) then 
          if (debug) write(*,*) me,'Posting receive from',prcid(i),rcv_pt
          call mpi_irecv(y%combuf(snd_pt),nesd,&
               & psb_mpi_c_dpk_,prcid(i),&
            & p2ptag, icomm,baseline_comm_handle%comid(i,2),iret)
        end if
        pnti   = pnti + nerv + nesd + 3
      end do

      if (debug) write(*,*) me,' Gather '
      !
      ! Then gather for sending.
      !    
      pnti   = 1
      snd_pt = 1
      do i=1, num_neighbors
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_

        idx_pt = rcv_pt
        call y%gth(idx_pt,nerv,idx)

        pnti   = pnti + nerv + nesd + 3
      end do

      !
      ! Then wait 
      !
      call y%device_wait()

      if (debug) write(*,*) me,' isend'
      !
      ! Then send
      !

      pnti   = 1
      snd_pt = 1
      rcv_pt = 1
      p2ptag = psb_double_swap_tag
      do i=1, num_neighbors
        proc_to_comm = idx%v(pnti+psb_proc_id_)
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_

        if ((nerv>0).and.(proc_to_comm /= me)) then 
          call mpi_isend(y%combuf(rcv_pt),nerv,&
               & psb_mpi_c_dpk_,prcid(i),&
            & p2ptag,icomm,baseline_comm_handle%comid(i,1),iret)
        end if

        if(iret /= mpi_success) then
          info=psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if

        pnti   = pnti + nerv + nesd + 3
      end do
    end if

    if (do_recv) then 
      if (debug) write(*,*) me,' do_Recv'
      if (.not.allocated(baseline_comm_handle%comid)) then 
        ! 
        ! No matching send? Something is wrong....
        !
        info=psb_err_mpi_error_
        call psb_errpush(info,name,m_err=(/-2/))
        goto 9999
      end if
      call psb_realloc(num_neighbors,prcid,info)

      if (debug) write(*,*) me,' wait'
      pnti   = 1
      p2ptag = psb_double_swap_tag
      do i=1, num_neighbors
        proc_to_comm = idx%v(pnti+psb_proc_id_)
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_

        if (proc_to_comm /= me)then 
          if (nerv>0) then 
            call mpi_wait(baseline_comm_handle%comid(i,1),p2pstat,iret)
            if(iret /= mpi_success) then
              info=psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          end if
          if (nesd>0) then 
            call mpi_wait(baseline_comm_handle%comid(i,2),p2pstat,iret)
            if(iret /= mpi_success) then
              info=psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          end if
        else if (proc_to_comm == me) then 
          if (nesd /= nerv) then 
            write(psb_err_unit,*) &
                 & 'Fatal error in swapdata: mismatch on self send',&
                 & nerv,nesd
          end if
          y%combuf(snd_pt:snd_pt+nesd-1) = y%combuf(rcv_pt:rcv_pt+nerv-1) 
        end if
        pnti   = pnti + nerv + nesd + 3
      end do

      if (debug) write(*,*) me,' scatter'      
      pnti   = 1
      snd_pt = 1
      rcv_pt = 1
      do i=1, num_neighbors
        proc_to_comm = idx%v(pnti+psb_proc_id_)
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        idx_pt = 1+pnti+psb_n_elem_recv_
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_

        if (debug) write(0,*)me,' Received from: ',prcid(i),&
             & y%combuf(snd_pt:snd_pt+nesd-1)        
        call y%sct(snd_pt,nesd,idx,beta)
        pnti   = pnti + nerv + nesd + 3
      end do
      !
      ! Waited for everybody, clean up
      !
      baseline_comm_handle%comid = mpi_request_null

      !
      ! Then wait for device
      !
      if (debug) write(*,*) me,' wait'
      call y%device_wait()
      if (debug) write(*,*) me,' free buffer'
      call y%maybe_free_buffer(info)
      if (info == 0) then
        if (allocated(y%comm_handle)) call psb_comm_free(y%comm_handle, info)
      end if
      if (info /= 0) then 
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
      if (debug) write(*,*) me,' done'
    end if


    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return

  end subroutine psi_ztran_baseline_vect




  subroutine psi_ztran_neighbor_topology_vect(ctxt,swap_status,beta,y,comm_indexes,&
    & num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
  use mpi
#endif
  implicit none
#ifdef PSB_MPI_H
  include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)             :: ctxt
    integer(psb_ipk_), intent(in)               :: swap_status
    complex(psb_dpk_), intent(in)                  :: beta
    class(psb_z_base_vect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)  :: comm_indexes
    integer(psb_ipk_), intent(in)               :: num_neighbors,total_send,total_recv
    class(psb_comm_handle_type), intent(inout)  :: comm_handle
    integer(psb_ipk_), intent(out)              :: info

    ! locals
    integer(psb_mpk_)                           :: icomm
    integer(psb_mpk_)                           :: np, me
    integer(psb_mpk_)                           :: iret, p2pstat(mpi_status_size)
    type(psb_comm_neighbor_handle), pointer     :: neighbor_comm_handle
    integer(psb_ipk_)                           :: err_act, topology_total_send, topology_total_recv, buffer_size
    logical                                     :: do_start, do_wait
    logical, parameter                          :: debug = .false.
    character(len=30)                           :: name


    info = psb_success_
    name = 'psi_ztran_neighbor_topology_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,me,np) 
    if (np == -1) then
      info=psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif

    icomm = ctxt%get_mpic()

    neighbor_comm_handle => null()
    select type(ch => comm_handle)
    type is(psb_comm_neighbor_handle)
      neighbor_comm_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected neighbor comm_handle in neighbor swaptran')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    ! ---------------------------------------------------------
    !  START phase: build topology (if needed), gather, post MPI
    ! ---------------------------------------------------------
    if (do_start) then
      if(debug) write(*,*) me,' nbr_vect: starting data exchange'
      ! Lazy initialization: build the topology on first call
      if (.not. neighbor_comm_handle%is_initialized) then
        if (debug) write(*,*) me,' nbr_vect: building topology'
        call neighbor_comm_handle%topology_init(comm_indexes%v, num_neighbors, total_send, total_recv, ctxt, icomm, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_, name, &
              & a_err='neighbor_topology_init')
          goto 9999
        end if
      end if

      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      ! Buffer layout:
      !   combuf(1 : total_send)                       = send area
      !   combuf(total_send+1 : total_send+total_recv) = recv area
      buffer_size = topology_total_send + topology_total_recv

      call y%new_buffer(buffer_size, info)
      if (info /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_, name)
        goto 9999
      end if
      neighbor_comm_handle%comm_request = mpi_request_null

        ! For transpose exchange: gather recv area first (we will send "recv" data)
        if (debug) write(*,*) me,' nbr_tran_vect: gathering recv data,', topology_total_recv,' elems'
        call y%gth(int(topology_total_recv,psb_mpk_), &
          & neighbor_comm_handle%recv_indexes, &
          & y%combuf(1:topology_total_recv))

        ! Wait for device (important for GPU subclasses)
        call y%device_wait()

        ! Post non-blocking neighborhood alltoallv swapping send/recv arrays
        if (debug) write(*,*) me,' nbr_tran_vect: posting MPI_Ineighbor_alltoallv (swapped)'
        call mpi_ineighbor_alltoallv( &
          & y%combuf(1),                           &  ! send buffer (recv_indexes gathered)
          & neighbor_comm_handle%recv_counts,       &
          & neighbor_comm_handle%recv_displs,       &
          & psb_mpi_c_dpk_,                        &
          & y%combuf(topology_total_recv + 1),     &  ! recv buffer (will contain send_indexes data)
          & neighbor_comm_handle%send_counts,       &
          & neighbor_comm_handle%send_displs,       &
          & psb_mpi_c_dpk_,                        &
          & neighbor_comm_handle%graph_comm,        &
          & neighbor_comm_handle%comm_request, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if

    end if ! do_start

    ! ---------------------------------------------------------
    !  WAIT phase: complete MPI, scatter received data
    ! ---------------------------------------------------------
    if (do_wait) then

      if (neighbor_comm_handle%comm_request == mpi_request_null) then
        ! No matching start? Something is wrong
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/-2/))
        goto 9999
      end if

      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      ! Wait for the non-blocking collective to complete
      if (debug) write(*,*) me,' nbr_vect: waiting on MPI request'
      call mpi_wait(neighbor_comm_handle%comm_request, p2pstat, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if

        ! For transpose exchange: scatter the data that correspond to peers' send area
        if (debug) write(*,*) me,' nbr_tran_vect: scattering send-index data,', topology_total_send,' elems'
        call y%sct(int(topology_total_send,psb_mpk_), &
          & neighbor_comm_handle%send_indexes, &
          & y%combuf(topology_total_recv+1:topology_total_recv+topology_total_send), &
          & beta)


      ! Clean up
      neighbor_comm_handle%comm_request = mpi_request_null
      call y%device_wait()
      call y%maybe_free_buffer(info)
      if (info /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_, name)
        goto 9999
      end if
      if (debug) write(*,*) me,' nbr_vect: done'

    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psi_ztran_neighbor_topology_vect




  !
  !
  !
  !
  ! Subroutine: psi_zswaptran_multivect
  !   Data exchange among processes.
  !
  !   Takes care of Y an encaspulated multivector.
  !   
  !   
  module subroutine psi_zswaptran_multivect(swap_status,beta,y,desc_a,info,data)

#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    integer(psb_ipk_), intent(in)                   :: swap_status
    complex(psb_dpk_), intent(in)                      :: beta
    class(psb_z_base_multivect_type), intent(inout) :: y
    type(psb_desc_type),target                      :: desc_a
    integer(psb_ipk_), intent(out)                  :: info
    integer(psb_ipk_), optional                     :: data

    ! locals
    type(psb_ctxt_type)                             :: ctxt
    integer(psb_mpk_)                               :: icomm
    integer(psb_ipk_)                               :: np, me, total_send, total_recv, num_neighbors, err_act, data_
    class(psb_i_base_vect_type), pointer            :: comm_indexes
    character(len=20)                               :: name

    integer(psb_ipk_)                               :: setflag


    info = psb_success_
    name = 'psi_zswaptran_multivect'
    call psb_erractionsave(err_act)

    ctxt = desc_a%get_context()
    icomm = ctxt%get_mpic()
    call psb_info(ctxt,me,np) 
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif

    if (.not.psb_is_asb_desc(desc_a)) then 
      info=psb_err_invalid_cd_state_
      call psb_errpush(info,name)
      goto 9999
    endif

    if (present(data)) then
      data_ = data
    else
      data_ = psb_comm_halo_
    end if

    call desc_a%get_list_p(data_,comm_indexes,num_neighbors,total_recv,total_send,info) 
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,a_err='psb_cd_get_list')
      goto 9999
    end if

    setflag = swap_status
    if (swap_status == psb_swap_start_) then
      setflag = psb_comm_status_start_
    else if (swap_status == psb_swap_wait_) then
      setflag = psb_comm_status_wait_
    else if ((iand(swap_status, psb_swap_send_) /= 0) .or. (iand(swap_status, psb_swap_recv_) /= 0) .or. &
      & (iand(swap_status, psb_swap_mpi_) /= 0) .or. (iand(swap_status, psb_swap_sync_) /= 0)) then
      setflag = psb_comm_status_sync_
    end if

    if ((setflag /= psb_comm_status_start_) .and. (setflag /= psb_comm_status_wait_) .and. &
      & (setflag /= psb_comm_status_sync_)) then
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Invalid swap_status')
      goto 9999
    end if

    if (.not. allocated(y%comm_handle)) then
      call psb_comm_set(psb_comm_isend_irecv_, y%comm_handle, info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_, name, a_err='init comm default baseline')
        goto 9999
      end if
    end if

    call y%comm_handle%set_swap_status(setflag, info)
    if (info /= psb_success_) then
      call psb_errpush(info,name,a_err='set_swap_status')
      goto 9999
    end if

    select case(y%comm_handle%comm_type)
    case(psb_comm_ineighbor_alltoallv_)
      call psi_ztran_neighbor_topology_multivect(ctxt,setflag,beta,y,comm_indexes,&
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='neighbor a2av tran')
        goto 9999
      end if
    case(psb_comm_persistent_ineighbor_alltoallv_)
      call psi_ztran_neighbor_persistent_topology_multivect(ctxt,setflag,beta,y,comm_indexes,&
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='persistent neighbor tran')
        goto 9999
      end if
    case(psb_comm_rma_pull_)
      call psi_ztran_rma_pull_multivect(ctxt,setflag,beta,y,comm_indexes,&
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='rma pull tran')
        goto 9999
      end if
    case(psb_comm_rma_push_)
      call psi_ztran_rma_push_multivect(ctxt,setflag,beta,y,comm_indexes,&
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='rma push tran')
        goto 9999
      end if
    case default
      call psi_ztran_baseline_multivect(ctxt,setflag,beta,y,comm_indexes,num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then
        call psb_errpush(info,name,a_err='baseline tran')
        goto 9999
      end if
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psi_zswaptran_multivect

  subroutine psi_ztran_baseline_multivect(ctxt,swap_status,beta,y,idx,&
      & num_neighbors,total_send,total_recv,comm_handle,info)

#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)                 :: ctxt
    integer(psb_ipk_), intent(in)                   :: swap_status
    integer(psb_ipk_), intent(out)                  :: info
    class(psb_z_base_multivect_type), intent(inout) :: y
    complex(psb_dpk_), intent(in)                      :: beta
    class(psb_i_base_vect_type), intent(inout)      :: idx
    integer(psb_ipk_), intent(in)                   :: num_neighbors,total_send,total_recv
    class(psb_comm_handle_type), intent(inout)      :: comm_handle

    ! locals
    integer(psb_mpk_)   :: np, me, nesd, nerv, n
    integer(psb_mpk_)   :: proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_), allocatable :: prcid(:)
    type(psb_comm_baseline_handle), pointer         :: baseline_comm_handle
    integer(psb_ipk_) :: err_act, i, idx_pt, total_send_, total_recv_,&
         & snd_pt, rcv_pt, pnti
    logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
         & albf,do_send,do_recv
    logical, parameter :: usersend=.false., debug=.false.
    character(len=20)  :: name

    info = psb_success_
    name = 'psi_ztran_vidx_multivect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,me,np) 
    if (np == -1) then
      info=psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif
    icomm = ctxt%get_mpic()

    baseline_comm_handle => null()
    select type(ch => comm_handle)
    type is(psb_comm_baseline_handle)
      baseline_comm_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected baseline comm_handle in baseline multivect swaptran')
      goto 9999
    end select

    n = y%get_ncols()

    do_send = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_recv = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    total_recv_ = total_recv * n
    total_send_ = total_send * n

    call idx%sync()

    if (debug) write(*,*) me,'Internal buffer'
    if (do_send) then 
      if (allocated(baseline_comm_handle%comid)) then 
        if (any(baseline_comm_handle%comid /= mpi_request_null)) then 
          ! 
          ! Unfinished communication? Something is wrong....
          !
          info=psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/-2/))
          goto 9999
        end if
      end if
      if (debug) write(*,*) me,'do_send start'
      call y%new_buffer(ione*size(idx%v),info)
      call psb_realloc(num_neighbors,2_psb_ipk_,baseline_comm_handle%comid,info)
      baseline_comm_handle%comid = mpi_request_null
      call psb_realloc(num_neighbors,prcid,info)
      ! First I post all the non blocking receives
      pnti   = 1
      snd_pt = total_recv_+1
      rcv_pt = 1
      p2ptag = psb_double_swap_tag
      do i=1, num_neighbors
        proc_to_comm = idx%v(pnti+psb_proc_id_)
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        prcid(i) = psb_get_mpi_rank(ctxt,proc_to_comm)      
        if ((nesd>0).and.(proc_to_comm /= me)) then 
          if (debug) write(*,*) me,'Posting receive from',prcid(i),snd_pt
          call mpi_irecv(y%combuf(snd_pt),n*nesd,&
               & psb_mpi_c_dpk_,prcid(i),&
            & p2ptag, icomm,baseline_comm_handle%comid(i,2),iret)
        end if
        rcv_pt = rcv_pt + n*nerv
        snd_pt = snd_pt + n*nesd
        pnti   = pnti + nerv + nesd + 3
      end do

      if (debug) write(*,*) me,' Gather '
      !
      ! Then gather for sending.
      !    
      pnti   = 1
      snd_pt = total_recv_+1
      rcv_pt = 1
      do i=1, num_neighbors
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        idx_pt = 1+pnti+psb_n_elem_recv_
        call y%gth(idx_pt,rcv_pt,nerv,idx)
        rcv_pt = rcv_pt + n*nerv
        snd_pt = snd_pt + n*nesd
        pnti   = pnti + nerv + nesd + 3
      end do

      !
      ! Then wait for device
      !
      call y%device_wait()

      if (debug) write(*,*) me,' isend'
      !
      ! Then send
      !

      pnti   = 1
      snd_pt = total_recv_+1
      rcv_pt = 1
      p2ptag = psb_double_swap_tag
      do i=1, num_neighbors
        proc_to_comm = idx%v(pnti+psb_proc_id_)
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        idx_pt = 1+pnti+psb_n_elem_recv_

        if ((nerv>0).and.(proc_to_comm /= me)) then 
          call mpi_isend(y%combuf(rcv_pt),n*nerv,&
               & psb_mpi_c_dpk_,prcid(i),&
            & p2ptag,icomm,baseline_comm_handle%comid(i,1),iret)
        end if

        if(iret /= mpi_success) then
          info=psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rcv_pt = rcv_pt + n*nerv
        snd_pt = snd_pt + n*nesd
        pnti   = pnti + nerv + nesd + 3
      end do
    end if

    if (do_recv) then 
      if (debug) write(*,*) me,' do_Recv'
      if (.not.allocated(baseline_comm_handle%comid)) then 
        ! 
        ! No matching send? Something is wrong....
        !
        info=psb_err_mpi_error_
        call psb_errpush(info,name,m_err=(/-2/))
        goto 9999
      end if
      call psb_realloc(num_neighbors,prcid,info)

      if (debug) write(*,*) me,' wait'
      pnti   = 1
      snd_pt = total_recv_+1
      rcv_pt = 1
      p2ptag = psb_double_swap_tag
      do i=1, num_neighbors
        proc_to_comm = idx%v(pnti+psb_proc_id_)
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        if (proc_to_comm /= me)then 
          if (nerv>0) then 
            call mpi_wait(baseline_comm_handle%comid(i,1),p2pstat,iret)
            if(iret /= mpi_success) then
              info=psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          end if
          if (nesd>0) then 
            call mpi_wait(baseline_comm_handle%comid(i,2),p2pstat,iret)
            if(iret /= mpi_success) then
              info=psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          end if
        else if (proc_to_comm == me) then 
          if (nesd /= nerv) then 
            write(psb_err_unit,*) &
                 & 'Fatal error in swapdata: mismatch on self send',&
                 & nerv,nesd
          end if
          y%combuf(snd_pt:snd_pt+n*nesd-1) = y%combuf(rcv_pt:rcv_pt+n*nerv-1) 
        end if
        rcv_pt = rcv_pt + n*nerv
        snd_pt = snd_pt + n*nesd
        pnti   = pnti + nerv + nesd + 3
      end do

      if (debug) write(*,*) me,' scatter'      
      pnti   = 1
      snd_pt = total_recv_+1
      rcv_pt = 1
      do i=1, num_neighbors
        proc_to_comm = idx%v(pnti+psb_proc_id_)
        nerv = idx%v(pnti+psb_n_elem_recv_)
        nesd = idx%v(pnti+nerv+psb_n_elem_send_)
        idx_pt = 1+pnti+nerv+psb_n_elem_send_

        if (debug) write(0,*)me,' Received from: ',prcid(i),&
             & y%combuf(snd_pt:snd_pt+n*nesd-1)        
        call y%sct(idx_pt,snd_pt,nesd,idx,beta)
        rcv_pt = rcv_pt + n*nerv
        snd_pt = snd_pt + n*nesd
        pnti   = pnti + nerv + nesd + 3
      end do


      !
      ! Waited for com, cleanup comid
      !
      baseline_comm_handle%comid = mpi_request_null

      !
      ! Then wait for device
      !
      if (debug) write(*,*) me,' wait'
      call y%device_wait()
      if (debug) write(*,*) me,' free buffer'
      call y%maybe_free_buffer(info)
      if (info == 0) then
        if (allocated(y%comm_handle)) call psb_comm_free(y%comm_handle, info)
      end if
      if (info /= 0) then 
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
      if (debug) write(*,*) me,' done'
    end if


    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return

  end subroutine psi_ztran_baseline_multivect


  subroutine psi_ztran_neighbor_topology_multivect(ctxt,swap_status,beta,y,comm_indexes,&
    & num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
  use mpi
#endif
  implicit none
#ifdef PSB_MPI_H
  include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)                 :: ctxt
    integer(psb_ipk_), intent(in)                   :: swap_status
    complex(psb_dpk_), intent(in)                      :: beta
    class(psb_z_base_multivect_type), intent(inout) :: y
    class(psb_i_base_vect_type), intent(inout)      :: comm_indexes
    integer(psb_ipk_), intent(in)                   :: num_neighbors,total_send,total_recv
    class(psb_comm_handle_type), intent(inout)      :: comm_handle
    integer(psb_ipk_), intent(out)                  :: info

    ! locals
    integer(psb_mpk_)                               :: icomm
    integer(psb_mpk_)                               :: np, me
    integer(psb_mpk_)                               :: iret, p2pstat(mpi_status_size)
    type(psb_comm_neighbor_handle), pointer         :: neighbor_comm_handle
    integer(psb_ipk_)                               :: err_act, topology_total_send, topology_total_recv, buffer_size
    logical                                         :: do_start, do_wait
    logical, parameter                              :: debug = .false.
    character(len=30)                               :: name


    info = psb_success_
    name = 'psi_ztran_neighbor_topology_multivect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,me,np) 
    if (np == -1) then
      info=psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif

    icomm = ctxt%get_mpic()

    neighbor_comm_handle => null()
    select type(ch => comm_handle)
    type is(psb_comm_neighbor_handle)
      neighbor_comm_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected neighbor comm_handle in neighbor multivect swaptran')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    ! ---------------------------------------------------------
    !  START phase: build topology (if needed), gather, post MPI
    ! ---------------------------------------------------------
    if (do_start) then
      if(debug) write(*,*) me,' nbr_vect: starting data exchange'
      ! Lazy initialization: build the topology on first call
      if (.not. neighbor_comm_handle%is_initialized) then
        if (debug) write(*,*) me,' nbr_vect: building topology'
        call neighbor_comm_handle%topology_init(comm_indexes%v, num_neighbors, total_send, total_recv, ctxt, icomm, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_, name, &
              & a_err='neighbor_topology_init')
          goto 9999
        end if
      end if

      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      ! Buffer layout:
      !   combuf(1 : total_send)                       = send area
      !   combuf(total_send+1 : total_send+total_recv) = recv area
      buffer_size = topology_total_send + topology_total_recv

      call y%new_buffer(buffer_size, info)
      if (info /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_, name)
        goto 9999
      end if
      neighbor_comm_handle%comm_request = mpi_request_null

        ! For transpose exchange: gather recv area first (we will send "recv" data)
        if (debug) write(*,*) me,' nbr_tran_vect: gathering recv data,', topology_total_recv,' elems'
        call y%gth(int(topology_total_recv,psb_mpk_), &
          & neighbor_comm_handle%recv_indexes, &
          & y%combuf(1:topology_total_recv))

        ! Wait for device (important for GPU subclasses)
        call y%device_wait()

        ! Post non-blocking neighborhood alltoallv swapping send/recv arrays
        if (debug) write(*,*) me,' nbr_tran_vect: posting MPI_Ineighbor_alltoallv (swapped)'
        call mpi_ineighbor_alltoallv( &
          & y%combuf(1),                           &  ! send buffer (recv_indexes gathered)
          & neighbor_comm_handle%recv_counts,       &
          & neighbor_comm_handle%recv_displs,       &
          & psb_mpi_c_dpk_,                        &
          & y%combuf(topology_total_recv + 1),     &  ! recv buffer (will contain send_indexes data)
          & neighbor_comm_handle%send_counts,       &
          & neighbor_comm_handle%send_displs,       &
          & psb_mpi_c_dpk_,                        &
          & neighbor_comm_handle%graph_comm,        &
          & neighbor_comm_handle%comm_request, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if

    end if ! do_start

    ! ---------------------------------------------------------
    !  WAIT phase: complete MPI, scatter received data
    ! ---------------------------------------------------------
    if (do_wait) then

      if (neighbor_comm_handle%comm_request == mpi_request_null) then
        ! No matching start? Something is wrong
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/-2/))
        goto 9999
      end if

      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      ! Wait for the non-blocking collective to complete
      if (debug) write(*,*) me,' nbr_vect: waiting on MPI request'
      call mpi_wait(neighbor_comm_handle%comm_request, p2pstat, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if

        ! For transpose exchange: scatter the data that correspond to peers' send area
        if (debug) write(*,*) me,' nbr_tran_vect: scattering send-index data,', topology_total_send,' elems'
        call y%sct(int(topology_total_send,psb_mpk_), &
          & neighbor_comm_handle%send_indexes, &
          & y%combuf(topology_total_recv+1:topology_total_recv+topology_total_send), &
          & beta)


      ! Clean up
      neighbor_comm_handle%comm_request = mpi_request_null
      call y%device_wait()
      call y%maybe_free_buffer(info)
      if (info /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_, name)
        goto 9999
      end if
      if (debug) write(*,*) me,' nbr_vect: done'

    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psi_ztran_neighbor_topology_multivect



  subroutine psi_ztran_neighbor_persistent_topology_vect(ctxt,swap_status,beta,y,comm_indexes,&
    & num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
  use mpi
#endif
  implicit none
#ifdef PSB_MPI_H
  include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)             :: ctxt
    integer(psb_ipk_), intent(in)               :: swap_status
    complex(psb_dpk_), intent(in)                  :: beta
    class(psb_z_base_vect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)  :: comm_indexes
    integer(psb_ipk_), intent(in)               :: num_neighbors,total_send,total_recv
    class(psb_comm_handle_type), intent(inout)  :: comm_handle
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_mpk_)                           :: icomm, np, my_rank, iret
    integer(psb_mpk_)                           :: p2pstat(mpi_status_size)
    type(psb_comm_neighbor_handle), pointer     :: neighbor_comm_handle
    integer(psb_ipk_)                           :: err_act, topology_total_send, topology_total_recv, buffer_size
    logical                                     :: do_start, do_wait
    logical, parameter                          :: debug = .false.
    character(len=30)                           :: name

    info = psb_success_
    name = 'psi_ztran_neighbor_persistent_topology_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info=psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif
    icomm = ctxt%get_mpic()

    neighbor_comm_handle => null()
    select type(ch => comm_handle)
    type is(psb_comm_neighbor_handle)
      neighbor_comm_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected neighbor comm_handle in persistent neighbor swaptran')
      goto 9999
    end select

    if (swap_status == psb_comm_status_unknown_) then
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='psb_comm_status_unknown_ not allowed in persistent neighbor swaptran')
      goto 9999
    end if

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    if (do_start) then
      if (neighbor_comm_handle%persistent_in_flight) then
        info = psb_err_mpi_error_
        call psb_errpush(info,name,a_err='Invalid START: persistent neighbor request already in flight')
        goto 9999
      end if
      if (.not. neighbor_comm_handle%is_initialized) then
        call neighbor_comm_handle%topology_init(comm_indexes%v, num_neighbors, total_send, total_recv, ctxt, icomm, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_, name, a_err='neighbor_topology_init')
          goto 9999
        end if
      end if
      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv
      buffer_size = topology_total_send + topology_total_recv

      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          if (neighbor_comm_handle%persistent_request_ready) then
            if (neighbor_comm_handle%persistent_request /= mpi_request_null) then
              call mpi_request_free(neighbor_comm_handle%persistent_request, iret)
            end if
            neighbor_comm_handle%persistent_request = mpi_request_null
            neighbor_comm_handle%persistent_request_ready = .false.
            neighbor_comm_handle%persistent_in_flight = .false.
            neighbor_comm_handle%persistent_buffer_size = 0
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= 0) then
            call psb_errpush(psb_err_alloc_dealloc_, name)
            goto 9999
          end if
        else if (size(y%combuf) < buffer_size) then
          if (neighbor_comm_handle%persistent_request_ready) then
            if (neighbor_comm_handle%persistent_request /= mpi_request_null) then
              call mpi_request_free(neighbor_comm_handle%persistent_request, iret)
            end if
            neighbor_comm_handle%persistent_request = mpi_request_null
            neighbor_comm_handle%persistent_request_ready = .false.
            neighbor_comm_handle%persistent_in_flight = .false.
            neighbor_comm_handle%persistent_buffer_size = 0
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= 0) then
            call psb_errpush(psb_err_alloc_dealloc_, name)
            goto 9999
          end if
        end if
      end if
      neighbor_comm_handle%comm_request = mpi_request_null

      if (buffer_size > 0) then
        ! Transpose: gather from recv_indexes (we "send" recv data)
        if (debug) write(*,*) my_rank,' tran_persistent_vect: gathering recv data,',topology_total_recv,' elems'
        call y%gth(int(topology_total_recv,psb_mpk_), &
          & neighbor_comm_handle%recv_indexes, &
          & y%combuf(1:topology_total_recv))
      else
        neighbor_comm_handle%persistent_in_flight = .false.
      end if

      call y%device_wait()

      if (.not. neighbor_comm_handle%persistent_request_ready) then
        if (buffer_size > 0) then
          ! Transpose: swap send/recv counts in alltoallv_init
          !   send = recv_indexes data with recv_counts/displs
          !   recv = into send_indexes area with send_counts/displs
          call mpi_neighbor_alltoallv_init( &
              & y%combuf(1),                             &
              & neighbor_comm_handle%recv_counts,        &
              & neighbor_comm_handle%recv_displs,        &
              & psb_mpi_c_dpk_,                              &
              & y%combuf(topology_total_recv + 1),       &
              & neighbor_comm_handle%send_counts,        &
              & neighbor_comm_handle%send_displs,        &
              & psb_mpi_c_dpk_,                              &
              & neighbor_comm_handle%graph_comm,         &
              & mpi_info_null,                           &
              & neighbor_comm_handle%persistent_request, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info, name, m_err=(/iret/))
            goto 9999
          end if
          neighbor_comm_handle%persistent_request_ready = .true.
          neighbor_comm_handle%persistent_buffer_size = buffer_size
        else
          neighbor_comm_handle%persistent_request_ready = .false.
          neighbor_comm_handle%persistent_buffer_size = 0
        end if
      end if

      if (buffer_size > 0) then
        call mpi_start(neighbor_comm_handle%persistent_request, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
        neighbor_comm_handle%persistent_in_flight = .true.
      else
        neighbor_comm_handle%persistent_in_flight = .false.
      end if
    end if ! do_start

    if (do_wait) then
      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      if ((topology_total_send + topology_total_recv) == 0) then
        neighbor_comm_handle%persistent_in_flight = .false.
      else
        if (.not. neighbor_comm_handle%persistent_in_flight) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, a_err='Invalid WAIT: no persistent neighbor request in flight')
          goto 9999
        end if
      end if

      if ((topology_total_send + topology_total_recv) > 0) then
        call mpi_wait(neighbor_comm_handle%persistent_request, p2pstat, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
        neighbor_comm_handle%persistent_in_flight = .false.

        ! Transpose: scatter to send_indexes
        if (debug) write(*,*) my_rank,' tran_persistent_vect: scattering to send_indexes,',topology_total_send,' elems'
        call y%sct(int(topology_total_send,psb_mpk_), &
          & neighbor_comm_handle%send_indexes, &
          & y%combuf(topology_total_recv+1:topology_total_recv+topology_total_send), &
          & beta)
      end if

      call y%device_wait()
    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_ztran_neighbor_persistent_topology_vect


  subroutine psi_ztran_rma_pull_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)               :: ctxt
    integer(psb_ipk_), intent(in)                 :: swap_status
    complex(psb_dpk_), intent(in)                    :: beta
    class(psb_z_base_vect_type), intent(inout)    :: y
    class(psb_i_base_vect_type), intent(inout)    :: comm_indexes
    integer(psb_ipk_), intent(in)                 :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)    :: comm_handle
    integer(psb_ipk_), intent(out)                :: info

    ! Effective sizes for transpose:
    !   eff_send = total_recv (we expose recv_indexes data)
    !   eff_recv = total_send (we GET into the send_indexes area)
    integer(psb_mpk_) :: np, my_rank, iret, element_bytes, icomm
    integer(psb_mpk_) :: proc_to_comm, prc_rank, recv_count, send_count, send_pos, recv_pos, list_pos
    integer(psb_mpk_) :: remote_base
    integer(kind=MPI_ADDRESS_KIND) :: remote_disp, exposed_bytes
    integer(psb_ipk_) :: err_act, neighbor_idx, buffer_size, eff_send, eff_recv
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    logical :: do_start, do_wait, layout_rebuild_needed
    type(psb_comm_rma_handle), pointer :: rma_handle
    character(len=30) :: name

    info = psb_success_
    name = 'psi_ztran_rma_pull_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    icomm = ctxt%get_mpic()
    eff_send = total_recv
    eff_recv = total_send

    select type(ch => comm_handle)
    type is(psb_comm_rma_handle)
      rma_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected RMA comm_handle for tran pull')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    if (do_start) then
      buffer_size = eff_send + eff_recv
      layout_rebuild_needed = (.not. rma_handle%layout_ready) .or. &
         & (rma_handle%layout_nnbr /= num_neighbors) .or. &
         & (rma_handle%layout_send /= eff_send) .or. &
         & (rma_handle%layout_recv /= eff_recv)

      if (layout_rebuild_needed) then
        if (allocated(peer_mpi_rank)) deallocate(peer_mpi_rank)
        if (num_neighbors > 0) then
          allocate(peer_mpi_rank(num_neighbors), stat=iret)
          if (iret /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info,name,a_err='RMA tran pull rank cache allocation')
            goto 9999
          end if
        end if
        list_pos = 1
        do neighbor_idx = 1, num_neighbors
          proc_to_comm = comm_indexes%v(list_pos+psb_proc_id_)
          peer_mpi_rank(neighbor_idx) = psb_get_mpi_rank(ctxt,proc_to_comm)
          recv_count = comm_indexes%v(list_pos+psb_n_elem_recv_)
          send_count = comm_indexes%v(list_pos+recv_count+psb_n_elem_send_)
          list_pos = list_pos + recv_count + send_count + 3
        end do
        call rma_handle%init_memory_buffer_layout_tran(info, comm_indexes%v, peer_mpi_rank, &
             & num_neighbors, total_send, total_recv, my_rank, icomm)
        if (info /= psb_success_) then
          call psb_errpush(info,name,a_err='RMA tran pull init_memory_buffer_layout_tran failure')
          goto 9999
        end if
      end if

      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
          end if
        else if (size(y%combuf) < size(comm_indexes%v)) then
          if (rma_handle%window_open) then
            call mpi_win_unlock_all(rma_handle%win, iret)
            rma_handle%window_open = .false.
          end if
          if (rma_handle%window_ready) then
            call mpi_win_free(rma_handle%win, iret)
            rma_handle%window_ready = .false.
            rma_handle%win = mpi_win_null
          end if
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
          end if
        end if
      end if

      if ((buffer_size > 0).and.(.not. rma_handle%window_ready)) then
        element_bytes = storage_size(y%combuf(1))/8
        exposed_bytes = int(size(y%combuf),kind=MPI_ADDRESS_KIND) * int(element_bytes,kind=MPI_ADDRESS_KIND)
        call mpi_win_create(y%combuf, exposed_bytes, element_bytes, &
             & mpi_info_null, ctxt%get_mpic(), rma_handle%win, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rma_handle%window_ready = .true.
      end if

      if (buffer_size > 0) then
        ! Transpose: gather recv_indexes data into combuf(1:eff_send)
        if (eff_send > 0) then
          call y%gth(int(eff_send,psb_mpk_), rma_handle%peer_send_indexes, y%combuf(1:eff_send))
        end if
        call y%device_wait()

        ! Pull from each peer's recv_indexes area with per-neighbor passive lock (neighbor-only sync).
        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          send_count = rma_handle%peer_send_counts(neighbor_idx)
          recv_count = rma_handle%peer_recv_counts(neighbor_idx)
          prc_rank   = rma_handle%peer_mpi_rank(neighbor_idx)
          send_pos   = rma_handle%peer_send_displs(neighbor_idx) + 1
          recv_pos   = eff_send + rma_handle%peer_recv_displs(neighbor_idx) + 1

          if (proc_to_comm /= my_rank) then
            remote_base = rma_handle%peer_remote_send_displs(neighbor_idx)
            if (remote_base < 1) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='Invalid remote metadata in RMA tran pull')
              goto 9999
            end if
            if (recv_count > 0) then
              call mpi_win_lock(MPI_LOCK_SHARED, prc_rank, 0, rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
              remote_disp = int(remote_base - 1, kind=MPI_ADDRESS_KIND)
              call mpi_get(y%combuf(recv_pos), recv_count, psb_mpi_r_dpk_, prc_rank, remote_disp, recv_count, psb_mpi_r_dpk_, &
                   & rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
              call mpi_win_unlock(prc_rank, rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
            end if
          else
            if (send_count /= recv_count) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA tran pull self-copy mismatch')
              goto 9999
            end if
            y%combuf(recv_pos:recv_pos+recv_count-1) = y%combuf(send_pos:send_pos+send_count-1)
          end if
        end do
      end if
    end if ! do_start

    ! WAIT phase: GETs already complete (per-neighbor unlock in START); scatter into send_indexes.
    if (do_wait) then
      ! Transpose: scatter to send_indexes (peer_recv_indexes in tran init = actual send_indexes)
      if (eff_recv > 0) then
        call y%sct(int(eff_recv,psb_mpk_), rma_handle%peer_recv_indexes, y%combuf(eff_send+1:eff_send+eff_recv), beta)
      end if
      call y%device_wait()
    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_ztran_rma_pull_vect


  subroutine psi_ztran_rma_push_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)               :: ctxt
    integer(psb_ipk_), intent(in)                 :: swap_status
    complex(psb_dpk_), intent(in)                    :: beta
    class(psb_z_base_vect_type), intent(inout)    :: y
    class(psb_i_base_vect_type), intent(inout)    :: comm_indexes
    integer(psb_ipk_), intent(in)                 :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)    :: comm_handle
    integer(psb_ipk_), intent(out)                :: info

    integer(psb_mpk_) :: np, my_rank, iret, element_bytes, icomm
    integer(psb_mpk_) :: proc_to_comm, prc_rank, recv_count, send_count, send_pos, recv_pos, list_pos
    integer(psb_mpk_) :: remote_base
    integer(kind=MPI_ADDRESS_KIND) :: remote_disp, exposed_bytes
    integer(psb_ipk_) :: err_act, neighbor_idx, buffer_size, eff_send, eff_recv
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    integer(psb_mpk_), parameter :: rma_push_notify_tag = 914_psb_mpk_
    logical :: do_start, do_wait, layout_rebuild_needed
    type(psb_comm_rma_handle), pointer :: rma_handle
    character(len=30) :: name

    info = psb_success_
    name = 'psi_ztran_rma_push_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    icomm = ctxt%get_mpic()
    eff_send = total_recv
    eff_recv = total_send

    select type(ch => comm_handle)
    type is(psb_comm_rma_handle)
      rma_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected RMA comm_handle for tran push')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    if (do_start) then
      buffer_size = eff_send + eff_recv
      layout_rebuild_needed = (.not. rma_handle%layout_ready) .or. &
         & (rma_handle%layout_nnbr /= num_neighbors) .or. &
         & (rma_handle%layout_send /= eff_send) .or. &
         & (rma_handle%layout_recv /= eff_recv)

      if (layout_rebuild_needed) then
        if (allocated(peer_mpi_rank)) deallocate(peer_mpi_rank)
        if (num_neighbors > 0) then
          allocate(peer_mpi_rank(num_neighbors), stat=iret)
          if (iret /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info,name,a_err='RMA tran push rank cache allocation')
            goto 9999
          end if
        end if
        list_pos = 1
        do neighbor_idx = 1, num_neighbors
          proc_to_comm = comm_indexes%v(list_pos+psb_proc_id_)
          peer_mpi_rank(neighbor_idx) = psb_get_mpi_rank(ctxt,proc_to_comm)
          recv_count = comm_indexes%v(list_pos+psb_n_elem_recv_)
          send_count = comm_indexes%v(list_pos+recv_count+psb_n_elem_send_)
          list_pos = list_pos + recv_count + send_count + 3
        end do
        call rma_handle%init_memory_buffer_layout_tran(info, comm_indexes%v, peer_mpi_rank, &
             & num_neighbors, total_send, total_recv, my_rank, icomm)
        if (info /= psb_success_) then
          call psb_errpush(info,name,a_err='RMA tran push init_memory_buffer_layout_tran failure')
          goto 9999
        end if
      end if

      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
          end if
        else if (size(y%combuf) < size(comm_indexes%v)) then
          if (rma_handle%window_open) then
            call mpi_win_unlock_all(rma_handle%win, iret)
            rma_handle%window_open = .false.
          end if
          if (rma_handle%window_ready) then
            call mpi_win_free(rma_handle%win, iret)
            rma_handle%window_ready = .false.
            rma_handle%win = mpi_win_null
          end if
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
          end if
        end if
      end if

      if ((buffer_size > 0).and.(.not. rma_handle%window_ready)) then
        element_bytes = storage_size(y%combuf(1))/8
        exposed_bytes = int(size(y%combuf),kind=MPI_ADDRESS_KIND) * int(element_bytes,kind=MPI_ADDRESS_KIND)
        call mpi_win_create(y%combuf, exposed_bytes, element_bytes, &
             & mpi_info_null, ctxt%get_mpic(), rma_handle%win, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rma_handle%window_ready = .true.
      end if

      if (buffer_size > 0) then
        ! Transpose: gather recv_indexes data into combuf(1:eff_send)
        if (eff_send > 0) then
          call y%gth(int(eff_send,psb_mpk_), rma_handle%peer_send_indexes, y%combuf(1:eff_send))
        end if
        call y%device_wait()

        ! Pre-post notification receives before opening the window.
        if (num_neighbors > 0) then
          rma_handle%notify_recv_reqs(1:num_neighbors) = MPI_REQUEST_NULL
          rma_handle%notify_send_reqs(1:num_neighbors) = MPI_REQUEST_NULL
        end if
        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          if (proc_to_comm /= my_rank) then
            prc_rank = rma_handle%peer_mpi_rank(neighbor_idx)
            call mpi_irecv(rma_handle%notify_buf(neighbor_idx), 1, psb_mpi_mpk_, prc_rank, &
                 & rma_push_notify_tag, icomm, rma_handle%notify_recv_reqs(neighbor_idx), iret)
            if (iret /= mpi_success) then
              info = psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          end if
        end do

        call mpi_win_lock_all(0, rma_handle%win, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rma_handle%window_open = .true.

        ! Push our recv_indexes data to each peer's send_indexes area; notify after flush.
        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          send_count = rma_handle%peer_send_counts(neighbor_idx)
          recv_count = rma_handle%peer_recv_counts(neighbor_idx)
          prc_rank   = rma_handle%peer_mpi_rank(neighbor_idx)
          send_pos   = rma_handle%peer_send_displs(neighbor_idx) + 1
          recv_pos   = eff_send + rma_handle%peer_recv_displs(neighbor_idx) + 1

          if (proc_to_comm /= my_rank) then
            remote_base = rma_handle%peer_remote_recv_displs(neighbor_idx)
            if (remote_base < 1) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='Invalid remote metadata in RMA tran push')
              goto 9999
            end if
            if (send_count > 0) then
              remote_disp = int(remote_base - 1, kind=MPI_ADDRESS_KIND)
              call mpi_put(y%combuf(send_pos), send_count, psb_mpi_r_dpk_, prc_rank, remote_disp, send_count, psb_mpi_r_dpk_, &
                   & rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
            end if
            call mpi_win_flush(prc_rank, rma_handle%win, iret)
            if (iret /= mpi_success) then
              info = psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
            call mpi_isend(rma_handle%notify_buf(neighbor_idx), 1, psb_mpi_mpk_, prc_rank, &
                 & rma_push_notify_tag, icomm, rma_handle%notify_send_reqs(neighbor_idx), iret)
            if (iret /= mpi_success) then
              info = psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          else
            if (send_count /= recv_count) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA tran push self-copy mismatch')
              goto 9999
            end if
            y%combuf(recv_pos:recv_pos+recv_count-1) = y%combuf(send_pos:send_pos+send_count-1)
          end if
        end do
      end if
    end if ! do_start

    ! WAIT phase: close epoch, wait for P2P notifications, then scatter into send_indexes.
    if (do_wait) then
      if (rma_handle%window_open) then
        call mpi_win_unlock_all(rma_handle%win, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rma_handle%window_open = .false.
      end if
      if (num_neighbors > 0) then
        call mpi_waitall(num_neighbors, rma_handle%notify_recv_reqs, MPI_STATUSES_IGNORE, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        call mpi_waitall(num_neighbors, rma_handle%notify_send_reqs, MPI_STATUSES_IGNORE, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
      end if
      ! Transpose: scatter to send_indexes (peer_recv_indexes in tran init = actual send_indexes)
      if (eff_recv > 0) then
        call y%sct(int(eff_recv,psb_mpk_), rma_handle%peer_recv_indexes, y%combuf(eff_send+1:eff_send+eff_recv), beta)
      end if
      call y%device_wait()
    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_ztran_rma_push_vect


  subroutine psi_ztran_neighbor_persistent_topology_multivect(ctxt,swap_status,beta,y,comm_indexes,&
    & num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
  use mpi
#endif
  implicit none
#ifdef PSB_MPI_H
  include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)                 :: ctxt
    integer(psb_ipk_), intent(in)                   :: swap_status
    complex(psb_dpk_), intent(in)                      :: beta
    class(psb_z_base_multivect_type), intent(inout) :: y
    class(psb_i_base_vect_type), intent(inout)      :: comm_indexes
    integer(psb_ipk_), intent(in)                   :: num_neighbors,total_send,total_recv
    class(psb_comm_handle_type), intent(inout)      :: comm_handle
    integer(psb_ipk_), intent(out)                  :: info

    integer(psb_mpk_)                               :: icomm, np, my_rank, iret, n
    integer(psb_mpk_)                               :: p2pstat(mpi_status_size)
    type(psb_comm_neighbor_handle), pointer         :: neighbor_comm_handle
    integer(psb_ipk_)                               :: err_act, topology_total_send, topology_total_recv, buffer_size
    integer(psb_mpk_)                               :: total_send_, total_recv_
    logical                                         :: do_start, do_wait
    logical, parameter                              :: debug = .false.
    character(len=30)                               :: name

    info = psb_success_
    name = 'psi_ztran_neighbor_persistent_topology_multivect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info=psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    endif
    icomm = ctxt%get_mpic()
    n = y%get_ncols()

    neighbor_comm_handle => null()
    select type(ch => comm_handle)
    type is(psb_comm_neighbor_handle)
      neighbor_comm_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected neighbor comm_handle in persistent neighbor multivect swaptran')
      goto 9999
    end select

    if (swap_status == psb_comm_status_unknown_) then
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='psb_comm_status_unknown_ not allowed in persistent neighbor multivect swaptran')
      goto 9999
    end if

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    if (do_start) then
      if (neighbor_comm_handle%persistent_in_flight) then
        info = psb_err_mpi_error_
        call psb_errpush(info,name,a_err='Invalid START: persistent neighbor request already in flight')
        goto 9999
      end if
      if (.not. neighbor_comm_handle%is_initialized) then
        call neighbor_comm_handle%topology_init(comm_indexes%v, num_neighbors, total_send, total_recv, ctxt, icomm, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_, name, a_err='neighbor_topology_init')
          goto 9999
        end if
      end if
      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv
      total_send_ = topology_total_send * n
      total_recv_ = topology_total_recv * n
      buffer_size = total_send_ + total_recv_

      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          if (neighbor_comm_handle%persistent_request_ready) then
            if (neighbor_comm_handle%persistent_request /= mpi_request_null) then
              call mpi_request_free(neighbor_comm_handle%persistent_request, iret)
            end if
            neighbor_comm_handle%persistent_request = mpi_request_null
            neighbor_comm_handle%persistent_request_ready = .false.
            neighbor_comm_handle%persistent_in_flight = .false.
            neighbor_comm_handle%persistent_buffer_size = 0
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= 0) then
            call psb_errpush(psb_err_alloc_dealloc_, name)
            goto 9999
          end if
        else if (size(y%combuf) < buffer_size) then
          if (neighbor_comm_handle%persistent_request_ready) then
            if (neighbor_comm_handle%persistent_request /= mpi_request_null) then
              call mpi_request_free(neighbor_comm_handle%persistent_request, iret)
            end if
            neighbor_comm_handle%persistent_request = mpi_request_null
            neighbor_comm_handle%persistent_request_ready = .false.
            neighbor_comm_handle%persistent_in_flight = .false.
            neighbor_comm_handle%persistent_buffer_size = 0
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= 0) then
            call psb_errpush(psb_err_alloc_dealloc_, name)
            goto 9999
          end if
        end if
      end if
      neighbor_comm_handle%comm_request = mpi_request_null

      if (buffer_size > 0) then
        ! Transpose: gather from recv_indexes
        if (debug) write(*,*) my_rank,' tran_persistent_mv: gathering recv data,',topology_total_recv,' elems'
        call y%gth(int(topology_total_recv,psb_mpk_), &
          & neighbor_comm_handle%recv_indexes, &
          & y%combuf(1:total_recv_))
      else
        neighbor_comm_handle%persistent_in_flight = .false.
      end if

      call y%device_wait()

      if (.not. neighbor_comm_handle%persistent_request_ready) then
        if (buffer_size > 0) then
          ! Transpose: swap send/recv in alltoallv_init
          call mpi_neighbor_alltoallv_init( &
              & y%combuf(1),                              &
              & n*neighbor_comm_handle%recv_counts,       &
              & n*neighbor_comm_handle%recv_displs,       &
              & psb_mpi_c_dpk_,                               &
              & y%combuf(total_recv_ + 1),                &
              & n*neighbor_comm_handle%send_counts,       &
              & n*neighbor_comm_handle%send_displs,       &
              & psb_mpi_c_dpk_,                               &
              & neighbor_comm_handle%graph_comm,          &
              & mpi_info_null,                            &
              & neighbor_comm_handle%persistent_request, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info, name, m_err=(/iret/))
            goto 9999
          end if
          neighbor_comm_handle%persistent_request_ready = .true.
          neighbor_comm_handle%persistent_buffer_size = buffer_size
        else
          neighbor_comm_handle%persistent_request_ready = .false.
          neighbor_comm_handle%persistent_buffer_size = 0
        end if
      end if

      if (buffer_size > 0) then
        call mpi_start(neighbor_comm_handle%persistent_request, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
        neighbor_comm_handle%persistent_in_flight = .true.
      else
        neighbor_comm_handle%persistent_in_flight = .false.
      end if
    end if ! do_start

    if (do_wait) then
      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv
      total_send_ = topology_total_send * n
      total_recv_ = topology_total_recv * n

      if ((topology_total_send + topology_total_recv) == 0) then
        neighbor_comm_handle%persistent_in_flight = .false.
      else
        if (.not. neighbor_comm_handle%persistent_in_flight) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, a_err='Invalid WAIT: no persistent neighbor request in flight')
          goto 9999
        end if
      end if

      if ((topology_total_send + topology_total_recv) > 0) then
        call mpi_wait(neighbor_comm_handle%persistent_request, p2pstat, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
        neighbor_comm_handle%persistent_in_flight = .false.

        ! Transpose: scatter to send_indexes
        if (debug) write(*,*) my_rank,' tran_persistent_mv: scattering to send_indexes,',topology_total_send,' elems'
        call y%sct(int(topology_total_send,psb_mpk_), &
          & neighbor_comm_handle%send_indexes, &
          & y%combuf(total_recv_+1:total_recv_+total_send_), &
          & beta)
      end if

      call y%device_wait()
    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_ztran_neighbor_persistent_topology_multivect


  subroutine psi_ztran_rma_pull_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)                   :: ctxt
    integer(psb_ipk_), intent(in)                     :: swap_status
    complex(psb_dpk_), intent(in)                        :: beta
    class(psb_z_base_multivect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)        :: comm_indexes
    integer(psb_ipk_), intent(in)                     :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)        :: comm_handle
    integer(psb_ipk_), intent(out)                    :: info

    integer(psb_mpk_) :: np, my_rank, iret, element_bytes, icomm, n
    integer(psb_mpk_) :: proc_to_comm, prc_rank, recv_count, send_count, send_pos, recv_pos, list_pos
    integer(psb_mpk_) :: remote_base
    integer(kind=MPI_ADDRESS_KIND) :: remote_disp, exposed_bytes
    integer(psb_ipk_) :: err_act, neighbor_idx, buffer_size, eff_send, eff_recv, total_eff_send_, total_eff_recv_
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    logical :: do_start, do_wait, layout_rebuild_needed
    type(psb_comm_rma_handle), pointer :: rma_handle
    character(len=30) :: name

    info = psb_success_
    name = 'psi_ztran_rma_pull_multivect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    icomm = ctxt%get_mpic()
    n = y%get_ncols()
    eff_send = total_recv
    eff_recv = total_send
    total_eff_send_ = eff_send * n
    total_eff_recv_ = eff_recv * n

    select type(ch => comm_handle)
    type is(psb_comm_rma_handle)
      rma_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected RMA comm_handle for tran pull multivect')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    if (do_start) then
      buffer_size = total_eff_send_ + total_eff_recv_
      layout_rebuild_needed = (.not. rma_handle%layout_ready) .or. &
         & (rma_handle%layout_nnbr /= num_neighbors) .or. &
         & (rma_handle%layout_send /= eff_send) .or. &
         & (rma_handle%layout_recv /= eff_recv)

      if (layout_rebuild_needed) then
        if (allocated(peer_mpi_rank)) deallocate(peer_mpi_rank)
        if (num_neighbors > 0) then
          allocate(peer_mpi_rank(num_neighbors), stat=iret)
          if (iret /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info,name,a_err='RMA tran pull multivect rank allocation')
            goto 9999
          end if
        end if
        list_pos = 1
        do neighbor_idx = 1, num_neighbors
          proc_to_comm = comm_indexes%v(list_pos+psb_proc_id_)
          peer_mpi_rank(neighbor_idx) = psb_get_mpi_rank(ctxt,proc_to_comm)
          recv_count = comm_indexes%v(list_pos+psb_n_elem_recv_)
          send_count = comm_indexes%v(list_pos+recv_count+psb_n_elem_send_)
          list_pos = list_pos + recv_count + send_count + 3
        end do
        call rma_handle%init_memory_buffer_layout_tran(info, comm_indexes%v, peer_mpi_rank, &
             & num_neighbors, total_send, total_recv, my_rank, icomm)
        if (info /= psb_success_) then
          call psb_errpush(info,name,a_err='RMA tran pull multivect init_memory_buffer_layout_tran')
          goto 9999
        end if
      end if

      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          call y%new_buffer(buffer_size, info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
          end if
        else if (size(y%combuf) < buffer_size) then
          if (rma_handle%window_open) then
            call mpi_win_unlock_all(rma_handle%win, iret)
            rma_handle%window_open = .false.
          end if
          if (rma_handle%window_ready) then
            call mpi_win_free(rma_handle%win, iret)
            rma_handle%window_ready = .false.
            rma_handle%win = mpi_win_null
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
          end if
        end if
      end if

      if ((buffer_size > 0).and.(.not. rma_handle%window_ready)) then
        element_bytes = storage_size(y%combuf(1))/8
        exposed_bytes = int(size(y%combuf),kind=MPI_ADDRESS_KIND) * int(element_bytes,kind=MPI_ADDRESS_KIND)
        call mpi_win_create(y%combuf, exposed_bytes, element_bytes, &
             & mpi_info_null, ctxt%get_mpic(), rma_handle%win, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rma_handle%window_ready = .true.
      end if

      if (buffer_size > 0) then
        if (eff_send > 0) then
          call y%gth(int(eff_send,psb_mpk_), rma_handle%peer_send_indexes, y%combuf(1:total_eff_send_))
        end if
        call y%device_wait()

        ! Pull from each peer's recv_indexes area with per-neighbor passive lock (neighbor-only sync).
        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          send_count = rma_handle%peer_send_counts(neighbor_idx)
          recv_count = rma_handle%peer_recv_counts(neighbor_idx)
          prc_rank   = rma_handle%peer_mpi_rank(neighbor_idx)
          send_pos   = rma_handle%peer_send_displs(neighbor_idx)*n + 1
          recv_pos   = total_eff_send_ + rma_handle%peer_recv_displs(neighbor_idx)*n + 1

          if (proc_to_comm /= my_rank) then
            remote_base = rma_handle%peer_remote_send_displs(neighbor_idx)
            if (remote_base < 1) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='Invalid remote metadata in RMA tran pull multivect')
              goto 9999
            end if
            if (recv_count > 0) then
              call mpi_win_lock(MPI_LOCK_SHARED, prc_rank, 0, rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
              remote_disp = int((remote_base - 1)*n, kind=MPI_ADDRESS_KIND)
              call mpi_get(y%combuf(recv_pos), recv_count*n, psb_mpi_r_dpk_, prc_rank, remote_disp, recv_count*n, psb_mpi_r_dpk_, &
                   & rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
              call mpi_win_unlock(prc_rank, rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
            end if
          else
            if (send_count /= recv_count) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA tran pull multivect self-copy mismatch')
              goto 9999
            end if
            y%combuf(recv_pos:recv_pos+recv_count*n-1) = y%combuf(send_pos:send_pos+send_count*n-1)
          end if
        end do
      end if
    end if ! do_start

    ! WAIT phase: GETs already complete (per-neighbor unlock in START); scatter into send_indexes.
    if (do_wait) then
      if (eff_recv > 0) then
        call y%sct(int(eff_recv,psb_mpk_), rma_handle%peer_recv_indexes, &
          & y%combuf(total_eff_send_+1:total_eff_send_+total_eff_recv_), beta)
      end if
      call y%device_wait()
    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_ztran_rma_pull_multivect


  subroutine psi_ztran_rma_push_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)                   :: ctxt
    integer(psb_ipk_), intent(in)                     :: swap_status
    complex(psb_dpk_), intent(in)                        :: beta
    class(psb_z_base_multivect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)        :: comm_indexes
    integer(psb_ipk_), intent(in)                     :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)        :: comm_handle
    integer(psb_ipk_), intent(out)                    :: info

    integer(psb_mpk_) :: np, my_rank, iret, element_bytes, icomm, n
    integer(psb_mpk_) :: proc_to_comm, prc_rank, recv_count, send_count, send_pos, recv_pos, list_pos
    integer(psb_mpk_) :: remote_base
    integer(kind=MPI_ADDRESS_KIND) :: remote_disp, exposed_bytes
    integer(psb_ipk_) :: err_act, neighbor_idx, buffer_size, eff_send, eff_recv, total_eff_send_, total_eff_recv_
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    integer(psb_mpk_), parameter :: rma_push_notify_tag = 914_psb_mpk_
    logical :: do_start, do_wait, layout_rebuild_needed
    type(psb_comm_rma_handle), pointer :: rma_handle
    character(len=30) :: name

    info = psb_success_
    name = 'psi_ztran_rma_push_multivect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    icomm = ctxt%get_mpic()
    n = y%get_ncols()
    eff_send = total_recv
    eff_recv = total_send
    total_eff_send_ = eff_send * n
    total_eff_recv_ = eff_recv * n

    select type(ch => comm_handle)
    type is(psb_comm_rma_handle)
      rma_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected RMA comm_handle for tran push multivect')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    if (do_start) then
      buffer_size = total_eff_send_ + total_eff_recv_
      layout_rebuild_needed = (.not. rma_handle%layout_ready) .or. &
         & (rma_handle%layout_nnbr /= num_neighbors) .or. &
         & (rma_handle%layout_send /= eff_send) .or. &
         & (rma_handle%layout_recv /= eff_recv)

      if (layout_rebuild_needed) then
        if (allocated(peer_mpi_rank)) deallocate(peer_mpi_rank)
        if (num_neighbors > 0) then
          allocate(peer_mpi_rank(num_neighbors), stat=iret)
          if (iret /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info,name,a_err='RMA tran push multivect rank allocation')
            goto 9999
          end if
        end if
        list_pos = 1
        do neighbor_idx = 1, num_neighbors
          proc_to_comm = comm_indexes%v(list_pos+psb_proc_id_)
          peer_mpi_rank(neighbor_idx) = psb_get_mpi_rank(ctxt,proc_to_comm)
          recv_count = comm_indexes%v(list_pos+psb_n_elem_recv_)
          send_count = comm_indexes%v(list_pos+recv_count+psb_n_elem_send_)
          list_pos = list_pos + recv_count + send_count + 3
        end do
        call rma_handle%init_memory_buffer_layout_tran(info, comm_indexes%v, peer_mpi_rank, &
             & num_neighbors, total_send, total_recv, my_rank, icomm)
        if (info /= psb_success_) then
          call psb_errpush(info,name,a_err='RMA tran push multivect init_memory_buffer_layout_tran')
          goto 9999
        end if
      end if

      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          call y%new_buffer(buffer_size, info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
          end if
        else if (size(y%combuf) < buffer_size) then
          if (rma_handle%window_open) then
            call mpi_win_unlock_all(rma_handle%win, iret)
            rma_handle%window_open = .false.
          end if
          if (rma_handle%window_ready) then
            call mpi_win_free(rma_handle%win, iret)
            rma_handle%window_ready = .false.
            rma_handle%win = mpi_win_null
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name)
            goto 9999
          end if
        end if
      end if

      if ((buffer_size > 0).and.(.not. rma_handle%window_ready)) then
        element_bytes = storage_size(y%combuf(1))/8
        exposed_bytes = int(size(y%combuf),kind=MPI_ADDRESS_KIND) * int(element_bytes,kind=MPI_ADDRESS_KIND)
        call mpi_win_create(y%combuf, exposed_bytes, element_bytes, &
             & mpi_info_null, ctxt%get_mpic(), rma_handle%win, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rma_handle%window_ready = .true.
      end if

      if (buffer_size > 0) then
        if (eff_send > 0) then
          call y%gth(int(eff_send,psb_mpk_), rma_handle%peer_send_indexes, y%combuf(1:total_eff_send_))
        end if
        call y%device_wait()

        ! Pre-post notification receives before opening the window.
        if (num_neighbors > 0) then
          rma_handle%notify_recv_reqs(1:num_neighbors) = MPI_REQUEST_NULL
          rma_handle%notify_send_reqs(1:num_neighbors) = MPI_REQUEST_NULL
        end if
        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          if (proc_to_comm /= my_rank) then
            prc_rank = rma_handle%peer_mpi_rank(neighbor_idx)
            call mpi_irecv(rma_handle%notify_buf(neighbor_idx), 1, psb_mpi_mpk_, prc_rank, &
                 & rma_push_notify_tag, icomm, rma_handle%notify_recv_reqs(neighbor_idx), iret)
            if (iret /= mpi_success) then
              info = psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          end if
        end do

        call mpi_win_lock_all(0, rma_handle%win, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rma_handle%window_open = .true.

        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          send_count = rma_handle%peer_send_counts(neighbor_idx)
          recv_count = rma_handle%peer_recv_counts(neighbor_idx)
          prc_rank   = rma_handle%peer_mpi_rank(neighbor_idx)
          send_pos   = rma_handle%peer_send_displs(neighbor_idx)*n + 1
          recv_pos   = total_eff_send_ + rma_handle%peer_recv_displs(neighbor_idx)*n + 1

          if (proc_to_comm /= my_rank) then
            remote_base = rma_handle%peer_remote_recv_displs(neighbor_idx)
            if (remote_base < 1) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='Invalid remote metadata in RMA tran push multivect')
              goto 9999
            end if
            if (send_count > 0) then
              remote_disp = int((remote_base - 1)*n, kind=MPI_ADDRESS_KIND)
              call mpi_put(y%combuf(send_pos), send_count*n, psb_mpi_r_dpk_, prc_rank, remote_disp, send_count*n, psb_mpi_r_dpk_, &
                   & rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
            end if
            call mpi_win_flush(prc_rank, rma_handle%win, iret)
            if (iret /= mpi_success) then
              info = psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
            call mpi_isend(rma_handle%notify_buf(neighbor_idx), 1, psb_mpi_mpk_, prc_rank, &
                 & rma_push_notify_tag, icomm, rma_handle%notify_send_reqs(neighbor_idx), iret)
            if (iret /= mpi_success) then
              info = psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          else
            if (send_count /= recv_count) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA tran push multivect self-copy mismatch')
              goto 9999
            end if
            y%combuf(recv_pos:recv_pos+recv_count*n-1) = y%combuf(send_pos:send_pos+send_count*n-1)
          end if
        end do
      end if
    end if ! do_start

    ! WAIT phase: close epoch, wait for P2P notifications, then scatter into send_indexes.
    if (do_wait) then
      if (rma_handle%window_open) then
        call mpi_win_unlock_all(rma_handle%win, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        rma_handle%window_open = .false.
      end if
      if (num_neighbors > 0) then
        call mpi_waitall(num_neighbors, rma_handle%notify_recv_reqs, MPI_STATUSES_IGNORE, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
        call mpi_waitall(num_neighbors, rma_handle%notify_send_reqs, MPI_STATUSES_IGNORE, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info,name,m_err=(/iret/))
          goto 9999
        end if
      end if
      if (eff_recv > 0) then
        call y%sct(int(eff_recv,psb_mpk_), rma_handle%peer_recv_indexes, &
          & y%combuf(total_eff_send_+1:total_eff_send_+total_eff_recv_), beta)
      end if
      call y%device_wait()
    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_ztran_rma_push_multivect


end submodule psi_z_swaptran_impl
