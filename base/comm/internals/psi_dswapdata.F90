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
! File: psi_dswapdata.F90
!
!
!
! Subroutine: psi_dswapdata_vect
!   Implements the data exchange among processes. Essentially this is doing 
!   a variable all-to-all data exchange (ALLTOALLV in MPI parlance), but 
!   it is capable of pruning empty exchanges, which are very likely in out 
!   application environment. All the variants have the same structure 
!   In all these subroutines X may be:    I    Integer
!                                         S    real(psb_spk_)
!                                         D    real(psb_dpk_)
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
! Arguments: 
!    swap_status     - integer                 Swap status selector.
!                                       It is interpreted as a communication status:
!                                         psb_comm_status_start_   -> START phase
!                                         psb_comm_status_wait_    -> WAIT phase
!                                         psb_comm_status_unknown_ -> START+WAIT
!                                       The communication scheme is selected from 
!                                       y%comm_handle%comm_type.
!
!
!    n        - integer                 Number of columns in Y               
!    beta     - real                    Choose overwrite or sum. 
!    y        - type(psb_@x@_vect_type) The data area                        
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                return code.
!    data     - integer                 which list is to be used to exchange data
!                                       default psb_comm_halo_
!                                       psb_comm_halo_    use halo_index
!                                       psb_comm_ext_     use ext_index 
!                                       psb_comm_ovrl_    use ovrl_index
!                                       psb_comm_mov_     use ovr_mst_idx
!
!
! 
submodule (psi_d_comm_v_mod)  psi_d_swapdata_impl
  use psb_desc_const_mod, only: psb_swap_start_, psb_swap_wait_
  use psb_base_mod
    use psb_comm_schemes_mod, only: psb_comm_handle_type, psb_comm_isend_irecv_, psb_comm_ineighbor_alltoallv_, &
      & psb_comm_persistent_ineighbor_alltoallv_, &
      & psb_comm_status_start_, psb_comm_status_wait_, psb_comm_status_unknown_
  use psb_comm_factory_mod, only: psb_comm_init, psb_comm_free
  use psb_comm_baseline_mod, only: psb_comm_baseline_handle, psb_comm_baseline_alloc_comid
  use psb_comm_neighbor_impl_mod, only: psb_comm_neighbor_handle
contains
  module subroutine psi_dswapdata_vect(swap_status,beta,y,desc_a,info,data)

#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    integer(psb_ipk_), intent(in)               :: swap_status
    class(psb_d_base_vect_type), intent(inout)  :: y
    real(psb_dpk_), intent(in)                  :: beta
    type(psb_desc_type), target                 :: desc_a ! TODO: should this be intent(in)?
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_), optional                 :: data

    ! locals
    type(psb_ctxt_type)                         :: ctxt
    integer(psb_ipk_)                           :: np, me, total_send, total_recv, num_neighbors, data_, err_act
    integer(psb_ipk_)                           :: setflag
    class(psb_i_base_vect_type), pointer        :: comm_indexes

    ! communication scheme/status selectors
    logical                                     :: baseline, neighbor_a2av

    ! error handling variables
    integer(psb_ipk_)                           :: err_act
    character(len=30)                           :: name

    info = psb_success_
    name = 'psi_dswapdata_vect'
    call psb_erractionsave(err_act)

    ctxt = desc_a%get_context()

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

    if(present(data)) then
      data_ = data
    else
      data_ = psb_comm_halo_
    end if

    call desc_a%get_list_p(data_,comm_indexes,num_neighbors,total_recv,total_send,info) 
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,a_err='desc_a%get_list_p')
      goto 9999
    end if

    ! Debug: report list sizes
    ! if(me == 0) then 
    !   write(psb_err_unit,*) me, 'DBG: get_list_p -> num_neighbors=', &
    !   & num_neighbors, ' total_send=', total_send, ' total_recv=', total_recv
    ! end if
    ! Accept both new comm-status enums and legacy descriptor bitfields.
    setflag = swap_status
    if (swap_status == psb_swap_start_) then
      setflag = psb_comm_status_start_
    else if (swap_status == psb_swap_wait_) then
      setflag = psb_comm_status_wait_
    else if (iand(swap_status, psb_swap_start_) /= 0 .and. iand(swap_status, psb_swap_wait_) /= 0) then
      setflag = psb_comm_status_unknown_
    end if

    if ((setflag /= psb_comm_status_start_) .and. (setflag /= psb_comm_status_wait_) .and. &
      & (setflag /= psb_comm_status_unknown_)) then
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Invalid swap_status swap_status')
      goto 9999
    end if

    if (.not. allocated(y%comm_handle)) then
      call psb_comm_init(psb_comm_isend_irecv_, y%comm_handle, info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_, name, a_err='init comm default baseline')
        goto 9999
      end if
    end if

    ! Debug: comm handle allocation and type
    ! if(me == 0) then
    !   write(psb_err_unit,*) me, 'CALLING'
    !   write(psb_err_unit,*) me, 'DBG: comm_handle allocated=', allocated(y%comm_handle)
    !   if (allocated(y%comm_handle)) then
    !     write(psb_err_unit,*) me, 'DBG: comm_handle%comm_type=', y%comm_handle%comm_type
    !   end if
    ! end if
    ! Set the normalized swap status on the comm handle
    call y%comm_handle%set_swap_status(setflag, info)
    if (info /= psb_success_) then
      call psb_errpush(info,name,a_err='set_swap_status')
      goto 9999
    end if

    ! if(me == 0) then
    !   write(psb_err_unit,*) me, 'DBG: after set_swap_status, info=', info
    ! end if

    baseline = .false.
    neighbor_a2av = .false.
    select case(y%comm_handle%comm_type)
    case(psb_comm_ineighbor_alltoallv_, psb_comm_persistent_ineighbor_alltoallv_)
      neighbor_a2av = .true.
    case default
      baseline = .true.
    end select

    ! if(me == 0) then
    !   write(psb_err_unit,*) me, 'DBG: selected baseline=', baseline, ' neighbor=', neighbor_a2av
    ! end if

    if (baseline) then 
      call psi_dswap_baseline_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='baseline swap')
        goto 9999
      end if
    else if (neighbor_a2av) then 
      call psi_dswap_neighbor_topology_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,&
        & total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='neighbor a2av swap')
        goto 9999
      end if
    else 
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Incompatible swap_status settings: neither baseline nor neighbor_a2av is true')
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_dswapdata_vect


 subroutine psi_dswap_baseline_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
  use mpi
#endif
  implicit none
#ifdef PSB_MPI_H
  include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)             :: ctxt
    integer(psb_ipk_), intent(in)               :: swap_status
    real(psb_dpk_), intent(in)                  :: beta
    class(psb_d_base_vect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)  :: comm_indexes
    integer(psb_ipk_), intent(in)               :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)  :: comm_handle
    integer(psb_ipk_), intent(out)              :: info

    ! locals
    integer(psb_mpk_)                           :: icomm
    integer(psb_mpk_) :: np, me
    integer(psb_mpk_) :: proc_to_comm, p2ptag, p2pstat(mpi_status_size),&
        & iret, nesd, nerv
    integer(psb_mpk_), allocatable :: prcid(:)
    type(psb_comm_baseline_handle), pointer     :: baseline_comm_handle
    integer(psb_ipk_) :: err_act, i, idx_pt, total_send_, total_recv_,&
        & snd_pt, rcv_pt, pnti, n
    logical :: do_send,do_recv
    logical, parameter :: usersend=.false., debug=.false.
    character(len=20)  :: name

    info = psb_success_
    name = 'psi_dswap_baseline_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,me,np) 
    if (np == -1) then
      info = psb_err_context_error_
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
      call psb_errpush(info,name,a_err='Expected baseline comm_handle in baseline swap')
      goto 9999
    end select

    n=1
    do_send = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_unknown_)
    do_recv = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_unknown_)

    total_recv_ = total_recv * n
    total_send_ = total_send * n
    call comm_indexes%sync()

    if (debug) write(*,*) me,'Internal buffer'
    if (do_send) then 
      if (allocated(baseline_comm_handle%comid)) then
        if (any(baseline_comm_handle%comid /= mpi_request_null)) then 
          ! 
          ! Unfinished communication? Something is wrong....
          !
          info=psb_err_mpi_error_
          call psb_errpush(info,name,a_err='Unfinished communication? Something is wrong....')
          goto 9999
        end if
      end if
      if (debug) write(*,*) me,'do_send start'
      call y%new_buffer(ione*size(comm_indexes%v),info)
      call psb_realloc(num_neighbors,2_psb_ipk_,baseline_comm_handle%comid,info)
      baseline_comm_handle%comid = mpi_request_null
      call psb_realloc(num_neighbors,prcid,info)
      ! First I post all the non blocking receives
      pnti   = 1
      do i=1, num_neighbors
        proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
        nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
        nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)

        rcv_pt = 1+pnti+psb_n_elem_recv_
        prcid(i) = psb_get_mpi_rank(ctxt,proc_to_comm)      
        if ((nerv>0).and.(proc_to_comm /= me)) then 
          if (debug) write(*,*) me,'Posting receive from',prcid(i),rcv_pt
          p2ptag = psb_double_swap_tag
            call mpi_irecv(y%combuf(rcv_pt),nerv,&
              & psb_mpi_r_dpk_,prcid(i),&
              & p2ptag, icomm,baseline_comm_handle%comid(i,2),iret)
        end if
        pnti   = pnti + nerv + nesd + 3
      end do
      if (debug) write(*,*) me,' Gather '
      !
      ! Then gather for sending.
      !    
      pnti   = 1
      do i=1, num_neighbors
        nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
        nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_
        idx_pt = snd_pt
        call y%gth(idx_pt,nesd,comm_indexes)
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
        proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
        nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
        nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_

        if ((nesd>0).and.(proc_to_comm /= me)) then 
            call mpi_isend(y%combuf(snd_pt),nesd,&
              & psb_mpi_r_dpk_,prcid(i),&
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
        proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
        nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
        nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_

        if (proc_to_comm /= me)then 
          if (nesd>0) then 
            call mpi_wait(baseline_comm_handle%comid(i,1),p2pstat,iret)
            if(iret /= mpi_success) then
              info=psb_err_mpi_error_
              call psb_errpush(info,name,m_err=(/iret/))
              goto 9999
            end if
          end if
          if (nerv>0) then 
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
          y%combuf(rcv_pt:rcv_pt+nerv-1) = y%combuf(snd_pt:snd_pt+nesd-1)
        end if
        pnti   = pnti + nerv + nesd + 3
      end do

      if (debug) write(*,*) me,' scatter'      
      pnti   = 1
      snd_pt = 1
      rcv_pt = 1
      do i=1, num_neighbors
        proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
        nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
        nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
        idx_pt = 1+pnti+psb_n_elem_recv_
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_

        if (debug) write(0,*)me,' Received from: ',prcid(i),&
            & y%combuf(rcv_pt:rcv_pt+nerv-1)        
        call y%sct(rcv_pt,nerv,comm_indexes,beta)
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
        if (allocated(y%comm_handle)) call y%comm_handle%free(info)
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
  end subroutine psi_dswap_baseline_vect



  subroutine psi_dswap_neighbor_topology_vect(ctxt,swap_status,beta,y,comm_indexes,&
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
    real(psb_dpk_), intent(in)                  :: beta
    class(psb_d_base_vect_type), intent(inout)  :: y
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
    name = 'psi_dswap_neighbor_topology_vect'
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
      call psb_errpush(info,name,a_err='Expected neighbor comm_handle in neighbor swap')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_unknown_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_unknown_)

    call comm_indexes%sync()

    ! ---------------------------------------------------------
    !  START phase: build topology (if needed), gather, post MPI
    ! ---------------------------------------------------------
    if (do_start) then
      if(debug) write(*,*) me,' nbr_vect: starting data exchange'
      if (.not. neighbor_comm_handle%is_initialized) then
        if (debug) write(*,*) me,' nbr_vect: building topology via handle'
        call neighbor_comm_handle%topology_init(comm_indexes%v, num_neighbors, total_send, total_recv, ctxt, icomm, info)
        if (info /= psb_success_) then
          call psb_errpush(psb_err_internal_error_, name, a_err='neighbor_topology_init')
          goto 9999
        end if
      end if
      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      ! Buffer layout:
      !   combuf(1 : total_send)                       = send area
      !   combuf(total_send+1 : total_send+total_recv) = recv area
      buffer_size = topology_total_send + topology_total_recv

      if (neighbor_comm_handle%use_persistent_buffers) then
        if ((.not.allocated(y%combuf)) .or. (size(y%combuf) < buffer_size)) then
          if (neighbor_comm_handle%persistent_request_ready) then
            if (neighbor_comm_handle%persistent_request /= mpi_request_null) then
              call mpi_request_free(neighbor_comm_handle%persistent_request, iret)
            end if
            neighbor_comm_handle%persistent_request = mpi_request_null
            neighbor_comm_handle%persistent_request_ready = .false.
            neighbor_comm_handle%persistent_buffer_size = 0
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= 0) then
            call psb_errpush(psb_err_alloc_dealloc_, name)
            goto 9999
          end if
        end if
      else
        call y%new_buffer(buffer_size, info)
        if (info /= 0) then
          call psb_errpush(psb_err_alloc_dealloc_, name)
          goto 9999
        end if
      end if
      neighbor_comm_handle%comm_request = mpi_request_null

      ! Gather send data into contiguous send buffer (polymorphic for GPU)
      if (debug) write(*,*) me,' nbr_vect: gathering send data,', topology_total_send,' elems'
      call y%gth(int(topology_total_send,psb_mpk_), &
          & neighbor_comm_handle%send_indexes, &
          & y%combuf(1:topology_total_send))

      ! Wait for device (important for GPU subclasses)
      call y%device_wait()

      if (neighbor_comm_handle%use_persistent_buffers) then
        ! Lazy persistent-init: build the request once, then reuse with START/WAIT.
        if (.not. neighbor_comm_handle%persistent_request_ready) then
#ifdef PSB_HAVE_MPI_NEIGHBOR_PERSISTENT
          if (debug) write(*,*) me,' nbr_vect: posting MPI_Neighbor_alltoallv_init'
          call mpi_neighbor_alltoallv_init( &
              & y%combuf(1),                          &  ! send buffer
              & neighbor_comm_handle%send_counts,     &
              & neighbor_comm_handle%send_displs,     &
              & psb_mpi_r_dpk_,                       &
              & y%combuf(topology_total_send + 1),    &  ! recv buffer
              & neighbor_comm_handle%recv_counts,     &
              & neighbor_comm_handle%recv_displs,     &
              & psb_mpi_r_dpk_,                       &
              & neighbor_comm_handle%graph_comm,      &
              & mpi_info_null,                        &
              & neighbor_comm_handle%persistent_request, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info, name, m_err=(/iret/))
            goto 9999
          end if
          neighbor_comm_handle%persistent_request_ready = .true.
          neighbor_comm_handle%persistent_buffer_size = buffer_size
#else
          ! Fallback when persistent neighborhood collectives are not available
          neighbor_comm_handle%persistent_request_ready = .false.
          neighbor_comm_handle%persistent_buffer_size = 0
#endif
        end if

#ifdef PSB_HAVE_MPI_NEIGHBOR_PERSISTENT
        call mpi_start(neighbor_comm_handle%persistent_request, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
#else
        call mpi_ineighbor_alltoallv( &
            & y%combuf(1),                        &  ! send buffer
            & neighbor_comm_handle%send_counts,     &
            & neighbor_comm_handle%send_displs,     &
            & psb_mpi_r_dpk_,                     &
            & y%combuf(topology_total_send + 1),       &  ! recv buffer
            & neighbor_comm_handle%recv_counts,     &
            & neighbor_comm_handle%recv_displs,     &
            & psb_mpi_r_dpk_,                     &
            & neighbor_comm_handle%graph_comm,      &
            & neighbor_comm_handle%comm_request, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
#endif
      else
        ! Post non-blocking neighborhood alltoallv
        if (debug) write(*,*) me,' nbr_vect: posting MPI_Ineighbor_alltoallv'
        call mpi_ineighbor_alltoallv( &
            & y%combuf(1),                        &  ! send buffer
            & neighbor_comm_handle%send_counts,     &
            & neighbor_comm_handle%send_displs,     &
            & psb_mpi_r_dpk_,                     &
            & y%combuf(topology_total_send + 1),       &  ! recv buffer
            & neighbor_comm_handle%recv_counts,     &
            & neighbor_comm_handle%recv_displs,     &
            & psb_mpi_r_dpk_,                     &
            & neighbor_comm_handle%graph_comm,      &
            & neighbor_comm_handle%comm_request, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
      end if

    end if ! do_start

    ! ---------------------------------------------------------
    !  WAIT phase: complete MPI, scatter received data
    ! ---------------------------------------------------------
    if (do_wait) then

      if (neighbor_comm_handle%use_persistent_buffers) then
#ifdef PSB_HAVE_MPI_NEIGHBOR_PERSISTENT
        if (.not. neighbor_comm_handle%persistent_request_ready) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/-2/))
          goto 9999
        end if
#else
        if (neighbor_comm_handle%comm_request == mpi_request_null) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/-2/))
          goto 9999
        end if
#endif
      else
        if (neighbor_comm_handle%comm_request == mpi_request_null) then
          write(psb_err_unit,*) me, 'DBG: neighbor WAIT but comm_request is NULL; is_initialized=', &
          & neighbor_comm_handle%is_initialized
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/-2/))
          goto 9999
        end if
      end if

      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      ! Wait for the non-blocking collective to complete
      if (debug) write(*,*) me,' nbr_vect: waiting on MPI request'
      if (neighbor_comm_handle%use_persistent_buffers) then
#ifdef PSB_HAVE_MPI_NEIGHBOR_PERSISTENT
        call mpi_wait(neighbor_comm_handle%persistent_request, p2pstat, iret)
#else
        call mpi_wait(neighbor_comm_handle%comm_request, p2pstat, iret)
#endif
      else
        call mpi_wait(neighbor_comm_handle%comm_request, p2pstat, iret)
      end if
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if

      ! Scatter received data to local vector positions (polymorphic for GPU)
      if (debug) write(*,*) me,' nbr_vect: scattering recv data,', topology_total_recv,' elems'
      call y%sct(int(topology_total_recv,psb_mpk_), &
          & neighbor_comm_handle%recv_indexes, &
          & y%combuf(topology_total_send+1:topology_total_send+topology_total_recv), &
          & beta)


      ! Clean up
      if ((.not. neighbor_comm_handle%use_persistent_buffers) .or. &
        & (neighbor_comm_handle%use_persistent_buffers .and. .not. neighbor_comm_handle%persistent_request_ready)) then
        neighbor_comm_handle%comm_request = mpi_request_null
      end if
      call y%device_wait()
      if (.not. neighbor_comm_handle%use_persistent_buffers) then
        call y%maybe_free_buffer(info)
        if (info /= 0) then
          call psb_errpush(psb_err_alloc_dealloc_, name)
          goto 9999
        end if
      end if
      if (debug) write(*,*) me,' nbr_vect: done'

    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psi_dswap_neighbor_topology_vect



  !
  !
  ! Subroutine: psi_dswapdata_multivect
  !   Data exchange among processes.
  !
  !   Takes care of Y an encaspulated multivector.
  !   
  !   
  module subroutine psi_dswapdata_multivect(swap_status,beta,y,desc_a,info,data)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    integer(psb_ipk_), intent(in)                   :: swap_status
    class(psb_d_base_multivect_type), intent(inout) :: y
    real(psb_dpk_), intent(in)                      :: beta
    type(psb_desc_type), target                     :: desc_a
    integer(psb_ipk_), intent(out)                  :: info
    integer(psb_ipk_), optional                     :: data

    ! communication scheme/status selectors
    logical                                         :: baseline, neighbor_a2av

    ! locals
    type(psb_ctxt_type)                             :: ctxt
    integer(psb_mpk_)                               :: icomm
    integer(psb_ipk_)                               :: np, me, total_send, total_recv, num_neighbors, data_, err_act
    class(psb_i_base_vect_type), pointer            :: comm_indexes
    character(len=30)                               :: name

    
    info = psb_success_
    name = 'psi_dswapdata_multivect'
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

    if(present(data)) then
      data_ = data
    else
      data_ = psb_comm_halo_
    end if

    call desc_a%get_list_p(data_,comm_indexes,num_neighbors,total_recv,total_send,info) 
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,a_err='psb_cd_get_list')
      goto 9999
    end if

    if ((swap_status /= psb_comm_status_start_) .and. (swap_status /= psb_comm_status_wait_) &
      & .and. (swap_status /= psb_comm_status_unknown_)) then
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Invalid swap_status swap_status')
      goto 9999
    end if

    if (.not. allocated(y%comm_handle)) then
      call psb_comm_init(psb_comm_isend_irecv_, y%comm_handle, info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_, name, a_err='init comm default baseline')
        goto 9999
      end if
    end if

    call y%comm_handle%set_swap_status(swap_status, info)
    if (info /= psb_success_) then
      call psb_errpush(info,name,a_err='set_swap_status')
      goto 9999
    end if

    baseline = .false.
    neighbor_a2av = .false.
    select case(y%comm_handle%comm_type)
    case(psb_comm_ineighbor_alltoallv_, psb_comm_persistent_ineighbor_alltoallv_)
      neighbor_a2av = .true.
    case default
      baseline = .true.
    end select

    if (baseline) then 
      call psi_dswap_baseline_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='baseline swap')
        goto 9999
      end if
    else if (neighbor_a2av) then 
       call psi_dswap_neighbor_topology_multivect(ctxt,swap_status,beta,y,comm_indexes, &
         & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='neighbor a2av swap')
        goto 9999
      end if
    else 
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Incompatible swap_status settings: neither baseline nor neighbor_a2av is true')
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psi_dswapdata_multivect




subroutine psi_dswap_baseline_multivect(ctxt,swap_status,beta,y,comm_indexes, &
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
  real(psb_dpk_), intent(in)                      :: beta
  class(psb_d_base_multivect_type), intent(inout) :: y
  class(psb_i_base_vect_type), intent(inout)      :: comm_indexes
  integer(psb_ipk_), intent(in)                   :: num_neighbors,total_send, total_recv
  class(psb_comm_handle_type), intent(inout)      :: comm_handle
  integer(psb_ipk_), intent(out)                  :: info

  ! locals
  integer(psb_mpk_)                           :: icomm
  integer(psb_mpk_)                           :: np, me, nesd, nerv, n
  integer(psb_mpk_)                           :: proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret
  integer(psb_mpk_), allocatable              :: prcid(:)
  type(psb_comm_baseline_handle), pointer     :: baseline_comm_handle
  integer(psb_ipk_)                           :: err_act, i, idx_pt, total_send_, total_recv_,&
       & snd_pt, rcv_pt, pnti
  logical                                     :: do_send,do_recv
  logical, parameter                          :: usersend=.false., debug=.false.
  character(len=20)                           :: name

  info = psb_success_
  name = 'psi_dswap_baseline_multivect'
  call psb_erractionsave(err_act)
  call psb_info(ctxt,me,np) 
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  n = y%get_ncols()

  baseline_comm_handle => null()
  select type(ch => comm_handle)
  type is(psb_comm_baseline_handle)
    baseline_comm_handle => ch
  class default
    info = psb_err_mpi_error_
    call psb_errpush(info,name,a_err='Expected baseline comm_handle in baseline multivect swap')
    goto 9999
  end select

  do_send = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_unknown_)
  do_recv = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_unknown_)

  total_recv_ = total_recv * n
  total_send_ = total_send * n

  call comm_indexes%sync()

  if (debug) write(*,*) me,'Internal buffer'
  if (do_send) then 
    if (allocated(baseline_comm_handle%comid)) then
      if (any(baseline_comm_handle%comid /= mpi_request_null)) then
        info = psb_err_mpi_error_
        call psb_errpush(info,name,m_err=(/-2/))
        goto 9999
      end if
    end if
    if (debug) write(*,*) me,'do_send start'
    call y%new_buffer(ione*size(comm_indexes%v),info)
    call psb_realloc(num_neighbors,2_psb_ipk_,baseline_comm_handle%comid,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
    baseline_comm_handle%comid = mpi_request_null
    call psb_realloc(num_neighbors,prcid,info)
    ! First I post all the non blocking receives
    pnti   = 1
    snd_pt = total_recv_+1
    rcv_pt = 1
    do i=1, num_neighbors
      proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
      nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
      nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
      prcid(i) = psb_get_mpi_rank(ctxt,proc_to_comm)      
      if ((nerv>0).and.(proc_to_comm /= me)) then 
        if (debug) write(*,*) me,'Posting receive from',prcid(i),rcv_pt
        p2ptag = psb_double_swap_tag
           call mpi_irecv(y%combuf(rcv_pt),n*nerv,&
             & psb_mpi_r_dpk_,prcid(i),&
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
      nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
      nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+nerv+psb_n_elem_send_
      call y%gth(idx_pt,snd_pt,nesd,comm_indexes)
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
      proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
      nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
      nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)

      if ((nesd>0).and.(proc_to_comm /= me)) then 
           call mpi_isend(y%combuf(snd_pt),n*nesd,&
             & psb_mpi_r_dpk_,prcid(i),&
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
      proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
      nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
      nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
      if (proc_to_comm /= me)then 
        if (nesd>0) then 
          call mpi_wait(baseline_comm_handle%comid(i,1),p2pstat,iret)
          if(iret /= mpi_success) then
            info=psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
        end if
        if (nerv>0) then 
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
        y%combuf(rcv_pt:rcv_pt+n*nerv-1) = y%combuf(snd_pt:snd_pt+n*nesd-1)
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
      proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
      nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
      nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_

      if (debug) write(0,*)me,' Received from: ',prcid(i),&
           & y%combuf(rcv_pt:rcv_pt+n*nerv-1)        
      call y%sct(idx_pt,rcv_pt,nerv,comm_indexes,beta)
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
    call y%free_buffer(info)
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
end subroutine psi_dswap_baseline_multivect



subroutine psi_dswap_neighbor_topology_multivect(ctxt,swap_status,beta,y,comm_indexes,&
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
  real(psb_dpk_), intent(in)                      :: beta
  class(psb_d_base_multivect_type), intent(inout) :: y
  class(psb_i_base_vect_type), intent(inout)      :: comm_indexes
  integer(psb_ipk_), intent(in)                   :: num_neighbors,total_send, total_recv
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
  name = 'psi_dswap_neighbor_topology_multivect'
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
    call psb_errpush(info,name,a_err='Expected neighbor comm_handle in neighbor multivect swap')
    goto 9999
  end select

  do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_unknown_)
  do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_unknown_)

  call comm_indexes%sync()

  ! ---------------------------------------------------------
  !  START phase: build topology (if needed), gather, post MPI
  ! ---------------------------------------------------------
  if (do_start) then
    if(debug) write(*,*) me,' nbr_vect: starting data exchange'
    if (.not. neighbor_comm_handle%is_initialized) then
      if (debug) write(*,*) me,' nbr_vect: building topology via handle'
      call neighbor_comm_handle%topology_init(comm_indexes%v, num_neighbors, total_send, total_recv, &
          & ctxt, icomm, info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_, name, a_err='neighbor_topology_init')
        goto 9999
      end if
    end if
    topology_total_send = neighbor_comm_handle%total_send
    topology_total_recv = neighbor_comm_handle%total_recv

    ! Buffer layout:
    !   combuf(1 : total_send)                       = send area
    !   combuf(total_send+1 : total_send+total_recv) = recv area
    buffer_size = topology_total_send + topology_total_recv

    if (neighbor_comm_handle%use_persistent_buffers) then
      if ((.not.allocated(y%combuf)) .or. (size(y%combuf) < buffer_size)) then
        if (neighbor_comm_handle%persistent_request_ready) then
          if (neighbor_comm_handle%persistent_request /= mpi_request_null) then
            call mpi_request_free(neighbor_comm_handle%persistent_request, iret)
          end if
          neighbor_comm_handle%persistent_request = mpi_request_null
          neighbor_comm_handle%persistent_request_ready = .false.
          neighbor_comm_handle%persistent_buffer_size = 0
        end if
        call y%new_buffer(buffer_size, info)
        if (info /= 0) then
          call psb_errpush(psb_err_alloc_dealloc_, name)
          goto 9999
        end if
      end if
    else
      call y%new_buffer(buffer_size, info)
      if (info /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_, name)
        goto 9999
      end if
    end if
    neighbor_comm_handle%comm_request = mpi_request_null

    ! Gather send data into contiguous send buffer (polymorphic for GPU)
    if (debug) write(*,*) me,' nbr_vect: gathering send data,', topology_total_send,' elems'
    call y%gth(int(topology_total_send,psb_mpk_), &
      & neighbor_comm_handle%send_indexes, &
      & y%combuf(1:topology_total_send))

    ! Wait for device (important for GPU subclasses)
    call y%device_wait()

    if (neighbor_comm_handle%use_persistent_buffers) then
      if (.not. neighbor_comm_handle%persistent_request_ready) then
#ifdef PSB_HAVE_MPI_NEIGHBOR_PERSISTENT
        if (debug) write(*,*) me,' nbr_vect: posting MPI_Neighbor_alltoallv_init'
        call mpi_neighbor_alltoallv_init( &
            & y%combuf(1),                          &  ! send buffer
            & neighbor_comm_handle%send_counts,     &
            & neighbor_comm_handle%send_displs,     &
            & psb_mpi_r_dpk_,                       &
            & y%combuf(topology_total_send + 1),    &  ! recv buffer
            & neighbor_comm_handle%recv_counts,     &
            & neighbor_comm_handle%recv_displs,     &
            & psb_mpi_r_dpk_,                       &
            & neighbor_comm_handle%graph_comm,      &
            & mpi_info_null,                        &
            & neighbor_comm_handle%persistent_request, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
        neighbor_comm_handle%persistent_request_ready = .true.
        neighbor_comm_handle%persistent_buffer_size = buffer_size
#else
        ! Fallback when persistent neighborhood collectives are not available
        neighbor_comm_handle%persistent_request_ready = .false.
        neighbor_comm_handle%persistent_buffer_size = 0
#endif
      end if

#ifdef PSB_HAVE_MPI_NEIGHBOR_PERSISTENT
      call mpi_start(neighbor_comm_handle%persistent_request, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if
#else
      call mpi_ineighbor_alltoallv( &
          & y%combuf(1),                        &  ! send buffer
          & neighbor_comm_handle%send_counts,     &
          & neighbor_comm_handle%send_displs,     &
          & psb_mpi_r_dpk_,                     &
          & y%combuf(topology_total_send + 1),       &  ! recv buffer
          & neighbor_comm_handle%recv_counts,     &
          & neighbor_comm_handle%recv_displs,     &
          & psb_mpi_r_dpk_,                     &
          & neighbor_comm_handle%graph_comm,      &
          & neighbor_comm_handle%comm_request, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if
#endif
    else
      ! Post non-blocking neighborhood alltoallv
      if (debug) write(*,*) me,' nbr_vect: posting MPI_Ineighbor_alltoallv'
      call mpi_ineighbor_alltoallv( &
          & y%combuf(1),                        &  ! send buffer
          & neighbor_comm_handle%send_counts,     &
          & neighbor_comm_handle%send_displs,     &
          & psb_mpi_r_dpk_,                     &
          & y%combuf(topology_total_send + 1),       &  ! recv buffer
          & neighbor_comm_handle%recv_counts,     &
          & neighbor_comm_handle%recv_displs,     &
          & psb_mpi_r_dpk_,                     &
          & neighbor_comm_handle%graph_comm,      &
          & neighbor_comm_handle%comm_request, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if
    end if

  end if ! do_start

  ! ---------------------------------------------------------
  !  WAIT phase: complete MPI, scatter received data
  ! ---------------------------------------------------------
  if (do_wait) then

    if (neighbor_comm_handle%use_persistent_buffers) then
#ifdef PSB_HAVE_MPI_NEIGHBOR_PERSISTENT
      if (.not. neighbor_comm_handle%persistent_request_ready) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/-2/))
        goto 9999
      end if
#else
      if (neighbor_comm_handle%comm_request == mpi_request_null) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/-2/))
        goto 9999
      end if
#endif
    else
      if (neighbor_comm_handle%comm_request == mpi_request_null) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/-2/))
        goto 9999
      end if
    end if

    topology_total_send = neighbor_comm_handle%total_send
    topology_total_recv = neighbor_comm_handle%total_recv

    ! Wait for the non-blocking collective to complete
    if (debug) write(*,*) me,' nbr_vect: waiting on MPI request'
    if (neighbor_comm_handle%use_persistent_buffers) then
#ifdef PSB_HAVE_MPI_NEIGHBOR_PERSISTENT
      call mpi_wait(neighbor_comm_handle%persistent_request, p2pstat, iret)
#else
      call mpi_wait(neighbor_comm_handle%comm_request, p2pstat, iret)
#endif
    else
      call mpi_wait(neighbor_comm_handle%comm_request, p2pstat, iret)
    end if
    if (iret /= mpi_success) then
      info = psb_err_mpi_error_
      call psb_errpush(info, name, m_err=(/iret/))
      goto 9999
    end if

    ! Scatter received data to local vector positions (polymorphic for GPU)
    if (debug) write(*,*) me,' nbr_vect: scattering recv data,', topology_total_recv,' elems'
    call y%sct(int(topology_total_recv,psb_mpk_), &
      & neighbor_comm_handle%recv_indexes, &
      & y%combuf(topology_total_send+1:topology_total_send+topology_total_recv), &
      & beta)


    ! Clean up
    if ((.not. neighbor_comm_handle%use_persistent_buffers) .or. &
      & (neighbor_comm_handle%use_persistent_buffers .and. .not. neighbor_comm_handle%persistent_request_ready)) then
      neighbor_comm_handle%comm_request = mpi_request_null
    end if
    call y%device_wait()
    if (.not. neighbor_comm_handle%use_persistent_buffers) then
      call y%maybe_free_buffer(info)
      if (info /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_, name)
        goto 9999
      end if
    end if
    if (debug) write(*,*) me,' nbr_vect: done'

  end if ! do_wait

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psi_dswap_neighbor_topology_multivect


end submodule psi_d_swapdata_impl

