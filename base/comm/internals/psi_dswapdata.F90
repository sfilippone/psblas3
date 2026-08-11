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
!         software without specific prior written permission.
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
!   so that special versions (i.e. GPU vectors can override them) 
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
  use psb_error_mod, only: psb_get_debug_level, psb_get_debug_unit, psb_debug_ext_
  use psb_comm_schemes_mod, only: psb_comm_isend_irecv_, psb_comm_ineighbor_alltoallv_, &
      & psb_comm_persistent_ineighbor_alltoallv_, psb_comm_rma_pull_, psb_comm_rma_push_, &
      & psb_comm_handle_type
  use psb_comm_rma_mod, only: psb_comm_rma_handle, psb_comm_rma_get_window
  use, intrinsic :: iso_c_binding, only: c_loc
  use psb_comm_factory_mod

contains

  module subroutine psi_dswapdata_vect(swap_status,beta,y,desc_a,info,data)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    integer(psb_ipk_), intent(in)                 :: swap_status
    class(psb_d_base_vect_type), intent(inout)    :: y
    real(psb_dpk_), intent(in)                    :: beta
    type(psb_desc_type), target                   :: desc_a ! TODO: should this be intent(in)?
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_), optional                   :: data

    ! locals
    type(psb_ctxt_type)                           :: ctxt
    integer(psb_ipk_)                             :: np, my_rank, total_send, total_recv, num_neighbors, data_
    class(psb_i_base_vect_type), pointer          :: comm_indexes

    ! communication scheme/status selectors
    logical                                       :: baseline, ineighbor_a2av, ineighbor_a2av_persistent

    ! error handling variables
    integer(psb_ipk_)                             :: err_act
    character(len=30)                             :: name
    logical                                       :: debug

    info = psb_success_
    name = 'psi_dswapdata_vect'
    call psb_erractionsave(err_act)

    debug = .false.

    ctxt = desc_a%get_context()

    call psb_info(ctxt,my_rank,np) 
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


    if( (swap_status /= psb_comm_status_start_).and.(swap_status /= psb_comm_status_wait_)&
      & .and.(swap_status /= psb_comm_status_sync_) ) then
        info = psb_err_mpi_error_
        call psb_errpush(info,name,a_err='Invalid swap_status swap_status')
        goto 9999
    end if

    if (.not. allocated(y%comm_handle)) then
      call psb_comm_set(desc_a%comm_type, y%comm_handle, info)
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

    if(debug .eqv. .true.) then 
      write(0,'(a,i0,a,i0)') 'DBG_DISPATCH comm_type=',y%comm_handle%comm_type,' desc_a%comm_type=',desc_a%comm_type
    endif
    select case(y%comm_handle%comm_type)
    case(psb_comm_isend_irecv_)
      call psi_dswap_baseline_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='baseline swap')
        goto 9999
      end if
    case(psb_comm_ineighbor_alltoallv_)
      call psi_dswap_neighbor_topology_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,&
        & total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='neighbor nonblocking swap')
        goto 9999
      end if
    case(psb_comm_persistent_ineighbor_alltoallv_)
      call psi_dswap_neighbor_persistent_topology_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,&
        & total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='neighbor persistent swap')
        goto 9999
      end if
    case(psb_comm_rma_pull_)
      call psi_dswap_rma_pull_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,&
           & y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='rma pull swap')
        goto 9999
      end if
    case(psb_comm_rma_push_)
      call psi_dswap_rma_push_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,&
           & y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='rma push swap')
        goto 9999
      end if
    case default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Incompatible swap_status settings: no valid communication mode selected')
      goto 9999
    end select

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
    integer(psb_mpk_) :: np, my_rank
    integer(psb_mpk_) :: proc_to_comm, p2ptag, p2pstat(mpi_status_size),&
        & iret, nesd, nerv
    integer(psb_mpk_), allocatable :: prcid(:)
    type(psb_comm_baseline_handle), pointer     :: baseline_comm_handle
    integer(psb_ipk_) :: err_act, i, idx_pt, total_send_, total_recv_,&
        & snd_pt, rcv_pt, pnti, n
    logical :: do_send,do_recv
    logical, parameter :: usersend=.false.
    logical                                     :: debug
    character(len=20)                           :: name

    info = psb_success_
    name = 'psi_dswap_baseline_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np) 
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

    if(swap_status == psb_comm_status_unknown_) then
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Invalid swap_status: psb_comm_status_unknown_ is not allowed in neighbor swap')
      goto 9999
    end if

    n=1
    do_send = (swap_status == psb_comm_status_start_).or.(swap_status == psb_comm_status_sync_) 
    do_recv = (swap_status == psb_comm_status_wait_).or.(swap_status == psb_comm_status_sync_)  

    total_recv_ = total_recv * n
    total_send_ = total_send * n
    call comm_indexes%sync()

    if (debug) write(*,*) my_rank,'Internal buffer'
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
        if ((nerv>0).and.(proc_to_comm /= my_rank)) then 
          if (debug) write(*,*) my_rank,'Posting receive from',prcid(i),rcv_pt
          p2ptag = psb_double_swap_tag
          call mpi_irecv(y%combuf(rcv_pt),nerv,&
            & psb_mpi_r_dpk_,prcid(i),&
            & p2ptag, icomm,baseline_comm_handle%comid(i,2),iret)
        end if
        pnti   = pnti + nerv + nesd + 3
      end do
      if (debug) write(*,*) my_rank,' Gather '
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
        if ((idx_pt < 1) .or. (nesd < 0) .or. (idx_pt+max(0,nesd)-1 > size(comm_indexes%v))) then
          info = psb_err_internal_error_
          call psb_errpush(info,name,a_err='baseline gather metadata out of bounds')
          goto 9999
        end if
        if ((idx_pt < 1) .or. (nesd < 0) .or. (idx_pt+max(0,nesd)-1 > size(y%combuf))) then
          info = psb_err_internal_error_
          call psb_errpush(info,name,a_err='baseline gather combuf bounds error')
          goto 9999
        end if
        call y%gth(idx_pt,nesd,comm_indexes)
        pnti   = pnti + nerv + nesd + 3
      end do

      !
      ! Then wait 
      !
      call y%device_wait()

      if (debug) write(*,*) my_rank,' isend'
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

        if ((nesd>0).and.(proc_to_comm /= my_rank)) then 
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
      if (debug) write(*,*) my_rank,' do_Recv'
      if (.not.allocated(baseline_comm_handle%comid)) then 
        ! 
        ! No matching send? Something is wrong....
        !
        info=psb_err_mpi_error_
        call psb_errpush(info,name,m_err=(/-2/))
        goto 9999
      end if
      call psb_realloc(num_neighbors,prcid,info)

      if (debug) write(*,*) my_rank,' wait'
      pnti   = 1
      p2ptag = psb_double_swap_tag
      do i=1, num_neighbors
        proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
        nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
        nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
        snd_pt = 1+pnti+nerv+psb_n_elem_send_
        rcv_pt = 1+pnti+psb_n_elem_recv_

        if (proc_to_comm /= my_rank)then 
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
        else if (proc_to_comm == my_rank) then 
          if (nesd /= nerv) then 
            write(psb_err_unit,*) &
                & 'Fatal error in swapdata: mismatch on self send',&
                & nerv,nesd
          end if
          y%combuf(rcv_pt:rcv_pt+nerv-1) = y%combuf(snd_pt:snd_pt+nesd-1)
        end if
        pnti = pnti + nerv + nesd + 3
      end do

      if (debug) write(*,*) my_rank,' scatter'      
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

        if ((idx_pt < 1) .or. (nerv < 0) .or. (idx_pt+max(0,nerv)-1 > size(comm_indexes%v))) then
          info = psb_err_internal_error_
          call psb_errpush(info,name,a_err='baseline scatter metadata out of bounds')
          goto 9999
        end if
        if ((rcv_pt < 1) .or. (nerv < 0) .or. (rcv_pt+max(0,nerv)-1 > size(y%combuf))) then
          info = psb_err_internal_error_
          call psb_errpush(info,name,a_err='baseline scatter combuf bounds error')
          goto 9999
        end if

        if (debug) write(*,*)my_rank,' Received from: ',prcid(i),&
            & y%combuf(rcv_pt:rcv_pt+nerv-1)        
        call y%sct(rcv_pt,nerv,comm_indexes,beta)
        pnti = pnti + nerv + nesd + 3
      end do

      !
      ! Waited for everybody, clean up
      !
      baseline_comm_handle%comid = mpi_request_null

      !
      ! Then wait for device
      !
      if (debug) write(*,*) my_rank,' wait'
      call y%device_wait()
      if (debug) write(*,*) my_rank,' free buffer'
      call y%maybe_free_buffer(info)
      if (info == 0) then
        if (allocated(y%comm_handle)) call y%comm_handle%free(info)
      end if
      if (info /= 0) then 
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
      if (debug) write(*,*) my_rank,' done'
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

    type(psb_ctxt_type), intent(in)               :: ctxt
    integer(psb_ipk_), intent(in)                 :: swap_status
    real(psb_dpk_), intent(in)                   :: beta
    class(psb_d_base_vect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)    :: comm_indexes
    integer(psb_ipk_), intent(in)                 :: num_neighbors,total_send,total_recv
    class(psb_comm_handle_type), intent(inout)    :: comm_handle
    integer(psb_ipk_), intent(out)                :: info

    ! locals
    integer(psb_mpk_)                             :: icomm
    integer(psb_mpk_)                             :: np, my_rank
    integer(psb_mpk_)                             :: iret, p2pstat(mpi_status_size)
    type(psb_comm_neighbor_handle), pointer       :: neighbor_comm_handle
    integer(psb_ipk_)                             :: err_act, topology_total_send, topology_total_recv, buffer_size
    logical                                       :: do_start, do_wait
    logical                                       :: debug
    character(len=30)                             :: name
    integer(psb_ipk_)                             :: i, nesd, nerv, snd_pt, rcv_pt, pnti, n

    info = psb_success_
    name = 'psi_dswap_neighbor_topology_vect'
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
      call psb_errpush(info,name,a_err='Expected neighbor comm_handle in neighbor swap')
      goto 9999
    end select

    if(swap_status == psb_comm_status_unknown_) then
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Invalid swap_status: psb_comm_status_unknown_ is not allowed in neighbor swap')
      goto 9999
    end if

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    ! ---------------------------------------------------------
    !  START phase: build topology (if needed), gather, post MPI
    ! ---------------------------------------------------------
    if (do_start) then
      if(debug) write(*,*) my_rank,' nbr_vect: starting data exchange (nonblocking)'
      if (.not. neighbor_comm_handle%is_initialized) then
        if (debug) write(*,*) my_rank,' nbr_vect: building topology via handle'
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

      if (buffer_size > 0) then
        call y%new_buffer(ione*size(comm_indexes%v), info)
        if (info /= 0) then
          call psb_errpush(psb_err_alloc_dealloc_, name)
          goto 9999
        end if
        neighbor_comm_handle%comm_request = mpi_request_null
        pnti   = 1
        do i=1, num_neighbors
          nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
          nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
          snd_pt = 1+pnti+nerv+psb_n_elem_send_
          rcv_pt = 1+pnti+psb_n_elem_recv_
          if ((snd_pt < 1) .or. (nesd < 0) .or. (snd_pt+max(0,nesd)-1 > size(comm_indexes%v))) then
            info = psb_err_internal_error_
            call psb_errpush(info,name,a_err='baseline gather metadata out of bounds')
            goto 9999
          end if
          if ((snd_pt < 1) .or. (nesd < 0) .or. (snd_pt+max(0,nesd)-1 > size(y%combuf))) then
            info = psb_err_internal_error_
            call psb_errpush(info,name,a_err='baseline gather combuf bounds error')
            goto 9999
          end if
          call y%gth(snd_pt,nesd,comm_indexes)
          pnti   = pnti + nerv + nesd + 3
        end do
      else
        neighbor_comm_handle%comm_request = mpi_request_null
      end if

      ! Wait for device (important for GPU subclasses)
      call y%device_wait()

      ! Post non-blocking neighborhood alltoallv
      if (debug) write(*,*) my_rank,' nbr_vect: posting MPI_Ineighbor_alltoallv'
      if (buffer_size > 0) then
        call mpi_ineighbor_alltoallv( &
          & y%combuf(1),                             &  ! send buffer (baseline/positional layout)
          & neighbor_comm_handle%send_counts,        &
          & neighbor_comm_handle%send_displs_v,      &  ! positional displacements (match gth layout)
          & psb_mpi_r_dpk_,                             &
          & y%combuf(1),                             &  ! recv buffer (baseline/positional layout)
          & neighbor_comm_handle%recv_counts,        &
          & neighbor_comm_handle%recv_displs_v,      &  ! positional displacements (match sct layout)
          & psb_mpi_r_dpk_,                             &
          & neighbor_comm_handle%graph_comm,         &
          & neighbor_comm_handle%comm_request, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
      else
        neighbor_comm_handle%comm_request = mpi_request_null
      end if

    end if ! do_start

    ! ---------------------------------------------------------
    !  WAIT phase: complete MPI, scatter received data
    ! ---------------------------------------------------------
    if (do_wait) then

      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      if ((topology_total_send + topology_total_recv) == 0) then
        ! Valid no-op exchange: nothing was posted in START and nothing to wait/scatter.
        neighbor_comm_handle%comm_request = mpi_request_null
      else
        if (neighbor_comm_handle%comm_request == mpi_request_null) then
          write(psb_err_unit,*) my_rank, 'DBG: neighbor WAIT but comm_request is NULL; is_initialized=', &
          & neighbor_comm_handle%is_initialized
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/-2/))
          goto 9999
        end if
      end if

      ! Only wait and scatter if there's data
      if ((topology_total_send + topology_total_recv) > 0) then
        ! Wait for the non-blocking collective to complete
        if (debug) write(*,*) my_rank,' nbr_vect: waiting on MPI request'
        call mpi_wait(neighbor_comm_handle%comm_request, p2pstat, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if

        ! Scatter received data to local vector positions (polymorphic for GPU)
        if (debug) write(*,*) my_rank,' nbr_vect: scattering recv data,', topology_total_recv,' elems'
        pnti   = 1
        snd_pt = 1
        rcv_pt = 1
        do i=1, num_neighbors
          nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
          nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
          snd_pt = 1+pnti+nerv+psb_n_elem_send_
          rcv_pt = 1+pnti+psb_n_elem_recv_

          if ((1+pnti+psb_n_elem_recv_ < 1) .or. (nerv < 0) .or. &
          & (1+pnti+psb_n_elem_recv_+max(0,nerv)-1 > size(comm_indexes%v))) then
            info = psb_err_internal_error_
            call psb_errpush(info,name,a_err='baseline scatter metadata out of bounds')
            goto 9999
          end if
          if ((rcv_pt < 1) .or. (nerv < 0) .or. (rcv_pt+max(0,nerv)-1 > size(y%combuf))) then
            info = psb_err_internal_error_
            call psb_errpush(info,name,a_err='baseline scatter combuf bounds error')
            goto 9999
          end if
          call y%sct(rcv_pt,nerv,comm_indexes,beta)
          pnti = pnti + nerv + nesd + 3
        end do
      else
        ! nothing to wait/scatter
      end if


      ! Clean up
      neighbor_comm_handle%comm_request = mpi_request_null
      call y%device_wait()
      call y%maybe_free_buffer(info)
      if (info /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_, name)
        goto 9999
      end if
      if (debug) write(*,*) my_rank,' nbr_vect: done'

    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psi_dswap_neighbor_topology_vect



  subroutine psi_dswap_neighbor_persistent_topology_vect(ctxt,swap_status,beta,y,comm_indexes,&
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
    integer(psb_mpk_)                           :: np, my_rank
    integer(psb_mpk_)                           :: iret, p2pstat(mpi_status_size)
    type(psb_comm_neighbor_handle), pointer     :: neighbor_comm_handle
    integer(psb_ipk_)                           :: err_act, topology_total_send, topology_total_recv, buffer_size
    logical                                     :: do_start, do_wait
    logical                                     :: debug
    character(len=30)                           :: name
    integer(psb_ipk_)                           :: i, nesd, nerv, snd_pt, rcv_pt, pnti, n

    info = psb_success_
    debug = .false.
    name = 'psi_dswap_neighbor_persistent_topology_vect'
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
      call psb_errpush(info,name,a_err='Expected neighbor comm_handle in persistent neighbor swap')
      goto 9999
    end select

    if(swap_status == psb_comm_status_unknown_) then
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Invalid swap_status: psb_comm_status_unknown_ is not allowed in neighbor swap')
      goto 9999
    end if

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    ! ---------------------------------------------------------
    !  START phase: build topology (if needed), gather, post MPI
    ! ---------------------------------------------------------
    if (do_start) then
      if(debug) write(*,*) my_rank,' nbr_vect: starting data exchange (persistent)'
      if (neighbor_comm_handle%persistent_in_flight) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, a_err='Invalid START: persistent neighbor request already in flight')
        goto 9999
      end if
      if (.not. neighbor_comm_handle%is_initialized) then
        if (debug) write(*,*) my_rank,' nbr_vect: building topology via handle'
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
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= 0) then
            call psb_errpush(psb_err_alloc_dealloc_, name)
            goto 9999
          end if
        else if (size(y%combuf) < ione*size(comm_indexes%v)) then
          if (neighbor_comm_handle%persistent_request_ready) then
            if (neighbor_comm_handle%persistent_request /= mpi_request_null) then
              call mpi_request_free(neighbor_comm_handle%persistent_request, iret)
            end if
            neighbor_comm_handle%persistent_request = mpi_request_null
            neighbor_comm_handle%persistent_request_ready = .false.
            neighbor_comm_handle%persistent_in_flight = .false.
            neighbor_comm_handle%persistent_buffer_size = 0
          end if
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= 0) then
            call psb_errpush(psb_err_alloc_dealloc_, name)
            goto 9999
          end if
        end if
      end if
      neighbor_comm_handle%comm_request = mpi_request_null

      if (buffer_size > 0) then
        ! Gather send data into combuf using baseline-style gth calls
        if (debug) write(*,*) my_rank,' nbr_vect: gathering send data (baseline layout),', topology_total_send,' elems'
        pnti = 1
        do i=1, num_neighbors
          nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
          nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
          snd_pt = 1+pnti+nerv+psb_n_elem_send_
          if ((snd_pt < 1) .or. (nesd < 0) .or. (snd_pt+max(0,nesd)-1 > size(comm_indexes%v))) then
            info = psb_err_internal_error_
            call psb_errpush(info,name,a_err='baseline gather metadata out of bounds')
            goto 9999
          end if
          if ((snd_pt < 1) .or. (nesd < 0) .or. (snd_pt+max(0,nesd)-1 > size(y%combuf))) then
            info = psb_err_internal_error_
            call psb_errpush(info,name,a_err='baseline gather combuf bounds error')
            goto 9999
          end if
          call y%gth(snd_pt,nesd,comm_indexes)
          pnti = pnti + nerv + nesd + 3
        end do
      else
        neighbor_comm_handle%persistent_in_flight = .false.
      end if

      ! Wait for device (important for GPU subclasses)
      call y%device_wait()

      ! Lazy persistent-init: build the request once, then reuse with START/WAIT.
      if (.not. neighbor_comm_handle%persistent_request_ready) then
        if (buffer_size > 0) then
          if (debug) write(*,*) my_rank,' nbr_vect: posting MPI_Neighbor_alltoallv_init'
          call mpi_neighbor_alltoallv_init( &
              & y%combuf(1),                            &  ! send buffer
              & neighbor_comm_handle%send_counts,       &
              & neighbor_comm_handle%send_displs_v,     &  ! positional (baseline layout)
              & psb_mpi_r_dpk_,                         &
              & y%combuf(1),                            &  ! recv buffer (baseline layout)
              & neighbor_comm_handle%recv_counts,       &
              & neighbor_comm_handle%recv_displs_v,     &  ! positional (baseline layout)
              & psb_mpi_r_dpk_,                         &
              & neighbor_comm_handle%graph_comm,        &
              & mpi_info_null,                          &
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

    ! ---------------------------------------------------------
    !  WAIT phase: complete MPI, scatter received data
    ! ---------------------------------------------------------
    if (do_wait) then

      topology_total_send = neighbor_comm_handle%total_send
      topology_total_recv = neighbor_comm_handle%total_recv

      if ((topology_total_send + topology_total_recv) == 0) then
        ! Valid no-op exchange: nothing was posted in START and nothing to wait/scatter.
        neighbor_comm_handle%persistent_in_flight = .false.
      else
        if (.not. neighbor_comm_handle%persistent_in_flight) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, a_err='Invalid WAIT: no persistent neighbor request in flight')
          goto 9999
        end if
      end if

      ! Only wait and scatter if there's data
      if ((topology_total_send + topology_total_recv) > 0) then
        ! Wait for the persistent collective to complete
        if (debug) write(*,*) my_rank,' nbr_vect: waiting on persistent MPI request'
        call mpi_wait(neighbor_comm_handle%persistent_request, p2pstat, iret)
        if (iret /= mpi_success) then
          info = psb_err_mpi_error_
          call psb_errpush(info, name, m_err=(/iret/))
          goto 9999
        end if
        neighbor_comm_handle%persistent_in_flight = .false.

        ! Scatter received data to local vector positions (baseline-style sct)
        if (debug) write(*,*) my_rank,' nbr_vect: scattering recv data (baseline layout),', topology_total_recv,' elems'
        pnti = 1
        do i=1, num_neighbors
          nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
          nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
          rcv_pt = 1+pnti+psb_n_elem_recv_
          if (nerv > 0) then
            call y%sct(rcv_pt,nerv,comm_indexes,beta)
          end if
          pnti = pnti + nerv + nesd + 3
        end do
      end if

      call y%device_wait()
      if (debug) write(*,*) my_rank,' nbr_vect: done'

    end if ! do_wait

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psi_dswap_neighbor_persistent_topology_vect


  subroutine psi_dswap_rma_pull_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)               :: ctxt
    integer(psb_ipk_), intent(in)                 :: swap_status
    real(psb_dpk_), intent(in)                    :: beta
    class(psb_d_base_vect_type), intent(inout), target    :: y
    class(psb_i_base_vect_type), intent(inout)    :: comm_indexes
    integer(psb_ipk_), intent(in)                 :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)    :: comm_handle
    integer(psb_ipk_), intent(out)                :: info

    integer(psb_mpk_) :: np, my_rank, iret, element_bytes, icomm
    integer(psb_mpk_) :: proc_to_comm, prc_rank, recv_count, send_count, send_pos, recv_pos, list_pos
    integer(psb_mpk_) :: remote_base
    integer(kind=MPI_ADDRESS_KIND) :: remote_disp, exposed_bytes
    integer(psb_ipk_) :: err_act, neighbor_idx, buffer_size
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    integer(psb_mpk_) :: comm_grp, nbr_grp, n_off, grp_idx
    integer(psb_mpk_), allocatable :: grp_ranks(:)
    logical :: do_start, do_wait, memory_buffer_layout_rebuild_needed
    type(psb_comm_rma_handle), pointer :: rma_handle
    character(len=30) :: name
    integer(psb_ipk_)                           :: i, nesd, nerv, snd_pt, rcv_pt, pnti, n


    info = psb_success_
    name = 'psi_dswap_rma_pull_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    icomm = ctxt%get_mpic()

    select type(ch => comm_handle)
    type is(psb_comm_rma_handle)
      rma_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected RMA comm_handle for pull mode')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    ! START phase: build the layout once, prepare the exposed buffer, then issue RMA.
    if (do_start) then
      buffer_size = total_send + total_recv

      memory_buffer_layout_rebuild_needed = (.not. rma_handle%layout_ready) .or. &
         & (rma_handle%layout_nnbr /= num_neighbors) .or. &
         & (rma_handle%layout_send /= total_send) .or. &
         & (rma_handle%layout_recv /= total_recv)

      if (memory_buffer_layout_rebuild_needed) then
        if (allocated(peer_mpi_rank)) deallocate(peer_mpi_rank)
        if (num_neighbors > 0) then
          allocate(peer_mpi_rank(num_neighbors), stat=iret)
          if (iret /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info,name,a_err='RMA pull rank cache allocation')
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

           call rma_handle%init_memory_buffer_layout(info, comm_indexes%v, peer_mpi_rank, &
             & num_neighbors, total_send, total_recv, my_rank, icomm)
        if (info /= psb_success_) then
          call psb_errpush(info,name,a_err='RMA pull init_memory_buffer_layout failure')
          goto 9999
        end if
      end if

      ! Expose y%combuf as the RMA window (created once, reused). On GPU vectors
      ! combuf is the pinned/registered buffer, so RMA uses pinned memory.
      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name); goto 9999
          end if
        else if (size(y%combuf) < size(comm_indexes%v)) then
          if (rma_handle%window_open) then
            call mpi_win_unlock_all(rma_handle%win, iret)
            rma_handle%window_open = .false.
          end if
          ! Dynamic window: detach the old buffer and keep the window. Detach
          ! is a local call, so growing the buffer no longer costs a collective.
          if (rma_handle%buf_attached) then
            call mpi_win_detach(rma_handle%win, y%combuf, iret)
            rma_handle%buf_attached = .false.
          end if
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name); goto 9999
          end if
        end if

        element_bytes = storage_size(y%combuf(1))/8
        if (.not. rma_handle%window_ready) then
          ! The window comes from the module pool, keyed by communicator: it is
          ! created on the first request of the run and reused by every handle
          ! afterwards. Creating it per handle, as before, meant paying a
          ! communicator-wide collective (~0.9 s at 448 ranks) each time.
          call psb_comm_rma_get_window(ctxt%get_mpic(), rma_handle%win, info)
          if (info /= psb_success_) then
            call psb_errpush(info,name); goto 9999
          end if
          rma_handle%window_ready = .true.
        end if

        if (.not. rma_handle%buf_attached) then
          exposed_bytes = int(size(y%combuf),kind=MPI_ADDRESS_KIND) * int(element_bytes,kind=MPI_ADDRESS_KIND)
          call mpi_win_attach(rma_handle%win, y%combuf, exposed_bytes, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
          call mpi_get_address(y%combuf, rma_handle%my_buf_addr, iret)
          ! Kept so that the handle can detach on free, when the buffer is no
          ! longer reachable from here.
          rma_handle%win_base  = c_loc(y%combuf(1))
          rma_handle%win_nelem = size(y%combuf)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
          ! Targets are addressed absolutely on a dynamic window, so neighbours
          ! must be told the new base before any Get/Put can reference it.
          call rma_handle%publish_buf_addr(ctxt%get_mpic(), info)
          if (info /= psb_success_) then
            call psb_errpush(info,name); goto 9999
          end if
          rma_handle%buf_attached = .true.
        end if
      end if

      if (buffer_size > 0) then
        if (total_send > 0) then
          call y%gth(int(total_send,psb_mpk_), rma_handle%peer_send_indexes, y%combuf(1:total_send))
        end if
        call y%device_wait()

        ! Open exposure (post) + access (start) epochs over the CACHED neighbor
        ! group (built once with the layout in the handle) -- no per-call build.
        if (rma_handle%nbr_grp /= mpi_group_null) then
          call mpi_win_post(rma_handle%nbr_grp, 0, rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
          call mpi_win_start(rma_handle%nbr_grp, 0, rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
        end if

        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          recv_count = rma_handle%peer_recv_counts(neighbor_idx)
          send_count = rma_handle%peer_send_counts(neighbor_idx)
          prc_rank = rma_handle%peer_mpi_rank(neighbor_idx)
          send_pos = rma_handle%peer_send_displs(neighbor_idx) + 1
          recv_pos = total_send + rma_handle%peer_recv_displs(neighbor_idx) + 1

          if (proc_to_comm /= my_rank) then
            remote_base = rma_handle%peer_remote_send_displs(neighbor_idx)
            if (remote_base < 1) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='Invalid remote metadata in RMA pull')
              goto 9999
            end if
            if ((recv_pos < 1) .or. (recv_count < 0) .or. (recv_pos+max(0,recv_count)-1 > size(y%combuf))) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA pull local receive bounds error')
              goto 9999
            end if

            if (recv_count > 0) then
              ! Dynamic window: disp_unit is 1 and the displacement is an
              ! absolute address in the target's space, so the offset has to be
              ! expressed in bytes rather than in window elements.
              remote_disp = rma_handle%peer_buf_addr(neighbor_idx) &
                   & + int(remote_base - 1, kind=MPI_ADDRESS_KIND) &
                   &   * int(element_bytes, kind=MPI_ADDRESS_KIND)
              call mpi_get(y%combuf(recv_pos), recv_count, psb_mpi_r_dpk_, prc_rank, remote_disp, recv_count, psb_mpi_r_dpk_, &
                   & rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
            end if
          else
            if (send_count /= recv_count) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA pull self-copy mismatch')
              goto 9999
            end if
            y%combuf(recv_pos:recv_pos+recv_count-1) = y%combuf(send_pos:send_pos+send_count-1)
          end if
        end do

        ! Close access epoch (my GETs done) then exposure epoch (peers done
        ! reading my window): received data is now in y%combuf. Cached group, no free.
        if (rma_handle%nbr_grp /= mpi_group_null) then
          call mpi_win_complete(rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
          call mpi_win_wait(rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
        end if
      end if
    end if

    ! WAIT phase: GETs already complete (per-neighbor unlock in START); scatter received data into Y.
    if (do_wait) then
      if (total_recv > 0) then
        call y%sct(int(total_recv,psb_mpk_), rma_handle%peer_recv_indexes, y%combuf(total_send+1:total_send+total_recv), beta)
      end if

      call y%device_wait()
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_dswap_rma_pull_vect


  subroutine psi_dswap_rma_push_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
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
    class(psb_d_base_vect_type), intent(inout), target  :: y
    class(psb_i_base_vect_type), intent(inout)  :: comm_indexes
    integer(psb_ipk_), intent(in)               :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)  :: comm_handle
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_mpk_) :: np, my_rank, iret, element_bytes, icomm
    integer(psb_mpk_) :: proc_to_comm, prc_rank, recv_count, send_count, send_pos, recv_pos, list_pos
    integer(psb_mpk_) :: remote_base
    integer(kind=MPI_ADDRESS_KIND) :: remote_disp, exposed_bytes
    integer(psb_ipk_) :: err_act, neighbor_idx, buffer_size
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    integer(psb_mpk_) :: comm_grp, nbr_grp, n_off, grp_idx
    integer(psb_mpk_), allocatable :: grp_ranks(:)
    logical :: do_start, do_wait, memory_buffer_layout_rebuild_needed
    type(psb_comm_rma_handle), pointer :: rma_handle
    character(len=30) :: name

    info = psb_success_
    name = 'psi_dswap_rma_push_vect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    icomm = ctxt%get_mpic()

    select type(ch => comm_handle)
    type is(psb_comm_rma_handle)
      rma_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected RMA comm_handle for push mode')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    ! START phase: identical layout handling, but the remote metadata describes the receive side.
    if (do_start) then
      buffer_size = total_send + total_recv

      memory_buffer_layout_rebuild_needed = (.not. rma_handle%layout_ready) .or. &
         & (rma_handle%layout_nnbr /= num_neighbors) .or. &
         & (rma_handle%layout_send /= total_send) .or. &
         & (rma_handle%layout_recv /= total_recv)

      if (memory_buffer_layout_rebuild_needed) then
        if (allocated(peer_mpi_rank)) deallocate(peer_mpi_rank)
        if (num_neighbors > 0) then
          allocate(peer_mpi_rank(num_neighbors), stat=iret)
          if (iret /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info,name,a_err='RMA put rank cache allocation')
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

           call rma_handle%init_memory_buffer_layout(info, comm_indexes%v, peer_mpi_rank, &
             & num_neighbors, total_send, total_recv, my_rank, icomm)
        if (info /= psb_success_) then
          call psb_errpush(info,name,a_err='RMA put ini_memory_buffer_layout')
          goto 9999
        end if
      end if

      ! Expose y%combuf as the RMA window (created once, reused). On GPU vectors
      ! combuf is the pinned/registered buffer, so RMA uses pinned memory.
      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name); goto 9999
          end if
        else if (size(y%combuf) < size(comm_indexes%v)) then
          if (rma_handle%window_open) then
            call mpi_win_unlock_all(rma_handle%win, iret)
            rma_handle%window_open = .false.
          end if
          ! Dynamic window: detach the old buffer and keep the window. Detach
          ! is a local call, so growing the buffer no longer costs a collective.
          if (rma_handle%buf_attached) then
            call mpi_win_detach(rma_handle%win, y%combuf, iret)
            rma_handle%buf_attached = .false.
          end if
          call y%new_buffer(ione*size(comm_indexes%v), info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name); goto 9999
          end if
        end if

        element_bytes = storage_size(y%combuf(1))/8
        if (.not. rma_handle%window_ready) then
          ! The window comes from the module pool, keyed by communicator: it is
          ! created on the first request of the run and reused by every handle
          ! afterwards. Creating it per handle, as before, meant paying a
          ! communicator-wide collective (~0.9 s at 448 ranks) each time.
          call psb_comm_rma_get_window(ctxt%get_mpic(), rma_handle%win, info)
          if (info /= psb_success_) then
            call psb_errpush(info,name); goto 9999
          end if
          rma_handle%window_ready = .true.
        end if

        if (.not. rma_handle%buf_attached) then
          exposed_bytes = int(size(y%combuf),kind=MPI_ADDRESS_KIND) * int(element_bytes,kind=MPI_ADDRESS_KIND)
          call mpi_win_attach(rma_handle%win, y%combuf, exposed_bytes, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
          call mpi_get_address(y%combuf, rma_handle%my_buf_addr, iret)
          ! Kept so that the handle can detach on free, when the buffer is no
          ! longer reachable from here.
          rma_handle%win_base  = c_loc(y%combuf(1))
          rma_handle%win_nelem = size(y%combuf)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
          ! Targets are addressed absolutely on a dynamic window, so neighbours
          ! must be told the new base before any Get/Put can reference it.
          call rma_handle%publish_buf_addr(ctxt%get_mpic(), info)
          if (info /= psb_success_) then
            call psb_errpush(info,name); goto 9999
          end if
          rma_handle%buf_attached = .true.
        end if
      end if

      if (buffer_size > 0) then
        if (total_send > 0) then
          call y%gth(int(total_send,psb_mpk_), rma_handle%peer_send_indexes, y%combuf(1:total_send))
        end if
        call y%device_wait()

        ! Open exposure (post) + access (start) epochs over the CACHED neighbor
        ! group (built once with the layout in the handle) -- no per-call build.
        if (rma_handle%nbr_grp /= mpi_group_null) then
          call mpi_win_post(rma_handle%nbr_grp, 0, rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
          call mpi_win_start(rma_handle%nbr_grp, 0, rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
        end if

        ! Issue every PUT (remote) or perform the self-copy (local).
        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          recv_count = rma_handle%peer_recv_counts(neighbor_idx)
          send_count = rma_handle%peer_send_counts(neighbor_idx)
          prc_rank = rma_handle%peer_mpi_rank(neighbor_idx)
          send_pos = rma_handle%peer_send_displs(neighbor_idx) + 1
          recv_pos = total_send + rma_handle%peer_recv_displs(neighbor_idx) + 1

          if (proc_to_comm /= my_rank) then
            remote_base = rma_handle%peer_remote_recv_displs(neighbor_idx)
            if (remote_base < 1) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='Invalid remote metadata in RMA push')
              goto 9999
            end if
            if ((send_pos < 1) .or. (send_count < 0) .or. (send_pos+max(0,send_count)-1 > size(y%combuf))) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA push local send bounds error')
              goto 9999
            end if

            if (send_count > 0) then
              ! Dynamic window: disp_unit is 1 and the displacement is an
              ! absolute address in the target's space, so the offset has to be
              ! expressed in bytes rather than in window elements.
              remote_disp = rma_handle%peer_buf_addr(neighbor_idx) &
                   & + int(remote_base - 1, kind=MPI_ADDRESS_KIND) &
                   &   * int(element_bytes, kind=MPI_ADDRESS_KIND)
              call mpi_put(y%combuf(send_pos), send_count, psb_mpi_r_dpk_, prc_rank, remote_disp, send_count, psb_mpi_r_dpk_, &
                   & rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
            end if
          else
            if (send_count /= recv_count) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA push self-copy mismatch')
              goto 9999
            end if
            y%combuf(recv_pos:recv_pos+recv_count-1) = y%combuf(send_pos:send_pos+send_count-1)
          end if
        end do

        ! Close access epoch (my PUTs done) then exposure epoch (peers' PUTs into
        ! me done): receive area now holds the pushed data. Cached group, no free.
        if (rma_handle%nbr_grp /= mpi_group_null) then
          call mpi_win_complete(rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
          call mpi_win_wait(rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
        end if
      end if
    end if

    ! WAIT phase: data already delivered by PSCW in START; just scatter into Y.
    if (do_wait) then
      if (total_recv > 0) then
        call y%sct(int(total_recv,psb_mpk_), rma_handle%peer_recv_indexes, y%combuf(total_send+1:total_send+total_recv), beta)
      end if

      call y%device_wait()
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_dswap_rma_push_vect



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

    integer(psb_ipk_), intent(in)                     :: swap_status
    class(psb_d_base_multivect_type), intent(inout) :: y
    real(psb_dpk_), intent(in)                       :: beta
    type(psb_desc_type), target                       :: desc_a
    integer(psb_ipk_), intent(out)                    :: info
    integer(psb_ipk_), optional                       :: data

    ! communication scheme/status selectors
    logical                                           :: baseline, ineighbor_a2av, ineighbor_a2av_persistent

    ! locals
    type(psb_ctxt_type)                               :: ctxt
    integer(psb_mpk_)                                 :: icomm
    integer(psb_ipk_)                                 :: np, my_rank, total_send, total_recv, num_neighbors, data_, err_act
      integer(psb_mpk_)                               :: n, total_send_, total_recv_
    class(psb_i_base_vect_type), pointer              :: comm_indexes
    character(len=30)                                 :: name

    
    info = psb_success_
    name = 'psi_dswapdata_multivect'
    call psb_erractionsave(err_act)

    ctxt = desc_a%get_context()
    icomm = ctxt%get_mpic()
    call psb_info(ctxt,my_rank,np) 
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
      call psb_comm_set(desc_a%comm_type, y%comm_handle, info)
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

    n = y%get_ncols()
    total_send_ = total_send * n
    total_recv_ = total_recv * n

    baseline = .false.
    ineighbor_a2av = .false.
    ineighbor_a2av_persistent = .false.
    select case(y%comm_handle%comm_type)
    case(psb_comm_ineighbor_alltoallv_)
      ineighbor_a2av = .true.
    case(psb_comm_persistent_ineighbor_alltoallv_)
      ineighbor_a2av_persistent = .true.
    case(psb_comm_rma_pull_)
      call psi_dswap_rma_pull_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,&
           & y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='rma pull swap')
        goto 9999
      end if
    case(psb_comm_rma_push_)
      call psi_dswap_rma_push_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,&
           & y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='rma push swap')
        goto 9999
      end if
    case default
      baseline = .true.
    end select

    if (baseline) then 
      call psi_dswap_baseline_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='baseline swap')
        goto 9999
      end if
    else if (ineighbor_a2av) then 
       call psi_dswap_neighbor_topology_multivect(ctxt,swap_status,beta,y,comm_indexes, &
         & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='neighbor nonblocking swap')
        goto 9999
      end if
    else if (ineighbor_a2av_persistent) then 
      call psi_dswap_neighbor_topology_multivect_persistent(ctxt,swap_status,beta,y,comm_indexes, &
        & num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='neighbor persistent swap')
        goto 9999
      end if
    else 
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Incompatible swap_status settings: no valid communication mode selected')
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
  integer(psb_mpk_)                           :: np, my_rank, nesd, nerv, n
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
  call psb_info(ctxt,my_rank,np) 
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

  if (debug) write(*,*) my_rank,'Internal buffer'
  if (do_send) then 
    if (allocated(baseline_comm_handle%comid)) then
      if (any(baseline_comm_handle%comid /= mpi_request_null)) then
        info = psb_err_mpi_error_
        call psb_errpush(info,name,m_err=(/-2/))
        goto 9999
      end if
    end if
    if (debug) write(*,*) my_rank,'do_send start'
    call y%new_buffer(total_send_+total_recv_,info)
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
      if ((nerv>0).and.(proc_to_comm /= my_rank)) then 
        if (debug) write(*,*) my_rank,'Posting receive from',prcid(i),rcv_pt
        p2ptag = psb_double_swap_tag
           call mpi_irecv(y%combuf(rcv_pt),n*nerv,&
             & psb_mpi_r_dpk_,prcid(i),&
             & p2ptag, icomm,baseline_comm_handle%comid(i,2),iret)
      end if
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3
    end do
    if (debug) write(*,*) my_rank,' Gather '
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

    if (debug) write(*,*) my_rank,' isend'
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

      if ((nesd>0).and.(proc_to_comm /= my_rank)) then 
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
    if (debug) write(*,*) my_rank,' do_Recv'
    if (.not.allocated(baseline_comm_handle%comid)) then 
      ! 
      ! No matching send? Something is wrong....
      !
      info=psb_err_mpi_error_
      call psb_errpush(info,name,m_err=(/-2/))
      goto 9999
    end if
    call psb_realloc(num_neighbors,prcid,info)

    if (debug) write(*,*) my_rank,' wait'
    pnti   = 1
    snd_pt = total_recv_+1
    rcv_pt = 1
    p2ptag = psb_double_swap_tag
    do i=1, num_neighbors
      proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
      nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
      nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
      if (proc_to_comm /= my_rank)then 
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
      else if (proc_to_comm == my_rank) then 
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

    if (debug) write(*,*) my_rank,' scatter'      
    pnti   = 1
    snd_pt = total_recv_+1
    rcv_pt = 1
    do i=1, num_neighbors
      proc_to_comm = comm_indexes%v(pnti+psb_proc_id_)
      nerv = comm_indexes%v(pnti+psb_n_elem_recv_)
      nesd = comm_indexes%v(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_

      if (debug) write(0,*)my_rank,' Received from: ',prcid(i),&
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
    if (debug) write(*,*) my_rank,' wait'
    call y%device_wait()
    if (debug) write(*,*) my_rank,' free buffer'
    call y%free_buffer(info)
    if (info == 0) then
      if (allocated(y%comm_handle)) call psb_comm_free(y%comm_handle, info)
    end if
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
    if (debug) write(*,*) my_rank,' done'
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
  integer(psb_mpk_)                               :: np, my_rank, n
  integer(psb_mpk_)                               :: iret, p2pstat(mpi_status_size)
  type(psb_comm_neighbor_handle), pointer         :: neighbor_comm_handle
  integer(psb_ipk_)                               :: err_act, topology_total_send, topology_total_recv, buffer_size
  integer(psb_mpk_)                               :: total_send_, total_recv_
  integer(psb_mpk_)                               :: total_send_, total_recv_
  logical                                         :: do_start, do_wait
  logical, parameter                              :: debug = .false.
  character(len=30)                               :: name


  info = psb_success_
  name = 'psi_dswap_neighbor_topology_multivect'
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
    if(debug) write(*,*) my_rank,' nbr_vect: starting data exchange (nonblocking)'
    if (.not. neighbor_comm_handle%is_initialized) then
      if (debug) write(*,*) my_rank,' nbr_vect: building topology via handle'
      call neighbor_comm_handle%topology_init(comm_indexes%v, num_neighbors, total_send, total_recv, &
          & ctxt, icomm, info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_, name, a_err='neighbor_topology_init')
        goto 9999
      end if
    end if
    topology_total_send = neighbor_comm_handle%total_send
    topology_total_recv = neighbor_comm_handle%total_recv
    total_send_ = topology_total_send * n
    total_recv_ = topology_total_recv * n
    total_send_ = topology_total_send * n
    total_recv_ = topology_total_recv * n

    ! Buffer layout:
    !   combuf(1 : total_send)                       = send area
    !   combuf(total_send+1 : total_send+total_recv) = recv area
    buffer_size = total_send_ + total_recv_

    if (buffer_size > 0) then
      call y%new_buffer(buffer_size, info)
      if (info /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_, name)
        goto 9999
      end if
    end if
    neighbor_comm_handle%comm_request = mpi_request_null

    ! Gather send data into contiguous send buffer (polymorphic for GPU)
    if (buffer_size > 0) then
      if (debug) write(*,*) my_rank,' nbr_vect: gathering send data,', topology_total_send,' elems'
      call y%gth(int(topology_total_send,psb_mpk_), &
        & neighbor_comm_handle%send_indexes, &
        & y%combuf(1:total_send_))
    end if

    ! Wait for device (important for GPU subclasses)
    call y%device_wait()

    ! Post non-blocking neighborhood alltoallv
    if (debug) write(*,*) my_rank,' nbr_vect: posting MPI_Ineighbor_alltoallv'
    if (buffer_size > 0) then
      call mpi_ineighbor_alltoallv( &
          & y%combuf(1),                        &  ! send buffer
          & neighbor_comm_handle%send_counts,     &
          & neighbor_comm_handle%send_displs,     &
          & psb_mpi_r_dpk_,                     &
          & y%combuf(total_send_ + 1),            &  ! recv buffer
          & n*neighbor_comm_handle%recv_counts,   &
          & n*neighbor_comm_handle%recv_displs,    &
          & psb_mpi_r_dpk_,                     &
          & neighbor_comm_handle%graph_comm,      &
          & neighbor_comm_handle%comm_request, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if
    else
      neighbor_comm_handle%comm_request = mpi_request_null
    end if

  end if ! do_start

  ! ---------------------------------------------------------
  !  WAIT phase: complete MPI, scatter received data
  ! ---------------------------------------------------------
  if (do_wait) then

    if ((topology_total_send + topology_total_recv) > 0) then
      if (neighbor_comm_handle%comm_request == mpi_request_null) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/-2/))
        goto 9999
      end if
    else
      neighbor_comm_handle%comm_request = mpi_request_null
    end if

    topology_total_send = neighbor_comm_handle%total_send
    topology_total_recv = neighbor_comm_handle%total_recv

    ! Wait for the non-blocking collective to complete
    if ((topology_total_send + topology_total_recv) > 0) then
      if (debug) write(*,*) my_rank,' nbr_vect: waiting on MPI request'
      call mpi_wait(neighbor_comm_handle%comm_request, p2pstat, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if
    end if

    ! Scatter received data to local vector positions (polymorphic for GPU)
    if ((topology_total_send + topology_total_recv) > 0) then
      if (debug) write(*,*) my_rank,' nbr_vect: scattering recv data,', topology_total_recv,' elems'
      call y%sct(int(topology_total_recv,psb_mpk_), &
        & neighbor_comm_handle%recv_indexes, &
        & y%combuf(total_send_+1:total_send_+total_recv_), &
        & beta)
    end if


    ! Clean up
    neighbor_comm_handle%comm_request = mpi_request_null
    call y%device_wait()
    call y%maybe_free_buffer(info)
    if (info /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_, name)
      goto 9999
    end if
    if (debug) write(*,*) my_rank,' nbr_vect: done'

  end if ! do_wait

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psi_dswap_neighbor_topology_multivect



subroutine psi_dswap_neighbor_topology_multivect_persistent(ctxt,swap_status,beta,y,comm_indexes,&
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
  integer(psb_mpk_)                               :: np, my_rank, n
  integer(psb_mpk_)                               :: iret, p2pstat(mpi_status_size)
  type(psb_comm_neighbor_handle), pointer         :: neighbor_comm_handle
  integer(psb_ipk_)                               :: err_act, topology_total_send, topology_total_recv, buffer_size
  integer(psb_mpk_)                               :: total_send_, total_recv_
  logical                                         :: do_start, do_wait
  logical, parameter                              :: debug = .false.
  character(len=30)                               :: name


  info = psb_success_
  name = 'psi_dswap_neighbor_topology_multivect_persistent'
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
    call psb_errpush(info,name,a_err='Expected neighbor comm_handle in persistent neighbor multivect swap')
    goto 9999
  end select

  do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_unknown_)
  do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_unknown_)

  call comm_indexes%sync()

  ! ---------------------------------------------------------
  !  START phase: build topology (if needed), gather, post MPI
  ! ---------------------------------------------------------
  if (do_start) then
    if(debug) write(*,*) my_rank,' nbr_vect: starting data exchange (persistent)'
    if (neighbor_comm_handle%persistent_in_flight) then
      info = psb_err_mpi_error_
      call psb_errpush(info, name, a_err='Invalid START: persistent neighbor request already in flight')
      goto 9999
    end if
    if (.not. neighbor_comm_handle%is_initialized) then
      if (debug) write(*,*) my_rank,' nbr_vect: building topology via handle'
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
      ! Gather send data into contiguous send buffer (polymorphic for GPU)
      if (debug) write(*,*) my_rank,' nbr_vect: gathering send data,', topology_total_send,' elems'
      call y%gth(int(topology_total_send,psb_mpk_), &
        & neighbor_comm_handle%send_indexes, &
        & y%combuf(1:total_send_))
    else
      neighbor_comm_handle%persistent_in_flight = .false.
    end if

    ! Wait for device (important for GPU subclasses)
    call y%device_wait()

    if (.not. neighbor_comm_handle%persistent_request_ready) then
      if (buffer_size > 0) then
        if (debug) write(*,*) my_rank,' nbr_vect: posting MPI_Neighbor_alltoallv_init'
        call mpi_neighbor_alltoallv_init( &
            & y%combuf(1),                          &  ! send buffer
            & n*neighbor_comm_handle%send_counts,   &
            & n*neighbor_comm_handle%send_displs,    &
            & psb_mpi_r_dpk_,                       &
            & y%combuf(total_send_ + 1),            &  ! recv buffer
            & n*neighbor_comm_handle%recv_counts,   &
            & n*neighbor_comm_handle%recv_displs,    &
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

  ! ---------------------------------------------------------
  !  WAIT phase: complete MPI, scatter received data
  ! ---------------------------------------------------------
  if (do_wait) then

    topology_total_send = neighbor_comm_handle%total_send
    topology_total_recv = neighbor_comm_handle%total_recv

    if ((topology_total_send + topology_total_recv) > 0) then
      if (.not. neighbor_comm_handle%persistent_in_flight) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, a_err='Invalid WAIT: no persistent neighbor request in flight')
        goto 9999
      end if

      ! Wait for the persistent collective to complete
      if (debug) write(*,*) my_rank,' nbr_vect: waiting on persistent MPI request'
      call mpi_wait(neighbor_comm_handle%persistent_request, p2pstat, iret)
      if (iret /= mpi_success) then
        info = psb_err_mpi_error_
        call psb_errpush(info, name, m_err=(/iret/))
        goto 9999
      end if
      neighbor_comm_handle%persistent_in_flight = .false.

      ! Scatter received data to local vector positions (polymorphic for GPU)
      if (debug) write(*,*) my_rank,' nbr_vect: scattering recv data,', topology_total_recv,' elems'
      call y%sct(int(topology_total_recv,psb_mpk_), &
        & neighbor_comm_handle%recv_indexes, &
        & y%combuf(total_send_+1:total_send_+total_recv_), &
        & beta)
    else
      neighbor_comm_handle%persistent_in_flight = .false.
    end if

    call y%device_wait()
    if (debug) write(*,*) my_rank,' nbr_vect: done'

  end if ! do_wait

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psi_dswap_neighbor_topology_multivect_persistent


  subroutine psi_dswap_rma_pull_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)               :: ctxt
    integer(psb_ipk_), intent(in)                 :: swap_status
    real(psb_dpk_), intent(in)                  :: beta
    class(psb_d_base_multivect_type), intent(inout), target :: y
    class(psb_i_base_vect_type), intent(inout)    :: comm_indexes
    integer(psb_ipk_), intent(in)                 :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)    :: comm_handle
    integer(psb_ipk_), intent(out)                :: info

    integer(psb_mpk_) :: np, my_rank, iret, element_bytes, icomm, n
    integer(psb_mpk_) :: proc_to_comm, prc_rank, recv_count, send_count, send_pos, recv_pos, list_pos
    integer(psb_mpk_) :: remote_base
    integer(kind=MPI_ADDRESS_KIND) :: remote_disp, exposed_bytes
    integer(psb_ipk_) :: err_act, neighbor_idx, buffer_size, total_send_, total_recv_
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    integer(psb_mpk_) :: comm_grp, nbr_grp, n_off, grp_idx
    integer(psb_mpk_), allocatable :: grp_ranks(:)
    logical :: do_start, do_wait, memory_buffer_layout_rebuild_needed
    type(psb_comm_rma_handle), pointer :: rma_handle
    character(len=30) :: name

    info = psb_success_
    name = 'psi_dswap_rma_pull_multivect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    icomm = ctxt%get_mpic()
    n = y%get_ncols()
    total_send_ = total_send * n
    total_recv_ = total_recv * n

    select type(ch => comm_handle)
    type is(psb_comm_rma_handle)
      rma_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected RMA comm_handle for pull mode')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    if (do_start) then
      buffer_size = total_send_ + total_recv_
      memory_buffer_layout_rebuild_needed = (.not. rma_handle%layout_ready) .or. &
         & (rma_handle%layout_nnbr /= num_neighbors) .or. &
         & (rma_handle%layout_send /= total_send) .or. &
         & (rma_handle%layout_recv /= total_recv)

      if (memory_buffer_layout_rebuild_needed) then
        if (allocated(peer_mpi_rank)) deallocate(peer_mpi_rank)
        if (num_neighbors > 0) then
          allocate(peer_mpi_rank(num_neighbors), stat=iret)
          if (iret /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info,name,a_err='RMA pull rank cache allocation')
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
        call rma_handle%init_memory_buffer_layout(info, comm_indexes%v, peer_mpi_rank, &
             & num_neighbors, total_send, total_recv, my_rank, icomm)
        if (info /= psb_success_) then
          call psb_errpush(info,name,a_err='RMA pull init_memory_buffer_layout failure')
          goto 9999
        end if
      end if

      ! Expose y%combuf as the RMA window (created once, reused). On GPU vectors
      ! combuf is the pinned/registered buffer, so RMA uses pinned memory.
      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          call y%new_buffer(buffer_size, info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name); goto 9999
          end if
        else if (size(y%combuf) < buffer_size) then
          if (rma_handle%window_open) then
            call mpi_win_unlock_all(rma_handle%win, iret)
            rma_handle%window_open = .false.
          end if
          ! Dynamic window: detach the old buffer and keep the window. Detach
          ! is a local call, so growing the buffer no longer costs a collective.
          if (rma_handle%buf_attached) then
            call mpi_win_detach(rma_handle%win, y%combuf, iret)
            rma_handle%buf_attached = .false.
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name); goto 9999
          end if
        end if

        element_bytes = storage_size(y%combuf(1))/8
        if (.not. rma_handle%window_ready) then
          ! The window comes from the module pool, keyed by communicator: it is
          ! created on the first request of the run and reused by every handle
          ! afterwards. Creating it per handle, as before, meant paying a
          ! communicator-wide collective (~0.9 s at 448 ranks) each time.
          call psb_comm_rma_get_window(ctxt%get_mpic(), rma_handle%win, info)
          if (info /= psb_success_) then
            call psb_errpush(info,name); goto 9999
          end if
          rma_handle%window_ready = .true.
        end if

        if (.not. rma_handle%buf_attached) then
          exposed_bytes = int(size(y%combuf),kind=MPI_ADDRESS_KIND) * int(element_bytes,kind=MPI_ADDRESS_KIND)
          call mpi_win_attach(rma_handle%win, y%combuf, exposed_bytes, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
          call mpi_get_address(y%combuf, rma_handle%my_buf_addr, iret)
          ! Kept so that the handle can detach on free, when the buffer is no
          ! longer reachable from here.
          rma_handle%win_base  = c_loc(y%combuf(1))
          rma_handle%win_nelem = size(y%combuf)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
          ! Targets are addressed absolutely on a dynamic window, so neighbours
          ! must be told the new base before any Get/Put can reference it.
          call rma_handle%publish_buf_addr(ctxt%get_mpic(), info)
          if (info /= psb_success_) then
            call psb_errpush(info,name); goto 9999
          end if
          rma_handle%buf_attached = .true.
        end if
      end if

      if (buffer_size > 0) then
        if (total_send > 0) then
          call y%gth(int(total_send,psb_mpk_), rma_handle%peer_send_indexes, y%combuf(1:total_send_))
        end if
        call y%device_wait()

        ! Open exposure (post) + access (start) epochs over the CACHED neighbor
        ! group (built once with the layout in the handle) -- no per-call build.
        if (rma_handle%nbr_grp /= mpi_group_null) then
          call mpi_win_post(rma_handle%nbr_grp, 0, rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
          call mpi_win_start(rma_handle%nbr_grp, 0, rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
        end if

        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          recv_count = rma_handle%peer_recv_counts(neighbor_idx)
          send_count = rma_handle%peer_send_counts(neighbor_idx)
          prc_rank = rma_handle%peer_mpi_rank(neighbor_idx)
          send_pos = rma_handle%peer_send_displs(neighbor_idx)*n + 1
          recv_pos = total_send_ + rma_handle%peer_recv_displs(neighbor_idx)*n + 1

          if (proc_to_comm /= my_rank) then
            remote_base = rma_handle%peer_remote_send_displs(neighbor_idx)
            if (remote_base < 1) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='Invalid remote metadata in RMA pull')
              goto 9999
            end if
            if (recv_count > 0) then
              ! Dynamic window: absolute target address, offset in bytes.
              remote_disp = rma_handle%peer_buf_addr(neighbor_idx) &
                   & + int((remote_base - 1) * n, kind=MPI_ADDRESS_KIND) &
                   &   * int(element_bytes, kind=MPI_ADDRESS_KIND)
              call mpi_get(y%combuf(recv_pos), recv_count*n, psb_mpi_r_dpk_, prc_rank, remote_disp, recv_count*n, psb_mpi_r_dpk_, &
                   & rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
            end if
          else
            if (send_count /= recv_count) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA pull self-copy mismatch')
              goto 9999
            end if
            y%combuf(recv_pos:recv_pos+recv_count*n-1) = y%combuf(send_pos:send_pos+send_count*n-1)
          end if
        end do

        ! Close access epoch (my GETs done) then exposure epoch (peers done
        ! reading my window): received data is now in y%combuf. Cached group, no free.
        if (rma_handle%nbr_grp /= mpi_group_null) then
          call mpi_win_complete(rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
          call mpi_win_wait(rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
        end if
      end if
    end if

    ! WAIT phase: GETs already complete (per-neighbor unlock in START); scatter received data into Y.
    if (do_wait) then
      if (total_recv > 0) then
        call y%sct(int(total_recv,psb_mpk_), rma_handle%peer_recv_indexes, y%combuf(total_send_+1:total_send_+total_recv_), beta)
      end if
      call y%device_wait()
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_dswap_rma_pull_multivect


  subroutine psi_dswap_rma_push_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,comm_handle,info)
#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)               :: ctxt
    integer(psb_ipk_), intent(in)                 :: swap_status
    real(psb_dpk_), intent(in)                  :: beta
    class(psb_d_base_multivect_type), intent(inout), target :: y
    class(psb_i_base_vect_type), intent(inout)    :: comm_indexes
    integer(psb_ipk_), intent(in)                 :: num_neighbors, total_send, total_recv
    class(psb_comm_handle_type), intent(inout)    :: comm_handle
    integer(psb_ipk_), intent(out)                :: info

    integer(psb_mpk_) :: np, my_rank, iret, element_bytes, icomm, n
    integer(psb_mpk_) :: proc_to_comm, prc_rank, recv_count, send_count, send_pos, recv_pos, list_pos
    integer(psb_mpk_) :: remote_base
    integer(kind=MPI_ADDRESS_KIND) :: remote_disp, exposed_bytes
    integer(psb_ipk_) :: err_act, neighbor_idx, buffer_size, total_send_, total_recv_
    integer(psb_ipk_), allocatable :: peer_mpi_rank(:)
    integer(psb_mpk_) :: comm_grp, nbr_grp, n_off, grp_idx
    integer(psb_mpk_), allocatable :: grp_ranks(:)
    logical :: do_start, do_wait, memory_buffer_layout_rebuild_needed
    type(psb_comm_rma_handle), pointer :: rma_handle
    character(len=30) :: name

    info = psb_success_
    name = 'psi_dswap_rma_push_multivect'
    call psb_erractionsave(err_act)
    call psb_info(ctxt,my_rank,np)
    if (np == -1) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    icomm = ctxt%get_mpic()
    n = y%get_ncols()
    total_send_ = total_send * n
    total_recv_ = total_recv * n

    select type(ch => comm_handle)
    type is(psb_comm_rma_handle)
      rma_handle => ch
    class default
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Expected RMA comm_handle for push mode')
      goto 9999
    end select

    do_start = (swap_status == psb_comm_status_start_) .or. (swap_status == psb_comm_status_sync_)
    do_wait  = (swap_status == psb_comm_status_wait_)  .or. (swap_status == psb_comm_status_sync_)

    call comm_indexes%sync()

    if (do_start) then
      buffer_size = total_send_ + total_recv_
      memory_buffer_layout_rebuild_needed = (.not. rma_handle%layout_ready) .or. &
         & (rma_handle%layout_nnbr /= num_neighbors) .or. &
         & (rma_handle%layout_send /= total_send) .or. &
         & (rma_handle%layout_recv /= total_recv)

      if (memory_buffer_layout_rebuild_needed) then
        if (allocated(peer_mpi_rank)) deallocate(peer_mpi_rank)
        if (num_neighbors > 0) then
          allocate(peer_mpi_rank(num_neighbors), stat=iret)
          if (iret /= 0) then
            info = psb_err_alloc_dealloc_
            call psb_errpush(info,name,a_err='RMA put rank cache allocation')
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
        call rma_handle%init_memory_buffer_layout(info, comm_indexes%v, peer_mpi_rank, &
             & num_neighbors, total_send, total_recv, my_rank, icomm)
        if (info /= psb_success_) then
          call psb_errpush(info,name,a_err='RMA put ini_memory_buffer_layout')
          goto 9999
        end if
      end if

      ! Expose y%combuf as the RMA window (created once, reused). On GPU vectors
      ! combuf is the pinned/registered buffer, so RMA uses pinned memory.
      if (buffer_size > 0) then
        if (.not. allocated(y%combuf)) then
          call y%new_buffer(buffer_size, info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name); goto 9999
          end if
        else if (size(y%combuf) < buffer_size) then
          if (rma_handle%window_open) then
            call mpi_win_unlock_all(rma_handle%win, iret)
            rma_handle%window_open = .false.
          end if
          ! Dynamic window: detach the old buffer and keep the window. Detach
          ! is a local call, so growing the buffer no longer costs a collective.
          if (rma_handle%buf_attached) then
            call mpi_win_detach(rma_handle%win, y%combuf, iret)
            rma_handle%buf_attached = .false.
          end if
          call y%new_buffer(buffer_size, info)
          if (info /= psb_success_) then
            call psb_errpush(psb_err_alloc_dealloc_,name); goto 9999
          end if
        end if

        element_bytes = storage_size(y%combuf(1))/8
        if (.not. rma_handle%window_ready) then
          ! The window comes from the module pool, keyed by communicator: it is
          ! created on the first request of the run and reused by every handle
          ! afterwards. Creating it per handle, as before, meant paying a
          ! communicator-wide collective (~0.9 s at 448 ranks) each time.
          call psb_comm_rma_get_window(ctxt%get_mpic(), rma_handle%win, info)
          if (info /= psb_success_) then
            call psb_errpush(info,name); goto 9999
          end if
          rma_handle%window_ready = .true.
        end if

        if (.not. rma_handle%buf_attached) then
          exposed_bytes = int(size(y%combuf),kind=MPI_ADDRESS_KIND) * int(element_bytes,kind=MPI_ADDRESS_KIND)
          call mpi_win_attach(rma_handle%win, y%combuf, exposed_bytes, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
          call mpi_get_address(y%combuf, rma_handle%my_buf_addr, iret)
          ! Kept so that the handle can detach on free, when the buffer is no
          ! longer reachable from here.
          rma_handle%win_base  = c_loc(y%combuf(1))
          rma_handle%win_nelem = size(y%combuf)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
          ! Targets are addressed absolutely on a dynamic window, so neighbours
          ! must be told the new base before any Get/Put can reference it.
          call rma_handle%publish_buf_addr(ctxt%get_mpic(), info)
          if (info /= psb_success_) then
            call psb_errpush(info,name); goto 9999
          end if
          rma_handle%buf_attached = .true.
        end if
      end if

      if (buffer_size > 0) then
        if (total_send > 0) then
          call y%gth(int(total_send,psb_mpk_), rma_handle%peer_send_indexes, y%combuf(1:total_send_))
        end if
        call y%device_wait()

        ! Open exposure (post) + access (start) epochs over the CACHED neighbor
        ! group (built once with the layout in the handle) -- no per-call build.
        if (rma_handle%nbr_grp /= mpi_group_null) then
          call mpi_win_post(rma_handle%nbr_grp, 0, rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
          call mpi_win_start(rma_handle%nbr_grp, 0, rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
        end if

        ! Issue every PUT (remote) or perform the self-copy (local).
        do neighbor_idx=1, num_neighbors
          proc_to_comm = rma_handle%peer_proc(neighbor_idx)
          recv_count = rma_handle%peer_recv_counts(neighbor_idx)
          send_count = rma_handle%peer_send_counts(neighbor_idx)
          prc_rank = rma_handle%peer_mpi_rank(neighbor_idx)
          send_pos = rma_handle%peer_send_displs(neighbor_idx)*n + 1
          recv_pos = total_send_ + rma_handle%peer_recv_displs(neighbor_idx)*n + 1

          if (proc_to_comm /= my_rank) then
            remote_base = rma_handle%peer_remote_recv_displs(neighbor_idx)
            if (remote_base < 1) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='Invalid remote metadata in RMA push')
              goto 9999
            end if
            if (send_count > 0) then
              ! Dynamic window: absolute target address, offset in bytes.
              remote_disp = rma_handle%peer_buf_addr(neighbor_idx) &
                   & + int((remote_base - 1) * n, kind=MPI_ADDRESS_KIND) &
                   &   * int(element_bytes, kind=MPI_ADDRESS_KIND)
              call mpi_put(y%combuf(send_pos), send_count*n, psb_mpi_r_dpk_, prc_rank, remote_disp, send_count*n, psb_mpi_r_dpk_, &
                   & rma_handle%win, iret)
              if (iret /= mpi_success) then
                info = psb_err_mpi_error_
                call psb_errpush(info,name,m_err=(/iret/))
                goto 9999
              end if
            end if
          else
            if (send_count /= recv_count) then
              info = psb_err_internal_error_
              call psb_errpush(info,name,a_err='RMA push self-copy mismatch')
              goto 9999
            end if
            y%combuf(recv_pos:recv_pos+recv_count*n-1) = y%combuf(send_pos:send_pos+send_count*n-1)
          end if
        end do

        ! Close access epoch (my PUTs done) then exposure epoch (peers' PUTs into
        ! me done): receive area now holds the pushed data. Cached group, no free.
        if (rma_handle%nbr_grp /= mpi_group_null) then
          call mpi_win_complete(rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
          call mpi_win_wait(rma_handle%win, iret)
          if (iret /= mpi_success) then
            info = psb_err_mpi_error_; call psb_errpush(info,name,m_err=(/iret/)); goto 9999
          end if
        end if
      end if
    end if

    ! WAIT phase: data already delivered by PSCW in START; just scatter into Y.
    if (do_wait) then
      if (total_recv > 0) then
        call y%sct(int(total_recv,psb_mpk_), rma_handle%peer_recv_indexes, y%combuf(total_send_+1:total_send_+total_recv_), beta)
      end if
      call y%device_wait()
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return
  end subroutine psi_dswap_rma_push_multivect


end submodule psi_d_swapdata_impl

