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
! File: psi_dswaptran.F90
!
! Subroutine: psi_dswaptran_vect
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
submodule (psi_d_comm_v_mod)  psi_d_swaptran_impl
  use psb_base_mod
  use psb_comm_schemes_mod, only: psb_comm_handle_type, psb_comm_isend_irecv_
  use psb_comm_factory_mod, only: psb_comm_init, psb_comm_free
  use psb_comm_baseline_mod, only: psb_comm_baseline_handle
  use psb_comm_neighbor_impl_mod, only: psb_comm_neighbor_handle
contains
  module subroutine psi_dswaptran_vect(swap_status,beta,y,desc_a,info,data)

#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    integer(psb_ipk_), intent(in)               :: swap_status
    real(psb_dpk_), intent(in)                  :: beta
    class(psb_d_base_vect_type), intent(inout)  :: y
    type(psb_desc_type),target                  :: desc_a
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_), optional                 :: data

    ! locals
    type(psb_ctxt_type)                         :: ctxt
    integer(psb_mpk_)                           :: icomm
    integer(psb_ipk_)                           :: np, me, total_send, total_recv, num_neighbors, err_act, data_
    class(psb_i_base_vect_type), pointer        :: comm_indexes
    character(len=20)  :: name

    ! local variables used to detect the communication scheme
    logical                                     :: swap_mpi, swap_sync, swap_send, swap_recv, swap_start, swap_wait
    logical                                     :: baseline, neighbor_a2av


    info = psb_success_
    name = 'psi_dswaptran_vect'
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

    swap_mpi    = iand(swap_status,psb_swap_mpi_) /= 0
    swap_sync   = iand(swap_status,psb_swap_sync_) /= 0
    swap_send   = iand(swap_status,psb_swap_send_) /= 0
    swap_recv   = iand(swap_status,psb_swap_recv_) /= 0
    swap_start  = iand(swap_status,psb_swap_start_) /= 0
    swap_wait   = iand(swap_status,psb_swap_wait_) /= 0

    baseline = swap_mpi .or. swap_send .or. swap_recv .or. swap_sync
    neighbor_a2av = swap_start .or. swap_wait 

    if( (baseline.eqv..true.).and.(neighbor_a2av.eqv..true.) ) then 
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Incompatible swap_status settings: both baseline and neighbor_a2av are true')
      goto 9999
    end if

    if (.not. allocated(y%comm_handle)) then
      call psb_comm_init(psb_comm_isend_irecv_, y%comm_handle, info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_, name, a_err='init comm default baseline')
        goto 9999
      end if
    end if

    if (baseline) then 
      call psi_dtran_baseline_vect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='baseline swap')
        goto 9999
      end if
    else if (neighbor_a2av) then 
      call psi_dtran_neighbor_topology_vect(ctxt,swap_status,beta,y,comm_indexes,&
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
  end subroutine psi_dswaptran_vect

  !
  !
  ! Subroutine: psi_dtran_baseline_vect
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
  module subroutine psi_dtran_baseline_vect(ctxt,swap_status,beta,y,idx,&
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
    name = 'psi_dtran_baseline_vect'
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
    swap_mpi  = iand(swap_status,psb_swap_mpi_) /= 0
    swap_sync = iand(swap_status,psb_swap_sync_) /= 0
    swap_send = iand(swap_status,psb_swap_send_) /= 0
    swap_recv = iand(swap_status,psb_swap_recv_) /= 0
    do_send = swap_mpi .or. swap_sync .or. swap_send
    do_recv = swap_mpi .or. swap_sync .or. swap_recv

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

  end subroutine psi_dtran_baseline_vect




  subroutine psi_dtran_neighbor_topology_vect(ctxt,swap_status,beta,y,comm_indexes,&
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
    name = 'psi_dtran_neighbor_topology_vect'
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

    do_start = iand(swap_status,psb_swap_start_) /= 0
    do_wait  = iand(swap_status,psb_swap_wait_)  /= 0

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
          & psb_mpi_r_dpk_,                        &
          & y%combuf(topology_total_recv + 1),     &  ! recv buffer (will contain send_indexes data)
          & neighbor_comm_handle%send_counts,       &
          & neighbor_comm_handle%send_displs,       &
          & psb_mpi_r_dpk_,                        &
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
  end subroutine psi_dtran_neighbor_topology_vect




  !
  !
  !
  !
  ! Subroutine: psi_dswaptran_multivect
  !   Data exchange among processes.
  !
  !   Takes care of Y an encaspulated multivector.
  !   
  !   
  module subroutine psi_dswaptran_multivect(swap_status,beta,y,desc_a,info,data)

#ifdef PSB_MPI_MOD
    use mpi
#endif
    implicit none
#ifdef PSB_MPI_H
    include 'mpif.h'
#endif

    integer(psb_ipk_), intent(in)                   :: swap_status
    real(psb_dpk_), intent(in)                      :: beta
    class(psb_d_base_multivect_type), intent(inout) :: y
    type(psb_desc_type),target                      :: desc_a
    integer(psb_ipk_), intent(out)                  :: info
    integer(psb_ipk_), optional                     :: data

    ! locals
    type(psb_ctxt_type)                             :: ctxt
    integer(psb_mpk_)                               :: icomm
    integer(psb_ipk_)                               :: np, me, total_send, total_recv, num_neighbors, err_act, data_
    class(psb_i_base_vect_type), pointer            :: comm_indexes
    character(len=20)                               :: name

    ! local variables used to detect the communication scheme
    logical                                         :: swap_mpi, swap_sync, swap_send, swap_recv, swap_start, swap_wait
    logical                                         :: baseline, neighbor_a2av


    info = psb_success_
    name = 'psi_dswaptran_multivect'
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

    swap_mpi    = iand(swap_status,psb_swap_mpi_) /= 0
    swap_sync   = iand(swap_status,psb_swap_sync_) /= 0
    swap_send   = iand(swap_status,psb_swap_send_) /= 0
    swap_recv   = iand(swap_status,psb_swap_recv_) /= 0
    swap_start  = iand(swap_status,psb_swap_start_) /= 0
    swap_wait   = iand(swap_status,psb_swap_wait_) /= 0

    baseline = swap_mpi .or. swap_send .or. swap_recv .or. swap_sync
    neighbor_a2av = swap_start .or. swap_wait 

    if( (baseline.eqv..true.).and.(neighbor_a2av.eqv..true.) ) then 
      info = psb_err_mpi_error_
      call psb_errpush(info,name,a_err='Incompatible swap_status settings: both baseline and neighbor_a2av are true')
      goto 9999
    end if

    if (.not. allocated(y%comm_handle)) then
      call psb_comm_init(psb_comm_isend_irecv_, y%comm_handle, info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_, name, a_err='init comm default baseline')
        goto 9999
      end if
    end if

    if (baseline) then 
      call psi_dtran_baseline_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,total_send,total_recv,y%comm_handle,info)
      if (info /= psb_success_) then 
        call psb_errpush(info,name,a_err='baseline swap')
        goto 9999
      end if
    else if (neighbor_a2av) then 
      call psi_dtran_neighbor_topology_multivect(ctxt,swap_status,beta,y,comm_indexes,num_neighbors,&
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
  end subroutine psi_dswaptran_multivect

  subroutine psi_dtran_baseline_multivect(ctxt,swap_status,beta,y,idx,&
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
    class(psb_d_base_multivect_type), intent(inout) :: y
    real(psb_dpk_), intent(in)                      :: beta
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
    name = 'psi_dtran_vidx_multivect'
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

    swap_mpi  = iand(swap_status,psb_swap_mpi_) /= 0
    swap_sync = iand(swap_status,psb_swap_sync_) /= 0
    swap_send = iand(swap_status,psb_swap_send_) /= 0
    swap_recv = iand(swap_status,psb_swap_recv_) /= 0
    do_send = swap_mpi .or. swap_sync .or. swap_send
    do_recv = swap_mpi .or. swap_sync .or. swap_recv

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

  end subroutine psi_dtran_baseline_multivect


  subroutine psi_dtran_neighbor_topology_multivect(ctxt,swap_status,beta,y,comm_indexes,&
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
    name = 'psi_dtran_neighbor_topology_multivect'
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

    do_start = iand(swap_status,psb_swap_start_) /= 0
    do_wait  = iand(swap_status,psb_swap_wait_)  /= 0

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
          & psb_mpi_r_dpk_,                        &
          & y%combuf(topology_total_recv + 1),     &  ! recv buffer (will contain send_indexes data)
          & neighbor_comm_handle%send_counts,       &
          & neighbor_comm_handle%send_displs,       &
          & psb_mpi_r_dpk_,                        &
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
  end subroutine psi_dtran_neighbor_topology_multivect



end submodule psi_d_swaptran_impl
