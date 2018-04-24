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
! File: psi_cswapdata.F90
!
! Subroutine: psi_cswapdatam
!   Does the data exchange among processes. Essentially this is doing 
!   a variable all-to-all data exchange (ALLTOALLV in MPI parlance), but 
!   it is capable of pruning empty exchanges, which are very likely in out 
!   application environment. All the variants have the same structure 
!   In all these subroutines X may be:    I    Integer
!                                         S    real(psb_spk_)
!                                         D    real(psb_dpk_)
!                                         C    complex(psb_spk_)
!                                         Z    complex(psb_dpk_)
!   Basically the operation is as follows: on each process, we identify 
!   sections SND(Y) and RCV(Y); then we do a send on (PACK(SND(Y)));
!   then we receive, and we do an update with Y = UNPACK(RCV(Y)) + BETA * Y 
!   but only on the elements involved in the UNPACK operation. 
!   Thus: for halo data exchange, the receive section is confined in the 
!   halo indices, and BETA=0, whereas for overlap exchange the receive section 
!   is scattered in the owned indices, and BETA=1.
! 
! Arguments: 
!    flag     - integer                 Choose the algorithm for data exchange: 
!                                       this is chosen through bit fields. 
!                                        swap_mpi  = iand(flag,psb_swap_mpi_)  /= 0
!                                        swap_sync = iand(flag,psb_swap_sync_) /= 0
!                                        swap_send = iand(flag,psb_swap_send_) /= 0
!                                        swap_recv = iand(flag,psb_swap_recv_) /= 0
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
!    beta     - X                       Choose overwrite or sum. 
!    y(:,:)   - X                       The data area                        
!    desc_a   - type(psb_desc_type).  The communication descriptor.        
!    work(:)  - X                       Buffer space. If not sufficient, will do 
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
!
!
! Subroutine: psi_cswapdata_vect
!   Data exchange among processes.
!
!   Takes care of Y an exanspulated vector.
!   
!   
! 
subroutine psi_cswapdata_vect(flag,beta,y,desc_a,work,info,data)

  use psi_mod, psb_protect_name => psi_cswapdata_vect
  use psb_c_base_vect_mod
  use psb_error_mod
  use psb_desc_mod
  use psb_penv_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer(psb_ipk_), intent(in)         :: flag
  integer(psb_ipk_), intent(out)        :: info
  class(psb_c_base_vect_type) :: y
  complex(psb_spk_)           :: beta
  complex(psb_spk_), target   :: work(:)
  type(psb_desc_type), target  :: desc_a
  integer(psb_ipk_), optional           :: data

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, icomm, idxs, idxr, totxch, data_, err_act
  class(psb_i_base_vect_type), pointer :: d_vidx
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_datav'
  call psb_erractionsave(err_act)

  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  call psb_info(ictxt,me,np) 
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

  call desc_a%get_list(data_,d_vidx,totxch,idxr,idxs,info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,a_err='psb_cd_get_list')
    goto 9999
  end if

  call psi_swapdata(ictxt,icomm,flag,beta,y,d_vidx,totxch,idxs,idxr,work,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

    return
end subroutine psi_cswapdata_vect


!
!
! Subroutine: psi_cswap_vidx_vect
!   Data exchange among processes.
!
!   Takes care of Y an exanspulated vector. Relies on the gather/scatter methods
!   of vectors. 
!   
!   The real workhorse: the outer routine will only choose the index list
!   this one takes the index list and does the actual exchange. 
!   
!   
! 
subroutine psi_cswap_vidx_vect(iictxt,iicomm,flag,beta,y,idx, &
     & totxch,totsnd,totrcv,work,info)

  use psi_mod, psb_protect_name => psi_cswap_vidx_vect
  use psb_error_mod
  use psb_realloc_mod
  use psb_desc_mod
  use psb_penv_mod
  use psb_c_base_vect_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer(psb_ipk_), intent(in)         :: iictxt,iicomm,flag
  integer(psb_ipk_), intent(out)        :: info
  class(psb_c_base_vect_type) :: y
  complex(psb_spk_)           :: beta
  complex(psb_spk_), target   :: work(:)
  class(psb_i_base_vect_type), intent(inout) :: idx
  integer(psb_ipk_), intent(in)              :: totxch,totsnd, totrcv

  ! locals
  integer(psb_mpk_) :: ictxt, icomm, np, me,&
       & proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret
  integer(psb_mpk_), allocatable :: prcid(:)
  integer(psb_ipk_) :: nesd, nerv,&
       & err_act, i, idx_pt, totsnd_, totrcv_,&
       & snd_pt, rcv_pt, pnti, n
  logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
       & albf,do_send,do_recv
  logical, parameter :: usersend=.false., debug=.false.
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_datav'
  call psb_erractionsave(err_act)
  ictxt = iictxt
  icomm = iicomm

  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  n=1
  swap_mpi  = iand(flag,psb_swap_mpi_) /= 0
  swap_sync = iand(flag,psb_swap_sync_) /= 0
  swap_send = iand(flag,psb_swap_send_) /= 0
  swap_recv = iand(flag,psb_swap_recv_) /= 0
  do_send = swap_mpi .or. swap_sync .or. swap_send
  do_recv = swap_mpi .or. swap_sync .or. swap_recv

  totrcv_ = totrcv * n
  totsnd_ = totsnd * n
  call idx%sync()

  if (debug) write(*,*) me,'Internal buffer'
  if (do_send) then 
    if (allocated(y%comid)) then
      if (any(y%comid /= mpi_request_null)) then 
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
    call y%new_comid(totxch,info)
    y%comid = mpi_request_null
    call psb_realloc(totxch,prcid,info)
    ! First I post all the non blocking receives
    pnti   = 1
    do i=1, totxch
      proc_to_comm = idx%v(pnti+psb_proc_id_)
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)

      rcv_pt = 1+pnti+psb_n_elem_recv_
      call psb_get_rank(prcid(i),ictxt,proc_to_comm)      
      if ((nerv>0).and.(proc_to_comm /= me)) then 
        if (debug) write(*,*) me,'Posting receive from',prcid(i),rcv_pt
        p2ptag = psb_complex_swap_tag
        call mpi_irecv(y%combuf(rcv_pt),nerv,&
             & psb_mpi_c_spk_,prcid(i),&
             & p2ptag, icomm,y%comid(i,2),iret)
      end if
      pnti   = pnti + nerv + nesd + 3
    end do
    if (debug) write(*,*) me,' Gather '
    !
    ! Then gather for sending.
    !    
    pnti   = 1
    do i=1, totxch
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)
      snd_pt = 1+pnti+nerv+psb_n_elem_send_
      rcv_pt = 1+pnti+psb_n_elem_recv_
      idx_pt = snd_pt
      call y%gth(idx_pt,nesd,idx)
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
    p2ptag = psb_complex_swap_tag
    do i=1, totxch
      proc_to_comm = idx%v(pnti+psb_proc_id_)
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)
      snd_pt = 1+pnti+nerv+psb_n_elem_send_
      rcv_pt = 1+pnti+psb_n_elem_recv_

      if ((nesd>0).and.(proc_to_comm /= me)) then 
        call mpi_isend(y%combuf(snd_pt),nesd,&
             & psb_mpi_c_spk_,prcid(i),&
             & p2ptag,icomm,y%comid(i,1),iret)
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
    if (.not.allocated(y%comid)) then 
      ! 
      ! No matching send? Something is wrong....
      !
      info=psb_err_mpi_error_
      call psb_errpush(info,name,m_err=(/-2/))
      goto 9999
    end if
    call psb_realloc(totxch,prcid,info)

    if (debug) write(*,*) me,' wait'
    pnti   = 1
    p2ptag = psb_complex_swap_tag
    do i=1, totxch
      proc_to_comm = idx%v(pnti+psb_proc_id_)
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)
      snd_pt = 1+pnti+nerv+psb_n_elem_send_
      rcv_pt = 1+pnti+psb_n_elem_recv_

      if (proc_to_comm /= me)then 
        if (nesd>0) then 
          call mpi_wait(y%comid(i,1),p2pstat,iret)
          if(iret /= mpi_success) then
            info=psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
        end if
        if (nerv>0) then 
          call mpi_wait(y%comid(i,2),p2pstat,iret)
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
    do i=1, totxch
      proc_to_comm = idx%v(pnti+psb_proc_id_)
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_
      snd_pt = 1+pnti+nerv+psb_n_elem_send_
      rcv_pt = 1+pnti+psb_n_elem_recv_

      if (debug) write(0,*)me,' Received from: ',prcid(i),&
           & y%combuf(rcv_pt:rcv_pt+nerv-1)        
      call y%sct(rcv_pt,nerv,idx,beta)
      pnti   = pnti + nerv + nesd + 3
    end do
    !
    ! Waited for everybody, clean up
    !
    y%comid = mpi_request_null

    !
    ! Then wait for device
    !
    if (debug) write(*,*) me,' wait'
    call y%device_wait()
    if (debug) write(*,*) me,' free buffer'
    call y%maybe_free_buffer(info)
    if (info == 0) call y%free_comid(info)
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
    if (debug) write(*,*) me,' done'
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(iictxt,err_act)

  return
end subroutine psi_cswap_vidx_vect

!
!
! Subroutine: psi_cswapdata_multivect
!   Data exchange among processes.
!
!   Takes care of Y an encaspulated vector.
!   
!   
subroutine psi_cswapdata_multivect(flag,beta,y,desc_a,work,info,data)

  use psi_mod, psb_protect_name => psi_cswapdata_multivect
  use psb_c_base_multivect_mod
  use psb_error_mod
  use psb_desc_mod
  use psb_penv_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer(psb_ipk_), intent(in)         :: flag
  integer(psb_ipk_), intent(out)        :: info
  class(psb_c_base_multivect_type) :: y
  complex(psb_spk_)           :: beta
  complex(psb_spk_), target   :: work(:)
  type(psb_desc_type), target  :: desc_a
  integer(psb_ipk_), optional           :: data

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, icomm, idxs, idxr, totxch, data_, err_act
  class(psb_i_base_vect_type), pointer :: d_vidx
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_datav'
  call psb_erractionsave(err_act)

  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  call psb_info(ictxt,me,np) 
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

  call desc_a%get_list(data_,d_vidx,totxch,idxr,idxs,info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,a_err='psb_cd_get_list')
    goto 9999
  end if

  call psi_swapdata(ictxt,icomm,flag,beta,y,d_vidx,totxch,idxs,idxr,work,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

    return
end subroutine psi_cswapdata_multivect


!
!
! Subroutine: psi_cswap_vidx_multivect
!   Data exchange among processes.
!
!   Takes care of Y an encapsulated multivector. Relies on the gather/scatter methods
!   of multivectors. 
!   
!   The real workhorse: the outer routine will only choose the index list
!   this one takes the index list and does the actual exchange. 
!   
!   
! 
subroutine psi_cswap_vidx_multivect(iictxt,iicomm,flag,beta,y,idx, &
     & totxch,totsnd,totrcv,work,info)

  use psi_mod, psb_protect_name => psi_cswap_vidx_multivect
  use psb_error_mod
  use psb_realloc_mod
  use psb_desc_mod
  use psb_penv_mod
  use psb_c_base_multivect_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer(psb_ipk_), intent(in)         :: iictxt,iicomm,flag
  integer(psb_ipk_), intent(out)        :: info
  class(psb_c_base_multivect_type) :: y
  complex(psb_spk_)         :: beta
  complex(psb_spk_), target :: work(:)
  class(psb_i_base_vect_type), intent(inout) :: idx
  integer(psb_ipk_), intent(in)              :: totxch,totsnd, totrcv

  ! locals
  integer(psb_mpk_) :: ictxt, icomm, np, me,&
       & proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret
  integer(psb_mpk_), allocatable :: prcid(:)
  integer(psb_ipk_) :: nesd, nerv,&
       & err_act, i, idx_pt, totsnd_, totrcv_,&
       & snd_pt, rcv_pt, pnti, n
  logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
       & albf,do_send,do_recv
  logical, parameter :: usersend=.false., debug=.false.
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_datav'
  call psb_erractionsave(err_act)
  ictxt = iictxt
  icomm = iicomm

  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  n = y%get_ncols()

  swap_mpi  = iand(flag,psb_swap_mpi_)  /= 0
  swap_sync = iand(flag,psb_swap_sync_) /= 0
  swap_send = iand(flag,psb_swap_send_) /= 0
  swap_recv = iand(flag,psb_swap_recv_) /= 0
  do_send = swap_mpi .or. swap_sync .or. swap_send
  do_recv = swap_mpi .or. swap_sync .or. swap_recv

  totrcv_ = totrcv * n
  totsnd_ = totsnd * n

  call idx%sync()

  if (debug) write(*,*) me,'Internal buffer'
  if (do_send) then 
    if (allocated(y%comid)) then 
      if (any(y%comid /= mpi_request_null)) then 
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
    call y%new_comid(totxch,info)
    y%comid = mpi_request_null
    call psb_realloc(totxch,prcid,info)
    ! First I post all the non blocking receives
    pnti   = 1
    snd_pt = totrcv_+1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx%v(pnti+psb_proc_id_)
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)
      call psb_get_rank(prcid(i),ictxt,proc_to_comm)      
      if ((nerv>0).and.(proc_to_comm /= me)) then 
        if (debug) write(*,*) me,'Posting receive from',prcid(i),rcv_pt
        p2ptag = psb_complex_swap_tag
        call mpi_irecv(y%combuf(rcv_pt),n*nerv,&
             & psb_mpi_c_spk_,prcid(i),&
             & p2ptag, icomm,y%comid(i,2),iret)
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
    snd_pt = totrcv_+1
    rcv_pt = 1
    do i=1, totxch
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+nerv+psb_n_elem_send_
      call y%gth(idx_pt,snd_pt,nesd,idx)
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
    snd_pt = totrcv_+1
    rcv_pt = 1
    p2ptag = psb_complex_swap_tag
    do i=1, totxch
      proc_to_comm = idx%v(pnti+psb_proc_id_)
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)

      if ((nesd>0).and.(proc_to_comm /= me)) then 
        call mpi_isend(y%combuf(snd_pt),n*nesd,&
             & psb_mpi_c_spk_,prcid(i),&
             & p2ptag,icomm,y%comid(i,1),iret)
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
    if (.not.allocated(y%comid)) then 
      ! 
      ! No matching send? Something is wrong....
      !
      info=psb_err_mpi_error_
      call psb_errpush(info,name,m_err=(/-2/))
      goto 9999
    end if
    call psb_realloc(totxch,prcid,info)

    if (debug) write(*,*) me,' wait'
    pnti   = 1
    snd_pt = totrcv_+1
    rcv_pt = 1
    p2ptag = psb_complex_swap_tag
    do i=1, totxch
      proc_to_comm = idx%v(pnti+psb_proc_id_)
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)
      if (proc_to_comm /= me)then 
        if (nesd>0) then 
          call mpi_wait(y%comid(i,1),p2pstat,iret)
          if(iret /= mpi_success) then
            info=psb_err_mpi_error_
            call psb_errpush(info,name,m_err=(/iret/))
            goto 9999
          end if
        end if
        if (nerv>0) then 
          call mpi_wait(y%comid(i,2),p2pstat,iret)
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
    snd_pt = totrcv_+1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx%v(pnti+psb_proc_id_)
      nerv = idx%v(pnti+psb_n_elem_recv_)
      nesd = idx%v(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_

      if (debug) write(0,*)me,' Received from: ',prcid(i),&
           & y%combuf(rcv_pt:rcv_pt+n*nerv-1)        
      call y%sct(idx_pt,rcv_pt,nerv,idx,beta)
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3
    end do
    !
    ! Waited for com, cleanup comid
    !
    y%comid = mpi_request_null
    
    !
    ! Then wait for device
    !
    if (debug) write(*,*) me,' wait'
    call y%device_wait()
    if (debug) write(*,*) me,' free buffer'
    call y%free_buffer(info)
    if (info == 0) call y%free_comid(info)
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
    if (debug) write(*,*) me,' done'
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(iictxt,err_act)

  return
end subroutine psi_cswap_vidx_multivect

