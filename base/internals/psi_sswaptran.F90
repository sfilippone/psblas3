!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
!
! File: psi_sswaptran.F90
!
! Subroutine: psi_sswaptranm
!   Does the data exchange among processes. This is similar to Xswapdata, but
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
!   sections SND(Y) and RCV(Y); then we do a SEND(PACK(SND(Y)));
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
subroutine psi_sswaptranm(flag,n,beta,y,desc_a,work,info,data)

  use psi_mod, psb_protect_name => psi_sswaptranm
  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer, intent(in)      :: flag, n
  integer, intent(out)     :: info
  real(psb_spk_)         :: y(:,:), beta
  real(psb_spk_), target :: work(:)
  type(psb_desc_type),target       :: desc_a
  integer, optional         :: data

  ! locals
  integer  :: ictxt, np, me, icomm, idxs, idxr, err_act, totxch, data_
  integer, pointer :: d_idx(:)
  integer :: int_err(5)
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_tran'
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

  call desc_a%get_list(data_,d_idx,totxch,idxr,idxs,info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,a_err='psb_cd_get_list')
    goto 9999
  end if

  call  psi_swaptran(ictxt,icomm,flag,n,beta,y,d_idx,totxch,idxs,idxr,work,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_sswaptranm

subroutine psi_stranidxm(ictxt,icomm,flag,n,beta,y,idx,totxch,totsnd,totrcv,work,info)

  use psi_mod, psb_protect_name => psi_stranidxm
  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer, intent(in)      :: ictxt,icomm,flag,n
  integer, intent(out)     :: info
  real(psb_spk_)         :: y(:,:), beta
  real(psb_spk_), target :: work(:)
  integer, intent(in)      :: idx(:),totxch,totsnd, totrcv

  ! locals
  integer  :: np, me, nesd, nerv,&
       & proc_to_comm, p2ptag, p2pstat(mpi_status_size),&
       & iret, err_act, i, idx_pt, totsnd_, totrcv_,&
       & snd_pt, rcv_pt, pnti, data_
  integer, allocatable, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, rvhd, sdhd
  integer :: int_err(5)
  logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
       & albf,do_send,do_recv
  logical, parameter :: usersend=.false.

  real(psb_spk_), pointer, dimension(:) :: sndbuf, rcvbuf
#ifdef HAVE_VOLATILE
  volatile :: sndbuf, rcvbuf
#endif
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_tran'
  call psb_erractionsave(err_act)

  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif


  swap_mpi  = iand(flag,psb_swap_mpi_) /= 0
  swap_sync = iand(flag,psb_swap_sync_) /= 0
  swap_send = iand(flag,psb_swap_send_) /= 0
  swap_recv = iand(flag,psb_swap_recv_) /= 0
  do_send = swap_mpi .or. swap_sync .or. swap_send
  do_recv = swap_mpi .or. swap_sync .or. swap_recv

  totrcv_ = totrcv * n
  totsnd_ = totsnd * n

  if (swap_mpi) then 
    allocate(sdsz(0:np-1), rvsz(0:np-1), bsdidx(0:np-1),&
         & brvidx(0:np-1), rvhd(0:np-1), sdhd(0:np-1), prcid(0:np-1),&
         & stat=info)    
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if

    rvhd(:) = mpi_request_null
    sdsz(:) = 0 
    rvsz(:) = 0 

    ! prepare info for communications

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)

      call psb_get_rank(prcid(proc_to_comm),ictxt,proc_to_comm)


      brvidx(proc_to_comm) = rcv_pt
      rvsz(proc_to_comm)   = n*nerv

      bsdidx(proc_to_comm) = snd_pt
      sdsz(proc_to_comm)   = n*nesd

      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3
    end do

  else
    allocate(rvhd(totxch),prcid(totxch),stat=info) 
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
  end if


  totrcv_ = max(totrcv_,1)
  totsnd_ = max(totsnd_,1)
  if((totrcv_+totsnd_) < size(work)) then
    sndbuf => work(1:totsnd_)
    rcvbuf => work(totsnd_+1:totsnd_+totrcv_)
    albf=.false.
  else
    allocate(sndbuf(totsnd_),rcvbuf(totrcv_), stat=info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
    albf=.true.
  end if

  if (do_send) then

    ! Pack send buffers
    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_

      call psi_gth(nerv,n,idx(idx_pt:idx_pt+nerv-1),&
           & y,rcvbuf(rcv_pt:rcv_pt+n*nerv-1))

      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3
    end do

  end if



  ! Case SWAP_MPI
  if (swap_mpi) then

    ! swap elements using mpi_alltoallv
    call mpi_alltoallv(rcvbuf,rvsz,brvidx,&
         & mpi_real,&
         & sndbuf,sdsz,bsdidx,mpi_real,icomm,iret)
    if(iret /= mpi_success) then
      int_err(1) = iret
      info=psb_err_mpi_error_
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

  else if (swap_sync) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)

      if (proc_to_comm  <  me) then
        if (nerv>0) call psb_snd(ictxt,&
             & rcvbuf(rcv_pt:rcv_pt+n*nerv-1), proc_to_comm)
        if (nesd>0) call psb_rcv(ictxt,&
             & sndbuf(snd_pt:snd_pt+n*nesd-1), proc_to_comm)
      else if (proc_to_comm  >  me) then
        if (nesd>0) call psb_rcv(ictxt,&
             & sndbuf(snd_pt:snd_pt+n*nesd-1), proc_to_comm)
        if (nerv>0) call psb_snd(ictxt,&
             & rcvbuf(rcv_pt:rcv_pt+n*nerv-1), proc_to_comm)
      else if (proc_to_comm == me) then 
        if (nesd /= nerv) then 
          write(psb_err_unit,*) &
               & 'Fatal error in swaptran: mismatch on self send', &
               & nerv,nesd
        end if
        sndbuf(snd_pt:snd_pt+n*nesd-1) = rcvbuf(rcv_pt:rcv_pt+n*nerv-1) 
      end if
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3

    end do


  else if (swap_send .and. swap_recv) then

    ! First I post all the non blocking receives
    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      call psb_get_rank(prcid(i),ictxt,proc_to_comm)      
      if ((nesd>0).and.(proc_to_comm /= me)) then 
        p2ptag = psb_real_swap_tag
        call mpi_irecv(sndbuf(snd_pt),n*nesd,&
             & mpi_real,prcid(i),&
             & p2ptag,icomm,rvhd(i),iret)
      end if
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3
    end do


    ! Then I post all the blocking sends
    if (usersend)  call mpi_barrier(icomm,info)

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)

      if ((nerv>0).and.(proc_to_comm /= me)) then 
        p2ptag = psb_real_swap_tag
        if (usersend) then 
          call mpi_rsend(rcvbuf(rcv_pt),n*nerv,&
               & mpi_real,prcid(i),&
               & p2ptag,icomm,iret)
        else
          call mpi_send(rcvbuf(rcv_pt),n*nerv,&
               & mpi_real,prcid(i),&
               & p2ptag,icomm,iret)
        end if

        if(iret /= mpi_success) then
          int_err(1) = iret
          info=psb_err_mpi_error_
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3

    end do


    pnti   = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)

      p2ptag = psb_real_swap_tag

      if ((proc_to_comm /= me).and.(nesd>0)) then
        call mpi_wait(rvhd(i),p2pstat,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=psb_err_mpi_error_
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      else if (proc_to_comm == me) then 
        if (nesd /= nerv) then 
          write(psb_err_unit,*) &
               & 'Fatal error in swaptran: mismatch on self send',&
               & nerv,nesd
        end if
        sndbuf(snd_pt:snd_pt+n*nesd-1) = rcvbuf(rcv_pt:rcv_pt+n*nerv-1) 
      end if
      pnti   = pnti + nerv + nesd + 3
    end do


  else if (swap_send) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      if (nerv>0) call psb_snd(ictxt,&
           & rcvbuf(rcv_pt:rcv_pt+n*nerv-1), proc_to_comm)
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  else if (swap_recv) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      if (nesd>0) call psb_rcv(ictxt,&
           & sndbuf(snd_pt:snd_pt+n*nesd-1), proc_to_comm)
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3
    end do

  end if

  if (do_recv) then 

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+nerv+psb_n_elem_send_
      call psi_sct(nesd,n,idx(idx_pt:idx_pt+nesd-1),&
           & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3
    end do

  end if


  if (swap_mpi) then 
    deallocate(sdsz,rvsz,bsdidx,brvidx,rvhd,prcid,sdhd,&
         & stat=info)
  else
    deallocate(rvhd,prcid,stat=info)
  end if
  if(info /= psb_success_) then
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if
  if(albf) deallocate(sndbuf,rcvbuf,stat=info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_stranidxm
!
!
! Subroutine: psi_sswaptranv
!   Does the data exchange among processes. This is similar to Xswapdata, but
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
!   sections SND(Y) and RCV(Y); then we do a SEND(PACK(SND(Y)));
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
!    y(:)     - X                       The data area                        
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
subroutine psi_sswaptranv(flag,beta,y,desc_a,work,info,data)

  use psi_mod, psb_protect_name => psi_sswaptranv
  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer, intent(in)      :: flag
  integer, intent(out)     :: info
  real(psb_spk_)         :: y(:), beta
  real(psb_spk_), target :: work(:)
  type(psb_desc_type),target  :: desc_a
  integer, optional    :: data

  ! locals
  integer  :: ictxt, np, me, icomm, idxs, idxr, totxch, err_act, data_
  integer, pointer :: d_idx(:)
  integer :: int_err(5)
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_tranv'
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

  if (present(data)) then
    data_ = data
  else
    data_ = psb_comm_halo_
  end if
  
  call desc_a%get_list(data_,d_idx,totxch,idxr,idxs,info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,a_err='psb_cd_get_list')
    goto 9999
  end if

  call  psi_swaptran(ictxt,icomm,flag,beta,y,d_idx,totxch,idxs,idxr,work,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_sswaptranv



subroutine psi_stranidxv(ictxt,icomm,flag,beta,y,idx,totxch,totsnd,totrcv,work,info)

  use psi_mod, psb_protect_name => psi_stranidxv
  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer, intent(in)      :: ictxt,icomm,flag
  integer, intent(out)     :: info
  real(psb_spk_)         :: y(:), beta
  real(psb_spk_), target :: work(:)
  integer, intent(in)      :: idx(:),totxch,totsnd, totrcv

  ! locals
  integer  :: np, me, nesd, nerv,&
       & proc_to_comm, p2ptag, p2pstat(mpi_status_size),&
       & iret, err_act, i, idx_pt, totsnd_, totrcv_,&
       & snd_pt, rcv_pt, pnti, data_, n 
  integer, allocatable, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, rvhd, sdhd
  integer :: int_err(5)
  logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
       & albf,do_send,do_recv
  logical, parameter :: usersend=.false.

  real(psb_spk_), pointer, dimension(:) :: sndbuf, rcvbuf
#ifdef HAVE_VOLATILE
  volatile :: sndbuf, rcvbuf
#endif
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_tran'
  call psb_erractionsave(err_act)

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

  if (swap_mpi) then 
    allocate(sdsz(0:np-1), rvsz(0:np-1), bsdidx(0:np-1),&
         & brvidx(0:np-1), rvhd(0:np-1), sdhd(0:np-1), prcid(0:np-1),&
         & stat=info)    
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if

    rvhd(:) = mpi_request_null
    sdsz(:) = 0 
    rvsz(:) = 0 

    ! prepare info for communications


    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      call psb_get_rank(prcid(proc_to_comm),ictxt,proc_to_comm)

      brvidx(proc_to_comm) = rcv_pt
      rvsz(proc_to_comm)   = nerv

      bsdidx(proc_to_comm) = snd_pt
      sdsz(proc_to_comm)   = nesd

      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  else
    allocate(rvhd(totxch),prcid(totxch),stat=info) 
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
  end if


  totrcv_ = max(totrcv_,1)
  totsnd_ = max(totsnd_,1)
  if((totrcv_+totsnd_) < size(work)) then
    sndbuf => work(1:totsnd_)
    rcvbuf => work(totsnd_+1:totsnd_+totrcv_)
    albf=.false.
  else
    allocate(sndbuf(totsnd_),rcvbuf(totrcv_), stat=info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
    albf=.true.
  end if


  if (do_send) then

    ! Pack send buffers
    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_

      call psi_gth(nerv,idx(idx_pt:idx_pt+nerv-1),&
           & y,rcvbuf(rcv_pt:rcv_pt+nerv-1))

      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3
    end do

  end if

  ! Case SWAP_MPI
  if (swap_mpi) then

    ! swap elements using mpi_alltoallv
    call mpi_alltoallv(rcvbuf,rvsz,brvidx,&
         & mpi_real,&
         & sndbuf,sdsz,bsdidx,mpi_real,icomm,iret)
    if(iret /= mpi_success) then
      int_err(1) = iret
      info=psb_err_mpi_error_
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

  else if (swap_sync) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)

      if (proc_to_comm  <  me) then
        if (nerv>0) call psb_snd(ictxt,&
             & rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
        if (nesd>0) call psb_rcv(ictxt,&
             & sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
      else if (proc_to_comm  >  me) then
        if (nesd>0) call psb_rcv(ictxt,&
             & sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
        if (nerv>0) call psb_snd(ictxt,&
             & rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
      else if (proc_to_comm ==  me) then
        if (nesd /= nerv) then 
          write(psb_err_unit,*) &
               & 'Fatal error in swaptran: mismatch on self send', &
               & nerv,nesd
        end if
        sndbuf(snd_pt:snd_pt+nesd-1) = rcvbuf(rcv_pt:rcv_pt+nerv-1) 
      end if
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do


  else if (swap_send .and. swap_recv) then

    ! First I post all the non blocking receives
    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      call psb_get_rank(prcid(i),ictxt,proc_to_comm)      
      if ((nesd>0).and.(proc_to_comm /= me)) then 
        p2ptag = psb_real_swap_tag
        call mpi_irecv(sndbuf(snd_pt),nesd,&
             & mpi_real,prcid(i),&
             & p2ptag,icomm,rvhd(i),iret)
      end if
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3
    end do


    ! Then I post all the blocking sends
    if (usersend)  call mpi_barrier(icomm,info)

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)

      if ((nerv>0).and.(proc_to_comm /= me)) then 
        p2ptag = psb_real_swap_tag
        if (usersend) then 
          call mpi_rsend(rcvbuf(rcv_pt),nerv,&
               & mpi_real,prcid(i),&
               & p2ptag, icomm,iret)
        else
          call mpi_send(rcvbuf(rcv_pt),nerv,&
               & mpi_real,prcid(i),&
               & p2ptag, icomm,iret)
        end if

        if(iret /= mpi_success) then
          int_err(1) = iret
          info=psb_err_mpi_error_
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do


    pnti   = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      p2ptag = psb_real_swap_tag

      if ((proc_to_comm /= me).and.(nesd>0)) then
        call mpi_wait(rvhd(i),p2pstat,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=psb_err_mpi_error_
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      else if (proc_to_comm ==  me) then
        if (nesd /= nerv) then 
          write(psb_err_unit,*) &
               & 'Fatal error in swaptran: mismatch on self send', &
               & nerv,nesd
        end if
        sndbuf(snd_pt:snd_pt+nesd-1) = rcvbuf(rcv_pt:rcv_pt+nerv-1) 
      end if
      pnti   = pnti + nerv + nesd + 3
    end do


  else if (swap_send) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      if (nerv>0) call psb_snd(ictxt,&
           & rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  else if (swap_recv) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      if (nesd>0) call psb_rcv(ictxt,&
           & sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  end if


  if (do_recv) then 

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+nerv+psb_n_elem_send_
      call psi_sct(nesd,idx(idx_pt:idx_pt+nesd-1),&
           & sndbuf(snd_pt:snd_pt+nesd-1),beta,y)
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3
    end do

  end if


  if (swap_mpi) then 
    deallocate(sdsz,rvsz,bsdidx,brvidx,rvhd,prcid,sdhd,&
         & stat=info)
  else
    deallocate(rvhd,prcid,stat=info)
  end if
  if(info /= psb_success_) then
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if
  if(albf) deallocate(sndbuf,rcvbuf,stat=info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_stranidxv


subroutine psi_sswaptran_vect(flag,beta,y,desc_a,work,info,data)

  use psi_mod, psb_protect_name => psi_sswaptran_vect
  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
  use psb_s_base_vect_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer, intent(in)         :: flag
  integer, intent(out)        :: info
  class(psb_s_base_vect_type) :: y
  real(psb_spk_)           :: beta
  real(psb_spk_), target   :: work(:)
  type(psb_desc_type),target  :: desc_a
  integer, optional    :: data

  ! locals
  integer  :: ictxt, np, me, icomm, idxs, idxr, totxch, err_act, data_
  integer, pointer :: d_idx(:)
  integer :: int_err(5)
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_tranv'
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

  if (present(data)) then
    data_ = data
  else
    data_ = psb_comm_halo_
  end if
  
  call desc_a%get_list(data_,d_idx,totxch,idxr,idxs,info) 
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name,a_err='psb_cd_get_list')
    goto 9999
  end if

  call  psi_swaptran(ictxt,icomm,flag,beta,y,d_idx,totxch,idxs,idxr,work,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_sswaptran_vect



subroutine psi_stranidx_vect(ictxt,icomm,flag,beta,y,idx,&
     & totxch,totsnd,totrcv,work,info)

  use psi_mod, psb_protect_name => psi_stranidx_vect
  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
  use psb_s_base_vect_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  integer, intent(in)      :: ictxt,icomm,flag
  integer, intent(out)     :: info
  class(psb_s_base_vect_type) :: y
  real(psb_spk_)         :: beta
  real(psb_spk_), target :: work(:)
  integer, intent(in)      :: idx(:),totxch,totsnd, totrcv

  ! locals
  integer  :: np, me, nesd, nerv,&
       & proc_to_comm, p2ptag, p2pstat(mpi_status_size),&
       & iret, err_act, i, idx_pt, totsnd_, totrcv_,&
       & snd_pt, rcv_pt, pnti, data_, n 
  integer, allocatable, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, rvhd, sdhd
  integer :: int_err(5)
  logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
       & albf,do_send,do_recv
  logical, parameter :: usersend=.false.

  real(psb_spk_), pointer, dimension(:) :: sndbuf, rcvbuf
#ifdef HAVE_VOLATILE
  volatile :: sndbuf, rcvbuf
#endif
  character(len=20)  :: name

  info=psb_success_
  name='psi_swap_tran'
  call psb_erractionsave(err_act)

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

  if (swap_mpi) then 
    allocate(sdsz(0:np-1), rvsz(0:np-1), bsdidx(0:np-1),&
         & brvidx(0:np-1), rvhd(0:np-1), sdhd(0:np-1), prcid(0:np-1),&
         & stat=info)    
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if

    rvhd(:) = mpi_request_null
    sdsz(:) = 0 
    rvsz(:) = 0 

    ! prepare info for communications


    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      call psb_get_rank(prcid(proc_to_comm),ictxt,proc_to_comm)

      brvidx(proc_to_comm) = rcv_pt
      rvsz(proc_to_comm)   = nerv

      bsdidx(proc_to_comm) = snd_pt
      sdsz(proc_to_comm)   = nesd

      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  else
    allocate(rvhd(totxch),prcid(totxch),stat=info) 
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
  end if


  totrcv_ = max(totrcv_,1)
  totsnd_ = max(totsnd_,1)
  if((totrcv_+totsnd_) < size(work)) then
    sndbuf => work(1:totsnd_)
    rcvbuf => work(totsnd_+1:totsnd_+totrcv_)
    albf=.false.
  else
    allocate(sndbuf(totsnd_),rcvbuf(totrcv_), stat=info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
    albf=.true.
  end if


  if (do_send) then

    ! Pack send buffers
    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_

      call y%gth(nerv,idx(idx_pt:idx_pt+nerv-1),&
           & rcvbuf(rcv_pt:rcv_pt+nerv-1))

      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3
    end do

  end if

  ! Case SWAP_MPI
  if (swap_mpi) then

    ! swap elements using mpi_alltoallv
    call mpi_alltoallv(rcvbuf,rvsz,brvidx,&
         & mpi_real,&
         & sndbuf,sdsz,bsdidx,mpi_real,icomm,iret)
    if(iret /= mpi_success) then
      int_err(1) = iret
      info=psb_err_mpi_error_
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

  else if (swap_sync) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)

      if (proc_to_comm  <  me) then
        if (nerv>0) call psb_snd(ictxt,&
             & rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
        if (nesd>0) call psb_rcv(ictxt,&
             & sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
      else if (proc_to_comm  >  me) then
        if (nesd>0) call psb_rcv(ictxt,&
             & sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
        if (nerv>0) call psb_snd(ictxt,&
             & rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
      else if (proc_to_comm ==  me) then
        if (nesd /= nerv) then 
          write(psb_err_unit,*) &
               & 'Fatal error in swaptran: mismatch on self send',& 
               & nerv,nesd
        end if
        sndbuf(snd_pt:snd_pt+nesd-1) = rcvbuf(rcv_pt:rcv_pt+nerv-1) 
      end if
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do


  else if (swap_send .and. swap_recv) then

    ! First I post all the non blocking receives
    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      call psb_get_rank(prcid(i),ictxt,proc_to_comm)      
      if ((nesd>0).and.(proc_to_comm /= me)) then 
        p2ptag = psb_real_swap_tag
        call mpi_irecv(sndbuf(snd_pt),nesd,&
             & mpi_real,prcid(i),&
             & p2ptag,icomm,rvhd(i),iret)
      end if
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3
    end do


    ! Then I post all the blocking sends
    if (usersend)  call mpi_barrier(icomm,info)

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)

      if ((nerv>0).and.(proc_to_comm /= me)) then 
        p2ptag = psb_real_swap_tag
        if (usersend) then 
          call mpi_rsend(rcvbuf(rcv_pt),nerv,&
               & mpi_real,prcid(i),&
               & p2ptag, icomm,iret)
        else
          call mpi_send(rcvbuf(rcv_pt),nerv,&
               & mpi_real,prcid(i),&
               & p2ptag, icomm,iret)
        end if

        if(iret /= mpi_success) then
          int_err(1) = iret
          info=psb_err_mpi_error_
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do


    pnti   = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      p2ptag = psb_real_swap_tag

      if ((proc_to_comm /= me).and.(nesd>0)) then
        call mpi_wait(rvhd(i),p2pstat,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=psb_err_mpi_error_
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      else if (proc_to_comm ==  me) then
        if (nesd /= nerv) then 
          write(psb_err_unit,*) &
               & 'Fatal error in swaptran: mismatch on self send', &
               & nerv,nesd
        end if
        sndbuf(snd_pt:snd_pt+nesd-1) = rcvbuf(rcv_pt:rcv_pt+nerv-1) 
      end if
      pnti   = pnti + nerv + nesd + 3
    end do


  else if (swap_send) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      if (nerv>0) call psb_snd(ictxt,&
           & rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  else if (swap_recv) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      if (nesd>0) call psb_rcv(ictxt,&
           & sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  end if


  if (do_recv) then 

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = idx(pnti+psb_proc_id_)
      nerv = idx(pnti+psb_n_elem_recv_)
      nesd = idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+nerv+psb_n_elem_send_
      call y%sct(nesd,idx(idx_pt:idx_pt+nesd-1),&
           & sndbuf(snd_pt:snd_pt+nesd-1),beta)
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3
    end do

  end if


  if (swap_mpi) then 
    deallocate(sdsz,rvsz,bsdidx,brvidx,rvhd,prcid,sdhd,&
         & stat=info)
  else
    deallocate(rvhd,prcid,stat=info)
  end if
  if(info /= psb_success_) then
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if
  if(albf) deallocate(sndbuf,rcvbuf,stat=info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_stranidx_vect



