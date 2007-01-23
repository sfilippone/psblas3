!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
subroutine psi_iswaptranm(flag,n,beta,y,desc_a,work,info,data)

  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
  use psi_gthsct_mod
  use mpi
  implicit none

  integer, intent(in)      :: flag, n
  integer, intent(out)     :: info
  integer          :: y(:,:), beta
  integer, target :: work(:)
  type(psb_desc_type),target   :: desc_a
  integer, optional        :: data

  ! locals
  integer  :: ictxt, np, me, point_to_proc, nesd, nerv,&
       & proc_to_comm, p2ptag, icomm, p2pstat(mpi_status_size),&
       & idxs, idxr, iret, err_act, totxch, ixrec, i, idx_pt,&
       & snd_pt, rcv_pt, pnti
  integer  :: krecvid, ksendid
  integer, allocatable, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, rvhd, sdhd
  integer, pointer :: d_idx(:)
  integer :: int_err(5)
  logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
       & albf,do_send,do_recv
  logical, parameter :: usersend=.false.

  integer, pointer, dimension(:) :: sndbuf, rcvbuf
  character(len=20)  :: name

  info = 0
  name='psi_swap_tran'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_asb_desc(desc_a)) then 
    info = 1122
    call psb_errpush(info,name)
    goto 9999
  endif

  icomm = desc_a%matrix_data(psb_mpi_c_)

  swap_mpi  = iand(flag,psb_swap_mpi_) /= 0
  swap_sync = iand(flag,psb_swap_sync_) /= 0
  swap_send = iand(flag,psb_swap_send_) /= 0
  swap_recv = iand(flag,psb_swap_recv_) /= 0
  do_send = swap_mpi .or. swap_sync .or. swap_send
  do_recv = swap_mpi .or. swap_sync .or. swap_recv

  if(present(data)) then
    if(data == psb_comm_halo_) then
      d_idx => desc_a%halo_index
      totxch = desc_a%matrix_data(psb_thal_xch_)
      idxr   = desc_a%matrix_data(psb_thal_rcv_)
      idxs   = desc_a%matrix_data(psb_thal_snd_)

    else if(data == psb_comm_ovr_) then
      d_idx => desc_a%ovrlap_index
      totxch = desc_a%matrix_data(psb_tovr_xch_)
      idxr   = desc_a%matrix_data(psb_tovr_rcv_)
      idxs   = desc_a%matrix_data(psb_tovr_snd_)
    else
      d_idx => desc_a%halo_index
      totxch = desc_a%matrix_data(psb_thal_xch_)
      idxr   = desc_a%matrix_data(psb_thal_rcv_)
      idxs   = desc_a%matrix_data(psb_thal_snd_)
    end if
  else
    d_idx => desc_a%halo_index
    totxch = desc_a%matrix_data(psb_thal_xch_)
    idxr   = desc_a%matrix_data(psb_thal_rcv_)
    idxs   = desc_a%matrix_data(psb_thal_snd_)
  end if
  idxr = idxr * n
  idxs = idxs * n

  if (swap_mpi) then 
    allocate(sdsz(0:np-1), rvsz(0:np-1), bsdidx(0:np-1),&
         & brvidx(0:np-1), rvhd(0:np-1), sdhd(0:np-1), prcid(0:np-1),&
         & stat=info)    
    if(info /= 0) then
      call psb_errpush(4000,name)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)

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
    if(info /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    end if
  end if


  if((idxr+idxs) < size(work)) then
    sndbuf => work(1:idxs)
    rcvbuf => work(idxs+1:idxs+idxr)
    albf=.false.
  else
    allocate(sndbuf(idxs),rcvbuf(idxr), stat=info)
    if(info /= 0) then
      call psb_errpush(4000,name)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_

      call psi_gth(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
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
         & mpi_integer,&
         & sndbuf,sdsz,bsdidx,mpi_integer,icomm,iret)
    if(iret /= mpi_success) then
      int_err(1) = iret
      info=400
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

  else if (swap_sync) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)

      if (proc_to_comm  <  me) then
        call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+n*nerv-1), proc_to_comm)
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+n*nesd-1), proc_to_comm)
      else if (proc_to_comm  >  me) then
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+n*nesd-1), proc_to_comm)
        call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+n*nerv-1), proc_to_comm)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      p2ptag = krecvid(ictxt,proc_to_comm,me)
      call psb_get_rank(prcid(i),ictxt,proc_to_comm)      
      call mpi_irecv(sndbuf(snd_pt),n*nesd,&
           & mpi_integer,prcid(i),&
           & p2ptag,icomm,rvhd(i),iret)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      p2ptag=ksendid(ictxt,proc_to_comm,me)

      if (usersend) then 
        call mpi_rsend(rcvbuf(rcv_pt),n*nerv,&
             & mpi_integer,prcid(i),&
             & p2ptag,icomm,iret)
      else
        call mpi_send(rcvbuf(rcv_pt),n*nerv,&
             & mpi_integer,prcid(i),&
             & p2ptag,icomm,iret)
      end if

      if(iret /= mpi_success) then
        int_err(1) = iret
        info=400
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
      end if
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3

    end do


    pnti   = 1
    do i=1, totxch
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)

      p2ptag = krecvid(ictxt,proc_to_comm,me)

      if (proc_to_comm /= me) then
        call mpi_wait(rvhd(i),p2pstat,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if
      pnti   = pnti + nerv + nesd + 3
    end do


  else if (swap_send) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+n*nerv-1), proc_to_comm)
      rcv_pt = rcv_pt + n*nerv
      snd_pt = snd_pt + n*nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  else if (swap_recv) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+n*nesd-1), proc_to_comm)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+nerv+psb_n_elem_send_
      call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
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
  if(info /= 0) then
    call psb_errpush(4000,name)
    goto 9999
  end if
  if(albf) deallocate(sndbuf,rcvbuf,stat=info)
  if(info /= 0) then
    call psb_errpush(4000,name)
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
end subroutine psi_iswaptranm


!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
subroutine psi_iswaptranv(flag,beta,y,desc_a,work,info,data)

  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
  use psi_gthsct_mod
  use mpi
  implicit none

  integer, intent(in)  :: flag
  integer, intent(out) :: info
  integer              :: y(:), beta
  integer, target :: work(:)
  type(psb_desc_type),target  :: desc_a
  integer, optional    :: data

  ! locals
  integer  :: ictxt, np, me, point_to_proc, nesd, nerv,&
       & proc_to_comm, p2ptag, icomm, p2pstat(mpi_status_size),&
       & idxs, idxr, iret, err_act, totxch, ixrec, i, &
       & idx_pt, snd_pt, rcv_pt, n, pnti

  integer, allocatable, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, rvhd, sdhd
  integer, pointer :: d_idx(:)
  integer  :: krecvid, ksendid
  integer :: int_err(5)
  logical :: swap_mpi, swap_sync, swap_send, swap_recv,&
       & albf,do_send,do_recv
  logical, parameter :: usersend=.false.

  integer, pointer, dimension(:) :: sndbuf, rcvbuf
  character(len=20)  :: name

  info = 0
  name='psi_swap_tranv'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_asb_desc(desc_a)) then 
    info = 1122
    call psb_errpush(info,name)
    goto 9999
  endif

  icomm = desc_a%matrix_data(psb_mpi_c_)
  n=1

  swap_mpi  = iand(flag,psb_swap_mpi_) /= 0
  swap_sync = iand(flag,psb_swap_sync_) /= 0
  swap_send = iand(flag,psb_swap_send_) /= 0
  swap_recv = iand(flag,psb_swap_recv_) /= 0
  do_send = swap_mpi .or. swap_sync .or. swap_send
  do_recv = swap_mpi .or. swap_sync .or. swap_recv

  if(present(data)) then
    if(data == psb_comm_halo_) then
      d_idx => desc_a%halo_index
      totxch = desc_a%matrix_data(psb_thal_xch_)
      idxr   = desc_a%matrix_data(psb_thal_rcv_)
      idxs   = desc_a%matrix_data(psb_thal_snd_)

    else if(data == psb_comm_ovr_) then
      d_idx => desc_a%ovrlap_index
      totxch = desc_a%matrix_data(psb_tovr_xch_)
      idxr   = desc_a%matrix_data(psb_tovr_rcv_)
      idxs   = desc_a%matrix_data(psb_tovr_snd_)
    else
      d_idx => desc_a%halo_index
      totxch = desc_a%matrix_data(psb_thal_xch_)
      idxr   = desc_a%matrix_data(psb_thal_rcv_)
      idxs   = desc_a%matrix_data(psb_thal_snd_)
    end if
  else
    d_idx => desc_a%halo_index
    totxch = desc_a%matrix_data(psb_thal_xch_)
    idxr   = desc_a%matrix_data(psb_thal_rcv_)
    idxs   = desc_a%matrix_data(psb_thal_snd_)
  end if


  idxr = idxr * n
  idxs = idxs * n

  if (swap_mpi) then 
    allocate(sdsz(0:np-1), rvsz(0:np-1), bsdidx(0:np-1),&
         & brvidx(0:np-1), rvhd(0:np-1), sdhd(0:np-1), prcid(0:np-1),&
         & stat=info)    
    if(info /= 0) then
      call psb_errpush(4000,name)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
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
    if(info /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    end if
  end if



  if((idxr+idxs) < size(work)) then
    sndbuf => work(1:idxs)
    rcvbuf => work(idxs+1:idxs+idxr)
    albf=.false.
  else
    allocate(sndbuf(idxs),rcvbuf(idxr), stat=info)
    if(info /= 0) then
      call psb_errpush(4000,name)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+psb_n_elem_recv_

      call psi_gth(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
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
         & mpi_integer,&
         & sndbuf,sdsz,bsdidx,mpi_integer,icomm,iret)
    if(iret /= mpi_success) then
      int_err(1) = iret
      info=400
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

  else if (swap_sync) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)

      if (proc_to_comm  <  me) then
        call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
      else if (proc_to_comm  >  me) then
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
        call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      p2ptag = krecvid(ictxt,proc_to_comm,me)
      call psb_get_rank(prcid(i),ictxt,proc_to_comm)      
      call mpi_irecv(sndbuf(snd_pt),nesd,&
           & mpi_integer,prcid(i),&
           & p2ptag,icomm,rvhd(i),iret)

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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      p2ptag=ksendid(ictxt,proc_to_comm,me)

      if (usersend) then 
        call mpi_rsend(rcvbuf(rcv_pt),nerv,&
             & mpi_integer,prcid(i),&
             & p2ptag, icomm,iret)
      else
        call mpi_send(rcvbuf(rcv_pt),nerv,&
             & mpi_integer,prcid(i),&
             & p2ptag, icomm,iret)
      end if

      if(iret /= mpi_success) then
        int_err(1) = iret
        info=400
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
      end if
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do


    pnti   = 1
    do i=1, totxch
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      p2ptag = krecvid(ictxt,proc_to_comm,me)

      if (proc_to_comm /= me) then
        call mpi_wait(rvhd(i),p2pstat,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if
      pnti   = pnti + nerv + nesd + 3
    end do


  else if (swap_send) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3

    end do

  else if (swap_recv) then

    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do i=1, totxch
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
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
      proc_to_comm = d_idx(pnti+psb_proc_id_)
      nerv = d_idx(pnti+psb_n_elem_recv_)
      nesd = d_idx(pnti+nerv+psb_n_elem_send_)
      idx_pt = 1+pnti+nerv+psb_n_elem_send_
      call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
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
  if(info /= 0) then
    call psb_errpush(4000,name)
    goto 9999
  end if
  if(albf) deallocate(sndbuf,rcvbuf,stat=info)
  if(info /= 0) then
    call psb_errpush(4000,name)
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
end subroutine psi_iswaptranv
