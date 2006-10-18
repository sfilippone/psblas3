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
subroutine psi_zswaptranm(flag,n,beta,y,desc_a,work,info,data)

  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
  use psi_gthsct_mod
  use mpi
  implicit none

  integer, intent(in)       :: flag, n
  integer, intent(out)      :: info
  complex(kind(1.d0))       :: y(:,:), beta
  complex(kind(1.d0)), target :: work(:)
  type(psb_desc_type)       :: desc_a
  integer, optional         :: data

  ! locals
  integer  :: ictxt, np, me, point_to_proc, nesd, nerv,&
       & proc_to_comm, p2ptag, icomm, p2pstat(mpi_status_size),&
       & idxs, idxr, iret, err_act, totxch, ixrec, i, idx_pt,&
       & snd_pt, rcv_pt

  integer, pointer, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, ptp, rvhd, d_idx
  integer :: int_err(5)
  integer  :: krecvid, ksendid
  logical :: swap_mpi, swap_sync, swap_send, swap_recv, all
  complex(kind(1.d0)), pointer, dimension(:) :: sndbuf, rcvbuf
  character(len=20)  :: name

  info = 0
  name='psi_zswaptranm'
  call psb_erractionsave(err_act)

  ictxt=desc_a%matrix_data(psb_ctxt_)
  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_get_mpicomm(ictxt,icomm)

  allocate(sdsz(0:np-1), rvsz(0:np-1), bsdidx(0:np-1),&
       & brvidx(0:np-1), rvhd(0:np-1), prcid(0:np-1),&
       & ptp(0:np-1), stat=info)
  if(info /= 0) then
    call psb_errpush(4000,name)
    goto 9999
  end if

  swap_mpi  = iand(flag,psb_swap_mpi_) /= 0
  swap_sync = iand(flag,psb_swap_sync_) /= 0
  swap_send = iand(flag,psb_swap_send_) /= 0
  swap_recv = iand(flag,psb_swap_recv_) /= 0

  if(present(data)) then
    if(data == psb_comm_halo_) then
      d_idx => desc_a%halo_index
    else if(data == psb_comm_ovr_) then
      d_idx => desc_a%ovrlap_index
    else
      d_idx => desc_a%halo_index
    end if
  else
    d_idx => desc_a%halo_index
  end if

  idxs = 1
  idxr = 1
  totxch = 0
  point_to_proc = 1
  rvhd(:) = mpi_request_null
  sdsz(:) = 0 
  rvsz(:) = 0 

  ! prepare info for communications
  proc_to_comm = d_idx(point_to_proc+psb_proc_id_)
  do while (proc_to_comm /= -1)
    if(proc_to_comm  /=  me) totxch = totxch+1
    nerv = d_idx(point_to_proc+psb_n_elem_recv_)
    nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

    call psb_get_rank(prcid(proc_to_comm),ictxt,proc_to_comm)
    ptp(proc_to_comm)   = point_to_proc

    brvidx(proc_to_comm) = idxr
    rvsz(proc_to_comm)   = n*nerv
    idxr                 = idxr+rvsz(proc_to_comm)

    bsdidx(proc_to_comm) = idxs
    sdsz(proc_to_comm)   = n*nesd
    idxs                 = idxs+sdsz(proc_to_comm)

    point_to_proc = point_to_proc+nerv+nesd+3
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
  end do

  if((idxr+idxs) < size(work)) then
    sndbuf => work(1:idxs)
    rcvbuf => work(idxs+1:idxs+idxr)
    all=.false.
  else
    allocate(sndbuf(idxs),rcvbuf(idxr), stat=info)
    if(info /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    end if
    all=.true.
  end if

  ! Case SWAP_MPI
  if(swap_mpi) then

    ! gather elements into sendbuffer for swapping
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+psb_elem_recv_
      rcv_pt = brvidx(proc_to_comm)

      call psi_gth(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
           & y,rcvbuf(rcv_pt:rcv_pt+nerv*n-1))

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    ! swap elements using mpi_alltoallv
    call mpi_alltoallv(rcvbuf,rvsz,brvidx,&
         & mpi_double_complex,sndbuf,sdsz,&
         & bsdidx,mpi_double_complex,icomm,iret)
    if(iret /= mpi_success) then
      int_err(1) = iret
      info=400
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    ! scatter elements from receivebuffer after swapping
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+nerv+psb_elem_send_
      snd_pt = bsdidx(proc_to_comm)
      call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
           & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_sync) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)
      if (proc_to_comm  <  me) then
        ! First I send
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_gth(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
             & y,rcvbuf(rcv_pt:rcv_pt+nerv*n-1))
        call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv*n-1), proc_to_comm)
        ! Then I receive
        snd_pt = brvidx(proc_to_comm)
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd*n-1), proc_to_comm)
      else if (proc_to_comm  >  me) then
        ! First I receive
        snd_pt = bsdidx(proc_to_comm)
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd*n-1), proc_to_comm)
        ! Then I send
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_gth(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
             & y,rcvbuf(rcv_pt:rcv_pt+nerv*n-1))
        call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv*n-1), proc_to_comm)
      else if (proc_to_comm  ==  me) then
        ! I send to myself
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = bsdidx(proc_to_comm)
        call psi_gth(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
             & y,rcvbuf(rcv_pt:rcv_pt+nerv*n-1))
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm /= me) then
        idx_pt = point_to_proc+nerv+psb_elem_send_
        snd_pt = bsdidx(proc_to_comm)
        call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
             & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)
      else
        idx_pt = point_to_proc+nerv+psb_elem_send_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
             & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_send .and. swap_recv) then

    ! First I post all the non blocking receives
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm /= me) then
        p2ptag = krecvid(ictxt,proc_to_comm,me)
        snd_pt = brvidx(proc_to_comm)
        call mpi_irecv(sndbuf(snd_pt),sdsz(proc_to_comm),&
             & mpi_double_complex,prcid(proc_to_comm),&
             & p2ptag, icomm,rvhd(proc_to_comm),iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    ! Then I post all the blocking sends
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+psb_elem_recv_
      rcv_pt = brvidx(proc_to_comm)
      call psi_gth(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
           & y,rcvbuf(rcv_pt:rcv_pt+nerv*n-1))

      if(proc_to_comm  /=  me) then
        p2ptag=ksendid(ictxt,proc_to_comm,me)
        call mpi_send(rcvbuf(rcv_pt),rvsz(proc_to_comm),&
             & mpi_double_complex,prcid(proc_to_comm),&
             & p2ptag,icomm,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if
      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    if(.false.) then
      do i=1, totxch
        call mpi_waitany(np,rvhd,ixrec,p2pstat,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if

        if (ixrec  /=  mpi_undefined) then
          ixrec=ixrec-1  ! mpi_waitany returns an 1 to np index
          point_to_proc = ptp(ixrec)
          proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
          nerv = d_idx(point_to_proc+psb_n_elem_recv_)
          nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

          idx_pt = point_to_proc+nerv+psb_elem_send_
          snd_pt = bsdidx(proc_to_comm)
          call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
               & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)
        else
          int_err(1) = ixrec
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end do

      point_to_proc = 1
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      do while (proc_to_comm  /=  -1)
        nerv = d_idx(point_to_proc+psb_n_elem_recv_)
        nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

        if(proc_to_comm  ==  me) then
          idx_pt = point_to_proc+nerv+psb_elem_send_
          rcv_pt = brvidx(proc_to_comm)
          call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
               & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)
        end if

        point_to_proc = point_to_proc+nerv+nesd+3
        proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      end do
    else

      point_to_proc = 1
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      do while (proc_to_comm  /=  -1)
        nerv = d_idx(point_to_proc+psb_n_elem_recv_)
        nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

        if(proc_to_comm /= me) then
          call mpi_wait(rvhd(proc_to_comm),p2pstat,iret)
          if(iret /= mpi_success) then
            int_err(1) = iret
            info=400
            call psb_errpush(info,name,i_err=int_err)
            goto 9999
          end if
          idx_pt = point_to_proc+nerv+psb_elem_send_
          snd_pt = bsdidx(proc_to_comm)
          call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
               & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)
        else
          idx_pt = point_to_proc+nerv+psb_elem_send_
          rcv_pt = brvidx(proc_to_comm)
          call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
               & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)

        end if

        point_to_proc = point_to_proc+nerv+nesd+3
        proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      end do

    end if


  else if (swap_send) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+psb_elem_recv_
      rcv_pt = brvidx(proc_to_comm)
      call psi_gth(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
           & y,rcvbuf(rcv_pt:rcv_pt+nerv*n-1))
      call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv*n-1), proc_to_comm)

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_recv) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm /= me) then
        snd_pt = bsdidx(proc_to_comm)
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd*n-1), proc_to_comm)
        idx_pt = point_to_proc+nerv+psb_elem_send_
        call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
             & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)
      else
        idx_pt = point_to_proc+nerv+psb_elem_send_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
             & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  end if

  deallocate(sdsz,rvsz,bsdidx,&
       & brvidx,rvhd,prcid,&
       & ptp,stat=info)
  if(all) deallocate(sndbuf,rcvbuf,stat=info)
  if(info /= 0) then
    call psb_errpush(4000,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == act_abort) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_zswaptranm






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

subroutine psi_zswaptranv(flag,beta,y,desc_a,work,info,data)

  use psb_error_mod
  use psb_descriptor_type
  use psb_penv_mod
  use psi_gthsct_mod
  use mpi
  implicit none

  integer, intent(in)  :: flag
  integer, intent(out) :: info
  complex(kind(1.d0))     :: y(:), beta
  complex(kind(1.d0)), target :: work(:)
  type(psb_desc_type)  :: desc_a
  integer, optional    :: data

  ! locals
  integer  :: ictxt, np, me, point_to_proc, nesd, nerv,&
       & proc_to_comm, p2ptag, icomm, p2pstat(mpi_status_size),&
       & idxs, idxr, iret, err_act, totxch, ixrec, i, idx_pt,&
       & snd_pt, rcv_pt, n

  integer, pointer, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, ptp, rvhd, d_idx
  integer :: int_err(5)
  integer  :: krecvid, ksendid
  logical :: swap_mpi, swap_sync, swap_send, swap_recv, all
  complex(kind(1.d0)), pointer, dimension(:) :: sndbuf, rcvbuf
  character(len=20)  :: name

  info = 0
  name='psi_zswaptranv'
  call psb_erractionsave(err_act)

  ictxt=desc_a%matrix_data(psb_ctxt_)
  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_get_mpicomm(ictxt,icomm)

  allocate(sdsz(0:np-1), rvsz(0:np-1), bsdidx(0:np-1),&
       & brvidx(0:np-1), rvhd(0:np-1), prcid(0:np-1),&
       & ptp(0:np-1), stat=info)
  if(info /= 0) then
    call psb_errpush(4000,name)
    goto 9999
  end if

  swap_mpi  = iand(flag,psb_swap_mpi_) /= 0
  swap_sync = iand(flag,psb_swap_sync_) /= 0
  swap_send = iand(flag,psb_swap_send_) /= 0
  swap_recv = iand(flag,psb_swap_recv_) /= 0

  if(present(data)) then
    if(data == psb_comm_halo_) then
      d_idx => desc_a%halo_index
    else if(data == psb_comm_ovr_) then
      d_idx => desc_a%ovrlap_index
    else
      d_idx => desc_a%halo_index
    end if
  else
    d_idx => desc_a%halo_index
  end if

  idxs = 1
  idxr = 1
  totxch = 0
  point_to_proc = 1
  rvhd(:) = mpi_request_null
  sdsz(:) = 0 
  rvsz(:) = 0 
  n=1

  ! prepare info for communications
  proc_to_comm = d_idx(point_to_proc+psb_proc_id_)
  do while (proc_to_comm /= -1)
    if(proc_to_comm  /=  me) totxch = totxch+1
    nerv = d_idx(point_to_proc+psb_n_elem_recv_)
    nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

    call psb_get_rank(prcid(proc_to_comm),ictxt,proc_to_comm)
    ptp(proc_to_comm)   = point_to_proc

    brvidx(proc_to_comm) = idxr
    rvsz(proc_to_comm)   = nerv
    idxr                 = idxr+rvsz(proc_to_comm)

    bsdidx(proc_to_comm) = idxs
    sdsz(proc_to_comm)   = nesd
    idxs                 = idxs+sdsz(proc_to_comm)

    point_to_proc = point_to_proc+nerv+nesd+3
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
  end do

  if((idxr+idxs) < size(work)) then
    sndbuf => work(1:idxs)
    rcvbuf => work(idxs+1:idxs+idxr)
    all=.false.
  else
    allocate(sndbuf(idxs),rcvbuf(idxr), stat=info)
    if(info /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    end if
    all=.true.
  end if

  ! Case SWAP_MPI
  if(swap_mpi) then

    ! gather elements into sendbuffer for swapping
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+psb_elem_recv_
      rcv_pt = brvidx(proc_to_comm)
      call psi_gth(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
           & y,rcvbuf(rcv_pt:rcv_pt+nerv-1))

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    ! swap elements using mpi_alltoallv
    call mpi_alltoallv(rcvbuf,rvsz,brvidx,&
         & mpi_double_complex,sndbuf,sdsz,&
         & bsdidx,mpi_double_complex,icomm,iret)
    if(iret /= mpi_success) then
      int_err(1) = iret
      info=400
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    ! scatter elements from receivebuffer after swapping
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+nerv+psb_elem_send_
      snd_pt = bsdidx(proc_to_comm)
      call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
           & sndbuf(snd_pt:snd_pt+nesd-1),beta,y)

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_sync) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)
      if (proc_to_comm  <  me) then
        ! First I send
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_gth(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
             & y,rcvbuf(rcv_pt:rcv_pt+nerv-1))
        call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
        ! Then I receive
        snd_pt = brvidx(proc_to_comm)
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
      else if (proc_to_comm  >  me) then
        ! First I receive
        snd_pt = bsdidx(proc_to_comm)
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
        ! Then I send
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_gth(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
             & y,rcvbuf(rcv_pt:rcv_pt+nerv-1))
        call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)
      else if (proc_to_comm  ==  me) then
        ! I send to myself
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = bsdidx(proc_to_comm)
        call psi_gth(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
             & y,rcvbuf(rcv_pt:rcv_pt+nerv-1))
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm /= me) then
        idx_pt = point_to_proc+nerv+psb_elem_send_
        snd_pt = bsdidx(proc_to_comm)
        call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
             & sndbuf(snd_pt:snd_pt+nesd-1),beta,y)
      else
        idx_pt = point_to_proc+nerv+psb_elem_send_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
             & rcvbuf(rcv_pt:rcv_pt+nerv-1),beta,y)
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_send .and. swap_recv) then

    ! First I post all the non blocking receives
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if (proc_to_comm /= me) then
        p2ptag = krecvid(ictxt,proc_to_comm,me)
        snd_pt = brvidx(proc_to_comm)
        call mpi_irecv(sndbuf(snd_pt),sdsz(proc_to_comm),&
             & mpi_double_complex,prcid(proc_to_comm),&
             & p2ptag, icomm,rvhd(proc_to_comm),iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    ! Then I post all the blocking sends
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+psb_elem_recv_
      rcv_pt = brvidx(proc_to_comm)
      call psi_gth(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
           & y,rcvbuf(rcv_pt:rcv_pt+nerv-1))

      if(proc_to_comm  /=  me) then
        p2ptag=ksendid(ictxt,proc_to_comm,me)
        call mpi_send(rcvbuf(rcv_pt),rvsz(proc_to_comm),&
             & mpi_double_complex,prcid(proc_to_comm),&
             & p2ptag,icomm,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end if
      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    if(.false.) then
      do i=1, totxch
        call mpi_waitany(np,rvhd,ixrec,p2pstat,iret)
        if(iret /= mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if

        if (ixrec  /=  mpi_undefined) then
          ixrec=ixrec-1  ! mpi_waitany returns an 1 to np index
          point_to_proc = ptp(ixrec)
          proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
          nerv = d_idx(point_to_proc+psb_n_elem_recv_)
          nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

          idx_pt = point_to_proc+nerv+psb_elem_send_
          rcv_pt = bsdidx(proc_to_comm)
          call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
               & rcvbuf(rcv_pt:rcv_pt+nerv-1),beta,y)
        else
          int_err(1) = ixrec
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end do

      point_to_proc = 1
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      do while (proc_to_comm  /=  -1)
        nerv = d_idx(point_to_proc+psb_n_elem_recv_)
        nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

        if(proc_to_comm  ==  me) then
          idx_pt = point_to_proc+nerv+psb_elem_send_
          rcv_pt = brvidx(proc_to_comm)
          call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
               & rcvbuf(rcv_pt:rcv_pt+nerv-1),beta,y)
        end if

        point_to_proc = point_to_proc+nerv+nesd+3
        proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      end do

    else

      point_to_proc = 1
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      do while (proc_to_comm  /=  -1)
        nerv = d_idx(point_to_proc+psb_n_elem_recv_)
        nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

        if(proc_to_comm /= me) then
          call mpi_wait(rvhd(proc_to_comm),p2pstat,iret)
          if(iret /= mpi_success) then
            int_err(1) = iret
            info=400
            call psb_errpush(info,name,i_err=int_err)
            goto 9999
          end if
          idx_pt = point_to_proc+nerv+psb_elem_send_
          snd_pt = bsdidx(proc_to_comm)
          call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
               & sndbuf(snd_pt:snd_pt+nesd-1),beta,y)
        else
          idx_pt = point_to_proc+nerv+psb_elem_send_
          rcv_pt = brvidx(proc_to_comm)
          call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
               & rcvbuf(rcv_pt:rcv_pt+nerv-1),beta,y)

        end if

        point_to_proc = point_to_proc+nerv+nesd+3
        proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      end do

    end if

  else if (swap_send) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+psb_elem_recv_
      rcv_pt = brvidx(proc_to_comm)
      call psi_gth(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
           & y,rcvbuf(rcv_pt:rcv_pt+nerv-1))
      call psb_snd(ictxt,rcvbuf(rcv_pt:rcv_pt+nerv-1), proc_to_comm)

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_recv) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm  /=  -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm /= me) then
        snd_pt = bsdidx(proc_to_comm)
        call psb_rcv(ictxt,sndbuf(snd_pt:snd_pt+nesd-1), proc_to_comm)
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
             & sndbuf(snd_pt:snd_pt+nesd-1),beta,y)
      else
        idx_pt = point_to_proc+nerv+psb_elem_send_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
             & rcvbuf(rcv_pt:rcv_pt+nerv-1),beta,y)
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  end if

  deallocate(sdsz,rvsz,bsdidx,&
       & brvidx,rvhd,prcid,&
       & ptp,stat=info)
  if ((info==0).and. (all)) then 
    write(0,*) me, 'Deallocating sndbuf? '
    deallocate(sndbuf,rcvbuf,stat=info)

  end if

  if(info /= 0) then
    write(0,*) 'Error on deallocate ',info
    call psb_errpush(4000,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == act_abort) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_zswaptranv
