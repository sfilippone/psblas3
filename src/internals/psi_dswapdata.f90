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
subroutine psi_dswapdatam(flag,n,beta,y,desc_a,work,info,data)

  use psb_error_mod
  use psb_descriptor_type
  use mpi
  implicit none

  integer, intent(in)      :: flag, n
  integer, intent(out)     :: info
  real(kind(1.d0))         :: y(:,:), beta
  real(kind(1.d0)), target :: work(:)
  type(psb_desc_type)      :: desc_a
  integer, optional        :: data

  ! locals
  integer  :: icontxt, nprow, npcol, myrow,&
       & mycol, point_to_proc, nesd, nerv,&
       & proc_to_comm, p2ptag, icomm, p2pstat(mpi_status_size),&
       & idxs, idxr, iret, errlen, ifcomm, rank,&
       & err_act, totxch, ixrec, i, lw, idx_pt,&
       & snd_pt, rcv_pt
  integer  :: blacs_pnum, krecvid, ksendid
  integer, pointer, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, ptp, rvhd, d_idx
  integer :: int_err(5)
  logical :: swap_mpi, swap_sync, swap_send, swap_recv, all
  real(kind(1.d0)), pointer, dimension(:) :: sndbuf, rcvbuf
  character(len=20)  :: name, ch_err

  interface psi_gth
    subroutine psi_dgthm(n,k,idx,x,y)
      integer :: n, k, idx(:)
      real(kind(1.d0)) :: x(:,:), y(:)
    end subroutine psi_dgthm
    subroutine psi_dgthv(n,idx,x,y)
      integer :: n, idx(:)
      real(kind(1.d0)) :: x(:), y(:)
    end subroutine psi_dgthv
    subroutine psi_igthm(n,k,idx,x,y)
      integer :: n, k, idx(:)
      integer :: x(:,:), y(:)
    end subroutine psi_igthm
    subroutine psi_igthv(n,idx,x,y)
      integer :: n, idx(:)
      integer :: x(:), y(:)
    end subroutine psi_igthv
  end interface

  interface psi_sct
    subroutine psi_dsctm(n,k,idx,x,beta,y)
      integer :: n, k, idx(:)
      real(kind(1.d0)) :: beta, x(:), y(:,:)
    end subroutine psi_dsctm
    subroutine psi_dsctv(n,idx,x,beta,y)
      integer :: n, idx(:)
      real(kind(1.d0)) :: beta, x(:), y(:)
    end subroutine psi_dsctv
    subroutine psi_isctm(n,k,idx,x,beta,y)
      integer :: n, k, idx(:)
      integer :: beta, x(:), y(:,:)
    end subroutine psi_isctm
    subroutine psi_isctv(n,idx,x,beta,y)
      integer :: n, idx(:)
      integer :: beta, x(:), y(:)
    end subroutine psi_isctv
  end interface

  info = 0
  name='psi_dswap_data'
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol) 
  if (nprow == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= 1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  call blacs_get(icontxt,10,icomm)

  allocate(sdsz(0:nprow-1), rvsz(0:nprow-1), bsdidx(0:nprow-1),&
       & brvidx(0:nprow-1), rvhd(0:nprow-1), prcid(0:nprow-1),&
       & ptp(0:nprow-1), stat=info)
  if(info.ne.0) then
    call psb_errpush(4000,name)
    goto 9999
  end if

  swap_mpi  = iand(flag,psb_swap_mpi_) .ne.0
  swap_sync = iand(flag,psb_swap_sync_).ne.0
  swap_send = iand(flag,psb_swap_send_).ne.0
  swap_recv = iand(flag,psb_swap_recv_).ne.0

  if(present(data)) then
    if(data.eq.psb_comm_halo_) then
      d_idx => desc_a%halo_index
    else if(data.eq.psb_comm_ovr_) then
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

  ! prepare info for communications
  proc_to_comm = d_idx(point_to_proc+psb_proc_id_)
  do while (proc_to_comm.ne.-1)
    if(proc_to_comm .ne. myrow) totxch = totxch+1
    nerv = d_idx(point_to_proc+psb_n_elem_recv_)
    nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

    prcid(proc_to_comm) = blacs_pnum(icontxt,proc_to_comm,mycol)
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

  if((idxr+idxs).lt.size(work)) then
    sndbuf => work(1:idxs)
    rcvbuf => work(idxs+1:idxs+idxr)
    all=.false.
  else
    allocate(sndbuf(idxs),rcvbuf(idxr), stat=info)
    if(info.ne.0) then
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
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+nerv+psb_elem_send_
      snd_pt = bsdidx(proc_to_comm)

      call psi_gth(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
           & y,sndbuf(snd_pt:snd_pt+nesd*n-1))

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    ! swap elements using mpi_alltoallv
    call mpi_alltoallv(sndbuf,sdsz,bsdidx,&
         & mpi_double_precision,rcvbuf,rvsz,&
         & brvidx,mpi_double_precision,icomm,iret)
    if(iret.ne.mpi_success) then
      int_err(1) = iret
      info=400
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    ! scatter elements from receivebuffer after swapping
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+psb_elem_recv_
      rcv_pt = brvidx(proc_to_comm)
      call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
           & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_sync) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)
      if (proc_to_comm .lt. myrow) then
        ! First I send
        idx_pt = point_to_proc+nerv+psb_elem_send_
        snd_pt = bsdidx(proc_to_comm)
        call psi_gth(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
             & y,sndbuf(snd_pt:snd_pt+nesd*n-1))
        call dgesd2d(icontxt,nesd,n,sndbuf(snd_pt),nesd,proc_to_comm,0)
        ! Then I receive
        rcv_pt = brvidx(proc_to_comm)
        call dgerv2d(icontxt,nerv,n,rcvbuf(rcv_pt),nerv,proc_to_comm,0)
      else if (proc_to_comm .gt. myrow) then
        ! First I receive
        rcv_pt = brvidx(proc_to_comm)
        call dgerv2d(icontxt,nerv,n,rcvbuf(rcv_pt),nerv,proc_to_comm,0)
        ! Then I send
        idx_pt = point_to_proc+nerv+psb_elem_send_
        snd_pt = bsdidx(proc_to_comm)
        call psi_gth(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
             & y,sndbuf(snd_pt:snd_pt+nesd*n-1))
        call dgesd2d(icontxt,nesd,n,sndbuf(snd_pt),nesd,proc_to_comm,0)
      else if (proc_to_comm .eq. myrow) then
        ! I send to myself
        idx_pt = point_to_proc+nerv+psb_elem_send_
        snd_pt = bsdidx(proc_to_comm)
        call psi_gth(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
             & y,sndbuf(snd_pt:snd_pt+nesd*n-1))
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm.ne.myrow) then
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
             & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)
      else
        idx_pt = point_to_proc+psb_elem_recv_
        snd_pt = bsdidx(proc_to_comm)
        call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
             & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_send .and. swap_recv) then
    ! First I post all the non blocking receives
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm.ne.myrow) then
        p2ptag = krecvid(icontxt,proc_to_comm,myrow)
        rcv_pt = brvidx(proc_to_comm)
        call mpi_irecv(rcvbuf(rcv_pt),rvsz(proc_to_comm),&
             & mpi_double_precision,prcid(proc_to_comm),&
             & p2ptag, icomm,rvhd(proc_to_comm),iret)
        if(iret.ne.mpi_success) then
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
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+nerv+psb_elem_send_
      snd_pt = bsdidx(proc_to_comm)

      call psi_gth(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
           & y,sndbuf(snd_pt:snd_pt+nesd*n-1))

      if(proc_to_comm .ne. myrow) then
        p2ptag=ksendid(icontxt,proc_to_comm,myrow)
        call mpi_send(sndbuf(snd_pt),sdsz(proc_to_comm),&
             & mpi_double_precision,prcid(proc_to_comm),&
             & p2ptag,icomm,iret)
        if(iret.ne.mpi_success) then
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
        call mpi_waitany(nprow,rvhd,ixrec,p2pstat,iret)
        if(iret.ne.mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if

        if (ixrec .ne. mpi_undefined) then
          ixrec=ixrec-1  ! mpi_waitany returns an 1 to nprow index
          point_to_proc = ptp(ixrec)
          proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
          nerv = d_idx(point_to_proc+psb_n_elem_recv_)
          nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

          idx_pt = point_to_proc+psb_elem_recv_
          rcv_pt = brvidx(proc_to_comm)
          call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
               & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)
        else
          int_err(1) = ixrec
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
      end do

      point_to_proc = 1
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      do while (proc_to_comm .ne. -1)
        nerv = d_idx(point_to_proc+psb_n_elem_recv_)
        nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

        if(proc_to_comm .eq. myrow) then
          idx_pt = point_to_proc+psb_elem_recv_
          snd_pt = bsdidx(proc_to_comm)
          call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
               & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)
        end if

        point_to_proc = point_to_proc+nerv+nesd+3
        proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      end do
    else

      point_to_proc = 1
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      do while (proc_to_comm .ne. -1)
        nerv = d_idx(point_to_proc+psb_n_elem_recv_)
        nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

        if(proc_to_comm.ne.myrow) then
          call mpi_wait(rvhd(proc_to_comm),p2pstat,iret)
          if(iret.ne.mpi_success) then
            int_err(1) = iret
            info=400
            call psb_errpush(info,name,i_err=int_err)
            goto 9999
          end if
          idx_pt = point_to_proc+psb_elem_recv_
          rcv_pt = brvidx(proc_to_comm)
          call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
               & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)
        else
          idx_pt = point_to_proc+psb_elem_recv_
          snd_pt = bsdidx(proc_to_comm)
          call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
               & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)

        end if

        point_to_proc = point_to_proc+nerv+nesd+3
        proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      end do

    end if


  else if (swap_send) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+nerv+psb_elem_send_
      snd_pt = bsdidx(proc_to_comm)
      call psi_gth(nesd,n,d_idx(idx_pt:idx_pt+nesd-1),&
           & y,sndbuf(snd_pt:snd_pt+nesd*n-1))
      call dgesd2d(icontxt,nesd,n,sndbuf(snd_pt),nesd,proc_to_comm,0)

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_recv) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm.ne.myrow) then
        rcv_pt = brvidx(proc_to_comm)
        call dgerv2d(icontxt,nerv,n,rcvbuf(rcv_pt),nerv,proc_to_comm,0)
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
             & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)
      else
        idx_pt = point_to_proc+psb_elem_recv_
        snd_pt = bsdidx(proc_to_comm)
        call psi_sct(nerv,n,d_idx(idx_pt:idx_pt+nerv-1),&
             & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  end if

  deallocate(sdsz,rvsz,bsdidx,&
       & brvidx,rvhd,prcid,&
       & ptp,stat=info)
  if(all) deallocate(sndbuf,rcvbuf,stat=info)
  if(info.ne.0) then
    call psb_errpush(4000,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error(icontxt)
    return
  end if
  return
end subroutine psi_dswapdatam


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
subroutine psi_dswapdatav(flag,beta,y,desc_a,work,info,data)

  use psb_error_mod
  use psb_descriptor_type
  use mpi
  implicit none

  integer, intent(in)      :: flag
  integer, intent(out)     :: info
  real(kind(1.d0))         :: y(:), beta
  real(kind(1.d0)), target :: work(:)
  type(psb_desc_type)      :: desc_a
  integer, optional        :: data

  ! locals
  integer  :: icontxt, nprow, npcol, myrow,&
       & mycol, point_to_proc, nesd, nerv,&
       & proc_to_comm, p2ptag, icomm, p2pstat(mpi_status_size),&
       & idxs, idxr, iret, errlen, ifcomm, rank,&
       & err_act, totxch, ixrec, i, lw, idx_pt,&
       & snd_pt, rcv_pt, n

  integer, pointer, dimension(:) :: bsdidx, brvidx,&
       & sdsz, rvsz, prcid, ptp, rvhd, d_idx
  integer  :: blacs_pnum, krecvid, ksendid
  integer :: int_err(5)
  logical :: swap_mpi, swap_sync, swap_send, swap_recv,all
  real(kind(1.d0)), pointer, dimension(:) :: sndbuf, rcvbuf
  character(len=20)  :: name, ch_err

  interface psi_gth
    subroutine psi_dgthm(n,k,idx,x,y)
      integer :: n, k, idx(:)
      real(kind(1.d0)) :: x(:,:), y(:)
    end subroutine psi_dgthm
    subroutine psi_dgthv(n,idx,x,y)
      integer :: n, idx(:)
      real(kind(1.d0)) :: x(:), y(:)
    end subroutine psi_dgthv
    subroutine psi_igthm(n,k,idx,x,y)
      integer :: n, k, idx(:)
      integer :: x(:,:), y(:)
    end subroutine psi_igthm
    subroutine psi_igthv(n,idx,x,y)
      integer :: n, idx(:)
      integer :: x(:), y(:)
    end subroutine psi_igthv
  end interface

  interface psi_sct
    subroutine psi_dsctm(n,k,idx,x,beta,y)
      integer :: n, k, idx(:)
      real(kind(1.d0)) :: beta, x(:), y(:,:)
    end subroutine psi_dsctm
    subroutine psi_dsctv(n,idx,x,beta,y)
      integer :: n, idx(:)
      real(kind(1.d0)) :: beta, x(:), y(:)
    end subroutine psi_dsctv
    subroutine psi_isctm(n,k,idx,x,beta,y)
      integer :: n, k, idx(:)
      integer :: beta, x(:), y(:,:)
    end subroutine psi_isctm
    subroutine psi_isctv(n,idx,x,beta,y)
      integer :: n, idx(:)
      integer :: beta, x(:), y(:)
    end subroutine psi_isctv
  end interface

  info = 0
  name='psi_dswap_datav'
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol) 
  if (nprow == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= 1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  call blacs_get(icontxt,10,icomm)


  allocate(sdsz(0:nprow-1), rvsz(0:nprow-1), bsdidx(0:nprow-1),&
       & brvidx(0:nprow-1), rvhd(0:nprow-1), prcid(0:nprow-1),&
       & ptp(0:nprow-1), stat=info)
  if(info.ne.0) then
    call psb_errpush(4000,name)
    goto 9999
  end if

  swap_mpi  = iand(flag,psb_swap_mpi_).ne.0
  swap_sync = iand(flag,psb_swap_sync_).ne.0
  swap_send = iand(flag,psb_swap_send_).ne.0
  swap_recv = iand(flag,psb_swap_recv_).ne.0

  if(present(data)) then
    if(data.eq.psb_comm_halo_) then
      d_idx => desc_a%halo_index
    else if(data.eq.psb_comm_ovr_) then
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
  n=1

  ! prepare info for communications
  proc_to_comm = d_idx(point_to_proc+psb_proc_id_)
  do while (proc_to_comm.ne.-1)
    if(proc_to_comm .ne. myrow) totxch = totxch+1
    nerv = d_idx(point_to_proc+psb_n_elem_recv_)
    nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

    prcid(proc_to_comm) = blacs_pnum(icontxt,proc_to_comm,mycol)
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

  if((idxr+idxs).lt.size(work)) then
    sndbuf => work(1:idxs)
    rcvbuf => work(idxs+1:idxs+idxr)
    all=.false.
  else
    allocate(sndbuf(idxs),rcvbuf(idxr), stat=info)
    if(info.ne.0) then
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
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+nerv+psb_elem_send_
      snd_pt = bsdidx(proc_to_comm)
      call psi_gth(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
           & y,sndbuf(snd_pt:snd_pt+nesd-1))

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    ! swap elements using mpi_alltoallv
    call mpi_alltoallv(sndbuf,sdsz,bsdidx,&
         & mpi_double_precision,rcvbuf,rvsz,&
         & brvidx,mpi_double_precision,icomm,iret)
    if(iret.ne.mpi_success) then
      int_err(1) = iret
      info=400
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    ! scatter elements from receivebuffer after swapping
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+psb_elem_recv_
      rcv_pt = brvidx(proc_to_comm)
      call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
           & rcvbuf(rcv_pt:rcv_pt+nerv-1),beta,y)

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_sync) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)
      if (proc_to_comm .lt. myrow) then
        ! First I send
        idx_pt = point_to_proc+nerv+psb_elem_send_
        snd_pt = bsdidx(proc_to_comm)
        call psi_gth(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
             & y,sndbuf(snd_pt:snd_pt+nesd-1))
        call dgesd2d(icontxt,nesd,1,sndbuf(snd_pt),nesd,proc_to_comm,0)
        ! Then I receive
        rcv_pt = brvidx(proc_to_comm)
        call dgerv2d(icontxt,nerv,1,rcvbuf(rcv_pt),nerv,proc_to_comm,0)
      else if (proc_to_comm .gt. myrow) then
        ! First I receive
        rcv_pt = brvidx(proc_to_comm)
        call dgerv2d(icontxt,nerv,1,rcvbuf(rcv_pt),nerv,proc_to_comm,0)
        ! Then I send
        idx_pt = point_to_proc+nerv+psb_elem_send_
        snd_pt = bsdidx(proc_to_comm)
        call psi_gth(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
             & y,sndbuf(snd_pt:snd_pt+nesd-1))
        call dgesd2d(icontxt,nesd,1,sndbuf(snd_pt),nesd,proc_to_comm,0)
      else if (proc_to_comm .eq. myrow) then
        ! I send to myself
        idx_pt = point_to_proc+nerv+psb_elem_send_
        snd_pt = bsdidx(proc_to_comm)
        call psi_gth(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
             & y,sndbuf(snd_pt:snd_pt+nesd-1))
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm.ne.myrow) then
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
             & rcvbuf(rcv_pt:rcv_pt+nerv-1),beta,y)
      else
        idx_pt = point_to_proc+psb_elem_recv_
        snd_pt = bsdidx(proc_to_comm)
        call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
             & sndbuf(snd_pt:snd_pt+nesd-1),beta,y)
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_send .and. swap_recv) then
    ! First I post all the non blocking receives
    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm.ne.myrow) then
        p2ptag = krecvid(icontxt,proc_to_comm,myrow)
        rcv_pt = brvidx(proc_to_comm)
        call mpi_irecv(rcvbuf(rcv_pt),rvsz(proc_to_comm),&
             & mpi_double_precision,prcid(proc_to_comm),&
             & p2ptag, icomm,rvhd(proc_to_comm),iret)
        if(iret.ne.mpi_success) then
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
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+nerv+psb_elem_send_
      snd_pt = bsdidx(proc_to_comm)
      call psi_gth(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
           & y,sndbuf(snd_pt:snd_pt+nesd-1))

      if(proc_to_comm .ne. myrow) then
        p2ptag=ksendid(icontxt,proc_to_comm,myrow)
        call mpi_send(sndbuf(snd_pt),sdsz(proc_to_comm),&
             & mpi_double_precision,prcid(proc_to_comm),&
             & p2ptag,icomm,iret)
        if(iret.ne.mpi_success) then
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
        call mpi_waitany(nprow,rvhd,ixrec,p2pstat,iret)
        if(iret.ne.mpi_success) then
          int_err(1) = iret
          info=400
          call psb_errpush(info,name,i_err=int_err)
          goto 9999
        end if
        if (ixrec .ne. mpi_undefined) then
          ixrec=ixrec-1  ! mpi_waitany returns an 1 to nprow index
          point_to_proc = ptp(ixrec)
          proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
          nerv = d_idx(point_to_proc+psb_n_elem_recv_)
          nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

          idx_pt = point_to_proc+psb_elem_recv_
          rcv_pt = brvidx(proc_to_comm)
          call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
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
      do while (proc_to_comm .ne. -1)
        nerv = d_idx(point_to_proc+psb_n_elem_recv_)
        nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

        if(proc_to_comm .eq. myrow) then
          idx_pt = point_to_proc+psb_elem_recv_
          snd_pt = bsdidx(proc_to_comm)
          call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
               & sndbuf(snd_pt:snd_pt+nesd-1),beta,y)
        end if

        point_to_proc = point_to_proc+nerv+nesd+3
        proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      end do

    else

      point_to_proc = 1
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      do while (proc_to_comm .ne. -1)
        nerv = d_idx(point_to_proc+psb_n_elem_recv_)
        nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

        if(proc_to_comm.ne.myrow) then
          call mpi_wait(rvhd(proc_to_comm),p2pstat,iret)
          if(iret.ne.mpi_success) then
            int_err(1) = iret
            info=400
            call psb_errpush(info,name,i_err=int_err)
            goto 9999
          end if
          idx_pt = point_to_proc+psb_elem_recv_
          rcv_pt = brvidx(proc_to_comm)
          call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
               & rcvbuf(rcv_pt:rcv_pt+n*nerv-1),beta,y)
        else
          idx_pt = point_to_proc+psb_elem_recv_
          snd_pt = bsdidx(proc_to_comm)
          call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
               & sndbuf(snd_pt:snd_pt+n*nesd-1),beta,y)

        end if

        point_to_proc = point_to_proc+nerv+nesd+3
        proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
      end do

    end if


  else if (swap_send) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      idx_pt = point_to_proc+nerv+psb_elem_send_
      snd_pt = bsdidx(proc_to_comm)
      call psi_gth(nesd,d_idx(idx_pt:idx_pt+nesd-1),&
           & y,sndbuf(snd_pt:snd_pt+nesd-1))
      call dgesd2d(icontxt,nesd,1,sndbuf(snd_pt),nesd,proc_to_comm,0)

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  else if (swap_recv) then

    point_to_proc = 1
    proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    do while (proc_to_comm .ne. -1)
      nerv = d_idx(point_to_proc+psb_n_elem_recv_)
      nesd = d_idx(point_to_proc+nerv+psb_n_elem_send_)

      if(proc_to_comm.ne.myrow) then
        rcv_pt = brvidx(proc_to_comm)
        call dgerv2d(icontxt,nerv,1,rcvbuf(rcv_pt),nerv,proc_to_comm,0)
        idx_pt = point_to_proc+psb_elem_recv_
        rcv_pt = brvidx(proc_to_comm)
        call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
             & rcvbuf(rcv_pt:rcv_pt+nerv-1),beta,y)
      else
        idx_pt = point_to_proc+psb_elem_recv_
        snd_pt = bsdidx(proc_to_comm)
        call psi_sct(nerv,d_idx(idx_pt:idx_pt+nerv-1),&
             & sndbuf(snd_pt:snd_pt+nesd-1),beta,y)
      end if

      point_to_proc = point_to_proc+nerv+nesd+3
      proc_to_comm  = d_idx(point_to_proc+psb_proc_id_)
    end do

  end if

  deallocate(sdsz,rvsz,bsdidx,&
       & brvidx,rvhd,prcid,&
       & ptp,stat=info)
  if(all) deallocate(sndbuf,rcvbuf,stat=info)
  if(info.ne.0) then
    call psb_errpush(4000,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error(icontxt)
    return
  end if
  return
end subroutine psi_dswapdatav
