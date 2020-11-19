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
!
! File: psi_adjcncy_fnd_owner.f90
!
! Subroutine: psi_adjcncy_fnd_owner
!   Figure out who owns  global indices. 
! 
! Arguments: 
!    idx(:)   - integer                 Required indices on the calling process.
!                                       Note: the indices should be unique!
!    iprc(:)  - integer(psb_ipk_), allocatable    Output: process identifiers for the corresponding
!                                       indices
!    adj(:)   - integer(psb_ipk_)       Input: list of topological neighbours for current process.
!
!    idxmap   - class(psb_indx_map).    The index map
!    info     - integer.                return code.
!
! This version takes on input a list of processes that are assumed to 
! be topological neighbours of the current one. Each process will send to all 
! of its neighbours the list of indices for which it is trying to find the
! owner, prepare its own answers, and collect answers from others.
! There are three possibile implementations: using mpi_alltoallv, using mpi_isend/irecv,
! using psb_snd/psb_rcv. The default is mpi_alltoallv.
! 
subroutine psi_adjcncy_fnd_owner(idx,iprc,adj,idxmap,info)
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psb_timers_mod
  use psb_indx_map_mod, psb_protect_name => psi_adjcncy_fnd_owner
#ifdef MPI_MOD
  use mpi
#endif

  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif
  integer(psb_lpk_), intent(in)   :: idx(:)
  integer(psb_ipk_), allocatable, intent(out) :: iprc(:)
  integer(psb_ipk_), intent(in)   :: adj(:)
  class(psb_indx_map), intent(in) :: idxmap
  integer(psb_ipk_), intent(out)  :: info


  integer(psb_lpk_), allocatable :: rmtidx(:)
  integer(psb_ipk_), allocatable :: tproc(:), lclidx(:)
  integer(psb_mpk_), allocatable :: hsz(:),hidx(:), sdidx(:), rvidx(:),&
       & sdsz(:), rvsz(:), sdhd(:), rvhd(:), p2pstat(:,:)
  integer(psb_mpk_) :: prc, p2ptag, iret
  integer(psb_mpk_) :: icomm, minfo
  integer(psb_ipk_) :: i,n_row,n_col,err_act,hsize,ip,isz,j, k,&
       & last_ih, last_j, nidx, nrecv, nadj
  integer(psb_lpk_) :: mglob, ih
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np,me
  logical, parameter  :: gettime=.true., debug=.false.
  integer(psb_mpk_)   :: xchg_alg 
  logical, parameter  :: do_timings=.false.
  integer(psb_ipk_), save  :: idx_phase1=-1, idx_phase2=-1, idx_phase3=-1
  integer(psb_ipk_), save  :: idx_phase11=-1, idx_phase12=-1, idx_phase13=-1
  real(psb_dpk_)      :: t0, t1, t2, t3, t4, tamx, tidx
  character(len=20)   :: name

  info = psb_success_
  name = 'psi_adjcncy_fnd_owner'
  call psb_erractionsave(err_act)

  ctxt   = idxmap%get_ctxt()
  icomm   = idxmap%get_mpic()
  mglob   = idxmap%get_gr()
  n_row   = idxmap%get_lr()
  n_col   = idxmap%get_lc()

  if ((do_timings).and.(idx_phase1==-1))       &
       & idx_phase1 = psb_get_timer_idx("ADJ_FND_OWN: phase1 ")
  if ((do_timings).and.(idx_phase2==-1))       &
       & idx_phase2 = psb_get_timer_idx("ADJ_FND_OWN: phase2")
  if ((do_timings).and.(idx_phase3==-1))       &
       & idx_phase3 = psb_get_timer_idx("ADJ_FND_OWN: phase3")
  if ((do_timings).and.(idx_phase11==-1))       &
       & idx_phase11 = psb_get_timer_idx("ADJ_FND_OWN: phase11 ")
  if ((do_timings).and.(idx_phase12==-1))       &
       & idx_phase12 = psb_get_timer_idx("ADJ_FND_OWN: phase12")
  if ((do_timings).and.(idx_phase13==-1))       &
       & idx_phase13 = psb_get_timer_idx("ADJ_FND_OWN: phase13")


  call psb_info(ctxt, me, np)

  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.(idxmap%is_valid())) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='invalid idxmap')
    goto 9999      
  end if

  if (gettime) then 
    t0 = psb_wtime()
  end if

  nadj = size(adj)
  nidx = size(idx)
  call psb_realloc(nidx,iprc,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_realloc')
    goto 9999      
  end if
  iprc = -1
  ! write(0,*) me,name,' Going through ',nidx,nadj
  xchg_alg = psi_get_adj_alg()
  select case(xchg_alg)
  case(psi_adj_fnd_a2av_) 
    if (do_timings) call psb_tic(idx_phase1)

    !
    ! First simple minded version with auxiliary arrays
    ! dimensioned on NP.
    ! Do the exchange with an alltoallv
    ! 
    !    
    Allocate(hidx(0:np),hsz(np),sdsz(0:np-1),rvsz(0:np-1), &
         & sdidx(0:np),rvidx(0:np),stat=info)
    !
    ! Same send buffer for everybody
    !
    sdidx(:) = 0
    !
    ! First, send sizes according to adjcncy list
    !
    if (do_timings) call psb_tic(idx_phase11)
    sdsz = 0 
    do j=1, nadj
      sdsz(adj(j)) = nidx
    end do
    !write(0,*)me,' Check on sizes into a2a:',adj(:),nadj,':',sdsz(:)

    call mpi_alltoall(sdsz,1,psb_mpi_mpk_,&
         & rvsz,1,psb_mpi_mpk_,icomm,minfo)
    if (do_timings) call psb_toc(idx_phase11)
    if (do_timings) call psb_tic(idx_phase12)
    rvidx(0) = 0
    do i=0, np-1
      rvidx(i+1) = rvidx(i) + rvsz(i)
    end do
    hsize = rvidx(np)
    
    ! write(0,*)me,' Check on sizes from a2a:',hsize,rvsz(:)
    !
    ! Second, allocate buffers and exchange data
    !
    if (do_timings) call psb_toc(idx_phase12)
    if (do_timings) call psb_tic(idx_phase13)
    Allocate(rmtidx(hsize),lclidx(max(hsize,nidx*nadj)),&
         & tproc(max(hsize,nidx)),stat=info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

    call mpi_alltoallv(idx,sdsz,sdidx,psb_mpi_lpk_,&
         & rmtidx,rvsz,rvidx,psb_mpi_lpk_,icomm,iret)
    if (do_timings) call psb_toc(idx_phase13)
    if (do_timings) call psb_toc(idx_phase1)
    if (do_timings) call psb_tic(idx_phase2)
    !
    ! Third, compute local answers
    !
    call idxmap%g2l(rmtidx(1:hsize),lclidx(1:hsize),info,owned=.true.)
    do i=1, hsize
      tproc(i) = -1
      if ((0 < lclidx(i)).and. (lclidx(i) <= n_row)) tproc(i) = me
    end do
    if (do_timings) call psb_toc(idx_phase2)
    if (do_timings) call psb_tic(idx_phase3)

    !
    ! Fourth, exchange the answers 
    !
    ! Adjust sdidx for reuse in receiving lclidx array 
    do i=0,np-1
      sdidx(i+1) = sdidx(i) + sdsz(i)
    end do
    call mpi_alltoallv(tproc,rvsz,rvidx,psb_mpi_ipk_,&
         & lclidx,sdsz,sdidx,psb_mpi_ipk_,icomm,iret)

    !
    ! Because IPRC has been initialized to -1, the MAX operation selects
    ! the answers. 
    !
    do i=0, np-1
      if (sdsz(i)>0) then
        ! Must be nidx == sdsz(i) 
        iprc(1:nidx) = max(iprc(1:nidx), lclidx(sdidx(i)+1:sdidx(i)+sdsz(i)))
      end if
    end do
    if (do_timings) call psb_toc(idx_phase3)

    if (debug) write(0,*) me,' End of adjcncy_fnd ',iprc(1:nidx)    

  case(psi_adj_fnd_irecv_)
    
    if (do_timings) call psb_tic(idx_phase1)
    !
    ! First simple minded version with auxiliary arrays
    ! dimensioned on NP.
    ! Could it be improved with a loop based on the maximum length
    ! of adj(:) ???
    !    
    Allocate(hidx(0:np),hsz(np),sdsz(0:np-1),rvsz(0:np-1),&
         & sdhd(0:np-1), rvhd(0:np-1), p2pstat(mpi_status_size,0:np-1),&
         & stat=info)
    if (do_timings) call psb_tic(idx_phase11)
    sdhd(:) = mpi_request_null
    rvhd(:) = mpi_request_null
    !
    ! First, send sizes according to adjcncy list
    !
    sdsz = 0 
    do j=1, nadj
      sdsz(adj(j)) = nidx
    end do
    !write(0,*)me,' Check on sizes into a2a:',adj(:),nadj,':',sdsz(:)

    call mpi_alltoall(sdsz,1,psb_mpi_mpk_,&
         & rvsz,1,psb_mpi_mpk_,icomm,minfo)
    hidx(0) = 0
    do i=0, np-1
      hidx(i+1) = hidx(i) + rvsz(i)
    end do
    hsize = hidx(np)
    ! write(0,*)me,' Check on sizes from a2a:',hsize,rvsz(:)
    !
    ! Second, allocate buffers and exchange data
    !
    Allocate(rmtidx(hsize),lclidx(max(hsize,nidx*nadj)),tproc(max(hsize,nidx)),stat=info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if
    do i = 0, np-1
      if (rvsz(i)>0) then
        ! write(0,*) me, ' First receive from ',i,rvsz(i)
        prc = psb_get_mpi_rank(ctxt,i)        
        p2ptag = psb_long_swap_tag
        !write(0,*) me, ' Posting first receive from ',i,rvsz(i),prc           
        call mpi_irecv(rmtidx(hidx(i)+1),rvsz(i),&
             & psb_mpi_lpk_,prc,&
             & p2ptag, icomm,rvhd(i),iret)
      end if
    end do
    if (do_timings) call psb_toc(idx_phase11)
    if (do_timings) call psb_tic(idx_phase12)    
    do j=1, nadj
      if (nidx > 0) then
        prc = psb_get_mpi_rank(ctxt,adj(j))        
        p2ptag = psb_long_swap_tag
        !write(0,*) me, ' First send to ',adj(j),nidx, prc
        call mpi_send(idx,nidx,&
             & psb_mpi_lpk_,prc,&
             & p2ptag, icomm,iret)
      end if
    end do
    if (do_timings) call psb_toc(idx_phase12)
    if (do_timings) call psb_tic(idx_phase13)
    call mpi_waitall(np,rvhd,p2pstat,iret)
    if (do_timings) call psb_toc(idx_phase13)
    if (do_timings) call psb_toc(idx_phase1)
    if (do_timings) call psb_tic(idx_phase2)

    !
    ! Third, compute local answers
    !
    call idxmap%g2l(rmtidx(1:hsize),lclidx(1:hsize),info,owned=.true.)
    do i=1, hsize
      tproc(i) = -1
      if ((0 < lclidx(i)).and. (lclidx(i) <= n_row)) tproc(i) = me
    end do
    if (do_timings) call psb_toc(idx_phase2)
    if (do_timings) call psb_tic(idx_phase3)
    !
    ! At this point we can reuse lclidx to receive messages
    !
    rvhd(:) = mpi_request_null
    do j=1, nadj
      !write(0,*) me, ' First send to ',adj(j),nidx
      if (nidx > 0) then
        prc = psb_get_mpi_rank(ctxt,adj(j))        
        p2ptag = psb_int_swap_tag
        !write(0,*) me, ' Posting second receive from ',adj(j),nidx, prc
        call mpi_irecv(lclidx((j-1)*nidx+1),nidx, &
             & psb_mpi_ipk_,prc,&
             & p2ptag, icomm,rvhd(j),iret)
      end if
    end do

    !
    ! Fourth, send data back; 
    !
    do i = 0, np-1
      if (rvsz(i)>0) then
        prc = psb_get_mpi_rank(ctxt,i)        
        p2ptag = psb_int_swap_tag
        !write(0,*) me, ' Second send to ',i,rvsz(i), prc
        call mpi_send(tproc(hidx(i)+1),rvsz(i),&
             & psb_mpi_ipk_,prc,&
             & p2ptag, icomm,iret)
      end if
    end do
    !
    ! Fifth: receive and combine. MAX works because default
    ! answer is -1. 
    !
    call mpi_waitall(np,rvhd,p2pstat,iret)    
    do j = 1, nadj
      iprc(1:nidx) = max(iprc(1:nidx), lclidx((j-1)*nidx+1:(j-1)*nidx+nidx))
    end do
    if (do_timings) call psb_toc(idx_phase3)
    if (debug) write(0,*) me,' End of adjcncy_fnd ',iprc(1:nidx)

  case(psi_adj_fnd_pbrcv_) 

    Allocate(hidx(0:np),hsz(np),&
         & sdsz(0:np-1),rvsz(0:np-1),stat=info)
    !
    ! First, send sizes according to adjcncy list
    !
    sdsz = 0 
    do j=1, nadj
      sdsz(adj(j)) = nidx
    end do
    !write(0,*)me,' Check on sizes into a2a:',adj(:),nadj,':',sdsz(:)

    call mpi_alltoall(sdsz,1,psb_mpi_mpk_,&
         & rvsz,1,psb_mpi_mpk_,icomm,minfo)
    hidx(0) = 0
    do i=0, np-1
      hidx(i+1) = hidx(i) + rvsz(i)
    end do
    hsize = hidx(np)
    ! write(0,*)me,' Check on sizes from a2a:',hsize,rvsz(:)
    !
    ! Second, allocate buffers and exchange data
    !
    Allocate(rmtidx(hsize),lclidx(hsize),tproc(max(hsize,nidx)),stat=info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if
    do j=1, nadj
      !write(0,*) me, ' First send to ',adj(j),nidx
      if (nidx > 0) call psb_snd(ctxt,idx(1:nidx),adj(j))
    end do
    do i = 0, np-1
      if (rvsz(i)>0) then
        ! write(0,*) me, ' First receive from ',i,rvsz(i)           
        call psb_rcv(ctxt,rmtidx(hidx(i)+1:hidx(i)+rvsz(i)),i)
      end if
    end do

    !
    ! Third, compute local answers
    !
    call idxmap%g2l(rmtidx(1:hsize),lclidx(1:hsize),info,owned=.true.)
    do i=1, hsize
      tproc(i) = -1
      if ((0 < lclidx(i)).and. (lclidx(i) <= n_row)) tproc(i) = me
    end do

    !
    ! Fourth, send data back; 
    !
    do i = 0, np-1
      if (rvsz(i)>0) then
        !write(0,*) me, ' Second send to ',i,rvsz(i)
        call psb_snd(ctxt,tproc(hidx(i)+1:hidx(i)+rvsz(i)),i)
      end if
    end do
    !
    ! Fifth: receive and combine. MAX works because default
    ! answer is -1. Reuse tproc
    !
    do j = 1, nadj
      !write(0,*) me, ' Second receive from ',adj(j), nidx          
      if (nidx > 0) call psb_rcv(ctxt,tproc(1:nidx),adj(j))
      iprc(1:nidx) = max(iprc(1:nidx), tproc(1:nidx))
    end do
  case default 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invalid exchange alg choice')
    goto 9999      
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psi_adjcncy_fnd_owner
