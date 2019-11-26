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
! File: psi_fnd_owner.f90
!
! Subroutine: psi_fnd_owner
!   Figure out who owns  global indices. 
! 
! Arguments: 
! 
subroutine psi_adjcncy_fnd_owner(idx,iprc,adj,idxmap,info)
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
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
  integer(psb_mpk_) :: icomm, minfo, iictxt
  integer(psb_ipk_) :: i,n_row,n_col,err_act,hsize,ip,isz,j, k,&
       & last_ih, last_j, nidx, nrecv, nadj
  integer(psb_lpk_) :: mglob, ih
  integer(psb_ipk_) :: ictxt,np,me
  logical, parameter  :: gettime=.false., new_impl=.true.
  logical, parameter  :: a2av_impl=.true., debug=.false.
  real(psb_dpk_)      :: t0, t1, t2, t3, t4, tamx, tidx
  character(len=20)   :: name

  info = psb_success_
  name = 'psi_adjcncy_fnd_owner'
  call psb_erractionsave(err_act)

  ictxt   = idxmap%get_ctxt()
  icomm   = idxmap%get_mpic()
  mglob   = idxmap%get_gr()
  n_row   = idxmap%get_lr()
  n_col   = idxmap%get_lc()
  iictxt = ictxt 

  call psb_info(ictxt, me, np)

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

  if (a2av_impl) then
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
    sdsz = 0 
    do j=1, nadj
      sdsz(adj(j)) = nidx
    end do
    !write(0,*)me,' Check on sizes into a2a:',adj(:),nadj,':',sdsz(:)

    call mpi_alltoall(sdsz,1,psb_mpi_mpk_,&
         & rvsz,1,psb_mpi_mpk_,icomm,minfo)
    rvidx(0) = 0
    do i=0, np-1
      rvidx(i+1) = rvidx(i) + rvsz(i)
    end do
    hsize = rvidx(np)
    ! write(0,*)me,' Check on sizes from a2a:',hsize,rvsz(:)
    !
    ! Second, allocate buffers and exchange data
    !
    Allocate(rmtidx(hsize),lclidx(max(hsize,nidx*nadj)),&
         & tproc(max(hsize,nidx)),stat=info)

    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

    call mpi_alltoallv(idx,sdsz,sdidx,psb_mpi_lpk_,&
         & rmtidx,rvsz,rvidx,psb_mpi_lpk_,icomm,iret)

    !
    ! Third, compute local answers
    !
    call idxmap%g2l(rmtidx(1:hsize),lclidx(1:hsize),info,owned=.true.)
    do i=1, hsize
      tproc(i) = -1
      if ((0 < lclidx(i)).and. (lclidx(i) <= n_row)) tproc(i) = me
    end do

    !
    ! Fourth, exchange the answers 
    !
    ! Adjust sdidx for reuse in receiving lclidx array 
    do i=0,np-1
      sdidx(i+1) = sdidx(i) + sdsz(i)
    end do
    call mpi_alltoallv(tproc,rvsz,rvidx,psb_mpi_ipk_,&
         & lclidx,sdsz,sdidx,psb_mpi_ipk_,icomm,iret)

    do i=0, np-1
      if (sdsz(i)>0) then
        ! Must be nidx == sdsz(i) 
        iprc(1:nidx) = max(iprc(1:nidx), lclidx(sdidx(i)+1:sdidx(i)+sdsz(i)))
      end if
    end do

    if (debug) write(0,*) me,' End of adjcncy_fnd ',iprc(1:nidx)    
  else
    if (new_impl) then
      !
      ! First simple minded version with auxiliary arrays
      ! dimensioned on NP.
      ! Could it be improved with a loop based on the maximum length
      ! of adj(:) ???
      !    
      Allocate(hidx(0:np),hsz(np),sdsz(0:np-1),rvsz(0:np-1),&
           & sdhd(0:np-1), rvhd(0:np-1), p2pstat(mpi_status_size,0:np-1),&
           & stat=info)
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
          prc = psb_get_rank(ictxt,i)        
          p2ptag = psb_long_swap_tag
          !write(0,*) me, ' Posting first receive from ',i,rvsz(i),prc           
          call mpi_irecv(rmtidx(hidx(i)+1),rvsz(i),&
               & psb_mpi_lpk_,prc,&
               & p2ptag, icomm,rvhd(i),iret)
        end if
      end do
      do j=1, nadj
        if (nidx > 0) then
          !call psb_snd(ictxt,idx(1:nidx),adj(j))
          prc = psb_get_rank(ictxt,adj(j))        
          p2ptag = psb_long_swap_tag
          !write(0,*) me, ' First send to ',adj(j),nidx, prc
          call mpi_send(idx,nidx,&
               & psb_mpi_lpk_,prc,&
               & p2ptag, icomm,iret)
        end if
      end do
!!$    do i = 0, np-1
!!$      if (rvsz(i)>0) then
!!$        ! write(0,*) me, ' First receive from ',i,rvsz(i)           
!!$        call psb_rcv(ictxt,rmtidx(hidx(i)+1:hidx(i)+rvsz(i)),i)
!!$      end if
!!$    end do
      call mpi_waitall(np,rvhd,p2pstat,iret)

      !
      ! Third, compute local answers
      !
      call idxmap%g2l(rmtidx(1:hsize),lclidx(1:hsize),info,owned=.true.)
      do i=1, hsize
        tproc(i) = -1
        if ((0 < lclidx(i)).and. (lclidx(i) <= n_row)) tproc(i) = me
      end do
      !
      ! At this point we can reuse lclidx to receive messages
      !
      rvhd(:) = mpi_request_null
      do j=1, nadj
        !write(0,*) me, ' First send to ',adj(j),nidx
        if (nidx > 0) then
          !call psb_snd(ictxt,idx(1:nidx),adj(j))
          prc = psb_get_rank(ictxt,adj(j))        
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
          !call psb_snd(ictxt,tproc(hidx(i)+1:hidx(i)+rvsz(i)),i)
          prc = psb_get_rank(ictxt,i)        
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
        !write(0,*) me, ' Second receive from ',adj(j), nidx          
        !if (nidx > 0) call psb_rcv(ictxt,tproc(1:nidx),adj(j))
        iprc(1:nidx) = max(iprc(1:nidx), lclidx((j-1)*nidx+1:(j-1)*nidx+nidx))
      end do
      if (debug) write(0,*) me,' End of adjcncy_fnd ',iprc(1:nidx)

    else

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
        if (nidx > 0) call psb_snd(ictxt,idx(1:nidx),adj(j))
      end do
      do i = 0, np-1
        if (rvsz(i)>0) then
          ! write(0,*) me, ' First receive from ',i,rvsz(i)           
          call psb_rcv(ictxt,rmtidx(hidx(i)+1:hidx(i)+rvsz(i)),i)
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
          call psb_snd(ictxt,tproc(hidx(i)+1:hidx(i)+rvsz(i)),i)
        end if
      end do
      !
      ! Fifth: receive and combine. MAX works because default
      ! answer is -1. Reuse tproc
      !
      do j = 1, nadj
        !write(0,*) me, ' Second receive from ',adj(j), nidx          
        if (nidx > 0) call psb_rcv(ictxt,tproc(1:nidx),adj(j))
        iprc(1:nidx) = max(iprc(1:nidx), tproc(1:nidx))
      end do
    end if
  end if
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psi_adjcncy_fnd_owner
