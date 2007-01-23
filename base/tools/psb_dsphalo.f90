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
! File: psb_dsphalo.f90
!
!*****************************************************************************
!*                                                                           *
!*  This routine does the retrieval of remote matrix rows.                   *
!*  Note that retrieval is done through GTBLK, therefore it should work      *
!*  for any format.                                                          *
!*  Currently the output is BLK%FIDA='CSR' but it would take little          *
!*  work to change that; the pieces are transferred in COO format            *
!*  thus we would only need a DCSDP at the end to exit in whatever format    *
!*  is needed.                                                               *
!*  But I'm feeling soooooo lazy today......                                 *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
Subroutine psb_dsphalo(a,desc_a,blk,info,rwcnv,clcnv,cliprow,outfmt,data)

  use psb_const_mod
  use psb_serial_mod
  use psb_descriptor_type
  use psb_realloc_mod
  use psb_tools_mod, only : psb_glob_to_loc, psb_loc_to_glob
  use psb_error_mod
  use psb_penv_mod
  use mpi
  Implicit None

  Type(psb_dspmat_type),Intent(in)    :: a
  Type(psb_dspmat_type),Intent(inout) :: blk
  Type(psb_desc_type),Intent(in), target  :: desc_a
  integer, intent(out)                :: info
  logical, optional, intent(in)       :: rwcnv,clcnv,cliprow
  character(len=5), optional          :: outfmt 
  integer, intent(in), optional       :: data
  !     ...local scalars....
  Integer    :: np,me,counter,proc,i, &
       &     n_el_send,k,n_el_recv,ictxt, idx, r, tot_elem,&
       &     n_elem, j, ipx,mat_recv, iszs, iszr,idxs,idxr,nz,nrmin,&
       &     data_
  Type(psb_dspmat_type)     :: tmp
  Integer :: l1, icomm, err_act
  Integer, allocatable  :: sdid(:,:), brvindx(:),rvid(:,:), &
       & rvsz(:), bsdindx(:),sdsz(:)
  integer, pointer :: idxv(:)
  logical :: rwcnv_,clcnv_,cliprow_
  character(len=5)  :: outfmt_
  Logical,Parameter :: debug=.false., debugprt=.false.
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6,t7,t8,t9
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=0
  name='psb_dsphalo'
  call psb_erractionsave(err_act)

  if(debug) write(0,*)'Inside DSPHALO'
  if (present(rwcnv)) then 
    rwcnv_ = rwcnv
  else
    rwcnv_ = .true.
  endif
  if (present(clcnv)) then 
    clcnv_ = clcnv
  else
    clcnv_ = .true.
  endif
  if (present(cliprow)) then 
    cliprow_ = cliprow
  else
    cliprow_ = .false.
  endif
  if (present(data)) then 
    data_ = data
  else
    data_ = psb_comm_halo_
  endif

  if (present(outfmt)) then 
    outfmt_ =  toupper(outfmt)
  else
    outfmt_ = 'CSR'
  endif

  ictxt = psb_cd_get_context(desc_a)

  Call psb_info(ictxt, me, np)

  t1 = mpi_wtime()
  Allocate(sdid(np,3),rvid(np,3),brvindx(np+1),&
       & rvsz(np),sdsz(np),bsdindx(np+1),stat=info)

  if (info /= 0) then
    info=4000
    call psb_errpush(info,name)
    goto 9999
  end if

  If (debug) Write(0,*)'dsphalo',me

  select case(data_) 
  case(psb_comm_halo_) 
    idxv => desc_a%halo_index

  case(psb_comm_ovr_) 
    idxv => desc_a%ovrlap_index

  case(psb_comm_ext_) 
    idxv => desc_a%ext_index

  case default
    call psb_errpush(4010,name,a_err='wrong Data selector')
    goto 9999
  end select


  l1  = 0

  sdsz(:)=0
  rvsz(:)=0
  ipx = 1
  brvindx(ipx) = 0
  bsdindx(ipx) = 0
  counter=1
  idx = 0
  idxs = 0
  idxr = 0
  blk%k = a%k
  blk%m = 0 
  ! For all rows in the halo descriptor, extract and send/receive.
  Do 
    proc=idxv(counter)
    if (proc == -1) exit
    n_el_recv = idxv(counter+psb_n_elem_recv_)
    counter   = counter+n_el_recv
    n_el_send = idxv(counter+psb_n_elem_send_)
    tot_elem = 0
    Do j=0,n_el_send-1
      idx = idxv(counter+psb_elem_send_+j)
      n_elem = psb_sp_get_nnz_row(idx,a)
      tot_elem = tot_elem+n_elem      
    Enddo
    sdsz(proc+1) = tot_elem

    blk%m = blk%m + n_el_recv

    counter   = counter+n_el_send+3
  Enddo
  call psb_get_mpicomm(ictxt,icomm)

  call mpi_alltoall(sdsz,1,mpi_integer,rvsz,1,mpi_integer,icomm,info)
  if (info /= 0) then
    info=4010
    ch_err='mpi_alltoall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  idxs = 0
  idxr = 0
  counter = 1
  Do 
    proc=idxv(counter)
    if (proc == -1) exit
    n_el_recv = idxv(counter+psb_n_elem_recv_)
    counter   = counter+n_el_recv
    n_el_send = idxv(counter+psb_n_elem_send_)

    bsdindx(proc+1) = idxs
    idxs = idxs + sdsz(proc+1)
    brvindx(proc+1) = idxr
    idxr = idxr + rvsz(proc+1)
    counter   = counter+n_el_send+3
  Enddo

  iszr=sum(rvsz)
  call psb_sp_reall(blk,max(iszr,1),info)
  if(debug)  write(0,*)me,'SPHALO Sizes:',size(blk%ia1),size(blk%ia2)
  if (info /= 0) then
    info=4010
    ch_err='psb_sp_reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  mat_recv = iszr
  iszs=sum(sdsz)
  call psb_nullify_sp(tmp)
  call psb_sp_all(0,0,tmp,max(iszs,1),info)
  tmp%fida='COO'
  call psb_sp_setifld(psb_spmat_asb_,psb_state_,tmp,info)

  t2 = mpi_wtime()

  l1  = 0
  ipx = 1
  counter=1
  idx = 0

  Do 
    proc=idxv(counter)
    if (proc == -1) exit 
    n_el_recv=idxv(counter+psb_n_elem_recv_)
    counter=counter+n_el_recv
    n_el_send=idxv(counter+psb_n_elem_send_)
    tot_elem=0

    Do j=0,n_el_send-1
      idx = idxv(counter+psb_elem_send_+j)
      n_elem = psb_sp_get_nnz_row(idx,a)

      call psb_sp_getblk(idx,a,tmp,info,append=.true.)
      if (info /= 0) then
        info=4010
        ch_err='psb_sp_getblk'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      tot_elem=tot_elem+n_elem
    Enddo

    ipx = ipx + 1 

    counter   = counter+n_el_send+3
  Enddo
  nz = tmp%infoa(psb_nnz_)

  if (rwcnv_) call psb_loc_to_glob(tmp%ia1(1:nz),desc_a,info,iact='I')
  if (clcnv_) call psb_loc_to_glob(tmp%ia2(1:nz),desc_a,info,iact='I')
  if (info /= 0) then
    info=4010
    ch_err='psb_loc_to_glob'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (debugprt) then 
    open(30+me)
    call psb_csprt(30+me,tmp,head='% SPHALO border SEND .')
    call flush(30+me)
    close(30+me)
  end if


  call mpi_alltoallv(tmp%aspk,sdsz,bsdindx,mpi_double_precision,&
       & blk%aspk,rvsz,brvindx,mpi_double_precision,icomm,info)
  call mpi_alltoallv(tmp%ia1,sdsz,bsdindx,mpi_integer,&
       & blk%ia1,rvsz,brvindx,mpi_integer,icomm,info)
  call mpi_alltoallv(tmp%ia2,sdsz,bsdindx,mpi_integer,&
       & blk%ia2,rvsz,brvindx,mpi_integer,icomm,info)
  if (info /= 0) then
    info=4010
    ch_err='mpi_alltoallv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  t3 = mpi_wtime()


  !
  ! Convert into local numbering 
  !
  if (rwcnv_) call psb_glob_to_loc(blk%ia1(1:iszr),desc_a,info,iact='I',owned=cliprow_)
  if (clcnv_) call psb_glob_to_loc(blk%ia2(1:iszr),desc_a,info,iact='I')

  if (info /= 0) then
    info=4010
    ch_err='psbglob_to_loc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (debugprt) then 
    blk%fida='COO'
    blk%infoa(psb_nnz_)=iszr
    open(40+me)
    call psb_csprt(40+me,blk,head='% SPHALO border .')
    call flush(40+me)
    close(40+me)
  end if
  l1  = 0
  blk%m=0
  nrmin=max(0,a%m)
  Do i=1,iszr
!!$      write(0,*) work5(i),work6(i)
    r=(blk%ia1(i))
    k=(blk%ia2(i))
    If ((r>nrmin).and.(k>0)) Then
      l1=l1+1
      blk%aspk(l1) = blk%aspk(i)
      blk%ia1(l1) = r
      blk%ia2(l1) = k
      blk%k = max(blk%k,k)
      blk%m = max(blk%m,r)
    End If
  Enddo
  blk%fida='COO'
  blk%infoa(psb_nnz_)=l1
  blk%m = blk%m - a%m
  if (debugprt) then 
    open(50+me)
    call psb_csprt(50+me,blk,head='% SPHALO border .')
    call flush(50+me)
    close(50+me)
    call psb_barrier(ictxt)
  end if
  t4 = mpi_wtime()

  if(debug) Write(0,*)me,'End first loop',counter,l1,blk%m

  !
  ! Combined sort & conversion to CSR. 
  !
  if(debug) write(0,*) me,'Calling ipcoo2csr from dsphalo ',blk%m,blk%k,l1,blk%ia2(2)

  select case(outfmt_)
  case ('CSR') 
    call psb_ipcoo2csr(blk,info,rwshr=.true.)
    if (info /= 0) then
      info=4010
      ch_err='psb_ipcoo2csr'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  case('COO')
    ! Do nothing! 
  case default
    write(0,*) 'Error in DSPHALO : invalid outfmt "',outfmt_,'"'
    info=4010
    ch_err='Bad outfmt'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end select
  t5 = mpi_wtime()



!!$  write(0,'(i3,1x,a,4(1x,i14))') me,'DSPHALO sizes:',iszr,iszs
!!$  write(0,'(i3,1x,a,4(1x,g14.5))') me,'DSPHALO timings:',t6-t2,t7-t6,t8-t7,t3-t8
!!$  write(0,'(i3,1x,a,4(1x,g14.5))') me,'DSPHALO timings:',t2-t1,t3-t2,t4-t3,t5-t4

  Deallocate(sdid,brvindx,rvid,bsdindx,rvsz,sdsz,stat=info)

  call psb_sp_free(tmp,info)
  if (info /= 0) then
    info=4010
    ch_err='psb_sp_free'
    call psb_errpush(info,name,a_err=ch_err)
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

End Subroutine psb_dsphalo
