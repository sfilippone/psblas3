!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
! File: psb_zsphalo.f90
!
! Subroutine: psb_zsphalo
!  This routine does the retrieval of remote matrix rows.                   
!  Note that retrieval is done through GTBLK, therefore it should work      
!  for any matrix format in A; as for the output, default is CSR. 
!    
! 
! Arguments: 
!    a        - type(psb_zspmat_type)   The local part of input matrix A
!    desc_a   - type(psb_desc_type).  The communication descriptor.
!    blck     - type(psb_zspmat_type)   The local part of output matrix BLCK
!    info     - integer.                Return code
!    rowcnv   - logical                 Should row/col indices be converted
!    colcnv   - logical                 to/from global numbering when sent/received?
!                                       default is .TRUE.
!    rowscale - logical                 Should row/col indices on output be remapped
!    colscale - logical                 from MIN:MAX  to 1:(MAX-MIN+1) ? 
!                                       default is .FALSE. 
!                                       (commmon use is ROWSCALE=.TRUE., COLSCALE=.FALSE.)
!    data     - integer                 Which index list in desc_a should be used to retrieve
!                                       rows, default psb_comm_halo_
!                                       psb_comm_halo_    use halo_index
!                                       psb_comm_ext_     use ext_index 
!                                       psb_comm_ovrl_  DISABLED for this routine.
!
Subroutine psb_zsphalo(a,desc_a,blk,info,rowcnv,colcnv,&
     &  rowscale,colscale,outfmt,data)

  use psb_const_mod
  use psb_serial_mod
  use psb_descriptor_type
  use psb_realloc_mod
  use psb_tools_mod, psb_protect_name => psb_zsphalo
  use psb_error_mod
  use psb_penv_mod
#ifdef MPI_MOD
    use mpi
#endif
  Implicit None
#ifdef MPI_H
    include 'mpif.h'
#endif

  Type(psb_zspmat_type),Intent(in)    :: a
  Type(psb_zspmat_type),Intent(inout) :: blk
  Type(psb_desc_type),Intent(in), target :: desc_a
  integer, intent(out)                :: info
  logical, optional, intent(in)       :: rowcnv,colcnv,rowscale,colscale
  character(len=5), optional          :: outfmt 
  integer, intent(in), optional       :: data
  !     ...local scalars....
  Integer    :: np,me,counter,proc,i, &
       &     n_el_send,k,n_el_recv,ictxt, idx, r, tot_elem,&
       &     n_elem, j, ipx,mat_recv, iszs, iszr,idxs,idxr,nz,&
       &     irmin,icmin,irmax,icmax,data_,ngtz
  Integer :: l1, icomm, err_act
  Integer, allocatable  :: sdid(:,:), brvindx(:),rvid(:,:), &
       & rvsz(:), bsdindx(:),sdsz(:), iasnd(:), jasnd(:)
  complex(psb_dpk_), allocatable :: valsnd(:)
  integer, pointer :: idxv(:)
  logical :: rowcnv_,colcnv_,rowscale_,colscale_
  character(len=5)  :: outfmt_
  integer           :: debug_level, debug_unit
  character(len=20) :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=0
  name='psb_zsphalo'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Start'

  if (present(rowcnv)) then 
    rowcnv_ = rowcnv
  else
    rowcnv_ = .true.
  endif
  if (present(colcnv)) then 
    colcnv_ = colcnv
  else
    colcnv_ = .true.
  endif
  if (present(rowscale)) then 
    rowscale_ = rowscale
  else
    rowscale_ = .false.
  endif
  if (present(colscale)) then 
    colscale_ = colscale
  else
    colscale_ = .false.
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
  icomm = psb_cd_get_mpic(desc_a)

  Call psb_info(ictxt, me, np)

  Allocate(sdid(np,3),rvid(np,3),brvindx(np+1),&
       & rvsz(np),sdsz(np),bsdindx(np+1),stat=info)

  if (info /= 0) then
    info=4000
    call psb_errpush(info,name)
    goto 9999
  end if

  If (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Data selector',data_
  select case(data_) 
  case(psb_comm_halo_) 
    idxv => desc_a%halo_index

  case(psb_comm_ext_) 
    idxv => desc_a%ext_index

! !$  case(psb_comm_ovr_) 
! !$    idxv => desc_a%ovrlap_index
! Do not accept OVRLAP_INDEX any longer. 
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
  call psb_sp_all(blk,max(iszr,1),info)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Sizes:',size(blk%ia1),size(blk%ia2),&
       & ' Send:',sdsz(:),' Receive:',rvsz(:)
  if (info /= 0) then
    info=4010
    ch_err='psb_sp_reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  mat_recv = iszr
  iszs=sum(sdsz)
  call psb_ensure_size(max(iszs,1),iasnd,info)
  if (info == 0) call psb_ensure_size(max(iszs,1),jasnd,info)
  if (info == 0) call psb_ensure_size(max(iszs,1),valsnd,info)

  l1  = 0
  ipx = 1
  counter=1
  idx = 0

  tot_elem=0
  Do 
    proc=idxv(counter)
    if (proc == -1) exit 
    n_el_recv=idxv(counter+psb_n_elem_recv_)
    counter=counter+n_el_recv
    n_el_send=idxv(counter+psb_n_elem_send_)

    Do j=0,n_el_send-1
      idx = idxv(counter+psb_elem_send_+j)
      n_elem = psb_sp_get_nnz_row(idx,a)
      
      call psb_sp_getrow(idx,a,ngtz,iasnd,jasnd,valsnd,info,&
           &  append=.true.,nzin=tot_elem)
      if (info /= 0) then
        info=4010
        ch_err='psb_sp_getrow'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      tot_elem=tot_elem+n_elem
    Enddo

    ipx = ipx + 1 

    counter   = counter+n_el_send+3
  Enddo
  nz = tot_elem

  if (rowcnv_) call psb_loc_to_glob(iasnd(1:nz),desc_a,info,iact='I')
  if (colcnv_) call psb_loc_to_glob(jasnd(1:nz),desc_a,info,iact='I')
  if (info /= 0) then
    info=4010
    ch_err='psb_loc_to_glob'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  call mpi_alltoallv(valsnd,sdsz,bsdindx,mpi_double_complex,&
       & blk%aspk,rvsz,brvindx,mpi_double_complex,icomm,info)
  call mpi_alltoallv(iasnd,sdsz,bsdindx,mpi_integer,&
       & blk%ia1,rvsz,brvindx,mpi_integer,icomm,info)
  call mpi_alltoallv(jasnd,sdsz,bsdindx,mpi_integer,&
       & blk%ia2,rvsz,brvindx,mpi_integer,icomm,info)
  if (info /= 0) then
    info=4010
    ch_err='mpi_alltoallv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  !
  ! Convert into local numbering 
  !
  if (rowcnv_) call psb_glob_to_loc(blk%ia1(1:iszr),desc_a,info,iact='I')
  if (colcnv_) call psb_glob_to_loc(blk%ia2(1:iszr),desc_a,info,iact='I')

  if (info /= 0) then
    info=4010
    ch_err='psbglob_to_loc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  l1  = 0
  blk%m=0
  !
  irmin = huge(irmin)
  icmin = huge(icmin)
  irmax = 0
  icmax = 0
  Do i=1,iszr
    r=(blk%ia1(i))
    k=(blk%ia2(i))
    ! Just in case some of the conversions were out-of-range
    If ((r>0).and.(k>0)) Then
      l1=l1+1
      blk%aspk(l1) = blk%aspk(i)
      blk%ia1(l1) = r 
      blk%ia2(l1) = k
      irmin = min(irmin,r)
      irmax = max(irmax,r)
      icmin = min(icmin,k)
      icmax = max(icmax,k)
    End If
  Enddo
  if (rowscale_) then 
    blk%m = irmax-irmin+1
    blk%ia1(1:l1) = blk%ia1(1:l1) - irmin + 1
  else
    blk%m = irmax
  end if
  if (colscale_) then 
    blk%k = icmax-icmin+1
    blk%ia2(1:l1) = blk%ia2(1:l1) - icmin + 1
  else
    blk%k = icmax
  end if

  blk%fida            = 'COO'
  blk%infoa(psb_nnz_) = l1
  call psb_ensure_size(1,blk%pl,info)
  if (info == 0) call psb_ensure_size(1,blk%pr,info)
  if (info /= 0) then 
    info=4010
    ch_err='psb_ensure_size'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  blk%pl = 0
  blk%pr = 0

  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),&
       & ': End data exchange',counter,l1,blk%m

  ! Do we expect any duplicates to appear???? 
  call psb_spcnv(blk,info,afmt=outfmt_,dupl=psb_dupl_add_)
  if (info /= 0) then
    info=4010
    ch_err='psb_spcnv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  Deallocate(sdid,brvindx,rvid,bsdindx,rvsz,sdsz,&
       & iasnd,jasnd,valsnd,stat=info)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

End Subroutine psb_zsphalo
