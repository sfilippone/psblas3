!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006, 2010, 2015, 2017
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: psb_ssphalo.f90
!
! Subroutine: psb_ssphalo
!  This routine does the retrieval of remote matrix rows.                   
!  Note that retrieval is done through GTBLK, therefore it should work      
!  for any matrix format in A; as for the output, default is CSR. 
!    
! 
! Arguments: 
!    a        - type(psb_sspmat_type)   The local part of input matrix A
!    desc_a   - type(psb_desc_type).  The communication descriptor.
!    blck     - type(psb_sspmat_type)   The local part of output matrix BLCK
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
Subroutine psb_ssphalo(a,desc_a,blk,info,rowcnv,colcnv,&
     &  rowscale,colscale,outfmt,data)
  use psb_base_mod, psb_protect_name => psb_ssphalo

#ifdef MPI_MOD
  use mpi
#endif
  Implicit None
#ifdef MPI_H
  include 'mpif.h'
#endif

  Type(psb_sspmat_type),Intent(in)    :: a
  Type(psb_sspmat_type),Intent(inout) :: blk
  Type(psb_desc_type),Intent(in), target :: desc_a
  integer(psb_ipk_), intent(out)                :: info
  logical, optional, intent(in)       :: rowcnv,colcnv,rowscale,colscale
  character(len=5), optional          :: outfmt 
  integer(psb_ipk_), intent(in), optional       :: data
  !     ...local scalars....
  integer(psb_ipk_) :: np,me,counter,proc,i, &
       &     n_el_send,k,n_el_recv,ictxt, idx, r, tot_elem,&
       &     n_elem, j, ipx,mat_recv, iszs, iszr,idxs,idxr,nz,&
       &     irmin,icmin,irmax,icmax,data_,ngtz,totxch,nxs, nxr,&
       &     l1, err_act
  integer(psb_mpik_) :: icomm, minfo
  integer(psb_mpik_), allocatable  :: brvindx(:), &
       & rvsz(:), bsdindx(:),sdsz(:)
  integer(psb_ipk_), allocatable  :: iasnd(:), jasnd(:)
  real(psb_spk_), allocatable :: valsnd(:)
  type(psb_s_coo_sparse_mat), allocatable :: acoo
  integer(psb_ipk_), pointer  :: idxv(:)
  class(psb_i_base_vect_type), pointer :: pdxv
  integer(psb_ipk_), allocatable :: ipdxv(:)
  logical           :: rowcnv_,colcnv_,rowscale_,colscale_
  character(len=5)  :: outfmt_
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20) :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='psb_ssphalo'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()

  Call psb_info(ictxt, me, np)

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
    outfmt_ =  psb_toupper(outfmt)
  else
    outfmt_ = 'CSR'
  endif

  Allocate(brvindx(np+1),&
       & rvsz(np),sdsz(np),bsdindx(np+1), acoo,stat=info)

  if (info /= psb_success_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  If (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Data selector',data_
  select case(data_) 
  case(psb_comm_halo_,psb_comm_ext_ ) 
    ! Do not accept OVRLAP_INDEX any longer. 
  case default
    call psb_errpush(psb_err_from_subroutine_,name,a_err='wrong Data selector')
    goto 9999
  end select


  sdsz(:)=0
  rvsz(:)=0
  l1  = 0
  ipx = 1
  brvindx(ipx) = 0
  bsdindx(ipx) = 0
  counter=1
  idx = 0
  idxs = 0
  idxr = 0

  call acoo%allocate(izero,a%get_ncols(),info)


  call desc_a%get_list(data_,pdxv,totxch,nxr,nxs,info)
  ipdxv = pdxv%get_vect()
  ! For all rows in the halo descriptor, extract and send/receive.
  Do 
    proc=ipdxv(counter)
    if (proc == -1) exit
    n_el_recv = ipdxv(counter+psb_n_elem_recv_)
    counter   = counter+n_el_recv
    n_el_send = ipdxv(counter+psb_n_elem_send_)
    tot_elem = 0
    Do j=0,n_el_send-1
      idx = ipdxv(counter+psb_elem_send_+j)
      n_elem = a%get_nz_row(idx)
      tot_elem = tot_elem+n_elem      
    Enddo
    sdsz(proc+1) = tot_elem
    call acoo%set_nrows(acoo%get_nrows() + n_el_recv)
    counter   = counter+n_el_send+3
  Enddo

  call mpi_alltoall(sdsz,1,psb_mpi_def_integer,& 
       & rvsz,1,psb_mpi_def_integer,icomm,minfo)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='mpi_alltoall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  idxs = 0
  idxr = 0
  counter = 1
  Do 
    proc=ipdxv(counter)
    if (proc == -1) exit
    n_el_recv = ipdxv(counter+psb_n_elem_recv_)
    counter   = counter+n_el_recv
    n_el_send = ipdxv(counter+psb_n_elem_send_)

    bsdindx(proc+1) = idxs
    idxs = idxs + sdsz(proc+1)
    brvindx(proc+1) = idxr
    idxr = idxr + rvsz(proc+1)
    counter   = counter+n_el_send+3
  Enddo

  iszr=sum(rvsz)
  call acoo%reallocate(max(iszr,1))
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Sizes:',acoo%get_size(),&
       & ' Send:',sdsz(:),' Receive:',rvsz(:)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_sp_reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  mat_recv = iszr
  iszs=sum(sdsz)
  call psb_ensure_size(max(iszs,1),iasnd,info)
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),jasnd,info)
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),valsnd,info)

  l1  = 0
  ipx = 1
  counter=1
  idx = 0

  tot_elem=0
  Do 
    proc=ipdxv(counter)
    if (proc == -1) exit 
    n_el_recv=ipdxv(counter+psb_n_elem_recv_)
    counter=counter+n_el_recv
    n_el_send=ipdxv(counter+psb_n_elem_send_)

    Do j=0,n_el_send-1
      idx = ipdxv(counter+psb_elem_send_+j)
      n_elem = a%get_nz_row(idx)
      call a%csget(idx,idx,ngtz,iasnd,jasnd,valsnd,info,&
           &  append=.true.,nzin=tot_elem)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
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
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_loc_to_glob'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  call mpi_alltoallv(valsnd,sdsz,bsdindx,psb_mpi_r_spk_,&
       & acoo%val,rvsz,brvindx,psb_mpi_r_spk_,icomm,minfo)
  call mpi_alltoallv(iasnd,sdsz,bsdindx,psb_mpi_ipk_integer,&
       & acoo%ia,rvsz,brvindx,psb_mpi_ipk_integer,icomm,minfo)
  call mpi_alltoallv(jasnd,sdsz,bsdindx,psb_mpi_ipk_integer,&
       & acoo%ja,rvsz,brvindx,psb_mpi_ipk_integer,icomm,minfo)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='mpi_alltoallv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  !
  ! Convert into local numbering 
  !
  if (rowcnv_) call psb_glob_to_loc(acoo%ia(1:iszr),desc_a,info,iact='I')
  if (colcnv_) call psb_glob_to_loc(acoo%ja(1:iszr),desc_a,info,iact='I')

  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psbglob_to_loc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  l1  = 0
  call acoo%set_nrows(izero)
  !
  irmin = huge(irmin)
  icmin = huge(icmin)
  irmax = 0
  icmax = 0
  Do i=1,iszr
    r=(acoo%ia(i))
    k=(acoo%ja(i))
    ! Just in case some of the conversions were out-of-range
    If ((r>0).and.(k>0)) Then
      l1=l1+1
      acoo%val(l1) = acoo%val(i)
      acoo%ia(l1)  = r 
      acoo%ja(l1)  = k
      irmin = min(irmin,r)
      irmax = max(irmax,r)
      icmin = min(icmin,k)
      icmax = max(icmax,k)
    End If
  Enddo
  if (rowscale_) then 
    call acoo%set_nrows(max(irmax-irmin+1,0))
    acoo%ia(1:l1) = acoo%ia(1:l1) - irmin + 1
  else    
    call acoo%set_nrows(irmax)
  end if
  if (colscale_) then 
    call acoo%set_ncols(max(icmax-icmin+1,0))
    acoo%ja(1:l1) = acoo%ja(1:l1) - icmin + 1
  else
    call acoo%set_ncols(icmax)
  end if

  call acoo%set_nzeros(l1)
  call acoo%set_sorted(.false.)

  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),&
       & ': End data exchange',counter,l1

  call move_alloc(acoo,blk%a)

  ! Do we expect any duplicates to appear???? 
  call blk%cscnv(info,type=outfmt_,dupl=psb_dupl_add_)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_spcnv'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  Deallocate(brvindx,bsdindx,rvsz,sdsz,&
       & iasnd,jasnd,valsnd,stat=info)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

End Subroutine psb_ssphalo
