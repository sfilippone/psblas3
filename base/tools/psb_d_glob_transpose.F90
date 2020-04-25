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
! File: psb_d_glob_transpose.f90
!
! Subroutine: psb_d_glob_transpose
! Version:    real
!
!  This file provides multiple related versions of a parallel
!  global transpose 
!
!         B = A^T
!
!  There are some variations of this routine, that are accounted for
!  in the workhorse  psb_ld_coo_glob_transpose
!  1. The row and column spaces can share the same descriptor
!     This essentially means that the descriptor relates to a
!     matrix with a symmetric pattern. Examples are: symmetric
!     matrices, matrices with symmetric pattern, lower or
!     upper halves of such matrices
!  2. The row and column index spaces are different
!     In this case you need to have two descriptor on input,
!     plus if you want the output to be distributed according
!     to the row descriptor, you will still need a new descriptor
!     because even if the row distribution is the same, the pattern
!     will be different.
!
!  This is handled in the workhorse by having one mandatory and
!  two optional descriptors:
!  1. If only the mandatory descriptor is present, then it is assumed that it
!     is both row and column descriptor, and that it is sufficient.
!  2. If two descriptors are available, then use the second
!  3. If the third output descriptor is available, then rebuild it
!     after the data exchange.
!
!  The main transpose algorithm works like this:
!  1. Compute sizes: any entry A(I,J) with J in the halo will have
!     to be sent to the process owning J, so walk through the
!     matrix and compute all the send sizes, then do an alltoall to figure
!     the receive sizes;
!  2. Adjust send bufffers;
!  3. Perform a local transpose;
!  4. Split the matrix: all local entries stay, all halo entries go into
!     the send buffers, and are converted to global numbering;
!  5. Do the all-to-all (see below for a discussion of the alternative
!     communication strategies)
!  6. The receive is in the extra section of the ACOO buffer; convert
!     the row indices to local numbering, and discard extra ones (there will
!     be some)
!  7. If desc_rx was requested, make sure to insert the (new) column indices
!  8. Cleanup and sort the output matrix
!  9. Copy back into AIN or ATRANS if requested.
!
!  There are three possible exchange algorithms:
!  1. Use MPI_Alltoallv
!  2. Use psb_simple_a2av
!  3. Use psb_simple_triad_a2av
!  Default choice is 3. The MPI variant has proved to be inefficient;
!  that is because it is not persistent, therefore you pay the initialization price
!  every time, and it is not optimized for a sparse communication pattern,
!  most MPI implementations assume that all communications are non-empty.
!  The PSB_SIMPLE variants reuse the same communicator, and go for a simplistic
!  sequence of sends/receive that is quite efficient for a sparse communication
!  pattern. To be refined/reviewed in the future to compare with neighbour
!  persistent collectives. 
!
!
subroutine psb_ld_coo_glob_transpose(ain,desc_r,info,atrans,desc_c,desc_rx)
#ifdef MPI_MOD
  use mpi
#endif
  use psb_base_mod, psb_protect_name => psb_ld_coo_glob_transpose
  Implicit None
#ifdef MPI_H
  include 'mpif.h'
#endif
  type(psb_ld_coo_sparse_mat), intent(inout) :: ain
  type(psb_desc_type), intent(inout), target   :: desc_r
  type(psb_ld_coo_sparse_mat), intent(out), optional :: atrans
  type(psb_desc_type), intent(inout), target, optional :: desc_c
  type(psb_desc_type), intent(out), optional   :: desc_rx
  integer(psb_ipk_), intent(out)               :: info

  !     ...local scalars....
  integer(psb_ipk_) :: ictxt, np,me
  integer(psb_ipk_) :: counter,proc, err_act, j
  integer(psb_lpk_) :: i,  k, idx, r, ipx,mat_recv, iszs, iszr,idxs,idxr,nz,&
       &     l1, nsnds, nrcvs, nr,nc,nzl, hlstart, nzt, nzd
  integer(psb_mpk_) :: icomm, minfo
  integer(psb_mpk_), allocatable  :: brvindx(:), &
       & rvsz(:), bsdindx(:), sdsz(:), tsdx(:), trvx(:) 
  integer(psb_ipk_), allocatable  :: halo_owner(:)
  integer(psb_lpk_), allocatable  :: iasnd(:), jasnd(:)
  real(psb_dpk_), allocatable :: valsnd(:)
  type(psb_ld_coo_sparse_mat), allocatable :: acoo
  logical           :: rowcnv_,colcnv_,rowscale_,colscale_,outcol_glob_
  type(psb_desc_type), pointer :: p_desc_c
  character(len=5)  :: outfmt_
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20) :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='mld_glob_transpose'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_r%get_context()
  icomm = desc_r%get_mpic()

  Call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Start'

  if (present(desc_c)) then
    p_desc_c => desc_c
  else
    p_desc_c => desc_r
  end if

  Allocate(brvindx(np+1),&
       & rvsz(np),sdsz(np),bsdindx(np+1), acoo,stat=info)

  if (info /= psb_success_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if


  sdsz(:)=0
  rvsz(:)=0
  l1  = 0
  brvindx(1) = 0
  bsdindx(1) = 0
  counter=1
  idx  = 0
  idxs = 0
  idxr = 0

  if (present(atrans)) then 
    call ain%cp_to_coo(acoo,info)
  else
    call ain%mv_to_coo(acoo,info)
  end if


  !
  ! Compute number of entries in the
  ! halo part, sorted by destination process
  !
  nr = desc_r%get_local_rows()
  nc = p_desc_c%get_local_cols()
  nzl = acoo%get_nzeros()
  hlstart = p_desc_c%get_local_rows()  
  do k=1, nzl
    j = acoo%ja(k)
    if (j > hlstart) then
      call p_desc_c%indxmap%fnd_halo_owner(j,proc,info)
      sdsz(proc+1) = sdsz(proc+1) +1
    end if
  end do

  !
  ! Exchange  sizes, so as to know sends/receives.
  ! This is different from the halo exchange because the
  ! number of entries was not  precomputed in the descriptor,
  ! which was vector-oriented and not matrix-entry-oriented
  ! 
  call mpi_alltoall(sdsz,1,psb_mpi_mpk_,& 
       & rvsz,1,psb_mpi_mpk_,icomm,minfo)

  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='mpi_alltoall')
    goto 9999
  end if
  nsnds = count(sdsz /= 0)
  nrcvs = count(rvsz /= 0)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done initial alltoall',nsnds,nrcvs

  idxs = 0
  idxr = 0
  counter = 1
  Do proc = 0, np-1 
    bsdindx(proc+1) = idxs
    idxs            = idxs + sdsz(proc+1)
    brvindx(proc+1) = idxr
    idxr            = idxr + rvsz(proc+1)
  Enddo

  tsdx = bsdindx
  trvx = brvindx

  iszr = sum(rvsz)
  iszs = sum(sdsz)

  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Sizes:',&
       & ' Send:',sdsz(:),' Receive:',rvsz(:)

  call psb_ensure_size(max(iszs,1),iasnd,info)
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),jasnd,info)
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),valsnd,info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='ensure_size')
    goto 9999
  end if

  !
  ! Now, transpose the matrix, then split between itself
  ! and the send buffers
  !
  call acoo%transp()
  if (acoo%get_nzeros()/= nzl) then
    write(0,*) me,'Something strange upon transpose: ',nzl,acoo%get_nzeros()
  end if
  nzl = acoo%get_nzeros()
  hlstart = p_desc_c%get_local_rows()


  nzd = 0
  do k = 1, nzl
    j = acoo%ia(k)
    if (j<=hlstart) then
      nzd = nzd + 1
      acoo%ia(nzd)  = acoo%ia(k) 
      acoo%ja(nzd)  = acoo%ja(k) 
      acoo%val(nzd) = acoo%val(k)
    else
      call p_desc_c%indxmap%fnd_halo_owner(j,proc,info)
      tsdx(proc+1)         = tsdx(proc+1) +1
      iasnd(tsdx(proc+1))  = acoo%ia(k)
      jasnd(tsdx(proc+1))  = acoo%ja(k)
      valsnd(tsdx(proc+1)) = acoo%val(k)
    end if
  end do
  call acoo%set_nzeros(nzd)
  !
  !  Put halo entries in global numbering
  !
  call desc_r%indxmap%l2gip(jasnd(1:iszs),info)
  call p_desc_c%indxmap%l2gip(iasnd(1:iszs),info)
  ! And exchange data.
  ! Normally we'll use our SIMPLE A2AV and not MPI, because
  ! the communication pattern is sparse, so ours is more
  ! efficient. Using ACOO for the receive buffers.
  nzl = acoo%get_nzeros()
  call acoo%ensure_size(nzl+iszr)

  select case(psb_get_sp_a2av_alg())
  case(psb_sp_a2av_smpl_triad_)
    call psb_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
         & acoo%val(nzl+1:nzl+iszr),acoo%ia(nzl+1:nzl+iszr),&
         & acoo%ja(nzl+1:nzl+iszr),rvsz,brvindx,ictxt,info)
  case(psb_sp_a2av_smpl_v_)
    call psb_simple_a2av(valsnd,sdsz,bsdindx,&
         & acoo%val(nzl+1:nzl+iszr),rvsz,brvindx,ictxt,info)
    if (info == psb_success_) call psb_simple_a2av(iasnd,sdsz,bsdindx,&
         & acoo%ia(nzl+1:nzl+iszr),rvsz,brvindx,ictxt,info)
    if (info == psb_success_) call psb_simple_a2av(jasnd,sdsz,bsdindx,&
         & acoo%ja(nzl+1:nzl+iszr),rvsz,brvindx,ictxt,info)
  case(psb_sp_a2av_mpi_)
    
    call mpi_alltoallv(valsnd,sdsz,bsdindx,psb_mpi_r_dpk_,&
         & acoo%val(nzl+1:nzl+iszr),rvsz,brvindx,psb_mpi_r_dpk_,icomm,minfo)
    if (minfo == mpi_success) &
         & call mpi_alltoallv(iasnd,sdsz,bsdindx,psb_mpi_lpk_,&
         & acoo%ia(nzl+1:nzl+iszr),rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
    if (minfo == mpi_success) &
         & call mpi_alltoallv(jasnd,sdsz,bsdindx,psb_mpi_lpk_,&
         & acoo%ja(nzl+1:nzl+iszr),rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
    if (minfo /= mpi_success) info = minfo
  case default
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='wrong A2AV alg selector')
    goto 9999
  end select

  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='alltoallv')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done alltoallv'

  if (present(desc_rx)) then 
    !
    ! Extend the appropriate descriptor; started as R but on
    ! transpose it now describes C
    ! 
    call desc_r%clone(desc_rx,info)
    call psb_cd_reinit(desc_rx,info)
    !
    ! Take out non-local rows
    !      
    call psb_glob_to_loc(acoo%ia(nzl+1:nzl+iszr),p_desc_c,info,iact='I',owned=.true.)
    call psb_coo_clean_negidx_inner(iszr,acoo%ia(nzl+1:nzl+iszr),acoo%ja(nzl+1:nzl+iszr),&
         & acoo%val(nzl+1:nzl+iszr),nzt,info)
    call desc_rx%g2lip_ins(acoo%ja(nzl+1:nzl+nzt),info)
    call psb_cdasb(desc_rx,info)
    nzl = nzl + nzt
    call acoo%set_nzeros(nzl)
    nzl = acoo%get_nzeros()
    call acoo%set_sorted(.false.)
    !
    ! Insert to extend descriptor
    !
    call acoo%set_nrows(p_desc_c%get_local_rows())
    call acoo%set_ncols(desc_rx%get_local_cols())
    !write(0,*) me,' Trans RX ',acoo%get_nrows(),acoo%get_ncols(),acoo%get_nzeros()
  else
    !
    !
    ! Take out non-local rows
    !
    call psb_glob_to_loc(acoo%ia(nzl+1:nzl+iszr),p_desc_c,info,iact='I',owned=.true.)
    call psb_glob_to_loc(acoo%ja(nzl+1:nzl+iszr),desc_r,info,iact='I')
    call psb_coo_clean_negidx_inner(iszr,acoo%ia(nzl+1:nzl+iszr),acoo%ja(nzl+1:nzl+iszr),&
         & acoo%val(nzl+1:nzl+iszr),nzt,info)

    nzl = nzl + nzt
    call acoo%set_nzeros(nzl)
    nzl = acoo%get_nzeros()
    call acoo%set_sorted(.false.)

    call acoo%set_nrows(p_desc_c%get_local_rows())
    call acoo%set_ncols(desc_r%get_local_cols())
    !write(0,*) me,' Trans R- ',acoo%get_nrows(),acoo%get_ncols(),acoo%get_nzeros()
  end if

!!$  write(0,*) me,' Sanity check after rx%g2l :',count(acoo%ja(1:nzl)<0)


  call acoo%fix(info)
  nzl = acoo%get_nzeros()

  if (present(atrans)) then
    call atrans%mv_from_coo(acoo,info)
  else
    call ain%mv_from_coo(acoo,info)
  end if
  Deallocate(brvindx,bsdindx,rvsz,sdsz,&
       & iasnd,jasnd,valsnd,stat=info)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_ld_coo_glob_transpose

subroutine psb_d_coo_glob_transpose(ain,desc_r,info,atrans,desc_c,desc_rx)
#ifdef MPI_MOD
  use mpi
#endif
  use psb_base_mod, psb_protect_name => psb_d_coo_glob_transpose
  Implicit None
#ifdef MPI_H
  include 'mpif.h'
#endif
  type(psb_d_coo_sparse_mat), intent(inout) :: ain
  type(psb_desc_type), intent(inout), target   :: desc_r
  type(psb_d_coo_sparse_mat), intent(out), optional :: atrans
  type(psb_desc_type), intent(inout), target, optional :: desc_c
  type(psb_desc_type), intent(out), optional   :: desc_rx
  integer(psb_ipk_), intent(out)               :: info
  
  integer(psb_ipk_) :: ictxt, np,me
  integer(psb_ipk_) :: counter,proc, err_act, j
  integer(psb_ipk_) :: i,  k, idx, r, ipx,mat_recv, iszs, iszr,idxs,idxr,nz,&
       &     l1, nsnds, nrcvs, nr,nc,nzl, hlstart, nzd
  integer(psb_lpk_) :: nzt, lszr
  integer(psb_mpk_) :: icomm, minfo
  integer(psb_mpk_), allocatable  :: brvindx(:), &
       & rvsz(:), bsdindx(:), sdsz(:), tsdx(:), trvx(:) 
  integer(psb_ipk_), allocatable  :: halo_owner(:)
  integer(psb_lpk_), allocatable  :: iasnd(:), jasnd(:),iarcv(:),jarcv(:)
  real(psb_dpk_), allocatable :: valsnd(:)
  type(psb_d_coo_sparse_mat), allocatable :: acoo
  logical           :: rowcnv_,colcnv_,rowscale_,colscale_,outcol_glob_
  type(psb_desc_type), pointer :: p_desc_c
  character(len=5)  :: outfmt_
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20) :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='mld_glob_transpose'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_r%get_context()
  icomm = desc_r%get_mpic()

  Call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Start'

  if (present(desc_c)) then
    p_desc_c => desc_c
  else
    p_desc_c => desc_r
  end if

  Allocate(brvindx(np+1),&
       & rvsz(np),sdsz(np),bsdindx(np+1), acoo,stat=info)

  if (info /= psb_success_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if


  sdsz(:)=0
  rvsz(:)=0
  l1  = 0
  brvindx(1) = 0
  bsdindx(1) = 0
  counter=1
  idx  = 0
  idxs = 0
  idxr = 0

  if (present(atrans)) then 
    call ain%cp_to_coo(acoo,info)
  else
    call ain%mv_to_coo(acoo,info)
  end if


  !
  ! Compute number of entries in the
  ! halo part, sorted by destination process
  !
  nr  = desc_r%get_local_rows()
  nc  = p_desc_c%get_local_cols()
  nzl = acoo%get_nzeros()
  hlstart = p_desc_c%get_local_rows()  
  do k=1, nzl
    j = acoo%ja(k)
    if (j > hlstart) then
      call p_desc_c%indxmap%fnd_halo_owner(j,proc,info)
      sdsz(proc+1) = sdsz(proc+1) +1
    end if
  end do

  !
  ! Exchange  sizes, so as to know sends/receives.
  ! This is different from the halo exchange because the
  ! number of entries was not  precomputed in the descriptor,
  ! which was vector-oriented and not matrix-entry-oriented
  ! 
  call mpi_alltoall(sdsz,1,psb_mpi_mpk_,& 
       & rvsz,1,psb_mpi_mpk_,icomm,minfo)

  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='mpi_alltoall')
    goto 9999
  end if
  nsnds = count(sdsz /= 0)
  nrcvs = count(rvsz /= 0)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done initial alltoall',nsnds,nrcvs

  idxs = 0
  idxr = 0
  counter = 1
  Do proc = 0, np-1 
    bsdindx(proc+1) = idxs
    idxs            = idxs + sdsz(proc+1)
    brvindx(proc+1) = idxr
    idxr            = idxr + rvsz(proc+1)
  Enddo

  tsdx = bsdindx
  trvx = brvindx

  iszr = sum(rvsz)
  iszs = sum(sdsz)

  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Sizes:',&
       & ' Send:',sdsz(:),' Receive:',rvsz(:)

  call psb_ensure_size(max(iszs,1),iasnd,info)
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),jasnd,info)
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),valsnd,info)
  if (info == psb_success_) call psb_ensure_size(max(iszr,1),iarcv,info)
  if (info == psb_success_) call psb_ensure_size(max(iszr,1),jarcv,info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='ensure_size')
    goto 9999
  end if

  !
  ! Now, transpose the matrix, then split between itself
  ! and the send buffers
  !
  call acoo%transp()
  if (acoo%get_nzeros()/= nzl) then
    write(0,*) me,'Something strange upon transpose: ',nzl,acoo%get_nzeros()
  end if
  nzl     = acoo%get_nzeros()
  hlstart = p_desc_c%get_local_rows()


  nzd = 0
  do k = 1, nzl
    j = acoo%ia(k)
    if (j<=hlstart) then
      nzd = nzd + 1
      acoo%ia(nzd)  = acoo%ia(k) 
      acoo%ja(nzd)  = acoo%ja(k) 
      acoo%val(nzd) = acoo%val(k)
    else
      call p_desc_c%indxmap%fnd_halo_owner(j,proc,info)
      tsdx(proc+1)         = tsdx(proc+1) +1
      iasnd(tsdx(proc+1))  = acoo%ia(k)
      jasnd(tsdx(proc+1))  = acoo%ja(k)
      valsnd(tsdx(proc+1)) = acoo%val(k)
    end if
  end do
  call acoo%set_nzeros(nzd)
  !
  !  Put halo entries in global numbering
  !
  call desc_r%indxmap%l2gip(jasnd(1:iszs),info)
  call p_desc_c%indxmap%l2gip(iasnd(1:iszs),info)
  ! And exchange data.
  ! Normally we'll use our SIMPLE A2AV and not MPI, because
  ! the communication pattern is sparse, so ours is more
  ! efficient. Using ACOO for the receive buffers.
  nzl = acoo%get_nzeros()
  call acoo%ensure_size(nzl+iszr)

  select case(psb_get_sp_a2av_alg())
  case(psb_sp_a2av_smpl_triad_)
    call psb_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
         & acoo%val(nzl+1:nzl+iszr),iarcv(1:iszr),&
         & jarcv(1:iszr),rvsz,brvindx,ictxt,info)
  case(psb_sp_a2av_smpl_v_)
    call psb_simple_a2av(valsnd,sdsz,bsdindx,&
         & acoo%val(nzl+1:nzl+iszr),rvsz,brvindx,ictxt,info)
    if (info == psb_success_) call psb_simple_a2av(iasnd,sdsz,bsdindx,&
         & iarcv(1:iszr),rvsz,brvindx,ictxt,info)
    if (info == psb_success_) call psb_simple_a2av(jasnd,sdsz,bsdindx,&
         & jarcv(1:iszr),rvsz,brvindx,ictxt,info)
  case(psb_sp_a2av_mpi_)
    
    call mpi_alltoallv(valsnd,sdsz,bsdindx,psb_mpi_r_dpk_,&
         & acoo%val(nzl+1:nzl+iszr),rvsz,brvindx,psb_mpi_r_dpk_,icomm,minfo)
    if (minfo == mpi_success) &
         & call mpi_alltoallv(iasnd,sdsz,bsdindx,psb_mpi_lpk_,&
         & iarcv(1:iszr),rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
    if (minfo == mpi_success) &
         & call mpi_alltoallv(jasnd,sdsz,bsdindx,psb_mpi_lpk_,&
         & jarcv(1:iszr),rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
    if (minfo /= mpi_success) info = minfo
  case default
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='wrong A2AV alg selector')
    goto 9999
  end select

  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='alltoallv')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done alltoallv'

  if (present(desc_rx)) then 
    !
    ! Extend the appropriate descriptor; started as R but on
    ! transpose it now describes C
    ! 
    call desc_r%clone(desc_rx,info)
    call psb_cd_reinit(desc_rx,info)
    !
    ! Take out non-local rows
    !      
    call psb_glob_to_loc(iarcv(1:iszr),p_desc_c,info,iact='I',owned=.true.)
    lszr = iszr 
    call psb_coo_clean_negidx_inner(lszr,iarcv(1:iszr),jarcv(1:iszr),&
         & acoo%val(nzl+1:nzl+iszr),nzt,info)
    call desc_rx%g2lip_ins(jarcv(1:nzt),info)
    call psb_cdasb(desc_rx,info)
    acoo%ia(nzl+1:nzl+nzt) = iarcv(1:nzt)
    acoo%ja(nzl+1:nzl+nzt) = jarcv(1:nzt)
    nzl = nzl + nzt
    call acoo%set_nzeros(nzl)
    nzl = acoo%get_nzeros()
    call acoo%set_sorted(.false.)
    !
    ! Insert to extend descriptor
    !
    call acoo%set_nrows(p_desc_c%get_local_rows())
    call acoo%set_ncols(desc_rx%get_local_cols())
    !write(0,*) me,' Trans RX ',acoo%get_nrows(),acoo%get_ncols(),acoo%get_nzeros()
  else
    !
    !
    ! Take out non-local rows
    !
    call psb_glob_to_loc(iarcv(1:iszr),p_desc_c,info,iact='I',owned=.true.)
    call psb_glob_to_loc(jarcv(1:iszr),desc_r,info,iact='I')
    lszr = iszr 
    call psb_coo_clean_negidx_inner(lszr,iarcv(1:iszr),jarcv(1:iszr),&
         & acoo%val(nzl+1:nzl+iszr),nzt,info)
    acoo%ia(nzl+1:nzl+nzt) = iarcv(1:nzt)
    acoo%ja(nzl+1:nzl+nzt) = jarcv(1:nzt)
    nzl = nzl + nzt
    call acoo%set_nzeros(nzl)
    nzl = acoo%get_nzeros()
    call acoo%set_sorted(.false.)

    call acoo%set_nrows(p_desc_c%get_local_rows())
    call acoo%set_ncols(desc_r%get_local_cols())
    !write(0,*) me,' Trans R- ',acoo%get_nrows(),acoo%get_ncols(),acoo%get_nzeros()
  end if

!!$  write(0,*) me,' Sanity check after rx%g2l :',count(acoo%ja(1:nzl)<0)


  call acoo%fix(info)
  nzl = acoo%get_nzeros()

  if (present(atrans)) then
    call atrans%mv_from_coo(acoo,info)
  else
    call ain%mv_from_coo(acoo,info)
  end if
  Deallocate(brvindx,bsdindx,rvsz,sdsz,&
       & iasnd,jasnd,valsnd,&
       & iarcv,jarcv,stat=info)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_d_coo_glob_transpose

subroutine psb_d_simple_glob_transpose_ip(ain,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_d_simple_glob_transpose_ip
  implicit none
  type(psb_dspmat_type), intent(inout)  :: ain
  type(psb_desc_type)           :: desc_a
  integer(psb_ipk_), intent(out) :: info

  !
  ! BEWARE: This routine works under the assumption
  ! that the same DESC_A works for both A and A^T, which
  ! essentially means that A has a symmetric pattern.
  !
  type(psb_d_coo_sparse_mat) :: tmpc1, tmpc2
  integer(psb_ipk_) :: nz1, nz2, nzh, nz
  integer(psb_ipk_) :: ictxt, me, np
  integer(psb_lpk_) :: i, j, k, nrow, ncol, nlz
  integer(psb_lpk_), allocatable :: ilv(:) 
  character(len=80) :: aname
  logical, parameter :: debug=.false., dump=.false., debug_sync=.false.

  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  if (debug_sync) then
    call psb_barrier(ictxt)
    if (me == 0) write(0,*) 'Start htranspose '
  end if


  call ain%mv_to(tmpc1)
  call psb_glob_transpose(tmpc1, desc_a,info,atrans=tmpc2)
  call ain%mv_from(tmpc2)

  if (dump) then
    block
      type(psb_ldspmat_type) :: aglb
      write(aname,'(a,i3.3,a)') 'atran-',me,'.mtx'    
      call ain%print(fname=aname,head='atrans ')
      call psb_gather(aglb,ain,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'atran.mtx'    
        call aglb%print(fname=aname,head='Test ')
      end if
    end block
  end if

end subroutine psb_d_simple_glob_transpose_ip

subroutine psb_d_simple_glob_transpose(ain,aout,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_d_simple_glob_transpose
  implicit none
  type(psb_dspmat_type), intent(in)  :: ain
  type(psb_dspmat_type), intent(out) :: aout
  type(psb_desc_type)           :: desc_a
  integer(psb_ipk_), intent(out) :: info

  !
  ! BEWARE: This routine works under the assumption
  ! that the same DESC_A works for both A and A^T, which
  ! essentially means that A has a symmetric pattern.
  !
  type(psb_d_coo_sparse_mat) :: tmpc1, tmpc2
  integer(psb_ipk_) :: nz1, nz2, nzh, nz
  integer(psb_ipk_) :: ictxt, me, np
  integer(psb_lpk_) :: i, j, k, nrow, ncol, nlz
  integer(psb_lpk_), allocatable :: ilv(:) 
  character(len=80) :: aname
  logical, parameter :: debug=.false., dump=.false., debug_sync=.false.

  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  if (debug_sync) then
    call psb_barrier(ictxt)
    if (me == 0) write(0,*) 'Start htranspose '
  end if


  call ain%cp_to(tmpc1)
  call psb_glob_transpose(tmpc1, desc_a,info,atrans=tmpc2)
  call aout%mv_from(tmpc2)

  if (dump) then
    block
      type(psb_ldspmat_type) :: aglb
      
      write(aname,'(a,i3.3,a)') 'atran-',me,'.mtx'    
      call aout%print(fname=aname,head='atrans ')
      call psb_gather(aglb,aout,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'atran.mtx'    
        call aglb%print(fname=aname,head='Test ')
      end if
    end block
  end if

end subroutine psb_d_simple_glob_transpose

subroutine psb_ld_simple_glob_transpose_ip(ain,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_ld_simple_glob_transpose_ip
  implicit none
  type(psb_ldspmat_type), intent(inout)  :: ain
  type(psb_desc_type)           :: desc_a
  integer(psb_ipk_), intent(out) :: info

  !
  ! BEWARE: This routine works under the assumption
  ! that the same DESC_A works for both A and A^T, which
  ! essentially means that A has a symmetric pattern.
  !
  type(psb_ld_coo_sparse_mat) :: tmpc1, tmpc2
  integer(psb_ipk_) :: nz1, nz2, nzh, nz
  integer(psb_ipk_) :: ictxt, me, np
  integer(psb_lpk_) :: i, j, k, nrow, ncol, nlz
  integer(psb_lpk_), allocatable :: ilv(:) 
  character(len=80) :: aname
  logical, parameter :: debug=.false., dump=.false., debug_sync=.false.

  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  if (debug_sync) then
    call psb_barrier(ictxt)
    if (me == 0) write(0,*) 'Start htranspose '
  end if


  call ain%mv_to(tmpc1)
  call psb_glob_transpose(tmpc1, desc_a,info,atrans=tmpc2)
  call ain%mv_from(tmpc2)

  if (dump) then
    block
      type(psb_ldspmat_type) :: aglb
      write(aname,'(a,i3.3,a)') 'atran-',me,'.mtx'    
      call ain%print(fname=aname,head='atrans ',iv=ilv)
      call psb_gather(aglb,ain,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'atran.mtx'    
        call aglb%print(fname=aname,head='Test ')
      end if
    end block
  end if

end subroutine psb_ld_simple_glob_transpose_ip

subroutine psb_ld_simple_glob_transpose(ain,aout,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_ld_simple_glob_transpose
  implicit none
  type(psb_ldspmat_type), intent(in)  :: ain
  type(psb_ldspmat_type), intent(out) :: aout
  type(psb_desc_type)           :: desc_a
  integer(psb_ipk_), intent(out) :: info

  !
  ! BEWARE: This routine works under the assumption
  ! that the same DESC_A works for both A and A^T, which
  ! essentially means that A has a symmetric pattern.
  !
  type(psb_ld_coo_sparse_mat) :: tmpc1, tmpc2
  integer(psb_ipk_) :: nz1, nz2, nzh, nz
  integer(psb_ipk_) :: ictxt, me, np
  integer(psb_lpk_) :: i, j, k, nrow, ncol, nlz
  integer(psb_lpk_), allocatable :: ilv(:) 
  character(len=80) :: aname
  logical, parameter :: debug=.false., dump=.false., debug_sync=.false.

  ictxt = desc_a%get_context()
  call psb_info(ictxt,me,np)

  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  if (debug_sync) then
    call psb_barrier(ictxt)
    if (me == 0) write(0,*) 'Start htranspose '
  end if


  call ain%cp_to(tmpc1)
  call psb_glob_transpose(tmpc1, desc_a,info,atrans=tmpc2)
  call aout%mv_from(tmpc2)

  if (dump) then
    block
      type(psb_ldspmat_type) :: aglb
      
      write(aname,'(a,i3.3,a)') 'atran-',me,'.mtx'    
      call aout%print(fname=aname,head='atrans ',iv=ilv)
      call psb_gather(aglb,aout,desc_a,info)
      if (me==psb_root_) then
        write(aname,'(a,i3.3,a)') 'atran.mtx'    
        call aglb%print(fname=aname,head='Test ')
      end if
    end block
  end if

end subroutine psb_ld_simple_glob_transpose


