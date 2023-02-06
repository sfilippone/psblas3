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
! File: psb_e_remote_vect.f90
!
! Subroutine: 
!  This routine does the retrieval of remote vector entries.                   
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
! Arguments: 
!    desc_a   - type(psb_desc_type).  The communication descriptor.
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
subroutine psb_e_remote_vect(n,v,iv,desc_a,x,ix, info)
  use psb_base_mod, psb_protect_name => psb_e_remote_vect

#ifdef MPI_MOD
  use mpi
#endif
  Implicit None
#ifdef MPI_H
  include 'mpif.h'
#endif
  integer(psb_ipk_), intent(in)  :: n
  integer(psb_epk_),   Intent(in)  :: v(:)
  integer(psb_lpk_), Intent(in)  :: iv(:)
  type(psb_desc_type),intent(in) :: desc_a
  integer(psb_epk_),   allocatable, intent(out) :: x(:)
  integer(psb_lpk_), allocatable, intent(out) :: ix(:)
  integer(psb_ipk_), intent(out) :: info
  !     ...local scalars....
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: counter, proc, i,  &
       &     j,  idxs,idxr, k, iszs, iszr
  integer(psb_ipk_) :: nrcvs, nsnds
  integer(psb_mpk_) :: icomm, minfo
  integer(psb_mpk_), allocatable  :: brvindx(:), &
       & rvsz(:), bsdindx(:), sdsz(:), sdsi(:), rvsi(:) 
  integer(psb_lpk_), allocatable  :: lsnd(:)
  integer(psb_epk_), allocatable :: valsnd(:)
  integer(psb_ipk_), allocatable :: iprc(:)
  integer(psb_ipk_) :: debug_level, debug_unit, err_act
  character(len=20) :: name, ch_err

  info=psb_success_
  name='psb_e_remote_vect'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt  = desc_a%get_context()
  icomm = desc_a%get_mpic()

  Call psb_info(ctxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Start'

  Allocate(rvsz(np),sdsz(np),sdsi(np),rvsi(np),brvindx(np+1),&
       & bsdindx(np+1), stat=info)

  if (info /= psb_success_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  call desc_a%indxmap%fnd_owner(iv(1:n),iprc,info)

  icomm   = desc_a%get_mpic()
  sdsz(:) = 0
  rvsz(:) = 0
  sdsi(:) = 0
  rvsi(:) = 0
  brvindx(:) = 0
  bsdindx(:) = 0
  counter = 1
  idxs    = 0
  idxr    = 0
  do  i=1,n
    if (iprc(i) >=0) then
      sdsz(iprc(i)+1) = sdsz(iprc(i)+1) +1
    else
      write(0,*)me,name,' Error from fnd_owner: ',iprc(i)
    end if
  end do
  call mpi_alltoall(sdsz,1,psb_mpi_mpk_,& 
       & rvsz,1,psb_mpi_mpk_,icomm,minfo)
  if (minfo /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='mpi_alltoall')
    goto 9999
  end if
  !write(0,*)me,name,' sdsz ',sdsz(:),' rvsz:',rvsz(:)        
  nsnds = count(sdsz /= 0)
  nrcvs = count(rvsz /= 0)
  idxs  = 0
  idxr  = 0
  counter = 1
  Do proc=0,np-1
    bsdindx(proc+1) = idxs
    idxs            = idxs + sdsz(proc+1)
    brvindx(proc+1) = idxr
    idxr     = idxr + rvsz(proc+1)
  Enddo

  iszs = sum(sdsz)
  iszr = sum(rvsz)
  call psb_realloc(iszs,lsnd,info)
  if (info == 0) call psb_realloc(iszs,valsnd,info)
  if (info == 0) call psb_realloc(iszr,x,info)
  if (info == 0) call psb_realloc(iszr,ix,info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='realloc')
    goto 9999
  end if
  do k=1, n
    proc           = iprc(k)
    sdsi(proc+1)   = sdsi(proc+1) + 1
    lsnd(bsdindx(proc+1)+sdsi(proc+1))   = iv(k)
    valsnd(bsdindx(proc+1)+sdsi(proc+1)) = v(k)
  end do
  do proc=0,np-1
    if (sdsi(proc+1) /= sdsz(proc+1)) &
         & write(0,*) me,name,'Send mismacth ',sdsi(proc+1),sdsz(proc+1)
  end do

  select case(psb_get_sp_a2av_alg())
  case(psb_sp_a2av_smpl_triad_,psb_sp_a2av_smpl_v_)
    call psb_simple_a2av(valsnd,sdsz,bsdindx,&
         & x,rvsz,brvindx,ctxt,info)
    if (info == psb_success_) call psb_simple_a2av(lsnd,sdsz,bsdindx,&
         & ix,rvsz,brvindx,ctxt,info)
  case(psb_sp_a2av_mpi_)

    call mpi_alltoallv(valsnd,sdsz,bsdindx,psb_mpi_epk_,&
         & x,rvsz,brvindx,psb_mpi_epk_,icomm,minfo)
    if (minfo == mpi_success) &
         & call mpi_alltoallv(lsnd,sdsz,bsdindx,psb_mpi_lpk_,&
         & ix,rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
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
  
  Deallocate(brvindx,bsdindx,rvsz,sdsz,&
       & lsnd,valsnd,stat=info)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

End Subroutine psb_e_remote_vect
