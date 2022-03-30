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
! File: psb_z_remote_mat.f90
!
! Subroutine: 
!  This routine does the retrieval of remote matrix rows.                   
!  Retrieval is done through GETROW, therefore it should work      
!  for any matrix format in A; as for the output, default is CSR.
!  
!  There is also a specialized version lz_CSR whose interface
!  is adapted for the needs of z_par_csr_spspmm. 
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
Subroutine psb_lz_remote_mat(a,desc_a,b,info)
  use psb_base_mod, psb_protect_name => psb_lz_remote_mat

#ifdef MPI_MOD
  use mpi
#endif
  Implicit None
#ifdef MPI_H
  include 'mpif.h'
#endif

  Type(psb_lz_coo_sparse_mat),Intent(inout)  :: a
  Type(psb_lz_coo_sparse_mat),Intent(inout)  :: b
  Type(psb_desc_type), Intent(inout)         :: desc_a
  integer(psb_ipk_), intent(out)                :: info

  !     ...local scalars....
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: counter, proc, i, n_el_send,n_el_recv, &
       &     n_elem, j, ipx,mat_recv, idxs,idxr,&
       &     data_,totxch,nxs, nxr, ncg
  integer(psb_lpk_) :: r, k, irmin, irmax, icmin, icmax, iszs, iszr, &
       & lidx, l1, lnr, lnc, lnnz, idx, ngtz, tot_elem
  integer(psb_lpk_) :: nz,nouth
  integer(psb_ipk_) :: nnp, nrcvs, nsnds
  integer(psb_mpk_) :: icomm, minfo
  integer(psb_mpk_), allocatable  :: brvindx(:), &
       & rvsz(:), bsdindx(:),sdsz(:), sdsi(:), rvsi(:) 
  integer(psb_lpk_), allocatable  :: iasnd(:), jasnd(:)
  complex(psb_dpk_), allocatable :: valsnd(:)
  type(psb_lz_coo_sparse_mat), allocatable :: acoo
  class(psb_i_base_vect_type), pointer :: pdxv
  integer(psb_ipk_), allocatable :: ladj(:), ila(:), iprc(:)
  logical           :: rowcnv_,colcnv_,rowscale_,colscale_
  character(len=5)  :: outfmt_
  integer(psb_ipk_) :: debug_level, debug_unit, err_act
  character(len=20) :: name, ch_err

  info=psb_success_
  name='psb_z_remote_mat'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc_a%get_context()
  icomm = desc_a%get_mpic()

  Call psb_info(ctxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Start'


  call b%free() 

  Allocate(rvsz(np),sdsz(np),sdsi(np),rvsi(np),brvindx(np+1),&
       & bsdindx(np+1), acoo,stat=info)

  if (info /= psb_success_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if


  nz = a%get_nzeros()
  allocate(ila(nz))
  !write(0,*) me,name,' size :',nz,size(ila)        
  call desc_a%g2l(a%ia(1:nz),ila(1:nz),info,owned=.false.)
  nouth = count(ila(1:nz)<0)
  !write(0,*) me,name,' Count out of halo :',nouth
  call psb_max(ctxt,nouth)
  if ((nouth/=0).and.(me==0)) &
       & write(0,*) 'Warning: would require reinit of DESC_A'

  call psi_graph_fnd_owner(a%ia(1:nz),iprc,ladj,desc_a%indxmap,info)
  call psb_msort_unique(ladj,nnp)
  !write(0,*) me,name,' Processes:',ladj(1:nnp)

  icomm = desc_a%get_mpic()
  sdsz(:)=0
  rvsz(:)=0
  sdsi(:)=0
  rvsi(:)=0
  ipx = 1
  brvindx(:) = 0
  bsdindx(:) = 0
  counter=1
  idx = 0
  idxs = 0
  idxr = 0
  do  i=1,nz
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
  idxs = 0
  idxr = 0
  counter = 1
  Do proc=0,np-1
    bsdindx(proc+1) = idxs
    idxs            = idxs + sdsz(proc+1)
    brvindx(proc+1) = idxr
    idxr     = idxr + rvsz(proc+1)
  Enddo

  iszs = sum(sdsz)
  iszr = sum(rvsz)
  call acoo%allocate(desc_a%get_global_rows(),desc_a%get_global_cols(),iszr)
  if (psb_errstatus_fatal()) then
    write(0,*) 'Error from acoo%allocate '
    info = 4010
    goto 9999
  end if
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Sizes:',acoo%get_size(),&
       & ' Send:',sdsz(:),' Receive:',rvsz(:)
  !write(debug_unit,*) me,' ',trim(name),': ',info
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),iasnd,info)
  !write(debug_unit,*) me,' ',trim(name),' iasnd: ',info
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),jasnd,info)
  !write(debug_unit,*) me,' ',trim(name),' jasnd: ',info
  if (info == psb_success_) call psb_ensure_size(max(iszs,1),valsnd,info)
  !write(debug_unit,*) me,' ',trim(name),' valsnd: ',info
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='ensure_size')
    goto 9999
  end if
  do k=1, nz
    proc = iprc(k)
    sdsi(proc+1) = sdsi(proc+1) + 1
    !rvsi(proc) = rvsi(proc) + 1 
    iasnd(bsdindx(proc+1)+sdsi(proc+1))  = a%ia(k)
    jasnd(bsdindx(proc+1)+sdsi(proc+1))  = a%ja(k)
    valsnd(bsdindx(proc+1)+sdsi(proc+1)) = a%val(k)
  end do
  do proc=0,np-1
    if (sdsi(proc+1) /= sdsz(proc+1)) &
         & write(0,*) me,name,'Send mismacth ',sdsi(proc+1),sdsz(proc+1)
  end do

  select case(psb_get_sp_a2av_alg())
  case(psb_sp_a2av_smpl_triad_)
    call psb_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
         & acoo%val,acoo%ia,acoo%ja,rvsz,brvindx,ctxt,info)
  case(psb_sp_a2av_smpl_v_)
    call psb_simple_a2av(valsnd,sdsz,bsdindx,&
         & acoo%val,rvsz,brvindx,ctxt,info)
    if (info == psb_success_) call psb_simple_a2av(iasnd,sdsz,bsdindx,&
         & acoo%ia,rvsz,brvindx,ctxt,info)
    if (info == psb_success_) call psb_simple_a2av(jasnd,sdsz,bsdindx,&
         & acoo%ja,rvsz,brvindx,ctxt,info)
  case(psb_sp_a2av_mpi_)

    call mpi_alltoallv(valsnd,sdsz,bsdindx,psb_mpi_c_dpk_,&
         & acoo%val,rvsz,brvindx,psb_mpi_c_dpk_,icomm,minfo)
    if (minfo == mpi_success) &
         & call mpi_alltoallv(iasnd,sdsz,bsdindx,psb_mpi_lpk_,&
         & acoo%ia,rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
    if (minfo == mpi_success) &
         & call mpi_alltoallv(jasnd,sdsz,bsdindx,psb_mpi_lpk_,&
         & acoo%ja,rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
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
  call acoo%set_nzeros(iszr)
  call acoo%mv_to_coo(b,info)
  
  Deallocate(brvindx,bsdindx,rvsz,sdsz,&
       & iasnd,jasnd,valsnd,stat=info)
  if (debug_level >= psb_debug_outer_)&
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

End Subroutine psb_lz_remote_mat


subroutine psb_z_remote_vect(v,desc_a, info, dupl)
  use psb_base_mod, psb_protect_name => psb_z_remote_vect

#ifdef MPI_MOD
  use mpi
#endif
  Implicit None
#ifdef MPI_H
  include 'mpif.h'
#endif

  type(psb_z_vect_type),Intent(inout)    :: v
  type(psb_desc_type),intent(in)         :: desc_a
  integer(psb_ipk_), intent(out)         :: info
  integer(psb_ipk_),optional, intent(in) :: dupl

  !     ...local scalars....
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: counter, proc, i, n_el_send,n_el_recv, &
       &     n_elem, j, ipx,mat_recv, idxs,idxr,&
       &     data_,totxch,nxs, nxr, ncg, dupl_
  integer(psb_lpk_) :: r, k, irmin, irmax, icmin, icmax, iszs, iszr, &
       & lidx, l1, lnr, lnc, lnnz, idx, ngtz, tot_elem
  integer(psb_lpk_) :: nz,nouth
  integer(psb_ipk_) :: nnp, nrcvs, nsnds
  integer(psb_mpk_) :: icomm, minfo
  integer(psb_mpk_), allocatable  :: brvindx(:), &
       & rvsz(:), bsdindx(:),sdsz(:), sdsi(:), rvsi(:) 
  integer(psb_lpk_), allocatable  :: iasnd(:), jasnd(:)
  complex(psb_dpk_), allocatable :: valsnd(:)
  integer(psb_ipk_), allocatable :: ladj(:), ila(:), iprc(:)
  logical           :: rowcnv_,colcnv_,rowscale_,colscale_
  character(len=5)  :: outfmt_
  integer(psb_ipk_) :: debug_level, debug_unit, err_act
  character(len=20) :: name, ch_err

  info=psb_success_
  name='psb_z_remote_vect'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if (present(dupl)) then 
    dupl_ = dupl
  else
    if (v%is_remote_build()) then
      dupl_ = psb_dupl_add_
    else
      dupl_ = psb_dupl_ovwrt_
    end if
  endif

  ctxt  = desc_a%get_context()
  icomm = desc_a%get_mpic()

  Call psb_info(ctxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Start'
  write(0,*) me, 'X_remote_vect implementation to be completed '

!!$  call b%free() 
!!$
!!$  Allocate(rvsz(np),sdsz(np),sdsi(np),rvsi(np),brvindx(np+1),&
!!$       & bsdindx(np+1), acoo,stat=info)
!!$
!!$  if (info /= psb_success_) then
!!$    info=psb_err_alloc_dealloc_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  end if
!!$
!!$
!!$  nz = a%get_nzeros()
!!$  allocate(ila(nz))
!!$  !write(0,*) me,name,' size :',nz,size(ila)        
!!$  call desc_a%g2l(a%ia(1:nz),ila(1:nz),info,owned=.false.)
!!$  nouth = count(ila(1:nz)<0)
!!$  !write(0,*) me,name,' Count out of halo :',nouth
!!$  call psb_max(ctxt,nouth)
!!$  if ((nouth/=0).and.(me==0)) &
!!$       & write(0,*) 'Warning: would require reinit of DESC_A'
!!$
!!$  call psi_graph_fnd_owner(a%ia(1:nz),iprc,ladj,desc_a%indxmap,info)
!!$  call psb_msort_unique(ladj,nnp)
!!$  !write(0,*) me,name,' Processes:',ladj(1:nnp)
!!$
!!$  icomm = desc_a%get_mpic()
!!$  sdsz(:)=0
!!$  rvsz(:)=0
!!$  sdsi(:)=0
!!$  rvsi(:)=0
!!$  ipx = 1
!!$  brvindx(:) = 0
!!$  bsdindx(:) = 0
!!$  counter=1
!!$  idx = 0
!!$  idxs = 0
!!$  idxr = 0
!!$  do  i=1,nz
!!$    if (iprc(i) >=0) then
!!$      sdsz(iprc(i)+1) = sdsz(iprc(i)+1) +1
!!$    else
!!$      write(0,*)me,name,' Error from fnd_owner: ',iprc(i)
!!$    end if
!!$  end do
!!$  call mpi_alltoall(sdsz,1,psb_mpi_mpk_,& 
!!$       & rvsz,1,psb_mpi_mpk_,icomm,minfo)
!!$  if (minfo /= psb_success_) then
!!$    info=psb_err_from_subroutine_
!!$    call psb_errpush(info,name,a_err='mpi_alltoall')
!!$    goto 9999
!!$  end if
!!$  !write(0,*)me,name,' sdsz ',sdsz(:),' rvsz:',rvsz(:)        
!!$  nsnds = count(sdsz /= 0)
!!$  nrcvs = count(rvsz /= 0)
!!$  idxs = 0
!!$  idxr = 0
!!$  counter = 1
!!$  Do proc=0,np-1
!!$    bsdindx(proc+1) = idxs
!!$    idxs            = idxs + sdsz(proc+1)
!!$    brvindx(proc+1) = idxr
!!$    idxr     = idxr + rvsz(proc+1)
!!$  Enddo
!!$
!!$  iszs = sum(sdsz)
!!$  iszr = sum(rvsz)
!!$  call acoo%allocate(desc_a%get_global_rows(),desc_a%get_global_cols(),iszr)
!!$  if (psb_errstatus_fatal()) then
!!$    write(0,*) 'Error from acoo%allocate '
!!$    info = 4010
!!$    goto 9999
!!$  end if
!!$  if (debug_level >= psb_debug_outer_)&
!!$       & write(debug_unit,*) me,' ',trim(name),': Sizes:',acoo%get_size(),&
!!$       & ' Send:',sdsz(:),' Receive:',rvsz(:)
!!$  !write(debug_unit,*) me,' ',trim(name),': ',info
!!$  if (info == psb_success_) call psb_ensure_size(max(iszs,1),iasnd,info)
!!$  !write(debug_unit,*) me,' ',trim(name),' iasnd: ',info
!!$  if (info == psb_success_) call psb_ensure_size(max(iszs,1),jasnd,info)
!!$  !write(debug_unit,*) me,' ',trim(name),' jasnd: ',info
!!$  if (info == psb_success_) call psb_ensure_size(max(iszs,1),valsnd,info)
!!$  !write(debug_unit,*) me,' ',trim(name),' valsnd: ',info
!!$  if (info /= psb_success_) then
!!$    info=psb_err_from_subroutine_
!!$    call psb_errpush(info,name,a_err='ensure_size')
!!$    goto 9999
!!$  end if
!!$  do k=1, nz
!!$    proc = iprc(k)
!!$    sdsi(proc+1) = sdsi(proc+1) + 1
!!$    !rvsi(proc) = rvsi(proc) + 1 
!!$    iasnd(bsdindx(proc+1)+sdsi(proc+1))  = a%ia(k)
!!$    jasnd(bsdindx(proc+1)+sdsi(proc+1))  = a%ja(k)
!!$    valsnd(bsdindx(proc+1)+sdsi(proc+1)) = a%val(k)
!!$  end do
!!$  do proc=0,np-1
!!$    if (sdsi(proc+1) /= sdsz(proc+1)) &
!!$         & write(0,*) me,name,'Send mismacth ',sdsi(proc+1),sdsz(proc+1)
!!$  end do
!!$
!!$  select case(psb_get_sp_a2av_alg())
!!$  case(psb_sp_a2av_smpl_triad_)
!!$    call psb_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
!!$         & acoo%val,acoo%ia,acoo%ja,rvsz,brvindx,ctxt,info)
!!$  case(psb_sp_a2av_smpl_v_)
!!$    call psb_simple_a2av(valsnd,sdsz,bsdindx,&
!!$         & acoo%val,rvsz,brvindx,ctxt,info)
!!$    if (info == psb_success_) call psb_simple_a2av(iasnd,sdsz,bsdindx,&
!!$         & acoo%ia,rvsz,brvindx,ctxt,info)
!!$    if (info == psb_success_) call psb_simple_a2av(jasnd,sdsz,bsdindx,&
!!$         & acoo%ja,rvsz,brvindx,ctxt,info)
!!$  case(psb_sp_a2av_mpi_)
!!$
!!$    call mpi_alltoallv(valsnd,sdsz,bsdindx,psb_mpi_c_dpk_,&
!!$         & acoo%val,rvsz,brvindx,psb_mpi_c_dpk_,icomm,minfo)
!!$    if (minfo == mpi_success) &
!!$         & call mpi_alltoallv(iasnd,sdsz,bsdindx,psb_mpi_lpk_,&
!!$         & acoo%ia,rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
!!$    if (minfo == mpi_success) &
!!$         & call mpi_alltoallv(jasnd,sdsz,bsdindx,psb_mpi_lpk_,&
!!$         & acoo%ja,rvsz,brvindx,psb_mpi_lpk_,icomm,minfo)
!!$    if (minfo /= mpi_success) info = minfo
!!$  case default
!!$    info = psb_err_internal_error_
!!$    call psb_errpush(info,name,a_err='wrong A2AV alg selector')
!!$    goto 9999
!!$  end select
!!$
!!$  if (info /= psb_success_) then
!!$    info=psb_err_from_subroutine_
!!$    call psb_errpush(info,name,a_err='alltoallv')
!!$    goto 9999
!!$  end if
!!$  call acoo%set_nzeros(iszr)
!!$  call acoo%mv_to_coo(b,info)
!!$  
!!$  Deallocate(brvindx,bsdindx,rvsz,sdsz,&
!!$       & iasnd,jasnd,valsnd,stat=info)
!!$  if (debug_level >= psb_debug_outer_)&
!!$       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

End Subroutine psb_z_remote_vect


