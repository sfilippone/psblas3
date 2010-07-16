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
!
!
! File: psi_desc_index.f90
!
! Subroutine: psi_desc_index
!    Converts a list of data exchanges from build format to assembled format. 
!    See below for a description of the formats.
!
! Arguments:
! desc_a       - type(psb_desc_type)   The descriptor; in this context only the index 
!                                       mapping parts are used.
! index_in(:)  - integer               The index list, build format  
! index_out(:) - integer, allocatable  The index list, assembled format
! glob_idx     - logical               Whether the input indices are in local or global
!                                      numbering; the global numbering is used when 
!                                      converting the overlap exchange lists.
! nxch         - integer               The number of data exchanges on the calling process
! nsnd         - integer               Total send buffer size       on the calling process
! nrcv         - integer               Total receive buffer size    on the calling process
!
! The format of the index lists. Copied from base/modules/psb_desc_type
!
!  7. The data exchange is based on lists of local indices to be exchanged; all the 
!     lists have the same format, as follows:
!     the list is  composed of variable dimension blocks, one for each process to 
!     communicate with; each block contains indices of local elements to be 
!     exchanged. We do choose the order of communications:  do not change 
!     the sequence of blocks unless you know what you're doing, or you'll 
!     risk a deadlock. NOTE: This is the format when the state is PSB_ASB_.
!     See below for BLD. The end-of-list is marked with a -1. 
!
!  notation        stored in		          explanation
!  --------------- --------------------------- -----------------------------------
!  process_id      index_v(p+proc_id_)      identifier of process with which 
!                                                data is  exchanged.
!  n_elements_recv index_v(p+n_elem_recv_)  number of elements to receive.
!  elements_recv   index_v(p+elem_recv_+i)  indexes of local elements to
!					          receive. these are stored in the
!					          array from location p+elem_recv_ to
!					          location p+elem_recv_+
!						  index_v(p+n_elem_recv_)-1.
!  n_elements_send index_v(p+n_elem_send_)  number of elements to send.
!  elements_send   index_v(p+elem_send_+i)  indexes of local elements to
!					          send. these are stored in the
!					          array from location p+elem_send_ to
!					          location p+elem_send_+
!						  index_v(p+n_elem_send_)-1.
!
!     This organization is valid for both halo and overlap indices; overlap entries
!     need to be updated to ensure that a variable at a given global index 
!     (assigned to multiple processes) has the same value. The way to resolve the 
!     issue is to exchange the data and then sum (or average) the values. See
!     psb_ovrl subroutine. 
!  
!  8. When the descriptor is in the BLD state the INDEX vectors contains only 
!     the indices to be received, organized as  a sequence 
!     of entries of the form (proc,N,(lx1,lx2,...,lxn)) with owning process,
!     number of indices (most often N=1), list of local indices. 
!     This is because we only know the list of halo indices to be received 
!     as we go about building the sparse matrix pattern, and we want the build 
!     phase to be loosely synchronized. Thus we record the indices we have to ask 
!     for, and at the time we call PSB_CDASB we match all the requests to figure 
!     out who should be sending what to whom.
!     However this implies that we know who owns the indices; if we are in the 
!     LARGE case (as described above) this is actually only true for the OVERLAP list 
!     that is filled in at CDALL time, and not for the HALO; thus the HALO list 
!     is rebuilt during the CDASB process (in the psi_ldsc_pre_halo subroutine). 
!
!
subroutine psi_desc_index(desc,index_in,dep_list,&
     & length_dl,nsnd,nrcv,desc_index,isglob_in,info)
  use psb_descriptor_type
  use psb_realloc_mod
  use psb_error_mod
  use psb_const_mod
#ifdef MPI_MOD
  use mpi
#endif
  use psb_penv_mod
  use psi_mod, psb_protect_name => psi_desc_index
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  !    ...array parameters.....
  type(psb_desc_type) :: desc
  integer         :: index_in(:),dep_list(:)
  integer,allocatable  :: desc_index(:)
  integer         :: length_dl,nsnd,nrcv,info
  logical         :: isglob_in
  !    ....local scalars...        
  integer :: j,me,np,i,proc
  !    ...parameters...
  integer :: ictxt
  integer, parameter  :: no_comm=-1
  !     ...local arrays..
  integer,allocatable  :: brvindx(:),rvsz(:),&
       & bsdindx(:),sdsz(:), sndbuf(:), rcvbuf(:)

  integer :: ihinsz,ntot,k,err_act,nidx,&
       & idxr, idxs, iszs, iszr, nesd, nerv, icomm

  logical,parameter :: usempi=.true.
  integer              :: debug_level, debug_unit
  character(len=20) :: name

  info = psb_success_
  name='psi_desc_index'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = psb_cd_get_context(desc)
  icomm = psb_cd_get_mpic(desc)
  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (debug_level >= psb_debug_inner_) then 
    write(debug_unit,*) me,' ',trim(name),': start'
    call psb_barrier(ictxt)
  endif

  ! 
  !     first, find out the sizes to be exchanged.
  !     note: things marked here as sndbuf/rcvbuf (for mpi) corresponds to things  
  !     to be received/sent (in the final psblas descriptor).
  !     be careful of the inversion
  !   
  allocate(sdsz(np),rvsz(np),bsdindx(np),brvindx(np),stat=info)
  if(info /= psb_success_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  sdsz(:) = 0
  rvsz(:) = 0
  bsdindx(:) = 0
  brvindx(:) = 0
  i = 1
  do 
    if (index_in(i) == -1) exit
    proc = index_in(i)
    i = i + 1 
    nerv = index_in(i)
    sdsz(proc+1) = sdsz(proc+1) + nerv
    i = i + nerv + 1 
  end do
  ihinsz=i
  call mpi_alltoall(sdsz,1,mpi_integer,rvsz,1,mpi_integer,icomm,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mpi_alltoall')
    goto 9999
  end if

  i    = 1
  idxs = 0
  idxr = 0
  do i=1, length_dl
    proc = dep_list(i)
    bsdindx(proc+1) = idxs
    idxs = idxs + sdsz(proc+1)
    brvindx(proc+1) = idxr
    idxr = idxr + rvsz(proc+1)
  end do
  iszs = sum(sdsz)
  iszr = sum(rvsz)
  nsnd = iszr
  nrcv = iszs  

  if ((iszs /= idxs).or.(iszr /= idxr)) then 
    write(psb_err_unit,*) me, trim(name),': Warning: strange results?', &
         & iszs,idxs,iszr,idxr
  end if
  if (debug_level >= psb_debug_inner_) then 
    write(debug_unit,*) me,' ',trim(name),': computed sizes ',iszr,iszs
    call psb_barrier(ictxt)
  endif

  ntot = (3*(count((sdsz>0).or.(rvsz>0)))+ iszs + iszr) + 1

  if (ntot > psb_size(desc_index)) then 
    call psb_realloc(ntot,desc_index,info) 
  endif
!!$  call psb_ensure_size(ntot,desc_index,info)

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_realloc')
    goto 9999
  end if

  if (debug_level >= psb_debug_inner_) then 
    write(debug_unit,*) me,' ',trim(name),': computed allocated workspace ',iszr,iszs
    call psb_barrier(ictxt)
  endif
  allocate(sndbuf(iszs),rcvbuf(iszr),stat=info)
  if(info /= psb_success_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  ! 
  !     Second build the lists of requests
  !   
  i = 1
  do 
    if (i > ihinsz) then 
!!$      write(psb_err_unit,*) me,' did not find index_in end??? ',i,ihinsz
      exit
    end if
    if (index_in(i) == -1) exit
    proc = index_in(i)
    i = i + 1 
    nerv = index_in(i)
    !  
    !   note that here bsdinx is zero-based, hence the following loop
    !          
    if (isglob_in) then 
      do j=1, nerv
        sndbuf(bsdindx(proc+1)+j) = (index_in(i+j))
      end do
    else
      call psb_map_l2g(index_in(i+1:i+nerv),&
           & sndbuf(bsdindx(proc+1)+1:bsdindx(proc+1)+nerv),&
           & desc%idxmap,info) 
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_map_l2g')
        goto 9999
      end if

    endif
    bsdindx(proc+1) = bsdindx(proc+1) + nerv
    i = i + nerv + 1 
  end do

  if (debug_level >= psb_debug_inner_) then 
    write(debug_unit,*) me,' ',trim(name),': prepared send buffer '
    call psb_barrier(ictxt)
  endif
  !
  !   now have to regenerate bsdindx
  !  
  idxs = 0
  idxr = 0
  do i=1, length_dl
    proc = dep_list(i)
    bsdindx(proc+1) = idxs
    idxs = idxs + sdsz(proc+1)
    brvindx(proc+1) = idxr
    idxr = idxr + rvsz(proc+1)
  end do

  call mpi_alltoallv(sndbuf,sdsz,bsdindx,mpi_integer,&
       & rcvbuf,rvsz,brvindx,mpi_integer,icomm,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='mpi_alltoallv')
    goto 9999
  end if

  !
  !  at this point we can finally build the output desc_index. beware
  !  of snd/rcv inversion. 
  !      
  i = 1
  do k = 1, length_dl
    proc = dep_list(k)
    desc_index(i) = proc
    i = i + 1 
    nerv = sdsz(proc+1) 
    desc_index(i) = nerv
    call psi_idx_cnv(nerv,sndbuf(bsdindx(proc+1)+1:bsdindx(proc+1)+nerv),&
         &  desc_index(i+1:i+nerv),desc,info)
    i = i + nerv + 1 
    nesd = rvsz(proc+1) 
    desc_index(i) = nesd
    call psi_idx_cnv(nesd,rcvbuf(brvindx(proc+1)+1:brvindx(proc+1)+nesd),&
         &  desc_index(i+1:i+nesd),desc,info)
    i = i + nesd + 1 
  end do
  desc_index(i) = - 1 

  deallocate(sdsz,rvsz,bsdindx,brvindx,sndbuf,rcvbuf,stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug_level >= psb_debug_inner_) then 
    write(debug_unit,*) me,' ',trim(name),': done'
    call psb_barrier(ictxt)
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_desc_index
