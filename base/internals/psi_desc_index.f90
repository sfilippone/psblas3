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
subroutine psi_desc_index(desc,index_in,dep_list,&
     & length_dl,nsnd,nrcv,desc_index,isglob_in,info)
  use psb_descriptor_type
  use psb_realloc_mod
  use psb_error_mod
  use psb_const_mod
  use mpi
  use psb_penv_mod
  use psi_mod, only : psi_idx_cnv
  implicit none

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
  integer :: no_comm,err
  parameter (no_comm=-1)
  !     ...local arrays..
  integer,allocatable  :: brvindx(:),rvsz(:),&
       & bsdindx(:),sdsz(:), sndbuf(:), rcvbuf(:)

  integer :: ihinsz,ntot,k,err_act,nidx,&
       & idxr, idxs, iszs, iszr, nesd, nerv, icomm

  logical,parameter :: debug=.false., usempi=.true.
  character(len=20)  :: name

  info = 0
  name='psi_desc_index'
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc)
  icomm = psb_cd_get_mpic(desc)
  call psb_info(ictxt,me,np) 
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (debug) then 
    write(0,*) me,'start desc_index'
    call psb_barrier(ictxt)
  endif

  ! 
  !     first, find out the total sizes to be exchanged.
  !     note: things marked here as sndbuf/rcvbuf (for mpi) corresponds to things  
  !     to be received/sent (in the final psblas descriptor).
  !     be careful of the inversion
  !   
  allocate(sdsz(np),rvsz(np),bsdindx(np),brvindx(np),stat=info)
  if(info /= 0) then
    info=4000
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
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='mpi_alltoall')
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
    write(0,*) 'strange results???', iszs,idxs,iszr,idxr
  end if
  if (debug) then 
    write(0,*) me,'computed sizes ',iszr,iszs
    call psb_barrier(ictxt)
  endif

  ntot = (3*(count((sdsz>0).or.(rvsz>0)))+ iszs + iszr) + 1
  if (allocated(desc_index)) then 
    nidx = size(desc_index)
  else
    nidx = 0 
  endif

  if (nidx < ntot) then 
    call psb_realloc(ntot,desc_index,info)
  endif
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='psb_realloc')
    goto 9999
  end if

  if (debug) then 
    write(0,*) me,'computed allocated workspace ',iszr,iszs
    call psb_barrier(ictxt)
  endif
  allocate(sndbuf(iszs),rcvbuf(iszr),stat=info)
  if(info /= 0) then
    info=4000
    call psb_errpush(info,name)
    goto 9999
  end if

  i = 1
  do 
    if (i > ihinsz) then 
      write(0,*) me,' did not find index_in end??? ',i,ihinsz
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
      
      do j=1, nerv
        sndbuf(bsdindx(proc+1)+j) = desc%loc_to_glob(index_in(i+j))
      end do
    endif
    bsdindx(proc+1) = bsdindx(proc+1) + nerv
    i = i + nerv + 1 
  end do

  if (debug) then 
    write(0,*) me,' prepared send buffer '
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
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='mpi_alltoallv')
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
!!$    do j=1, nerv
!!$      desc_index(i+j) = glob_to_loc(sndbuf(bsdindx(proc+1)+j))
!!$    end do
    i = i + nerv + 1 
    nesd = rvsz(proc+1) 
    desc_index(i) = nesd
    call psi_idx_cnv(nesd,rcvbuf(brvindx(proc+1)+1:brvindx(proc+1)+nesd),&
         &  desc_index(i+1:i+nesd),desc,info)
!!$    do j=1, nesd
!!$      desc_index(i+j) = glob_to_loc(rcvbuf(brvindx(proc+1)+j))
!!$    end do
    i = i + nesd + 1 
  end do
  desc_index(i) = - 1 

  deallocate(sdsz,rvsz,bsdindx,brvindx,sndbuf,rcvbuf,stat=info)
  if (info /= 0) then 
    info=4000
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug) then 
    write(0,*) me,'end desc_index'
    call psb_barrier(ictxt)
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_desc_index
