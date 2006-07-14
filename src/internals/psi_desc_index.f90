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
subroutine psi_desc_index(desc_data,index_in,dep_list,&
     & length_dl,loc_to_glob,glob_to_loc,desc_index,&
     & isglob_in,info)

  use psb_realloc_mod
  use psb_error_mod
  use psb_const_mod
  use mpi
  use psb_penv_mod
  implicit none

  !c     ...array parameters.....
  integer         :: desc_data(:),index_in(:),dep_list(:)
  integer         :: loc_to_glob(:),glob_to_loc(:)
  integer,pointer :: desc_index(:)
  integer         :: length_dl, info
  logical         :: isglob_in
  !c     ....local scalars...        
  integer :: j,me,np,i,proc,dim
  !c     ...parameters...
  integer :: ictxt
  integer :: no_comm,err
  parameter (no_comm=-1)
  !c     ...local arrays..
  integer :: int_err(5)
  integer,pointer :: brvindx(:),rvsz(:),&
       & bsdindx(:),sdsz(:), sndbuf(:), rcvbuf(:)

  integer :: ihinsz,ntot,k,err_act,&
       & idxr, idxs, iszs, iszr, nesd, nerv, icomm, iret

  logical,parameter :: debug=.false., usempi=.true.
  character(len=20)    :: name, ch_err

  info = 0
  name='psi_desc_index'
  call psb_erractionsave(err_act)

  !c     if mode == 1 then we can use glob_to_loc array
  !c     else we can't utilize it
  ictxt=desc_data(psb_ctxt_)
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

  call psb_get_mpicomm(ictxt,icomm)
  !c 
  !c     first, find out the total sizes to be exchanged.
  !c     note: things marked here as sndbuf/rcvbuf (for mpi) corresponds to things  
  !c     to be received/sent (in the final psblas descriptor).
  !c     be careful of the inversion
  !c   
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
  if ((iszs /= idxs).or.(iszr /= idxr)) then 
    write(0,*) 'strange results???', iszs,idxs,iszr,idxr
  end if
  if (debug) then 
    write(0,*) me,'computed sizes ',iszr,iszs
    call psb_barrier(ictxt)
  endif

  ntot = (3*(max(count(sdsz>0),count(rvsz>0)))+ iszs + iszr) + 1
  if (size(desc_index) < ntot) then 
    !c$$$          write(0,*)  'potential error on desc_index :',
    !c$$$     +      length_dh, size(desc_index),ntot
    write(0,*) 'calling irealloc psi_desc_index ',ntot
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
    !     c
    !     c    note that here bsdinx is zero-based, hence the following loop
    !     c          
    if (isglob_in) then 
      do j=1, nerv
        sndbuf(bsdindx(proc+1)+j) = (index_in(i+j))
      end do
    else
      do j=1, nerv
        sndbuf(bsdindx(proc+1)+j) = loc_to_glob(index_in(i+j))
      end do
    endif
    bsdindx(proc+1) = bsdindx(proc+1) + nerv
    i = i + nerv + 1 
  end do

  if (debug) then 
    write(0,*) me,' prepared send buffer '
    call psb_barrier(ictxt)
  endif
  !c
  !c     now have to regenerate bsdindx
  !c   
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

  !c
  !c   at this point we can finally build the output desc_index. beware
  !c   of snd/rcv inversion. 
  !c       
  i = 1
  do k = 1, length_dl
    proc = dep_list(k)
    desc_index(i) = proc
    i = i + 1 
    nerv = sdsz(proc+1) 
    desc_index(i) = nerv
    do j=1, nerv
      desc_index(i+j) = glob_to_loc(sndbuf(bsdindx(proc+1)+j))
    end do
    i = i + nerv + 1 
    nesd = rvsz(proc+1) 
    desc_index(i) = nesd
    do j=1, nesd
      desc_index(i+j) = glob_to_loc(rcvbuf(brvindx(proc+1)+j))
    end do
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
  if (err_act.eq.act_abort) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_desc_index
