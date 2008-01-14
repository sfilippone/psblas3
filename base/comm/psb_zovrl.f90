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
!
! File:  psb_zovrl.f90
!
! Subroutine: psb_zovrlm
!   This subroutine performs the exchange of the overlap elements in a 
!    distributed dense matrix between all the processes.
!
! Arguments:
!   x(:,:)      -  complex                   The local part of the dense matrix.
!   desc_a      -  type(psb_desc_type).    The communication descriptor.
!   info        -  integer.                  Return code.
!   jx          -  integer(optional).        The starting column of the global matrix
!   ik          -  integer(optional).        The number of columns to gather. 
!   work        -  real(optional).           A work area.
!   update      -  integer(optional).        Type of update:
!                                            psb_none_   do nothing
!                                            psb_sum_    sum of overlaps
!                                            psb_avg_    average of overlaps
!   mode        -  integer(optional).        Choose the algorithm for data exchange: 
!                                       this is chosen through bit fields. 
!                                       - swap_mpi  = iand(flag,psb_swap_mpi_)  /= 0
!                                       - swap_sync = iand(flag,psb_swap_sync_) /= 0
!                                       - swap_send = iand(flag,psb_swap_send_) /= 0
!                                       - swap_recv = iand(flag,psb_swap_recv_) /= 0
!                                       - if (swap_mpi):  use underlying MPI_ALLTOALLV.
!                                       - if (swap_sync): use PSB_SND and PSB_RCV in 
!                                                       synchronized pairs
!                                       - if (swap_send .and. swap_recv): use mpi_irecv 
!                                                       and mpi_send
!                                       - if (swap_send): use psb_snd (but need another 
!                                                       call with swap_recv to complete)
!                                       - if (swap_recv): use psb_rcv (completing a 
!                                                       previous call with swap_send)
!
!
subroutine  psb_zovrlm(x,desc_a,info,jx,ik,work,update,mode)
  use psb_descriptor_type
  use psb_const_mod
  use psi_mod
  use psb_realloc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(kind(1.d0)), intent(inout), target   :: x(:,:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  complex(kind(1.d0)), optional, target        :: work(:)
  integer, intent(in), optional             :: update,jx,ik,mode

  ! locals
  integer                  :: ictxt, np, me, &
       & err_act, m, n, iix, jjx, ix, ijx, nrow, ncol, k, maxk, update_,&
       & mode_, err, liwork
  complex(kind(1.d0)),pointer :: iwork(:), xp(:,:)
  logical                  :: do_swap
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_zovrlm'
  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  if (present(jx)) then
    ijx = jx
  else
    ijx = 1
  endif

  m = psb_cd_get_global_rows(desc_a)
  n = psb_cd_get_global_cols(desc_a)
  nrow = psb_cd_get_local_rows(desc_a)
  ncol = psb_cd_get_local_cols(desc_a)

  maxk=size(x,2)-ijx+1

  if(present(ik)) then
    if(ik > maxk) then
      k=maxk
    else
      k=ik
    end if
  else
    k = maxk
  end if

  if (present(update)) then 
    update_ = update
  else
    update_ = psb_avg_
  endif

  if (present(mode)) then 
    mode_ = mode
  else
    mode_ = IOR(psb_swap_send_,psb_swap_recv_)
  endif
  do_swap = (mode_ /= 0)

  ! check vector correctness
  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if(info /= 0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix /= 1) then
    info=3040
    call psb_errpush(info,name)
  end if

  err=info
  call psb_errcomm(ictxt,err)
  if(err /= 0) goto 9999

  ! check for presence/size of a work area
  liwork=ncol
  if (present(work)) then
    if(size(work) >= liwork) then
      aliw=.false.
    else
      aliw=.true.
    end if
  else
    aliw=.true.
  end if
  if (aliw) then 
    allocate(iwork(liwork),stat=info)
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if
  else
    iwork => work    
  end if
  ! exchange overlap elements
  if(do_swap) then
    xp => x(iix:size(x,1),jjx:jjx+k-1)
    call psi_swapdata(mode_,k,zone,xp,&
         & desc_a,iwork,info,data=psb_comm_ovr_)
  end if
  if (info == 0) call psi_ovrl_upd(xp,desc_a,update_,info)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='Inner updates')
    goto 9999
  end if

  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zovrlm
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
!
! Subroutine: psb_zovrlv
!   This subroutine performs the exchange of the overlap elements in a 
!    distributed dense vector between all the processes.
!
! Arguments:
!   x(:)        -  complex                   The local part of the dense vector.
!   desc_a      -  type(psb_desc_type).    The communication descriptor.
!   info        -  integer.                  Return code.
!   work        -  real(optional).           A work area.
!   update      -  integer(optional).        Type of update:
!                                            psb_none_   do nothing
!                                            psb_sum_    sum of overlaps
!                                            psb_avg_    average of overlaps
!   mode        -  integer(optional).        Choose the algorithm for data exchange: 
!                                       this is chosen through bit fields. 
!                                       - swap_mpi  = iand(flag,psb_swap_mpi_)  /= 0
!                                       - swap_sync = iand(flag,psb_swap_sync_) /= 0
!                                       - swap_send = iand(flag,psb_swap_send_) /= 0
!                                       - swap_recv = iand(flag,psb_swap_recv_) /= 0
!                                       - if (swap_mpi):  use underlying MPI_ALLTOALLV.
!                                       - if (swap_sync): use PSB_SND and PSB_RCV in 
!                                                       synchronized pairs
!                                       - if (swap_send .and. swap_recv): use mpi_irecv 
!                                                       and mpi_send
!                                       - if (swap_send): use psb_snd (but need another 
!                                                       call with swap_recv to complete)
!                                       - if (swap_recv): use psb_rcv (completing a 
!                                                       previous call with swap_send)
!
!
subroutine  psb_zovrlv(x,desc_a,info,work,update,mode)
  use psb_descriptor_type
  use psi_mod
  use psb_const_mod
  use psb_realloc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(kind(1.d0)), intent(inout), target :: x(:)
  type(psb_desc_type), intent(in)            :: desc_a
  integer, intent(out)                       :: info
  complex(kind(1.d0)), optional, target      :: work(:)
  integer, intent(in), optional              :: update,mode

  ! locals
  integer                  :: ictxt, np, me, &
       & err_act, m, n, iix, jjx, ix, ijx, nrow, ncol, k, update_,&
       & mode_, err, liwork
  complex(kind(1.d0)),pointer :: iwork(:)
  logical                  :: do_swap
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_zovrlv'
  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  ijx = 1

  m = psb_cd_get_global_rows(desc_a)
  n = psb_cd_get_global_cols(desc_a)
  nrow = psb_cd_get_local_rows(desc_a)
  ncol = psb_cd_get_local_cols(desc_a)

  k = 1

  if (present(update)) then 
    update_ = update
  else
    update_ = psb_avg_
  endif

  if (present(mode)) then 
    mode_ = mode
  else
    mode_ = IOR(psb_swap_send_,psb_swap_recv_)
  endif
  do_swap = (mode_ /= 0)

  ! check vector correctness
  call psb_chkvect(m,1,size(x),ix,ijx,desc_a,info,iix,jjx)
  if(info /= 0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix /= 1) then
    info=3040
    call psb_errpush(info,name)
  end if

  err=info
  call psb_errcomm(ictxt,err)
  if(err /= 0) goto 9999

  ! check for presence/size of a work area
  liwork=ncol
  if (present(work)) then
    if(size(work) >= liwork) then
      aliw=.false.
    else
      aliw=.true.
    end if
  else
    aliw=.true.
  end if
  if (aliw) then 
    allocate(iwork(liwork),stat=info)
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if
  else
    iwork => work    
  end if

  ! exchange overlap elements
  if (do_swap) then
    call psi_swapdata(mode_,zone,x(:),&
         & desc_a,iwork,info,data=psb_comm_ovr_)
  end if
  if (info == 0) call psi_ovrl_upd(x,desc_a,update_,info)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='Inner updates')
    goto 9999
  end if
  
  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zovrlv
