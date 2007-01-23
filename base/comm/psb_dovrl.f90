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
! File:  psb_dovrl.f90
!
! Subroutine: psb_dovrlm
!   This subroutine performs the exchange of the overlap elements in a 
!   distributed dense matrix between all the processes.
!
! Parameters:
!   x           -  real,dimension(:,:).         The local part of the dense matrix.
!   desc_a      -  type(<psb_desc_type>).       The communication descriptor.
!   info        -  integer.                     A return code.
!   jx          -  integer(optional).           The starting column of the global matrix
!   ik          -  integer(optional).           The number of columns to gather. 
!   work        -  real(optional).              A work area.
!   update -  integer(optional).            ???.
!
subroutine  psb_dovrlm(x,desc_a,info,jx,ik,work,update)
  use psb_descriptor_type
  use psb_const_mod
  use psi_mod
  use psb_realloc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(inout), target   :: x(:,:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), optional, target        :: work(:)
  integer, intent(in), optional             :: update,jx,ik

  ! locals
  integer                  :: int_err(5), ictxt, np, me, &
       & err_act, m, n, iix, jjx, ix, ijx, nrow, ncol, k, maxk, iupdate,&
       & imode, err, liwork, i
  real(kind(1.d0)),pointer :: iwork(:), xp(:,:)
  logical                  :: do_update
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_dovrlm'
  if(psb_get_errstatus().ne.0) return 
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
    if(ik.gt.maxk) then
      k=maxk
    else
      k=ik
    end if
  else
    k = maxk
  end if

  if (present(update)) then 
    iupdate = update
  else
    iupdate = psb_avg_
  endif

  do_update = (iupdate /= psb_none_)
  imode = IOR(psb_swap_send_,psb_swap_recv_)

  ! check vector correctness
  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix.ne.1) then
    info=3040
    call psb_errpush(info,name)
  end if

  err=info
  call psb_errcomm(ictxt,err)
  if(err.ne.0) goto 9999

  ! check for presence/size of a work area
  liwork=ncol
  if (present(work)) then
    if(size(work).ge.liwork) then
      iwork => work
      aliw=.false.
    else
      aliw=.true.
      allocate(iwork(liwork),stat=info)
      if(info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if
  else
    aliw=.true.
    allocate(iwork(liwork),stat=info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if
    
  ! exchange overlap elements
  if(do_update) then
    xp => x(iix:size(x,1),jjx:jjx+k-1)
    call psi_swapdata(imode,k,1.d0,xp,&
         & desc_a,iwork,info,data=psb_comm_ovr_)
  end if

  if(info.ne.0) then
    call psb_errpush(4010,name,a_err='psi_swapdata')
    goto 9999
  end if

  i=1
  ! switch on update type
  select case (iupdate)
  case(psb_square_root_)
    do while(desc_a%ovrlap_elem(i).ne.-ione)
      x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:) =&
           & x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:)/&
           & sqrt(real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_)))
      i = i+2
    end do
  case(psb_avg_)
    do while(desc_a%ovrlap_elem(i).ne.-ione)
      x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:) =&
           & x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:)/&
           & real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_))
      i = i+2
    end do
  case(psb_sum_)
    ! do nothing
  case default 
    ! wrong value for choice argument
    info = 70
    int_err=(/10,iupdate,0,0,0/)
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end select

  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dovrlm






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
! Subroutine: psb_dovrlv
!   This subroutine performs the exchange of the overlap elements in a 
!   distributed dense vector between all the processes.
!
! Parameters:
!   x           -  real,dimension(:).          The local part of the dense vector.
!   desc_a      -  type(<psb_desc_type>).        The communication descriptor.
!   info        -  integer.                      Eventually returns an error code.
!   work        -  real(optional).               A working area.
!   update -  integer(optional).            ???.
!
subroutine  psb_dovrlv(x,desc_a,info,work,update)
  use psb_descriptor_type
  use psi_mod
  use psb_const_mod
  use psb_realloc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(inout), target   :: x(:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), optional, target        :: work(:)
  integer, intent(in), optional             :: update

  ! locals
  integer                  :: int_err(5), ictxt, np, me, &
       & err_act, m, n, iix, jjx, ix, ijx, nrow, ncol, k, iupdate,&
       & imode, err, liwork, i
  real(kind(1.d0)),pointer :: iwork(:)
  logical                  :: do_update
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_dovrlv'
  if(psb_get_errstatus().ne.0) return 
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
    iupdate = update
  else
    iupdate = psb_avg_
  endif

  do_update = (iupdate /= psb_none_)
  imode = IOR(psb_swap_send_,psb_swap_recv_)

  ! check vector correctness
  call psb_chkvect(m,1,size(x),ix,ijx,desc_a,info,iix,jjx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix.ne.1) then
    info=3040
    call psb_errpush(info,name)
  end if

  err=info
  call psb_errcomm(ictxt,err)
  if(err.ne.0) goto 9999

  ! check for presence/size of a work area
  liwork=ncol
  if (present(work)) then
    if(size(work).ge.liwork) then
      iwork => work
      aliw=.false.
    else
      aliw=.true.
      allocate(iwork(liwork),stat=info)
      if(info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if
  else
    aliw=.true.
    allocate(iwork(liwork),stat=info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if

  ! exchange overlap elements
  if(do_update) then
    call psi_swapdata(imode,1.d0,x(iix:size(x)),&
         & desc_a,iwork,info,data=psb_comm_ovr_)
  end if

  if(info.ne.0) then
    call psb_errpush(4010,name,a_err='PSI_SwapData')
    goto 9999
  end if

  i=1
  ! switch on update type
  select case (iupdate)
  case(psb_square_root_)
    do while(desc_a%ovrlap_elem(i).ne.-ione)
      x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_)) =&
           & x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_))/&
           & sqrt(real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_)))
      i = i+2
    end do
  case(psb_avg_)
    do while(desc_a%ovrlap_elem(i).ne.-ione)
      x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_)) =&
           & x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_))/&
           & real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_))
      i = i+2
    end do
  case(psb_sum_)
    ! do nothing
  case default 
    ! wrong value for choice argument
    info = 70
    int_err=(/10,iupdate,0,0,0/)
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end select

  if (aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dovrlv
