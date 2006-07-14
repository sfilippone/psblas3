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
! File:  psb_dhalo.f90
!
! Subroutine: psb_dhalom
!   This subroutine performs the exchange of the halo elements in a distributed dense matrix between all the processes.
!
! Parameters:
!   x         -  real,dimension(:,:).          The local part of the dense matrix.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   alpha     -  real(optional).               ???.
!   jx        -  integer(optional).            The starting column of the global matrix. 
!   ik        -  integer(optional).            The number of columns to gather. 
!   work      -  real(optional).               A working area.
!   tran      -  character(optional).          ???.
!   mode      -  integer(optional).
!
subroutine  psb_dhalom(x,desc_a,info,alpha,jx,ik,work,tran,mode)
  use psb_descriptor_type
  use psb_const_mod
  use psi_mod
  use psb_check_mod
  use psb_realloc_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(inout), target   :: x(:,:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), intent(in), optional    :: alpha
  real(kind(1.d0)), optional, target        :: work(:)
  integer, intent(in), optional             :: mode,jx,ik
  character, intent(in), optional           :: tran

  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, m, n, iix, jjx, ix, ijx, k, maxk, nrow, imode, i,&
       & err, liwork, ncol
  real(kind(1.d0)),pointer :: iwork(:), xp(:,:)
  character                :: ltran
  character(len=20)        :: name, ch_err

  name='psb_dhalom'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=desc_a%matrix_data(psb_ctxt_)

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

  m = desc_a%matrix_data(psb_m_)
  n = desc_a%matrix_data(psb_n_)
  nrow = desc_a%matrix_data(psb_n_row_)

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

  if (present(tran)) then     
    ltran = tran
  else
    ltran = 'N'
  endif
  if (present(mode)) then 
    imode = mode
  else
    imode = IOR(psb_swap_send_,psb_swap_recv_)
  endif

  ! check vector correctness
  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
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

  if(present(alpha)) then
    if(alpha.ne.1.d0) then
      do i=0, k-1
        call dscal(nrow,alpha,x(1,jjx+i),1)
      end do
    end if
  end if

  liwork=nrow
  if (present(work)) then
    if(size(work).ge.liwork) then
      iwork => work
    else
      call psb_realloc(liwork,iwork,info)
      if(info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if
  else
    call psb_realloc(liwork,iwork,info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if

  ! exchange halo elements
  xp => x(iix:size(x,1),jjx:jjx+k-1)
  if(ltran.eq.'N') then
    call psi_swapdata(imode,k,0.d0,xp,&
         & desc_a,iwork,info,data=psb_comm_halo_)
!!$     call PSI_dSwapData(imode,k,0.d0,x(1,jjx),&
!!$          & size(x,1),desc_a%matrix_data,&
!!$          & desc_a%halo_index,iwork,liwork,info)
  else if((ltran.eq.'T').or.(ltran.eq.'H')) then
    call psi_swaptran(imode,k,1.d0,xp,&
         &desc_a,iwork,info)
!!$     call PSI_dSwapTran(imode,k,1.d0,x(1,jjx),&
!!$          & size(x,1),desc_a%matrix_data,&
!!$          & desc_a%halo_index,iwork,liwork,info)
  end if

  if(info.ne.0) then
    ch_err='PSI_dSwap...'
    call psb_errpush(4010,name,a_err=ch_err)
    goto 9999
  end if

  if(.not.present(work)) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dhalom




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
! Subroutine: psb_dhalov
!   This subroutine performs the exchange of the halo elements in a distributed dense vector between all the processes.
!
! Parameters:
!   x         -  real,dimension(:).            The local part of the dense vector.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   alpha     -  real(optional).               ???.
!   work      -  real(optional).               A working area.
!   tran      -  character(optional).          ???.
!   mode      -  integer(optional).
!
subroutine  psb_dhalov(x,desc_a,info,alpha,work,tran,mode)
  use psb_descriptor_type
  use psb_const_mod
  use psi_mod
  use psb_check_mod
  use psb_realloc_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(inout)           :: x(:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), intent(in), optional    :: alpha
  real(kind(1.d0)), target, optional        :: work(:)
  integer, intent(in), optional             :: mode
  character, intent(in), optional           :: tran

  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, m, n, iix, jjx, ix, ijx, k, maxk, nrow, imode, i,&
       & err, liwork, ncol
  real(kind(1.d0)),pointer :: iwork(:)
  character                :: ltran
  character(len=20)        :: name, ch_err

  name='psb_dhalov'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  ijx = 1

  m = desc_a%matrix_data(psb_m_)
  n = desc_a%matrix_data(psb_n_)
  nrow = desc_a%matrix_data(psb_n_row_)

  if (present(tran)) then     
    ltran = tran
  else
    ltran = 'N'
  endif
  if (present(mode)) then 
    imode = mode
  else
    imode = IOR(psb_swap_send_,psb_swap_recv_)
  endif

  ! check vector correctness
  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
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

  if(present(alpha)) then
    if(alpha.ne.1.d0) then
      call dscal(nrow,alpha,x,ione)
    end if
  end if

  liwork=nrow
  if (present(work)) then
    if(size(work).ge.liwork) then
      iwork => work
    else
      call psb_realloc(liwork,iwork,info)
      if(info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if
  else
    call psb_realloc(liwork,iwork,info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if

  ! exchange halo elements
  if(ltran.eq.'N') then
    call psi_swapdata(imode,0.d0,x(iix:size(x)),&
         & desc_a,iwork,info,data=psb_comm_halo_)
  else if((ltran.eq.'T').or.(ltran.eq.'H')) then
    call psi_swaptran(imode,1.d0,x(iix:size(x)),&
         & desc_a,iwork,info)
  end if

  if(info.ne.0) then
    ch_err='PSI_dSwap...'
    call psb_errpush(4010,name,a_err=ch_err)
    goto 9999
  end if

  if(.not.present(work)) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dhalov



