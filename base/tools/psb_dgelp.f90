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
! File: psb_dgelp.f90
!
! Subroutine: psb_dgelp
!    ???????????
!
! Parameters:
! trans    - character. 
! iperm    - integer.
! x        - real, dimension(:,:).
! info     - integer.                 Eventually returns an error code.
subroutine psb_dgelp(trans,iperm,x,desc_a,info)
  !....assembly dense matrix x .....
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_psblas_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  type(psb_desc_type), intent(in)      ::  desc_a
  real(kind(1.d0)), intent(inout)      ::  x(:,:)
  integer, intent(inout)               ::  iperm(:),info
  character, intent(in)                :: trans

  ! local variables
  integer                  :: ictxt,np, me,nrow,ncol
  real(kind(1.d0)),pointer :: dtemp(:)
  integer                  :: int_err(5), i1sz, i2sz, dectype, i, err_act
  real(kind(1.d0)),parameter    :: one=1
  logical, parameter :: debug=.false.

  interface dgelp
    subroutine dgelp(trans,m,n,p,b,ldb,work,lwork,ierror)
      integer, intent(in)  :: ldb, m, n, lwork
      integer, intent(out) :: ierror
      character, intent(in) :: trans
      double precision, intent(inout) ::  b(ldb,*), work(*)
      integer, intent(in)  :: p(*)
    end subroutine dgelp
  end interface

  interface isaperm

    logical function isaperm(n,ip)
      integer, intent(in)    :: n   
      integer, intent(inout) :: ip(*)
    end function isaperm
  end interface

  character(len=20)   :: name, ch_err
  name = 'psb_dgelp'

  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt   = psb_cd_get_context(desc_a)
  dectype = psb_cd_get_dectype(desc_a)
  nrow    = psb_cd_get_local_rows(desc_a)
  ncol    = psb_cd_get_local_cols(desc_a)
  i1sz    = size(x,dim=1)
  i2sz    = size(x,dim=2)

  call psb_info(ictxt, me, np)

  if (debug) write(*,*) 'asb start: ',np,me,&
       & psb_cd_get_dectype(desc_a)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (.not.psb_is_asb_desc(desc_a)) then
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif


  if (.not.isaperm(i1sz,iperm)) then
    info = 70
    int_err(1) = 1      
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  if (debug) write(*,*) 'asb: ',i1sz,i2sz,nrow,ncol
  allocate(dtemp(i1sz),stat=info)

  call dgelp(trans,i1sz,i2sz,iperm,x,i1sz,dtemp,i1sz,info)
  if(info /= 0) then
    info=4010
    ch_err='dgelp'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  deallocate(dtemp)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == act_ret) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_dgelp



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
! Subroutine: psb_dgelpv
!    ???????????
!
! Parameters:
! trans    - character. 
! iperm    - integer.
! x        - real, dimension(:).
! info     - integer.                 Eventually returns an error code.
subroutine psb_dgelpv(trans,iperm,x,desc_a,info)
  !....assembly dense matrix x .....
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_psblas_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  type(psb_desc_type), intent(in)    ::  desc_a
  real(kind(1.d0)), intent(inout)    ::  x(:)
  integer, intent(inout)             ::  iperm(:), info
  character, intent(in)              ::  trans

  ! local variables
  integer :: ictxt,np,me
  integer :: int_err(5), i1sz,nrow,ncol,dectype, err_act
  real(kind(1.d0)),pointer ::  dtemp(:)
  real(kind(1.d0)),parameter    :: one=1
  logical, parameter :: debug=.false.

  interface dgelp
    subroutine dgelp(trans,m,n,p,b,ldb,work,lwork,ierror)
      integer, intent(in)  :: ldb, m, n, lwork
      integer, intent(out) :: ierror
      character, intent(in) :: trans
      double precision, intent(inout) ::  b(*), work(*)
      integer, intent(in)  :: p(*)
    end subroutine dgelp
  end interface

  interface isaperm

    logical function isaperm(n,ip)
      integer, intent(in)    :: n   
      integer, intent(inout) :: ip(*)
    end function isaperm
  end interface

  character(len=20)   :: name, ch_err
  name = 'psb_dgelpv'

  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  i1sz = size(x)

  ictxt   = psb_cd_get_context(desc_a)
  dectype = psb_cd_get_dectype(desc_a)
  nrow    = psb_cd_get_local_rows(desc_a)
  ncol    = psb_cd_get_local_cols(desc_a)

  call psb_info(ictxt, me, np)

  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (.not.psb_is_asb_desc(desc_a)) then
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif

  if (debug) write(0,*) 'calling isaperm ',i1sz,size(iperm),trans

  if (.not.isaperm(i1sz,iperm)) then
    info = 70
    int_err(1) = 1      
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  allocate(dtemp(i1sz),stat=info)

  call dgelp(trans,i1sz,1,iperm,x,i1sz,dtemp,i1sz,info)
  if(info /= 0) then
    info=4010
    ch_err='dgelp'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  deallocate(dtemp)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == act_ret) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_dgelpv

