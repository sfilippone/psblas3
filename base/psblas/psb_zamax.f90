!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: psb_zamax.f90
!
! Function: psb_zamax
!    Searches the absolute max of X.
!
!    normi := max(abs(sub(X)(i))  
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Arguments:
!    x(:,:) -  complex               The input vector.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    jx     -  integer(optional).    The column offset.
!
function psb_zamax (x,desc_a, info, jx)
  use psb_penv_mod 
  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  complex(psb_dpk_), intent(in)      :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  integer, optional, intent(in)     :: jx
  real(psb_dpk_)                  :: psb_zamax

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, iix, jjx, ix, ijx, m, imax, izamax
  real(psb_dpk_)         :: amax
  character(len=20)        :: name, ch_err
  double complex    ::         zdum
  double precision  ::  cabs1
  cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )

  name='psb_zamax'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  amax=0.d0

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999 
  endif
  
  ix = 1
  if (present(jx)) then
     ijx = jx
  else
     ijx = 1
  endif

  m = desc_a%get_global_rows()

  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
     info=psb_err_from_subroutine_
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (iix /= 1) then
     info=psb_err_ix_n1_iy_n1_unsupported_
     call psb_errpush(info,name)
     goto 9999
  end if

  ! compute local max
  if ((desc_a%get_local_rows() > 0).and.(m /= 0)) then
     imax=izamax(desc_a%get_local_rows()-iix+1,x(iix,jjx),1)
     amax=cabs1(x(iix+imax-1,jjx))
  else 
    amax = dzero
  end if
  
  ! compute global max
  call psb_amx(ictxt, amax)

  psb_zamax=amax

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end function psb_zamax




!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Function: psb_zamaxv
!    Searches the absolute max of X.
!
!    normi := max(abs(X(i))  
!
! Arguments:
!    x(:)   -  complex               The input vector.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!
function psb_zamaxv (x,desc_a, info)
  use psb_penv_mod
  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  complex(psb_dpk_), intent(in)   :: x(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  real(psb_dpk_)                  :: psb_zamaxv

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, iix, jjx, jx, ix, m, imax, izamax
  real(psb_dpk_)         :: amax
  complex(psb_dpk_)      :: cmax
  character(len=20)        :: name, ch_err
  double complex    ::         zdum
  double precision  ::  cabs1
  cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )

  name='psb_zamaxv'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  amax=0.d0

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  
  ix = 1
  jx = 1

  m = desc_a%get_global_rows()

  call psb_chkvect(m,1,size(x,1),ix,jx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
     info=psb_err_from_subroutine_
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (iix /= 1) then
     info=psb_err_ix_n1_iy_n1_unsupported_
     call psb_errpush(info,name)
     goto 9999
  end if

  ! compute local max
  if ((desc_a%get_local_rows() > 0).and.(m /= 0)) then
     imax=izamax(desc_a%get_local_rows()-iix+1,x(iix),1)
     cmax=(x(iix+imax-1))
     amax=cabs1(cmax)
  else 
    amax = dzero
  end if
  
  ! compute global max
  call psb_amx(ictxt, amax)

  psb_zamaxv=amax

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end function psb_zamaxv

!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Subroutine: psb_zamaxvs
!    Searches the absolute max of X.
!
!    normi := max(abs(sub(X)(i))  
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Arguments:
!    res    -  real                 The result.
!    x(:,:) -  complex              The input vector.
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!    jx     -  integer(optional).   The column offset.
!
subroutine psb_zamaxvs(res,x,desc_a, info)
  use psb_base_mod, psb_protect_name => psb_zamaxvs
  implicit none

  complex(psb_dpk_), intent(in)   :: x(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  real(psb_dpk_), intent(out)     :: res

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, iix, jjx, ix, ijx, m, imax, izamax
  real(psb_dpk_)         :: amax
  character(len=20)        :: name, ch_err
  complex(psb_dpk_)      :: cmax
  double complex    ::         zdum
  double precision  ::  cabs1
  cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )

  name='psb_zamaxvs'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  amax=0.d0

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  ijx=1

  m = desc_a%get_global_rows()
  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
     info=psb_err_from_subroutine_
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (iix /= 1) then
     info=psb_err_ix_n1_iy_n1_unsupported_
     call psb_errpush(info,name)
     goto 9999
  end if

  ! compute local max
  if ((desc_a%get_local_rows() > 0).and.(m /= 0)) then
     imax=izamax(desc_a%get_local_rows()-iix+1,x(iix),1)
     cmax=(x(iix+imax-1))
     amax=cabs1(cmax)
  else 
    amax = dzero
  end if
  
  ! compute global max
  call psb_amx(ictxt, amax)

  res = amax

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end subroutine psb_zamaxvs


!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Subroutine: psb_zmamaxs
!    Searches the absolute max of X.
!
!    normi := max(abs(X(i))  
!
! Arguments:
!    res(:) -  real.                The result.
!    x(:,:) -  complex              The input vector.
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!
subroutine psb_zmamaxs(res,x,desc_a, info,jx)
  use psb_base_mod, psb_protect_name => psb_zmamaxs
  implicit none

  complex(psb_dpk_), intent(in)   :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  integer, optional, intent(in)     :: jx
  real(psb_dpk_), intent(out)     :: res(:)

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, iix, jjx, ix, ijx, m, imax, i, k, izamax
  real(psb_dpk_)         :: amax
  character(len=20)        :: name, ch_err
  complex(psb_dpk_)      :: cmax
  double complex    ::         zdum
  double precision  ::  cabs1
  cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )

  name='psb_zmamaxs'
  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  amax=0.d0

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  
  ix = 1
  if (present(jx)) then
     ijx = jx
  else
     ijx = 1
  endif

  m = desc_a%get_global_rows()
  k  = min(size(x,2),size(res,1))

  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
     info=psb_err_from_subroutine_
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (iix /= 1) then
     info=psb_err_ix_n1_iy_n1_unsupported_
     call psb_errpush(info,name)
     goto 9999
  end if

  ! compute local max
  if ((desc_a%get_local_rows() > 0).and.(m /= 0)) then
     do i=1,k
        imax=izamax(desc_a%get_local_rows()-iix+1,x(iix,jjx+i-1),1)
        cmax=(x(iix+imax-1,jjx+i-1))
        res(i)=cabs1(cmax)
     end do
  else 
    amax = dzero
  end if
  
  ! compute global max
  call psb_amx(ictxt, res(1:k))

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end subroutine psb_zmamaxs
