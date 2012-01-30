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
! File: psb_damax.f90
!
! Function: psb_damax
!    Searches the absolute max of X.
!
!    normi := max(abs(sub(X)(i))  
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Arguments:
!    x(:,:) -  real               The input vector.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    jx     -  integer(optional).    The column offset.
!
function psb_damax(x,desc_a, info, jx) result(res)
  use psb_base_mod, psb_protect_name => psb_damax

  implicit none

  real(psb_dpk_), intent(in)    :: x(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_), optional, intent(in)    :: jx
  real(psb_dpk_)                   :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ix, ijx, m, ldx
  character(len=20)      :: name, ch_err

  name='psb_damax'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = desc_a%get_context()

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
  ldx = size(x,1)

  call psb_chkvect(m,ione,ldx,ix,ijx,desc_a,info,iix,jjx)
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
    res = psb_amax(desc_a%get_local_rows()-iix+1,x(:,jjx))
  else 
    res = dzero
  end if
  
  ! compute global max
  call psb_amx(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end function psb_damax




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
! Function: psb_damaxv
!    Searches the absolute max of X.
!
!    normi := max(abs(X(i))  
!
! Arguments:
!    x(:)   -  real               The input vector.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!
function psb_damaxv (x,desc_a, info) result(res)
  use psb_base_mod, psb_protect_name => psb_damaxv

  implicit none

  real(psb_dpk_), intent(in)   :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info
  real(psb_dpk_)                  :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, jx, ix, m, ldx

  character(len=20)        :: name, ch_err

  name='psb_damaxv'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


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
  ldx = size(x,1)

  call psb_chkvect(m,ione,ldx,ix,jx,desc_a,info,iix,jjx)
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
    res = psb_amax(desc_a%get_local_rows()-iix+1,x)
  else 
    res = dzero
  end if
  
  ! compute global max
  call psb_amx(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end function psb_damaxv


function psb_damax_vect(x, desc_a, info) result(res)
  use psb_penv_mod
  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_d_vect_mod
  implicit none

  real(psb_dpk_)                        :: res
  type(psb_d_vect_type), intent (inout) :: x
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)         :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, jx, ix, m
  character(len=20)      :: name, ch_err

  name='psb_damaxv'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  jx = 1

  m = desc_a%get_global_rows()
  call psb_chkvect(m,ione,x%get_nrows(),ix,jx,desc_a,info,iix,jjx)
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
    res = x%amax(desc_a%get_local_rows())
  else 
    res = dzero
  end if

  ! compute global max
  call psb_amx(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end function psb_damax_vect


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
! Subroutine: psb_damaxvs
!    Searches the absolute max of X.
!
!    normi := max(abs(sub(X)(i))  
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Arguments:
!    res    -  real                 The result.
!    x(:,:) -  real              The input vector.
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!    jx     -  integer(optional).   The column offset.
!
subroutine psb_damaxvs(res,x,desc_a, info)
  use psb_base_mod, psb_protect_name => psb_damaxvs

  implicit none

  real(psb_dpk_), intent(in)   :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info
  real(psb_dpk_), intent(out)      :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ix, ijx, m, ldx
  character(len=20)      :: name, ch_err

  name='psb_damaxvs'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)


  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  ijx=1

  m = desc_a%get_global_rows()
  ldx=size(x,1)
  call psb_chkvect(m,ione,ldx,ix,ijx,desc_a,info,iix,jjx)
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
    res = psb_amax(desc_a%get_local_rows()-iix+1,x)
  else 
    res = dzero
  end if
  
  ! compute global max
  call psb_amx(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end subroutine psb_damaxvs


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
! Subroutine: psb_dmamaxs
!    Searches the absolute max of X.
!
!    normi := max(abs(X(i))  
!
! Arguments:
!    res(:) -  real.                The result.
!    x(:,:) -  real              The input vector.
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!
subroutine psb_dmamaxs(res,x,desc_a, info,jx)
  use psb_base_mod, psb_protect_name => psb_dmamaxs

  implicit none

  real(psb_dpk_), intent(in)   :: x(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), optional, intent(in)   :: jx
  real(psb_dpk_), intent(out)     :: res(:)

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ix, ijx, m, ldx, i, k
  character(len=20)        :: name, ch_err

  name='psb_dmamaxs'
  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

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
  ldx = size(x,1)
  call psb_chkvect(m,ione,ldx,ix,ijx,desc_a,info,iix,jjx)
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

  res(1:k) = dzero
  ! compute local max
  if ((desc_a%get_local_rows() > 0).and.(m /= 0)) then
    do i=1,k
      res(i) = psb_amax(desc_a%get_local_rows()-iix+1,x(:,jjx+i-1))
    end do
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
end subroutine psb_dmamaxs
