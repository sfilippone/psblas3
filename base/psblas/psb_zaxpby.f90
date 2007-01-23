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
! File: psb_daxpby.f90
!
! Subroutine: psb_daxpby
!    Adds one distributed matrix to another,
!
!    sub( Y ) := beta * sub( Y ) + alpha * sub( X )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Parameters:
!    alpha  -  real.                      The scalar used to multiply each component of sub( X ).
!    x      -  real,dimension(:,:).       The input vector containing the entries of sub( X ).
!    beta   -  real.                      The scalar used to multiply each component of sub( Y ).
!    y      -  real,dimension(:,:).       The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset for sub( X ).
!    jy     -  integer(optional).         The column offset for sub( Y ).
!
subroutine  psb_zaxpby(alpha, x, beta,y,desc_a,info, n, jx, jy)
  use psb_descriptor_type
  use psb_check_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none                    

  integer, intent(in), optional   :: n, jx, jy
  integer, intent(out)            :: info
  type(psb_desc_type), intent(in) :: desc_a
  complex(kind(1.D0)), intent(in)    :: alpha, beta
  complex(kind(1.D0)), intent(in)    :: x(:,:)
  complex(kind(1.D0)), intent(inout) :: y(:,:)

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, iix, jjx, ix, iy, ijx, ijy, m, iiy, in, jjy
  character(len=20)        :: name, ch_err

  name='psb_dgeaxpby'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -ione) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = ione
  if (present(jx)) then
    ijx = jx
  else
    ijx = ione
  endif

  iy = ione
  if (present(jy)) then
    ijy = jy
  else
    ijy = ione
  endif

  if (present(n)) then
    if(((ijx+n).le.size(x,2)).and.&
         & ((ijy+n).le.size(y,2))) then 
      in = n
    else
      in = min(size(x,2),size(y,2))
    end if
  else
    in = min(size(x,2),size(y,2))
  endif

  if(ijx.ne.ijy) then
    info=3050
    call psb_errpush(info,name)
    goto 9999
  end if

  m = psb_cd_get_global_rows(desc_a)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if (info == 0) &
       & call psb_chkvect(m,ione,size(y,1),iy,ijy,desc_a,info,iiy,jjy)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((iix.ne.ione).or.(iiy.ne.ione)) then
    info=3040
    call psb_errpush(info,name)
    goto 9999
  end if

  if ((in.ne.0)) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      call zaxpby(psb_cd_get_local_cols(desc_a),in,&
           & alpha,x(iix,jjx),size(x,1),beta,&
           & y(iiy,jjy),size(y,1),info)
    end if
  end if

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zaxpby





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
! Subroutine: psb_dgeaxpbyv
!    Adds one distributed matrix to another,
!
!    Y := beta * Y + alpha * X
!
! Parameters:
!    alpha  -  real.                      The scalar used to multiply each component of X.
!    x      -  real,dimension(:).         The input vector containing the entries of X.
!    beta   -  real.                      The scalar used to multiply each component of Y.
!    y      -  real,dimension(:).         The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
subroutine  psb_zaxpbyv(alpha, x, beta,y,desc_a,info)
  use psb_descriptor_type
  use psb_const_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none                    

  integer, intent(out)            :: info
  type(psb_desc_type), intent(in) :: desc_a
  complex(kind(1.D0)), intent(in)    :: alpha, beta
  complex(kind(1.D0)), intent(in)    :: x(:)
  complex(kind(1.D0)), intent(inout) :: y(:)

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, n, iix, jjx, ix, iy, ijx, m, iiy, in, jjy
  character(len=20)        :: name, ch_err
  logical, parameter :: debug=.false.

  name='psb_dgeaxpby'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)
  if (np == -ione) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = psb_cd_get_global_rows(desc_a)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x),ix,ione,desc_a,info,iix,jjx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect 1'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,ione,size(y),iy,ione,desc_a,info,iiy,jjy)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect 2'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((iix.ne.ione).or.(iiy.ne.ione)) then
    info=3040
    call psb_errpush(info,name)
  end if

  if(psb_cd_get_local_rows(desc_a).gt.0) then
    call zaxpby(psb_cd_get_local_cols(desc_a),ione,&
         & alpha,x,size(x),beta,&
         & y,size(y),info)
  end if

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zaxpbyv
