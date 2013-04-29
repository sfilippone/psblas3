!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
! File: psb_zaxpby.f90

subroutine psb_zaxpby_vect(alpha, x, beta, y,&
     & desc_a, info)
  use psb_base_mod, psb_protect_name => psb_zaxpby_vect
  implicit none                    
  type(psb_z_vect_type), intent (inout) ::  x
  type(psb_z_vect_type), intent (inout) ::  y
  complex(psb_dpk_), intent (in)        :: alpha, beta
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)                  :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ix, iy, m, iiy, jjy
  character(len=20)        :: name, ch_err

  name='psb_zgeaxpby'
  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(y%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m,ione,x%get_nrows(),ix,ione,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 1'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,ione,y%get_nrows(),iy,ione,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 2'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((iix /= ione).or.(iiy /= ione)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
  end if

  if(desc_a%get_local_rows() > 0) then
    call y%axpby(desc_a%get_local_rows(),&
         & alpha,x,beta,info)
  end if

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_zaxpby_vect

!
! Subroutine: psb_zaxpby
!    Adds one distributed matrix to another,
!
!    sub( Y ) := beta * sub( Y ) + alpha * sub( X )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Arguments:
!    alpha  -  complex,input        The scalar used to multiply each component of X
!    x(:,:) -  complex,input        The input vector containing the entries of X
!    beta   -  complex,input        The scalar used to multiply each component of Y
!    y(:,:) -  complex,inout        The input vector Y
!    desc_a -  type(psb_desc_type)  The communication descriptor.
!    info   -  integer              Return code
!    jx     -  integer(optional)    The column offset for X 
!    jy     -  integer(optional)    The column offset for Y 
!
subroutine  psb_zaxpby(alpha, x, beta,y,desc_a,info, n, jx, jy)
  use psb_base_mod, psb_protect_name => psb_zaxpby

  implicit none                    

  integer(psb_ipk_), intent(in), optional   :: n, jx, jy
  integer(psb_ipk_), intent(out)            :: info
  type(psb_desc_type), intent(in) :: desc_a
  complex(psb_dpk_), intent(in)    :: alpha, beta
  complex(psb_dpk_), intent(in)    :: x(:,:)
  complex(psb_dpk_), intent(inout) :: y(:,:)

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ix, iy, ijx, ijy, m, iiy, in, jjy, &
       & lldx, lldy
  character(len=20)        :: name, ch_err

  name='psb_geaxpby'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -ione) then
    info = psb_err_context_error_
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
    if(((ijx+n) <= size(x,2)).and.&
         & ((ijy+n) <= size(y,2))) then 
      in = n
    else
      in = min(size(x,2),size(y,2))
    end if
  else
    in = min(size(x,2),size(y,2))
  endif

  if(ijx /= ijy) then
    info=3050
    call psb_errpush(info,name)
    goto 9999
  end if

  m = desc_a%get_global_rows()
  lldx = size(x,1)
  lldy = size(y,1)
  ! check vector correctness
  call psb_chkvect(m,ione,lldx,ix,ijx,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,ione,lldy,iy,ijy,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((iix /= ione).or.(iiy /= ione)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if ((in /= 0)) then
    if(desc_a%get_local_rows() > 0) then
      call caxpby(desc_a%get_local_cols(),in,&
           & alpha,x(iix:,jjx),lldx,beta,&
           & y(iiy:,jjy),lldy,info)
    end if
  end if

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zaxpby





!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
! Subroutine: psb_zaxpbyv
!    Adds one distributed matrix to another,
!
!    Y := beta * Y + alpha * X
!
! Arguments:
!    alpha  -  complex,input        The scalar used to multiply each component of X
!    x(:)   -  complex,input        The input vector containing the entries of X
!    beta   -  complex,input        The scalar used to multiply each component of Y
!    y(:)   -  complex,inout        The input vector Y
!    desc_a -  type(psb_desc_type)  The communication descriptor.
!    info   -  integer              Return code
!
!
subroutine  psb_zaxpbyv(alpha, x, beta,y,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_zaxpbyv
  implicit none                    

  integer(psb_ipk_), intent(out)            :: info
  type(psb_desc_type), intent(in) :: desc_a
  complex(psb_dpk_), intent(in)    :: alpha, beta
  complex(psb_dpk_), intent(in)    :: x(:)
  complex(psb_dpk_), intent(inout) :: y(:)

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ix, iy, m, iiy, jjy, &
       & lldx, lldy
  character(len=20)        :: name, ch_err
  logical, parameter :: debug=.false.

  name='psb_geaxpby'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()
  lldx = size(x,1)
  lldy = size(y,1)
  ! check vector correctness
  call psb_chkvect(m,ione,lldx,ix,ione,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 1'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,ione,lldy,iy,ione,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 2'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((iix /= ione).or.(iiy /= ione)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
  end if

  if(desc_a%get_local_rows() > 0) then
    call caxpby(desc_a%get_local_cols(),ione,&
         & alpha,x,lldx,beta,&
         & y,lldy,info)
  end if

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zaxpbyv
