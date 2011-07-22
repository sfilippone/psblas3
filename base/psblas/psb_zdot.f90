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
! File: psb_zdot.f90
!
! Function: psb_zdot
!    psb_zdot forms the dot product of two distributed vectors,
!
!    dot := sub( X )**C * sub( Y )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Arguments:
!    x(:,:) -  complex               The input vector containing the entries of sub( X ).
!    y(:,:) -  complex               The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    jx     -  integer(optional).    The column offset for sub( X ).
!    jy     -  integer(optional).    The column offset for sub( Y ).
!
function psb_zdot(x, y,desc_a, info, jx, jy)  
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(psb_dpk_), intent(in)     :: x(:,:), y(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(in), optional    :: jx, jy
  integer, intent(out)             :: info
  complex(psb_dpk_)              :: psb_zdot

  ! locals
  integer                  :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, ijx, iy, ijy, iiy, jjy, i, m
  complex(psb_dpk_)         :: dot_local
  complex(psb_dpk_)         :: zdotc
  character(len=20)        :: name, ch_err

  name='psb_zdot'
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

  if(ijx /= ijy) then
    info=3050
    call psb_errpush(info,name)
    goto 9999
  end if

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,ione,size(y,1),iy,ijy,desc_a,info,iiy,jjy)
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

  if(m /= 0) then
    if(desc_a%get_local_rows() > 0) then
      dot_local = zdotc(desc_a%get_local_rows(),&
           & x(iix:,jjx),ione,y(iiy:,jjy),ione)
      ! adjust dot_local because overlapped elements are computed more than once
      do i=1,size(desc_a%ovrlap_elem,1)
        idx  = desc_a%ovrlap_elem(i,1)
        ndm  = desc_a%ovrlap_elem(i,2)
        dot_local = dot_local - (real(ndm-1)/real(ndm))*(conjg(x(idx,jjx))*y(idx,jjy))
      end do
    else
      dot_local=0.d0
    end if
  else
    dot_local=0.d0
  end if

  ! compute global sum
  call psb_sum(ictxt, dot_local)

  psb_zdot = dot_local

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end function psb_zdot




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
! Function: psb_zdotv
!    psb_zdotv forms the dot product of two distributed vectors,
!
!    dot := X**C * Y
!
! Arguments:
!    x(:)   -  complex               The input vector containing the entries of X.
!    y(:)   -  complex               The input vector containing the entries of Y.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!
function psb_zdotv(x, y,desc_a, info)  
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(psb_dpk_), intent(in)  :: x(:), y(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  complex(psb_dpk_)              :: psb_zdotv

  ! locals
  integer                  :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, jx, iy, jy, iiy, jjy, i, m
  complex(psb_dpk_)         :: dot_local
  complex(psb_dpk_)         :: zdotc
  character(len=20)        :: name, ch_err

  name='psb_zdot'
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
  jx = ione
  jy = ione
  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,jx,desc_a,info,iix,jjx)
  if (info == psb_success_)&
       & call psb_chkvect(m,ione,size(y,1),iy,jy,desc_a,info,iiy,jjy)
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

  if(m /= 0) then
    if(desc_a%get_local_rows() > 0) then
      dot_local = zdotc(desc_a%get_local_rows(),&
           & x,ione,y,ione)
      ! adjust dot_local because overlapped elements are computed more than once
      do i=1,size(desc_a%ovrlap_elem,1)
        idx  = desc_a%ovrlap_elem(i,1)
        ndm  = desc_a%ovrlap_elem(i,2)
        dot_local = dot_local - (real(ndm-1)/real(ndm))*(conjg(x(idx))*y(idx))
      end do
    else
      dot_local=0.d0
    end if
  else
    dot_local=0.d0
  end if

  ! compute global sum
  call psb_sum(ictxt, dot_local)

  psb_zdotv = dot_local

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end function psb_zdotv



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
! Subroutine: psb_zdotvs
!    psb_zdotvs forms the dot product of two distributed vectors,
!
!    res := X**C * Y
!
! Arguments:
!    res    -  complex.             The result.
!    x(:)   -  complex              The input vector containing the entries of X.
!    y(:)   -  complex              The input vector containing the entries of Y.
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!
subroutine psb_zdotvs(res, x, y,desc_a, info)  
  use psb_base_mod, psb_protect_name => psb_zdotvs
  implicit none

  complex(psb_dpk_), intent(in)     :: x(:), y(:)
  complex(psb_dpk_), intent(out)    :: res
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer                  :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, iy, iiy, jjy, i, m
  complex(psb_dpk_)         :: dot_local
  complex(psb_dpk_)         :: zdotc
  character(len=20)        :: name, ch_err

  name='psb_zdot'
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
  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ix,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,ione,size(y,1),iy,iy,desc_a,info,iiy,jjy)
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

  if(m /= 0) then
    if(desc_a%get_local_rows() > 0) then
      dot_local = zdotc(desc_a%get_local_rows(),&
           & x,ione,y,ione)
      ! adjust dot_local because overlapped elements are computed more than once
      do i=1,size(desc_a%ovrlap_elem,1)
        idx  = desc_a%ovrlap_elem(i,1)
        ndm  = desc_a%ovrlap_elem(i,2)
        dot_local = dot_local - (real(ndm-1)/real(ndm))*(conjg(x(idx))*y(idx))
      end do
    else
      dot_local=0.d0
    end if
  else
    dot_local=0.d0
  end if

  ! compute global sum
  call psb_sum(ictxt, dot_local)

  res = dot_local

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zdotvs




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
! Subroutine: psb_zmdots
!    psb_zmdots forms the dot product of multiple distributed vectors,
!
!    res(i) := ( X(:,i) )**C * ( Y(:,i) )
!
! Arguments:
!    res(:) -  complex.             The result.
!    x(:)   -  complex              The input vector containing the entries of sub( X ).
!    y(:)   -  complex              The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!
subroutine psb_zmdots(res, x, y, desc_a, info)  
  use psb_base_mod, psb_protect_name => psb_zmdots
  implicit none

  complex(psb_dpk_), intent(in)     :: x(:,:), y(:,:)
  complex(psb_dpk_), intent(out)    :: res(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer                  :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, iy, iiy, jjy, i, m, j, k
  complex(psb_dpk_),allocatable  :: dot_local(:)
  complex(psb_dpk_)         :: zdotc
  character(len=20)        :: name, ch_err

  name='psb_zmdots'
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

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ix,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,ione,size(y,1),iy,iy,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((ix /= ione).or.(iy /= ione)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  k = min(size(x,2),size(y,2))
  allocate(dot_local(k))

  if(m /= 0) then
    if(desc_a%get_local_rows() > 0) then
      do j=1,k
        dot_local(j) = zdotc(desc_a%get_local_rows(),&
             & x(1:,j),ione,y(1:,j),ione)
        ! adjust dot_local because overlapped elements are computed more than once
      end do
      do i=1,size(desc_a%ovrlap_elem,1)
        idx  = desc_a%ovrlap_elem(i,1)
        ndm  = desc_a%ovrlap_elem(i,2)
        dot_local(1:k) = dot_local(1:k) - (real(ndm-1)/real(ndm))*(conjg(x(idx,1:k))*y(idx,1:k))
      end do
    else
      dot_local(:)=0.d0
    end if
  else
    dot_local(:)=0.d0
  end if

  ! compute global sum
  call psb_sum(ictxt, dot_local(1:k))

  res(1:k) = dot_local(1:k)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zmdots
