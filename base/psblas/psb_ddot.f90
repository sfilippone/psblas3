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
! File: psb_ddot.f90
!
! Function: psb_ddot
!    psb_ddot forms the dot product of two distributed vectors,
!
!    dot := sub( X )**T * sub( Y )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Parameters:
!    x      -  real,dimension(:,:).       The input vector containing the entries of sub( X ).
!    y      -  real,dimension(:,:).       The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset for sub( X ).
!    jy     -  integer(optional).         The column offset for sub( Y ).
!
function psb_ddot(x, y,desc_a, info, jx, jy)  
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:,:), y(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(in), optional    :: jx, jy
  integer, intent(out)             :: info
  real(kind(1.D0))                 :: psb_ddot

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, n, iix, jjx, ix, ijx, iy, ijy, iiy, jjy, i, m, j, k
  real(kind(1.D0))         :: dot_local
  real(kind(1.d0))         :: ddot
  character(len=20)        :: name, ch_err

  name='psb_ddot'
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

  if(m.ne.0) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      dot_local = ddot(psb_cd_get_local_rows(desc_a),&
           & x(iix,jjx),ione,y(iiy,jjy),ione)
      ! adjust dot_local because overlapped elements are computed more than once
      i=1
      do while (desc_a%ovrlap_elem(i).ne.-ione)
        dot_local = dot_local -&
             & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
             & x(iix+desc_a%ovrlap_elem(i)-1,jjx)*&
             & y(iiy+desc_a%ovrlap_elem(i)-1,jjy)
        i = i+2
      end do
    else
      dot_local=0.d0
    end if
  else
    dot_local=0.d0
  end if

  ! compute global sum
  call psb_sum(ictxt, dot_local)

  psb_ddot = dot_local

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end function psb_ddot




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
! Function: psb_ddotv
!    psb_ddot forms the dot product of two distributed vectors,
!
!    dot := X**T * Y
!
! Parameters:
!    x      -  real,dimension(:).         The input vector containing the entries of X.
!    y      -  real,dimension(:).         The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
function psb_ddotv(x, y,desc_a, info)  
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:), y(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  real(kind(1.D0))                 :: psb_ddotv

  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, n, iix, jjx, ix, jx, iy, jy, iiy, jjy, i, m, j, k
  real(kind(1.D0))         :: dot_local
  real(kind(1.d0))         :: ddot
  character(len=20)        :: name, ch_err

  name='psb_ddot'
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
  jx = ione
  jy = ione
  m = psb_cd_get_global_rows(desc_a)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,jx,desc_a,info,iix,jjx)
  if (info == 0) &
       & call psb_chkvect(m,ione,size(y,1),iy,jy,desc_a,info,iiy,jjy)
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

  if(m.ne.0) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      dot_local = ddot(psb_cd_get_local_rows(desc_a),&
           & x,ione,y,ione)
      ! adjust dot_local because overlapped elements are computed more than once
      i=1
      do while (desc_a%ovrlap_elem(i).ne.-ione)
        dot_local = dot_local -&
             & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
             & x(desc_a%ovrlap_elem(i))*&
             & y(desc_a%ovrlap_elem(i))
        i = i+2
      end do
    else
      dot_local=0.d0
    end if
  else
    dot_local=0.d0
  end if

  ! compute global sum
  call psb_sum(ictxt, dot_local)

  psb_ddotv = dot_local

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end function psb_ddotv



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
! Subroutine: psb_ddotvs
!    psb_ddot forms the dot product of two distributed vectors,
!
!    dot := X**T * Y
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:).         The input vector containing the entries of X.
!    y      -  real,dimension(:).         The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
subroutine psb_ddotvs(res, x, y,desc_a, info)  
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:), y(:)
  real(kind(1.d0)), intent(out)    :: res
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, n, iix, jjx, ix, ijx, iy, ijy, iiy, jjy, i, m, j, k
  real(kind(1.D0))         :: dot_local
  real(kind(1.d0))         :: ddot
  character(len=20)        :: name, ch_err

  name='psb_ddot'
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
  call psb_chkvect(m,ione,size(x,1),ix,ix,desc_a,info,iix,jjx)
  if (info == 0) &
       & call psb_chkvect(m,ione,size(y,1),iy,iy,desc_a,info,iiy,jjy)
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

  if(m.ne.0) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      dot_local = ddot(psb_cd_get_local_rows(desc_a),&
           & x,ione,y,ione)
      ! adjust dot_local because overlapped elements are computed more than once
      i=1
      do while (desc_a%ovrlap_elem(i).ne.-ione)
        dot_local = dot_local -&
             & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
             & x(desc_a%ovrlap_elem(i))*&
             & y(desc_a%ovrlap_elem(i))
        i = i+2
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

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_ddotvs




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
! Subroutine: psb_dmdots
!    psb_ddot forms the dot product of two distributed vectors,
!
!    dot := sub( X )**T * sub( Y )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:,:).       The input vector containing the entries of sub( X ).
!    y      -  real,dimension(:,:).       The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
subroutine psb_dmdots(res, x, y, desc_a, info)  
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:,:), y(:,:)
  real(kind(1.d0)), intent(out)    :: res(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, n, iix, jjx, ix, ijx, iy, ijy, iiy, jjy, i, m, j, k
  real(kind(1.d0)),allocatable  :: dot_local(:)
  real(kind(1.d0))         :: ddot
  character(len=20)        :: name, ch_err

  name='psb_dmdots'
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
  call psb_chkvect(m,ione,size(x,1),ix,ix,desc_a,info,iix,jjx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,ione,size(y,1),iy,iy,desc_a,info,iiy,jjy)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((ix.ne.ione).or.(iy.ne.ione)) then
    info=3040
    call psb_errpush(info,name)
    goto 9999
  end if

  k = min(size(x,2),size(y,2))
  allocate(dot_local(k))

  if(m.ne.0) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      do j=1,k
        dot_local(j) = ddot(psb_cd_get_local_rows(desc_a),&
             & x(1,j),ione,y(1,j),ione)
        ! adjust dot_local because overlapped elements are computed more than once
        i=1
        do while (desc_a%ovrlap_elem(i).ne.-ione)
          dot_local(j) = dot_local(j) -&
               & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
               & x(desc_a%ovrlap_elem(i)-1,j)*&
               & y(desc_a%ovrlap_elem(i)-1,j)
          i = i+2
        end do
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

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dmdots


subroutine psb_ddot2v(res, x, y,w,z,desc_a, info)  
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:), y(:),w(:), z(:)
  real(kind(1.d0)), intent(out)    :: res(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, n, iix, jjx, ix, ijx, iy, ijy, iiy, jjy, i, m, j, k
  real(kind(1.D0))         :: dot_local(2)
  real(kind(1.d0))         :: ddot
  character(len=20)        :: name, ch_err

  name='psb_ddot'
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
  call psb_chkvect(m,ione,size(x,1),ix,ix,desc_a,info,iix,jjx)
  if (info == 0) &
       & call psb_chkvect(m,ione,size(y,1),iy,iy,desc_a,info,iiy,jjy)
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

  if(m.ne.0) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      dot_local(1) = ddot(psb_cd_get_local_rows(desc_a),&
           & x,ione,y,ione)
      dot_local(2) = ddot(psb_cd_get_local_rows(desc_a),&
           & w,ione,z,ione)
      ! adjust dot_local because overlapped elements are computed more than once
      i=1
      do while (desc_a%ovrlap_elem(i).ne.-ione)
        dot_local(1) = dot_local(1) -&
             & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
             & x(desc_a%ovrlap_elem(i))*&
             & y(desc_a%ovrlap_elem(i))
        dot_local(2) = dot_local(2) -&
             & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
             & w(desc_a%ovrlap_elem(i))*&
             & z(desc_a%ovrlap_elem(i))
        i = i+2
      end do
    else
      dot_local=0.d0
    end if
  else
    dot_local=0.d0
  end if

  ! compute global sum
  call psb_sum(ictxt, dot_local)

  res(1:2) = dot_local(1:2)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_ddot2v

