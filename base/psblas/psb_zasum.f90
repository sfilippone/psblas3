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
! File: psb_zasum.f90
!
! Function: psb_zasum 
!    Computes norm1 of X
!
!    norm1 := sum(sub( X )(i))
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Parameters:
!    x      -  real,dimension(:,:).       The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset.
!
function psb_zasum (x,desc_a, info, jx)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(kind(1.d0)), intent(in)      :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  integer, optional, intent(in)     :: jx
  real(kind(1.d0))                  :: psb_zasum

  ! locals
  integer                  :: int_err(5), ictxt, np, npcol, me, mycol,&
       & err_act, n, iix, jjx, temp(2), ix, ijx, m, i
  real(kind(1.d0))         :: asum, dzasum
  character(len=20)        :: name, ch_err
  complex(kind(1.d0))      :: cmax
  double complex    ::         zdum
  double precision  ::  cabs1
  cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )

  name='psb_zasum'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  asum=0.d0

  ictxt=psb_cd_get_context(desc_a)

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

  ! check vector correctness
  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iix.ne.1) then
    info=3040
    call psb_errpush(info,name)
    goto 9999
  end if

  ! compute local max
  if ((m.ne.0)) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      asum=dzasum(psb_cd_get_local_rows(desc_a)-iix+1,x(iix,jjx),ione)

      ! adjust asum because overlapped elements are computed more than once
      i=1
      do while (desc_a%ovrlap_elem(i).ne.-ione)
        cmax = x(desc_a%ovrlap_elem(i)-iix+1,jjx)
        asum = asum -&
             & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
             &  cabs1(cmax)
        i = i+2
      end do

      ! compute global sum
      call psb_sum(ictxt, asum)

    else
      asum=0.d0
      ! compute global sum
      call psb_sum(ictxt, asum)
    end if
  else
     asum=0.d0
  end if
  

  psb_zasum=asum

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end function psb_zasum


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
! Function: psb_zasumv 
!    Computes norm1 of X
!
!    norm1 := sum(X(i))
!
! Parameters:
!    x      -  real,dimension(:).       The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
function psb_zasumv (x,desc_a, info)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(kind(1.d0)), intent(in)      :: x(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  real(kind(1.d0))                  :: psb_zasumv

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, n, iix, jjx, jx, ix, ijx, m, i
  real(kind(1.d0))         :: asum, dzasum
  character(len=20)        :: name, ch_err
  complex(kind(1.d0))      :: cmax
  double complex    ::         zdum
  double precision  ::  cabs1
  cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )

  name='psb_zasumv'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  asum=0.d0

  ictxt=psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  jx=1

  m = psb_cd_get_global_rows(desc_a)

  ! check vector correctness
  call psb_chkvect(m,1,size(x),ix,jx,desc_a,info,iix,jjx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iix.ne.1) then
    info=3040
    call psb_errpush(info,name)
    goto 9999
  end if

  ! compute local max
  if ((m.ne.0)) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      asum=dzasum(psb_cd_get_local_rows(desc_a),x,ione)

      ! adjust asum because overlapped elements are computed more than once
      i=1
      do while (desc_a%ovrlap_elem(i).ne.-ione)
        cmax = x(desc_a%ovrlap_elem(i))
        asum = asum -&
             & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
             & cabs1(cmax)
        i = i+2
      end do

      ! compute global sum
      call psb_sum(ictxt, asum)

    else
      asum=0.d0
      ! compute global sum
      call psb_sum(ictxt, asum)
    end if
  else
    asum=0.d0
  end if

  psb_zasumv=asum

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end function psb_zasumv


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
! Subroutine: psb_zasum vs
!    Computes norm1 of X
!
!    norm1 := sum(X(i))
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:).         The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset.
!
subroutine psb_zasumvs (res,x,desc_a, info)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(kind(1.d0)), intent(in)      :: x(:)
  real(kind(1.d0)), intent(out)     :: res
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, n, iix, jjx, ix, jx, ijx, m, i
  real(kind(1.d0))         :: asum, dzasum
  character(len=20)        :: name, ch_err
  complex(kind(1.d0))      :: cmax
  double complex    ::         zdum
  double precision  ::  cabs1
  cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )

  name='psb_zasumvs'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  asum=0.d0

  ictxt=psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  jx = 1

  m = psb_cd_get_global_rows(desc_a)
  ! check vector correctness
  call psb_chkvect(m,1,size(x),ix,jx,desc_a,info,iix,jjx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iix.ne.1) then
    info=3040
    call psb_errpush(info,name)
    goto 9999
  end if

  ! compute local max
  if ((m.ne.0)) then
    if(psb_cd_get_local_rows(desc_a).gt.0) then
      asum=dzasum(psb_cd_get_local_rows(desc_a),x,ione)

      ! adjust asum because overlapped elements are computed more than once
      i=1
      do while (desc_a%ovrlap_elem(i).ne.-ione)
        cmax = x(desc_a%ovrlap_elem(i))
        asum = asum -&
             & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
             & cabs1(cmax)
        i = i+2
      end do

      ! compute global sum
      call psb_sum(ictxt,asum)

    else
      asum=0.d0
      ! compute global sum
      call psb_sum(ictxt, asum)
    end if
  else
    asum=0.d0
  end if


  res = asum

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zasumvs
