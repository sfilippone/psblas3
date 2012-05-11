!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
! File: psb_sasum.f90
!
! Function: psb_sasum 
!    Computes norm1 of X
!
!    norm1 := sum(sub( X )(i))
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Arguments:
!    x(:,:) -  real                The input vector.
!    desc_a -  type(psb_desc_type).   The communication descriptor.
!    info   -  integer.               Return code
!    jx     -  integer(optional).     The column offset.
!
function psb_sasum (x,desc_a, info, jx) result(res)
  use psb_base_mod, psb_protect_name => psb_sasum

  implicit none

  real(psb_spk_), intent(in)   :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_), optional, intent(in)     :: jx
  real(psb_spk_)                  :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, &
       & err_act, iix, jjx, ix, ijx, m, i, idx, ndm, ldx
  character(len=20)        :: name, ch_err

  name='psb_sasum'
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
  if (present(jx)) then
    ijx = jx
  else
    ijx = 1
  endif

  m = desc_a%get_global_rows()
  ldx = size(x,1)
  ! check vector correctness
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
  if(desc_a%get_local_rows() > 0) then
    res = psb_asum(desc_a%get_local_rows()-iix+1,x(:,jjx))

    ! adjust res because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx  = desc_a%ovrlap_elem(i,1)
      ndm  = desc_a%ovrlap_elem(i,2)
      res = res - (real(ndm-1)/real(ndm))*psb_nrm1(x(idx,jjx))
    end do

  else
    res = szero
  end if
  ! compute global sum
  call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end function psb_sasum


function psb_sasum_vect(x, desc_a, info) result(res)
  use psb_base_mod, psb_protect_name => psb_sasum_vect
  implicit none

  real(psb_spk_)                        :: res
  type(psb_s_vect_type), intent (inout) :: x
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, jx, ix, m, imax
  character(len=20)        :: name, ch_err

  name='psb_sasumv'
  if (psb_errstatus_fatal()) return 
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
    res = x%asum(desc_a%get_local_rows())
  else 
    res = szero
  end if

  ! compute global sum
  call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end function psb_sasum_vect



!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
! Function: psb_sasumv 
!    Computes norm1 of X
!
!    norm1 := sum(X(i))
!
! Arguments:
!    x(:)   -  real               The input vector.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!
function psb_sasumv(x,desc_a, info) result(res)
  use psb_base_mod, psb_protect_name => psb_sasumv

  implicit none

  real(psb_spk_), intent(in)   :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info
  real(psb_spk_)                  :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, jx, ix, m, i, idx, ndm, ldx
  character(len=20)        :: name, ch_err

  name='psb_sasumv'
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
  jx=1

  m = desc_a%get_global_rows()
  ldx = size(x,1)
  ! check vector correctness
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
  if(desc_a%get_local_rows() > 0) then
    res = psb_asum(desc_a%get_local_rows(),x)

    ! adjust asum because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      res = res - (real(ndm-1)/real(ndm))*psb_nrm1(x(idx))
    end do

  else
    res = szero
  end if

  ! compute global sum
  call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end function psb_sasumv


!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
! Subroutine: psb_sasumvs
!    Computes norm1 of X
!
!    norm1 := sum(X(i))
!
! Arguments:
!    res    -  real.                 The result.
!    x(:)   -  real               The input vector.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    jx     -  integer(optional).    The column offset.
!
subroutine psb_sasumvs(res,x,desc_a, info)
  use psb_base_mod, psb_protect_name => psb_sasumvs

  implicit none

  real(psb_spk_), intent(in)    :: x(:)
  real(psb_spk_), intent(out)      :: res
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)   :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ix, jx, m, i, idx, ndm, ldx
  character(len=20)        :: name, ch_err

  name='psb_sasumvs'
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
  ! check vector correctness
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
  if(desc_a%get_local_rows() > 0) then
    res = psb_asum(desc_a%get_local_rows(),x)

    ! adjust asum because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      res = res - (real(ndm-1)/real(ndm))*psb_nrm1(x(idx))
    end do

  else
    res = szero
  end if

  ! compute global sum
  call psb_sum(ictxt,res)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_sasumvs
