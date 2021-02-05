!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
! File: psb_cnrm2.f90
!
! Function: psb_cnrm2
!    Computes the norm2 of a distributed vector,
!
!    norm2 := sqrt ( sub( X )**C * sub( X ) )
!
!    where sub( X ) denotes X(:,JX).
!
! Arguments:
!    x(:,:) -  complex              The input vector containing the entries of sub( X ).
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!    jx     -  integer(optional).   The column offset for sub( X ).
!    global -  logical(optional)    Whether to perform the global reduction, default: .true.
!
function psb_cnrm2(x, desc_a, info, jx,global)  result(res)
  use psb_desc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(psb_spk_), intent(in)      ::  x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer(psb_ipk_), intent(in), optional     :: jx
  integer(psb_ipk_), intent(out)              :: info
  real(psb_spk_)                  :: res
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, ndim, i, id, idx, ndm, ldx
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
  logical :: global_
  real(psb_spk_)         :: scnrm2, dd
  character(len=20)      :: name, ch_err

  name='psb_cnrm2'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = 1
  if (present(jx)) then
    ijx = jx
  else
    ijx = 1
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  m = desc_a%get_global_rows()
  ldx = size(x,1)
  call psb_chkvect(m,lone,ldx,ix,ijx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix /= 1) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (desc_a%get_local_rows() > 0) then
    ndim = desc_a%get_local_rows()
    res  = scnrm2( int(ndim,kind=psb_mpk_), x(iix:,jjx), int(ione,kind=psb_mpk_) )

    ! adjust  because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      dd  = real(ndm-1)/real(ndm)
      res = res * sqrt(sone - dd*(abs(x(idx,jjx))/res)**2)
    end do
  else
    res = szero
  end if

  if (global_) call psb_nrm2(ctxt,res)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end function psb_cnrm2



!!$
!!$              Parallel Sparse BLAS  version 3.5
!!$    (C) Copyright 2006-2018
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari
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
! Function: psb_cnrm2
!    Computes the norm2 of a distributed vector,
!
!    norm2 := sqrt ( X**C * X)
!
! Arguments:
!    x(:)   -  complex               The input vector containing the entries of X.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    global -  logical(optional)    Whether to perform the global reduction, default: .true.
!
function psb_cnrm2v(x, desc_a, info,global)  result(res)
  use psb_desc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(psb_spk_), intent(in)    :: x(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)   :: info
  real(psb_spk_)                   :: res
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, ndim, i, id, idx, ndm, ldx
  integer(psb_lpk_) :: ix, jx, iy, ijy, m
  real(psb_spk_)         :: scnrm2, dd
  logical :: global_
  character(len=20)        :: name, ch_err

  name='psb_cnrm2v'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if
  ix = 1
  jx=1
  m = desc_a%get_global_rows()
  ldx = size(x,1)
  call psb_chkvect(m,lone,ldx,ix,jx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix /= 1) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (desc_a%get_local_rows() > 0) then
    ndim = desc_a%get_local_rows()
    res  = scnrm2( int(ndim,kind=psb_mpk_), x, int(ione,kind=psb_mpk_) )
    ! adjust  because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      dd  = real(ndm-1)/real(ndm)
      res = res * sqrt(sone - dd*(abs(x(idx))/res)**2)
    end do
  else
    res = szero
  end if

  if (global_) call psb_nrm2(ctxt,res)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end function psb_cnrm2v


! Function: psb_cnrm2_vect
!    Computes the norm2 of a distributed vector,
!
!    norm2 := sqrt ( X**C * X)
!
! Arguments:
!    x      -  type(psb_c_vect_type) The input vector containing the entries of X.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    global -  logical(optional)    Whether to perform the global reduction, default: .true.
!
function psb_cnrm2_vect(x, desc_a, info,global)  result(res)
  use psb_desc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_c_vect_mod
  implicit none

  real(psb_spk_)                        :: res
  type(psb_c_vect_type), intent (inout) :: x
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, ndim, i, id, idx, ndm, ldx
  integer(psb_lpk_) :: ix, jx, iy, ijy, m
  logical :: global_
  real(psb_spk_)         :: snrm2, dd
  character(len=20)      :: name, ch_err

  name='psb_cnrm2v'
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  info=psb_success_

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  ix = 1
  jx = 1
  m  = desc_a%get_global_rows()
  ldx = x%get_nrows()
  call psb_chkvect(m,lone,ldx,ix,jx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix /= 1) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (desc_a%get_local_rows() > 0) then
    ndim = desc_a%get_local_rows()
    res  = x%nrm2(ndim)
    ! adjust  because overlapped elements are computed more than once
    if (size(desc_a%ovrlap_elem,1)>0) then
      if (x%is_dev()) call x%sync()
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        dd  = dble(ndm-1)/dble(ndm)
        res = res - sqrt(cone - dd*(abs(x%v%v(idx))/res)**2)
      end do
    end if
  else
    res = szero
  end if

  if (global_) call psb_nrm2(ctxt,res)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end function psb_cnrm2_vect

! Function: psb_cnrm2_weight_vect
!    Computes the weighted norm2 of a distributed vector,
!
!    norm2 := sqrt ( (w.*X)**C * (w.*X))
!
! Arguments:
!    x      -  type(psb_c_vect_type) The input vector containing the entries of X.
!    w      -  type(psb_c_vect_type) The input vector containing the entries of W.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    global -  logical(optional)    Whether to perform the global reduction, default: .true.
!
function psb_cnrm2_weight_vect(x,w, desc_a, info,global,aux)  result(res)
  use psb_desc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_c_vect_mod
  implicit none

  real(psb_spk_)                        :: res
  type(psb_c_vect_type), intent (inout) :: x
  type(psb_c_vect_type), intent (inout) :: w
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  logical, intent(in), optional        :: global
  type(psb_c_vect_type), intent(inout), optional :: aux

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, ndim, i, id, idx, ndm, ldx
  integer(psb_lpk_) :: ix, jx, iy, ijy, m
  logical :: global_
  real(psb_spk_)         :: snrm2, dd
  character(len=20)      :: name, ch_err

  name='psb_cnrm2v_weight'
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  info=psb_success_

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  ix = 1
  jx = 1
  m  = desc_a%get_global_rows()
  ldx = x%get_nrows()
  call psb_chkvect(m,lone,ldx,ix,jx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix /= 1) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (desc_a%get_local_rows() > 0) then
    ndim = desc_a%get_local_rows()
    res  = x%nrm2(ndim,w,aux)
    ! adjust  because overlapped elements are computed more than once
    if (size(desc_a%ovrlap_elem,1)>0) then
      if (x%is_dev()) call x%sync()
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        dd  = dble(ndm-1)/dble(ndm)
        res = res - sqrt(cone - dd*(abs(x%v%v(idx))/res)**2)
      end do
    end if
  else
    res = szero
  end if

  if (global_) call psb_nrm2(ctxt,res)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end function psb_cnrm2_weight_vect

! Function: psb_cnrm2_weight_vect
!    Computes the weighted norm2 of a distributed vector with respect to a mask
!    contained in the vector id.
!
!    norm2 := sqrt ( (w(id > 0).*X(id > 0))**C * (w(id > 0).*X(id > 0)))
!
! Arguments:
!    x      -  type(psb_c_vect_type) The input vector containing the entries of X.
!    w      -  type(psb_c_vect_type) The input vector containing the entries of W.
!    id     -  type(psb_c_vect_type) The inpute vector containing the mask
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    global -  logical(optional)    Whether to perform the global reduction, default: .true.
!
function psb_cnrm2_weightmask_vect(x,w,idv, desc_a, info,global, aux)  result(res)
  use psb_desc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_c_vect_mod
  implicit none

  real(psb_spk_)                        :: res
  type(psb_c_vect_type), intent (inout) :: x
  type(psb_c_vect_type), intent (inout) :: w
  type(psb_c_vect_type), intent (inout) :: idv
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  logical, intent(in), optional        :: global
  type(psb_c_vect_type), intent(inout), optional :: aux

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, ndim, i, id, idx, ndm, ldx
  integer(psb_lpk_) :: ix, jx, iy, ijy, m
  logical :: global_
  real(psb_spk_)         :: snrm2, dd
  character(len=20)      :: name, ch_err

  name='psb_cnrm2v_weightmask'
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  info=psb_success_

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  ix = 1
  jx = 1
  m  = desc_a%get_global_rows()
  ldx = x%get_nrows()
  call psb_chkvect(m,lone,ldx,ix,jx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix /= 1) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (desc_a%get_local_rows() > 0) then
    ndim = desc_a%get_local_rows()
    res  = x%nrm2(ndim,w,idv,info,aux)
    ! adjust  because overlapped elements are computed more than once
    if (size(desc_a%ovrlap_elem,1)>0) then
      if (x%is_dev()) call x%sync()
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        dd  = dble(ndm-1)/dble(ndm)
        res = res - sqrt(cone - dd*(abs(x%v%v(idx))/res)**2)
      end do
    end if
  else
    res = szero
  end if

  if (global_) call psb_nrm2(ctxt,res)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end function psb_cnrm2_weightmask_vect

!!$
!!$              Parallel Sparse BLAS  version 3.5
!!$    (C) Copyright 2006-2018
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari
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
! Subroutine: psb_cnrm2vs
!    Computes the norm2 of a distributed vector, subroutine version
!
!    norm2 := sqrt ( X**C * X)
!
! Arguments:
!    res    -  real                  The result.
!    x(:)   -  complex               The input vector containing the entries of X.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    global -  logical(optional)    Whether to perform the global reduction, default: .true.
!
subroutine psb_cnrm2vs(res, x, desc_a, info,global)
  use psb_desc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(psb_spk_), intent(in)    :: x(:)
  real(psb_spk_), intent(out)      :: res
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)   :: info
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, ndim, i, id, idx, ndm, ldx
  integer(psb_lpk_) :: ix, jx, iy, ijy, m
  logical :: global_
  real(psb_spk_)         :: nrm2, scnrm2, dd
  character(len=20)        :: name, ch_err

  name='psb_cnrm2'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info=psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  ix = 1
  jx = 1
  m = desc_a%get_global_rows()
  ldx = size(x,1)
  call psb_chkvect(m,lone,ldx,ix,jx,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix /= 1) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (desc_a%get_local_rows() > 0) then
    ndim = desc_a%get_local_rows()
    res  = scnrm2( int(ndim,kind=psb_mpk_), x, int(ione,kind=psb_mpk_) )

    ! adjust  because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      dd  = real(ndm-1)/real(ndm)
      res = res * sqrt(sone - dd*(abs(x(idx))/res)**2)
    end do
  else
    res = szero
  end if

  if (global_) call psb_nrm2(ctxt,res)


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psb_cnrm2vs
