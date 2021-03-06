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
! File: psb_zdot.f90
!
! Function: psb_zdot_vect
!    psb_zdot computes the dot product of two distributed vectors,
!
!    dot := ( X )**C * ( Y )
!
!
! Arguments:
!    x      -  type(psb_z_vect_type) The input vector containing the entries of sub( X ).
!    y      -  type(psb_z_vect_type) The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    global -  logical(optional)     Whether to perform the global sum, default: .true.
!
!  Note: from a functional point of view, X and Y are input, but here
!        they are declared INOUT because of the sync() methods. 
!
!
function psb_zdot_vect(x, y, desc_a,info,global) result(res)
  use psb_desc_mod
  use psb_z_base_mat_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_z_vect_mod
  use psb_z_psblas_mod, psb_protect_name => psb_zdot_vect
  implicit none 
  complex(psb_dpk_)                    :: res
  type(psb_z_vect_type), intent(inout) :: x, y
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)       :: info
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, idx, ndm,&
       & err_act, iix, jjx, iiy, jjy, i, nr
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
  logical :: global_
  character(len=20)      :: name, ch_err

  name='psb_zdot_vect'
  res = zzero
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()
  call psb_info(ctxt, me, np)
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

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  ix = ione
  ijx = ione

  iy = ione
  ijy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m,lone,x%get_nrows(),ix,ijx,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,lone,y%get_nrows(),iy,ijy,desc_a,info,iiy,jjy)
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

  nr = desc_a%get_local_rows() 
  if(nr > 0) then
    res = x%dot(nr,y)
    ! FIXME
    ! adjust dot_local because overlapped elements are computed more than once
    if (size(desc_a%ovrlap_elem,1)>0) then
      if (x%is_dev()) call x%sync()
      if (y%is_dev()) call y%sync()
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        res = res - (real(ndm-1)/real(ndm))*(x%v%v(idx)*y%v%v(idx))
      end do
    end if
  else
    res = zzero
  end if

  ! compute global sum
  if (global_) call psb_sum(ctxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

  return

end function psb_zdot_vect
!
! Function: psb_zdot
!    psb_zdot computes the dot product of two distributed vectors,
!
!    dot := sub( X )**C * sub( Y )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Arguments:
!    x(:,:) -  complex                The input vector containing the entries of sub( X ).
!    y(:,:) -  complex                The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    jx     -  integer(optional).    The column offset for sub( X ).
!    jy     -  integer(optional).    The column offset for sub( Y ).
!    global -  logical(optional)     Whether to perform the global sum, default: .true.
!
function psb_zdot(x, y,desc_a, info, jx, jy,global)  result(res)
  use psb_base_mod, psb_protect_name => psb_zdot
  implicit none

  complex(psb_dpk_), intent(in)    :: x(:,:), y(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(in), optional    :: jx, jy
  integer(psb_ipk_), intent(out)   :: info
  complex(psb_dpk_)              :: res
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, idx, ndm,&
       & err_act, iix, jjx, iiy, jjy, i, nr, lldx, lldy
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
  complex(psb_dpk_)        :: zdotc
  logical :: global_
  character(len=20)        :: name, ch_err

  name='psb_zdot'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()
  call psb_info(ctxt, me, np)
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

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  m = desc_a%get_global_rows()
  lldx = size(x,1)
  lldy = size(y,1)

  ! check vector correctness
  call psb_chkvect(m,lone,lldx,ix,ijx,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,lone,lldy,iy,ijy,desc_a,info,iiy,jjy)
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

  nr = desc_a%get_local_rows() 
  if(nr > 0) then
    res = zdotc(int(nr,kind=psb_mpk_), x(iix:,jjx),1,y(iiy:,jjy),1)
    ! adjust dot_local because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx  = desc_a%ovrlap_elem(i,1)
      ndm  = desc_a%ovrlap_elem(i,2)
      res = res - (real(ndm-1)/real(ndm))*(conjg(x(idx,jjx))*y(idx,jjy))
    end do
  else
    res = zzero
  end if

  ! compute global sum
  if (global_) call psb_sum(ctxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

  return
end function psb_zdot




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
! Function: psb_zdotv
!    psb_zdotv computes the dot product of two distributed vectors,
!
!    dot := X**C * Y
!
! Arguments:
!    x(:)   -  complex               The input vector containing the entries of X.
!    y(:)   -  complex               The input vector containing the entries of Y.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    global -  logical(optional)     Whether to perform the global sum, default: .true.
!
function psb_zdotv(x, y,desc_a, info,global)  result(res)
  use psb_base_mod, psb_protect_name => psb_zdotv
  implicit none

  complex(psb_dpk_), intent(in)   :: x(:), y(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info
  complex(psb_dpk_)              :: res
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, idx, ndm,&
       & err_act, iix, jjx, iiy, jjy, i, nr, lldx, lldy
  integer(psb_lpk_) :: ix, jx, iy, jy, m
  logical :: global_
  complex(psb_dpk_)         :: zdotc
  character(len=20)        :: name, ch_err

  name='psb_zdot'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  ix = ione
  iy = ione
  jx = ione
  jy = ione
  m = desc_a%get_global_rows()
  lldx = size(x,1)
  lldy = size(y,1)
  ! check vector correctness
  call psb_chkvect(m,lone,lldx,ix,jx,desc_a,info,iix,jjx)
  if (info == psb_success_)&
       & call psb_chkvect(m,lone,lldy,iy,jy,desc_a,info,iiy,jjy)
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

  nr = desc_a%get_local_rows() 
  if(nr > 0) then
    res = zdotc(int(nr,kind=psb_mpk_), x,1,y,1)
    ! adjust res because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx  = desc_a%ovrlap_elem(i,1)
      ndm  = desc_a%ovrlap_elem(i,2)
      res = res - (real(ndm-1)/real(ndm))*(conjg(x(idx))*y(idx))
    end do
  else
    res = zzero
  end if

  ! compute global sum
  if (global_) call psb_sum(ctxt, res)


  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

  return
end function psb_zdotv



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
! Subroutine: psb_zdotvs
!    psb_zdotvs computes the dot product of two distributed vectors,
!
!    res := X**C * Y
!
! Arguments:
!    res    -  complex.             The result.
!    x(:)   -  complex              The input vector containing the entries of X.
!    y(:)   -  complex              The input vector containing the entries of Y.
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!    global -  logical(optional)     Whether to perform the global sum, default: .true.
!
subroutine psb_zdotvs(res, x, y,desc_a, info,global)  
  use psb_base_mod, psb_protect_name => psb_zdotvs
  implicit none

  complex(psb_dpk_), intent(in)    :: x(:), y(:)
  complex(psb_dpk_), intent(out)   :: res
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)   :: info
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, idx, ndm,&
       & err_act, iix, jjx, iiy, jjy, i,nr, lldx, lldy
  integer(psb_lpk_) :: ix, jx, iy, jy, m
  logical :: global_
  complex(psb_dpk_)        :: zdotc
  character(len=20)        :: name, ch_err

  name='psb_zdot'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  ix = ione
  iy = ione
  m = desc_a%get_global_rows()
  lldx = size(x,1)
  lldy = size(y,1)
  ! check vector correctness
  call psb_chkvect(m,lone,lldx,ix,ix,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,lone,lldy,iy,iy,desc_a,info,iiy,jjy)
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

  nr = desc_a%get_local_rows() 
  if(nr > 0) then
    res = zdotc(int(nr,kind=psb_mpk_), x,1,y,1)
    ! adjust res because overlapped elements are computed more than once
    do i=1,size(desc_a%ovrlap_elem,1)
      idx  = desc_a%ovrlap_elem(i,1)
      ndm  = desc_a%ovrlap_elem(i,2)
      res = res - (real(ndm-1)/real(ndm))*(conjg(x(idx))*y(idx))
    end do
  else
    res = zzero
  end if

  ! compute global sum
  if (global_) call psb_sum(ctxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psb_zdotvs




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
! Subroutine: psb_zmdots
!    psb_zmdots computes the dot product of multiple distributed vectors,
!
!    res(i) := ( X(:,i) )**C * ( Y(:,i) )
!
! Arguments:
!    res(:) -  complex.             The result.
!    x(:)   -  complex              The input vector containing the entries of sub( X ).
!    y(:)   -  complex              The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type). The communication descriptor.
!    info   -  integer.             Return code
!    global -  logical(optional)     Whether to perform the global sum, default: .true.
!
subroutine psb_zmdots(res, x, y, desc_a, info,global)  
  use psb_base_mod, psb_protect_name => psb_zmdots
  implicit none

  complex(psb_dpk_), intent(in)    :: x(:,:), y(:,:)
  complex(psb_dpk_), intent(out)   :: res(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)   :: info
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, idx, ndm,&
       & err_act, iix, jjx, iiy, jjy, i, j, k, nr, lldx, lldy
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
  logical :: global_
  complex(psb_dpk_)        :: zdotc
  character(len=20)        :: name, ch_err

  name='psb_zmdots'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if
  ix = ione
  iy = ione

  m = desc_a%get_global_rows()
  lldx = size(x,1)
  lldy = size(y,1)

  ! check vector correctness
  call psb_chkvect(m,lone,lldx,ix,ix,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,lone,lldy,iy,iy,desc_a,info,iiy,jjy)
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

  nr = desc_a%get_local_rows() 
  if(nr > 0) then
    do j=1,k
      res(j) = zdotc(int(nr,kind=psb_mpk_),x(1:,j),1,y(1:,j),1)
      ! adjust res because overlapped elements are computed more than once
    end do
    do i=1,size(desc_a%ovrlap_elem,1)
      idx  = desc_a%ovrlap_elem(i,1)
      ndm  = desc_a%ovrlap_elem(i,2)
      res(1:k) = res(1:k) - &
           & (real(ndm-1)/real(ndm))*(conjg(x(idx,1:k))*y(idx,1:k))
    end do
  else
    res(:) = zzero
  end if


  ! compute global sum
  if (global_) call psb_sum(ctxt, res(1:k))

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psb_zmdots
