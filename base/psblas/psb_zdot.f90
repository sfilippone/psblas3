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
function psb_zdot_vect(x, y, desc_a,info) result(res)
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

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, ijx, iy, ijy, iiy, jjy, i, m, nr
  character(len=20)      :: name, ch_err

  name='psb_sdot'
  res = zzero
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
  ijx = ione

  iy = ione
  ijy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m,ione,x%get_nrows(),ix,ijx,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,ione,y%get_nrows(),iy,ijy,desc_a,info,iiy,jjy)
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
!!$      do i=1,size(desc_a%ovrlap_elem,1)
!!$        idx  = desc_a%ovrlap_elem(i,1)
!!$        ndm  = desc_a%ovrlap_elem(i,2)
!!$        dot_local = dot_local - (real(ndm-1)/real(ndm))*(x(idx)*y(idx))
!!$      end do
  else
    res = zzero
  end if

  ! compute global sum
  call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return

end function psb_zdot_vect

function psb_zdot(x, y,desc_a, info, jx, jy)  result(res)
  use psb_base_mod, psb_protect_name => psb_zdot
  implicit none

  complex(psb_dpk_), intent(in)    :: x(:,:), y(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(in), optional    :: jx, jy
  integer(psb_ipk_), intent(out)   :: info
  complex(psb_dpk_)              :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, ijx, iy, ijy, iiy, jjy, i, m, nr, &
       & lldx, lldy
  complex(psb_dpk_)        :: zdotc
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

  nr = desc_a%get_local_rows() 
  if(nr > 0) then
    res = zdotc(int(nr,kind=psb_mpik_), x(iix:,jjx),1,y(iiy:,jjy),1)
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
  call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

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
function psb_zdotv(x, y,desc_a, info)  result(res)
  use psb_base_mod, psb_protect_name => psb_zdotv
  implicit none

  complex(psb_dpk_), intent(in)   :: x(:), y(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info
  complex(psb_dpk_)              :: res

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, jx, iy, jy, iiy, jjy, i, m, nr, &
       & lldx, lldy
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
  lldx = size(x,1)
  lldy = size(y,1)
  ! check vector correctness
  call psb_chkvect(m,ione,lldx,ix,jx,desc_a,info,iix,jjx)
  if (info == psb_success_)&
       & call psb_chkvect(m,ione,lldy,iy,jy,desc_a,info,iiy,jjy)
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
    res = zdotc(int(nr,kind=psb_mpik_), x,1,y,1)
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
  call psb_sum(ictxt, res)


  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

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

  complex(psb_dpk_), intent(in)    :: x(:), y(:)
  complex(psb_dpk_), intent(out)   :: res
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)   :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, iy, iiy, jjy, i, m,nr, &
       & lldx, lldy
  complex(psb_dpk_)        :: zdotc
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
  lldx = size(x,1)
  lldy = size(y,1)
  ! check vector correctness
  call psb_chkvect(m,ione,lldx,ix,ix,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,ione,lldy,iy,iy,desc_a,info,iiy,jjy)
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
    res = zdotc(int(nr,kind=psb_mpik_), x,1,y,1)
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
  call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

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

  complex(psb_dpk_), intent(in)    :: x(:,:), y(:,:)
  complex(psb_dpk_), intent(out)   :: res(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)   :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, idx, ndm,&
       & err_act, iix, jjx, ix, iy, iiy, jjy, i, m, j, k, nr, &
       & lldx, lldy
  complex(psb_dpk_)        :: zdotc
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
  lldx = size(x,1)
  lldy = size(y,1)

  ! check vector correctness
  call psb_chkvect(m,ione,lldx,ix,ix,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,ione,lldy,iy,iy,desc_a,info,iiy,jjy)
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
      res(j) = zdotc(int(nr,kind=psb_mpik_),x(1:,j),1,y(1:,j),1)
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
  call psb_sum(ictxt, res(1:k))

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_zmdots
