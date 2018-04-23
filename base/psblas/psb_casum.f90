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
! File: psb_casum.f90
!
! Function: psb_casum 
!    Computes norm1 of X
!
!    norm1 := sum(sub( X )(i))
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Arguments:
!    x(:,:) -  complex                The input vector.
!    desc_a -  type(psb_desc_type).   The communication descriptor.
!    info   -  integer.               Return code
!    jx     -  integer(optional).     The column offset.
!
function psb_casum (x,desc_a, info, jx,global) result(res)
  use psb_base_mod, psb_protect_name => psb_casum

  implicit none

  complex(psb_spk_), intent(in)   :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_), optional, intent(in)     :: jx
  real(psb_spk_)                  :: res
  logical, intent(in), optional        :: global

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, &
       & err_act, iix, jjx, ix, ijx, m, i, idx, ndm, ldx
  logical :: global_
  character(len=20)        :: name, ch_err

  name='psb_casum'
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

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

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
  if (global_) call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end function psb_casum


function psb_casum_vect(x, desc_a, info,global) result(res)
  use psb_base_mod, psb_protect_name => psb_casum_vect
  implicit none

  real(psb_spk_)                        :: res
  type(psb_c_vect_type), intent (inout) :: x
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  logical, intent(in), optional        :: global

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, jx, ix, m, imax, i, idx, ndm  
  logical :: global_
  character(len=20)        :: name, ch_err

  name='psb_casumv'
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

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

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
    if (size(desc_a%ovrlap_elem,1)>0) then
      if (x%is_dev()) call x%sync()
      ! adjust res because overlapped elements are computed more than once
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        res = res - (real(ndm-1)/real(ndm))*abs(x%v%v(idx))
      end do
    end if
  else 
    res = szero
  end if

  ! compute global sum
  if (global_) call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return

end function psb_casum_vect



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
! Function: psb_casumv 
!    Computes norm1 of X
!
!    norm1 := sum(X(i))
!
! Arguments:
!    x(:)   -  complex               The input vector.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!
function psb_casumv(x,desc_a, info,global) result(res)
  use psb_base_mod, psb_protect_name => psb_casumv

  implicit none

  complex(psb_spk_), intent(in)   :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info
  real(psb_spk_)                  :: res
  logical, intent(in), optional        :: global

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, jx, ix, m, i, idx, ndm, ldx
  logical :: global_
  character(len=20)        :: name, ch_err

  name='psb_casumv'
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

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

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
  if (global_) call psb_sum(ictxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end function psb_casumv


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
! Subroutine: psb_casumvs
!    Computes norm1 of X
!
!    norm1 := sum(X(i))
!
! Arguments:
!    res    -  real.                 The result.
!    x(:)   -  complex               The input vector.
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    jx     -  integer(optional).    The column offset.
!
subroutine psb_casumvs(res,x,desc_a, info,global)
  use psb_base_mod, psb_protect_name => psb_casumvs

  implicit none

  complex(psb_spk_), intent(in)    :: x(:)
  real(psb_spk_), intent(out)      :: res
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)   :: info
  logical, intent(in), optional        :: global

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ix, jx, m, i, idx, ndm, ldx
  logical :: global_
  character(len=20)        :: name, ch_err

  name='psb_casumvs'
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

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

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
  if (global_) call psb_sum(ictxt,res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_casumvs
