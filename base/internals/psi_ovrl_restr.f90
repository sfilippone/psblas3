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
subroutine  psi_sovrl_restrr1(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_sovrl_restrr1

  implicit none

  real(psb_spk_), intent(inout)  :: x(:)
  real(psb_spk_)                 :: xs(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_sovrl_restrr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx    = desc_a%ovrlap_elem(i,1)
    x(idx) = xs(i) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_sovrl_restrr1

subroutine  psi_sovrl_restrr2(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_sovrl_restrr2

  implicit none

  real(psb_spk_), intent(inout)  :: x(:,:)
  real(psb_spk_)                 :: xs(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_sovrl_restrr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (size(x,2) /= size(xs,2)) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name, a_err='Mismacth columns X vs XS')
    goto 9999
  endif


  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx      = desc_a%ovrlap_elem(i,1)
    x(idx,:) = xs(i,:) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_sovrl_restrr2


subroutine  psi_dovrl_restrr1(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_restrr1

  implicit none

  real(psb_dpk_), intent(inout)  :: x(:)
  real(psb_dpk_)                 :: xs(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_dovrl_restrr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx    = desc_a%ovrlap_elem(i,1)
    x(idx) = xs(i) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_dovrl_restrr1


subroutine  psi_dovrl_restrr2(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_restrr2

  implicit none

  real(psb_dpk_), intent(inout)  :: x(:,:)
  real(psb_dpk_)                 :: xs(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_dovrl_restrr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (size(x,2) /= size(xs,2)) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name, a_err='Mismacth columns X vs XS')
    goto 9999
  endif


  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx      = desc_a%ovrlap_elem(i,1)
    x(idx,:) = xs(i,:) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_dovrl_restrr2


subroutine  psi_covrl_restrr1(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_covrl_restrr1

  implicit none

  complex(psb_spk_), intent(inout)  :: x(:)
  complex(psb_spk_)                 :: xs(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_covrl_restrr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx    = desc_a%ovrlap_elem(i,1)
    x(idx) = xs(i) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_covrl_restrr1

subroutine  psi_covrl_restrr2(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_covrl_restrr2

  implicit none

  complex(psb_spk_), intent(inout)  :: x(:,:)
  complex(psb_spk_)                 :: xs(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_covrl_restrr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (size(x,2) /= size(xs,2)) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name, a_err='Mismacth columns X vs XS')
    goto 9999
  endif


  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx      = desc_a%ovrlap_elem(i,1)
    x(idx,:) = xs(i,:) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_covrl_restrr2


subroutine  psi_zovrl_restrr1(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_zovrl_restrr1

  implicit none

  complex(psb_dpk_), intent(inout)  :: x(:)
  complex(psb_dpk_)                 :: xs(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_zovrl_restrr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx    = desc_a%ovrlap_elem(i,1)
    x(idx) = xs(i) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_zovrl_restrr1

subroutine  psi_zovrl_restrr2(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_zovrl_restrr2

  implicit none

  complex(psb_dpk_), intent(inout)  :: x(:,:)
  complex(psb_dpk_)                 :: xs(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_zovrl_restrr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (size(x,2) /= size(xs,2)) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name, a_err='Mismacth columns X vs XS')
    goto 9999
  endif


  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx      = desc_a%ovrlap_elem(i,1)
    x(idx,:) = xs(i,:) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_zovrl_restrr2


subroutine  psi_iovrl_restrr1(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_iovrl_restrr1

  implicit none

  integer, intent(inout)  :: x(:)
  integer                 :: xs(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_iovrl_restrr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx    = desc_a%ovrlap_elem(i,1)
    x(idx) = xs(i) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_iovrl_restrr1

subroutine  psi_iovrl_restrr2(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_iovrl_restrr2

  implicit none

  integer, intent(inout)  :: x(:,:)
  integer                 :: xs(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_iovrl_restrr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (size(x,2) /= size(xs,2)) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name, a_err='Mismacth columns X vs XS')
    goto 9999
  endif


  isz = size(desc_a%ovrlap_elem,1)

  do i=1, isz
    idx      = desc_a%ovrlap_elem(i,1)
    x(idx,:) = xs(i,:) 
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_iovrl_restrr2
