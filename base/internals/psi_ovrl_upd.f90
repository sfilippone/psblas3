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

subroutine  psi_sovrl_updr1(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_sovrl_updr1

  implicit none

  real(psb_spk_), intent(inout), target :: x(:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(in)                     :: update
  integer, intent(out)                    :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_sovrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
  case(psb_square_root_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/sqrt(real(ndm))
    end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx) = szero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_sovrl_updr1


subroutine  psi_sovrl_updr2(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_sovrl_updr2

  implicit none

  real(psb_spk_), intent(inout), target :: x(:,:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(in)                     :: update
  integer, intent(out)                    :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_sovrl_updr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
  case(psb_square_root_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/sqrt(real(ndm))
    end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx,:) = szero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_sovrl_updr2

subroutine  psi_dovrl_updr1(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_updr1

  implicit none

  real(psb_dpk_), intent(inout), target :: x(:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(in)                     :: update
  integer, intent(out)                    :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_dovrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
  case(psb_square_root_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/sqrt(real(ndm))
    end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx) = dzero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_dovrl_updr1


subroutine  psi_dovrl_updr2(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_updr2

  implicit none

  real(psb_dpk_), intent(inout), target :: x(:,:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(in)                     :: update
  integer, intent(out)                    :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_dovrl_updr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
  case(psb_square_root_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/sqrt(real(ndm))
    end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx,:) = dzero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_dovrl_updr2

subroutine  psi_covrl_updr1(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_covrl_updr1

  implicit none

  complex(psb_spk_), intent(inout), target :: x(:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(in)                     :: update
  integer, intent(out)                    :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_covrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
  case(psb_square_root_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/sqrt(real(ndm))
    end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx) = czero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_covrl_updr1


subroutine  psi_covrl_updr2(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_covrl_updr2

  implicit none

  complex(psb_spk_), intent(inout), target :: x(:,:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(in)                     :: update
  integer, intent(out)                    :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_covrl_updr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
  case(psb_square_root_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/sqrt(real(ndm))
    end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx,:) = czero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_covrl_updr2

subroutine  psi_zovrl_updr1(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_zovrl_updr1

  implicit none

  complex(psb_dpk_), intent(inout), target :: x(:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(in)                     :: update
  integer, intent(out)                    :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_zovrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
  case(psb_square_root_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/sqrt(real(ndm))
    end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx) = zzero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_zovrl_updr1


subroutine  psi_zovrl_updr2(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_zovrl_updr2

  implicit none

  complex(psb_dpk_), intent(inout), target :: x(:,:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(in)                     :: update
  integer, intent(out)                    :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_zovrl_updr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
  case(psb_square_root_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/sqrt(real(ndm))
    end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx,:) = zzero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_zovrl_updr2

subroutine  psi_iovrl_updr1(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_iovrl_updr1

  implicit none

  integer, intent(inout), target   :: x(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(in)              :: update
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_iovrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
    ! Square root does not make sense here
!!$    case(psb_square_root_)
!!$      do i=1,size(desc_a%ovrlap_elem,1)
!!$        idx = desc_a%ovrlap_elem(i,1)
!!$        ndm = desc_a%ovrlap_elem(i,2)
!!$        x(idx) = x(idx)/sqrt(real(ndm))
!!$      end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx) = x(idx)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx) = izero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_iovrl_updr1


subroutine  psi_iovrl_updr2(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_iovrl_updr2

  implicit none

  integer, intent(inout), target   :: x(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(in)              :: update
  integer, intent(out)             :: info

  ! locals
  integer           :: ictxt, np, me, err_act, i, idx, ndm
  character(len=20) :: name, ch_err

  name='psi_iovrl_updr2'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! switch on update type
  select case (update)
    ! Square root does not make sense here
!!$    case(psb_square_root_)
!!$      do i=1,size(desc_a%ovrlap_elem,1)
!!$        idx = desc_a%ovrlap_elem(i,1)
!!$        ndm = desc_a%ovrlap_elem(i,2)
!!$        x(idx,:) = x(idx,:)/sqrt(real(ndm))
!!$      end do
  case(psb_avg_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      ndm = desc_a%ovrlap_elem(i,2)
      x(idx,:) = x(idx,:)/real(ndm)
    end do
  case(psb_setzero_)
    do i=1,size(desc_a%ovrlap_elem,1)
      idx = desc_a%ovrlap_elem(i,1)
      if (me /= desc_a%ovrlap_elem(i,3))&
           & x(idx,:) = izero
    end do
  case(psb_sum_)
    ! do nothing

  case default 
    ! wrong value for choice argument
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psi_iovrl_updr2


subroutine  psi_sovrl_upd_vect(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_sovrl_upd_vect
  use psb_realloc_mod
  use psb_s_base_vect_mod

  implicit none

  class(psb_s_base_vect_type)     :: x
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(in)             :: update
  integer, intent(out)            :: info

  ! locals
  real(psb_spk_), allocatable :: xs(:)
  integer           :: ictxt, np, me, err_act, i, idx, ndm, nx
  character(len=20) :: name, ch_err


  name='psi_sovrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  nx = size(desc_a%ovrlap_elem,1)
  call psb_realloc(nx,xs,info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_Dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (update /= psb_sum_) then 
    call x%gth(nx,desc_a%ovrlap_elem(:,1),xs)
    ! switch on update type

    select case (update)
    case(psb_square_root_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i) = xs(i)/sqrt(real(ndm))
      end do
    case(psb_avg_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i) = xs(i)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,nx
        if (me /= desc_a%ovrlap_elem(i,3))&
             & xs(i) = szero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = psb_err_iarg_invalid_value_
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select
    call x%sct(nx,desc_a%ovrlap_elem(:,1),xs,szero)
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
end subroutine psi_sovrl_upd_vect

subroutine  psi_dovrl_upd_vect(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_upd_vect
  use psb_realloc_mod
  use psb_d_base_vect_mod

  implicit none

  class(psb_d_base_vect_type)     :: x
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(in)             :: update
  integer, intent(out)            :: info

  ! locals
  real(psb_dpk_), allocatable :: xs(:)
  integer           :: ictxt, np, me, err_act, i, idx, ndm, nx
  character(len=20) :: name, ch_err


  name='psi_dovrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  nx = size(desc_a%ovrlap_elem,1)
  call psb_realloc(nx,xs,info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_Dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (update /= psb_sum_) then 
    call x%gth(nx,desc_a%ovrlap_elem(:,1),xs)
    ! switch on update type

    select case (update)
    case(psb_square_root_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i) = xs(i)/sqrt(dble(ndm))
      end do
    case(psb_avg_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i) = xs(i)/dble(ndm)
      end do
    case(psb_setzero_)
      do i=1,nx
        if (me /= desc_a%ovrlap_elem(i,3))&
             & xs(i) = dzero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = psb_err_iarg_invalid_value_
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select
    call x%sct(nx,desc_a%ovrlap_elem(:,1),xs,dzero)
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
end subroutine psi_dovrl_upd_vect


subroutine  psi_covrl_upd_vect(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_covrl_upd_vect
  use psb_realloc_mod
  use psb_c_base_vect_mod

  implicit none

  class(psb_c_base_vect_type)     :: x
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(in)             :: update
  integer, intent(out)            :: info

  ! locals
  complex(psb_spk_), allocatable :: xs(:)
  integer           :: ictxt, np, me, err_act, i, idx, ndm, nx
  character(len=20) :: name, ch_err


  name='psi_covrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  nx = size(desc_a%ovrlap_elem,1)
  call psb_realloc(nx,xs,info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_Dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (update /= psb_sum_) then 
    call x%gth(nx,desc_a%ovrlap_elem(:,1),xs)
    ! switch on update type

    select case (update)
    case(psb_square_root_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i) = xs(i)/sqrt(real(ndm))
      end do
    case(psb_avg_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i) = xs(i)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,nx
        if (me /= desc_a%ovrlap_elem(i,3))&
             & xs(i) = szero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = psb_err_iarg_invalid_value_
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select
    call x%sct(nx,desc_a%ovrlap_elem(:,1),xs,czero)
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
end subroutine psi_covrl_upd_vect

subroutine  psi_zovrl_upd_vect(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_zovrl_upd_vect
  use psb_realloc_mod
  use psb_z_base_vect_mod

  implicit none

  class(psb_z_base_vect_type)     :: x
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(in)             :: update
  integer, intent(out)            :: info

  ! locals
  complex(psb_dpk_), allocatable :: xs(:)
  integer           :: ictxt, np, me, err_act, i, idx, ndm, nx
  character(len=20) :: name, ch_err


  name='psi_zovrl_updr1'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  nx = size(desc_a%ovrlap_elem,1)
  call psb_realloc(nx,xs,info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_Dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (update /= psb_sum_) then 
    call x%gth(nx,desc_a%ovrlap_elem(:,1),xs)
    ! switch on update type

    select case (update)
    case(psb_square_root_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i) = xs(i)/sqrt(dble(ndm))
      end do
    case(psb_avg_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i) = xs(i)/dble(ndm)
      end do
    case(psb_setzero_)
      do i=1,nx
        if (me /= desc_a%ovrlap_elem(i,3))&
             & xs(i) = dzero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = psb_err_iarg_invalid_value_
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select
    call x%sct(nx,desc_a%ovrlap_elem(:,1),xs,zzero)
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
end subroutine psi_zovrl_upd_vect


