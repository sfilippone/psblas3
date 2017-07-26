!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006, 2010, 2015, 2017
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
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
subroutine  psi_dovrl_updr1(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_updr1

  implicit none

  real(psb_dpk_), intent(inout), target :: x(:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer(psb_ipk_), intent(in)                     :: update
  integer(psb_ipk_), intent(out)                    :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, err_act, i, idx, ndm
  integer(psb_ipk_) :: ierr(5)
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
    ierr(1) = 3; ierr(2)=update;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_dovrl_updr1


subroutine  psi_dovrl_updr2(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_updr2

  implicit none

  real(psb_dpk_), intent(inout), target :: x(:,:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer(psb_ipk_), intent(in)                     :: update
  integer(psb_ipk_), intent(out)                    :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, err_act, i, idx, ndm
  integer(psb_ipk_) :: ierr(5)
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
    ierr(1) = 3; ierr(2)=update;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_dovrl_updr2


subroutine  psi_dovrl_upd_vect(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_upd_vect
  use psb_realloc_mod
  use psb_d_base_vect_mod

  implicit none

  class(psb_d_base_vect_type)     :: x
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(in)             :: update
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  real(psb_dpk_), allocatable :: xs(:)
  integer(psb_ipk_) :: ictxt, np, me, err_act, i, idx, ndm, nx
  integer(psb_ipk_) :: ierr(5)
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
             & xs(i) = dzero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = psb_err_iarg_invalid_value_
      ierr(1) = 3; ierr(2)=update;
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end select
    call x%sct(nx,desc_a%ovrlap_elem(:,1),xs,dzero)
  end if

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_dovrl_upd_vect

subroutine  psi_dovrl_upd_multivect(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_dovrl_upd_multivect
  use psb_realloc_mod
  use psb_d_base_vect_mod

  implicit none

  class(psb_d_base_multivect_type)     :: x
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(in)             :: update
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  real(psb_dpk_), allocatable :: xs(:,:)
  integer(psb_ipk_) :: ictxt, np, me, err_act, i, idx, ndm, nx, nc
  integer(psb_ipk_) :: ierr(5)
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
  nc = x%get_ncols()
  call psb_realloc(nx,nc,xs,info)
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
        xs(i,:) = xs(i,:)/sqrt(real(ndm))
      end do
    case(psb_avg_)
      do i=1,nx
        ndm = desc_a%ovrlap_elem(i,2)
        xs(i,:) = xs(i,:)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,nx
        if (me /= desc_a%ovrlap_elem(i,3))&
             & xs(i,:) = dzero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = psb_err_iarg_invalid_value_
      ierr(1) = 3; ierr(2)=update;
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end select
    call x%sct(nx,desc_a%ovrlap_elem(:,1),xs,dzero)
  end if

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_dovrl_upd_multivect
