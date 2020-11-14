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
! Subroutine: psi_zovrl_update
!   These subroutines update the overlap region of a vector; they are  used
!   for the transpose  matrix-vector product when there is a nonempty overlap,
!   or for the application of Additive Schwarz preconditioners.                                           
!    
!    
subroutine  psi_zovrl_updr1(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_zovrl_updr1

  implicit none

  complex(psb_dpk_), intent(inout), target :: x(:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer(psb_ipk_), intent(in)                     :: update
  integer(psb_ipk_), intent(out)                    :: info

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me, err_act, i, idx, ndm
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name='psi_zovrl_updr1'
  info = psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
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
    ierr(1) = 3; ierr(2)=update;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_zovrl_updr1

subroutine  psi_zovrl_updr2(x,desc_a,update,info)
  use psi_mod, psi_protect_name =>   psi_zovrl_updr2

  implicit none

  complex(psb_dpk_), intent(inout), target :: x(:,:)
  type(psb_desc_type), intent(in)         :: desc_a
  integer(psb_ipk_), intent(in)                     :: update
  integer(psb_ipk_), intent(out)                    :: info

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me, err_act, i, idx, ndm
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name='psi_zovrl_updr2'
  info = psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
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
    ierr(1) = 3; ierr(2)=update;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end select

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_zovrl_updr2
