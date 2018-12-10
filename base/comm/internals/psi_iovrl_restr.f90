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

subroutine  psi_iovrl_restr_vect(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_iovrl_restr_vect
  use psb_i_base_vect_mod

  implicit none

  class(psb_i_base_vect_type)     :: x
  integer(psb_ipk_)               :: xs(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, err_act, i, idx, isz
  character(len=20) :: name, ch_err

  name='psi_iovrl_restr_vect'
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

  isz = size(desc_a%ovrlap_elem,1)

  call x%sct(isz,desc_a%ovrlap_elem(:,1),xs,izero)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_iovrl_restr_vect


subroutine  psi_iovrl_restr_multivect(x,xs,desc_a,info)
  use psi_mod, psi_protect_name =>   psi_iovrl_restr_multivect
  use psb_i_base_vect_mod

  implicit none

  class(psb_i_base_multivect_type)     :: x
  integer(psb_ipk_)               :: xs(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, err_act, i, idx, isz, nc
  character(len=20) :: name, ch_err

  name='psi_iovrl_restr_mv'
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

  isz = size(desc_a%ovrlap_elem,1)
  call x%sct(isz,desc_a%ovrlap_elem(:,1),xs,izero)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_iovrl_restr_multivect


