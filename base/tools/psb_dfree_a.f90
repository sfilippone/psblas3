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
! File: psb_dfree.f90
!
! Subroutine: psb_dfree
!    frees a dense matrix structure
! 
! Arguments: 
!    x(:,:)   - real, allocatable          The dense matrix to be freed.
!    desc_a   - type(psb_desc_type).        The communication descriptor.
!    info     - integer.                      Return code
subroutine psb_dfree(x, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_dfree
  implicit none

  !....parameters...
  real(psb_dpk_),allocatable, intent(inout)    :: x(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  !...locals....
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me, err_act
  character(len=20)   :: name

  name='psb_dfree'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  
  if (.not.psb_is_ok_desc(desc_a)) then
    info=psb_err_forgot_spall_
    call psb_errpush(info,name)
    return
  end if

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.allocated(x)) then
    info=psb_err_forgot_spall_
    call psb_errpush(info,name)
    goto 9999
  end if

  !deallocate x
  deallocate(x,stat=info)
  if (info /= psb_no_err_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_dfree



! Subroutine: psb_dfreev
!    frees a dense matrix structure
! 
! Arguments: 
!    x(:)     - real, allocatable        The dense matrix to be freed.
!    desc_a   - type(psb_desc_type).      The communication descriptor.
!    info     - integer.                    Return code
subroutine psb_dfreev(x, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_dfreev
  implicit none
  !....parameters...
  real(psb_dpk_),allocatable, intent(inout)    :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  !...locals....
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me, err_act
  character(len=20)   :: name

  name='psb_dfreev'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if (.not.psb_is_ok_desc(desc_a)) then
    info=psb_err_forgot_spall_
    call psb_errpush(info,name)
    goto 9999
  end if
  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999

  endif

  if (.not.allocated(x)) then
    info=psb_err_forgot_spall_
    call psb_errpush(info,name)
    goto 9999
  end if

  !deallocate x
  deallocate(x,stat=info)
  if (info /= psb_no_err_) then
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_dfreev
