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
! File: psb_smlt_vect

subroutine psb_smlt_vect(x,y,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_smlt_vect
  implicit none
  type(psb_s_vect_type), intent (inout) :: x
  type(psb_s_vect_type), intent (inout) :: y
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
  character(len=20)        :: name, ch_err

  name='psb_s_mlt_vect'
  if (psb_errstatus_fatal()) return
  info=psb_success_
  call psb_erractionsave(err_act)

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


  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m,lone,x%get_nrows(),ix,lone,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 1'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,lone,y%get_nrows(),iy,lone,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 2'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(desc_a%get_local_rows() > 0) then
    call y%mlt(x,info)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_smlt_vect

!
! Subroutine: psb_smlt_vect2
!

subroutine psb_smlt_vect2(alpha,x,y,beta,z,desc_a,info,conjgx, conjgy)
  use psb_base_mod, psb_protect_name => psb_smlt_vect2
  implicit none
  real(psb_spk_), intent(in)        :: alpha,beta
  type(psb_s_vect_type), intent (inout) :: x
  type(psb_s_vect_type), intent (inout) :: y
  type(psb_s_vect_type), intent (inout) :: z
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  character(len=1), intent(in), optional :: conjgx, conjgy

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_) :: ix, ijx, iy, ijy, iz, ijz, m
  character(len=20)        :: name, ch_err

  name='psb_s_mlt_vect2'
  if (psb_errstatus_fatal()) return
  info=psb_success_
  call psb_erractionsave(err_act)

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
  if (.not.allocated(z%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  ix = ione
  iy = ione
  iz = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m,lone,x%get_nrows(),ix,lone,desc_a,info,iix,jjx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 1'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,lone,y%get_nrows(),iy,lone,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 2'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_chkvect(m,lone,z%get_nrows(),iz,lone,desc_a,info,iiz,jjz)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 3'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(desc_a%get_local_rows() > 0) then
    call z%mlt(alpha,x,y,beta,info,conjgx,conjgy)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_smlt_vect2

subroutine psb_smlt_mltvec_va(x,a,v,desc,info)
  use psb_base_mod, psb_protect_name => psb_smlt_mltvec_va
  implicit none
  type(psb_s_multivect_type), intent(inout) :: x
  real(psb_spk_), dimension(:), allocatable, intent(inout) :: a
  type(psb_s_vect_type), intent(inout) :: v
  type(psb_desc_type), intent(inout) :: desc
  integer(psb_ipk_), intent(inout) :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_) :: ix, ijx, iy, ijy, iz, ijz, m
  character(len=20)        :: name, ch_err

  name='psb_s_mlt_mltvec_va'
  if (psb_errstatus_fatal()) return
  info=psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc%get_ctxt()

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
  if (.not.allocated(v%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(a)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if(desc%get_local_rows() > 0) then
    call x%mlt(a,v,info)
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_smlt_mltvec_va
