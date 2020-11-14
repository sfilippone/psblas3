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
! File: psb_sdiv_vect

subroutine psb_sdiv_vect(x,y,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_sdiv_vect
  implicit none
  type(psb_s_vect_type), intent (inout) :: x
  type(psb_s_vect_type), intent (inout) :: y
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
  character(len=20)        :: name, ch_err

  name='psb_s_div_vect'
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
    call x%div(y,info)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sdiv_vect

subroutine psb_sdiv_vect2(x,y,z,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_sdiv_vect2
  implicit none
  type(psb_s_vect_type), intent (inout) :: x
  type(psb_s_vect_type), intent (inout) :: y
  type(psb_s_vect_type), intent (inout) :: z
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_) :: ix, ijx, iy, ijy, iz, ijz, m
  character(len=20)        :: name, ch_err

  name='psb_s_div_vect2'
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
    ch_err='psb_chkvect 2'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(desc_a%get_local_rows() > 0) then
    call z%div(x,y,info)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sdiv_vect2

subroutine psb_sdiv_vect_check(x,y,desc_a,info,flag)
  use psb_base_mod, psb_protect_name => psb_sdiv_vect_check
  implicit none
  type(psb_s_vect_type), intent (inout) :: x
  type(psb_s_vect_type), intent (inout) :: y
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  logical, intent(in)                   :: flag

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
  character(len=20)        :: name, ch_err

  name='psb_s_div_vect_check'
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
    call x%div(y,info,flag)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sdiv_vect_check

subroutine psb_sdiv_vect2_check(x,y,z,desc_a,info,flag)
  use psb_base_mod, psb_protect_name => psb_sdiv_vect2_check
  implicit none
  type(psb_s_vect_type), intent (inout) :: x
  type(psb_s_vect_type), intent (inout) :: y
  type(psb_s_vect_type), intent (inout) :: z
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  logical, intent(in)                   :: flag

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_) :: ix, ijx, iy, ijy, iz, ijz, m
  character(len=20)        :: name, ch_err

  name='psb_s_div_vect2_check'
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
    ch_err='psb_chkvect 2'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(desc_a%get_local_rows() > 0) then
    call z%div(x,y,info,flag)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sdiv_vect2_check

function psb_sminquotient_vect(x,y,desc_a,info,global) result(res)
  use psb_penv_mod
  use psb_serial_mod
  use psb_desc_mod
  use psb_check_mod
  use psb_error_mod
  use psb_s_vect_mod
  implicit none

  real(psb_spk_)                        :: res
  type(psb_s_vect_type), intent (inout) :: x
  type(psb_s_vect_type), intent (inout) :: y
  type(psb_desc_type), intent (in)        :: desc_a
  integer(psb_ipk_), intent(out)          :: info
  logical, intent(in), optional           :: global

  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx
  integer(psb_lpk_) :: ix, jx, iy, ijy, m
  logical :: global_
  character(len=20)      :: name, ch_err

  name='psb_sminquotient_vect'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

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
  call psb_chkvect(m,lone,x%get_nrows(),ix,jx,desc_a,info,iix,jjx)
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
    res = x%minquotient(y,info)
  else
    res = szero
  end if

  ! compute global min
  if (global_) call psb_min(ictxt, res)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end function psb_sminquotient_vect
