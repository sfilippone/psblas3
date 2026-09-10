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
!         software without specific prior written permission.
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
! File: psb_daxpby.f90

!
! Subroutine: psb_daxpby_vect
!    Adds one distributed vector to another, 
!
!    Y := beta * Y + alpha * X
!
! Arguments:
!    alpha  - real, input           The scalar used to multiply each component of X
!    x      - type(psb_d_vect_type) The input vector containing the entries of X
!    beta   - real, input           The scalar used to multiply each component of Y
!    y      - type(psb_d_vect_type)  The input/output vector Y
!    desc_a - type(psb_desc_type)    The communication descriptor.
!    info   - integer                Return code
!
!  Note: from a functional point of view, X is input, but here
!        it's declared INOUT because of the sync() methods.
!
subroutine psb_daxpby_vect(alpha, x, beta, y, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_vect
  implicit none
  real(psb_dpk_), intent(in)           :: alpha, beta
  type(psb_d_vect_type), intent(inout) :: x, y
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)       :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_dgeaxpby'
  if(psb_errstatus_fatal()) return
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif
  if(.not. allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif
  if(.not. allocated(y%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call y%axpby(desc_a%get_local_rows(), alpha, x, beta, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_vect


! Subroutines: psb_daxpby_multivect_*
!    Adds one distributed multivector/vector to another multivector/vector, 
!
!    Y(:) = alpha * X(:, j) + beta * Y(:)           (psb_daxpby_extract_c) 

!    Y(:, :) := beta * Y(:, :) + alpha * X          (psb_daxpby_mv_v_full)
!    Y(:, j) := beta * Y(:, j) + alpha * X          (psb_daxpby_mv_v_idxs)
!    Y(:, :) := beta * Y(:, :) + alpha * X(:, :)    (psb_daxpby_mv_m_full)
!    Y(:, j) := beta * Y(:, j) + alpha * X(:, k)    (psb_daxpby_mv_m_idxs)
!    Z(:, :) := beta * Y(:, :) + alpha * X(:, :)    (psb_daxpby_mv_m_full_out)
!
!    Z(:, k) := gamma * Z(:, k) + beta * Y + alpha * X                (psb_daxpby_mv_vv)
!    Z(:, k) := gamma * Z(:, k) + beta * Y(:, j) + alpha * X          (psb_daxpby_mv_mv)
!    Z(:, k) := gamma * Z(:, k) + beta * Y(:, j) + alpha * X(:, i)    (psb_daxpby_mv_mm)
!
!    Z(:, :) := gamma * Z(:, :) + beta * Y(:, :) + alpha * X(:, :)    (psb_daxpby_mv_mm_full)
!
! Arguments: ....
!
subroutine psb_daxpby_extract_c(alpha, x, idx_x, beta, y, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_extract_c
  implicit none
  real(psb_dpk_), intent(in)                :: alpha, beta
  type(psb_d_multivect_type), intent(inout) :: x
  integer(psb_ipk_), intent(in)             :: idx_x
  type(psb_d_vect_type), intent(inout)      :: y
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpy_extract_c'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if(.not. allocated(x%v)) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif
  if(.not. allocated(y%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call x%axpby(desc_a%get_local_rows(), alpha, idx_x, beta, y, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_extract_c

subroutine psb_daxpby_mv_v_full(alpha, x, beta, y, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_v_full
  implicit none
  real(psb_dpk_), intent(in)                :: alpha, beta
  type(psb_d_vect_type), intent(inout)      :: x
  type(psb_d_multivect_type), intent(inout) :: y
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_v_full'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if(.not. allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif
  if(.not. allocated(y%v)) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call y%axpby(desc_a%get_local_rows(), alpha, x, beta, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_v_full

subroutine psb_daxpby_mv_v_idxs(alpha, x, beta, y, idx_y, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_v_idxs
  implicit none
  type(psb_d_vect_type), intent(inout)      :: x
  type(psb_d_multivect_type), intent(inout) :: y
  integer(psb_ipk_), intent(in)             :: idx_y
  real(psb_dpk_), intent(in)                :: alpha, beta
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_v_idxs'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if(.not. allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif
  if(.not. allocated(y%v)) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    call y%axpby(desc_a%get_local_rows(), alpha, x, beta, idx_y, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_v_idxs

subroutine psb_daxpby_mv_m_full(alpha, x, beta, y, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_m_full
  implicit none
  type(psb_d_multivect_type), intent(inout) :: x, y
  real(psb_dpk_), intent(in)                :: alpha, beta
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_m_full'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call y%axpby(desc_a%get_local_rows(), alpha, x, beta, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_m_full

subroutine psb_daxpby_mv_m_idxs(alpha, x, idx_x, beta, y, idx_y, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_m_idxs
  implicit none
  type(psb_d_multivect_type), intent(inout) :: x, y
  integer(psb_ipk_), intent(in)             :: idx_x, idx_y
  real(psb_dpk_), intent(in)                :: alpha, beta
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_m_idxs'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if(.not. allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif
  if(.not. allocated(y%v)) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call y%axpby(desc_a%get_local_rows(), alpha, x, idx_x, beta, idx_y, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_m_idxs

subroutine psb_daxpby_mv_m_full_out(alpha, x, beta, y, z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_m_full_out
  implicit none
  type(psb_d_multivect_type), intent(inout) :: x, y, z
  real(psb_dpk_), intent(in)                :: alpha, beta
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, iz, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_mm_full'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione
  iz = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iz, lone, desc_a, info, iiz, jjz)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 3'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione) .or. (iiz /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) then
    call z%axpby(desc_a%get_local_rows(), alpha, x, beta, y, info)
  end if

  call psb_erractionrestore(err_act)
  return
9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_m_full_out

subroutine psb_daxpby_mv_vv(alpha, x, beta, y, gamma, z, idx_z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_vv
  implicit none
  type(psb_d_vect_type), intent(inout)      :: x, y
  type(psb_d_multivect_type), intent(inout) :: z
  integer(psb_ipk_), intent(in)             :: idx_z
  real(psb_dpk_), intent(in)                :: alpha, beta, gamma
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, iz, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_vv'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif
  if(.not. allocated(z%v)) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione
  iz = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iz, lone, desc_a, info, iiz, jjz)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 3'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione) .or. (iiz /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call z%axpby(desc_a%get_local_rows(), alpha, x, beta, y, gamma, idx_z, info)

  call psb_erractionrestore(err_act)
  return
9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_vv

subroutine psb_daxpby_mv_mv(alpha, x, beta, y, idx_y, gamma, z, idx_z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_mv
  implicit none
  type(psb_d_vect_type), intent(inout)      :: x
  type(psb_d_multivect_type), intent(inout) :: y, z
  integer(psb_ipk_), intent(in)             :: idx_y, idx_z
  real(psb_dpk_), intent(in)                :: alpha, beta, gamma
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, iz, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_mv'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if(.not. allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif
  if((.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione
  iz = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iz, lone, desc_a, info, iiz, jjz)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 3'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione) .or. (iiz /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) & 
    & call z%axpby(desc_a%get_local_rows(), alpha, x, beta, y, idx_y, gamma, idx_z, info)

  call psb_erractionrestore(err_act)
  return
9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_mv

subroutine psb_daxpby_mv_mm_idxs(alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_mm_idxs
  implicit none
  type(psb_d_multivect_type), intent(inout) :: x, y, z
  integer(psb_ipk_), intent(in)             :: idx_x, idx_y, idx_z
  real(psb_dpk_), intent(in)                :: alpha, beta, gamma
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, iz, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_mm_idxs'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione
  iz = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iz, lone, desc_a, info, iiz, jjz)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 3'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione) .or. (iiz /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call z%axpby(desc_a%get_local_rows(), alpha, x, idx_x, beta, y, idx_y, gamma, idx_z, info)

  call psb_erractionrestore(err_act)
  return
9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_mm_idxs

subroutine psb_daxpby_mv_mm_full(alpha, x, beta, y, gamma, z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_mm_full
  implicit none
  type(psb_d_multivect_type), intent(inout) :: x, y, z
  real(psb_dpk_), intent(in)                :: alpha, beta, gamma
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_) :: ix, ijx, iy, ijy, iz, m
  character(len=20) :: name, ch_err

  name = 'psb_daxpby_mv_mm_full'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione
  iz = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iz, lone, desc_a, info, iiz, jjz)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 3'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione) .or. (iiz /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call z%axpby(desc_a%get_local_rows(), alpha, x, beta, y, gamma, info)

  call psb_erractionrestore(err_act)
  return
9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_mm_full

subroutine psb_daxpby_mv_mm_out(alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, w, idx_w, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_mm_out
  implicit none
  type(psb_d_multivect_type), intent(inout) :: x, y, z, w
  integer(psb_ipk_), intent(in)             :: idx_x, idx_y, idx_z, idx_w
  real(psb_dpk_), intent(in)                :: alpha, beta, gamma
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy, iiz, jjz, iiw, jjw
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, iz, iw, ijw, m
  character(len=20)   :: name, ch_err

  name = 'psb_daxpby_mv_mm'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) & 
        .or. (.not. allocated(z%v)) .or. (.not. allocated(w%v))) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione
  iz = ione
  iw = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iz, lone, desc_a, info, iiz, jjz)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 3'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, w%get_nrows(), iw, lone, desc_a, info, iiw, jjw)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 3'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione) .or. (iiz /= ione) .or. (iiw /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) then
    call w%axpby(desc_a%get_local_rows(), alpha, x, idx_x, beta, y, idx_y, gamma, z, idx_z, idx_w, info)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_mm_out

subroutine psb_daxpby_mv_cspan1D(x, coeff, y, desc_a, info, upd_flag)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_cspan1D
  implicit none
  type(psb_d_multivect_type), intent(inout) :: x
  real(psb_dpk_), intent(in)                :: coeff(:)
  type(psb_d_vect_type), intent(inout)      :: y
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info
  logical, intent(in), optional             :: upd_flag

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err
  logical             :: upd_flag_

  upd_flag_ = .false.
  if(present(upd_flag)) upd_flag_ = upd_flag
  
  name = 'psb_daxpby_mv_cspan1D'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if(.not. allocated(x%v)) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif
  if(.not. allocated(y%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) then
    call x%axpby(desc_a%get_local_rows(), coeff, y, info, upd_flag_)
  end if

  call psb_erractionrestore(err_act)
  return
9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_cspan1D

subroutine psb_daxpby_mv_cspan2D(x, coeff, y, desc_a, info, upd_flag)
  use psb_base_mod, psb_protect_name => psb_daxpby_mv_cspan2D
  implicit none
  type(psb_d_multivect_type), intent(inout) :: x
  real(psb_dpk_), intent(in)                :: coeff(:, :)
  type(psb_d_multivect_type), intent(inout) :: y
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)            :: info
  logical, intent(in), optional             :: upd_flag

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err
  logical             :: upd_flag_ 
  
  upd_flag_ = .false.
  if(present(upd_flag)) upd_flag_ = upd_flag

  name = 'psb_daxpby_mv_cspan2D'
  if(psb_errstatus_fatal()) return
  
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v))) then
    info = psb_err_invalid_mvect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) & 
    & call x%axpby(desc_a%get_local_rows(), coeff, y, info, upd_flag_)

  call psb_erractionrestore(err_act)
  return
9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_mv_cspan2D


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
!         software without specific prior written permission.
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
! File: psb_daxpby.f90

!
! Subroutine: psb_daxpby_vect_out
!    Adds one distributed vector to another, 
!
!    Z := beta * Y + alpha * X
!
! Arguments:
!    alpha  - real, input            The scalar used to multiply each component of X
!    x      - type(psb_d_vect_type) The input vector containing the entries of X
!    beta   - real, input            The scalar used to multiply each component of Y
!    y      - type(psb_d_vect_type) The input vector Y
!    z      - type(psb_d_vect_type) The output vector Z
!    desc_a - type(psb_desc_type)   The communication descriptor.
!    info   - integer               Return code
!
!  Note: from a functional point of view, X is input, but here
!        it's declared INOUT because of the sync() methods.
!
subroutine psb_daxpby_vect_out(alpha, x, beta, y, z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpby_vect_out
  implicit none
  real(psb_dpk_), intent(in)            :: alpha, beta
  type(psb_d_vect_type), intent(inout)  :: x, y, z
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy, iiz, jjz
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, iz, ijz, m
  character(len=20)   :: name, ch_err

  name = 'psb_dgeaxpby'
  if(psb_errstatus_fatal()) return
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione
  iz = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, y%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iz, lone, desc_a, info, iiz, jjz)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 3'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione) .or. (iiz /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call z%axpby(desc_a%get_local_rows(), alpha, x, beta, y, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby_vect_out

!
! Subroutine: psb_daxpby
!    Adds one distributed matrix to another, 
!
!    sub( Y ) := beta * sub( Y ) + alpha * sub( X )
!
!    where sub( X ) denotes X(:, JX)
!
!    sub( Y ) denotes Y(:, JY).
!
! Arguments:
!    alpha   - real, input          The scalar used to multiply each component of X
!    x(:, :) - real, input          The input vector containing the entries of X
!    beta    - real, input          The scalar used to multiply each component of Y
!    y(:, :) - real, inout          The input vector Y
!    desc_a  - type(psb_desc_type)  The communication descriptor.
!    info    - integer              Return code
!    jx      - integer, optional    The column offset for X
!    jy      - integer, optional    The column offset for Y
!
subroutine psb_daxpby(alpha, x, beta, y, desc_a, info, n, jx, jy)
  use psb_base_mod, psb_protect_name => psb_daxpby
  use psi_d_serial_mod
  implicit none
  real(psb_dpk_), intent(in)      :: alpha, beta
  real(psb_dpk_), intent(in)      :: x(:, :)
  real(psb_dpk_), intent(inout)   :: y(:, :)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info
  integer(psb_ipk_), intent(in), optional :: n, jx, jy

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, in, jjy, lldx, lldy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  if(present(jx)) then
    ijx = jx
  else
    ijx = ione
  endif

  iy = ione
  if(present(jy)) then
    ijy = jy
  else
    ijy = ione
  endif

  if(present(n)) then
    if(((ijx + n) <= size(x, 2)) .and. ((ijy + n) <= size(y, 2))) then
      in = n
    else
      in = min(size(x, 2), size(y, 2))
    end if
  else
    in = min(size(x, 2), size(y, 2))
  endif

  if(ijx /= ijy) then
    info = 3050
    call psb_errpush(info, name)
    goto 9999
  end if

  m = desc_a%get_global_rows()
  lldx = size(x, 1)
  lldy = size(y, 1)

  ! check vector correctness
  call psb_chkvect(m, lone, lldx, ix, ijx, desc_a, info, iix, jjx)
  if(info == psb_success_) call psb_chkvect(m, lone, lldy, iy, ijy, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
    goto 9999
  end if

  if((in /= 0)) then
    if(desc_a%get_local_rows() > 0) then
      call psi_daxpby(desc_a%get_local_cols(), in, alpha, x(iix:, jjx:), beta, y(iiy:, jjy:), info)
    end if
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpby


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
!!$       software without specific prior written permission.
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
! Subroutine: psb_daxpbyv
!    Adds one distributed vector to another, 
!
!    Y := beta * Y + alpha * X
!
! Arguments:
!    alpha  - real, input         The scalar used to multiply each component of X
!    x(:)   - real, input         The input vector containing the entries of X
!    beta   - real, input         The scalar used to multiply each component of Y
!    y(:)   - real, inout         The input vector Y
!    desc_a - type(psb_desc_type) The communication descriptor.
!    info   - integer             Return code
!
subroutine psb_daxpbyv(alpha, x, beta, y, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpbyv
  implicit none
  real(psb_dpk_), intent(in)      :: alpha, beta
  real(psb_dpk_), intent(in)      :: x(:)
  real(psb_dpk_), intent(inout)   :: y(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy, lldx, lldy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err
  logical, parameter  :: debug = .false.

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_
    goto 9999
  end if

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()
  lldx = size(x, 1)
  lldy = size(y, 1)

  ! check vector correctness
  call psb_chkvect(m, lone, lldx, ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, lldy, iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call psb_geaxpby(desc_a%get_local_cols(), alpha, x, beta, y, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpbyv

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
!!$       software without specific prior written permission.
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
! Subroutine: psb_daxpbyvout
!    Adds one distributed vector to another, 
!
!    Z := beta * Y + alpha * X
!
! Arguments:
!    alpha  - real, input         The scalar used to multiply each component of X
!    x(:)   - real, input         The input vector containing the entries of X
!    beta   - real, input         The scalar used to multiply each component of Y
!    y(:)   - real, input         The input vector Y containing the entries of Y
!    Z(:)   - real, inout         The output vector Z
!    desc_a - type(psb_desc_type) The communication descriptor.
!    info   - integer             Return code
!
subroutine psb_daxpbyvout(alpha, x, beta, y, z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daxpbyvout
  implicit none
  real(psb_dpk_), intent(in)      :: alpha, beta
  real(psb_dpk_), intent(in)      :: x(:), y(:)
  real(psb_dpk_), intent(inout)   :: z(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)  :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, err_act, iix, jjx, iiy, jjy, iiz, jjz, &
                      & lldx, lldy, lldz
  integer(psb_lpk_) :: ix, ijx, iy, ijy, iz, ijz, m
  character(len=20) :: name, ch_err
  logical, parameter :: debug = .false.

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_
    goto 9999
  end if

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione
  iz = ione

  m = desc_a%get_global_rows()
  lldx = size(x, 1)
  lldy = size(y, 1)
  lldz = size(z, 1)

  ! check vector correctness
  call psb_chkvect(m, lone, lldx, ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, lldy, iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, lldz, iz, lone, desc_a, info, iiz, jjz)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if((iix /= ione) .or. (iiy /= ione) .or. (iiz /= ione)) then
    info = psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info, name)
  end if

  if(desc_a%get_local_rows() > 0) &
    & call psb_geaxpby(desc_a%get_local_cols(), alpha, x, beta, y, z, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daxpbyvout

!
! Subroutine: psb_daddconst_vect
!    Adds one distributed vector to another, 
!
!    Z(i) := X(i) + b
!
! Arguments:
!    x      - type(psb_d_vect_type) The input vector containing the entries of X
!    b      - real, input         The scalar used to add each component of X
!    z      - type(psb_d_vect_type) The input/output vector Z
!    desc_a - type(psb_desc_type)   The communication descriptor.
!    info   - integer               Return code
!
subroutine psb_daddconst_vect(x, b, z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_daddconst_vect
  implicit none
  type(psb_d_vect_type), intent(inout)  :: x
  type(psb_d_vect_type), intent(inout)  :: z
  real(psb_dpk_), intent(in)            :: b
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me, err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_d_addconst_vect'
  if(psb_errstatus_fatal()) return
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(z%v))) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%get_global_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if(desc_a%get_local_rows() > 0) call z%addconst(x, b, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_daddconst_vect

subroutine psb_d_upd_xyz_vect(alpha, beta, gamma, delta, x, y, z, desc_a, info)
  use psb_base_mod, psb_protect_name => psb_d_upd_xyz_vect
  implicit none 
  type(psb_d_vect_type), intent(inout)  :: x
  type(psb_d_vect_type), intent(inout)  :: y
  type(psb_d_vect_type), intent(inout)  :: z
  real(psb_dpk_), intent(in)            :: alpha, beta, gamma, delta
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)       :: info
  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np, me,  err_act, iix, jjx, iiy, jjy, nr
  integer(psb_lpk_)   :: ix, ijx, iy, ijy, m
  character(len=20)   :: name, ch_err

  name = 'psb_d_addconst_vect'
  if(psb_errstatus_fatal()) return
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if(np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info, name)
    goto 9999
  endif

  if((.not. allocated(x%v)) .or. (.not. allocated(y%v)) .or. (.not. allocated(z%v))) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m  = desc_a%get_global_rows()
  nr = desc_a%get_local_rows()

  ! check vector correctness
  call psb_chkvect(m, lone, x%get_nrows(), ix, lone, desc_a, info, iix, jjx)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 1'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  call psb_chkvect(m, lone, z%get_nrows(), iy, lone, desc_a, info, iiy, jjy)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err = 'psb_chkvect 2'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  if(desc_a%get_local_rows() > 0) &
    & call z%upd_xyz(nr, alpha, beta, gamma, delta, x, y, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_d_upd_xyz_vect
