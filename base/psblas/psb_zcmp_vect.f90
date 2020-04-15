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
! File: psb_zcmp_vect

subroutine psb_zcmp_vect(x,c,z,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_zcmp_vect
  implicit none
  type(psb_z_vect_type), intent (inout) :: x
  type(psb_z_vect_type), intent (inout) :: z
  real(psb_dpk_), intent(in)             :: c
  type(psb_desc_type), intent (in)        :: desc_a
  integer(psb_ipk_), intent(out)          :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, iiy, jjy
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
  character(len=20)        :: name, ch_err

  name='psb_z_cmp_vect'
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
  if (.not.allocated(z%v)) then
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
  call psb_chkvect(m,lone,z%get_nrows(),iy,lone,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect 2'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(desc_a%get_local_rows() > 0) then
    call z%acmp(x,c,info)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_zcmp_vect

subroutine psb_zcmp_spmatval(a,val,tol,desc_a,res,info)
  use psb_base_mod, psb_protect_name => psb_zcmp_spmatval
  implicit none
  type(psb_zspmat_type), intent(inout)  :: a
  complex(psb_dpk_), intent(in)             :: val
  real(psb_dpk_), intent(in)            :: tol
  type(psb_desc_type), intent (in)        :: desc_a
  integer(psb_ipk_), intent(out)          :: info
  logical, intent(out)                    :: res

  ! Local
  integer(psb_ipk_) :: ictxt, np, me
  integer(psb_ipk_) :: err_act
  character(len=20) :: name, ch_err
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_zcmp_spmatval'
  info=psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt=desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.((desc_a%get_local_rows() == a%get_nrows())&
    .and.(desc_a%get_local_cols() == a%get_ncols()))) then
      res = .false.
  else
    res = a%spcmp(val,tol,info)
  end if

  call psb_lallreduceand(ictxt,res)

  call psb_erractionrestore(err_act)
  if (debug_level >= psb_debug_comp_) then
    call psb_barrier(ictxt)
    write(debug_unit,*) me,' ',trim(name),' Returning '
  endif
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_zcmp_spmatval

subroutine psb_zcmp_spmat(a,b,tol,desc_a,res,info)
  use psb_base_mod, psb_protect_name => psb_zcmp_spmat
  implicit none
  type(psb_zspmat_type), intent(inout)  :: a
  type(psb_zspmat_type), intent(inout)  :: b
  real(psb_dpk_), intent(in)            :: tol
  type(psb_desc_type), intent (in)        :: desc_a
  integer(psb_ipk_), intent(out)          :: info
  logical, intent(out)                    :: res

  ! Local
  integer(psb_ipk_) :: ictxt, np, me
  integer(psb_ipk_) :: err_act
  character(len=20) :: name, ch_err
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_zcmp_spmatval'
  info=psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt=desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.((desc_a%get_local_rows() == a%get_nrows())&
    .and.(desc_a%get_local_rows() == b%get_nrows())&
    .and.(desc_a%get_local_cols() == a%get_ncols())&
      .and.(desc_a%get_local_cols() == b%get_ncols()))) then
      res = .false.
  else
    res = a%spcmp(b,tol,info)
  end if

  call psb_lallreduceand(ictxt,res)


  call psb_erractionrestore(err_act)
  if (debug_level >= psb_debug_comp_) then
    call psb_barrier(ictxt)
    write(debug_unit,*) me,' ',trim(name),' Returning '
  endif
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_zcmp_spmat
