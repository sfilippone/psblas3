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
!
! File: psb_sallc.f90
!
! Function: psb_salloc_vect
!    Allocates dense vector for PSBLAS routines. 
!    The descriptor may be in either the build or assembled state.
! 
! Arguments: 
!    x      - the vector to be allocated.
!    desc_a - the communication descriptor.
!    info   - Return code
subroutine psb_salloc_vect(x, desc_a,info, dupl, bldmode)
  use psb_base_mod, psb_protect_name => psb_salloc_vect
  use psi_mod
  implicit none

  !....parameters...
  type(psb_s_vect_type), intent(out)  :: x
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_),intent(out)             :: info
  integer(psb_ipk_), optional, intent(in) :: dupl, bldmode

  !locals
  integer(psb_ipk_) :: np,me,nr,i,err_act
  integer(psb_ipk_) :: dupl_, bldmode_, nrmt_
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info=psb_success_
  if (psb_errstatus_fatal()) return 
  name='psb_geall'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check m and n parameters....
  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ! As this is a rank-1 array, optional parameter N is actually ignored.

  !....allocate x .....
  if (psb_is_asb_desc(desc_a).or.psb_is_upd_desc(desc_a)) then
    nr = max(1,desc_a%get_local_cols())
  else if (psb_is_bld_desc(desc_a)) then
    nr = max(1,desc_a%get_local_rows())
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='Invalid desc_a')
    goto 9999
  endif

  allocate(psb_s_base_vect_type :: x%v, stat=info) 
  if (info == 0) call x%all(nr,info)
  if (psb_errstatus_fatal()) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/nr/),a_err='real(psb_spk_)')
    goto 9999
  endif
  call x%zero()

  if (present(bldmode)) then
    bldmode_ = bldmode
  else
    bldmode_ = psb_matbld_noremote_
  end if
  if (present(dupl)) then
    dupl_ = dupl 
  else
    dupl_ = psb_dupl_def_
  end if
  call x%set_dupl(dupl_)
  call x%set_remote_build(bldmode_)
  call x%set_nrmv(0)
  if (x%is_remote_build()) then
    nrmt_ = max(100,(desc_a%get_local_cols()-desc_a%get_local_rows()))
    call psb_ensure_size(nrmt_,x%rmtv,info)
    call psb_ensure_size(nrmt_,x%rmidx,info)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_salloc_vect

! Function: psb_salloc_vect_r2
!    Allocates a vector of dense vectors for PSBLAS routines. 
!    The descriptor may be in either the build or assembled state.
! 
! Arguments: 
!    x      - the vector to be allocated.
!    desc_a - the communication descriptor.
!    info   - Return code
!    n      - optional number of columns.
!    lb     - optional lower bound on column indices

subroutine psb_salloc_vect_r2(x, desc_a,info,n,lb, dupl, bldmode)
  use psb_base_mod, psb_protect_name => psb_salloc_vect_r2
  use psi_mod
  implicit none

  !....parameters...
  type(psb_s_vect_type), allocatable, intent(out)  :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_),intent(out)             :: info
  integer(psb_ipk_), optional, intent(in)   :: n,lb
  integer(psb_ipk_), optional, intent(in) :: dupl, bldmode

  !locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me,nr,i,err_act, n_, lb_
  integer(psb_ipk_) :: dupl_, bldmode_, nrmt_
  integer(psb_ipk_) :: exch(1)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)  :: name

  info=psb_success_
  if (psb_errstatus_fatal()) return 
  name='psb_geall'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check m and n parameters....
  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(n)) then 
    n_ = n
  else
    n_ = 1
  endif
  if (present(lb)) then 
    lb_ = lb
  else
    lb_ = 1
  endif

  !global check on n parameters
  if (me == psb_root_) then
    exch(1)=n_
    call psb_bcast(ctxt,exch(1),root=psb_root_)
  else
    call psb_bcast(ctxt,exch(1),root=psb_root_)
    if (exch(1) /= n_) then
      info=psb_err_parm_differs_among_procs_
      call psb_errpush(info,name,i_err=(/ione/))
      goto 9999
    endif
  endif
  ! As this is a rank-1 array, optional parameter N is actually ignored.

  !....allocate x .....
  if (desc_a%is_asb().or.desc_a%is_upd()) then
    nr = max(1,desc_a%get_local_cols())
  else if (desc_a%is_bld()) then
    nr = max(1,desc_a%get_local_rows())
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='Invalid desc_a')
    goto 9999
  endif

  allocate(x(lb_:lb_+n_-1), stat=info)
  if (info == 0) then 
    do i=lb_, lb_+n_-1
      allocate(psb_s_base_vect_type :: x(i)%v, stat=info) 
      if (info == 0) call x(i)%all(nr,info)
      if (info == 0) call x(i)%zero()
      if (info /= 0) exit
    end do
  end if

  if (present(bldmode)) then
    bldmode_ = bldmode
  else
    bldmode_ = psb_matbld_noremote_
  end if
  if (present(dupl)) then
    dupl_ = dupl 
  else
    dupl_ = psb_dupl_def_
  end if
  
  do i=lb_, lb_+n_-1
    call x(i)%set_dupl(dupl_)
    call x(i)%set_remote_build(bldmode_)
    if (x(i)%is_remote_build()) then
      nrmt_ = max(100,(desc_a%get_local_cols()-desc_a%get_local_rows()))
      allocate(x(i)%rmtv(nrmt_))
    end if
  end do
  if (psb_errstatus_fatal()) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/nr/),a_err='real(psb_spk_)')
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_salloc_vect_r2


subroutine psb_salloc_multivect(x, desc_a,info,n, dupl, bldmode)
  use psb_base_mod, psb_protect_name => psb_salloc_multivect
  use psi_mod
  implicit none

  !....parameters...
  type(psb_s_multivect_type), intent(out)  :: x
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_),intent(out)             :: info
  integer(psb_ipk_), optional, intent(in)   :: n
  integer(psb_ipk_), optional, intent(in) :: dupl, bldmode

  !locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me,nr,i,err_act, n_, lb_
  integer(psb_ipk_) :: dupl_, bldmode_, nrmt_
  integer(psb_ipk_) :: exch(1)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)  :: name

  info=psb_success_
  if (psb_errstatus_fatal()) return 
  name='psb_geall'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check m and n parameters....
  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(n)) then 
    n_ = n
  else
    n_ = 1
  endif

  !global check on n parameters
  if (me == psb_root_) then
    exch(1)=n_
    call psb_bcast(ctxt,exch(1),root=psb_root_)
  else
    call psb_bcast(ctxt,exch(1),root=psb_root_)
    if (exch(1) /= n_) then
      info=psb_err_parm_differs_among_procs_
      call psb_errpush(info,name,i_err=(/ione/))
      goto 9999
    endif
  endif
  ! As this is a rank-1 array, optional parameter N is actually ignored.

  !....allocate x .....
  if (desc_a%is_asb().or.desc_a%is_upd()) then
    nr = max(1,desc_a%get_local_cols())
  else if (desc_a%is_bld()) then
    nr = max(1,desc_a%get_local_rows())
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='Invalid desc_a')
    goto 9999
  endif

  allocate(psb_s_base_multivect_type :: x%v, stat=info) 
  if (info == 0) call x%all(nr,n_,info)
  if (info == 0) call x%zero()
  
  if (psb_errstatus_fatal()) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/nr/),a_err='real(psb_spk_)')
    goto 9999
  endif
  if (present(bldmode)) then
    bldmode_ = bldmode
  else
    bldmode_ = psb_matbld_noremote_
  end if
  if (present(dupl)) then
    dupl_ = dupl 
  else
    dupl_ = psb_dupl_def_
  end if
  call x%set_dupl(dupl_)
  call x%set_remote_build(bldmode_)
  if (x%is_remote_build()) then
    nrmt_ = max(100,(desc_a%get_local_cols()-desc_a%get_local_rows()))
    allocate(x%rmtv(nrmt_,n_))
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_salloc_multivect
