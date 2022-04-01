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
! File: psb_zasb.f90
!
! Subroutine: psb_zasb_vect, _vect_r2 and _multivect
!    Assembles a dense vector, an array of vectors or a multivector
!    for PSBLAS routines. 
!    Since the allocation may have been called with the desciptor 
!    in the build state we make sure that X has a number of rows 
!    allowing for the halo indices, reallocating if necessary. 
!    We also call the halo routine for good measure.
!    However, sometimes we need to create temporary vectors whose contents
!    will be initialized through some subsequent call to geaxpby or similar;
!    for this situation we provide the SCRATCH flag. 
!    
! 
! Arguments: 
!    x       - type(psb_z_vect_type) The matrix to be assembled.
!    desc_a  - type(psb_desc_type).    The communication descriptor.
!    info    - integer.                return code
!    mold    - type(psb_z_base_vect_type), optional   A mold for the inner storage format
!    scratch - logical, optional       If true, allocate without checking/zeroing contents.
!                                      default: .false.
!    
subroutine psb_zasb_vect(x, desc_a, info, mold, scratch)
  use psb_base_mod, psb_protect_name => psb_zasb_vect
  implicit none

  type(psb_desc_type), intent(in)      ::  desc_a
  type(psb_z_vect_type), intent(inout) ::  x
  integer(psb_ipk_), intent(out)                 ::  info
  class(psb_z_base_vect_type), intent(in), optional :: mold
  logical, intent(in), optional        :: scratch

  ! local variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: i1sz,nrow,ncol, err_act, dupl_
  logical :: scratch_
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name,ch_err

  info = psb_success_
  name = 'psb_zgeasb_v'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt       = desc_a%get_context()
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  scratch_ = .false.
  if (present(scratch)) scratch_ = scratch
  call psb_info(ctxt, me, np)
  dupl_ = x%get_dupl()
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  else   if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': sizes: ',nrow,ncol

  if (scratch_) then 
    call x%free(info)
    call x%bld(ncol,mold=mold)
  else
    if (x%is_remote_build()) call psb_z_remote_vect(x,desc_a,info)
    call x%asb(ncol,info)
    ! ..update halo elements..
    call psb_halo(x,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_halo')
      goto 9999
    end if
    call x%cnv(mold)
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_zasb_vect


subroutine psb_zasb_vect_r2(x, desc_a, info, mold, scratch)
  use psb_base_mod, psb_protect_name => psb_zasb_vect_r2
  implicit none

  type(psb_desc_type), intent(in)      ::  desc_a
  type(psb_z_vect_type), intent(inout) ::  x(:)
  integer(psb_ipk_), intent(out)                 ::  info
  class(psb_z_base_vect_type), intent(in), optional :: mold
  logical, intent(in), optional        :: scratch

  ! local variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me, i, n 
  integer(psb_ipk_) :: i1sz,nrow,ncol, err_act, dupl_
  logical :: scratch_
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name,ch_err

  info = psb_success_
  name = 'psb_zgeasb_v'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt       = desc_a%get_context()
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  scratch_ = .false.
  if (present(scratch)) scratch_ = scratch
  call psb_info(ctxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  else   if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  n    = size(x)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': sizes: ',nrow,ncol

  if (scratch_) then 
    do i=1,n
      call x(i)%free(info)
      call x(i)%bld(ncol,mold=mold)
    end do

  else
    do i=1, n
      dupl_ = x(i)%get_dupl()
      call x(i)%asb(ncol,info)
      if (info /= 0) exit
      ! ..update halo elements..
      call psb_halo(x(i),desc_a,info)
      if (info /= 0) exit
      call x(i)%cnv(mold)
    end do
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_halo')
      goto 9999
    end if
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_zasb_vect_r2


subroutine psb_zasb_multivect(x, desc_a, info, mold, scratch,n)
  use psb_base_mod, psb_protect_name => psb_zasb_multivect
  implicit none

  type(psb_desc_type), intent(in)      ::  desc_a
  type(psb_z_multivect_type), intent(inout) ::  x
  integer(psb_ipk_), intent(out)                 ::  info
  class(psb_z_base_multivect_type), intent(in), optional :: mold
  integer(psb_ipk_), optional, intent(in)   :: n    
  logical, intent(in), optional        :: scratch

  ! local variables
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: i1sz,nrow,ncol, err_act, n_, dupl_
  logical :: scratch_
  
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name,ch_err

  info = psb_success_
  if (psb_errstatus_fatal()) return 

  name = 'psb_zgeasb'

  ctxt       = desc_a%get_context()
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  scratch_ = .false.
  if (present(scratch)) scratch_ = scratch

  if (present(n)) then 
    n_ = n
  else
    if (allocated(x%v)) then
      n_ = x%v%get_ncols()
    else
      n_ = 1
    end if
  endif

  call psb_info(ctxt, me, np)

  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  else   if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': sizes: ',nrow,ncol

  dupl_ = x%get_dupl()
  if (scratch_) then 
    call x%free(info)
    call x%bld(ncol,n_,mold=mold)
  else
    call x%asb(ncol,n_,info)
    ! ..update halo elements..
    call psb_halo(x,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_halo')
      goto 9999
    end if
    if (present(mold)) then 
      call x%cnv(mold)
    end if
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_zasb_multivect

