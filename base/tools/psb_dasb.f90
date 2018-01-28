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
! File: psb_dasb.f90
!
! Subroutine: psb_dasb
!    Assembles a dense matrix for PSBLAS routines
!    Since the allocation may have been called with the desciptor 
!    in the build state we make sure that X has a number of rows 
!    allowing for the halo indices, reallocating if necessary. 
!    We also call the halo routine for good measure.
! 
! Arguments: 
!    x(:,:)  - real, allocatable    The matrix to be assembled.
!    desc_a  - type(psb_desc_type).  The communication descriptor.
!    info    - integer.                return code
subroutine psb_dasb(x, desc_a, info, scratch)
  use psb_base_mod, psb_protect_name => psb_dasb
  implicit none

  type(psb_desc_type), intent(in) ::  desc_a
  real(psb_dpk_), allocatable, intent(inout) ::  x(:,:)
  integer(psb_ipk_), intent(out)            ::  info
  logical, intent(in), optional        :: scratch

  ! local variables
  integer(psb_ipk_) :: ictxt,np,me,nrow,ncol, err_act
  integer(psb_ipk_) :: i1sz, i2sz
  integer(psb_ipk_) :: debug_level, debug_unit
  logical :: scratch_
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='psb_dgeasb_m'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  scratch_ = .false.
  if (present(scratch)) scratch_ = scratch

  if (.not.psb_is_ok_desc(desc_a)) then
    info=psb_err_input_matrix_unassembled_
    call psb_errpush(info,name)
    goto 9999
  endif
  ictxt   = desc_a%get_context()

  call psb_info(ictxt, me, np)


  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': start: ',np,&
       & desc_a%get_dectype()
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  else if (.not.psb_is_asb_desc(desc_a)) then
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),' error '
    info = psb_err_input_matrix_unassembled_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! check size
  ictxt = desc_a%get_context()
  nrow  = desc_a%get_local_rows()
  ncol  = desc_a%get_local_cols()
  i1sz = size(x,dim=1)
  i2sz = size(x,dim=2)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': ',i1sz,i2sz,nrow,ncol

  if (i1sz < ncol) then
    call psb_realloc(ncol,i2sz,x,info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_realloc')
      goto 9999
    endif
  endif

  if (.not.scratch_) then 
    ! ..update halo elements..
    call psb_halo(x,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_halo'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_dasb


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
!!$       software without specific written permission.
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
! Subroutine: psb_dasb
!    Assembles a dense matrix for PSBLAS routines
!    Since the allocation may have been called with the desciptor 
!    in the build state we make sure that X has a number of rows 
!    allowing for the halo indices, reallocating if necessary. 
!    We also call the halo routine for good measure.
! 
! Arguments: 
!    x(:)    - real, allocatable    The matrix to be assembled.
!    desc_a  - type(psb_desc_type).  The communication descriptor.
!    info    - integer.                Return  code
subroutine psb_dasbv(x, desc_a, info, scratch)
  use psb_base_mod, psb_protect_name => psb_dasbv
  implicit none

  type(psb_desc_type), intent(in)                 ::  desc_a
  real(psb_dpk_), allocatable, intent(inout) ::  x(:)
  integer(psb_ipk_), intent(out)        ::  info
  logical, intent(in), optional        :: scratch

  ! local variables
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: int_err(5), i1sz,nrow,ncol, err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  logical :: scratch_
  character(len=20)    :: name,ch_err

  info = psb_success_
  int_err(1) = 0
  name = 'psb_dgeasb_v'

  ictxt   = desc_a%get_context()
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  scratch_ = .false.
  if (present(scratch)) scratch_ = scratch

  call psb_info(ictxt, me, np)

  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  else if (.not.psb_is_asb_desc(desc_a)) then
    info = psb_err_input_matrix_unassembled_
    call psb_errpush(info,name)
    goto 9999
  endif

  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': sizes: ',nrow,ncol
  i1sz = size(x)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': sizes ',i1sz,ncol
  if (i1sz < ncol) then
    call psb_realloc(ncol,x,info)
    if (info /= psb_success_) then           
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_realloc')
      goto 9999
    endif
  endif

  if (.not.scratch_) then 
    ! ..update halo elements..
    call psb_halo(x,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='f90_pshalo'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_dasbv


subroutine psb_dasb_vect(x, desc_a, info, mold, scratch)
  use psb_base_mod, psb_protect_name => psb_dasb_vect
  implicit none

  type(psb_desc_type), intent(in)      ::  desc_a
  type(psb_d_vect_type), intent(inout) ::  x
  integer(psb_ipk_), intent(out)                 ::  info
  class(psb_d_base_vect_type), intent(in), optional :: mold
  logical, intent(in), optional        :: scratch

  ! local variables
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: int_err(5), i1sz,nrow,ncol, err_act
  logical :: scratch_
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name,ch_err

  info = psb_success_
  if (psb_errstatus_fatal()) return 

  int_err(1) = 0
  name = 'psb_dgeasb_v'

  ictxt       = desc_a%get_context()
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  scratch_ = .false.
  if (present(scratch)) scratch_ = scratch
  call psb_info(ictxt, me, np)

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
    call x%asb(ncol,info)
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

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_dasb_vect


subroutine psb_dasb_vect_r2(x, desc_a, info, mold, scratch)
  use psb_base_mod, psb_protect_name => psb_dasb_vect_r2
  implicit none

  type(psb_desc_type), intent(in)      ::  desc_a
  type(psb_d_vect_type), intent(inout) ::  x(:)
  integer(psb_ipk_), intent(out)                 ::  info
  class(psb_d_base_vect_type), intent(in), optional :: mold
  logical, intent(in), optional        :: scratch

  ! local variables
  integer(psb_ipk_) :: ictxt,np,me, i, n 
  integer(psb_ipk_) :: int_err(5), i1sz,nrow,ncol, err_act
  logical :: scratch_
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name,ch_err

  info = psb_success_
  if (psb_errstatus_fatal()) return 

  int_err(1) = 0
  name = 'psb_dgeasb_v'

  ictxt       = desc_a%get_context()
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  scratch_ = .false.
  if (present(scratch)) scratch_ = scratch
  call psb_info(ictxt, me, np)

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
      call x(i)%asb(ncol,info)
      if (info /= 0) exit
      ! ..update halo elements..
      call psb_halo(x(i),desc_a,info)
      if (info /= 0) exit
      if (present(mold)) then 
        call x(i)%cnv(mold)
      end if
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

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_dasb_vect_r2


subroutine psb_dasb_multivect(x, desc_a, info, mold, scratch,n)
  use psb_base_mod, psb_protect_name => psb_dasb_multivect
  implicit none

  type(psb_desc_type), intent(in)      ::  desc_a
  type(psb_d_multivect_type), intent(inout) ::  x
  integer(psb_ipk_), intent(out)                 ::  info
  class(psb_d_base_multivect_type), intent(in), optional :: mold
  integer(psb_ipk_), optional, intent(in)   :: n    
  logical, intent(in), optional        :: scratch

  ! local variables
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: int_err(5), i1sz,nrow,ncol, err_act, n_
  logical :: scratch_
  
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name,ch_err

  info = psb_success_
  if (psb_errstatus_fatal()) return 

  int_err(1) = 0
  name = 'psb_dgeasb'

  ictxt       = desc_a%get_context()
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

  call psb_info(ictxt, me, np)

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

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_dasb_multivect

