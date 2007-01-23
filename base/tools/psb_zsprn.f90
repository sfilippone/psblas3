!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
! File: psb_zsprn.f90
!
! Subroutine: psb_zsprn
!    Reinit sparse matrix structure for psblas routines.
! 
! Parameters: 
!    a        - type(<psb_zspmat_type>).          The sparse matrix to be reinitiated.      
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code.
!
Subroutine psb_zsprn(a, desc_a,info,clear)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  Implicit None

  !....Parameters...
  Type(psb_desc_type), intent(in)      :: desc_a
  Type(psb_zspmat_type), intent(inout) :: a
  integer, intent(out)                 :: info
  logical, intent(in), optional        :: clear

  !locals
  Integer             :: ictxt, np,me,err,err_act
  logical, parameter  :: debug=.false.
  integer             :: int_err(5)
  character(len=20)   :: name
  logical             :: clear_

  info = 0
  err  = 0
  int_err(1)=0
  name = 'psb_zsprn'
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (debug) &
       &write(*,*) 'starting spalloc ',ictxt,np,me

  
  if (psb_is_bld_desc(desc_a)) then
    ! Should do nothing, we are called redundantly
    return
  endif
  if (.not.psb_is_asb_desc(desc_a)) then
    info=590     
    call psb_errpush(info,name)
    goto 9999
  endif
  if (present(clear)) then 
    clear_ = clear
  else
    clear_ = .true.
  end if

  call psb_sp_reinit(a,info,clear=clear_)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_zsprn
