!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: psb_loc_to_glob.f90
!
! Subroutine: psb_loc_to_glob2v
!    Performs local to global index translation. If an index is out of range
!    a negative value is returned (see also iact).
! 
! Arguments: 
!    x(:)     - integer                   Array containing the indices to be translated.
!    y(:)     - integer                   Array containing the translated indices.
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
!    iact     - character, optional       A character defining the behaviour on 
!                                         an out of range index 
!                                         'I'gnore, 'W'arning, 'A'bort
!
subroutine psb_loc_to_glob2v(x,y,desc_a,info,iact)
  use psb_base_mod, psb_protect_name => psb_loc_to_glob2v
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in)    ::  desc_a
  integer(psb_ipk_), intent(in)                ::  x(:)  
  integer(psb_ipk_), intent(out)               ::  y(:)  
  integer(psb_ipk_), intent(out)               ::  info
  character, intent(in), optional    ::  iact

  !....locals....
  integer(psb_ipk_) ::  n, i, tmp
  character                          ::  act
  integer(psb_ipk_) ::  int_err(5), err_act
  integer(psb_ipk_), parameter                 ::  zero=0
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='psb_loc_to_glob2'
  call psb_erractionsave(err_act)
  if (.not.desc_a%is_valid()) then 
    info =  psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  if (present(iact)) then
    act=iact
  else
    act='I'
  endif
  act=psb_toupper(act)

  call desc_a%indxmap%l2g(x,y,info) 

  if (info /= psb_success_) then
    select case(act)
    case('E','I')
      ! do nothing, silently.
      info = psb_success_
    case('W')
      write(psb_err_unit,'("Error ",i5," in subroutine loc_to_glob")') info
      info = psb_success_
    case('A')
      call psb_errpush(info,name)
      goto 9999
    end select
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_loc_to_glob2v


!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Subroutine: psb_loc_to_glob1v
!    Performs local to global index translation. If an index is out of range
!    a negative value is returned (see also iact).
! 
! Arguments: 
!    x(:)     - integer                   Array containing the indices to be translated.
!                                         Overwritten on output. 
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
!    iact     - character, optional       A character defining the behaviour on 
!                                         an out of range index 
!                                         'I'gnore, 'W'arning, 'A'bort
!
subroutine psb_loc_to_glob1v(x,desc_a,info,iact)
  use psb_base_mod, psb_protect_name => psb_loc_to_glob1v
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in)    ::  desc_a
  integer(psb_ipk_), intent(inout)             ::  x(:)  
  integer(psb_ipk_), intent(out)               ::  info
  character, intent(in), optional    ::  iact

  !....locals....
  integer(psb_ipk_) ::  n ,i, tmp, err_act
  character                          ::  act
  integer(psb_ipk_) ::  int_err(5)
  integer(psb_ipk_), parameter                 ::  zero=0
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='psb_loc_to_glob'
  call psb_erractionsave(err_act)
  if (.not.desc_a%is_valid()) then 
    info =  psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iact)) then
    act=iact
  else
    act='I'
  endif
  act = psb_toupper(act)

  call desc_a%indxmap%l2gip(x,info) 

  if (info /= psb_success_) then
    select case(act)
    case('E','I')
      ! do nothing, silently.
      info = psb_success_
    case('W')
      write(psb_err_unit,'("Error ",i5," in subroutine loc_to_glob")') info
      info = psb_success_
    case('A')
      call psb_errpush(info,name)
      goto 9999
    end select
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_loc_to_glob1v

subroutine psb_loc_to_glob2s(x,y,desc_a,info,iact)
  use psb_desc_mod
  use psb_tools_mod, psb_protect_name => psb_loc_to_glob2s
  implicit none 
  type(psb_desc_type), intent(in)    ::  desc_a
  integer(psb_ipk_),intent(in)                 ::  x
  integer(psb_ipk_),intent(out)                ::  y  
  integer(psb_ipk_), intent(out)               ::  info
  character, intent(in), optional    ::  iact

  integer(psb_ipk_) :: iv1(1), iv2(1)

  iv1(1) = x
  call psb_loc_to_glob(iv1,iv2,desc_a,info,iact)
  y      = iv2(1)
end subroutine psb_loc_to_glob2s

subroutine psb_loc_to_glob1s(x,desc_a,info,iact)
  use psb_desc_mod
  use psb_tools_mod, psb_protect_name => psb_loc_to_glob1s
  implicit none 
  type(psb_desc_type), intent(in)    ::  desc_a
  integer(psb_ipk_),intent(inout)              ::  x  
  integer(psb_ipk_), intent(out)               ::  info
  character, intent(in), optional    ::  iact
  integer(psb_ipk_) :: iv1(1)

  iv1(1) = x
  call psb_loc_to_glob(iv1,desc_a,info,iact)
  x      = iv1(1)

end subroutine psb_loc_to_glob1s

