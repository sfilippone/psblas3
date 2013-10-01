!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
! File: psb_glob_to_loc.f90
!
! Subroutine: psb_glob_to_loc2v
!    Performs global to local index translation. If an index does not belong
!    to the current process, a negative value is returned (see also iact).
! 
! Arguments: 
!    x(:)     - integer                   Array containing the indices to be translated.
!    y(:)     - integer                   Array containing the translated indices.
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
!    iact     - character, optional       A character defining the behaviour on 
!                                         an index not belonging to the calling process
!                                         'I'gnore, 'W'arning, 'A'bort
!    owned    - logical, optional         When .true. limits the input to indices strictly
!                                         owned by the process, i.e. excludes halo.
!
subroutine psb_glob_to_loc2v(x,y,desc_a,info,iact,owned)
  use psb_base_mod, psb_protect_name => psb_glob_to_loc2v
  use psi_mod
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in) ::  desc_a
  integer(psb_ipk_), intent(in)                ::  x(:)  
  integer(psb_ipk_), intent(out)               ::  y(:), info
  character, intent(in), optional    ::  iact
  logical, intent(in),  optional     :: owned

  !....locals....
  integer(psb_ipk_) :: n, ictxt, iam, np 
  character                          :: act
  integer(psb_ipk_) :: int_err(5), err_act
  logical                            :: owned_
  integer(psb_ipk_), parameter                 :: zero=0
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'glob_to_loc'
  call psb_erractionsave(err_act)
  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ictxt = desc_a%get_context()
  call psb_info(ictxt,iam,np)


  if (present(iact)) then
    act=iact
  else
    act='I'
  endif
  act = psb_toupper(act)
  if (present(iact)) then
    act=iact
  else
    act='I'
  endif

  act = psb_toupper(act)

  n = size(x)
  call desc_a%indxmap%g2l(x(1:n),y(1:n),info,owned=owned)    

  select case(act)
  case('E','I')
    call psb_erractionrestore(err_act)
    return
  case('W')
    if ((info /= psb_success_).or.(count(y(1:n)<0) >0)) then
      write(psb_err_unit,'("Error ",i5," in subroutine glob_to_loc")') info
    end if
  case('A')
    if ((info /= psb_success_).or.(count(y(1:n)<0) >0)) then
      call psb_errpush(info,name)
      goto 9999
    end if
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error()
  end if
  return


end subroutine psb_glob_to_loc2v


!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
! Subroutine: psb_glob_to_loc1v
!    Performs global to local index translation. If an index does not belong
!    to the current process, a negative value is returned (see also iact).
! 
! Arguments: 
!    x(:)     - integer                   Array containing the indices to be translated.
!                                         overwritten on output with the result.  
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
!    iact     - character, optional       A character defining the behaviour on 
!                                         an index not belonging to the calling process
!                                         'I'gnore, 'W'arning, 'A'bort
!    owned    - logical, optional         When .true. limits the input to indices strictly
!                                         owned by the process, i.e. excludes halo.
!
subroutine psb_glob_to_loc1v(x,desc_a,info,iact,owned)
  use psb_base_mod, psb_protect_name => psb_glob_to_loc1v
  use psi_mod
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(inout)           :: x(:)  
  integer(psb_ipk_), intent(out)             :: info
  logical, intent(in),  optional     :: owned
  character, intent(in), optional  :: iact

  !....locals....
  integer(psb_ipk_) :: n
  character                        :: act
  integer(psb_ipk_) :: err_act
  logical                          :: owned_
  integer(psb_ipk_), parameter               :: zero=0
  character(len=20)   :: name
  integer(psb_ipk_) :: ictxt, iam, np

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'glob_to_loc'

  call psb_erractionsave(err_act)
  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ictxt = desc_a%get_context()
  call psb_info(ictxt,iam,np)


  if (present(iact)) then
    act=iact
  else
    act='I'
  endif

  act = psb_toupper(act)

  call desc_a%indxmap%g2lip(x,info,owned=owned)    

  select case(act)
  case('E','I')
    call psb_erractionrestore(err_act)
    return
  case('W')
    if ((info /= psb_success_).or.(count(x(1:n)<0) >0)) then
      write(psb_err_unit,'("Error ",i5," in subroutine glob_to_loc")') info
    end if
  case('A')
    if ((info /= psb_success_).or.(count(x(1:n)<0) >0)) then
      write(psb_err_unit,*) count(x(1:n)<0)
      call psb_errpush(info,name)
      goto 9999
    end if
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error()
  end if
  return

end subroutine psb_glob_to_loc1v

subroutine psb_glob_to_loc2s(x,y,desc_a,info,iact,owned)
  use psb_base_mod, psb_protect_name => psb_glob_to_loc2s
  implicit none 
  type(psb_desc_type), intent(in)    ::  desc_a
  integer(psb_ipk_),intent(in)                 ::  x
  integer(psb_ipk_),intent(out)                ::  y  
  integer(psb_ipk_), intent(out)               ::  info
  character, intent(in), optional    ::  iact
  logical, intent(in),  optional     :: owned

  integer(psb_ipk_) :: iv1(1), iv2(1)

  iv1(1) = x
  call psb_glob_to_loc(iv1,iv2,desc_a,info,iact,owned)
  y      = iv2(1)
end subroutine psb_glob_to_loc2s

subroutine psb_glob_to_loc1s(x,desc_a,info,iact,owned)
  use psb_base_mod, psb_protect_name => psb_glob_to_loc1s
  implicit none 
  type(psb_desc_type), intent(in)    ::  desc_a
  integer(psb_ipk_),intent(inout)              ::  x  
  integer(psb_ipk_), intent(out)               ::  info
  character, intent(in), optional    ::  iact
  logical, intent(in),  optional     :: owned
  integer(psb_ipk_) :: iv1(1)

  iv1(1) = x
  call psb_glob_to_loc(iv1,desc_a,info,iact,owned)
  x      = iv1(1)

end subroutine psb_glob_to_loc1s

