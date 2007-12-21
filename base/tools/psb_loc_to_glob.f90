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
! File: psb_loc_to_glob.f90
!
! Subroutine: psb_loc_to_glob2
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
subroutine psb_loc_to_glob2(x,y,desc_a,info,iact)

  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  use psb_string_mod
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in)    ::  desc_a
  integer, intent(in)                ::  x(:)  
  integer, intent(out)               ::  y(:)  
  integer, intent(out)               ::  info
  character, intent(in), optional    ::  iact

  !....locals....
  integer                            ::  n, i, tmp
  character                          ::  act
  integer                            ::  int_err(5), err_act
  integer, parameter                 ::  zero=0
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=0
  name='psb_loc_to_glob2'
  call psb_erractionsave(err_act)

  if (present(iact)) then
    act=iact
  else
    act='I'
  endif
  act=toupper(act)

  n=size(x)
  do i=1,n
    if ((x(i).gt.psb_cd_get_local_cols(desc_a)).or.&
         &  (x(i).le.zero)) then
      info=140
      int_err(1)=tmp
      int_err(2)=psb_cd_get_local_cols(desc_a)  
      exit
    else
      tmp=desc_a%loc_to_glob(x(i))
      if((tmp.gt.zero).or.(tmp.le.psb_cd_get_global_rows(desc_a))) then
        y(i)=tmp
      else
        info = 140
        int_err(1)=tmp
        int_err(2)=psb_cd_get_local_cols(desc_a)
        exit
      end if
    end if
  enddo

  if (info /= 0) then
    select case(act)
    case('E','I')
      ! do nothing
    case('W')
      write(0,'("Error ",i5," in subroutine glob_to_loc")') info
    case('A')
      call psb_errpush(info,name)
      goto 9999
    end select
  endif

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

end subroutine psb_loc_to_glob2


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
! Subroutine: psb_loc_to_glob
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
subroutine psb_loc_to_glob(x,desc_a,info,iact)

  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  use psb_string_mod
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in)    ::  desc_a
  integer, intent(inout)             ::  x(:)  
  integer, intent(out)               ::  info
  character, intent(in), optional    ::  iact

  !....locals....
  integer                            ::  n ,i, tmp, err_act
  character                          ::  act
  integer                            ::  int_err(5)
  integer, parameter                 ::  zero=0
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=0
  name='psb_loc_to_glob'
  call psb_erractionsave(err_act)

  if (present(iact)) then
    act=iact
  else
    act='I'
  endif
  act = toupper(act)

  n=size(x)
  do i=1,n
    if ((x(i).gt.psb_cd_get_local_cols(desc_a)).or.&
         &  (x(i).le.zero)) then
      info=140
      int_err(1)=x(i)
      int_err(2)=psb_cd_get_local_cols(desc_a)  
      exit
    else
      tmp=desc_a%loc_to_glob(x(i))
      if((tmp.gt.zero).or.(tmp.le.psb_cd_get_global_rows(desc_a))) then
        x(i)=tmp
      else
        info = 140
        exit
      end if
    end if
  enddo

  if (info /= 0) then
    select case(act)
    case('E','I')
!!$      call psb_erractionrestore(err_act)
!!$      return
    case('W')
      write(0,'("Error ",i5," in subroutine glob_to_loc")') info
    case('A')
      call psb_errpush(info,name)
      goto 9999
    end select
  endif

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

end subroutine psb_loc_to_glob

