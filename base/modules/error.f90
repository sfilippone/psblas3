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
!
!  Wrapper subroutines to provide error tools to F77 and C code
!

subroutine FCpsb_errcomm(ictxt, err)
  use psb_error_mod
  integer, intent(in)   :: ictxt
  integer, intent(inout):: err

  call psb_errcomm(ictxt, err)

end subroutine FCpsb_errcomm

subroutine FCpsb_errpush(err_c, r_name, i_err)
  use psb_error_mod
  implicit none
  
  integer, intent(in)              ::  err_c
  character(len=20), intent(in)    ::  r_name
  integer                          ::  i_err(5)

  call psb_errpush(err_c, r_name, i_err)
  
end subroutine FCpsb_errpush



subroutine FCpsb_serror()
  use psb_error_mod
  implicit none

  call psb_error()

end subroutine FCpsb_serror





subroutine FCpsb_perror(ictxt)
  use psb_error_mod
  implicit none

  integer, intent(in)   :: ictxt

  call psb_error(ictxt)

end subroutine FCpsb_perror





function FCpsb_get_errstatus()
  use psb_error_mod
  implicit none

  integer :: FCpsb_get_errstatus

  FCpsb_get_errstatus = psb_get_errstatus()

end function FCpsb_get_errstatus





subroutine FCpsb_get_errverbosity(v)
  use psb_error_mod
  implicit none

  integer, intent(out)   :: v

  v = psb_get_errverbosity()

end subroutine FCpsb_get_errverbosity




subroutine FCpsb_set_errverbosity(v)
  use psb_error_mod
  implicit none

  integer, intent(inout)   :: v

  call psb_set_errverbosity(v)

end subroutine FCpsb_set_errverbosity





subroutine FCpsb_erractionsave(err_act)
  use psb_error_mod
  implicit none

  integer, intent(out) :: err_act

  call psb_erractionsave(err_act)

end subroutine FCpsb_erractionsave


subroutine FCpsb_get_erraction(err_act)
  use psb_error_mod
  implicit none
  integer, intent(out) :: err_act 

  call psb_get_erraction(err_act)
end subroutine FCpsb_get_erraction



subroutine FCpsb_erractionrestore(err_act)
  use psb_error_mod
  implicit none

  integer, intent(in) :: err_act

  call psb_erractionrestore(err_act)

end subroutine FCpsb_erractionrestore






