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
module psb_string_mod

  public tolower, toupper, touppers
  interface tolower
    module procedure tolowerc
  end interface

  interface toupper
    module procedure toupperc
  end interface

  interface touppers
    module procedure sub_toupperc
  end interface

  private
  character(len=*), parameter   :: lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter   :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

contains 

  function  tolowerc(string)
    character(len=*), intent(in)  :: string
    character(len=len(string))    :: tolowerc
    integer  :: i,k

    do i=1,len(string)
      k = index(ucase,string(i:i))
      if (k /=0 ) then 
        tolowerc(i:i) = lcase(k:k)
      else          
        tolowerc(i:i) = string(i:i)
      end if
    enddo
  end function tolowerc

  function  toupperc(string)
    character(len=*), intent(in)  :: string
    character(len=len(string))    :: toupperc
    integer  :: i,k

    do i=1,len(string)
      k = index(lcase,string(i:i))
      if (k /=0 ) then 
        toupperc(i:i) = ucase(k:k)
      else          
        toupperc(i:i) = string(i:i)
      end if
    enddo
  end function toupperc

  subroutine   sub_toupperc(string,strout)
    character(len=*), intent(in)  :: string
    character(len=*), intent(out)  :: strout
    integer  :: i,k

    do i=1,len(string)
      k = index(lcase,string(i:i))
      if (k /=0 ) then 
        strout(i:i) = ucase(k:k)
      else          
        strout(i:i) = string(i:i)
      end if
    enddo
  end subroutine sub_toupperc

end module psb_string_mod
