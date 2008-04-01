!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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

  public psb_tolower, psb_toupper, psb_touppers
  interface psb_tolower
    module procedure psb_tolowerc
  end interface

  interface psb_toupper
    module procedure psb_toupperc
  end interface

  interface psb_touppers
    module procedure psb_sub_toupperc
  end interface

  private lcase, ucase, upper1c, lower1c
  character(len=*), parameter   :: lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter   :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

contains 

  function  psb_tolowerc(string)
    character(len=*), intent(in)  :: string
    character(len=len(string))    :: psb_tolowerc
    integer  :: i,k

    do i=1,len(string)
      psb_tolowerc(i:i) = lower1c(string(i:i))
!!$      k = index(ucase,string(i:i))
!!$      if (k /=0 ) then 
!!$        psb_tolowerc(i:i) = lcase(k:k)
!!$      else          
!!$        psb_tolowerc(i:i) = string(i:i)
!!$      end if
    enddo
  end function psb_tolowerc

  function  psb_toupperc(string)
    character(len=*), intent(in)  :: string
    character(len=len(string))    :: psb_toupperc
    integer  :: i,k

    do i=1,len(string)
      psb_toupperc(i:i) = upper1c(string(i:i))
!!$      k = index(lcase,string(i:i))
!!$      if (k /=0 ) then 
!!$        psb_toupperc(i:i) = ucase(k:k)
!!$      else          
!!$        psb_toupperc(i:i) = string(i:i)
!!$      end if
    enddo
  end function psb_toupperc

  subroutine   psb_sub_toupperc(string,strout)
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
  end subroutine psb_sub_toupperc


  function  lower1c(ch)
    character(len=1), intent(in) :: ch
    character(len=1)             :: lower1c

    select case(ch) 
    case ('A')
      lower1c = 'a'
    case ('B')
      lower1c = 'b'
    case ('C')
      lower1c = 'c'
    case ('D')
      lower1c = 'd'
    case ('E')
      lower1c = 'e'
    case ('F')
      lower1c = 'f'
    case ('G')
      lower1c = 'g'
    case ('H')
      lower1c = 'h'
    case ('I')
      lower1c = 'i'
    case ('J')
      lower1c = 'j'
    case ('K')
      lower1c = 'k'
    case ('L')
      lower1c = 'l'
    case ('M')
      lower1c = 'm'
    case ('N')
      lower1c = 'n'
    case ('O')
      lower1c = 'o'
    case ('P')
      lower1c = 'p'
    case ('Q')
      lower1c = 'q'
    case ('R')
      lower1c = 'r'
    case ('S')
      lower1c = 's'
    case ('T')
      lower1c = 't'
    case ('U')
      lower1c = 'u'
    case ('V')
      lower1c = 'v'
    case ('W')
      lower1c = 'w'
    case ('X')
      lower1c = 'x'
    case ('Y')
      lower1c = 'y'
    case ('Z')
      lower1c = 'z'
    case default
      lower1c = ch 
    end select
  end function lower1c

  function  upper1c(ch)
    character(len=1), intent(in) :: ch
    character(len=1)             :: upper1c

    select case(ch) 
    case ('a')
      upper1c = 'A'
    case ('b')
      upper1c = 'B'
    case ('c')
      upper1c = 'C'
    case ('d')
      upper1c = 'D'
    case ('e')
      upper1c = 'E'
    case ('f')
      upper1c = 'F'
    case ('g')
      upper1c = 'G'
    case ('h')
      upper1c = 'H'
    case ('i')
      upper1c = 'I'
    case ('j')
      upper1c = 'J'
    case ('k')
      upper1c = 'K'
    case ('l')
      upper1c = 'L'
    case ('m')
      upper1c = 'M'
    case ('n')
      upper1c = 'N'
    case ('o')
      upper1c = 'O'
    case ('p')
      upper1c = 'P'
    case ('q')
      upper1c = 'Q'
    case ('r')
      upper1c = 'R'
    case ('s')
      upper1c = 'S'
    case ('t')
      upper1c = 'T'
    case ('u')
      upper1c = 'U'
    case ('v')
      upper1c = 'V'
    case ('w')
      upper1c = 'W'
    case ('x')
      upper1c = 'X'
    case ('y')
      upper1c = 'Y'
    case ('z')
      upper1c = 'Z'
    case default
      upper1c = ch 
    end select
  end function upper1c


end module psb_string_mod
