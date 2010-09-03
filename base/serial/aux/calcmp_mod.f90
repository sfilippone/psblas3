!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
module calcmp_mod
  use psb_const_mod
  interface operator(<)
    module procedure callt
  end interface
  interface operator(<=)
    module procedure calle
  end interface
  interface operator(>)
    module procedure calgt
  end interface
  interface operator(>=)
    module procedure calge
  end interface

contains

  function callt(a,b)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: callt
    
    callt = (abs(real(a))<abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and.(abs(aimag(a))<abs(aimag(b))))
  end function callt
  function calle(a,b)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: calle
    
    calle = (abs(real(a))<abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and.(abs(aimag(a))<=abs(aimag(b))))
  end function calle

  function calgt(a,b)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: calgt
    
    calgt = (abs(real(a))>abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and.(abs(aimag(a))>abs(aimag(b))))
  end function calgt
  function calge(a,b)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: calge
    
    calge = (abs(real(a))>abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and.(abs(aimag(a))>=abs(aimag(b))))
  end function calge

end module calcmp_mod
  
