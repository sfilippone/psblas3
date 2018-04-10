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
module psi_alcx_mod
  use psb_const_mod
  interface operator(<)
    module procedure psi_callt, psi_zallt
  end interface
  interface operator(<=)
    module procedure psi_calle, psi_zalle
  end interface
  interface operator(>)
    module procedure psi_calgt, psi_zalgt
  end interface
  interface operator(>=)
    module procedure psi_calge, psi_zalge
  end interface

contains

  function psi_callt(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(real(a))<abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and. &
         & (abs(aimag(a))<abs(aimag(b))))
  end function psi_callt

  function psi_calle(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(real(a))<abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and. &
         & (abs(aimag(a))<=abs(aimag(b))))
  end function psi_calle

  function psi_calgt(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(real(a))>abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and. &
         & (abs(aimag(a))>abs(aimag(b))))
  end function psi_calgt

  function psi_calge(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(real(a))>abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and. &
         & (abs(aimag(a))>=abs(aimag(b))))
  end function psi_calge

  function psi_zallt(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(real(a))<abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and. &
         & (abs(aimag(a))<abs(aimag(b))))
  end function psi_zallt

  function psi_zalle(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(real(a))<abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and. &
         & (abs(aimag(a))<=abs(aimag(b))))
  end function psi_zalle

  function psi_zalgt(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(real(a))>abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and. &
         & (abs(aimag(a))>abs(aimag(b))))
  end function psi_zalgt

  function psi_zalge(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(real(a))>abs(real(b))).or. &
         & ((abs(real(a)) == abs(real(b))).and. &
         & (abs(aimag(a))>=abs(aimag(b))))
  end function psi_zalge

end module psi_alcx_mod
  
