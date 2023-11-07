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
module psi_acx_mod
  use psb_const_mod
  interface operator(<)
    module procedure psi_calt, psi_zalt
  end interface
  interface operator(<=)
    module procedure psi_cale, psi_zale
  end interface
  interface operator(>)
    module procedure psi_cagt, psi_zagt
  end interface
  interface operator(>=)
    module procedure psi_cage, psi_zage
  end interface

contains

  function psi_calt(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res

    res = (abs(a) < abs(b))
  end function psi_calt

  function psi_cale(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(a) <= abs(b))
  end function psi_cale

  function psi_cagt(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(a) > abs(b))
  end function psi_cagt

  function psi_cage(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(a) >= abs(b))
  end function psi_cage

  function psi_zalt(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res

    res = (abs(a) < abs(b))
  end function psi_zalt

  function psi_zale(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(a) <= abs(b))
  end function psi_zale

  function psi_zagt(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(a) > abs(b))
  end function psi_zagt

  function psi_zage(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (abs(a) >= abs(b))
  end function psi_zage

end module psi_acx_mod
