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
module psi_lcx_mod
  use psb_const_mod
  interface operator(<)
    module procedure psi_cllt, psi_zllt
  end interface
  interface operator(<=)
    module procedure psi_clle, psi_zlle
  end interface
  interface operator(>)
    module procedure psi_clgt, psi_zlgt
  end interface
  interface operator(>=)
    module procedure psi_clge, psi_zlge
  end interface

contains

  function psi_cllt(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (real(a)<real(b)).or. &
         & ((real(a) == real(b)).and.(aimag(a)<aimag(b)))
  end function psi_cllt

  function psi_clle(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (real(a)<real(b)).or. &
         & ((real(a) == real(b)).and.(aimag(a)<=aimag(b)))
  end function psi_clle

  function psi_clgt(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (real(a)>real(b)).or. &
         & ((real(a) == real(b)).and.(aimag(a)>aimag(b)))
  end function psi_clgt

  function psi_clge(a,b) result(res)
    use psb_const_mod
    complex(psb_spk_), intent(in) :: a,b
    logical :: res
    
    res = (real(a)>real(b)).or. &
         & ((real(a) == real(b)).and.(aimag(a)>=aimag(b)))
  end function psi_clge

  function psi_zllt(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (real(a)<real(b)).or. &
         & ((real(a) == real(b)).and.(aimag(a)<aimag(b)))
  end function psi_zllt

  function psi_zlle(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (real(a)<real(b)).or. &
         & ((real(a) == real(b)).and.(aimag(a)<=aimag(b)))
  end function psi_zlle

  function psi_zlgt(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (real(a)>real(b)).or. &
         & ((real(a) == real(b)).and.(aimag(a)>aimag(b)))
  end function psi_zlgt

  function psi_zlge(a,b) result(res)
    use psb_const_mod
    complex(psb_dpk_), intent(in) :: a,b
    logical :: res
    
    res = (real(a)>real(b)).or. &
         & ((real(a) == real(b)).and.(aimag(a)>=aimag(b)))
  end function psi_zlge

end module psi_lcx_mod

