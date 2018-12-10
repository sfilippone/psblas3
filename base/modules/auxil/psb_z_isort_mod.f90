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
!
!  Sorting routines
!  References:
!  D. Knuth
!  The Art of Computer Programming, vol. 3
!  Addison-Wesley
!  
!  Aho, Hopcroft, Ullman
!  Data Structures and Algorithms
!  Addison-Wesley
!
module psb_z_isort_mod
  use psb_const_mod

  interface psb_isort
    subroutine psb_zisort(x,ix,dir,flag)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_zisort
  end interface psb_isort



  interface 
    subroutine psi_zlisrx_up(n,x,ix)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zlisrx_up
    subroutine psi_zlisrx_dw(n,x,ix)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zlisrx_dw
    subroutine psi_zlisr_up(n,x)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zlisr_up
    subroutine psi_zlisr_dw(n,x)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zlisr_dw
    subroutine psi_zalisrx_up(n,x,ix)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zalisrx_up
    subroutine psi_zalisrx_dw(n,x,ix)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zalisrx_dw
    subroutine psi_zalisr_up(n,x)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zalisr_up
    subroutine psi_zalisr_dw(n,x)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zalisr_dw
    subroutine psi_zaisrx_up(n,x,ix)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zaisrx_up
    subroutine psi_zaisrx_dw(n,x,ix)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zaisrx_dw
    subroutine psi_zaisr_up(n,x)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zaisr_up
    subroutine psi_zaisr_dw(n,x)
      import 
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_zaisr_dw
  end interface


end module psb_z_isort_mod
