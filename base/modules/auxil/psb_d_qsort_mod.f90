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
module psb_d_qsort_mod
  use psb_const_mod

  interface psb_bsrch
    function  psb_dbsrch(key,n,v,dir,find) result(ipos)
      import 
      integer(psb_ipk_) :: ipos, n
      real(psb_dpk_) :: key
      real(psb_dpk_) :: v(:)
      integer(psb_ipk_), optional :: dir, find
    end function psb_dbsrch
  end interface psb_bsrch

  interface psb_ssrch
    function psb_dssrch(key,n,v) result(ipos)
      import 
      implicit none
      integer(psb_ipk_) :: ipos, n
      real(psb_dpk_) :: key
      real(psb_dpk_) :: v(:)
    end function psb_dssrch
  end interface psb_ssrch

  interface psb_qsort
    subroutine psb_dqsort(x,ix,dir,flag)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_dqsort
  end interface psb_qsort
  
  interface 
    subroutine psi_dqsrx_up(n,x,ix)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_dqsrx_up
    subroutine psi_dqsrx_dw(n,x,ix)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_dqsrx_dw
    subroutine psi_dqsr_up(n,x)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_dqsr_up
    subroutine psi_dqsr_dw(n,x)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_dqsr_dw
    subroutine psi_daqsrx_up(n,x,ix)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_daqsrx_up
    subroutine psi_daqsrx_dw(n,x,ix)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_daqsrx_dw
    subroutine psi_daqsr_up(n,x)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_daqsr_up
    subroutine psi_daqsr_dw(n,x)
      import 
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_daqsr_dw
  end interface

end module psb_d_qsort_mod
