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
module psb_e_hsort_mod
  use psb_const_mod

  interface psb_hsort
    subroutine psb_ehsort(x,ix,dir,flag)
      import 
      integer(psb_epk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_epk_), optional, intent(inout) :: ix(:)
    end subroutine psb_ehsort
  end interface psb_hsort


  interface psi_insert_heap
    subroutine psi_e_insert_heap(key,last,heap,dir,info)
      import 
      implicit none 

      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction

      integer(psb_epk_), intent(in)     :: key
      integer(psb_epk_), intent(inout)  :: heap(:)
      integer(psb_epk_), intent(in)     :: dir
      integer(psb_epk_), intent(inout)  :: last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_e_insert_heap
  end interface psi_insert_heap

  interface  psi_idx_insert_heap
    subroutine psi_e_idx_insert_heap(key,index,last,heap,idxs,dir,info)
      import 
      implicit none 

      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction

      integer(psb_epk_), intent(in)     :: key
      integer(psb_epk_), intent(inout)  :: heap(:)
      integer(psb_epk_), intent(in)     :: index
      integer(psb_epk_), intent(in)     :: dir
      integer(psb_epk_), intent(inout)  :: idxs(:)
      integer(psb_epk_), intent(inout)  :: last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_e_idx_insert_heap
  end interface psi_idx_insert_heap


  interface  psi_heap_get_first
    subroutine psi_e_heap_get_first(key,last,heap,dir,info)
      import 
      implicit none 
      integer(psb_epk_), intent(inout)   :: key
      integer(psb_epk_), intent(inout) :: last
      integer(psb_epk_), intent(in)     :: dir
      integer(psb_epk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_e_heap_get_first
  end interface psi_heap_get_first

  interface psi_idx_heap_get_first
    subroutine psi_e_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import
      integer(psb_epk_), intent(inout)    :: key
      integer(psb_epk_), intent(out)    :: index
      integer(psb_epk_), intent(inout)    :: heap(:)
      integer(psb_epk_), intent(in)     :: dir
      integer(psb_epk_), intent(inout)  :: last
      integer(psb_epk_), intent(inout)  :: idxs(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_e_idx_heap_get_first
  end interface psi_idx_heap_get_first


end module psb_e_hsort_mod
