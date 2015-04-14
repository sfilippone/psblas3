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
!
!  The merge-sort and quicksort routines are implemented in the
!  serial/aux directory
!  References:
!  D. Knuth
!  The Art of Computer Programming, vol. 3
!  Addison-Wesley
!  
!  Aho, Hopcroft, Ullman
!  Data Structures and Algorithms
!  Addison-Wesley
!
module psb_sort_mod
  use psb_const_mod

  ! 
  !  The up/down constant are defined in pairs having 
  !  opposite values. We make use of this fact in the heapsort routine.
  !
  integer(psb_ipk_), parameter :: psb_sort_up_      = 1, psb_sort_down_     = -1
  integer(psb_ipk_), parameter :: psb_lsort_up_     = 2, psb_lsort_down_    = -2
  integer(psb_ipk_), parameter :: psb_asort_up_     = 3, psb_asort_down_    = -3
  integer(psb_ipk_), parameter :: psb_alsort_up_    = 4, psb_alsort_down_   = -4
  integer(psb_ipk_), parameter :: psb_sort_ovw_idx_ = 0, psb_sort_keep_idx_ =  1
  integer(psb_ipk_), parameter :: psb_heap_resize   = 200

  type psb_int_heap
    integer(psb_ipk_) :: last, dir
    integer(psb_ipk_), allocatable :: keys(:)
  end type psb_int_heap
  type psb_int_idx_heap
    integer(psb_ipk_) :: last, dir
    integer(psb_ipk_), allocatable :: keys(:)
    integer(psb_ipk_), allocatable :: idxs(:)
  end type psb_int_idx_heap
  type psb_sreal_idx_heap
    integer(psb_ipk_) :: last, dir
    real(psb_spk_), allocatable    :: keys(:)
    integer(psb_ipk_), allocatable :: idxs(:)
  end type psb_sreal_idx_heap
  type psb_dreal_idx_heap
    integer(psb_ipk_) :: last, dir
    real(psb_dpk_), allocatable    :: keys(:)
    integer(psb_ipk_), allocatable :: idxs(:)
  end type psb_dreal_idx_heap
  type psb_scomplex_idx_heap
    integer(psb_ipk_) :: last, dir
    complex(psb_spk_), allocatable :: keys(:)
    integer(psb_ipk_), allocatable :: idxs(:)
  end type psb_scomplex_idx_heap
  type psb_dcomplex_idx_heap
    integer(psb_ipk_) :: last, dir
    complex(psb_dpk_), allocatable :: keys(:)
    integer(psb_ipk_), allocatable :: idxs(:)
  end type psb_dcomplex_idx_heap


  interface psb_iblsrch
    function  psb_iblsrch(key,n,v) result(ipos)
      import :: psb_ipk_
      integer(psb_ipk_) :: ipos, key, n
      integer(psb_ipk_) :: v(:)
    end function psb_iblsrch
  end interface

  interface psb_ibsrch
    function  psb_ibsrch(key,n,v) result(ipos)
      import :: psb_ipk_
      integer(psb_ipk_) :: ipos, key, n
      integer(psb_ipk_) :: v(:)
    end function psb_ibsrch
  end interface

  interface psb_issrch
    function psb_issrch(key,n,v) result(ipos)
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_) :: ipos, key, n
      integer(psb_ipk_) :: v(:)
    end function psb_issrch
  end interface

  interface psb_isaperm
    logical function psb_isaperm(n,eip)               
      import :: psb_ipk_
      integer(psb_ipk_), intent(in) :: n                             
      integer(psb_ipk_), intent(in) :: eip(n)
    end function psb_isaperm
  end interface


  interface psb_msort
    subroutine imsort(x,ix,dir,flag)
      import :: psb_ipk_
      integer(psb_ipk_), intent(inout)           :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine imsort
    subroutine smsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine smsort
    subroutine dmsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine dmsort
    subroutine camsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      complex(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine camsort
    subroutine zamsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine zamsort
  end interface


  interface psb_msort_unique
    subroutine imsort_u(x,nout,dir)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(inout)           :: x(:) 
      integer(psb_ipk_), intent(out)             :: nout
      integer(psb_ipk_), optional, intent(in)    :: dir
    end subroutine imsort_u
  end interface

  interface psb_qsort
    subroutine iqsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(inout)           :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine iqsort
    subroutine sqsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine sqsort
    subroutine dqsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine dqsort
    subroutine cqsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      complex(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine cqsort
    subroutine zqsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      complex(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine zqsort
  end interface
  

  interface psb_hsort
!!$    module procedure ihsort, shsort, chsort, dhsort, zhsort
    subroutine ihsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(inout)           :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine ihsort
    subroutine shsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      real(psb_spk_), intent(inout)    :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine shsort
    subroutine dhsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      real(psb_dpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine dhsort
    subroutine chsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      complex(psb_spk_), intent(inout) :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine chsort
    subroutine zhsort(x,ix,dir,flag)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      complex(psb_dpk_), intent(inout) :: x(:) 
      integer(psb_ipk_), optional, intent(in)      :: dir, flag
      integer(psb_ipk_), optional, intent(inout)   :: ix(:)
    end subroutine zhsort
  end interface


  interface psb_howmany_heap
    function  psb_howmany_int_heap(heap)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_int_heap
      type(psb_int_heap), intent(in) :: heap
      integer(psb_ipk_) :: psb_howmany_int_heap
    end function psb_howmany_int_heap
    function  psb_howmany_sreal_idx_heap(heap)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_sreal_idx_heap
      type(psb_sreal_idx_heap), intent(in) :: heap
      integer(psb_ipk_) :: psb_howmany_sreal_idx_heap
    end function psb_howmany_sreal_idx_heap
    function  psb_howmany_dreal_idx_heap(heap)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_dreal_idx_heap
      type(psb_dreal_idx_heap), intent(in) :: heap
      integer(psb_ipk_) :: psb_howmany_dreal_idx_heap
    end function psb_howmany_dreal_idx_heap
    function  psb_howmany_int_idx_heap(heap)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_int_idx_heap
      type(psb_int_idx_heap), intent(in) :: heap
      integer(psb_ipk_) :: psb_howmany_int_idx_heap
    end function psb_howmany_int_idx_heap
    function  psb_howmany_scomplex_idx_heap(heap)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_scomplex_idx_heap
      type(psb_scomplex_idx_heap), intent(in) :: heap
      integer(psb_ipk_) :: psb_howmany_scomplex_idx_heap
    end function psb_howmany_scomplex_idx_heap
    function  psb_howmany_dcomplex_idx_heap(heap)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_dcomplex_idx_heap
      type(psb_dcomplex_idx_heap), intent(in) :: heap
      integer(psb_ipk_) :: psb_howmany_dcomplex_idx_heap
    end function psb_howmany_dcomplex_idx_heap
  end interface
 

  interface psb_init_heap
    subroutine psb_init_int_heap(heap,info,dir)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_int_heap
      type(psb_int_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: dir
    end subroutine psb_init_int_heap
    subroutine psb_init_sreal_idx_heap(heap,info,dir)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_sreal_idx_heap
      type(psb_sreal_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: dir
    end subroutine psb_init_sreal_idx_heap
    subroutine psb_init_int_idx_heap(heap,info,dir)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_int_idx_heap
      type(psb_int_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: dir
    end subroutine psb_init_int_idx_heap
    subroutine psb_init_scomplex_idx_heap(heap,info,dir)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_scomplex_idx_heap
      type(psb_scomplex_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: dir
    end subroutine psb_init_scomplex_idx_heap
    subroutine psb_init_dcomplex_idx_heap(heap,info,dir)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_dcomplex_idx_heap
      type(psb_dcomplex_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: dir
    end subroutine psb_init_dcomplex_idx_heap
    subroutine psb_init_dreal_idx_heap(heap,info,dir)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_dreal_idx_heap
      type(psb_dreal_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: dir
    end subroutine psb_init_dreal_idx_heap
  end interface


  interface psb_dump_heap
    subroutine psb_dump_int_heap(iout,heap,info)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_int_heap
      type(psb_int_heap), intent(in) :: heap
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in)            :: iout
    end subroutine psb_dump_int_heap
    subroutine psb_dump_sreal_idx_heap(iout,heap,info)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_sreal_idx_heap
      type(psb_sreal_idx_heap), intent(in) :: heap
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in)            :: iout
    end subroutine psb_dump_sreal_idx_heap
    subroutine psb_dump_dreal_idx_heap(iout,heap,info)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_dreal_idx_heap
      type(psb_dreal_idx_heap), intent(in) :: heap
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in)            :: iout
    end subroutine psb_dump_dreal_idx_heap
    subroutine psb_dump_int_idx_heap(iout,heap,info)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_int_idx_heap
      type(psb_int_idx_heap), intent(in) :: heap
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in)            :: iout
    end subroutine psb_dump_int_idx_heap
    subroutine psb_dump_scomplex_idx_heap(iout,heap,info)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_scomplex_idx_heap
      type(psb_scomplex_idx_heap), intent(in) :: heap
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in)            :: iout
    end subroutine psb_dump_scomplex_idx_heap
    subroutine psb_dump_dcomplex_idx_heap(iout,heap,info)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      import :: psb_dcomplex_idx_heap
      type(psb_dcomplex_idx_heap), intent(in) :: heap
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), intent(in)            :: iout
    end subroutine psb_dump_dcomplex_idx_heap
  end interface


  interface psb_insert_heap
    subroutine psb_insert_int_heap(key,heap,info)
      import :: psb_int_heap, psb_ipk_
      integer(psb_ipk_), intent(in)               :: key
      type(psb_int_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_insert_int_heap
    subroutine psb_insert_int_idx_heap(key,index,heap,info)
      import :: psb_dpk_, psb_int_idx_heap, psb_ipk_
      integer(psb_ipk_), intent(in)                   :: key
      integer(psb_ipk_), intent(in)                   :: index
      type(psb_int_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)                  :: info
    end subroutine psb_insert_int_idx_heap
    subroutine psb_insert_sreal_idx_heap(key,index,heap,info)
      import :: psb_spk_, psb_sreal_idx_heap, psb_ipk_
      real(psb_spk_), intent(in)      :: key
      integer(psb_ipk_), intent(in)               :: index
      type(psb_sreal_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_insert_sreal_idx_heap
    subroutine psb_insert_dreal_idx_heap(key,index,heap,info)
      import :: psb_dpk_, psb_dreal_idx_heap, psb_ipk_
      real(psb_dpk_), intent(in)      :: key
      integer(psb_ipk_), intent(in)               :: index
      type(psb_dreal_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_insert_dreal_idx_heap
    subroutine psb_insert_scomplex_idx_heap(key,index,heap,info)
      import :: psb_spk_, psb_scomplex_idx_heap, psb_ipk_
      complex(psb_spk_), intent(in)              :: key
      integer(psb_ipk_), intent(in)                        :: index
      type(psb_scomplex_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_insert_scomplex_idx_heap
    subroutine psb_insert_dcomplex_idx_heap(key,index,heap,info)
      import :: psb_dpk_, psb_dcomplex_idx_heap, psb_ipk_
      complex(psb_dpk_), intent(in)            :: key
      integer(psb_ipk_), intent(in)                        :: index
      type(psb_dcomplex_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)                       :: info
    end subroutine psb_insert_dcomplex_idx_heap
  end interface

  interface psb_heap_get_first
    subroutine psb_int_heap_get_first(key,heap,info)
      import :: psb_int_heap, psb_ipk_
      type(psb_int_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)              :: key,info
    end subroutine psb_int_heap_get_first
    subroutine psb_int_idx_heap_get_first(key,index,heap,info)
      import :: psb_int_idx_heap, psb_ipk_
      type(psb_int_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)                  :: index,info
      integer(psb_ipk_), intent(out)                  :: key
    end subroutine psb_int_idx_heap_get_first
    subroutine psb_sreal_idx_heap_get_first(key,index,heap,info)
      import :: psb_spk_, psb_sreal_idx_heap, psb_ipk_
      type(psb_sreal_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)              :: index,info
      real(psb_spk_), intent(out)     :: key
    end subroutine psb_sreal_idx_heap_get_first
    subroutine psb_dreal_idx_heap_get_first(key,index,heap,info)
      import :: psb_dpk_, psb_dreal_idx_heap, psb_ipk_
      type(psb_dreal_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)              :: index,info
      real(psb_dpk_), intent(out)     :: key
    end subroutine psb_dreal_idx_heap_get_first
    subroutine psb_scomplex_idx_heap_get_first(key,index,heap,info)
      import :: psb_spk_, psb_scomplex_idx_heap, psb_ipk_
      type(psb_scomplex_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)                       :: index,info
      complex(psb_spk_), intent(out)           :: key
    end subroutine psb_scomplex_idx_heap_get_first
    
    subroutine psb_dcomplex_idx_heap_get_first(key,index,heap,info)
      import :: psb_dpk_, psb_dcomplex_idx_heap, psb_ipk_
      type(psb_dcomplex_idx_heap), intent(inout) :: heap
      integer(psb_ipk_), intent(out)                       :: index,info
      complex(psb_dpk_), intent(out)           :: key
    end subroutine psb_dcomplex_idx_heap_get_first
  end interface

  interface 
    subroutine psi_insert_int_heap(key,last,heap,dir,info)
      import :: psb_ipk_
      implicit none 
      
      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction
      
      integer(psb_ipk_), intent(in)     :: key,dir
      integer(psb_ipk_), intent(inout)  :: heap(:),last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_insert_int_heap
  end interface
  
  interface 
    subroutine psi_int_heap_get_first(key,last,heap,dir,info)
      import :: psb_ipk_
      implicit none 
     
      integer(psb_ipk_), intent(inout)  :: key,last
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_int_heap_get_first
  end interface
  
  interface 
    subroutine psi_insert_real_heap(key,last,heap,dir,info)
      import :: psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)    :: key
      integer(psb_ipk_), intent(in)           :: dir
      real(psb_spk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(inout)        :: last
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psi_insert_real_heap
  end interface
  
  interface 
    subroutine psi_real_heap_get_first(key,last,heap,dir,info)
      import :: psb_spk_, psb_ipk_
      real(psb_spk_), intent(inout) :: key
      integer(psb_ipk_), intent(inout)        :: last
      integer(psb_ipk_), intent(in)           :: dir
      real(psb_spk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)          :: info
    end subroutine psi_real_heap_get_first
  end interface
  
  interface 
    subroutine psi_insert_double_heap(key,last,heap,dir,info)
      import :: psb_dpk_, psb_ipk_
      real(psb_dpk_), intent(in)    :: key
      integer(psb_ipk_), intent(in)             :: dir
      real(psb_dpk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(inout)          :: last
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_insert_double_heap
  end interface
  
  interface 
    subroutine psi_double_heap_get_first(key,last,heap,dir,info)
      import :: psb_dpk_, psb_ipk_
      real(psb_dpk_), intent(inout) :: key
      integer(psb_ipk_), intent(inout)          :: last
      integer(psb_ipk_), intent(in)             :: dir
      real(psb_dpk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_double_heap_get_first
  end interface
  
  interface 
    subroutine psi_insert_scomplex_heap(key,last,heap,dir,info)
      import :: psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)    :: key
      integer(psb_ipk_), intent(in)              :: dir
      complex(psb_spk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(inout)           :: last
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_insert_scomplex_heap
  end interface
  
  interface 
    subroutine psi_scomplex_heap_get_first(key,last,heap,dir,info)
      import :: psb_spk_, psb_ipk_
      complex(psb_spk_), intent(inout) :: key
      integer(psb_ipk_), intent(inout)           :: last
      integer(psb_ipk_), intent(in)              :: dir
      complex(psb_spk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_scomplex_heap_get_first
  end interface
  
  interface 
    subroutine psi_insert_dcomplex_heap(key,last,heap,dir,info)
      import :: psb_dpk_, psb_ipk_
      complex(psb_dpk_), intent(in)    :: key
      integer(psb_ipk_), intent(in)                :: dir
      complex(psb_dpk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(inout)             :: last
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psi_insert_dcomplex_heap
  end interface
  
  interface 
    subroutine psi_dcomplex_heap_get_first(key,last,heap,dir,info)
      import :: psb_dpk_, psb_ipk_
      complex(psb_dpk_), intent(inout) :: key
      integer(psb_ipk_), intent(inout)             :: last
      integer(psb_ipk_), intent(in)                :: dir
      complex(psb_dpk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psi_dcomplex_heap_get_first
  end interface

  interface 
    subroutine psi_insert_int_idx_heap(key,index,last,heap,idxs,dir,info)
      import :: psb_ipk_
      integer(psb_ipk_), intent(in)     :: key
      integer(psb_ipk_), intent(in)     :: index,dir
      integer(psb_ipk_), intent(inout)  :: heap(:),last
      integer(psb_ipk_), intent(inout)  :: idxs(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_insert_int_idx_heap
  end interface
  
  interface 
    subroutine psi_int_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import :: psb_ipk_
      integer(psb_ipk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)   :: index,info
      integer(psb_ipk_), intent(inout) :: last,idxs(:)
      integer(psb_ipk_), intent(in)    :: dir
      integer(psb_ipk_), intent(out)   :: key
    end subroutine psi_int_idx_heap_get_first
  end interface
  
  interface 
    subroutine psi_insert_sreal_idx_heap(key,index,last,heap,idxs,dir,info)
      import :: psb_spk_, psb_ipk_
      real(psb_spk_), intent(in)     :: key
      integer(psb_ipk_), intent(in)            :: index,dir
      real(psb_spk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(inout)         :: idxs(:),last
      integer(psb_ipk_), intent(out)           :: info
    end subroutine psi_insert_sreal_idx_heap
  end interface
  
  interface 
    subroutine psi_sreal_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import :: psb_spk_, psb_ipk_
      real(psb_spk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)          :: index,info
      integer(psb_ipk_), intent(inout)        :: last,idxs(:)
      integer(psb_ipk_), intent(in)           :: dir
      real(psb_spk_), intent(out)   :: key
    end subroutine psi_sreal_idx_heap_get_first
  end interface

  interface 
    subroutine psi_insert_dreal_idx_heap(key,index,last,heap,idxs,dir,info)
      import :: psb_dpk_, psb_ipk_
      real(psb_dpk_), intent(in)     :: key
      integer(psb_ipk_), intent(in)              :: index,dir
      real(psb_dpk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(inout)           :: idxs(:),last
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_insert_dreal_idx_heap
  end interface
  
  interface 
    subroutine psi_dreal_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import :: psb_dpk_, psb_ipk_
      real(psb_dpk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)            :: index,info
      integer(psb_ipk_), intent(inout)          :: last,idxs(:)
      integer(psb_ipk_), intent(in)             :: dir
      real(psb_dpk_), intent(out)   :: key
    end subroutine psi_dreal_idx_heap_get_first
  end interface

  interface 
    subroutine psi_insert_scomplex_idx_heap(key,index,last,heap,idxs,dir,info)
      import :: psb_spk_, psb_ipk_
      complex(psb_spk_), intent(in)    :: key
      integer(psb_ipk_), intent(in)              :: index,dir
      complex(psb_spk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(inout)           :: idxs(:),last
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_insert_scomplex_idx_heap
  end interface

  interface 
    subroutine psi_scomplex_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import :: psb_spk_, psb_ipk_
      complex(psb_spk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)             :: index,info
      integer(psb_ipk_), intent(inout)           :: last,idxs(:)
      integer(psb_ipk_), intent(in)              :: dir
      complex(psb_spk_), intent(out)   :: key
    end subroutine psi_scomplex_idx_heap_get_first
  end interface

  interface 
    subroutine psi_insert_dcomplex_idx_heap(key,index,last,heap,idxs,dir,info)
      import :: psb_dpk_, psb_ipk_
      complex(psb_dpk_), intent(in)    :: key
      integer(psb_ipk_), intent(in)                :: index,dir
      complex(psb_dpk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(inout)             :: idxs(:),last
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psi_insert_dcomplex_idx_heap
  end interface

  interface 
    subroutine psi_dcomplex_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import :: psb_dpk_, psb_ipk_
      complex(psb_dpk_), intent(inout) :: heap(:)
      integer(psb_ipk_), intent(out)               :: index,info
      integer(psb_ipk_), intent(inout)             :: last,idxs(:)
      integer(psb_ipk_), intent(in)                :: dir
      complex(psb_dpk_), intent(out)   :: key
    end subroutine psi_dcomplex_idx_heap_get_first
  end interface
  

  interface psb_free_heap
    module procedure psb_free_int_heap, psb_free_int_idx_heap,&
         & psb_free_sreal_idx_heap, psb_free_scomplex_idx_heap, &
         & psb_free_dreal_idx_heap, psb_free_dcomplex_idx_heap
  end interface

contains

  subroutine psb_free_int_heap(heap,info)
    implicit none 
    type(psb_int_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)

  end subroutine psb_free_int_heap

  subroutine psb_free_sreal_idx_heap(heap,info)
    implicit none 
    type(psb_sreal_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)
    if ((info == psb_success_).and.(allocated(heap%idxs))) deallocate(heap%idxs,stat=info)

  end subroutine psb_free_sreal_idx_heap

  subroutine psb_free_dreal_idx_heap(heap,info)
    implicit none 
    type(psb_dreal_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)
    if ((info == psb_success_).and.(allocated(heap%idxs))) deallocate(heap%idxs,stat=info)

  end subroutine psb_free_dreal_idx_heap

  subroutine psb_free_int_idx_heap(heap,info)
    implicit none 
    type(psb_int_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)
    if ((info == psb_success_).and.(allocated(heap%idxs))) deallocate(heap%idxs,stat=info)

  end subroutine psb_free_int_idx_heap

  subroutine psb_free_scomplex_idx_heap(heap,info)
    implicit none 
    type(psb_scomplex_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)
    if ((info == psb_success_).and.(allocated(heap%idxs))) deallocate(heap%idxs,stat=info)

  end subroutine psb_free_scomplex_idx_heap

  subroutine psb_free_dcomplex_idx_heap(heap,info)
    implicit none 
    type(psb_dcomplex_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)
    if ((info == psb_success_).and.(allocated(heap%idxs))) deallocate(heap%idxs,stat=info)

  end subroutine psb_free_dcomplex_idx_heap

end module psb_sort_mod
