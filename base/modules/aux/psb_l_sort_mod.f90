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
module psb_l_sort_mod
  use psb_const_mod

  interface psb_isaperm
    logical function psb_lisaperm(n,eip)               
      import 
      integer(psb_lpk_), intent(in) :: n                             
      integer(psb_lpk_), intent(in) :: eip(n)
    end function psb_lisaperm
  end interface psb_isaperm

  interface psb_msort_unique
    subroutine psb_lmsort_u(x,nout,dir)
      import 
      integer(psb_lpk_), intent(inout)           :: x(:) 
      integer(psb_lpk_), intent(out)             :: nout
      integer(psb_ipk_), optional, intent(in)    :: dir
    end subroutine psb_lmsort_u
  end interface psb_msort_unique

  type psb_l_heap
    integer(psb_ipk_) :: last, dir
    integer(psb_lpk_), allocatable    :: keys(:)
  contains
    procedure, pass(heap) :: init       => psb_l_init_heap
    procedure, pass(heap) :: howmany    => psb_l_howmany
    procedure, pass(heap) :: insert     => psb_l_insert_heap
    procedure, pass(heap) :: get_first  => psb_l_heap_get_first
    procedure, pass(heap) :: dump       => psb_l_dump_heap
    procedure, pass(heap) :: free       => psb_l_free_heap   
  end type psb_l_heap

  type psb_l_idx_heap
    integer(psb_ipk_) :: last, dir
    integer(psb_lpk_), allocatable    :: keys(:)
    integer(psb_lpk_), allocatable :: idxs(:)
  contains
    procedure, pass(heap) :: init       => psb_l_idx_init_heap
    procedure, pass(heap) :: howmany    => psb_l_idx_howmany
    procedure, pass(heap) :: insert     => psb_l_idx_insert_heap
    procedure, pass(heap) :: get_first  => psb_l_idx_heap_get_first
    procedure, pass(heap) :: dump       => psb_l_idx_dump_heap
    procedure, pass(heap) :: free       => psb_l_idx_free_heap   
  end type psb_l_idx_heap


  interface psb_msort
    subroutine psb_lmsort(x,ix,dir,flag)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_lpk_), optional, intent(inout) :: ix(:)
    end subroutine psb_lmsort
  end interface psb_msort


  interface psb_bsrch
    function  psb_lbsrch(key,n,v) result(ipos)
      import 
      integer(psb_ipk_) :: ipos, n
      integer(psb_lpk_) :: key
      integer(psb_lpk_) :: v(:)
    end function psb_lbsrch
  end interface psb_bsrch

  interface psb_ssrch
    function psb_lssrch(key,n,v) result(ipos)
      import 
      implicit none
      integer(psb_ipk_) :: ipos, n
      integer(psb_lpk_) :: key
      integer(psb_lpk_) :: v(:)
    end function psb_lssrch
  end interface psb_ssrch

  interface 
    subroutine psi_l_msort_up(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      integer(psb_lpk_)  ::  k(n)
      integer(psb_lpk_) :: l(0:n+1)
    end subroutine psi_l_msort_up
    subroutine psi_l_msort_dw(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      integer(psb_lpk_)  ::  k(n)
      integer(psb_lpk_) :: l(0:n+1)
    end subroutine psi_l_msort_dw
  end interface
  interface 
    subroutine psi_l_amsort_up(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      integer(psb_lpk_)  ::  k(n)
      integer(psb_lpk_) :: l(0:n+1)
    end subroutine psi_l_amsort_up
    subroutine psi_l_amsort_dw(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      integer(psb_lpk_)  ::  k(n)
      integer(psb_lpk_) :: l(0:n+1)
    end subroutine psi_l_amsort_dw
  end interface
  
  
  interface psb_qsort
    subroutine psb_lqsort(x,ix,dir,flag)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_lpk_), optional, intent(inout) :: ix(:)
    end subroutine psb_lqsort
  end interface psb_qsort
  
  interface psb_isort
    subroutine psb_lisort(x,ix,dir,flag)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_lpk_), optional, intent(inout) :: ix(:)
    end subroutine psb_lisort
  end interface psb_isort


  interface psb_hsort
    subroutine psb_lhsort(x,ix,dir,flag)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_lpk_), optional, intent(inout) :: ix(:)
    end subroutine psb_lhsort
  end interface psb_hsort


  interface 
    subroutine psi_l_insert_heap(key,last,heap,dir,info)
      import 
      implicit none 

      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction

      integer(psb_lpk_), intent(in)     :: key
      integer(psb_lpk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_l_insert_heap
  end interface

  interface 
    subroutine psi_l_idx_insert_heap(key,index,last,heap,idxs,dir,info)
      import 
      implicit none 

      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction

      integer(psb_lpk_), intent(in)     :: key
      integer(psb_lpk_), intent(inout)  :: heap(:)
      integer(psb_lpk_), intent(in)     :: index
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_lpk_), intent(inout)  :: idxs(:)
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_l_idx_insert_heap
  end interface


  interface 
    subroutine psi_l_heap_get_first(key,last,heap,dir,info)
      import 
      implicit none 
      integer(psb_lpk_), intent(inout)  :: key
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_lpk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_l_heap_get_first
  end interface

  interface 
    subroutine psi_l_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import
      integer(psb_lpk_), intent(inout)    :: key
      integer(psb_lpk_), intent(out)    :: index
      integer(psb_lpk_), intent(inout)    :: heap(:)
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_lpk_), intent(inout)  :: idxs(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_l_idx_heap_get_first
  end interface

  interface 
    subroutine psi_lisrx_up(n,x,ix)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(inout) :: ix(:)
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_lisrx_up
    subroutine psi_lisrx_dw(n,x,ix)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(inout) :: ix(:)
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_lisrx_dw
    subroutine psi_lisr_up(n,x)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_lisr_up
    subroutine psi_lisr_dw(n,x)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_lisr_dw
    subroutine psi_laisrx_up(n,x,ix)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(inout) :: ix(:)
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_laisrx_up
    subroutine psi_laisrx_dw(n,x,ix)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(inout) :: ix(:)
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_laisrx_dw
    subroutine psi_laisr_up(n,x)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_laisr_up
    subroutine psi_laisr_dw(n,x)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_laisr_dw
  end interface

  interface 
    subroutine psi_lqsrx_up(n,x,ix)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(inout) :: ix(:)
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_lqsrx_up
    subroutine psi_lqsrx_dw(n,x,ix)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(inout) :: ix(:)
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_lqsrx_dw
    subroutine psi_lqsr_up(n,x)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_lqsr_up
    subroutine psi_lqsr_dw(n,x)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_lqsr_dw
    subroutine psi_laqsrx_up(n,x,ix)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(inout) :: ix(:)
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_laqsrx_up
    subroutine psi_laqsrx_dw(n,x,ix)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(inout) :: ix(:)
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_laqsrx_dw
    subroutine psi_laqsr_up(n,x)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_laqsr_up
    subroutine psi_laqsr_dw(n,x)
      import 
      integer(psb_lpk_), intent(inout)  :: x(:) 
      integer(psb_lpk_), intent(in)   :: n
    end subroutine psi_laqsr_dw
  end interface

contains

  subroutine psb_l_init_heap(heap,info,dir)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 
    class(psb_l_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)            :: info
    integer(psb_ipk_), intent(in), optional   :: dir

    info = psb_success_
    heap%last=0
    if (present(dir)) then 
      heap%dir = dir
    else
      heap%dir = psb_sort_up_
    endif
    select case(heap%dir) 
    case (psb_sort_up_,psb_sort_down_,psb_asort_up_,psb_asort_down_)
      ! ok, do nothing
    case default
      write(psb_err_unit,*) 'Invalid direction, defaulting to psb_sort_up_'
      heap%dir = psb_sort_up_
    end select
    call psb_ensure_size(psb_heap_resize,heap%keys,info)

    return
  end subroutine psb_l_init_heap


  function psb_l_howmany(heap) result(res)
    implicit none 
    class(psb_l_heap), intent(in) :: heap
    integer(psb_ipk_) :: res
    res  = heap%last
  end function psb_l_howmany

  subroutine psb_l_insert_heap(key,heap,info)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 

    integer(psb_lpk_), intent(in)              :: key
    class(psb_l_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)                       :: info

    info = psb_success_
    if (heap%last < 0) then 
      write(psb_err_unit,*) 'Invalid last in heap ',heap%last
      info = heap%last
      return
    endif

    call psb_ensure_size(heap%last+1,heap%keys,info,addsz=psb_heap_resize)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) 'Memory allocation failure in heap_insert'
      info = -5
      return
    end if
    call psi_l_insert_heap(key,&
         & heap%last,heap%keys,heap%dir,info)

    return
  end subroutine psb_l_insert_heap

  subroutine psb_l_heap_get_first(key,heap,info)
    implicit none 

    class(psb_l_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)     :: info
    integer(psb_lpk_), intent(out)       :: key


    info = psb_success_

    call psi_l_heap_get_first(key,&
         & heap%last,heap%keys,heap%dir,info)

    return
  end subroutine psb_l_heap_get_first

  subroutine psb_l_dump_heap(iout,heap,info)

    implicit none 
    class(psb_l_heap), intent(in) :: heap
    integer(psb_ipk_), intent(out)    :: info
    integer(psb_ipk_), intent(in)     :: iout

    info = psb_success_
    if (iout < 0) then
      write(psb_err_unit,*) 'Invalid file '
      info =-1
      return
    end if

    write(iout,*) 'Heap direction ',heap%dir
    write(iout,*) 'Heap size      ',heap%last
    if ((heap%last > 0).and.((.not.allocated(heap%keys)).or.&
         & (size(heap%keys)<heap%last))) then
      write(iout,*) 'Inconsistent size/allocation status!!'
    else
      write(iout,*) heap%keys(1:heap%last)
    end if
  end subroutine psb_l_dump_heap

  subroutine psb_l_free_heap(heap,info)
    implicit none 
    class(psb_l_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)

  end subroutine psb_l_free_heap

  subroutine psb_l_idx_init_heap(heap,info,dir)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 
    class(psb_l_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)            :: info
    integer(psb_ipk_), intent(in), optional   :: dir

    info = psb_success_
    heap%last=0
    if (present(dir)) then 
      heap%dir = dir
    else
      heap%dir = psb_sort_up_
    endif
    select case(heap%dir) 
    case (psb_sort_up_,psb_sort_down_,psb_asort_up_,psb_asort_down_)
      ! ok, do nothing
    case default
      write(psb_err_unit,*) 'Invalid direction, defaulting to psb_sort_up_'
      heap%dir = psb_sort_up_
    end select

    call psb_ensure_size(psb_heap_resize,heap%keys,info)
    call psb_ensure_size(psb_heap_resize,heap%idxs,info)
    return
  end subroutine psb_l_idx_init_heap


  function psb_l_idx_howmany(heap) result(res)
    implicit none 
    class(psb_l_idx_heap), intent(in) :: heap
    integer(psb_ipk_) :: res
    res  = heap%last
  end function psb_l_idx_howmany

  subroutine psb_l_idx_insert_heap(key,index,heap,info)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 

    integer(psb_lpk_), intent(in)              :: key
    integer(psb_lpk_), intent(in)                        :: index
    class(psb_l_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)                       :: info

    info = psb_success_
    if (heap%last < 0) then 
      write(psb_err_unit,*) 'Invalid last in heap ',heap%last
      info = heap%last
      return
    endif

    call psb_ensure_size(heap%last+1,heap%keys,info,addsz=psb_heap_resize)
    if (info == psb_success_) &
         & call psb_ensure_size(heap%last+1,heap%idxs,info,addsz=psb_heap_resize)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) 'Memory allocation failure in heap_insert'
      info = -5
      return
    end if
    call psi_l_idx_insert_heap(key,index,&
         & heap%last,heap%keys,heap%idxs,heap%dir,info)

    return
  end subroutine psb_l_idx_insert_heap

  subroutine psb_l_idx_heap_get_first(key,index,heap,info)
    implicit none 

    class(psb_l_idx_heap), intent(inout) :: heap
    integer(psb_lpk_), intent(out)       :: index
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_lpk_), intent(out)           :: key


    info = psb_success_

    call psi_l_idx_heap_get_first(key,index,&
         & heap%last,heap%keys,heap%idxs,heap%dir,info)

    return
  end subroutine psb_l_idx_heap_get_first

  subroutine psb_l_idx_dump_heap(iout,heap,info)

    implicit none 
    class(psb_l_idx_heap), intent(in) :: heap
    integer(psb_ipk_), intent(out)    :: info
    integer(psb_ipk_), intent(in)     :: iout

    info = psb_success_
    if (iout < 0) then
      write(psb_err_unit,*) 'Invalid file '
      info =-1
      return
    end if

    write(iout,*) 'Heap direction ',heap%dir
    write(iout,*) 'Heap size      ',heap%last
    if ((heap%last > 0).and.((.not.allocated(heap%keys)).or.&
         & (size(heap%keys)<heap%last))) then
      write(iout,*) 'Inconsistent size/allocation status!!'
    else    if ((heap%last > 0).and.((.not.allocated(heap%idxs)).or.&
         & (size(heap%idxs)<heap%last))) then
      write(iout,*) 'Inconsistent size/allocation status!!'
    else
      write(iout,*) heap%keys(1:heap%last)
      write(iout,*) heap%idxs(1:heap%last)
    end if
  end subroutine psb_l_idx_dump_heap

  subroutine psb_l_idx_free_heap(heap,info)
    implicit none 
    class(psb_l_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)
    if ((info == psb_success_).and.(allocated(heap%idxs))) deallocate(heap%idxs,stat=info)

  end subroutine psb_l_idx_free_heap

end module psb_l_sort_mod
