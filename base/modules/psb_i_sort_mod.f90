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
module psb_i_sort_mod
  use psb_const_mod

  interface psb_iblsrch
    function  psb_iblsrch(key,n,v) result(ipos)
      import :: psb_ipk_
      integer(psb_ipk_) :: ipos, key, n
      integer(psb_ipk_) :: v(:)
    end function psb_iblsrch
  end interface psb_iblsrch

  interface psb_ibsrch
    function  psb_ibsrch(key,n,v) result(ipos)
      import :: psb_ipk_
      integer(psb_ipk_) :: ipos, key, n
      integer(psb_ipk_) :: v(:)
    end function psb_ibsrch
  end interface psb_ibsrch

  interface psb_issrch
    function psb_issrch(key,n,v) result(ipos)
      import :: psb_ipk_
      implicit none
      integer(psb_ipk_) :: ipos, key, n
      integer(psb_ipk_) :: v(:)
    end function psb_issrch
  end interface psb_issrch

  interface psb_isaperm
    logical function psb_isaperm(n,eip)               
      import :: psb_ipk_
      integer(psb_ipk_), intent(in) :: n                             
      integer(psb_ipk_), intent(in) :: eip(n)
    end function psb_isaperm
  end interface psb_isaperm

  interface psb_msort_unique
    subroutine psb_imsort_u(x,nout,dir)
      import :: psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(inout)           :: x(:) 
      integer(psb_ipk_), intent(out)             :: nout
      integer(psb_ipk_), optional, intent(in)    :: dir
    end subroutine psb_imsort_u
  end interface psb_msort_unique

  type psb_i_heap
    integer(psb_ipk_) :: last, dir
    integer(psb_ipk_), allocatable    :: keys(:)
  contains
    procedure, pass(heap) :: init       => psb_i_init_heap
    procedure, pass(heap) :: howmany    => psb_i_howmany
    procedure, pass(heap) :: insert     => psb_i_insert_heap
    procedure, pass(heap) :: get_first  => psb_i_heap_get_first
    procedure, pass(heap) :: dump       => psb_i_dump_heap
    procedure, pass(heap) :: free       => psb_i_free_heap   
  end type psb_i_heap

  type psb_i_idx_heap
    integer(psb_ipk_) :: last, dir
    integer(psb_ipk_), allocatable    :: keys(:)
    integer(psb_ipk_), allocatable :: idxs(:)
  contains
    procedure, pass(heap) :: init       => psb_i_idx_init_heap
    procedure, pass(heap) :: howmany    => psb_i_idx_howmany
    procedure, pass(heap) :: insert     => psb_i_idx_insert_heap
    procedure, pass(heap) :: get_first  => psb_i_idx_heap_get_first
    procedure, pass(heap) :: dump       => psb_i_idx_dump_heap
    procedure, pass(heap) :: free       => psb_i_idx_free_heap   
  end type psb_i_idx_heap


  interface psb_msort
    subroutine psb_imsort(x,ix,dir,flag)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_imsort
  end interface psb_msort

  interface 
    subroutine psi_i_msort_up(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      integer(psb_ipk_)  ::  k(n)
      integer(psb_ipk_) :: l(0:n+1)
    end subroutine psi_i_msort_up
    subroutine psi_i_msort_dw(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      integer(psb_ipk_)  ::  k(n)
      integer(psb_ipk_) :: l(0:n+1)
    end subroutine psi_i_msort_dw
  end interface
  interface 
    subroutine psi_i_amsort_up(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      integer(psb_ipk_)  ::  k(n)
      integer(psb_ipk_) :: l(0:n+1)
    end subroutine psi_i_amsort_up
    subroutine psi_i_amsort_dw(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      integer(psb_ipk_)  ::  k(n)
      integer(psb_ipk_) :: l(0:n+1)
    end subroutine psi_i_amsort_dw
  end interface
  
  
  interface psb_qsort
    subroutine psb_iqsort(x,ix,dir,flag)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_iqsort
  end interface psb_qsort
  
  interface psb_isort
    subroutine psb_iisort(x,ix,dir,flag)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_iisort
  end interface psb_isort


  interface psb_hsort
    subroutine psb_ihsort(x,ix,dir,flag)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_ihsort
  end interface psb_hsort


!!$  interface !psb_howmany_heap
!!$    module procedure psb_i_howmany,  psb_i_idx_howmany
!!$  end interface 
!!$
!!$
!!$  interface !psb_init_heap
!!$    module procedure psb_i_init_heap, psb_i_idx_init_heap
!!$  end interface 
!!$
!!$
!!$  interface !psb_dump_heap
!!$    module procedure psb_i_dump_heap, psb_dump_i_idx_heap
!!$  end interface 
!!$
!!$
!!$  interface !psb_insert_heap
!!$    module procedure psb_i_insert_heap,  psb_i_idx_insert_heap
!!$  end interface 
!!$
!!$  interface !psb_heap_get_first
!!$    module procedure psb_i_heap_get_first, psb_i_idx_heap_get_first
!!$  end interface 
!!$  
!!$  interface !psb_free_heap
!!$    module procedure psb_free_i_heap, psb_free_i_idx_heap
!!$  end interface 

  interface 
    subroutine psi_i_insert_heap(key,last,heap,dir,info)
      import 
      implicit none 

      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction

      integer(psb_ipk_), intent(in)     :: key
      integer(psb_ipk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_i_insert_heap
  end interface

  interface 
    subroutine psi_i_idx_insert_heap(key,index,last,heap,idxs,dir,info)
      import 
      implicit none 

      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction

      integer(psb_ipk_), intent(in)     :: key
      integer(psb_ipk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(in)     :: index
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: idxs(:)
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_i_idx_insert_heap
  end interface


  interface 
    subroutine psi_i_heap_get_first(key,last,heap,dir,info)
      import 
      implicit none 
      integer(psb_ipk_), intent(inout)  :: key
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_i_heap_get_first
  end interface

  interface 
    subroutine psi_i_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import
      integer(psb_ipk_), intent(inout)    :: key
      integer(psb_ipk_), intent(out)    :: index
      integer(psb_ipk_), intent(inout)    :: heap(:)
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(inout)  :: idxs(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_i_idx_heap_get_first
  end interface

  interface 
    subroutine psi_iisrx_up(n,x,ix)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iisrx_up
    subroutine psi_iisrx_dw(n,x,ix)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iisrx_dw
    subroutine psi_iisr_up(n,x)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iisr_up
    subroutine psi_iisr_dw(n,x)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iisr_dw
    subroutine psi_iaisrx_up(n,x,ix)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iaisrx_up
    subroutine psi_iaisrx_dw(n,x,ix)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iaisrx_dw
    subroutine psi_iaisr_up(n,x)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iaisr_up
    subroutine psi_iaisr_dw(n,x)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iaisr_dw
  end interface

  interface 
    subroutine psi_iqsrx_up(n,x,ix)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iqsrx_up
    subroutine psi_iqsrx_dw(n,x,ix)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iqsrx_dw
    subroutine psi_iqsr_up(n,x)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iqsr_up
    subroutine psi_iqsr_dw(n,x)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iqsr_dw
    subroutine psi_iaqsrx_up(n,x,ix)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iaqsrx_up
    subroutine psi_iaqsrx_dw(n,x,ix)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iaqsrx_dw
    subroutine psi_iaqsr_up(n,x)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iaqsr_up
    subroutine psi_iaqsr_dw(n,x)
      import 
      integer(psb_ipk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_iaqsr_dw
  end interface

contains

  subroutine psb_i_init_heap(heap,info,dir)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 
    class(psb_i_heap), intent(inout) :: heap
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
  end subroutine psb_i_init_heap


  function psb_i_howmany(heap) result(res)
    implicit none 
    class(psb_i_heap), intent(in) :: heap
    integer(psb_ipk_) :: res
    res  = heap%last
  end function psb_i_howmany

  subroutine psb_i_insert_heap(key,heap,info)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 

    integer(psb_ipk_), intent(in)              :: key
    class(psb_i_heap), intent(inout) :: heap
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
    call psi_i_insert_heap(key,&
         & heap%last,heap%keys,heap%dir,info)

    return
  end subroutine psb_i_insert_heap

  subroutine psb_i_heap_get_first(key,heap,info)
    implicit none 

    class(psb_i_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)     :: info
    integer(psb_ipk_), intent(out)       :: key


    info = psb_success_

    call psi_i_heap_get_first(key,&
         & heap%last,heap%keys,heap%dir,info)

    return
  end subroutine psb_i_heap_get_first

  subroutine psb_i_dump_heap(iout,heap,info)

    implicit none 
    class(psb_i_heap), intent(in) :: heap
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
  end subroutine psb_i_dump_heap

  subroutine psb_i_free_heap(heap,info)
    implicit none 
    class(psb_i_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)

  end subroutine psb_i_free_heap

  subroutine psb_i_idx_init_heap(heap,info,dir)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 
    class(psb_i_idx_heap), intent(inout) :: heap
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
  end subroutine psb_i_idx_init_heap


  function psb_i_idx_howmany(heap) result(res)
    implicit none 
    class(psb_i_idx_heap), intent(in) :: heap
    integer(psb_ipk_) :: res
    res  = heap%last
  end function psb_i_idx_howmany

  subroutine psb_i_idx_insert_heap(key,index,heap,info)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 

    integer(psb_ipk_), intent(in)              :: key
    integer(psb_ipk_), intent(in)                        :: index
    class(psb_i_idx_heap), intent(inout) :: heap
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
    call psi_i_idx_insert_heap(key,index,&
         & heap%last,heap%keys,heap%idxs,heap%dir,info)

    return
  end subroutine psb_i_idx_insert_heap

  subroutine psb_i_idx_heap_get_first(key,index,heap,info)
    implicit none 

    class(psb_i_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)       :: index,info
    integer(psb_ipk_), intent(out)           :: key


    info = psb_success_

    call psi_i_idx_heap_get_first(key,index,&
         & heap%last,heap%keys,heap%idxs,heap%dir,info)

    return
  end subroutine psb_i_idx_heap_get_first

  subroutine psb_i_idx_dump_heap(iout,heap,info)

    implicit none 
    class(psb_i_idx_heap), intent(in) :: heap
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
  end subroutine psb_i_idx_dump_heap

  subroutine psb_i_idx_free_heap(heap,info)
    implicit none 
    class(psb_i_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)
    if ((info == psb_success_).and.(allocated(heap%idxs))) deallocate(heap%idxs,stat=info)

  end subroutine psb_i_idx_free_heap

end module psb_i_sort_mod
