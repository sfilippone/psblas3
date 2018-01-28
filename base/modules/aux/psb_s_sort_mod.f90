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
module psb_s_sort_mod
  use psb_const_mod


  type psb_s_heap
    integer(psb_ipk_) :: last, dir
    real(psb_spk_), allocatable    :: keys(:)
  contains
    procedure, pass(heap) :: init       => psb_s_init_heap
    procedure, pass(heap) :: howmany    => psb_s_howmany
    procedure, pass(heap) :: insert     => psb_s_insert_heap
    procedure, pass(heap) :: get_first  => psb_s_heap_get_first
    procedure, pass(heap) :: dump       => psb_s_dump_heap
    procedure, pass(heap) :: free       => psb_s_free_heap   
  end type psb_s_heap

  type psb_s_idx_heap
    integer(psb_ipk_) :: last, dir
    real(psb_spk_), allocatable    :: keys(:)
    integer(psb_ipk_), allocatable :: idxs(:)
  contains
    procedure, pass(heap) :: init       => psb_s_idx_init_heap
    procedure, pass(heap) :: howmany    => psb_s_idx_howmany
    procedure, pass(heap) :: insert     => psb_s_idx_insert_heap
    procedure, pass(heap) :: get_first  => psb_s_idx_heap_get_first
    procedure, pass(heap) :: dump       => psb_s_idx_dump_heap
    procedure, pass(heap) :: free       => psb_s_idx_free_heap   
  end type psb_s_idx_heap


  interface psb_msort
    subroutine psb_smsort(x,ix,dir,flag)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_smsort
  end interface psb_msort

  interface 
    subroutine psi_s_msort_up(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      real(psb_spk_)  ::  k(n)
      integer(psb_ipk_) :: l(0:n+1)
    end subroutine psi_s_msort_up
    subroutine psi_s_msort_dw(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      real(psb_spk_)  ::  k(n)
      integer(psb_ipk_) :: l(0:n+1)
    end subroutine psi_s_msort_dw
  end interface
  interface 
    subroutine psi_s_amsort_up(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      real(psb_spk_)  ::  k(n)
      integer(psb_ipk_) :: l(0:n+1)
    end subroutine psi_s_amsort_up
    subroutine psi_s_amsort_dw(n,k,l,iret)
      import
      implicit none
      integer(psb_ipk_) :: n, iret
      real(psb_spk_)  ::  k(n)
      integer(psb_ipk_) :: l(0:n+1)
    end subroutine psi_s_amsort_dw
  end interface
  
  
  interface psb_qsort
    subroutine psb_sqsort(x,ix,dir,flag)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_sqsort
  end interface psb_qsort
  
  interface psb_isort
    subroutine psb_sisort(x,ix,dir,flag)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_sisort
  end interface psb_isort


  interface psb_hsort
    subroutine psb_shsort(x,ix,dir,flag)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), optional, intent(in)    :: dir, flag
      integer(psb_ipk_), optional, intent(inout) :: ix(:)
    end subroutine psb_shsort
  end interface psb_hsort


  interface 
    subroutine psi_s_insert_heap(key,last,heap,dir,info)
      import 
      implicit none 

      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction

      real(psb_spk_), intent(in)     :: key
      real(psb_spk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_s_insert_heap
  end interface

  interface 
    subroutine psi_s_idx_insert_heap(key,index,last,heap,idxs,dir,info)
      import 
      implicit none 

      !  
      ! Input: 
      !   key:  the new value
      !   last: pointer to the last occupied element in heap
      !   heap: the heap
      !   dir:  sorting direction

      real(psb_spk_), intent(in)     :: key
      real(psb_spk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(in)     :: index
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: idxs(:)
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_s_idx_insert_heap
  end interface


  interface 
    subroutine psi_s_heap_get_first(key,last,heap,dir,info)
      import 
      implicit none 
      real(psb_spk_), intent(inout)  :: key
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(in)     :: dir
      real(psb_spk_), intent(inout)  :: heap(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_s_heap_get_first
  end interface

  interface 
    subroutine psi_s_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
      import
      real(psb_spk_), intent(inout)    :: key
      integer(psb_ipk_), intent(out)    :: index
      real(psb_spk_), intent(inout)    :: heap(:)
      integer(psb_ipk_), intent(in)     :: dir
      integer(psb_ipk_), intent(inout)  :: last
      integer(psb_ipk_), intent(inout)  :: idxs(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psi_s_idx_heap_get_first
  end interface

  interface 
    subroutine psi_sisrx_up(n,x,ix)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_sisrx_up
    subroutine psi_sisrx_dw(n,x,ix)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_sisrx_dw
    subroutine psi_sisr_up(n,x)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_sisr_up
    subroutine psi_sisr_dw(n,x)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_sisr_dw
    subroutine psi_saisrx_up(n,x,ix)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_saisrx_up
    subroutine psi_saisrx_dw(n,x,ix)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_saisrx_dw
    subroutine psi_saisr_up(n,x)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_saisr_up
    subroutine psi_saisr_dw(n,x)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_saisr_dw
  end interface

  interface 
    subroutine psi_sqsrx_up(n,x,ix)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_sqsrx_up
    subroutine psi_sqsrx_dw(n,x,ix)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_sqsrx_dw
    subroutine psi_sqsr_up(n,x)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_sqsr_up
    subroutine psi_sqsr_dw(n,x)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_sqsr_dw
    subroutine psi_saqsrx_up(n,x,ix)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_saqsrx_up
    subroutine psi_saqsrx_dw(n,x,ix)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(inout) :: ix(:)
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_saqsrx_dw
    subroutine psi_saqsr_up(n,x)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_saqsr_up
    subroutine psi_saqsr_dw(n,x)
      import 
      real(psb_spk_), intent(inout)  :: x(:) 
      integer(psb_ipk_), intent(in)   :: n
    end subroutine psi_saqsr_dw
  end interface

contains

  subroutine psb_s_init_heap(heap,info,dir)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 
    class(psb_s_heap), intent(inout) :: heap
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
  end subroutine psb_s_init_heap


  function psb_s_howmany(heap) result(res)
    implicit none 
    class(psb_s_heap), intent(in) :: heap
    integer(psb_ipk_) :: res
    res  = heap%last
  end function psb_s_howmany

  subroutine psb_s_insert_heap(key,heap,info)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 

    real(psb_spk_), intent(in)              :: key
    class(psb_s_heap), intent(inout) :: heap
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
    call psi_s_insert_heap(key,&
         & heap%last,heap%keys,heap%dir,info)

    return
  end subroutine psb_s_insert_heap

  subroutine psb_s_heap_get_first(key,heap,info)
    implicit none 

    class(psb_s_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)     :: info
    real(psb_spk_), intent(out)       :: key


    info = psb_success_

    call psi_s_heap_get_first(key,&
         & heap%last,heap%keys,heap%dir,info)

    return
  end subroutine psb_s_heap_get_first

  subroutine psb_s_dump_heap(iout,heap,info)

    implicit none 
    class(psb_s_heap), intent(in) :: heap
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
  end subroutine psb_s_dump_heap

  subroutine psb_s_free_heap(heap,info)
    implicit none 
    class(psb_s_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)

  end subroutine psb_s_free_heap

  subroutine psb_s_idx_init_heap(heap,info,dir)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 
    class(psb_s_idx_heap), intent(inout) :: heap
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
  end subroutine psb_s_idx_init_heap


  function psb_s_idx_howmany(heap) result(res)
    implicit none 
    class(psb_s_idx_heap), intent(in) :: heap
    integer(psb_ipk_) :: res
    res  = heap%last
  end function psb_s_idx_howmany

  subroutine psb_s_idx_insert_heap(key,index,heap,info)
    use psb_realloc_mod, only : psb_ensure_size
    implicit none 

    real(psb_spk_), intent(in)              :: key
    integer(psb_ipk_), intent(in)                        :: index
    class(psb_s_idx_heap), intent(inout) :: heap
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
    call psi_s_idx_insert_heap(key,index,&
         & heap%last,heap%keys,heap%idxs,heap%dir,info)

    return
  end subroutine psb_s_idx_insert_heap

  subroutine psb_s_idx_heap_get_first(key,index,heap,info)
    implicit none 

    class(psb_s_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)       :: index,info
    real(psb_spk_), intent(out)           :: key


    info = psb_success_

    call psi_s_idx_heap_get_first(key,index,&
         & heap%last,heap%keys,heap%idxs,heap%dir,info)

    return
  end subroutine psb_s_idx_heap_get_first

  subroutine psb_s_idx_dump_heap(iout,heap,info)

    implicit none 
    class(psb_s_idx_heap), intent(in) :: heap
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
  end subroutine psb_s_idx_dump_heap

  subroutine psb_s_idx_free_heap(heap,info)
    implicit none 
    class(psb_s_idx_heap), intent(inout) :: heap
    integer(psb_ipk_), intent(out)           :: info

    info=psb_success_
    if (allocated(heap%keys)) deallocate(heap%keys,stat=info)
    if ((info == psb_success_).and.(allocated(heap%idxs))) deallocate(heap%idxs,stat=info)

  end subroutine psb_s_idx_free_heap

end module psb_s_sort_mod
