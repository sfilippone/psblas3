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
subroutine psb_dhsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => psb_dhsort
  use psb_error_mod
  implicit none 
  real(psb_dpk_), intent(inout)           :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: flag_, n, i, err_act,info
  integer(psb_ipk_) :: dir_, l
  real(psb_dpk_) :: key
  integer(psb_ipk_) :: index

  integer(psb_ipk_)  :: ierr(5)
  character(len=20)  :: name

  name='psb_hsort'
  call psb_erractionsave(err_act)

  if (present(flag)) then 
    flag_ = flag
  else 
    flag_ = psb_sort_ovw_idx_
  end if
  select case(flag_) 
  case( psb_sort_ovw_idx_, psb_sort_keep_idx_)
    ! OK keep going
  case default
    ierr(1) = 4; ierr(2) = flag_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  if (present(dir)) then 
    dir_ = dir
  else
    dir_= psb_sort_up_
  end if

  select case(dir_)
  case(psb_sort_up_,psb_sort_down_)
    ! OK
  case (psb_asort_up_,psb_asort_down_) 
    ! OK    
  case default
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  n = size(x)

  !
  ! Dirty trick to sort with heaps: if we want 
  ! to sort in place upwards, first we set up a heap so that
  ! we can easily get the LARGEST element, then we take it out 
  ! and put it in the last entry, and so on. 
  ! So,  we invert dir_
  !
  dir_ = -dir_ 

  if (present(ix)) then 
    if (size(ix) < n) then 
      ierr(1) = 2; ierr(2) = size(ix); 
      call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
      goto 9999
    end if
    if (flag_ == psb_sort_ovw_idx_) then 
      do i=1, n
        ix(i) = i
      end do
    end if
    l = 0
    do i=1, n 
      key   = x(i)
      index = ix(i)
      call psi_idx_insert_heap(key,index,l,x,ix,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! '
      end if
    end do
    do i=n, 2, -1 
      call psi_idx_heap_get_first(key,index,l,x,ix,dir_,info)
      if (l /= i-1) then 
        write(psb_err_unit,*) 'Mismatch while pulling out of heap ',l,i
      end if
      x(i)  = key
      ix(i) = index
    end do
  else if (.not.present(ix)) then 
    l = 0
    do i=1, n 
      key   = x(i)
      call psi_insert_heap(key,l,x,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! ',l,i
      end if
    end do
    do i=n, 2, -1 
      call psi_d_heap_get_first(key,l,x,dir_,info)
      if (l /= i-1) then 
        write(psb_err_unit,*) 'Mismatch while pulling out of heap ',l,i
      end if
      x(i)  = key
    end do
  end if


  return

9999 call psb_error_handler(err_act)

  return
end subroutine psb_dhsort



!
! These are packaged so that they can be used to implement 
! a heapsort.
!
!
!   Programming note:
!   In the implementation of the heap_get_first function
!   we have code like this
!
!      if ( ( heap(2*i) < heap(2*i+1) ) .or.&
!           & (2*i == last)) then 
!        j = 2*i
!      else
!        j = 2*i + 1
!      end if
!
!   It looks like the 2*i+1 could overflow the array, but this
!   is not true because there is a guard statement
!       if (i>last/2) exit
!   and because last has just been reduced by 1 when defining the return value,
!   therefore 2*i+1 may be greater than the current value of last,
!   but cannot be greater than the value of last when the routine was entered
!   hence it is safe.
!
!
!

subroutine psi_d_insert_heap(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_d_insert_heap
  implicit none 

  !  
  ! Input: 
  !   key:  the new value
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   dir:  sorting direction

  real(psb_dpk_), intent(in)    :: key
  integer(psb_ipk_), intent(in)           :: dir
  real(psb_dpk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(inout)        :: last
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_) :: i, i2
  real(psb_dpk_)   :: temp

  info = psb_success_
  if (last < 0) then 
    write(psb_err_unit,*) 'Invalid last in heap ',last
    info = last
    return
  endif
  last    = last + 1
  if (last > size(heap)) then 
    write(psb_err_unit,*) 'out of bounds '
    info = -1
    return
  end if
  i       = last
  heap(i) = key

  select case(dir)
  case (psb_sort_up_)

    do 
      if (i<=1) exit
      i2 = i/2
      if (heap(i) < heap(i2)) then 
        temp     = heap(i)
        heap(i)  = heap(i2)
        heap(i2) = temp
        i        = i2
      else
        exit
      end if
    end do


  case (psb_sort_down_)

    do 
      if (i<=1) exit
      i2 = i/2
      if (heap(i) > heap(i2)) then 
        temp     = heap(i)
        heap(i)  = heap(i2)
        heap(i2) = temp
        i        = i2
      else
        exit
      end if
    end do

  case (psb_asort_up_)

    do 
      if (i<=1) exit
      i2 = i/2
      if (abs(heap(i)) < abs(heap(i2))) then 
        temp     = heap(i)
        heap(i)  = heap(i2)
        heap(i2) = temp
        i        = i2
      else
        exit
      end if
    end do


  case (psb_asort_down_)

    do 
      if (i<=1) exit
      i2 = i/2
      if (abs(heap(i)) > abs(heap(i2))) then 
        temp     = heap(i)
        heap(i)  = heap(i2)
        heap(i2) = temp
        i        = i2
      else
        exit
      end if
    end do


  case default
    write(psb_err_unit,*) 'Invalid direction in heap ',dir
  end select

  return
end subroutine psi_d_insert_heap


subroutine psi_d_heap_get_first(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_d_heap_get_first
  implicit none 

  real(psb_dpk_), intent(inout)     :: key
  integer(psb_ipk_), intent(inout)  :: last
  integer(psb_ipk_), intent(in)     :: dir
  real(psb_dpk_), intent(inout)     :: heap(:)
  integer(psb_ipk_), intent(out)    :: info

  integer(psb_ipk_) :: i, j
  real(psb_dpk_)   :: temp


  info = psb_success_
  if (last <= 0) then 
    key  = 0
    info = -1
    return
  endif

  key     = heap(1)
  heap(1) = heap(last)
  last    = last - 1

  select case(dir)
  case (psb_sort_up_)

    i = 1
    do 
      if (i > (last/2)) exit
      if ( (heap(2*i) < heap(2*i+1)) .or.&
           & (2*i == last)) then 
        j = 2*i
      else
        j = 2*i + 1
      end if

      if (heap(i) > heap(j)) then 
        temp     = heap(i)
        heap(i)  = heap(j)
        heap(j)  = temp
        i        = j 
      else
        exit
      end if
    end do


  case (psb_sort_down_)

    i = 1
    do 
      if (i > (last/2)) exit
      if ( (heap(2*i) > heap(2*i+1)) .or.&
           & (2*i == last)) then 
        j = 2*i
      else
        j = 2*i + 1
      end if

      if (heap(i) < heap(j)) then 
        temp     = heap(i)
        heap(i)  = heap(j)
        heap(j)  = temp
        i        = j 
      else
        exit
      end if
    end do

  case (psb_asort_up_)

    i = 1
    do 
      if (i > (last/2)) exit
      if ( (abs(heap(2*i)) < abs(heap(2*i+1))) .or.&
           & (2*i == last)) then 
        j = 2*i
      else
        j = 2*i + 1
      end if

      if (abs(heap(i)) > abs(heap(j))) then 
        temp     = heap(i)
        heap(i)  = heap(j)
        heap(j)  = temp
        i        = j 
      else
        exit
      end if
    end do


  case (psb_asort_down_)

    i = 1
    do 
      if (i > (last/2)) exit
      if ( (abs(heap(2*i)) > abs(heap(2*i+1))) .or.&
           & (2*i == last)) then 
        j = 2*i
      else
        j = 2*i + 1
      end if

      if (abs(heap(i)) < abs(heap(j))) then 
        temp     = heap(i)
        heap(i)  = heap(j)
        heap(j)  = temp
        i        = j 
      else
        exit
      end if
    end do

  case default
    write(psb_err_unit,*) 'Invalid direction in heap ',dir
  end select

  return
end subroutine psi_d_heap_get_first


subroutine psi_d_idx_insert_heap(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_d_idx_insert_heap

  implicit none 
  !  
  ! Input: 
  !   key:  the new value
  !   index: the new index
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   idxs: the indices
  !   dir:  sorting direction

  real(psb_dpk_), intent(in)       :: key
  integer(psb_ipk_), intent(in)     :: index,dir
  real(psb_dpk_), intent(inout)    :: heap(:)
  integer(psb_ipk_), intent(inout)  :: idxs(:),last
  integer(psb_ipk_), intent(out)    :: info
  integer(psb_ipk_) :: i, i2, itemp
  real(psb_dpk_)   :: temp 

  info = psb_success_
  if (last < 0) then 
    write(psb_err_unit,*) 'Invalid last in heap ',last
    info = last
    return
  endif

  last    = last + 1
  if (last > size(heap)) then 
    write(psb_err_unit,*) 'out of bounds '
    info = -1
    return
  end if

  i       = last
  heap(i) = key
  idxs(i) = index

  select case(dir)
  case (psb_sort_up_)

    do 
      if (i<=1) exit
      i2 = i/2
      if (heap(i) < heap(i2)) then 
        itemp    = idxs(i)
        idxs(i)  = idxs(i2)
        idxs(i2) = itemp
        temp     = heap(i)
        heap(i)  = heap(i2)
        heap(i2) = temp 
        i        = i2
      else
        exit
      end if
    end do


  case (psb_sort_down_)

    do 
      if (i<=1) exit
      i2 = i/2
      if (heap(i) > heap(i2)) then 
        itemp    = idxs(i)
        idxs(i)  = idxs(i2)
        idxs(i2) = itemp
        temp     = heap(i)
        heap(i)  = heap(i2)
        heap(i2) = temp 
        i        = i2
      else
        exit
      end if
    end do

  case (psb_asort_up_)

    do 
      if (i<=1) exit
      i2 = i/2
      if (abs(heap(i)) < abs(heap(i2))) then 
        itemp    = idxs(i)
        idxs(i)  = idxs(i2)
        idxs(i2) = itemp
        temp     = heap(i)
        heap(i)  = heap(i2)
        heap(i2) = temp
        i        = i2
      else
        exit
      end if
    end do


  case (psb_asort_down_)

    do 
      if (i<=1) exit
      i2 = i/2
      if (abs(heap(i)) > abs(heap(i2))) then 
        itemp    = idxs(i)
        idxs(i)  = idxs(i2)
        idxs(i2) = itemp
        temp     = heap(i)
        heap(i)  = heap(i2)
        heap(i2) = temp 
        i        = i2
      else
        exit
      end if
    end do


  case default
    write(psb_err_unit,*) 'Invalid direction in heap ',dir
  end select

  return
end subroutine psi_d_idx_insert_heap

subroutine psi_d_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_d_idx_heap_get_first
  implicit none 

  real(psb_dpk_), intent(inout)    :: key
  real(psb_dpk_), intent(inout)    :: heap(:)
  integer(psb_ipk_), intent(out)    :: index
  integer(psb_ipk_), intent(out)    :: info
  integer(psb_ipk_), intent(inout)  :: last,idxs(:)
  integer(psb_ipk_), intent(in)     :: dir

  integer(psb_ipk_) :: i, j,itemp
  real(psb_dpk_)   :: temp

  info = psb_success_
  if (last <= 0) then 
    key   = 0
    index = 0
    info  = -1
    return
  endif

  key     = heap(1)
  index   = idxs(1)
  heap(1) = heap(last)
  idxs(1) = idxs(last)
  last    = last - 1

  select case(dir)
  case (psb_sort_up_)

    i = 1
    do 
      if (i > (last/2)) exit
      if ( (heap(2*i) < heap(2*i+1)) .or.&
           & (2*i == last)) then 
        j = 2*i
      else
        j = 2*i + 1
      end if

      if (heap(i) > heap(j)) then 
        itemp    = idxs(i)
        idxs(i)  = idxs(j)
        idxs(j)  = itemp
        temp     = heap(i)
        heap(i)  = heap(j)
        heap(j)  = temp 
        i        = j 
      else
        exit
      end if
    end do


  case (psb_sort_down_)

    i = 1
    do 
      if (i > (last/2)) exit
      if ( (heap(2*i) > heap(2*i+1)) .or.&
           & (2*i == last)) then 
        j = 2*i
      else
        j = 2*i + 1
      end if

      if (heap(i) < heap(j)) then 
        itemp    = idxs(i)
        idxs(i)  = idxs(j)
        idxs(j)  = itemp
        temp     = heap(i)
        heap(i)  = heap(j)
        heap(j)  = temp 
        i        = j 
      else
        exit
      end if
    end do

  case (psb_asort_up_)

    i = 1
    do 
      if (i > (last/2)) exit
      if ( (abs(heap(2*i)) < abs(heap(2*i+1))) .or.&
           & (2*i == last)) then 
        j = 2*i
      else
        j = 2*i + 1
      end if

      if (abs(heap(i)) > abs(heap(j))) then 
        itemp    = idxs(i)
        idxs(i)  = idxs(j)
        idxs(j)  = itemp
        temp     = heap(i)
        heap(i)  = heap(j)
        heap(j)  = temp 
        i        = j 
      else
        exit
      end if
    end do


  case (psb_asort_down_)

    i = 1
    do 
      if (i > (last/2)) exit
      if ( (abs(heap(2*i)) > abs(heap(2*i+1))) .or.&
           & (2*i == last)) then 
        j = 2*i
      else
        j = 2*i + 1
      end if

      if (abs(heap(i)) < abs(heap(j))) then 
        itemp    = idxs(i)
        idxs(i)  = idxs(j)
        idxs(j)  = itemp
        temp     = heap(i)
        heap(i)  = heap(j)
        heap(j)  = temp 
        i        = j 
      else
        exit
      end if
    end do

  case default
    write(psb_err_unit,*) 'Invalid direction in heap ',dir
  end select

  return
end subroutine psi_d_idx_heap_get_first




