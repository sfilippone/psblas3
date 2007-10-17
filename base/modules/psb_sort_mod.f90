!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
module psb_sort_mod


  integer, parameter :: psb_sort_up_=1,  psb_sort_down_=-1
  integer, parameter :: psb_lsort_up_=2, psb_lsort_down_=-2
  integer, parameter :: psb_asort_up_=3, psb_asort_down_=-3
  integer, parameter :: psb_alsort_up_=4, psb_alsort_down_=-4
  integer, parameter :: psb_sort_ovw_idx_=0, psb_sort_keep_idx_=1
  integer, parameter :: psb_heap_resize=200

  type psb_int_heap
    integer              :: last, dir
    integer, allocatable :: keys(:)
  end type psb_int_heap
  type psb_int_idx_heap
    integer              :: last, dir
    integer, allocatable :: keys(:)
    integer, allocatable :: idxs(:)
  end type psb_int_idx_heap
  type psb_double_idx_heap
    integer              :: last, dir
    real(kind(1.d0)), allocatable :: keys(:)
    integer, allocatable          :: idxs(:)
  end type psb_double_idx_heap
  type psb_dcomplex_idx_heap
    integer              :: last, dir
    complex(kind(1.d0)), allocatable :: keys(:)
    integer, allocatable             :: idxs(:)
  end type psb_dcomplex_idx_heap

  interface psb_msort
    module procedure imsort
  end interface

  interface psb_msort_unique
    module procedure imsort_u
  end interface

  interface psb_qsort
    module procedure iqsort, dqsort, zqsort
  end interface

  interface psb_init_heap
    module procedure psb_init_int_heap, psb_init_int_idx_heap,&
         & psb_init_double_idx_heap, psb_init_dcomplex_idx_heap
  end interface

  interface psb_dump_heap
    module procedure psb_dump_int_heap, psb_dump_int_idx_heap,&
         & psb_dump_double_idx_heap, psb_dump_dcomplex_idx_heap
  end interface

  interface psb_howmany_heap
    module procedure psb_howmany_int_heap, psb_howmany_int_idx_heap,&
         & psb_howmany_double_idx_heap, psb_howmany_dcomplex_idx_heap
  end interface

  interface psb_insert_heap
    module procedure psb_insert_int_heap, psb_insert_int_idx_heap,&
         & psb_insert_double_idx_heap, psb_insert_dcomplex_idx_heap
  end interface

  interface psb_heap_get_first
    module procedure psb_int_heap_get_first, psb_int_idx_heap_get_first,&
         & psb_double_idx_heap_get_first, psb_dcomplex_idx_heap_get_first
  end interface


contains

  subroutine imsort(x,ix,dir,flag)
    use psb_error_mod
    implicit none 
    integer, intent(inout)           :: x(:) 
    integer, optional, intent(in)    :: dir, flag
    integer, optional, intent(inout) :: ix(:)
    
    integer  :: dir_, flag_, n, err_act
    
    character(len=20)  :: name

    name='psb_msort'
    call psb_erractionsave(err_act)

    if (present(dir)) then 
      dir_ = dir
    else
      dir_= psb_sort_up_
    end if
    select case(dir_) 
    case( psb_sort_up_, psb_sort_down_)
      ! OK keep going
    case default
      call psb_errpush(30,name,i_err=(/3,dir_,0,0,0/))
      goto 9999
    end select
      
    n = size(x)
 
    if (present(ix)) then 
      if (size(ix) < n) then 
        call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
        goto 9999
      end if
      if (present(flag)) then 
        flag_ = flag
      else 
        flag_ = psb_sort_ovw_idx_
      end if
      select case(flag_) 
      case( psb_sort_ovw_idx_, psb_sort_keep_idx_)
        ! OK keep going
      case default
        call psb_errpush(30,name,i_err=(/4,flag_,0,0,0/))
        goto 9999
      end select

      call imsrx(n,x,ix,dir_,flag_)
    else
      call imsr(n,x,dir_)
    end if

9999 continue 
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
  end subroutine imsort

  subroutine imsort_u(x,nout,dir)
    use psb_error_mod
    implicit none 
    integer, intent(inout)           :: x(:) 
    integer, intent(out)             :: nout
    integer, optional, intent(in)    :: dir
    
    integer  :: dir_, flag_, n, err_act
    
    character(len=20)  :: name

    name='psb_msort_u'
    call psb_erractionsave(err_act)

    if (present(dir)) then 
      dir_ = dir
    else
      dir_= psb_sort_up_
    end if
    select case(dir_) 
    case( psb_sort_up_, psb_sort_down_)
      ! OK keep going
    case default
      call psb_errpush(30,name,i_err=(/3,dir_,0,0,0/))
      goto 9999
    end select
      
    n = size(x)
 
    call imsru(n,x,dir_,nout)
      

9999 continue 
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
  end subroutine imsort_u


  subroutine iqsort(x,ix,dir,flag)
    use psb_error_mod
    implicit none 
    integer, intent(inout)           :: x(:) 
    integer, optional, intent(in)    :: dir, flag
    integer, optional, intent(inout) :: ix(:)
    
    integer  :: dir_, flag_, n, err_act
    
    character(len=20)  :: name

    name='psb_qsort'
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
      call psb_errpush(30,name,i_err=(/4,flag_,0,0,0/))
      goto 9999
    end select
    
    if (present(dir)) then 
      dir_ = dir
    else
      dir_= psb_sort_up_
    end if

    n = size(x)

    select case(dir_) 
    case( psb_sort_up_, psb_sort_down_)
      if (present(ix)) then 
        if (size(ix) < n) then 
          call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
          goto 9999
        end if
        
        call isrx(n,x,ix,dir_,flag_)
      else
        call isr(n,x,dir_)
      end if
      
    case( psb_asort_up_, psb_asort_down_)
      ! OK keep going
      if (present(ix)) then 
        if (size(ix) < n) then 
          call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
          goto 9999
        end if
        
        call iasrx(n,x,ix,dir_,flag_)
      else
        call iasr(n,x,dir_)
      end if
      
    case default
      call psb_errpush(30,name,i_err=(/3,dir_,0,0,0/))
      goto 9999
    end select
      
    

9999 continue 
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
  end subroutine iqsort


  subroutine dqsort(x,ix,dir,flag)
    use psb_error_mod
    implicit none 
    real(kind(1.d0)), intent(inout)  :: x(:) 
    integer, optional, intent(in)    :: dir, flag
    integer, optional, intent(inout) :: ix(:)
    
    integer  :: dir_, flag_, n, err_act
    
    character(len=20)  :: name

    name='psb_qsort'
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
      call psb_errpush(30,name,i_err=(/4,flag_,0,0,0/))
      goto 9999
    end select
    
    if (present(dir)) then 
      dir_ = dir
    else
      dir_= psb_sort_up_
    end if

    n = size(x)

    select case(dir_) 
    case( psb_sort_up_, psb_sort_down_)
      if (present(ix)) then 
        if (size(ix) < n) then 
          call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
          goto 9999
        end if
        
        call dsrx(n,x,ix,dir_,flag_)
      else
        call dsr(n,x,dir_)
      end if
      
    case( psb_asort_up_, psb_asort_down_)
      ! OK keep going
      if (present(ix)) then 
        if (size(ix) < n) then 
          call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
          goto 9999
        end if
        
        call dasrx(n,x,ix,dir_,flag_)
      else
        call dasr(n,x,dir_)
      end if
      
    case default
      call psb_errpush(30,name,i_err=(/3,dir_,0,0,0/))
      goto 9999
    end select
      
    

9999 continue 
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
  end subroutine dqsort


  subroutine zqsort(x,ix,dir,flag)
    use psb_error_mod
    implicit none 
    complex(kind(1.d0)), intent(inout)  :: x(:) 
    integer, optional, intent(in)    :: dir, flag
    integer, optional, intent(inout) :: ix(:)
    
    integer  :: dir_, flag_, n, err_act
    
    character(len=20)  :: name

    name='psb_qsort'
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
      call psb_errpush(30,name,i_err=(/4,flag_,0,0,0/))
      goto 9999
    end select
    
    if (present(dir)) then 
      dir_ = dir
    else
      dir_= psb_lsort_up_
    end if

    n = size(x)

    select case(dir_) 
    case( psb_lsort_up_, psb_lsort_down_)
      if (present(ix)) then 
        if (size(ix) < n) then 
          call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
          goto 9999
        end if
        
        call zlsrx(n,x,ix,dir_,flag_)
      else
        call zlsr(n,x,dir_)
      end if
      
    case( psb_alsort_up_, psb_alsort_down_)
      ! OK keep going
      if (present(ix)) then 
        if (size(ix) < n) then 
          call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
          goto 9999
        end if
        
        call zalsrx(n,x,ix,dir_,flag_)
      else
        call zalsr(n,x,dir_)
      end if

    case( psb_asort_up_, psb_asort_down_)
      ! OK keep going
      if (present(ix)) then 
        if (size(ix) < n) then 
          call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
          goto 9999
        end if
        
        call zasrx(n,x,ix,dir_,flag_)
      else
        call zasr(n,x,dir_)
      end if
      
    case default
      call psb_errpush(30,name,i_err=(/3,dir_,0,0,0/))
      goto 9999
    end select
      
    

9999 continue 
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
  end subroutine zqsort



  function  psb_howmany_int_heap(heap)
    implicit none 
    type(psb_int_heap), intent(in) :: heap
    integer :: psb_howmany_int_heap
    psb_howmany_int_heap = heap%last
  end function psb_howmany_int_heap

  subroutine psb_init_int_heap(heap,info,dir)
    use psb_realloc_mod
    implicit none 
    type(psb_int_heap), intent(inout) :: heap
    integer, intent(out)            :: info
    integer, intent(in), optional   :: dir
    
    info = 0
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
      write(0,*) 'Invalid direction, defaulting to psb_sort_up_'
      heap%dir = psb_sort_up_
    end select
    
    call psb_ensure_size(psb_heap_resize,heap%keys,info)
    return
  end subroutine psb_init_int_heap

  subroutine psb_dump_int_heap(iout,heap,info)
    implicit none 
    type(psb_int_heap), intent(in) :: heap
    integer, intent(out)           :: info
    integer, intent(in)            :: iout
    
    info = 0
    if (iout < 0) then
      write(0,*) 'Invalid file '
      info =-1
      return
    end if
    
    write(iout,*) 'Heap direction ',heap%dir
    write(iout,*) 'Heap size      ',heap%last
    if ((heap%last > 0).and.((.not.allocated(heap%keys)).or.&
         & (size(heap%keys)<heap%last))) then
      write(iout,*) 'Inconsistent size/allocation status!!'
    else
      if (heap%last > 0) then 
        write(iout,*) heap%keys(1:heap%last)
      end if
    end if
  end subroutine psb_dump_int_heap

  subroutine psb_insert_int_heap(key,heap,info)
    use psb_realloc_mod
    implicit none 

    integer, intent(in)               :: key
    type(psb_int_heap), intent(inout) :: heap
    integer, intent(out)              :: info
    integer                           :: i, i2
    integer                           :: temp
    info = 0
    if (heap%last < 0) then 
      write(0,*) 'Invalid last in heap ',heap%last
      info = heap%last
      return
    endif

    heap%last = heap%last + 1
    call psb_ensure_size(heap%last,heap%keys,info,addsz=psb_heap_resize)
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in heap_insert'
      info = -5
      return
    end if

    i            = heap%last
    heap%keys(i) = key
    
    select case(heap%dir)
    case (psb_sort_up_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (heap%keys(i) < heap%keys(i2)) then 
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp
          i             = i2
        else
          exit
        end if
      end do
      
          
    case (psb_sort_down_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (heap%keys(i) > heap%keys(i2)) then 
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp
          i             = i2
        else
          exit
        end if
      end do

    case (psb_asort_up_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (abs(heap%keys(i)) < abs(heap%keys(i2))) then 
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp
          i             = i2
        else
          exit
        end if
      end do
      
          
    case (psb_asort_down_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (abs(heap%keys(i)) > abs(heap%keys(i2))) then 
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp
          i             = i2
        else
          exit
        end if
      end do
      

    case default
      write(0,*) 'Invalid direction in heap ',heap%dir
    end select

    return
  end subroutine psb_insert_int_heap


  subroutine psb_int_heap_get_first(key,heap,info)
    implicit none 

    type(psb_int_heap), intent(inout) :: heap
    integer, intent(out)              :: key,info
    
    integer                           :: i, i2, last,j
    integer                           :: temp

    
    info = 0
    if (heap%last <= 0) then 
      key  = 0
      info = -1
      return
    endif

    key          = heap%keys(1)
    heap%keys(1) = heap%keys(heap%last)
    heap%last    = heap%last - 1
    last         = heap%last

    select case(heap%dir)
    case (psb_sort_up_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (heap%keys(2*i) < heap%keys(2*i+1)) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (heap%keys(i) > heap%keys(j)) then 
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp
          i             = j 
        else
          exit
        end if
      end do
      
          
    case (psb_sort_down_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (heap%keys(2*i) > heap%keys(2*i+1)) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (heap%keys(i) < heap%keys(j)) then 
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp
          i             = j 
        else
          exit
        end if
      end do

    case (psb_asort_up_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (abs(heap%keys(2*i)) < abs(heap%keys(2*i+1))) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (abs(heap%keys(i)) > abs(heap%keys(j))) then 
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp
          i             = j 
        else
          exit
        end if
      end do
      
          
    case (psb_asort_down_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (abs(heap%keys(2*i)) > abs(heap%keys(2*i+1))) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if
        
        if (abs(heap%keys(i)) < abs(heap%keys(j))) then 
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp
          i             = j 
        else
          exit
        end if
      end do

    case default
      write(0,*) 'Invalid direction in heap ',heap%dir
    end select

    return
  end subroutine psb_int_heap_get_first



  function  psb_howmany_double_idx_heap(heap)
    implicit none 
    type(psb_double_idx_heap), intent(in) :: heap
    integer :: psb_howmany_double_idx_heap
    psb_howmany_double_idx_heap = heap%last
  end function psb_howmany_double_idx_heap

  subroutine psb_init_double_idx_heap(heap,info,dir)
    use psb_realloc_mod
    implicit none 
    type(psb_double_idx_heap), intent(inout) :: heap
    integer, intent(out)            :: info
    integer, intent(in), optional   :: dir
    
    info = 0
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
      write(0,*) 'Invalid direction, defaulting to psb_sort_up_'
      heap%dir = psb_sort_up_
    end select
    
    call psb_ensure_size(psb_heap_resize,heap%keys,info)
    call psb_ensure_size(psb_heap_resize,heap%idxs,info)
    return
  end subroutine psb_init_double_idx_heap

  subroutine psb_dump_double_idx_heap(iout,heap,info)
    implicit none 
    type(psb_double_idx_heap), intent(in) :: heap
    integer, intent(out)           :: info
    integer, intent(in)            :: iout
    
    info = 0
    if (iout < 0) then
      write(0,*) 'Invalid file '
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
      if (heap%last > 0) then 
        write(iout,*) heap%keys(1:heap%last)
        write(iout,*) heap%idxs(1:heap%last)
      end if
    end if
  end subroutine psb_dump_double_idx_heap

  subroutine psb_insert_double_idx_heap(key,index,heap,info)
    use psb_realloc_mod
    implicit none 

    real(kind(1.d0)), intent(in)      :: key
    integer, intent(in)               :: index
    type(psb_double_idx_heap), intent(inout) :: heap
    integer, intent(out)              :: info
    integer                           :: i, i2, itemp
    real(kind(1.d0))                  :: temp 
    info = 0
    if (heap%last < 0) then 
      write(0,*) 'Invalid last in heap ',heap%last
      info = heap%last
      return
    endif

    heap%last = heap%last + 1
    call psb_ensure_size(heap%last,heap%keys,info,addsz=psb_heap_resize)
    if (info == 0) &
         & call psb_ensure_size(heap%last,heap%idxs,info,addsz=psb_heap_resize)
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in heap_insert'
      info = -5
      return
    end if
    
    i            = heap%last
    heap%keys(i) = key
    heap%idxs(i) = index

    select case(heap%dir)
    case (psb_sort_up_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (heap%keys(i) < heap%keys(i2)) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp 
          i             = i2
        else
          exit
        end if
      end do
      
          
    case (psb_sort_down_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (heap%keys(i) > heap%keys(i2)) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp 
          i             = i2
        else
          exit
        end if
      end do
      
    case (psb_asort_up_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (abs(heap%keys(i)) < abs(heap%keys(i2))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp
          i             = i2
        else
          exit
        end if
      end do
      
          
    case (psb_asort_down_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (abs(heap%keys(i)) > abs(heap%keys(i2))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp 
          i             = i2
        else
          exit
        end if
      end do
      

    case default
      write(0,*) 'Invalid direction in heap ',heap%dir
    end select

    return
  end subroutine psb_insert_double_idx_heap

  subroutine psb_double_idx_heap_get_first(key,index,heap,info)
    implicit none 

    type(psb_double_idx_heap), intent(inout) :: heap
    integer, intent(out)              :: index,info
    real(kind(1.d0)), intent(out)     :: key
    
    integer                           :: i, i2, last,j,itemp
    real(kind(1.d0))                  :: temp

    
    info = 0
    if (heap%last <= 0) then 
      key  = 0
      info = -1
      return
    endif

    key          = heap%keys(1)
    index        = heap%idxs(1)
    heap%keys(1) = heap%keys(heap%last)
    heap%idxs(1) = heap%idxs(heap%last)
    heap%last    = heap%last - 1
    last         = heap%last

    select case(heap%dir)
    case (psb_sort_up_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (heap%keys(2*i) < heap%keys(2*i+1)) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (heap%keys(i) > heap%keys(j)) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do
      
          
    case (psb_sort_down_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (heap%keys(2*i) > heap%keys(2*i+1)) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (heap%keys(i) < heap%keys(j)) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do

    case (psb_asort_up_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (abs(heap%keys(2*i)) < abs(heap%keys(2*i+1))) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (abs(heap%keys(i)) > abs(heap%keys(j))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do
      
          
    case (psb_asort_down_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (abs(heap%keys(2*i)) > abs(heap%keys(2*i+1))) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (abs(heap%keys(i)) < abs(heap%keys(j))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do

    case default
      write(0,*) 'Invalid direction in heap ',heap%dir
    end select

    return
  end subroutine psb_double_idx_heap_get_first

  function  psb_howmany_int_idx_heap(heap)
    implicit none 
    type(psb_int_idx_heap), intent(in) :: heap
    integer :: psb_howmany_int_idx_heap
    psb_howmany_int_idx_heap = heap%last
  end function psb_howmany_int_idx_heap

  subroutine psb_init_int_idx_heap(heap,info,dir)
    use psb_realloc_mod
    implicit none 
    type(psb_int_idx_heap), intent(inout) :: heap
    integer, intent(out)            :: info
    integer, intent(in), optional   :: dir
    
    info = 0
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
      write(0,*) 'Invalid direction, defaulting to psb_sort_up_'
      heap%dir = psb_sort_up_
    end select
    
    call psb_ensure_size(psb_heap_resize,heap%keys,info)
    call psb_ensure_size(psb_heap_resize,heap%idxs,info)
    return
  end subroutine psb_init_int_idx_heap

  subroutine psb_dump_int_idx_heap(iout,heap,info)
    implicit none 
    type(psb_int_idx_heap), intent(in) :: heap
    integer, intent(out)           :: info
    integer, intent(in)            :: iout
    
    info = 0
    if (iout < 0) then
      write(0,*) 'Invalid file '
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
      if (heap%last > 0) then 
        write(iout,*) heap%keys(1:heap%last)
        write(iout,*) heap%idxs(1:heap%last)
      end if
    end if
  end subroutine psb_dump_int_idx_heap

  subroutine psb_insert_int_idx_heap(key,index,heap,info)
    use psb_realloc_mod
    implicit none 

    integer, intent(in)                   :: key
    integer, intent(in)                   :: index
    type(psb_int_idx_heap), intent(inout) :: heap
    integer, intent(out)                  :: info
    integer                               :: i, i2, itemp
    integer                               :: temp 
    info = 0
    if (heap%last < 0) then 
      write(0,*) 'Invalid last in heap ',heap%last
      info = heap%last
      return
    endif

    heap%last = heap%last + 1
    call psb_ensure_size(heap%last,heap%keys,info,addsz=psb_heap_resize)
    if (info == 0) &
         & call psb_ensure_size(heap%last,heap%idxs,info,addsz=psb_heap_resize)
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in heap_insert'
      info = -5
      return
    end if
    
    i            = heap%last
    heap%keys(i) = key
    heap%idxs(i) = index

    select case(heap%dir)
    case (psb_sort_up_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (heap%keys(i) < heap%keys(i2)) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp 
          i             = i2
        else
          exit
        end if
      end do
      
          
    case (psb_sort_down_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (heap%keys(i) > heap%keys(i2)) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp 
          i             = i2
        else
          exit
        end if
      end do
      
    case (psb_asort_up_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (abs(heap%keys(i)) < abs(heap%keys(i2))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp
          i             = i2
        else
          exit
        end if
      end do
      
          
    case (psb_asort_down_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (abs(heap%keys(i)) > abs(heap%keys(i2))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp 
          i             = i2
        else
          exit
        end if
      end do
      

    case default
      write(0,*) 'Invalid direction in heap ',heap%dir
    end select

    return
  end subroutine psb_insert_int_idx_heap

  subroutine psb_int_idx_heap_get_first(key,index,heap,info)
    implicit none 

    type(psb_int_idx_heap), intent(inout) :: heap
    integer, intent(out)                  :: index,info
    integer, intent(out)                  :: key
    
    integer                               :: i, i2, last,j,itemp
    integer                               :: temp

    
    info = 0
    if (heap%last <= 0) then 
      key  = 0
      info = -1
      return
    endif

    key          = heap%keys(1)
    index        = heap%idxs(1)
    heap%keys(1) = heap%keys(heap%last)
    heap%idxs(1) = heap%idxs(heap%last)
    heap%last    = heap%last - 1
    last         = heap%last

    select case(heap%dir)
    case (psb_sort_up_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (heap%keys(2*i) < heap%keys(2*i+1)) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (heap%keys(i) > heap%keys(j)) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do
      
          
    case (psb_sort_down_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (heap%keys(2*i) > heap%keys(2*i+1)) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (heap%keys(i) < heap%keys(j)) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do

    case (psb_asort_up_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (abs(heap%keys(2*i)) < abs(heap%keys(2*i+1))) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (abs(heap%keys(i)) > abs(heap%keys(j))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do
      
          
    case (psb_asort_down_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (abs(heap%keys(2*i)) > abs(heap%keys(2*i+1))) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (abs(heap%keys(i)) < abs(heap%keys(j))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do

    case default
      write(0,*) 'Invalid direction in heap ',heap%dir
    end select

    return
  end subroutine psb_int_idx_heap_get_first



  function  psb_howmany_dcomplex_idx_heap(heap)
    implicit none 
    type(psb_dcomplex_idx_heap), intent(in) :: heap
    integer :: psb_howmany_dcomplex_idx_heap
    psb_howmany_dcomplex_idx_heap = heap%last
  end function psb_howmany_dcomplex_idx_heap

  subroutine psb_init_dcomplex_idx_heap(heap,info,dir)
    use psb_realloc_mod
    implicit none 
    type(psb_dcomplex_idx_heap), intent(inout) :: heap
    integer, intent(out)            :: info
    integer, intent(in), optional   :: dir
    
    info = 0
    heap%last=0
    if (present(dir)) then 
      heap%dir = dir
    else
      heap%dir = psb_sort_up_
    endif
    select case(heap%dir) 
!!$    case (psb_sort_up_,psb_sort_down_,psb_asort_up_,psb_asort_down_)
    case (psb_asort_up_,psb_asort_down_)
      ! ok, do nothing
    case default
      write(0,*) 'Invalid direction, defaulting to psb_sort_up_'
      heap%dir = psb_asort_up_
    end select
    
    call psb_ensure_size(psb_heap_resize,heap%keys,info)
    call psb_ensure_size(psb_heap_resize,heap%idxs,info)
    return
  end subroutine psb_init_dcomplex_idx_heap

  subroutine psb_dump_dcomplex_idx_heap(iout,heap,info)
    implicit none 
    type(psb_dcomplex_idx_heap), intent(in) :: heap
    integer, intent(out)           :: info
    integer, intent(in)            :: iout
    
    info = 0
    if (iout < 0) then
      write(0,*) 'Invalid file '
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
      if (heap%last > 0) then 
        write(iout,*) heap%keys(1:heap%last)
        write(iout,*) heap%idxs(1:heap%last)
      end if
    end if
  end subroutine psb_dump_dcomplex_idx_heap

  subroutine psb_insert_dcomplex_idx_heap(key,index,heap,info)
    use psb_realloc_mod
    implicit none 

    complex(kind(1.d0)), intent(in)            :: key
    integer, intent(in)                        :: index
    type(psb_dcomplex_idx_heap), intent(inout) :: heap
    integer, intent(out)                       :: info
    integer                                    :: i, i2, itemp
    complex(kind(1.d0))                        :: temp 
    info = 0
    if (heap%last < 0) then 
      write(0,*) 'Invalid last in heap ',heap%last
      info = heap%last
      return
    endif

    heap%last = heap%last + 1
    call psb_ensure_size(heap%last,heap%keys,info,addsz=psb_heap_resize)
    if (info == 0) &
         & call psb_ensure_size(heap%last,heap%idxs,info,addsz=psb_heap_resize)
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in heap_insert'
      info = -5
      return
    end if
    
    i            = heap%last
    heap%keys(i) = key
    heap%idxs(i) = index

    select case(heap%dir)
!!$    case (psb_sort_up_)
!!$
!!$      do 
!!$        if (i<=1) exit
!!$        i2 = i/2
!!$        if (heap%keys(i) < heap%keys(i2)) then 
!!$          itemp         = heap%idxs(i)
!!$          heap%idxs(i)  = heap%idxs(i2)
!!$          heap%idxs(i2) = itemp
!!$          temp          = heap%keys(i)
!!$          heap%keys(i)  = heap%keys(i2)
!!$          heap%keys(i2) = temp 
!!$          i             = i2
!!$        else
!!$          exit
!!$        end if
!!$      end do
!!$      
!!$          
!!$    case (psb_sort_down_)
!!$
!!$      do 
!!$        if (i<=1) exit
!!$        i2 = i/2
!!$        if (heap%keys(i) > heap%keys(i2)) then 
!!$          itemp         = heap%idxs(i)
!!$          heap%idxs(i)  = heap%idxs(i2)
!!$          heap%idxs(i2) = itemp
!!$          temp          = heap%keys(i)
!!$          heap%keys(i)  = heap%keys(i2)
!!$          heap%keys(i2) = temp 
!!$          i             = i2
!!$        else
!!$          exit
!!$        end if
!!$      end do
      
    case (psb_asort_up_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (abs(heap%keys(i)) < abs(heap%keys(i2))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp
          i             = i2
        else
          exit
        end if
      end do
      
          
    case (psb_asort_down_)

      do 
        if (i<=1) exit
        i2 = i/2
        if (abs(heap%keys(i)) > abs(heap%keys(i2))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(i2)
          heap%idxs(i2) = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(i2)
          heap%keys(i2) = temp 
          i             = i2
        else
          exit
        end if
      end do
      

    case default
      write(0,*) 'Invalid direction in heap ',heap%dir
    end select

    return
  end subroutine psb_insert_dcomplex_idx_heap

  subroutine psb_dcomplex_idx_heap_get_first(key,index,heap,info)
    implicit none 

    type(psb_dcomplex_idx_heap), intent(inout) :: heap
    integer, intent(out)                       :: index,info
    complex(kind(1.d0)), intent(out)           :: key
    
    integer                                    :: i, i2, last,j,itemp
    complex(kind(1.d0))                        :: temp

    
    info = 0
    if (heap%last <= 0) then 
      key  = 0
      info = -1
      return
    endif

    key          = heap%keys(1)
    index        = heap%idxs(1)
    heap%keys(1) = heap%keys(heap%last)
    heap%idxs(1) = heap%idxs(heap%last)
    heap%last    = heap%last - 1
    last         = heap%last

    select case(heap%dir)
!!$    case (psb_sort_up_)
!!$
!!$      i = 1
!!$      do 
!!$        if (i > (last/2)) exit
!!$        if ( (heap%keys(2*i) < heap%keys(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap%keys(i) > heap%keys(j)) then 
!!$          itemp         = heap%idxs(i)
!!$          heap%idxs(i)  = heap%idxs(j)
!!$          heap%idxs(j)  = itemp
!!$          temp          = heap%keys(i)
!!$          heap%keys(i)  = heap%keys(j)
!!$          heap%keys(j)  = temp 
!!$          i             = j 
!!$        else
!!$          exit
!!$        end if
!!$      end do
!!$      
!!$          
!!$    case (psb_sort_down_)
!!$
!!$      i = 1
!!$      do 
!!$        if (i > (last/2)) exit
!!$        if ( (heap%keys(2*i) > heap%keys(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap%keys(i) < heap%keys(j)) then 
!!$          itemp         = heap%idxs(i)
!!$          heap%idxs(i)  = heap%idxs(j)
!!$          heap%idxs(j)  = itemp
!!$          temp          = heap%keys(i)
!!$          heap%keys(i)  = heap%keys(j)
!!$          heap%keys(j)  = temp 
!!$          i             = j 
!!$        else
!!$          exit
!!$        end if
!!$      end do

    case (psb_asort_up_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (abs(heap%keys(2*i)) < abs(heap%keys(2*i+1))) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (abs(heap%keys(i)) > abs(heap%keys(j))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do
      
          
    case (psb_asort_down_)

      i = 1
      do 
        if (i > (last/2)) exit
        if ( (abs(heap%keys(2*i)) > abs(heap%keys(2*i+1))) .or.&
             & (2*i == last)) then 
          j = 2*i
        else
          j = 2*i + 1
        end if

        if (abs(heap%keys(i)) < abs(heap%keys(j))) then 
          itemp         = heap%idxs(i)
          heap%idxs(i)  = heap%idxs(j)
          heap%idxs(j)  = itemp
          temp          = heap%keys(i)
          heap%keys(i)  = heap%keys(j)
          heap%keys(j)  = temp 
          i             = j 
        else
          exit
        end if
      end do

    case default
      write(0,*) 'Invalid direction in heap ',heap%dir
    end select

    return
  end subroutine psb_dcomplex_idx_heap_get_first


end module psb_sort_mod
