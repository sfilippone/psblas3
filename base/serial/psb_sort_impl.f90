!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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


logical function psb_isaperm(n,eip)               
  use psb_sort_mod, psb_protect_name => psb_isaperm
  implicit none

  integer(psb_ipk_), intent(in) :: n                                                      
  integer(psb_ipk_), intent(in) :: eip(n)
  integer(psb_ipk_), allocatable :: ip(:)
  integer(psb_ipk_) :: i,j,m, info


  psb_isaperm = .true.
  if (n <= 0) return
  allocate(ip(n), stat=info) 
  if (info /= psb_success_) return
  !
  !   sanity check first 
  !     
  do i=1, n 
    ip(i) = eip(i)
    if ((ip(i) < 1).or.(ip(i) > n)) then
      write(psb_err_unit,*) 'Out of bounds in isaperm' ,ip(i), n
      psb_isaperm = .false.
      return
    endif
  enddo

  !
  ! now work through the cycles, by marking each successive item as negative.
  ! no cycle should intersect with any other, hence the  >= 1 check. 
  !
  do m = 1, n    
    i = ip(m) 
    if (i < 0) then      
      ip(m) = -i          
    else if (i /= m) then 
      j     = ip(i)               
      ip(i) = -j          
      i     = j
      do while ((j >= 1).and.(j /= m))
        j     = ip(i)               
        ip(i) = -j 
        i     = j               
      enddo
      ip(m) = abs(ip(m))
      if (j /= m) then 
        psb_isaperm = .false.
        goto 9999
      endif
    end if
  enddo
9999 continue 

  return                                                                    
end function psb_isaperm

function  psb_iblsrch(key,n,v) result(ipos)
  use psb_sort_mod, psb_protect_name => psb_iblsrch
  implicit none
  integer(psb_ipk_) :: ipos, key, n
  integer(psb_ipk_) :: v(:)

  integer(psb_ipk_) :: lb, ub, m

  if (n < 5) then 
    ! don't bother with binary search for very
    ! small vectors
    ipos = 0
    do
      if (ipos == n) return
      if (key < v(ipos+1)) return 
      ipos = ipos + 1 
    end do
  else
    lb = 1 
    ub = n
    ipos = -1 
    
    do while (lb <= ub) 
      m = (lb+ub)/2
      if (key==v(m))  then
        ipos = m 
        return
      else if (key < v(m))  then
        ub = m-1
      else 
        lb = m + 1
      end if
    enddo
    if (v(ub) > key) then
!!$      write(0,*) 'Check: ',ub,v(ub),key
      ub = ub - 1 
    end if
    ipos = ub 
  endif
  return
end function psb_iblsrch

function  psb_ibsrch(key,n,v) result(ipos)
  use psb_sort_mod, psb_protect_name => psb_ibsrch
  implicit none
  integer(psb_ipk_) :: ipos, key, n
  integer(psb_ipk_) :: v(:)

  integer(psb_ipk_) :: lb, ub, m

  lb = 1 
  ub = n
  ipos = -1 

  do while (lb.le.ub) 
    m = (lb+ub)/2
    if (key.eq.v(m))  then
      ipos = m 
      lb   = ub + 1
    else if (key < v(m))  then
      ub = m-1
    else 
      lb = m + 1
    end if
  enddo
  return
end function psb_ibsrch

function psb_issrch(key,n,v) result(ipos)
  use psb_sort_mod, psb_protect_name => psb_issrch
  implicit none
  integer(psb_ipk_) :: ipos, key, n
  integer(psb_ipk_) :: v(:)

  integer(psb_ipk_) :: i

  ipos = -1 
  do i=1,n
    if (key.eq.v(i))  then
      ipos = i
      return
    end if
  enddo
  
  return
end function psb_issrch


subroutine imsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => imsort
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(inout)           :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  n = size(x)

  if (present(ix)) then 
    if (size(ix) < n) then 
      ierr(1) = 2; ierr(2) = size(ix); 
      call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
      ierr(1) = 4; ierr(2) = flag_; 
      call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
      goto 9999
    end select

    call imsrx(n,x,ix,dir_,flag_)
  else
    call imsr(n,x,dir_)
  end if

9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine imsort


subroutine smsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => smsort
  use psb_error_mod
  implicit none 
  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  n = size(x)

  if (present(ix)) then 
    if (size(ix) < n) then 
      ierr(1) = 2; ierr(2) = size(ix); 
      call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
      ierr(1) = 4; ierr(2) = flag_; 
      call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
      goto 9999
    end select

    call smsrx(n,x,ix,dir_,flag_)
  else
    call smsr(n,x,dir_)
  end if

9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine smsort

subroutine dmsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => dmsort
  use psb_error_mod
  implicit none 
  real(psb_dpk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  n = size(x)

  if (present(ix)) then 
    if (size(ix) < n) then 
      ierr(1) = 2; ierr(2) = size(ix); 
      call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
      ierr(1) = 4; ierr(2) = flag_; 
      call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
      goto 9999
    end select

    call dmsrx(n,x,ix,dir_,flag_)
  else
    call dmsr(n,x,dir_)
  end if

9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine dmsort

subroutine camsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => camsort
  use psb_error_mod
  implicit none 
  complex(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
  character(len=20)  :: name

  name='psb_msort'
  call psb_erractionsave(err_act)

  if (present(dir)) then 
    dir_ = dir
  else
    dir_= psb_asort_up_
  end if
  select case(dir_) 
  case( psb_asort_up_, psb_asort_down_)
    ! OK keep going
  case default
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  n = size(x)

  if (present(ix)) then 
    if (size(ix) < n) then 
      ierr(1) = 2; ierr(2) = size(ix); 
      call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
      ierr(1) = 4; ierr(2) = flag_; 
      call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
      goto 9999
    end select

    call camsrx(n,x,ix,dir_,flag_)
  else
    call camsr(n,x,dir_)
  end if

9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine camsort

subroutine zamsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => zamsort
  use psb_error_mod
  implicit none 
  complex(psb_dpk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
  character(len=20)  :: name

  name='psb_msort'
  call psb_erractionsave(err_act)

  if (present(dir)) then 
    dir_ = dir
  else
    dir_= psb_asort_up_
  end if
  select case(dir_) 
  case( psb_asort_up_, psb_asort_down_)
    ! OK keep going
  case default
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  n = size(x)

  if (present(ix)) then 
    if (size(ix) < n) then 
      ierr(1) = 2; ierr(2) = size(ix); 
      call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
      ierr(1) = 4; ierr(2) = flag_; 
      call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
      goto 9999
    end select

    call zamsrx(n,x,ix,dir_,flag_)
  else
    call zamsr(n,x,dir_)
  end if

9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine zamsort


subroutine imsort_u(x,nout,dir)
  use psb_sort_mod, psb_protect_name => imsort_u
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(inout)           :: x(:) 
  integer(psb_ipk_), intent(out)             :: nout
  integer(psb_ipk_), optional, intent(in)    :: dir

  integer(psb_ipk_) :: dir_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  n = size(x)

  call imsru(n,x,dir_,nout)


9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine imsort_u


subroutine iqsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => iqsort
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(inout)           :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 4; ierr(2) = flag_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if

      call iasrx(n,x,ix,dir_,flag_)
    else
      call iasr(n,x,dir_)
    end if

  case default
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select



9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine iqsort


subroutine sqsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => sqsort
  use psb_error_mod
  implicit none 
  real(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 4; ierr(2) = flag_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if

      call ssrx(n,x,ix,dir_,flag_)
    else
      call ssr(n,x,dir_)
    end if

  case( psb_asort_up_, psb_asort_down_)
    ! OK keep going
    if (present(ix)) then 
      if (size(ix) < n) then 
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if

      call sasrx(n,x,ix,dir_,flag_)
    else
      call sasr(n,x,dir_)
    end if

  case default
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select



9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine sqsort

subroutine dqsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => dqsort
  use psb_error_mod
  implicit none 
  real(psb_dpk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 4; ierr(2) = flag_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if

      call dasrx(n,x,ix,dir_,flag_)
    else
      call dasr(n,x,dir_)
    end if

  case default
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select



9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine dqsort


subroutine cqsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => cqsort
  use psb_error_mod
  implicit none 
  complex(psb_spk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 4; ierr(2) = flag_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if

      call clsrx(n,x,ix,dir_,flag_)
    else
      call clsr(n,x,dir_)
    end if

  case( psb_alsort_up_, psb_alsort_down_)
    ! OK keep going
    if (present(ix)) then 
      if (size(ix) < n) then 
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if

      call calsrx(n,x,ix,dir_,flag_)
    else
      call calsr(n,x,dir_)
    end if

  case( psb_asort_up_, psb_asort_down_)
    ! OK keep going
    if (present(ix)) then 
      if (size(ix) < n) then 
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if

      call casrx(n,x,ix,dir_,flag_)
    else
      call casr(n,x,dir_)
    end if

  case default
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select



9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine cqsort


subroutine zqsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => zqsort
  use psb_error_mod
  implicit none 
  complex(psb_dpk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, err_act

  integer(psb_ipk_)  :: ierr(5)
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
    ierr(1) = 4; ierr(2) = flag_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
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
        ierr(1) = 2; ierr(2) = size(ix); 
        call psb_errpush(psb_err_input_asize_invalid_i_,name,i_err=ierr)
        goto 9999
      end if

      call zasrx(n,x,ix,dir_,flag_)
    else
      call zasr(n,x,dir_)
    end if

  case default
    ierr(1) = 3; ierr(2) = dir_; 
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select



9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine zqsort




subroutine ihsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => ihsort
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(inout)           :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, i, l, err_act,info
  integer(psb_ipk_) :: key
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
  case(psb_sort_up_,psb_sort_down_,psb_asort_up_,psb_asort_down_) 
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
  ! So,  we invert dir_!
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
      call psi_insert_int_idx_heap(key,index,l,x,ix,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! '
      end if
    end do
    do i=n, 2, -1 
      call psi_int_idx_heap_get_first(key,index,l,x,ix,dir_,info)
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
      call psi_insert_int_heap(key,l,x,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! ',l,i
      end if
    end do
    do i=n, 2, -1 
      call psi_int_heap_get_first(key,l,x,dir_,info)
      if (l /= i-1) then 
        write(psb_err_unit,*) 'Mismatch while pulling out of heap ',l,i
      end if
      x(i)  = key
    end do
  end if


9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine ihsort


subroutine shsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => shsort
  use psb_error_mod
  implicit none 
  real(psb_spk_), intent(inout)    :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, i, l, err_act,info
  real(psb_spk_) :: key
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
  case(psb_sort_up_,psb_sort_down_,psb_asort_up_,psb_asort_down_) 
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
  ! So,  we invert dir_!
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
      call psi_insert_sreal_idx_heap(key,index,l,x,ix,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! '
      end if
    end do
    do i=n, 2, -1 
      call psi_sreal_idx_heap_get_first(key,index,l,x,ix,dir_,info)
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
      call psi_insert_real_heap(key,l,x,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! ',l,i
      end if
    end do
    do i=n, 2, -1 
      call psi_real_heap_get_first(key,l,x,dir_,info)
      if (l /= i-1) then 
        write(psb_err_unit,*) 'Mismatch while pulling out of heap ',l,i
      end if
      x(i)  = key
    end do
  end if


9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine shsort


subroutine dhsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => dhsort
  use psb_error_mod
  implicit none 
  real(psb_dpk_), intent(inout)  :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, i, l, err_act,info
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
  case(psb_sort_up_,psb_sort_down_,psb_asort_up_,psb_asort_down_) 
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
  ! So,  we invert dir_!
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
      call psi_insert_dreal_idx_heap(key,index,l,x,ix,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! '
      end if
    end do
    do i=n, 2, -1 
      call psi_dreal_idx_heap_get_first(key,index,l,x,ix,dir_,info)
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
      call psi_insert_double_heap(key,l,x,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! ',l,i
      end if
    end do
    do i=n, 2, -1 
      call psi_double_heap_get_first(key,l,x,dir_,info)
      if (l /= i-1) then 
        write(psb_err_unit,*) 'Mismatch while pulling out of heap ',l,i
      end if
      x(i)  = key
    end do
  end if


9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine dhsort


subroutine chsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => chsort
  use psb_error_mod
  implicit none 
  complex(psb_spk_), intent(inout) :: x(:) 
  integer(psb_ipk_), optional, intent(in)    :: dir, flag
  integer(psb_ipk_), optional, intent(inout) :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, i, l, err_act,info
  complex(psb_spk_) :: key
  integer(psb_ipk_) :: index

  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name

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
    dir_= psb_asort_up_
  end if

  select case(dir_)
  case(psb_asort_up_,psb_asort_down_) 
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
  ! So,  we invert dir_!
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
      call psi_insert_scomplex_idx_heap(key,index,l,x,ix,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! '
      end if
    end do
    do i=n, 2, -1 
      call psi_scomplex_idx_heap_get_first(key,index,l,x,ix,dir_,info)
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
      call psi_insert_scomplex_heap(key,l,x,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! ',l,i
      end if
    end do
    do i=n, 2, -1 
      call psi_scomplex_heap_get_first(key,l,x,dir_,info)
      if (l /= i-1) then 
        write(psb_err_unit,*) 'Mismatch while pulling out of heap ',l,i
      end if
      x(i)  = key
    end do
  end if


9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine chsort


subroutine zhsort(x,ix,dir,flag)
  use psb_sort_mod, psb_protect_name => zhsort
  use psb_error_mod
  implicit none 
  complex(psb_dpk_), intent(inout) :: x(:) 
  integer(psb_ipk_), optional, intent(in)      :: dir, flag
  integer(psb_ipk_), optional, intent(inout)   :: ix(:)

  integer(psb_ipk_) :: dir_, flag_, n, i, l, err_act,info
  complex(psb_dpk_) :: key
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
    dir_= psb_asort_up_
  end if

  select case(dir_)
  case(psb_asort_up_,psb_asort_down_) 
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
  ! So,  we invert dir_!
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
      call psi_insert_dcomplex_idx_heap(key,index,l,x,ix,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! '
      end if
    end do
    do i=n, 2, -1 
      call psi_dcomplex_idx_heap_get_first(key,index,l,x,ix,dir_,info)
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
      call psi_insert_dcomplex_heap(key,l,x,dir_,info)
      if (l /= i) then 
        write(psb_err_unit,*) 'Mismatch while heapifying ! ',l,i
      end if
    end do
    do i=n, 2, -1 
      call psi_dcomplex_heap_get_first(key,l,x,dir_,info)
      if (l /= i-1) then 
        write(psb_err_unit,*) 'Mismatch while pulling out of heap ',l,i
      end if
      x(i)  = key
    end do
  end if


9999 continue 
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine zhsort


function  psb_howmany_int_heap(heap)
  use psb_sort_mod, psb_protect_name => psb_howmany_int_heap
  implicit none 
  type(psb_int_heap), intent(in) :: heap
  integer(psb_ipk_) :: psb_howmany_int_heap
  psb_howmany_int_heap = heap%last
end function psb_howmany_int_heap

subroutine psb_init_int_heap(heap,info,dir)
  use psb_sort_mod, psb_protect_name => psb_init_int_heap
  use psb_realloc_mod
  implicit none 
  type(psb_int_heap), intent(inout) :: heap
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
end subroutine psb_init_int_heap

subroutine psb_dump_int_heap(iout,heap,info)
  use psb_sort_mod, psb_protect_name => psb_dump_int_heap
  implicit none 
  type(psb_int_heap), intent(in) :: heap
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), intent(in)            :: iout

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
    if (heap%last > 0) then 
      write(iout,*) heap%keys(1:heap%last)
    end if
  end if
end subroutine psb_dump_int_heap

subroutine psb_insert_int_heap(key,heap,info)
  use psb_sort_mod, psb_protect_name => psb_insert_int_heap
  use psb_realloc_mod
  implicit none 

  integer(psb_ipk_), intent(in)               :: key
  type(psb_int_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)              :: info

  info = psb_success_
  if (heap%last < 0) then 
    write(psb_err_unit,*) 'Invalid last in heap ',heap%last
    info = heap%last
    return
  endif

  heap%last = heap%last 
  call psb_ensure_size(heap%last+1,heap%keys,info,addsz=psb_heap_resize)
  if (info /= psb_success_) then 
    write(psb_err_unit,*) 'Memory allocation failure in heap_insert'
    info = -5
    return
  end if
  call  psi_insert_int_heap(key,heap%last,heap%keys,heap%dir,info)

  return
end subroutine psb_insert_int_heap


subroutine psb_int_heap_get_first(key,heap,info)
  use psb_sort_mod, psb_protect_name => psb_int_heap_get_first
  implicit none 

  type(psb_int_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)              :: key,info

  info = psb_success_

  call psi_int_heap_get_first(key,heap%last,heap%keys,heap%dir,info)

  return
end subroutine psb_int_heap_get_first


function  psb_howmany_sreal_idx_heap(heap)
  use psb_sort_mod, psb_protect_name => psb_howmany_sreal_idx_heap
  implicit none 
  type(psb_sreal_idx_heap), intent(in) :: heap
  integer(psb_ipk_) :: psb_howmany_sreal_idx_heap
  psb_howmany_sreal_idx_heap = heap%last
end function psb_howmany_sreal_idx_heap

subroutine psb_init_sreal_idx_heap(heap,info,dir)
  use psb_sort_mod, psb_protect_name => psb_init_sreal_idx_heap
  use psb_realloc_mod
  implicit none 
  type(psb_sreal_idx_heap), intent(inout) :: heap
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
end subroutine psb_init_sreal_idx_heap

subroutine psb_dump_sreal_idx_heap(iout,heap,info)
  use psb_sort_mod, psb_protect_name => psb_dump_sreal_idx_heap
  implicit none 
  type(psb_sreal_idx_heap), intent(in) :: heap
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), intent(in)            :: iout

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
    if (heap%last > 0) then 
      write(iout,*) heap%keys(1:heap%last)
      write(iout,*) heap%idxs(1:heap%last)
    end if
  end if
end subroutine psb_dump_sreal_idx_heap

subroutine psb_insert_sreal_idx_heap(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_insert_sreal_idx_heap
  use psb_realloc_mod
  implicit none 

  real(psb_spk_), intent(in)      :: key
  integer(psb_ipk_), intent(in)               :: index
  type(psb_sreal_idx_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)              :: info

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

  call psi_insert_sreal_idx_heap(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_insert_sreal_idx_heap

subroutine psb_sreal_idx_heap_get_first(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_sreal_idx_heap_get_first
  implicit none 

  type(psb_sreal_idx_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)              :: index,info
  real(psb_spk_), intent(out)     :: key

  info = psb_success_

  call psi_sreal_idx_heap_get_first(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_sreal_idx_heap_get_first


function  psb_howmany_dreal_idx_heap(heap)
  use psb_sort_mod, psb_protect_name => psb_howmany_dreal_idx_heap
  implicit none 
  type(psb_dreal_idx_heap), intent(in) :: heap
  integer(psb_ipk_) :: psb_howmany_dreal_idx_heap
  psb_howmany_dreal_idx_heap = heap%last
end function psb_howmany_dreal_idx_heap

subroutine psb_init_dreal_idx_heap(heap,info,dir)
  use psb_sort_mod, psb_protect_name => psb_init_dreal_idx_heap
  use psb_realloc_mod
  implicit none 
  type(psb_dreal_idx_heap), intent(inout) :: heap
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
end subroutine psb_init_dreal_idx_heap

subroutine psb_dump_dreal_idx_heap(iout,heap,info)
  use psb_sort_mod, psb_protect_name => psb_dump_dreal_idx_heap
  implicit none 
  type(psb_dreal_idx_heap), intent(in) :: heap
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), intent(in)            :: iout

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
    if (heap%last > 0) then 
      write(iout,*) heap%keys(1:heap%last)
      write(iout,*) heap%idxs(1:heap%last)
    end if
  end if
end subroutine psb_dump_dreal_idx_heap

subroutine psb_insert_dreal_idx_heap(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_insert_dreal_idx_heap
  use psb_realloc_mod
  implicit none 

  real(psb_dpk_), intent(in)      :: key
  integer(psb_ipk_), intent(in)               :: index
  type(psb_dreal_idx_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)              :: info

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

  call psi_insert_dreal_idx_heap(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_insert_dreal_idx_heap

subroutine psb_dreal_idx_heap_get_first(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_dreal_idx_heap_get_first
  implicit none 

  type(psb_dreal_idx_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)              :: index,info
  real(psb_dpk_), intent(out)     :: key

  info = psb_success_

  call psi_dreal_idx_heap_get_first(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_dreal_idx_heap_get_first

function  psb_howmany_int_idx_heap(heap)
  use psb_sort_mod, psb_protect_name => psb_howmany_int_idx_heap
  implicit none 
  type(psb_int_idx_heap), intent(in) :: heap
  integer(psb_ipk_) :: psb_howmany_int_idx_heap
  psb_howmany_int_idx_heap = heap%last
end function psb_howmany_int_idx_heap

subroutine psb_init_int_idx_heap(heap,info,dir)
  use psb_sort_mod, psb_protect_name => psb_init_int_idx_heap
  use psb_realloc_mod
  implicit none 
  type(psb_int_idx_heap), intent(inout) :: heap
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
end subroutine psb_init_int_idx_heap

subroutine psb_dump_int_idx_heap(iout,heap,info)
  use psb_sort_mod, psb_protect_name => psb_dump_int_idx_heap
  implicit none 
  type(psb_int_idx_heap), intent(in) :: heap
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), intent(in)            :: iout

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
    if (heap%last > 0) then 
      write(iout,*) heap%keys(1:heap%last)
      write(iout,*) heap%idxs(1:heap%last)
    end if
  end if
end subroutine psb_dump_int_idx_heap

subroutine psb_insert_int_idx_heap(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_insert_int_idx_heap
  use psb_realloc_mod
  implicit none 

  integer(psb_ipk_), intent(in)                   :: key
  integer(psb_ipk_), intent(in)                   :: index
  type(psb_int_idx_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)                  :: info

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

  call psi_insert_int_idx_heap(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_insert_int_idx_heap

subroutine psb_int_idx_heap_get_first(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_int_idx_heap_get_first
  implicit none 

  type(psb_int_idx_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)                  :: index,info
  integer(psb_ipk_), intent(out)                  :: key

  info = psb_success_

  call psi_int_idx_heap_get_first(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_int_idx_heap_get_first



function  psb_howmany_scomplex_idx_heap(heap)
  use psb_sort_mod, psb_protect_name => psb_howmany_scomplex_idx_heap
  implicit none 
  type(psb_scomplex_idx_heap), intent(in) :: heap
  integer(psb_ipk_) :: psb_howmany_scomplex_idx_heap
  psb_howmany_scomplex_idx_heap = heap%last
end function psb_howmany_scomplex_idx_heap

subroutine psb_init_scomplex_idx_heap(heap,info,dir)
  use psb_sort_mod, psb_protect_name => psb_init_scomplex_idx_heap
  use psb_realloc_mod
  implicit none 
  type(psb_scomplex_idx_heap), intent(inout) :: heap
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
!!$    case (psb_sort_up_,psb_sort_down_,psb_asort_up_,psb_asort_down_)
  case (psb_asort_up_,psb_asort_down_)
    ! ok, do nothing
  case default
    write(psb_err_unit,*) 'Invalid direction, defaulting to psb_sort_up_'
    heap%dir = psb_asort_up_
  end select

  call psb_ensure_size(psb_heap_resize,heap%keys,info)
  call psb_ensure_size(psb_heap_resize,heap%idxs,info)
  return
end subroutine psb_init_scomplex_idx_heap

subroutine psb_dump_scomplex_idx_heap(iout,heap,info)
  use psb_sort_mod, psb_protect_name => psb_dump_scomplex_idx_heap
  implicit none 
  type(psb_scomplex_idx_heap), intent(in) :: heap
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), intent(in)            :: iout

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
    if (heap%last > 0) then 
      write(iout,*) heap%keys(1:heap%last)
      write(iout,*) heap%idxs(1:heap%last)
    end if
  end if
end subroutine psb_dump_scomplex_idx_heap

subroutine psb_insert_scomplex_idx_heap(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_insert_scomplex_idx_heap
  use psb_realloc_mod
  implicit none 

  complex(psb_spk_), intent(in)              :: key
  integer(psb_ipk_), intent(in)                        :: index
  type(psb_scomplex_idx_heap), intent(inout) :: heap
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
  call psi_insert_scomplex_idx_heap(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_insert_scomplex_idx_heap

subroutine psb_scomplex_idx_heap_get_first(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_scomplex_idx_heap_get_first
  implicit none 

  type(psb_scomplex_idx_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)                       :: index,info
  complex(psb_spk_), intent(out)           :: key


  info = psb_success_

  call psi_scomplex_idx_heap_get_first(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_scomplex_idx_heap_get_first



function  psb_howmany_dcomplex_idx_heap(heap)
  use psb_sort_mod, psb_protect_name => psb_howmany_dcomplex_idx_heap
  implicit none 
  type(psb_dcomplex_idx_heap), intent(in) :: heap
  integer(psb_ipk_) :: psb_howmany_dcomplex_idx_heap
  psb_howmany_dcomplex_idx_heap = heap%last
end function psb_howmany_dcomplex_idx_heap

subroutine psb_init_dcomplex_idx_heap(heap,info,dir)
  use psb_sort_mod, psb_protect_name => psb_init_dcomplex_idx_heap
  use psb_realloc_mod
  implicit none 
  type(psb_dcomplex_idx_heap), intent(inout) :: heap
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
!!$    case (psb_sort_up_,psb_sort_down_,psb_asort_up_,psb_asort_down_)
  case (psb_asort_up_,psb_asort_down_)
    ! ok, do nothing
  case default
    write(psb_err_unit,*) 'Invalid direction, defaulting to psb_sort_up_'
    heap%dir = psb_asort_up_
  end select

  call psb_ensure_size(psb_heap_resize,heap%keys,info)
  call psb_ensure_size(psb_heap_resize,heap%idxs,info)
  return
end subroutine psb_init_dcomplex_idx_heap

subroutine psb_dump_dcomplex_idx_heap(iout,heap,info)
  use psb_sort_mod, psb_protect_name => psb_dump_dcomplex_idx_heap
  implicit none 
  type(psb_dcomplex_idx_heap), intent(in) :: heap
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), intent(in)            :: iout

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
    if (heap%last > 0) then 
      write(iout,*) heap%keys(1:heap%last)
      write(iout,*) heap%idxs(1:heap%last)
    end if
  end if
end subroutine psb_dump_dcomplex_idx_heap

subroutine psb_insert_dcomplex_idx_heap(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_insert_dcomplex_idx_heap
  use psb_realloc_mod
  implicit none 

  complex(psb_dpk_), intent(in)            :: key
  integer(psb_ipk_), intent(in)                        :: index
  type(psb_dcomplex_idx_heap), intent(inout) :: heap
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
  call psi_insert_dcomplex_idx_heap(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_insert_dcomplex_idx_heap

subroutine psb_dcomplex_idx_heap_get_first(key,index,heap,info)
  use psb_sort_mod, psb_protect_name => psb_dcomplex_idx_heap_get_first
  implicit none 

  type(psb_dcomplex_idx_heap), intent(inout) :: heap
  integer(psb_ipk_), intent(out)                       :: index,info
  complex(psb_dpk_), intent(out)           :: key


  info = psb_success_

  call psi_dcomplex_idx_heap_get_first(key,index,&
       & heap%last,heap%keys,heap%idxs,heap%dir,info)

  return
end subroutine psb_dcomplex_idx_heap_get_first



!
! These are packaged so that they can be used to implement 
! a heapsort, should the need arise
!


subroutine psi_insert_int_heap(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_int_heap
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
  integer(psb_ipk_) :: i, i2
  integer(psb_ipk_) :: temp

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
end subroutine psi_insert_int_heap



!
!
!   Programming note:
!   In the implementation of the various get_first heap function
!   we often have the following code (or similar)
!
!      if ( (abs(heap(2*i)) < abs(heap(2*i+1))) .or.&
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

subroutine psi_int_heap_get_first(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_int_heap_get_first
  implicit none 

  integer(psb_ipk_), intent(inout)  :: key,last
  integer(psb_ipk_), intent(in)     :: dir
  integer(psb_ipk_), intent(inout)  :: heap(:)
  integer(psb_ipk_), intent(out)    :: info

  integer(psb_ipk_) :: i, j
  integer(psb_ipk_) :: temp


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
end subroutine psi_int_heap_get_first



subroutine psi_insert_real_heap(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_real_heap
  implicit none 

  !  
  ! Input: 
  !   key:  the new value
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   dir:  sorting direction

  real(psb_spk_), intent(in)    :: key
  integer(psb_ipk_), intent(in)           :: dir
  real(psb_spk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(inout)        :: last
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_) :: i, i2
  real(psb_spk_)                :: temp

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
end subroutine psi_insert_real_heap


subroutine psi_real_heap_get_first(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_real_heap_get_first
  implicit none 

  real(psb_spk_), intent(inout) :: key
  integer(psb_ipk_), intent(inout)        :: last
  integer(psb_ipk_), intent(in)           :: dir
  real(psb_spk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)          :: info

  integer(psb_ipk_) :: i, j
  real(psb_spk_)        :: temp


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
end subroutine psi_real_heap_get_first


subroutine psi_insert_double_heap(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_double_heap
  implicit none 

  !  
  ! Input: 
  !   key:  the new value
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   dir:  sorting direction

  real(psb_dpk_), intent(in)    :: key
  integer(psb_ipk_), intent(in)             :: dir
  real(psb_dpk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(inout)          :: last
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_) :: i, i2
  real(psb_dpk_)                :: temp

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
end subroutine psi_insert_double_heap


subroutine psi_double_heap_get_first(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_double_heap_get_first
  implicit none 

  real(psb_dpk_), intent(inout) :: key
  integer(psb_ipk_), intent(inout)          :: last
  integer(psb_ipk_), intent(in)             :: dir
  real(psb_dpk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: i, j
  real(psb_dpk_)        :: temp


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
end subroutine psi_double_heap_get_first




subroutine psi_insert_scomplex_heap(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_scomplex_heap
  implicit none 

  !  
  ! Input: 
  !   key:  the new value
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   dir:  sorting direction

  complex(psb_spk_), intent(in)    :: key
  integer(psb_ipk_), intent(in)              :: dir
  complex(psb_spk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(inout)           :: last
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_) :: i, i2
  complex(psb_spk_)                :: temp

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
!!$    case (psb_sort_up_)
!!$
!!$      do 
!!$        if (i<=1) exit
!!$        i2 = i/2
!!$        if (heap(i) < heap(i2)) then 
!!$          temp     = heap(i)
!!$          heap(i)  = heap(i2)
!!$          heap(i2) = temp
!!$          i        = i2
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
!!$        if (heap(i) > heap(i2)) then 
!!$          temp     = heap(i)
!!$          heap(i)  = heap(i2)
!!$          heap(i2) = temp
!!$          i        = i2
!!$        else
!!$          exit
!!$        end if
!!$      end do

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
end subroutine psi_insert_scomplex_heap


subroutine psi_scomplex_heap_get_first(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_scomplex_heap_get_first
  implicit none 

  complex(psb_spk_), intent(inout) :: key
  integer(psb_ipk_), intent(inout)           :: last
  integer(psb_ipk_), intent(in)              :: dir
  complex(psb_spk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)             :: info

  integer(psb_ipk_) :: i, j
  complex(psb_spk_)        :: temp


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
!!$    case (psb_sort_up_)
!!$
!!$      i = 1
!!$      do 
!!$        if (i > (last/2)) exit
!!$        if ( (heap(2*i) < heap(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap(i) > heap(j)) then 
!!$          temp     = heap(i)
!!$          heap(i)  = heap(j)
!!$          heap(j)  = temp
!!$          i        = j 
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
!!$        if ( (heap(2*i) > heap(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap(i) < heap(j)) then 
!!$          temp     = heap(i)
!!$          heap(i)  = heap(j)
!!$          heap(j)  = temp
!!$          i        = j 
!!$        else
!!$          exit
!!$        end if
!!$      end do

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
end subroutine psi_scomplex_heap_get_first


subroutine psi_insert_dcomplex_heap(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_dcomplex_heap
  implicit none 

  !  
  ! Input: 
  !   key:  the new value
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   dir:  sorting direction

  complex(psb_dpk_), intent(in)    :: key
  integer(psb_ipk_), intent(in)                :: dir
  complex(psb_dpk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(inout)             :: last
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_) :: i, i2
  complex(psb_dpk_)                :: temp

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
!!$    case (psb_sort_up_)
!!$
!!$      do 
!!$        if (i<=1) exit
!!$        i2 = i/2
!!$        if (heap(i) < heap(i2)) then 
!!$          temp     = heap(i)
!!$          heap(i)  = heap(i2)
!!$          heap(i2) = temp
!!$          i        = i2
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
!!$        if (heap(i) > heap(i2)) then 
!!$          temp     = heap(i)
!!$          heap(i)  = heap(i2)
!!$          heap(i2) = temp
!!$          i        = i2
!!$        else
!!$          exit
!!$        end if
!!$      end do

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
end subroutine psi_insert_dcomplex_heap


subroutine psi_dcomplex_heap_get_first(key,last,heap,dir,info)
  use psb_sort_mod, psb_protect_name => psi_dcomplex_heap_get_first
  implicit none 

  complex(psb_dpk_), intent(inout) :: key
  integer(psb_ipk_), intent(inout)             :: last
  integer(psb_ipk_), intent(in)                :: dir
  complex(psb_dpk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)               :: info

  integer(psb_ipk_) :: i, j
  complex(psb_dpk_)        :: temp


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
!!$    case (psb_sort_up_)
!!$
!!$      i = 1
!!$      do 
!!$        if (i > (last/2)) exit
!!$        if ( (heap(2*i) < heap(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap(i) > heap(j)) then 
!!$          temp     = heap(i)
!!$          heap(i)  = heap(j)
!!$          heap(j)  = temp
!!$          i        = j 
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
!!$        if ( (heap(2*i) > heap(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap(i) < heap(j)) then 
!!$          temp     = heap(i)
!!$          heap(i)  = heap(j)
!!$          heap(j)  = temp
!!$          i        = j 
!!$        else
!!$          exit
!!$        end if
!!$      end do

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
end subroutine psi_dcomplex_heap_get_first




subroutine psi_insert_int_idx_heap(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_int_idx_heap

  implicit none 
  !  
  ! Input: 
  !   key:  the new value
  !   index: the new index
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   idxs: the indices
  !   dir:  sorting direction

  integer(psb_ipk_), intent(in)     :: key
  integer(psb_ipk_), intent(in)     :: index,dir
  integer(psb_ipk_), intent(inout)  :: heap(:),last
  integer(psb_ipk_), intent(inout)  :: idxs(:)
  integer(psb_ipk_), intent(out)    :: info
  integer(psb_ipk_) :: i, i2, itemp
  integer(psb_ipk_) :: temp 

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
end subroutine psi_insert_int_idx_heap

subroutine psi_int_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_int_idx_heap_get_first
  implicit none 

  integer(psb_ipk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)   :: index,info
  integer(psb_ipk_), intent(inout) :: last,idxs(:)
  integer(psb_ipk_), intent(in)    :: dir
  integer(psb_ipk_), intent(out)   :: key

  integer(psb_ipk_) :: i, j,itemp
  integer(psb_ipk_) :: temp

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
end subroutine psi_int_idx_heap_get_first

subroutine psi_insert_sreal_idx_heap(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_sreal_idx_heap

  implicit none 
  !  
  ! Input: 
  !   key:  the new value
  !   index: the new index
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   idxs: the indices
  !   dir:  sorting direction

  real(psb_spk_), intent(in)     :: key
  integer(psb_ipk_), intent(in)            :: index,dir
  real(psb_spk_), intent(inout)  :: heap(:)
  integer(psb_ipk_), intent(inout)         :: idxs(:),last
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_) :: i, i2, itemp
  real(psb_spk_)                 :: temp 

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
end subroutine psi_insert_sreal_idx_heap

subroutine psi_sreal_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_sreal_idx_heap_get_first
  implicit none 

  real(psb_spk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)          :: index,info
  integer(psb_ipk_), intent(inout)        :: last,idxs(:)
  integer(psb_ipk_), intent(in)           :: dir
  real(psb_spk_), intent(out)   :: key

  integer(psb_ipk_) :: i, j,itemp
  real(psb_spk_)                :: temp

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
end subroutine psi_sreal_idx_heap_get_first


subroutine psi_insert_dreal_idx_heap(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_dreal_idx_heap

  implicit none 
  !  
  ! Input: 
  !   key:  the new value
  !   index: the new index
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   idxs: the indices
  !   dir:  sorting direction

  real(psb_dpk_), intent(in)     :: key
  integer(psb_ipk_), intent(in)              :: index,dir
  real(psb_dpk_), intent(inout)  :: heap(:)
  integer(psb_ipk_), intent(inout)           :: idxs(:),last
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_) :: i, i2, itemp
  real(psb_dpk_)                 :: temp 

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
end subroutine psi_insert_dreal_idx_heap

subroutine psi_dreal_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_dreal_idx_heap_get_first
  implicit none 

  real(psb_dpk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)            :: index,info
  integer(psb_ipk_), intent(inout)          :: last,idxs(:)
  integer(psb_ipk_), intent(in)             :: dir
  real(psb_dpk_), intent(out)   :: key

  integer(psb_ipk_) :: i, j,itemp
  real(psb_dpk_)                :: temp

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
end subroutine psi_dreal_idx_heap_get_first


subroutine psi_insert_scomplex_idx_heap(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_scomplex_idx_heap

  implicit none 
  !  
  ! Input: 
  !   key:  the new value
  !   index: the new index
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   idxs: the indices
  !   dir:  sorting direction

  complex(psb_spk_), intent(in)    :: key
  integer(psb_ipk_), intent(in)              :: index,dir
  complex(psb_spk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(inout)           :: idxs(:),last
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_) :: i, i2, itemp
  complex(psb_spk_)                :: temp 

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
!!$    case (psb_sort_up_)
!!$
!!$      do 
!!$        if (i<=1) exit
!!$        i2 = i/2
!!$        if (heap(i) < heap(i2)) then 
!!$          itemp    = idxs(i)
!!$          idxs(i)  = idxs(i2)
!!$          idxs(i2) = itemp
!!$          temp     = heap(i)
!!$          heap(i)  = heap(i2)
!!$          heap(i2) = temp 
!!$          i        = i2
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
!!$        if (heap(i) > heap(i2)) then 
!!$          itemp    = idxs(i)
!!$          idxs(i)  = idxs(i2)
!!$          idxs(i2) = itemp
!!$          temp     = heap(i)
!!$          heap(i)  = heap(i2)
!!$          heap(i2) = temp 
!!$          i        = i2
!!$        else
!!$          exit
!!$        end if
!!$      end do

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
end subroutine psi_insert_scomplex_idx_heap

subroutine psi_scomplex_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_scomplex_idx_heap_get_first
  implicit none 

  complex(psb_spk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)             :: index,info
  integer(psb_ipk_), intent(inout)           :: last,idxs(:)
  integer(psb_ipk_), intent(in)              :: dir
  complex(psb_spk_), intent(out)   :: key

  integer(psb_ipk_) :: i, j, itemp
  complex(psb_spk_)                :: temp

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
!!$    case (psb_sort_up_)
!!$
!!$      i = 1
!!$      do 
!!$        if (i > (last/2)) exit
!!$        if ( (heap(2*i) < heap(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap(i) > heap(j)) then 
!!$          itemp    = idxs(i)
!!$          idxs(i)  = idxs(j)
!!$          idxs(j)  = itemp
!!$          temp     = heap(i)
!!$          heap(i)  = heap(j)
!!$          heap(j)  = temp 
!!$          i        = j 
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
!!$        if ( (heap(2*i) > heap(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap(i) < heap(j)) then 
!!$          itemp    = idxs(i)
!!$          idxs(i)  = idxs(j)
!!$          idxs(j)  = itemp
!!$          temp     = heap(i)
!!$          heap(i)  = heap(j)
!!$          heap(j)  = temp 
!!$          i        = j 
!!$        else
!!$          exit
!!$        end if
!!$      end do

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
end subroutine psi_scomplex_idx_heap_get_first


subroutine psi_insert_dcomplex_idx_heap(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_insert_dcomplex_idx_heap

  implicit none 
  !  
  ! Input: 
  !   key:  the new value
  !   index: the new index
  !   last: pointer to the last occupied element in heap
  !   heap: the heap
  !   idxs: the indices
  !   dir:  sorting direction

  complex(psb_dpk_), intent(in)    :: key
  integer(psb_ipk_), intent(in)                :: index,dir
  complex(psb_dpk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(inout)             :: idxs(:),last
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_) :: i, i2, itemp
  complex(psb_dpk_)                :: temp 

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
!!$    case (psb_sort_up_)
!!$
!!$      do 
!!$        if (i<=1) exit
!!$        i2 = i/2
!!$        if (heap(i) < heap(i2)) then 
!!$          itemp    = idxs(i)
!!$          idxs(i)  = idxs(i2)
!!$          idxs(i2) = itemp
!!$          temp     = heap(i)
!!$          heap(i)  = heap(i2)
!!$          heap(i2) = temp 
!!$          i        = i2
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
!!$        if (heap(i) > heap(i2)) then 
!!$          itemp    = idxs(i)
!!$          idxs(i)  = idxs(i2)
!!$          idxs(i2) = itemp
!!$          temp     = heap(i)
!!$          heap(i)  = heap(i2)
!!$          heap(i2) = temp 
!!$          i        = i2
!!$        else
!!$          exit
!!$        end if
!!$      end do

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
end subroutine psi_insert_dcomplex_idx_heap

subroutine psi_dcomplex_idx_heap_get_first(key,index,last,heap,idxs,dir,info)
  use psb_sort_mod, psb_protect_name => psi_dcomplex_idx_heap_get_first
  implicit none 

  complex(psb_dpk_), intent(inout) :: heap(:)
  integer(psb_ipk_), intent(out)               :: index,info
  integer(psb_ipk_), intent(inout)             :: last,idxs(:)
  integer(psb_ipk_), intent(in)                :: dir
  complex(psb_dpk_), intent(out)   :: key

  integer(psb_ipk_) :: i, j, itemp
  complex(psb_dpk_)                :: temp

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
!!$    case (psb_sort_up_)
!!$
!!$      i = 1
!!$      do 
!!$        if (i > (last/2)) exit
!!$        if ( (heap(2*i) < heap(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap(i) > heap(j)) then 
!!$          itemp    = idxs(i)
!!$          idxs(i)  = idxs(j)
!!$          idxs(j)  = itemp
!!$          temp     = heap(i)
!!$          heap(i)  = heap(j)
!!$          heap(j)  = temp 
!!$          i        = j 
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
!!$        if ( (heap(2*i) > heap(2*i+1)) .or.&
!!$             & (2*i == last)) then 
!!$          j = 2*i
!!$        else
!!$          j = 2*i + 1
!!$        end if
!!$
!!$        if (heap(i) < heap(j)) then 
!!$          itemp    = idxs(i)
!!$          idxs(i)  = idxs(j)
!!$          idxs(j)  = itemp
!!$          temp     = heap(i)
!!$          heap(i)  = heap(j)
!!$          heap(j)  = temp 
!!$          i        = j 
!!$        else
!!$          exit
!!$        end if
!!$      end do

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
end subroutine psi_dcomplex_idx_heap_get_first



