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
!         software without specific prior written permission.
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
subroutine psi_e_exscanv(n, x, info, shift)
  use psi_e_serial_mod, psb_protect_name => psi_e_exscanv
  use psb_const_mod
  use psb_error_mod
#if defined(PSB_OPENMP)
  use omp_lib
#endif
  implicit none
  integer(psb_ipk_), intent(in)   :: n
  integer(psb_epk_), intent(inout)   :: x(:)
  integer(psb_ipk_), intent(out)  :: info
  integer(psb_epk_), intent(in), optional :: shift
  
  integer(psb_epk_) :: shift_, tp, ts
  integer(psb_ipk_) :: i
  logical is_nested, is_parallel

  info = psb_success_
  
  shift_ = ezero
  if(present(shift)) shift_ = shift
    
#if defined(PSB_OPENMP)
  is_parallel = omp_in_parallel()
  if(is_parallel) then 
    call inner_e_exscan()
  else
    !$OMP PARALLEL default(shared) 
    call inner_e_exscan()
    !$OMP END PARALLEL
  end if
#else
  tp = shift_
  do i = 1, n
    ts = x(i)
    x(i) = tp
    tp = tp + ts
  end do

#endif
#if defined(PSB_OPENMP)
contains
  subroutine inner_e_exscan()
    ! Note: all these variables are private, but SUMB should *really* be
    ! a pointer. The semantics of COPYPRIVATE is that the POINTER is copied
    ! so effectively we are recovering a SHARED SUMB which is what
    ! we need in this case. If it was an ALLOCATABLE, then it would be the contents
    ! that would get copied, and the SHARED effect would  no longer be there.
    ! Simple parallel version of EXSCAN
    integer(psb_ipk_)       :: i, ithread, nthreads, idxstart, idxend, wrk
    integer(psb_epk_), pointer :: sumb(:)
    integer(psb_epk_)          :: tp, ts

    nthreads = omp_get_num_threads()
    ithread = omp_get_thread_num()
    !$OMP SINGLE
    allocate(sumb(nthreads+1))
    sumb(:) = 0
    !$OMP END SINGLE COPYPRIVATE(sumb)

    wrk = (n)/nthreads
    if(ithread < MOD((n), nthreads)) then
      wrk = wrk + 1
      idxstart = ithread*wrk + 1
    else
      idxstart = ithread*wrk + MOD((n), nthreads) + 1
    end if

    idxend = min(idxstart + wrk - 1, n )
    tp = ezero
    if(idxstart<=idxend) then
      do i = idxstart, idxend
        ts = x(i)
        x(i) = tp 
        tp = tp + ts 
      end do
    end if
    sumb(ithread+2) = tp 
    !$OMP BARRIER
    
    !$OMP SINGLE
    do i = 2, nthreads+1
      sumb(i) = sumb(i) + sumb(i-1)
    end do
    !$OMP END SINGLE      

    !$OMP BARRIER

    !$OMP DO SCHEDULE(STATIC)
    do i = 1, n
      x(i) = x(i) + sumb(ithread+1) + shift_ 
    end do
    !$OMP END DO
    !$OMP SINGLE
    deallocate(sumb)
    !$OMP END SINGLE
  end subroutine inner_e_exscan
#endif
end subroutine psi_e_exscanv

subroutine psb_m_egelp(trans, iperm, x, info)
  use psb_serial_mod, psb_protect_name => psb_m_egelp
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_epk_), intent(inout)   :: x(:, :)
  integer(psb_mpk_), intent(in)   :: iperm(:)
  integer(psb_ipk_), intent(out)  :: info
  character, intent(in)           :: trans

  ! local variables
  integer(psb_epk_), allocatable     :: temp(:)
  integer(psb_ipk_), allocatable  :: itemp(:)
  integer(psb_ipk_)               :: int_err(5), i1sz, i2sz, err_act, i, j
  integer(psb_epk_), parameter  :: one = 1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20) :: name
  name = 'psb_egelp'

  if(psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  i1sz = size(x, dim=1)
  i2sz = size(x, dim=2)

  if(debug_level >= psb_debug_serial_) write(debug_unit, *) trim(name), ': size', i1sz, i2sz

  allocate(temp(i1sz), itemp(size(iperm)), stat=info)
  if(info /= psb_success_) then
    info = 2040
    call psb_errpush(info, name)
    goto 9999
  end if
  itemp(:) = iperm(:) 

  if(.not. psb_isaperm(i1sz, itemp)) then
    info = psb_err_iarg_invalid_value_
    int_err(1) = 1      
    call psb_errpush(info, name, i_err = int_err)
    goto 9999
  endif
  select case(psb_toupper(trans))
    case('N') 
      do j = 1, i2sz
        do i = 1, i1sz
          temp(i) = x(itemp(i), j)
        end do
        do i = 1, i1sz
          x(i, j) = temp(i) 
        end do
      end do

    case('T')
      do j = 1, i2sz
        do i = 1, i1sz
          temp(itemp(i)) = x(i, j)
        end do
        do i = 1, i1sz
          x(i, j) = temp(i) 
        end do
      end do

    case default
      info = psb_err_from_subroutine_
      call psb_errpush(info, name, a_err='egelp')
  end select

  deallocate(temp, itemp)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_m_egelp

!!$ 
!!$              Parallel Sparse BLAS  version 3.5
!!$    (C) Copyright 2006-2018
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
!!$       software without specific prior written permission.
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
!
! Subroutine: psb_egelpv
!             Apply a left permutation to a dense matrix
!
! Arguments:
! trans    - character. 
! iperm    - integer.
! x        - real, dimension(:).
! info     - integer.                 Return code.
subroutine psb_m_egelpv(trans, iperm, x, info)
  use psb_serial_mod, psb_protect_name => psb_m_egelpv
  use psb_const_mod
  use psb_error_mod
  implicit none
  !Arguments
  character, intent(in)           :: trans
  integer(psb_mpk_), intent(in)   :: iperm(:)
  integer(psb_epk_), intent(inout)   :: x(:)
  integer(psb_ipk_), intent(out)  :: info

  !Local variables
  integer(psb_epk_), allocatable    :: temp(:)
  integer(psb_ipk_), allocatable  :: itemp(:)
  integer(psb_ipk_)           :: int_err(5), i1sz, err_act, i
  integer(psb_epk_), parameter  :: one = 1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20) :: name
  name = 'psb_egelpv'

  if(psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  i1sz = min(size(x), size(iperm))

  if(debug_level >= psb_debug_serial_) write(debug_unit, *)  trim(name), ': size', i1sz

  allocate(temp(i1sz), itemp(size(iperm)), stat=info)
  if(info /= psb_success_) then
    info = 2040
    call psb_errpush(info, name)
    goto 9999
  end if
  itemp(:) = iperm(:) 

  if(.not. psb_isaperm(i1sz, itemp)) then
    info = psb_err_iarg_invalid_value_
    int_err(1) = 1      
    call psb_errpush(info, name, i_err = int_err)
    goto 9999
  endif

  select case(psb_toupper(trans))
    case('N') 
      do i = 1, i1sz
        temp(i) = x(itemp(i))
      end do
      do i = 1, i1sz
        x(i) = temp(i) 
      end do

    case('T')
      do i = 1, i1sz
        temp(itemp(i)) = x(i)
      end do
      do i = 1, i1sz
        x(i) = temp(i) 
      end do

    case default
      info = psb_err_from_subroutine_
      call psb_errpush(info, name, a_err='egelp')
  end select

  deallocate(temp, itemp)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_m_egelpv

subroutine psb_e_egelp(trans, iperm, x, info)
  use psb_serial_mod, psb_protect_name => psb_e_egelp
  use psb_const_mod
  use psb_error_mod
  implicit none
  character, intent(in)           :: trans
  integer(psb_epk_), intent(in)   :: iperm(:)
  integer(psb_epk_), intent(inout)   :: x(:, :)
  integer(psb_ipk_), intent(out)  :: info

  ! local variables
  integer(psb_epk_), allocatable     :: temp(:)
  integer(psb_epk_), allocatable  :: itemp(:)
  integer(psb_ipk_)           :: int_err(5), err_act
  integer(psb_epk_)           :: i1sz, i2sz, i, j
  integer(psb_epk_), parameter  :: one = 1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20) :: name
  name = 'psb_egelp'

  if(psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  i1sz = size(x, dim=1)
  i2sz = size(x, dim=2)

  if(debug_level >= psb_debug_serial_) write(debug_unit, *)  trim(name), ': size', i1sz, i2sz

  allocate(temp(i1sz), itemp(size(iperm)), stat=info)
  if(info /= psb_success_) then
    info = 2040
    call psb_errpush(info, name)
    goto 9999
  end if
  itemp(:) = iperm(:) 

  if(.not. psb_isaperm(i1sz, itemp)) then
    info = psb_err_iarg_invalid_value_
    int_err(1) = 1      
    call psb_errpush(info, name, i_err = int_err)
    goto 9999
  endif
  select case(psb_toupper(trans))
    case('N') 
      do j = 1, i2sz
        do i = 1, i1sz
          temp(i) = x(itemp(i), j)
        end do
        do i = 1, i1sz
          x(i, j) = temp(i) 
        end do
      end do

    case('T')
      do j = 1, i2sz
        do i = 1, i1sz
          temp(itemp(i)) = x(i, j)
        end do
        do i = 1, i1sz
          x(i, j) = temp(i) 
        end do
      end do

    case default
      info = psb_err_from_subroutine_
      call psb_errpush(info, name, a_err='egelp')
  end select

  deallocate(temp, itemp)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_e_egelp

!!$ 
!!$              Parallel Sparse BLAS  version 3.5
!!$    (C) Copyright 2006-2018
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
!!$       software without specific prior written permission.
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
!
! Subroutine: psb_egelpv
!             Apply a left permutation to a dense matrix
!
! Arguments:
! trans    - character. 
! iperm    - integer.
! x        - real, dimension(:).
! info     - integer.                 Return code.
subroutine psb_e_egelpv(trans, iperm, x, info)
  use psb_serial_mod, psb_protect_name => psb_e_egelpv
  use psb_const_mod
  use psb_error_mod
  implicit none
  character, intent(in)           :: trans
  integer(psb_epk_), intent(in)   :: iperm(:)
  integer(psb_epk_), intent(inout)   :: x(:)
  integer(psb_ipk_), intent(out)  :: info

  ! local variables
  integer(psb_epk_), allocatable    :: temp(:)
  integer(psb_epk_), allocatable  :: itemp(:)
  integer(psb_epk_) :: i1sz, i
  integer(psb_ipk_) :: int_err(5), err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20) :: name

  name = 'psb_egelp'
  if(psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  i1sz = min(size(x), size(iperm))

  if(debug_level >= psb_debug_serial_) write(debug_unit, *)  trim(name), ': size', i1sz

  allocate(temp(i1sz), itemp(size(iperm)), stat=info)
  if(info /= psb_success_) then
    info = 2040
    call psb_errpush(info, name)
    goto 9999
  end if
  itemp(:) = iperm(:) 

  if(.not. psb_isaperm(i1sz, itemp)) then
    info = psb_err_iarg_invalid_value_
    int_err(1) = 1      
    call psb_errpush(info, name, i_err = int_err)
    goto 9999
  endif

  select case(psb_toupper(trans))
    case('N') 
      do i = 1, i1sz
        temp(i) = x(itemp(i))
      end do
      do i = 1, i1sz
        x(i) = temp(i) 
      end do

    case('T')
      do i = 1, i1sz
        temp(itemp(i)) = x(i)
      end do
      do i = 1, i1sz
        x(i) = temp(i) 
      end do

    case default
      info = psb_err_from_subroutine_
      call psb_errpush(info, name, a_err='egelp')
  end select

  deallocate(temp, itemp)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_e_egelpv

subroutine psi_eaxpby(m, n, alpha, x, beta, y, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m, n
  integer(psb_epk_), intent(in)      :: x(:, :)
  integer(psb_epk_), intent(inout)   :: y(:, :)
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_ipk_), intent(out)  :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly, i
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_
    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)
  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 4; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if((m > 0) .and. (n > 0)) call eaxpby(m, n, alpha, x, lx, beta, y, ly, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_eaxpby

subroutine psi_eaxpby2(m, n, alpha, x, beta, y, z, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m, n
  integer(psb_epk_), intent(in)      :: x(:, :)
  integer(psb_epk_), intent(in)      :: y(:, :)
  integer(psb_epk_), intent(inout)   :: z(:, :)
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly, lz, i
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_
    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if((m > 0) .and. (n > 0)) call eaxpbyv2(m, n, alpha, x, lx, beta, y, ly, z, lz, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_eaxpby2

subroutine psi_eaxpby3(m, n, alpha, x, beta, y, gamma, z, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m, n
  integer(psb_epk_), intent(in)      :: x(:, :)
  integer(psb_epk_), intent(in)      :: y(:, :)
  integer(psb_epk_), intent(inout)   :: z(:, :)
  integer(psb_epk_), intent(in)      :: alpha, beta, gamma
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly, lz, i
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_
    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if((m > 0) .and. (n > 0)) call eaxpbyv3(m, n, alpha, x, lx, beta, y, ly, gamma, z, lz, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_eaxpby3

subroutine psi_eaxpbyv(m, alpha, x, beta, y, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)      :: m
  integer(psb_epk_), intent(in)       :: x(:)
  integer(psb_epk_), intent(inout)    :: y(:)
  integer(psb_epk_), intent(in)       :: alpha, beta
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly
  integer(psb_ipk_) :: ierr(5)
  integer(psb_ipk_) :: i
  character(len=20) :: name, ch_err

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  lx = size(x, 1)
  ly = size(y, 1)
  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  ! if(m>0) call eaxpby(m, ione, alpha, x, lx, beta, y, ly, info)

  if(alpha.eq.ezero) then
    if(beta.eq.ezero) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = ezero
      enddo
    else if(beta.eq.eone) then
      !
      !        Do nothing!
      !

    else if(beta.eq.-eone) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i = 1, m
        y(i) =  beta*y(i)
      enddo
    endif

  else if(alpha.eq.eone) then

    if(beta.eq.ezero) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = x(i)
      enddo
    else if(beta.eq.eone) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = x(i) + y(i)
      enddo

    else if(beta.eq.-eone) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = x(i) + beta*y(i)
      enddo
    endif

  else if(alpha.eq.-eone) then

    if(beta.eq.ezero) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -x(i)
      enddo
    else if(beta.eq.eone) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -x(i) + y(i)
      enddo
    else if(beta.eq.-eone) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -x(i) + beta*y(i)
      enddo
    endif

  else

    if(beta.eq.ezero) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = alpha*x(i)
      enddo
    else if(beta.eq.eone) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = alpha*x(i) + y(i)
      enddo
    else if(beta.eq.-eone) then
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = alpha*x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = alpha*x(i) + beta*y(i)
      enddo
    endif
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_eaxpbyv

subroutine psi_eaxpbyv2(m, alpha, x, beta, y, z, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(in)      :: y(:)
  integer(psb_epk_), intent(inout)   :: z(:)
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_ipk_), intent(out)  :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly, lz, i
  integer(psb_ipk_) :: ierr(5)
  character(len=20)        :: name, ch_err

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)
  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(alpha.eq.ezero) then
    if(beta.eq.ezero) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = ezero
      enddo
    else if(beta.eq.eone) then
      !
      !        Do nothing!
      !

    else if(beta.eq.-eone) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) =  beta*y(i)
      enddo
    endif

  else if(alpha.eq.eone) then

    if(beta.eq.ezero) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = x(i)
      enddo
    else if(beta.eq.eone) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = x(i) + y(i)
      enddo

    else if(beta.eq.-eone) then
      !$omp parallel do private(i)
      do i = 1, m
          Z(i) = x(i) - y(i)
        enddo
    else
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = x(i) + beta*y(i)
      enddo
    endif

  else if(alpha.eq.-eone) then

    if(beta.eq.ezero) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = -x(i)
      enddo
    else if(beta.eq.eone) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = -x(i) + y(i)
      enddo

    else if(beta.eq.-eone) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = -x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = -x(i) + beta*y(i)
      enddo
    endif

  else

    if(beta.eq.ezero) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = alpha*x(i)
      enddo
    else if(beta.eq.eone) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = alpha*x(i) + y(i)
      enddo

    else if(beta.eq.-eone) then
      !$omp parallel do private(i)
      do i = 1, m
        Z(i) = alpha*x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i = 1, m
          Z(i) = alpha*x(i) + beta*y(i)
      enddo
    endif
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_eaxpbyv2

subroutine psi_eaxpbyv3(m, alpha, x, beta, y, gamma, z, info)
  use psi_e_serial_mod, psb_protect_name => psi_eaxpbyv3
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(in)      :: y(:)
  integer(psb_epk_), intent(inout)   :: z(:)
  integer(psb_epk_), intent(in)      :: alpha, beta, gamma
  integer(psb_ipk_), intent(out)  :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly, lz, i, code
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)
  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  ! Get the op-code based on the values of alpha, beta, gamma
  code = get_axpbylike_code(alpha, beta, gamma)

  select case (code)
    case( 0) ! (alpha, beta, gamma) = ( *,  *,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + beta*y(i) + gamma*z(i)
      end do
    case( 1) ! (alpha, beta, gamma) = ( 1,  *,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + beta*y(i) + gamma*z(i)
      end do
    case( 2) ! (alpha, beta, gamma) = ( 0,  *,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = beta*y(i) + gamma*z(i)
      end do
    case( 3) ! (alpha, beta, gamma) = (-1,  *,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + beta*y(i) + gamma*z(i)
      end do
    case( 4) ! (alpha, beta, gamma) = ( *,  1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + y(i) + gamma*z(i)
      end do
    case( 5) ! (alpha, beta, gamma) = ( 1,  1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + y(i) + gamma*z(i)
      end do
    case( 6) ! (alpha, beta, gamma) = ( 0,  1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = y(i) + gamma*z(i)
      end do
    case( 7) ! (alpha, beta, gamma) = (-1,  1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + y(i) + gamma*z(i)
      end do
    case( 8) ! (alpha, beta, gamma) = ( *,  0,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + gamma*z(i)
      end do
    case( 9) ! (alpha, beta, gamma) = ( 1,  0,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + gamma*z(i)
      end do
    case(10) ! (alpha, beta, gamma) = ( 0,  0,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = gamma*z(i)
      end do
    case(11) ! (alpha, beta, gamma) = (-1,  0,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + gamma*z(i)
      end do
    case(12) ! (alpha, beta, gamma) = ( *, -1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) - y(i) + gamma*z(i)
      end do
    case(13) ! (alpha, beta, gamma) = ( 1, -1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) - y(i) + gamma*z(i)
      end do
    case(14) ! (alpha, beta, gamma) = ( 0, -1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -y(i) + gamma*z(i)
      end do
    case(15) ! (alpha, beta, gamma) = (-1, -1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) - y(i) + gamma*z(i)
      end do
    case(16) ! (alpha, beta, gamma) = ( *,  *,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + beta*y(i) + z(i)
      end do
    case(17) ! (alpha, beta, gamma) = ( 1,  *,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + beta*y(i) + z(i)
      end do
    case(18) ! (alpha, beta, gamma) = ( 0,  *,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = beta*y(i) + z(i)
      end do
    case(19) ! (alpha, beta, gamma) = (-1,  *,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + beta*y(i) + z(i)
      end do
    case(20) ! (alpha, beta, gamma) = ( *,  1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + y(i) + z(i)
      end do
    case(21) ! (alpha, beta, gamma) = ( 1,  1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + y(i) + z(i)
      end do
    case(22) ! (alpha, beta, gamma) = ( 0,  1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = y(i) + z(i)
      end do
    case(23) ! (alpha, beta, gamma) = (-1,  1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + y(i) + z(i)
      end do
    case(24) ! (alpha, beta, gamma) = ( *,  0,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + z(i)
      end do
    case(25) ! (alpha, beta, gamma) = ( 1,  0,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + z(i)
      end do
    case(26) ! (alpha, beta, gamma) = ( 0,  0,  1)
      ! empty case: z(i) = z(i)
    case(27) ! (alpha, beta, gamma) = (-1,  0,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + z(i)
      end do
    case(28) ! (alpha, beta, gamma) = ( *, -1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) - y(i) + z(i)
      end do
    case(29) ! (alpha, beta, gamma) = ( 1, -1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) - y(i) + z(i)
      end do
    case(30) ! (alpha, beta, gamma) = ( 0, -1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -y(i) + z(i)
      end do
    case(31) ! (alpha, beta, gamma) = (-1, -1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) - y(i) + z(i)
      end do
    case(32) ! (alpha, beta, gamma) = ( *,  *,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + beta*y(i)
      end do
    case(33) ! (alpha, beta, gamma) = ( 1,  *,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + beta*y(i)
      end do
    case(34) ! (alpha, beta, gamma) = ( 0,  *,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = beta*y(i)
      end do
    case(35) ! (alpha, beta, gamma) = ( -1,  *,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + beta*y(i)
      end do
    case(36) ! (alpha, beta, gamma) = ( *,  1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + y(i)
      end do
    case(37) ! (alpha, beta, gamma) = ( 1,  1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + y(i)
      end do
    case(38) ! (alpha, beta, gamma) = ( 0,  1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = y(i)
      end do
    case(39) ! (alpha, beta, gamma) = (-1,  1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -y(i)
      end do
    case(40) ! (alpha, beta, gamma) = ( *,  0,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i)
      end do
    case(41) ! (alpha, beta, gamma) = ( 1,  0,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i)
      end do
    case(42) ! (alpha, beta, gamma) = ( 0,  0,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = ezero
      end do
    case(43) ! (alpha, beta, gamma) = (-1,  0,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i)
      end do
    case(44) ! (alpha, beta, gamma) = ( *, -1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) - y(i)
      end do
    case(45) ! (alpha, beta, gamma) = ( 1, -1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) - y(i)
      end do
    case(46) ! (alpha, beta, gamma) = ( 0, -1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -y(i)
      end do
    case(47) ! (alpha, beta, gamma) = (-1, -1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) - y(i)
      end do
    case(48) ! (alpha, beta, gamma) = ( *,  *, -1) 
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + beta*y(i) - z(i)
      end do
    case(49) ! (alpha, beta, gamma) = ( 1,  *, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + beta*y(i) - z(i)
      end do
    case(50) ! (alpha, beta, gamma) = ( 0,  *, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = beta*y(i) - z(i)
      end do
    case(51) ! (alpha, beta, gamma) = (-1,  *, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + beta*y(i) - z(i)
      end do
    case(52) ! (alpha, beta, gamma) = ( *,  1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) + y(i) - z(i)
      end do
    case(53) ! (alpha, beta, gamma) = ( 1,  1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) + y(i) - z(i)
      end do
    case(54) ! (alpha, beta, gamma) = ( 0,  1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = y(i) - z(i)
      end do
    case(55) ! (alpha, beta, gamma) = (-1,  1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) + y(i) - z(i)
      end do
    case(56) ! (alpha, beta, gamma) = ( *,  0, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) - z(i)
      end do
    case(57) ! (alpha, beta, gamma) = ( 1,  0, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) - z(i)
      end do
    case(58) ! (alpha, beta, gamma) = ( 0,  0, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -z(i)
      end do
    case(59) ! (alpha, beta, gamma) = (-1,  0, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) - z(i)
      end do
    case(60) ! (alpha, beta, gamma) = ( *, -1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = alpha*x(i) - y(i) - z(i)
      end do
    case(61) ! (alpha, beta, gamma) = ( 1, -1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = x(i) - y(i) - z(i)
      end do
    case(62) ! (alpha, beta, gamma) = ( 0, -1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -y(i) - z(i)
      end do
    case(63) ! (alpha, beta, gamma) = (-1, -1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          z(i) = -x(i) - y(i) - z(i)
      end do
    case default
      info = psb_err_internal_error_
      call psb_errpush(info, name)
      goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_eaxpbyv3

subroutine psi_eaxpbyv3_out(m, alpha, x, beta, y, gamma, z, w, info)
  use psi_e_serial_mod, psb_protect_name => psi_eaxpbyv3_out
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(in)      :: y(:)
  integer(psb_epk_), intent(in)      :: z(:)
  integer(psb_epk_), intent(inout)   :: w(:)
  integer(psb_epk_), intent(in)      :: alpha, beta, gamma
  integer(psb_ipk_), intent(out)  :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly, lz, lw, i, code
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name = 'psb_geaxpby'
  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)
  lw = size(w, 1)
  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 7; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  if(lw < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 8; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  ! Get the op-code based on the values of alpha, beta, gamma
  code = get_axpbylike_code(alpha, beta, gamma)

  select case (code)
    case( 0) ! (alpha, beta, gamma) = ( *,  *,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + beta*y(i) + gamma*z(i)
      end do
    case( 1) ! (alpha, beta, gamma) = ( 1,  *,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + beta*y(i) + gamma*z(i)
      end do
    case( 2) ! (alpha, beta, gamma) = ( 0,  *,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = beta*y(i) + gamma*z(i)
      end do
    case( 3) ! (alpha, beta, gamma) = (-1,  *,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + beta*y(i) + gamma*z(i)
      end do
    case( 4) ! (alpha, beta, gamma) = ( *,  1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + y(i) + gamma*z(i)
      end do
    case( 5) ! (alpha, beta, gamma) = ( 1,  1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + y(i) + gamma*z(i)
      end do
    case( 6) ! (alpha, beta, gamma) = ( 0,  1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = y(i) + gamma*z(i)
      end do
    case( 7) ! (alpha, beta, gamma) = (-1,  1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + y(i) + gamma*z(i)
      end do
    case( 8) ! (alpha, beta, gamma) = ( *,  0,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + gamma*z(i)
      end do
    case( 9) ! (alpha, beta, gamma) = ( 1,  0,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + gamma*z(i)
      end do
    case(10) ! (alpha, beta, gamma) = ( 0,  0,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = gamma*z(i)
      end do
    case(11) ! (alpha, beta, gamma) = (-1,  0,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + gamma*z(i)
      end do
    case(12) ! (alpha, beta, gamma) = ( *, -1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) - y(i) + gamma*z(i)
      end do
    case(13) ! (alpha, beta, gamma) = ( 1, -1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) - y(i) + gamma*z(i)
      end do
    case(14) ! (alpha, beta, gamma) = ( 0, -1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -y(i) + gamma*z(i)
      end do
    case(15) ! (alpha, beta, gamma) = (-1, -1,  *)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) - y(i) + gamma*z(i)
      end do
    case(16) ! (alpha, beta, gamma) = ( *,  *,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + beta*y(i) + z(i)
      end do
    case(17) ! (alpha, beta, gamma) = ( 1,  *,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + beta*y(i) + z(i)
      end do
    case(18) ! (alpha, beta, gamma) = ( 0,  *,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = beta*y(i) + z(i)
      end do
    case(19) ! (alpha, beta, gamma) = (-1,  *,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + beta*y(i) + z(i)
      end do
    case(20) ! (alpha, beta, gamma) = ( *,  1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + y(i) + z(i)
      end do
    case(21) ! (alpha, beta, gamma) = ( 1,  1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + y(i) + z(i)
      end do
    case(22) ! (alpha, beta, gamma) = ( 0,  1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = y(i) + z(i)
      end do
    case(23) ! (alpha, beta, gamma) = (-1,  1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + y(i) + z(i)
      end do
    case(24) ! (alpha, beta, gamma) = ( *,  0,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + z(i)
      end do
    case(25) ! (alpha, beta, gamma) = ( 1,  0,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + z(i)
      end do
    case(26) ! (alpha, beta, gamma) = ( 0,  0,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = z(i)
      end do
    case(27) ! (alpha, beta, gamma) = (-1,  0,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + z(i)
      end do
    case(28) ! (alpha, beta, gamma) = ( *, -1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) - y(i) + z(i)
      end do
    case(29) ! (alpha, beta, gamma) = ( 1, -1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) - y(i) + z(i)
      end do
    case(30) ! (alpha, beta, gamma) = ( 0, -1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -y(i) + z(i)
      end do
    case(31) ! (alpha, beta, gamma) = (-1, -1,  1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) - y(i) + z(i)
      end do
    case(32) ! (alpha, beta, gamma) = ( *,  *,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + beta*y(i)
      end do
    case(33) ! (alpha, beta, gamma) = ( 1,  *,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + beta*y(i)
      end do
    case(34) ! (alpha, beta, gamma) = ( 0,  *,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = beta*y(i)
      end do
    case(35) ! (alpha, beta, gamma) = ( -1,  *,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + beta*y(i)
      end do
    case(36) ! (alpha, beta, gamma) = ( *,  1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + y(i)
      end do
    case(37) ! (alpha, beta, gamma) = ( 1,  1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + y(i)
      end do
    case(38) ! (alpha, beta, gamma) = ( 0,  1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = y(i)
      end do
    case(39) ! (alpha, beta, gamma) = (-1,  1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -y(i)
      end do
    case(40) ! (alpha, beta, gamma) = ( *,  0,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i)
      end do
    case(41) ! (alpha, beta, gamma) = ( 1,  0,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i)
      end do
    case(42) ! (alpha, beta, gamma) = ( 0,  0,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = ezero
      end do
    case(43) ! (alpha, beta, gamma) = (-1,  0,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i)
      end do
    case(44) ! (alpha, beta, gamma) = ( *, -1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) - y(i)
      end do
    case(45) ! (alpha, beta, gamma) = ( 1, -1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) - y(i)
      end do
    case(46) ! (alpha, beta, gamma) = ( 0, -1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -y(i)
      end do
    case(47) ! (alpha, beta, gamma) = (-1, -1,  0)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) - y(i)
      end do
    case(48) ! (alpha, beta, gamma) = ( *,  *, -1) 
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + beta*y(i) - z(i)
      end do
    case(49) ! (alpha, beta, gamma) = ( 1,  *, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + beta*y(i) - z(i)
      end do
    case(50) ! (alpha, beta, gamma) = ( 0,  *, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = beta*y(i) - z(i)
      end do
    case(51) ! (alpha, beta, gamma) = (-1,  *, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + beta*y(i) - z(i)
      end do
    case(52) ! (alpha, beta, gamma) = ( *,  1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) + y(i) - z(i)
      end do
    case(53) ! (alpha, beta, gamma) = ( 1,  1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) + y(i) - z(i)
      end do
    case(54) ! (alpha, beta, gamma) = ( 0,  1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = y(i) - z(i)
      end do
    case(55) ! (alpha, beta, gamma) = (-1,  1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) + y(i) - z(i)
      end do
    case(56) ! (alpha, beta, gamma) = ( *,  0, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) - z(i)
      end do
    case(57) ! (alpha, beta, gamma) = ( 1,  0, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) - z(i)
      end do
    case(58) ! (alpha, beta, gamma) = ( 0,  0, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -z(i)
      end do
    case(59) ! (alpha, beta, gamma) = (-1,  0, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) - z(i)
      end do
    case(60) ! (alpha, beta, gamma) = ( *, -1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = alpha*x(i) - y(i) - z(i)
      end do
    case(61) ! (alpha, beta, gamma) = ( 1, -1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = x(i) - y(i) - z(i)
      end do
    case(62) ! (alpha, beta, gamma) = ( 0, -1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -y(i) - z(i)
      end do
    case(63) ! (alpha, beta, gamma) = (-1, -1, -1)
      !$omp parallel do private(i)
      do i = 1, m
          w(i) = -x(i) - y(i) - z(i)
      end do
    case default
      info = psb_err_internal_error_
      call psb_errpush(info, name)
      goto 9999
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_eaxpbyv3_out


subroutine psi_emlt(m, n, alpha, x, y, beta, info)
  use psi_e_serial_mod, psb_protect_name => psi_emlt
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m, n
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_epk_), intent(in)      :: x(:, :)
  integer(psb_epk_), intent(inout)   :: y(:, :)
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: lx, ly, i, j, code
  integer(psb_ipk_) :: ierr(5), err_act
  character name*20
  name = 'emlt'

  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 4; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  endif

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  
  ! Make a separate subroutine like in the axpby case?
  code = get_axpbylike_code(alpha, beta)
  select case (code)
    case( 0) ! (alpha, beta) = ( *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = alpha*x(i, j)*y(i, j) + beta*y(i, j)
        end do
      end do
    case( 1) ! (alpha, beta) = ( 1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = x(i, j)*y(i, j) + beta*y(i, j)
        end do
      end do
    case( 2) ! (alpha, beta) = ( 0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = beta*y(i, j)
        end do
      end do
    case( 3) ! (alpha, beta) = (-1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -x(i, j)*y(i, j) + beta*y(i, j)
        end do
      end do
    case( 4) ! (alpha, beta) = ( *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = alpha*x(i, j)*y(i, j) + y(i, j)
        end do
      end do
    case( 5) ! (alpha, beta) = ( 1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = x(i, j)*y(i, j) + y(i, j)
        end do
      end do
    case( 6) ! (alpha, beta) = ( 0,  1)
      ! empty case: y(i, j) = y(i, j)
    case( 7) ! (alpha, beta) = (-1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -x(i, j)*y(i, j) + y(i, j)
        end do
      end do
    case( 8) ! (alpha, beta) = ( *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = alpha*x(i, j)*y(i, j)
        end do
      end do
    case( 9) ! (alpha, beta) = ( 1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = x(i, j)*y(i, j)
        end do
      end do
    case(10) ! (alpha, beta) = ( 0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = 0
        end do
      end do
    case(11) ! (alpha, beta) = (-1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -x(i, j)*y(i, j)
        end do
      end do
    case(12) ! (alpha, beta) = ( *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = alpha*x(i, j)*y(i, j) - y(i, j)
        end do
      end do
    case(13) ! (alpha, beta) = ( 1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = x(i, j)*y(i, j) - y(i, j)
        end do
      end do
    case(14) ! (alpha, beta) = ( 0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -y(i, j)
        end do
      end do
    case(15) ! (alpha, beta) = (-1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -x(i, j)*y(i, j) - y(i, j)
        end do
      end do
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_emlt 

subroutine psi_emlt2(m, n, alpha, x, y, beta, z, info)
  use psi_e_serial_mod, psb_protect_name => psi_emlt2
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m, n
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_epk_), intent(in)      :: x(:, :)
  integer(psb_epk_), intent(in)      :: y(:, :)
  integer(psb_epk_), intent(inout)   :: z(:, :)
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: lx, ly, lz, i, j, code
  integer(psb_ipk_) :: ierr(5), err_act
  character name*20
  name = 'emlt'

  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 4; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  endif

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  
  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 7; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  ! Make a separate subroutine like in the axpby case?
  code = get_axpbylike_code(alpha, beta)
  select case (code)
    case( 0) ! (alpha, beta) = ( *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j)*y(i, j) + beta*z(i, j)
        end do
      end do
    case( 1) ! (alpha, beta) = ( 1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j)*y(i, j) + beta*z(i, j)
        end do
      end do
    case( 2) ! (alpha, beta) = ( 0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = beta*z(i, j)
        end do
      end do
    case( 3) ! (alpha, beta) = (-1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j)*y(i, j) + beta*z(i, j)
        end do
      end do
    case( 4) ! (alpha, beta) = ( *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j)*y(i, j) + z(i, j)
        end do
      end do
    case( 5) ! (alpha, beta) = ( 1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j)*y(i, j) + z(i, j)
        end do
      end do
    case( 6) ! (alpha, beta) = ( 0,  1)
      ! empty case: z(i, j) = z(i, j)
    case( 7) ! (alpha, beta) = (-1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j)*y(i, j) + z(i, j)
        end do
      end do
    case( 8) ! (alpha, beta) = ( *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j)*y(i, j)
        end do
      end do
    case( 9) ! (alpha, beta) = ( 1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j)*y(i, j)
        end do
      end do
    case(10) ! (alpha, beta) = ( 0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = 0
        end do
      end do
    case(11) ! (alpha, beta) = (-1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j)*y(i, j)
        end do
      end do
    case(12) ! (alpha, beta) = ( *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j)*y(i, j) - z(i, j)
        end do
      end do
    case(13) ! (alpha, beta) = ( 1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j)*y(i, j) - z(i, j)
        end do
      end do
    case(14) ! (alpha, beta) = ( 0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -y(i, j)
        end do
      end do
    case(15) ! (alpha, beta) = (-1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j)*y(i, j) - z(i, j)
        end do
      end do
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_emlt2 

subroutine psi_emltv(m, alpha, x, y, beta, info)
  use psi_e_serial_mod, psb_protect_name => psi_emltv
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(inout)   :: y(:)
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: lx, ly, i, code
  integer(psb_ipk_) :: ierr(5), err_act
  character name*20
  name = 'emlt'

  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  endif

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 4; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  
  ! Make a separate subroutine like in the axpby case?
  code = get_axpbylike_code(alpha, beta)
  select case (code)
    case( 0) ! (alpha, beta) = ( *,  *)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = alpha*x(i)*y(i) + beta*y(i)
      end do
    case( 1) ! (alpha, beta) = ( 1,  *)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = x(i)*y(i) + beta*y(i)
      end do
    case( 2) ! (alpha, beta) = ( 0,  *)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = beta*y(i)
      end do
    case( 3) ! (alpha, beta) = (-1,  *)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -x(i)*y(i) + beta*y(i)
      end do
    case( 4) ! (alpha, beta) = ( *,  1)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = alpha*x(i)*y(i) + y(i)
      end do
    case( 5) ! (alpha, beta) = ( 1,  1)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = x(i)*y(i) + y(i)
      end do
    case( 6) ! (alpha, beta) = ( 0,  1)
      ! empty case: y(i) = y(i)
    case( 7) ! (alpha, beta) = (-1,  1)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -x(i)*y(i) + y(i)
      end do
    case( 8) ! (alpha, beta) = ( *,  0)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = alpha*x(i)*y(i)
      end do
    case( 9) ! (alpha, beta) = ( 1,  0)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = x(i)*y(i)
      end do
    case(10) ! (alpha, beta) = ( 0,  0)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = 0
      end do
    case(11) ! (alpha, beta) = (-1,  0)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -x(i)*y(i)
      end do
    case(12) ! (alpha, beta) = ( *, -1)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = alpha*x(i)*y(i) - y(i)
      end do
    case(13) ! (alpha, beta) = ( 1, -1)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = x(i)*y(i) - y(i)
      end do
    case(14) ! (alpha, beta) = ( 0, -1)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -y(i)
      end do
    case(15) ! (alpha, beta) = (-1, -1)
      !$omp parallel do private(i)
      do i = 1, m
        y(i) = -x(i)*y(i) - y(i)
      end do
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_emltv 

subroutine psi_emltv2(m, alpha, x, y, beta, z, info)
  use psi_e_serial_mod, psb_protect_name => psi_emltv2
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(in)      :: y(:)
  integer(psb_epk_), intent(inout)   :: z(:)
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: lx, ly, lz, i, j, code
  integer(psb_ipk_) :: ierr(5), err_act
  character name*20
  name = 'emlt'

  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  
  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  endif

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  
  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  ! Make a separate subroutine like in the axpby case?
  code = get_axpbylike_code(alpha, beta)
  select case (code)
    case( 0) ! (alpha, beta) = ( *,  *)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = alpha*x(i)*y(i) + beta*z(i)
      end do
    case( 1) ! (alpha, beta) = ( 1,  *)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = x(i)*y(i) + beta*z(i)
      end do
    case( 2) ! (alpha, beta) = ( 0,  *)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = beta*z(i)
      end do
    case( 3) ! (alpha, beta) = (-1,  *)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = -x(i)*y(i) + beta*z(i)
      end do
    case( 4) ! (alpha, beta) = ( *,  1)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = alpha*x(i)*y(i) + z(i)
      end do
    case( 5) ! (alpha, beta) = ( 1,  1)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = x(i)*y(i) + z(i)
      end do
    case( 6) ! (alpha, beta) = ( 0,  1)
      ! empty case: z(i) = z(i)
    case( 7) ! (alpha, beta) = (-1,  1)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = -x(i)*y(i) + z(i)
      end do
    case( 8) ! (alpha, beta) = ( *,  0)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = alpha*x(i)*y(i)
      end do
    case( 9) ! (alpha, beta) = ( 1,  0)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = x(i)*y(i)
      end do
    case(10) ! (alpha, beta) = ( 0,  0)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = 0
      end do
    case(11) ! (alpha, beta) = (-1,  0)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = -x(i)*y(i)
      end do
    case(12) ! (alpha, beta) = ( *, -1)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = alpha*x(i)*y(i) - z(i)
      end do
    case(13) ! (alpha, beta) = ( 1, -1)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = x(i)*y(i) - z(i)
      end do
    case(14) ! (alpha, beta) = ( 0, -1)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = -y(i)
      end do
    case(15) ! (alpha, beta) = (-1, -1)
      !$omp parallel do private(i)
      do i = 1, m
        z(i) = -x(i)*y(i) - z(i)
      end do
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_emltv2

subroutine psi_emltx(m, n, alpha, x, y, beta, info)
  use psi_e_serial_mod, psb_protect_name => psi_emltx
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m, n
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(inout)   :: y(:, :)
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: lx, ly, i, j, code
  integer(psb_ipk_) :: ierr(5), err_act
  character name*20
  name = 'emlt'

  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 4; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  endif

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  
  ! Make a separate subroutine like in the axpby case?
  code = get_axpbylike_code(alpha, beta)
  select case (code)
    case( 0) ! (alpha, beta) = ( *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = alpha*x(i)*y(i, j) + beta*y(i, j)
        end do
      end do
    case( 1) ! (alpha, beta) = ( 1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = x(i)*y(i, j) + beta*y(i, j)
        end do
      end do
    case( 2) ! (alpha, beta) = ( 0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = beta*y(i, j)
        end do
      end do
    case( 3) ! (alpha, beta) = (-1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -x(i)*y(i, j) + beta*y(i, j)
        end do
      end do
    case( 4) ! (alpha, beta) = ( *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = alpha*x(i)*y(i, j) + y(i, j)
        end do
      end do
    case( 5) ! (alpha, beta) = ( 1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = x(i)*y(i, j) + y(i, j)
        end do
      end do
    case( 6) ! (alpha, beta) = ( 0,  1)
      ! empty case: y(i, j) = y(i, j)
    case( 7) ! (alpha, beta) = (-1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -x(i)*y(i, j) + y(i, j)
        end do
      end do
    case( 8) ! (alpha, beta) = ( *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = alpha*x(i)*y(i, j)
        end do
      end do
    case( 9) ! (alpha, beta) = ( 1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = x(i)*y(i, j)
        end do
      end do
    case(10) ! (alpha, beta) = ( 0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = 0
        end do
      end do
    case(11) ! (alpha, beta) = (-1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -x(i)*y(i, j)
        end do
      end do
    case(12) ! (alpha, beta) = ( *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = alpha*x(i)*y(i, j) - y(i, j)
        end do
      end do
    case(13) ! (alpha, beta) = ( 1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = x(i)*y(i, j) - y(i, j)
        end do
      end do
    case(14) ! (alpha, beta) = ( 0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -y(i, j)
        end do
      end do
    case(15) ! (alpha, beta) = (-1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            y(i, j) = -x(i)*y(i, j) - y(i, j)
        end do
      end do
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_emltx

subroutine psi_emltx2(m, n, alpha, x, y, beta, z, info)
  use psi_e_serial_mod, psb_protect_name => psi_emltx2
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m, n
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(in)      :: y(:, :)
  integer(psb_epk_), intent(inout)   :: z(:, :)
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: lx, ly, lz, i, j, code
  integer(psb_ipk_) :: ierr(5), err_act
  character name*20
  name = 'emlt'

  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 4; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  endif

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  
  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 7; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  ! Make a separate subroutine like in the axpby case?
  code = get_axpbylike_code(alpha, beta)
  select case (code)
    case( 0) ! (alpha, beta) = ( *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i)*y(i, j) + beta*z(i, j)
        end do
      end do
    case( 1) ! (alpha, beta) = ( 1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i)*y(i, j) + beta*z(i, j)
        end do
      end do
    case( 2) ! (alpha, beta) = ( 0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = beta*z(i, j)
        end do
      end do
    case( 3) ! (alpha, beta) = (-1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i)*y(i, j) + beta*z(i, j)
        end do
      end do
    case( 4) ! (alpha, beta) = ( *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i)*y(i, j) + z(i, j)
        end do
      end do
    case( 5) ! (alpha, beta) = ( 1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i)*y(i, j) + z(i, j)
        end do
      end do
    case( 6) ! (alpha, beta) = ( 0,  1)
      ! empty case: z(i, j) = z(i, j)
    case( 7) ! (alpha, beta) = (-1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i)*y(i, j) + z(i, j)
        end do
      end do
    case( 8) ! (alpha, beta) = ( *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i)*y(i, j)
        end do
      end do
    case( 9) ! (alpha, beta) = ( 1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i)*y(i, j)
        end do
      end do
    case(10) ! (alpha, beta) = ( 0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = 0
        end do
      end do
    case(11) ! (alpha, beta) = (-1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i)*y(i, j)
        end do
      end do
    case(12) ! (alpha, beta) = ( *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i)*y(i, j) - z(i, j)
        end do
      end do
    case(13) ! (alpha, beta) = ( 1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i)*y(i, j) - z(i, j)
        end do
      end do
    case(14) ! (alpha, beta) = ( 0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -y(i, j)
        end do
      end do
    case(15) ! (alpha, beta) = (-1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i)*y(i, j) - z(i, j)
        end do
      end do
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_emltx2

subroutine psi_emlte2(m, n, alpha, x, y, beta, z, info)
  use psi_e_serial_mod, psb_protect_name => psi_emlte2
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m, n
  integer(psb_epk_), intent(in)      :: alpha, beta
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(in)      :: y(:)
  integer(psb_epk_), intent(inout)   :: z(:, :)
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: lx, ly, lz, i, j, code
  integer(psb_ipk_) :: ierr(5), err_act
  character name*20
  name = 'emlt'

  info = psb_success_
  call psb_erractionsave(err_act)
  if(psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if(m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  if(n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  lx = size(x, 1)
  ly = size(y, 1)
  lz = size(z, 1)

  if(lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 4; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  endif

  if(ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if
  
  if(lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 7; ierr(2) = m
    call psb_errpush(info, name, i_err = ierr)
    goto 9999
  end if

  ! Make a separate subroutine like in the axpby case?
  code = get_axpbylike_code(alpha, beta)
  select case (code)
    case( 0) ! (alpha, beta) = ( *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i)*y(i) + beta*z(i, j)
        end do
      end do
    case( 1) ! (alpha, beta) = ( 1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i)*y(i) + beta*z(i, j)
        end do
      end do
    case( 2) ! (alpha, beta) = ( 0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = beta*z(i, j)
        end do
      end do
    case( 3) ! (alpha, beta) = (-1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i)*y(i) + beta*z(i, j)
        end do
      end do
    case( 4) ! (alpha, beta) = ( *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i)*y(i) + z(i, j)
        end do
      end do
    case( 5) ! (alpha, beta) = ( 1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i)*y(i) + z(i, j)
        end do
      end do
    case( 6) ! (alpha, beta) = ( 0,  1)
      ! empty case: z(i, j) = z(i, j)
    case( 7) ! (alpha, beta) = (-1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i)*y(i) + z(i, j)
        end do
      end do
    case( 8) ! (alpha, beta) = ( *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i)*y(i)
        end do
      end do
    case( 9) ! (alpha, beta) = ( 1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i)*y(i)
        end do
      end do
    case(10) ! (alpha, beta) = ( 0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = 0
        end do
      end do
    case(11) ! (alpha, beta) = (-1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i)*y(i)
        end do
      end do
    case(12) ! (alpha, beta) = ( *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i)*y(i) - z(i, j)
        end do
      end do
    case(13) ! (alpha, beta) = ( 1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i)*y(i) - z(i, j)
        end do
      end do
    case(14) ! (alpha, beta) = ( 0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -y(i)
        end do
      end do
    case(15) ! (alpha, beta) = (-1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i)*y(i) - z(i, j)
        end do
      end do
  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psi_emlte2 

subroutine psi_egthmv(n, k, idx, alpha, x, beta, y)
  use psb_const_mod
  implicit none
  integer(psb_mpk_) :: n, k
  integer(psb_ipk_) :: idx(:)
  integer(psb_epk_)    :: alpha, x(:, :), beta, y(:)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  if(beta == ezero) then
    if(alpha == ezero) then
      pt = 0
      do j = 1, k
        do i = 1, n
          pt = pt + 1 
          y(pt) = ezero
        end do
      end do
    else if(alpha == eone) then
      pt = 0
      do j = 1, k
        do i = 1, n
          pt = pt + 1 
          y(pt) = x(idx(i), j)
        end do
      end do
    else if(alpha == -eone) then
      pt = 0
      do j = 1, k
        do i = 1, n
          pt = pt + 1 
          y(pt) = -x(idx(i), j)
        end do
      end do
    else
      pt = 0
      do j = 1, k
        do i = 1, n
          pt = pt + 1 
          y(pt) = alpha*x(idx(i), j)
        end do
      end do
    end if
  else
    if(beta == eone) then
      ! Do nothing
    else if(beta == -eone) then
      y(1:n*k) = -y(1:n*k)
    else
      y(1:n*k) = beta*y(1:n*k)
    end if

    if(alpha == ezero) then
      ! do nothing
    else if(alpha == eone) then
      pt = 0
      do j = 1, k
        do i = 1, n
          pt = pt + 1 
          y(pt) = y(pt) + x(idx(i), j)
        end do
      end do
    else if(alpha == -eone) then
      pt = 0
      do j = 1, k
        do i = 1, n
          pt = pt + 1 
          y(pt) = y(pt) - x(idx(i), j)
        end do
      end do
    else
      pt = 0
      do j = 1, k
        do i = 1, n
          pt = pt + 1 
          y(pt) = y(pt) + alpha*x(idx(i), j)
        end do
      end do
    end if
  end if
end subroutine psi_egthmv

subroutine psi_egthv(n, idx, alpha, x, beta, y)
  use psb_const_mod
  implicit none
  integer(psb_mpk_) :: n
  integer(psb_ipk_) :: idx(:)
  integer(psb_epk_)    :: alpha, x(:), beta, y(:)

  ! Locals
  integer(psb_ipk_) :: i

  if(beta == ezero) then
    if(alpha == ezero) then
      do i = 1, n
        y(i) = ezero
      end do
    else if(alpha == eone) then
      do i = 1, n
        y(i) = x(idx(i))
      end do
    else if(alpha == -eone) then
      do i = 1, n
        y(i) = -x(idx(i))
      end do
    else
      do i = 1, n
        y(i) = alpha*x(idx(i))
      end do
    end if
  else
    if(beta == eone) then
      ! Do nothing
    else if(beta == -eone) then
      y(1:n) = -y(1:n)
    else
      y(1:n) = beta*y(1:n)
    end if

    if(alpha == ezero) then
      ! do nothing
    else if(alpha == eone) then
      do i = 1, n
        y(i) = y(i) + x(idx(i))
      end do
    else if(alpha == -eone) then
      do i = 1, n
        y(i) = y(i) - x(idx(i))
      end do
    else
      do i = 1, n
        y(i) = y(i) + alpha*x(idx(i))
      end do
    end if
  end if
end subroutine psi_egthv

subroutine psi_egthzmm(n, k, idx, x, y)
  use psb_const_mod
  implicit none
  integer(psb_mpk_) :: n, k
  integer(psb_ipk_) :: idx(:)
  integer(psb_epk_)    :: x(:, :), y(:, :)

  ! Locals
  integer(psb_ipk_) :: i

  do i = 1, n
    y(i, 1:k) = x(idx(i), 1:k)
  end do
end subroutine psi_egthzmm

subroutine psi_egthzmv(n, k, idx, x, y)
  use psb_const_mod
  implicit none
  integer(psb_mpk_) :: n, k
  integer(psb_ipk_) :: idx(:)
  integer(psb_epk_)    :: x(:, :), y(:)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  pt = 0
  do j = 1, k
    do i = 1, n
      pt = pt + 1 
      y(pt) = x(idx(i), j)
    end do
  end do
end subroutine psi_egthzmv

subroutine psi_egthzv(n, idx, x, y)
  use psb_const_mod
  implicit none
  integer(psb_mpk_) :: n
  integer(psb_ipk_) :: idx(:)
  integer(psb_epk_)    :: x(:), y(:)

  ! Locals
  integer(psb_ipk_) :: i

  do i = 1, n
    y(i) = x(idx(i))
  end do
end subroutine psi_egthzv

subroutine psi_esctmm(n, k, idx, x, beta, y)
  use psb_const_mod
  implicit none
  integer(psb_mpk_) :: n, k
  integer(psb_ipk_) :: idx(:)
  integer(psb_epk_)    :: x(:, :), beta, y(:, :)

  ! Locals
  integer(psb_ipk_) :: i, j

  if(beta == ezero) then
    do i = 1, n
      y(idx(i), 1:k) = x(i, 1:k)
    end do
  else if(beta == eone) then
    do i = 1, n
      y(idx(i), 1:k) = y(idx(i), 1:k)+x(i, 1:k)
    end do
  else
    do i = 1, n
      y(idx(i), 1:k) = beta*y(idx(i), 1:k)+x(i, 1:k)
    end do
  end if
end subroutine psi_esctmm

subroutine psi_esctmv(n, k, idx, x, beta, y)
  use psb_const_mod
  implicit none
  integer(psb_mpk_) :: n, k
  integer(psb_ipk_) :: idx(:)
  integer(psb_epk_)    :: x(:), beta, y(:, :)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  if(beta == ezero) then
    pt = 0
    do j = 1, k
      do i = 1, n
        pt = pt + 1 
        y(idx(i), j) = x(pt)
      end do
    end do
  else if(beta == eone) then
    pt = 0
    do j = 1, k
      do i = 1, n
        pt = pt + 1 
        y(idx(i), j) = y(idx(i), j)+x(pt)
      end do
    end do
  else
    pt = 0
    do j = 1, k
      do i = 1, n
        pt = pt + 1 
        y(idx(i), j) = beta*y(idx(i), j)+x(pt)
      end do
    end do
  end if
end subroutine psi_esctmv

subroutine psi_esctv(n, idx, x, beta, y)
  use psb_const_mod
  implicit none
  integer(psb_mpk_) :: n
  integer(psb_ipk_) :: idx(:)
  integer(psb_epk_) :: beta, x(:), y(:)

  ! Locals
  integer(psb_ipk_) :: i

  if(beta == ezero) then
    do i = 1, n
      y(idx(i)) = x(i)
    end do
  else if(beta == eone) then
    do i = 1, n
      y(idx(i)) = y(idx(i))+x(i)
    end do
  else
    do i = 1, n
      y(idx(i)) = beta*y(idx(i))+x(i)
    end do
  end if
end subroutine psi_esctv

subroutine eaxpby(m, n, alpha, X, lldx, beta, Y, lldy, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_) :: n, m, lldx, lldy, info
  integer(psb_epk_) X(lldx, *), Y(lldy, *)
  integer(psb_epk_) alpha, beta
  integer(psb_ipk_) :: i, j
  integer(psb_ipk_) :: int_err(5)
  character name*20
  name = 'eaxpby'

  !
  !     Error handling
  !
  info = psb_success_
  if(m .lt. 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(n .lt. 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = n
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(lldx .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 5
    int_err(2) = 1
    int_err(3) = lldx
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(lldy .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 8
    int_err(2) = 1
    int_err(3) = lldy
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif

  if(alpha.eq.ezero) then
    if(beta.eq.ezero) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = ezero
        enddo
      enddo
    else if(beta.eq.eone) then
      !
      !        Do nothing!
      !

    else if(beta.eq.-eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = - y(i, j)
        enddo
      enddo
    else
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) =  beta*y(i, j)
        enddo
      enddo
    endif

  else if(alpha.eq.eone) then

    if(beta.eq.ezero) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = x(i, j)
        enddo
      enddo
    else if(beta.eq.eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = x(i, j) + y(i, j)
        enddo
      enddo

    else if(beta.eq.-eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = x(i, j) - y(i, j)
        enddo
      enddo
    else
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = x(i, j) + beta*y(i, j)
        enddo
      enddo
    endif

  else if(alpha.eq.-eone) then

    if(beta.eq.ezero) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = -x(i, j)
        enddo
      enddo
    else if(beta.eq.eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = -x(i, j) + y(i, j)
        enddo
      enddo

    else if(beta.eq.-eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = -x(i, j) - y(i, j)
        enddo
      enddo
    else
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = -x(i, j) + beta*y(i, j)
        enddo
      enddo
    endif

  else

    if(beta.eq.ezero) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = alpha*x(i, j)
        enddo
      enddo
    else if(beta.eq.eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = alpha*x(i, j) + y(i, j)
        enddo
      enddo

    else if(beta.eq.-eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = alpha*x(i, j) - y(i, j)
        enddo
      enddo
    else
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          y(i, j) = alpha*x(i, j) + beta*y(i, j)
        enddo
      enddo
    endif
  endif
  return

9999 continue
  call fcpsb_serror()
  return
end subroutine eaxpby

subroutine eaxpbyv2(m, n, alpha, X, lldx, beta, Y, lldy, Z, lldz, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_) :: n, m, lldx, lldy, lldz, info
  integer(psb_epk_)    :: X(lldx, *), Y(lldy, *), Z(lldy, *)
  integer(psb_epk_)    :: alpha, beta
  integer(psb_ipk_) :: i, j
  integer(psb_ipk_) :: int_err(5)
  character name*20
  name = 'eaxpby'

  !
  !     Error handling
  !
  info = psb_success_
  if(m .lt. 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(n .lt. 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = n
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(lldx .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 5
    int_err(2) = 1
    int_err(3) = lldx
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(lldy .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 8
    int_err(2) = 1
    int_err(3) = lldy
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(lldz .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 8
    int_err(2) = 1
    int_err(3) = lldz
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif

  if(alpha.eq.ezero) then
    if(beta.eq.ezero) then
      do j = 1, n
        do i = 1, m
          Z(i, j) = ezero
        enddo
      enddo
    else if(beta.eq.eone) then
      !
      !        Do nothing!
      !

    else if(beta.eq.-eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = - y(i, j)
        enddo
      enddo
    else
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) =  beta*y(i, j)
        enddo
      enddo
    endif

  else if(alpha.eq.eone) then

    if(beta.eq.ezero) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = x(i, j)
        enddo
      enddo
    else if(beta.eq.eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = x(i, j) + y(i, j)
        enddo
      enddo

    else if(beta.eq.-eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = x(i, j) - y(i, j)
        enddo
      enddo
    else
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = x(i, j) + beta*y(i, j)
        enddo
      enddo
    endif

  else if(alpha.eq.-eone) then

    if(beta.eq.ezero) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = -x(i, j)
        enddo
      enddo
    else if(beta.eq.eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = -x(i, j) + y(i, j)
        enddo
      enddo

    else if(beta.eq.-eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = -x(i, j) - y(i, j)
        enddo
      enddo
    else
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = -x(i, j) + beta*y(i, j)
        enddo
      enddo
    endif

  else

    if(beta.eq.ezero) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = alpha*x(i, j)
        enddo
      enddo
    else if(beta.eq.eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = alpha*x(i, j) + y(i, j)
        enddo
      enddo

    else if(beta.eq.-eone) then
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = alpha*x(i, j) - y(i, j)
        enddo
      enddo
    else
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
          Z(i, j) = alpha*x(i, j) + beta*y(i, j)
        enddo
      enddo
    endif
  endif
  return

9999 continue
  call fcpsb_serror()
  return
end subroutine eaxpbyv2

subroutine eaxpbyv3(m, n, alpha, X, lldx, beta, Y, lldy, gamma, Z, lldz, info)
  use psi_e_serial_mod, psb_protect_name => psi_eaxpbyv3
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_) :: n, m, lldx, lldy, lldz, info
  integer(psb_epk_)    :: X(lldx, *), Y(lldy, *), Z(lldy, *)
  integer(psb_epk_)    :: alpha, beta, gamma
  integer(psb_ipk_) :: i, j, code
  integer(psb_ipk_) :: int_err(5)
  character name*20
  name = 'eaxpby'

  !
  !     Error handling
  !
  info = psb_success_
  if(m .lt. 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif

  if(n .lt. 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = n
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif

  if(lldx .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 5
    int_err(2) = 1
    int_err(3) = lldx
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif

  if(lldy .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 8
    int_err(2) = 1
    int_err(3) = lldy
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif

  if(lldz .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 8
    int_err(2) = 1
    int_err(3) = lldz
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif

  ! Get the op-code based on the values of alpha, beta, gamma
  code = get_axpbylike_code(alpha, beta, gamma)
  select case (code)
    case( 0) ! (alpha, beta, gamma) = ( *,  *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + beta*y(i, j) + gamma*z(i, j)
        end do
      end do
    case( 1) ! (alpha, beta, gamma) = ( 1,  *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + beta*y(i, j) + gamma*z(i, j)
        end do
      end do
    case( 2) ! (alpha, beta, gamma) = ( 0,  *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = beta*y(i, j) + gamma*z(i, j)
        end do
      end do
    case( 3) ! (alpha, beta, gamma) = (-1,  *,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + beta*y(i, j) + gamma*z(i, j)
        end do
      end do
    case( 4) ! (alpha, beta, gamma) = ( *,  1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + y(i, j) + gamma*z(i, j)
        end do
      end do
    case( 5) ! (alpha, beta, gamma) = ( 1,  1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + y(i, j) + gamma*z(i, j)
        end do
      end do
    case( 6) ! (alpha, beta, gamma) = ( 0,  1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = y(i, j) + gamma*z(i, j)
        end do
      end do
    case( 7) ! (alpha, beta, gamma) = (-1,  1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + y(i, j) + gamma*z(i, j)
        end do
      end do
    case( 8) ! (alpha, beta, gamma) = ( *,  0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + gamma*z(i, j)
        end do
      end do
    case( 9) ! (alpha, beta, gamma) = ( 1,  0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + gamma*z(i, j)
        end do
      end do
    case(10) ! (alpha, beta, gamma) = ( 0,  0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = gamma*z(i, j)
        end do
      end do
    case(11) ! (alpha, beta, gamma) = (-1,  0,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + gamma*z(i, j)
        end do
      end do
    case(12) ! (alpha, beta, gamma) = ( *, -1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) - y(i, j) + gamma*z(i, j)
        end do
      end do
    case(13) ! (alpha, beta, gamma) = ( 1, -1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) - y(i, j) + gamma*z(i, j)
        end do
      end do
    case(14) ! (alpha, beta, gamma) = ( 0, -1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -y(i, j) + gamma*z(i, j)
        end do
      end do
    case(15) ! (alpha, beta, gamma) = (-1, -1,  *)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) - y(i, j) + gamma*z(i, j)
        end do
      end do
    case(16) ! (alpha, beta, gamma) = ( *,  *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + beta*y(i, j) + z(i, j)
        end do
      end do
    case(17) ! (alpha, beta, gamma) = ( 1,  *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + beta*y(i, j) + z(i, j)
        end do
      end do
    case(18) ! (alpha, beta, gamma) = ( 0,  *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = beta*y(i, j) + z(i, j)
        end do
      end do
    case(19) ! (alpha, beta, gamma) = (-1,  *,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + beta*y(i, j) + z(i, j)
        end do
      end do
    case(20) ! (alpha, beta, gamma) = ( *,  1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + y(i, j) + z(i, j)
        end do
      end do
    case(21) ! (alpha, beta, gamma) = ( 1,  1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + y(i, j) + z(i, j)
        end do
      end do
    case(22) ! (alpha, beta, gamma) = ( 0,  1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = y(i, j) + z(i, j)
        end do
      end do
    case(23) ! (alpha, beta, gamma) = (-1,  1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + y(i, j) + z(i, j)
        end do
      end do
    case(24) ! (alpha, beta, gamma) = ( *,  0,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + z(i, j)
        end do
      end do
    case(25) ! (alpha, beta, gamma) = ( 1,  0,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + z(i, j)
        end do
      end do
    case(26) ! (alpha, beta, gamma) = ( 0,  0,  1)
      ! empty case: z(i) = z(i)
    case(27) ! (alpha, beta, gamma) = (-1,  0,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + z(i, j)
        end do
      end do
    case(28) ! (alpha, beta, gamma) = ( *, -1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) - y(i, j) + z(i, j)
        end do
      end do
    case(29) ! (alpha, beta, gamma) = ( 1, -1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) - y(i, j) + z(i, j)
        end do
      end do
    case(30) ! (alpha, beta, gamma) = ( 0, -1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -y(i, j) + z(i, j)
        end do
      end do
    case(31) ! (alpha, beta, gamma) = (-1, -1,  1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) - y(i, j) + z(i, j)
        end do
      end do
    case(32) ! (alpha, beta, gamma) = ( *,  *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + beta*y(i, j)
        end do
      end do
    case(33) ! (alpha, beta, gamma) = ( 1,  *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + beta*y(i, j)
        end do
      end do
    case(34) ! (alpha, beta, gamma) = ( 0,  *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = beta*y(i, j)
        end do
      end do
    case(35) ! (alpha, beta, gamma) = ( -1,  *,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + beta*y(i, j)
        end do
      end do
    case(36) ! (alpha, beta, gamma) = ( *,  1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + y(i, j)
        end do
      end do
    case(37) ! (alpha, beta, gamma) = ( 1,  1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + y(i, j)
        end do
      end do
    case(38) ! (alpha, beta, gamma) = ( 0,  1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = y(i, j)
        end do
      end do
    case(39) ! (alpha, beta, gamma) = (-1,  1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -y(i, j)
        end do
      end do
    case(40) ! (alpha, beta, gamma) = ( *,  0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j)
        end do
      end do
    case(41) ! (alpha, beta, gamma) = ( 1,  0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j)
        end do
      end do
    case(42) ! (alpha, beta, gamma) = ( 0,  0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = ezero
        end do
      end do
    case(43) ! (alpha, beta, gamma) = (-1,  0,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j)
        end do
      end do
    case(44) ! (alpha, beta, gamma) = ( *, -1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) - y(i, j)
        end do
      end do
    case(45) ! (alpha, beta, gamma) = ( 1, -1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) - y(i, j)
        end do
      end do
    case(46) ! (alpha, beta, gamma) = ( 0, -1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -y(i, j)
        end do
      end do
    case(47) ! (alpha, beta, gamma) = (-1, -1,  0)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) - y(i, j)
        end do
      end do
    case(48) ! (alpha, beta, gamma) = ( *,  *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + beta*y(i, j) - z(i, j)
        end do
      end do
    case(49) ! (alpha, beta, gamma) = ( 1,  *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + beta*y(i, j) - z(i, j)
        end do
      end do
    case(50) ! (alpha, beta, gamma) = ( 0,  *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = beta*y(i, j) - z(i, j)
        end do
      end do
    case(51) ! (alpha, beta, gamma) = (-1,  *, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + beta*y(i, j) - z(i, j)
        end do
      end do
    case(52) ! (alpha, beta, gamma) = ( *,  1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) + y(i, j) - z(i, j)
        end do
      end do
    case(53) ! (alpha, beta, gamma) = ( 1,  1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) + y(i, j) - z(i, j)
        end do
      end do
    case(54) ! (alpha, beta, gamma) = ( 0,  1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = y(i, j) - z(i, j)
        end do
      end do
    case(55) ! (alpha, beta, gamma) = (-1,  1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) + y(i, j) - z(i, j)
        end do
      end do
    case(56) ! (alpha, beta, gamma) = ( *,  0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) - z(i, j)
        end do
      end do
    case(57) ! (alpha, beta, gamma) = ( 1,  0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) - z(i, j)
        end do
      end do
    case(58) ! (alpha, beta, gamma) = ( 0,  0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -z(i, j)
        end do
      end do
    case(59) ! (alpha, beta, gamma) = (-1,  0, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) - z(i, j)
        end do
      end do
    case(60) ! (alpha, beta, gamma) = ( *, -1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = alpha*x(i, j) - y(i, j) - z(i, j)
        end do
      end do
    case(61) ! (alpha, beta, gamma) = ( 1, -1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = x(i, j) - y(i, j) - z(i, j)
        end do
      end do
    case(62) ! (alpha, beta, gamma) = ( 0, -1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -y(i, j) - z(i, j)
        end do
      end do
    case(63) ! (alpha, beta, gamma) = (-1, -1, -1)
      do j = 1, n
        !$omp parallel do private(i)
        do i = 1, m
            z(i, j) = -x(i, j) - y(i, j) - z(i, j)
        end do
      end do
    case default
      info = psb_err_internal_error_
      call psb_errpush(info, name)
      goto 9999
  end select
  return

9999 continue
  call fcpsb_serror()
  return
end subroutine eaxpbyv3

subroutine psi_e_upd_xyz(m, alpha, beta, gamma, delta, x, y, z, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)   :: m
  integer(psb_epk_), intent(in)      :: x(:)
  integer(psb_epk_), intent(inout)   :: y(:)
  integer(psb_epk_), intent(inout)   :: z(:)
  integer(psb_epk_), intent(in)      :: alpha, beta, gamma, delta
  integer(psb_ipk_), intent(out)  :: info

  integer(psb_ipk_) :: i
  integer(psb_ipk_) :: int_err(5)
  character name*20
  name = 'e_upd_xyz'

  info = psb_success_
  if(m .lt. 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(size(x) .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 6
    int_err(2) = 1
    int_err(3) = size(x)
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(size(y) .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 7
    int_err(2) = 1
    int_err(3) = size(y)
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(size(z) .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 8
    int_err(2) = 1
    int_err(3) = size(z)
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif
 
  if(beta == ezero) then
    if(gamma == ezero) then
      if(alpha == ezero) then
        if(delta == ezero) then
          !  a 0   b 0 g 0 d 0 
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = ezero
            z(i) = ezero
          end do
        else if(delta /= ezero) then
          !  a 0   b 0 g 0 d n 
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = ezero
            z(i) = delta*z(i)
          end do
        end if
      else if(alpha /= ezero) then
        if(delta == ezero) then
          !  a n   b 0 g 0 d 0 
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = alpha*x(i)
            z(i) = ezero
          end do
        else if(delta /= ezero) then
          !  a n   b 0 g 0 d n 
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = alpha*x(i)
            z(i) = delta*z(i)
          end do
        end if
      end if

    else  if(gamma /= ezero) then

      if(alpha == ezero) then
      
        if(delta == ezero) then
          !  a 0   b 0 g n d 0
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = ezero
            z(i) = ezero  ! gamma*y(i)
          end do
          
        else if(delta /= ezero) then
          !  a 0   b 0 g n d n
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = ezero
            z(i) = delta*z(i)
          end do
        end if

      else if(alpha /= ezero) then
        
        if(delta == ezero) then
          !  a n   b 0 g n d 0
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = alpha*x(i)
            z(i) = gamma*y(i)
          end do
          
        else if(delta /= ezero) then
          !  a n   b 0 g n d n
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = alpha*x(i)
            z(i) = gamma*y(i)+delta*z(i)
          end do
          
        end if
      end if
    end if

  else  if(beta /= ezero) then
    
    if(gamma == ezero) then
      if(alpha == ezero) then
        if(delta == ezero) then
          !  a 0   b n g 0 d 0
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = beta*y(i)
            z(i) = ezero
          end do
          
        else  if(delta /= ezero) then
          !  a 0   b n g 0 d n
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = beta*y(i)
            z(i) = delta*z(i)
          end do
          
        end if

      else  if(alpha /= ezero) then
        if(delta == ezero) then
          !  a n  b n g 0 d 0
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = alpha*x(i)+beta*y(i)
            z(i) = ezero
          end do
          
        else if(delta /= ezero) then
          !  a n  b n g 0 d n
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = alpha*x(i)+beta*y(i)
            z(i) = delta*z(i)
          end do
          
        end if

      end if
    else  if(gamma /= ezero) then
      if(alpha == ezero) then
        if(delta == ezero) then
          !  a 0  b n g n d 0
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = beta*y(i)
            z(i) = gamma*y(i)
          end do
          
        else if(delta /= ezero) then
          !  a 0  b n g n d n
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = beta*y(i)
            z(i) = gamma*y(i)+delta*z(i)
          end do

        end if

      else if(alpha /= ezero) then
        if(delta == ezero) then
          !  a n b n g n d 0
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = alpha*x(i)+beta*y(i)
            z(i) = gamma*y(i)
          end do
          
        else if(delta /= ezero) then
          !  a n b n g n d n
          !$omp parallel do private(i)
          do i = 1, m
            y(i) = alpha*x(i)+beta*y(i)
            z(i) = gamma*y(i)+delta*z(i)
          end do
          
        end if
      end if
    end if
  end if

  return

9999 continue
  call fcpsb_serror()
  return
end subroutine psi_e_upd_xyz

subroutine psi_exyzw(m, a, b, c, d, e, f, x, y, z, w, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)      :: m
  integer(psb_epk_), intent(in)       :: x(:)
  integer(psb_epk_), intent(inout)    :: y(:)
  integer(psb_epk_), intent(inout)    :: z(:)
  integer(psb_epk_), intent(inout)    :: w(:)
  integer(psb_epk_), intent(in)       :: a, b, c, d, e, f
  integer(psb_ipk_), intent(out)     :: info

  integer(psb_ipk_) :: i
  integer(psb_ipk_) :: int_err(5)
  character name*20
  name = 'e_xyzw'

  info = psb_success_
  if(m .lt. 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(size(x) .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 6
    int_err(2) = 1
    int_err(3) = size(x)
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(size(y) .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 7
    int_err(2) = 1
    int_err(3) = size(y)
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  else if(size(z) .lt. max(1, m)) then
    info = psb_err_iarg_not_gtia_ii_
    int_err(1) = 8
    int_err(2) = 1
    int_err(3) = size(z)
    int_err(4) = m
    call fcpsb_errpush(info, name, int_err)
    goto 9999
  endif

  if((a == ezero) .or. (b == ezero) .or. &
       & (c == ezero) .or. (d == ezero) .or. &
       & (e == ezero) .or. (f == ezero)) then
    write(0, *) 'XYZW assumes  a, b, c, d, e, f are all nonzero'
  else
    !$omp parallel do private(i)
    do i = 1, m
      y(i) = a*x(i) + b*y(i)
      z(i) = c*y(i) + d*z(i)
      w(i) = e*z(i) + f*w(i)
    end do
  end if
  return

9999 continue
  call fcpsb_serror()
  return
end subroutine psi_exyzw
