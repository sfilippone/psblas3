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
subroutine psi_z_exscanv(n,x,info,shift)
  use psi_z_serial_mod, psb_protect_name => psi_z_exscanv
  use psb_const_mod
  use psb_error_mod
#if defined(OPENMP)
  use omp_lib
#endif
  implicit none
  integer(psb_ipk_), intent(in)      :: n
  complex(psb_dpk_), intent (inout)    :: x(:)
  integer(psb_ipk_), intent(out)     :: info
  complex(psb_dpk_), intent(in), optional :: shift
  
  complex(psb_dpk_) :: shift_, tp, ts
  integer(psb_ipk_) :: i
  logical is_nested, is_parallel
  
  if (present(shift)) then
    shift_ = shift
  else
    shift_ = zzero
  end if
    
#if defined(OPENMP)
  is_parallel = omp_in_parallel()
  if (is_parallel) then 
    call inner_z_exscan()
  else
    !$OMP PARALLEL default(shared) 
    call inner_z_exscan()
    !$OMP END PARALLEL
  end if
#else
  tp = shift_
  do i=1,n
    ts = x(i)
    x(i) = tp
    tp = tp + ts
  end do

#endif
#if defined(OPENMP)
contains
  subroutine inner_z_exscan()
    ! Note: all these variables are private, but SUMB should *really* be
    ! a pointer. The semantics of COPYPRIVATE is that the POINTER is copied
    ! so effectively we are recovering a SHARED SUMB which is what
    ! we need in this case. If it was an ALLOCATABLE, then it would be the contents
    ! that would get copied, and the SHARED effect would  no longer be there.
    ! Simple parallel version of EXSCAN
    integer(psb_ipk_)  :: i,ithread,nthreads,idxstart,idxend,wrk
    complex(psb_dpk_), pointer :: sumb(:)
    complex(psb_dpk_)   :: tp, ts

    nthreads = omp_get_num_threads()
    ithread = omp_get_thread_num()
    !$OMP SINGLE
    allocate(sumb(nthreads+1))
    sumb(:) = 0
    !$OMP END SINGLE COPYPRIVATE(sumb)

    wrk = (n)/nthreads
    if (ithread < MOD((n),nthreads)) then
      wrk = wrk + 1
      idxstart = ithread*wrk + 1
    else
      idxstart = ithread*wrk + MOD((n),nthreads) + 1
    end if

    idxend = min(idxstart + wrk - 1,n )
    tp = zzero
    if (idxstart<=idxend) then
      do i=idxstart,idxend
        ts = x(i)
        x(i) = tp 
        tp = tp + ts 
      end do
    end if
    sumb(ithread+2) = tp 
    !$OMP BARRIER
    
    !$OMP SINGLE
    do i=2,nthreads+1
      sumb(i) = sumb(i) + sumb(i-1)
    end do
    !$OMP END SINGLE      

    !$OMP BARRIER

    !$OMP DO SCHEDULE(STATIC)
    do i=1,n
      x(i) = x(i) + sumb(ithread+1) + shift_ 
    end do
    !$OMP END DO
    !$OMP SINGLE
    deallocate(sumb)
    !$OMP END SINGLE
  end subroutine inner_z_exscan
#endif
end subroutine psi_z_exscanv

subroutine psb_m_zgelp(trans,iperm,x,info)
  use psb_serial_mod, psb_protect_name => psb_m_zgelp
  use psb_const_mod
  use psb_error_mod
  implicit none

  complex(psb_dpk_), intent(inout) ::  x(:,:)
  integer(psb_mpk_), intent(in)           ::  iperm(:)
  integer(psb_ipk_), intent(out)          ::  info
  character, intent(in)         :: trans

  ! local variables
  complex(psb_dpk_),allocatable :: temp(:)
  integer(psb_ipk_) :: int_err(5), i1sz, i2sz, err_act,i,j
  integer(psb_ipk_), allocatable          :: itemp(:)
  complex(psb_dpk_),parameter   :: one=1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name
  name = 'psb_zgelp'

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  i1sz    = size(x,dim=1)
  i2sz    = size(x,dim=2)

  if (debug_level >= psb_debug_serial_)&
       & write(debug_unit,*)  trim(name),': size',i1sz,i2sz

  allocate(temp(i1sz),itemp(size(iperm)),stat=info)
  if (info /= psb_success_) then
    info=2040
    call psb_errpush(info,name)
    goto 9999
  end if
  itemp(:) = iperm(:) 

  if (.not.psb_isaperm(i1sz,itemp)) then
    info=psb_err_iarg_invalid_value_
    int_err(1) = 1      
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif
  select case( psb_toupper(trans))
  case('N') 
    do j=1,i2sz
      do i=1,i1sz
        temp(i) = x(itemp(i),j)
      end do
      do i=1,i1sz
        x(i,j) = temp(i) 
      end do
    end do
  case('T')
    do j=1,i2sz
      do i=1,i1sz
        temp(itemp(i)) = x(i,j)
      end do
      do i=1,i1sz
        x(i,j) = temp(i) 
      end do
    end do
  case default
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='zgelp')
  end select

  deallocate(temp,itemp)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_m_zgelp



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
!
! Subroutine: psb_zgelpv
!             Apply a left permutation to a dense matrix
!
! Arguments:
! trans    - character. 
! iperm    - integer.
! x        - real, dimension(:).
! info     - integer.                 Return code.
subroutine psb_m_zgelpv(trans,iperm,x,info)
  use psb_serial_mod, psb_protect_name => psb_m_zgelpv
  use psb_const_mod
  use psb_error_mod
  implicit none

  complex(psb_dpk_), intent(inout) ::  x(:)
  integer(psb_mpk_), intent(in)           ::  iperm(:)
  integer(psb_ipk_), intent(out)          ::  info
  character, intent(in)         ::  trans

  ! local variables
  integer(psb_ipk_) :: int_err(5), i1sz, err_act, i
  complex(psb_dpk_),allocatable ::  temp(:)
  integer(psb_ipk_), allocatable          :: itemp(:)
  complex(psb_dpk_),parameter   :: one=1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name
  name = 'psb_zgelpv'

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  i1sz = min(size(x),size(iperm))

  if (debug_level >= psb_debug_serial_)&
       & write(debug_unit,*)  trim(name),': size',i1sz
  allocate(temp(i1sz),itemp(size(iperm)),stat=info)
  if (info /= psb_success_) then
    info=2040
    call psb_errpush(info,name)
    goto 9999
  end if
  itemp(:) = iperm(:) 

  if (.not.psb_isaperm(i1sz,itemp)) then
    info=psb_err_iarg_invalid_value_
    int_err(1) = 1      
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  select case( psb_toupper(trans))
  case('N') 
    do i=1,i1sz
      temp(i) = x(itemp(i))
    end do
    do i=1,i1sz
      x(i) = temp(i) 
    end do
  case('T')
    do i=1,i1sz
      temp(itemp(i)) = x(i)
    end do
    do i=1,i1sz
      x(i) = temp(i) 
    end do
  case default
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='zgelp')
  end select

  deallocate(temp,itemp)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_m_zgelpv

subroutine psb_e_zgelp(trans,iperm,x,info)
  use psb_serial_mod, psb_protect_name => psb_e_zgelp
  use psb_const_mod
  use psb_error_mod
  implicit none

  complex(psb_dpk_), intent(inout) ::  x(:,:)
  integer(psb_epk_), intent(in)           ::  iperm(:)
  integer(psb_ipk_), intent(out)          ::  info
  character, intent(in)         :: trans

  ! local variables
  complex(psb_dpk_),allocatable :: temp(:)
  integer(psb_ipk_) :: int_err(5), err_act
  integer(psb_epk_) :: i1sz, i2sz, i, j
  integer(psb_epk_), allocatable          :: itemp(:)
  complex(psb_dpk_),parameter   :: one=1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name
  name = 'psb_zgelp'

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  i1sz    = size(x,dim=1)
  i2sz    = size(x,dim=2)

  if (debug_level >= psb_debug_serial_)&
       & write(debug_unit,*)  trim(name),': size',i1sz,i2sz

  allocate(temp(i1sz),itemp(size(iperm)),stat=info)
  if (info /= psb_success_) then
    info=2040
    call psb_errpush(info,name)
    goto 9999
  end if
  itemp(:) = iperm(:) 

  if (.not.psb_isaperm(i1sz,itemp)) then
    info=psb_err_iarg_invalid_value_
    int_err(1) = 1      
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif
  select case( psb_toupper(trans))
  case('N') 
    do j=1,i2sz
      do i=1,i1sz
        temp(i) = x(itemp(i),j)
      end do
      do i=1,i1sz
        x(i,j) = temp(i) 
      end do
    end do
  case('T')
    do j=1,i2sz
      do i=1,i1sz
        temp(itemp(i)) = x(i,j)
      end do
      do i=1,i1sz
        x(i,j) = temp(i) 
      end do
    end do
  case default
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='zgelp')
  end select

  deallocate(temp,itemp)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_e_zgelp



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
!
! Subroutine: psb_zgelpv
!             Apply a left permutation to a dense matrix
!
! Arguments:
! trans    - character. 
! iperm    - integer.
! x        - real, dimension(:).
! info     - integer.                 Return code.
subroutine psb_e_zgelpv(trans,iperm,x,info)
  use psb_serial_mod, psb_protect_name => psb_e_zgelpv
  use psb_const_mod
  use psb_error_mod
  implicit none

  complex(psb_dpk_), intent(inout) ::  x(:)
  integer(psb_epk_), intent(in)           ::  iperm(:)
  integer(psb_ipk_), intent(out)          ::  info
  character, intent(in)         ::  trans

  ! local variables
  integer(psb_ipk_) :: int_err(5), err_act
  complex(psb_dpk_),allocatable ::  temp(:)
  integer(psb_epk_) :: i1sz, i
  integer(psb_epk_), allocatable          :: itemp(:)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20) :: name

  name = 'psb_zgelp'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  i1sz = min(size(x),size(iperm))

  if (debug_level >= psb_debug_serial_)&
       & write(debug_unit,*)  trim(name),': size',i1sz
  allocate(temp(i1sz),itemp(size(iperm)),stat=info)
  if (info /= psb_success_) then
    info=2040
    call psb_errpush(info,name)
    goto 9999
  end if
  itemp(:) = iperm(:) 

  if (.not.psb_isaperm(i1sz,itemp)) then
    info=psb_err_iarg_invalid_value_
    int_err(1) = 1      
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  select case( psb_toupper(trans))
  case('N') 
    do i=1,i1sz
      temp(i) = x(itemp(i))
    end do
    do i=1,i1sz
      x(i) = temp(i) 
    end do
  case('T')
    do i=1,i1sz
      temp(itemp(i)) = x(i)
    end do
    do i=1,i1sz
      x(i) = temp(i) 
    end do
  case default
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='zgelp')
  end select

  deallocate(temp,itemp)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_e_zgelpv

subroutine psi_zaxpby(m,n,alpha, x, beta, y, info)

  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)      :: m, n
  complex(psb_dpk_), intent (in)       ::  x(:,:)
  complex(psb_dpk_), intent (inout)    ::  y(:,:)
  complex(psb_dpk_), intent (in)       ::  alpha, beta
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly, i
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name='psb_geaxpby'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if (m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (n < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = n
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  lx = size(x,1)
  ly = size(y,1)
  if (lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 4; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 6; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if ((m>0).and.(n>0)) call zaxpby(m,n,alpha,x,lx,beta,y,ly,info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine psi_zaxpby

subroutine psi_zaxpbyv(m,alpha, x, beta, y, info)

  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)      :: m
  complex(psb_dpk_), intent (in)       ::  x(:)
  complex(psb_dpk_), intent (inout)    ::  y(:)
  complex(psb_dpk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly
  integer(psb_ipk_) :: ierr(5)
  integer(psb_ipk_) :: i
  character(len=20) :: name, ch_err

  name='psb_geaxpby'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if (m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  lx = size(x,1)
  ly = size(y,1)
  if (lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

!  if (m>0) call zaxpby(m,ione,alpha,x,lx,beta,y,ly,info)

  if (alpha.eq.zzero) then
    if (beta.eq.zzero) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = zzero
      enddo
    else if (beta.eq.zone) then
      !
      !        Do nothing!
      !

    else if (beta.eq.-zone) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i=1,m
        y(i) =  beta*y(i)
      enddo
    endif

  else if (alpha.eq.zone) then

    if (beta.eq.zzero) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = x(i)
      enddo
    else if (beta.eq.zone) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = x(i) + y(i)
      enddo

    else if (beta.eq.-zone) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i=1,m
        y(i) = x(i) + beta*y(i)
      enddo
    endif

  else if (alpha.eq.-zone) then

    if (beta.eq.zzero) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = -x(i)
      enddo
    else if (beta.eq.zone) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = -x(i) + y(i)
      enddo
    else if (beta.eq.-zone) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = -x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i=1,m
        y(i) = -x(i) + beta*y(i)
      enddo
    endif

  else

    if (beta.eq.zzero) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = alpha*x(i)
      enddo
    else if (beta.eq.zone) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = alpha*x(i) + y(i)
      enddo
    else if (beta.eq.-zone) then
      !$omp parallel do private(i)
      do i=1,m
        y(i) = alpha*x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i=1,m
        y(i) = alpha*x(i) + beta*y(i)
      enddo
    endif

  endif


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psi_zaxpbyv

subroutine psi_zaxpbyv2(m,alpha, x, beta, y, z, info)

  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)      :: m
  complex(psb_dpk_), intent (in)       ::  x(:)
  complex(psb_dpk_), intent (in)       ::  y(:)
  complex(psb_dpk_), intent (inout)    ::  z(:)
  complex(psb_dpk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly, lz, i
  integer(psb_ipk_) :: ierr(5)
  character(len=20)        :: name, ch_err

  name='psb_geaxpby'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  if (m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  lx = size(x,1)
  ly = size(y,1)
  lz = size(z,1)
  if (lx < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (ly < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (lz < m) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (alpha.eq.zzero) then
    if (beta.eq.zzero) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = zzero
      enddo
    else if (beta.eq.zone) then
      !
      !        Do nothing!
      !

    else if (beta.eq.-zone) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i=1,m
        Z(i) =  beta*y(i)
      enddo
    endif

  else if (alpha.eq.zone) then

    if (beta.eq.zzero) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = x(i)
      enddo
    else if (beta.eq.zone) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = x(i) + y(i)
      enddo

    else if (beta.eq.-zone) then
      !$omp parallel do private(i)
      do i=1,m
          Z(i) = x(i) - y(i)
        enddo
    else
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = x(i) + beta*y(i)
      enddo
    endif

  else if (alpha.eq.-zone) then

    if (beta.eq.zzero) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = -x(i)
      enddo
    else if (beta.eq.zone) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = -x(i) + y(i)
      enddo

    else if (beta.eq.-zone) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = -x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = -x(i) + beta*y(i)
      enddo
    endif

  else

    if (beta.eq.zzero) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = alpha*x(i)
      enddo
    else if (beta.eq.zone) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = alpha*x(i) + y(i)
      enddo

    else if (beta.eq.-zone) then
      !$omp parallel do private(i)
      do i=1,m
        Z(i) = alpha*x(i) - y(i)
      enddo
    else
      !$omp parallel do private(i)
      do i=1,m
          Z(i) = alpha*x(i) + beta*y(i)
      enddo
    endif

  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psi_zaxpbyv2

subroutine psi_zgthmv(n,k,idx,alpha,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  complex(psb_dpk_) :: x(:,:), y(:),alpha,beta

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  if (beta == zzero) then
    if (alpha == zzero) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = zzero
        end do
      end do
    else if (alpha == zone) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = x(idx(i),j)
        end do
      end do
    else if (alpha == -zone) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = -x(idx(i),j)
        end do
      end do
    else
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = alpha*x(idx(i),j)
        end do
      end do
    end if
  else
    if (beta == zone) then
      ! Do nothing
    else if (beta == -zone) then
      y(1:n*k) = -y(1:n*k)
    else
      y(1:n*k) = beta*y(1:n*k)
    end if

    if (alpha == zzero) then
      ! do nothing
    else if (alpha == zone) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = y(pt) + x(idx(i),j)
        end do
      end do
    else if (alpha == -zone) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = y(pt) - x(idx(i),j)
        end do
      end do
    else
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = y(pt) + alpha*x(idx(i),j)
        end do
      end do
    end if
  end if

end subroutine psi_zgthmv

subroutine psi_zgthv(n,idx,alpha,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, idx(:)
  complex(psb_dpk_) :: x(:), y(:),alpha,beta

  ! Locals
  integer(psb_ipk_) :: i
  if (beta == zzero) then
    if (alpha == zzero) then
      do i=1,n
        y(i) = zzero
      end do
    else if (alpha == zone) then
      do i=1,n
        y(i) = x(idx(i))
      end do
    else if (alpha == -zone) then
      do i=1,n
        y(i) = -x(idx(i))
      end do
    else
      do i=1,n
        y(i) = alpha*x(idx(i))
      end do
    end if
  else
    if (beta == zone) then
      ! Do nothing
    else if (beta == -zone) then
      y(1:n) = -y(1:n)
    else
      y(1:n) = beta*y(1:n)
    end if

    if (alpha == zzero) then
      ! do nothing
    else if (alpha == zone) then
      do i=1,n
        y(i) = y(i) + x(idx(i))
      end do
    else if (alpha == -zone) then
      do i=1,n
        y(i) = y(i) - x(idx(i))
      end do
    else
      do i=1,n
        y(i) = y(i) + alpha*x(idx(i))
      end do
    end if
  end if

end subroutine psi_zgthv

subroutine psi_zgthzmm(n,k,idx,x,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  complex(psb_dpk_) :: x(:,:), y(:,:)

  ! Locals
  integer(psb_ipk_) :: i


  do i=1,n
    y(i,1:k)=x(idx(i),1:k)
  end do

end subroutine psi_zgthzmm

subroutine psi_zgthzmv(n,k,idx,x,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  complex(psb_dpk_) :: x(:,:), y(:)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  pt=0
  do j=1,k
    do i=1,n
      pt=pt+1
      y(pt)=x(idx(i),j)
    end do
  end do

end subroutine psi_zgthzmv

subroutine psi_zgthzv(n,idx,x,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, idx(:)
  complex(psb_dpk_) :: x(:), y(:)

  ! Locals
  integer(psb_ipk_) :: i

  do i=1,n
    y(i)=x(idx(i))
  end do

end subroutine psi_zgthzv

subroutine psi_zsctmm(n,k,idx,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  complex(psb_dpk_) :: beta, x(:,:), y(:,:)

  ! Locals
  integer(psb_ipk_) :: i, j

  if (beta == zzero) then
    do i=1,n
      y(idx(i),1:k) = x(i,1:k)
    end do
  else if (beta == zone) then
    do i=1,n
      y(idx(i),1:k) = y(idx(i),1:k)+x(i,1:k)
    end do
  else
    do i=1,n
      y(idx(i),1:k) = beta*y(idx(i),1:k)+x(i,1:k)
    end do
  end if
end subroutine psi_zsctmm

subroutine psi_zsctmv(n,k,idx,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  complex(psb_dpk_) :: beta, x(:), y(:,:)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  if (beta == zzero) then
    pt=0
    do j=1,k
      do i=1,n
        pt=pt+1
        y(idx(i),j) = x(pt)
      end do
    end do
  else if (beta == zone) then
    pt=0
    do j=1,k
      do i=1,n
        pt=pt+1
        y(idx(i),j) = y(idx(i),j)+x(pt)
      end do
    end do
  else
    pt=0
    do j=1,k
      do i=1,n
        pt=pt+1
        y(idx(i),j) = beta*y(idx(i),j)+x(pt)
      end do
    end do
  end if
end subroutine psi_zsctmv

subroutine psi_zsctv(n,idx,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, idx(:)
  complex(psb_dpk_) :: beta, x(:), y(:)

  ! Locals
  integer(psb_ipk_) :: i

  if (beta == zzero) then
    do i=1,n
      y(idx(i)) = x(i)
    end do
  else if (beta == zone) then
    do i=1,n
      y(idx(i)) = y(idx(i))+x(i)
    end do
  else
    do i=1,n
      y(idx(i)) = beta*y(idx(i))+x(i)
    end do
  end if
end subroutine psi_zsctv

subroutine  zaxpby(m, n, alpha, X, lldx, beta, Y, lldy, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_) :: n, m, lldx, lldy, info
  complex(psb_dpk_) X(lldx,*), Y(lldy,*)
  complex(psb_dpk_) alpha, beta
  integer(psb_ipk_) :: i, j
  integer(psb_ipk_) :: int_err(5)
  character  name*20
  name='zaxpby'


  !
  !     Error handling
  !
  info = psb_success_
  if (m.lt.0) then
    info=psb_err_iarg_neg_
    int_err(1)=1
    int_err(2)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (n.lt.0) then
    info=psb_err_iarg_neg_
    int_err(1)=1
    int_err(2)=n
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (lldx.lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=5
    int_err(2)=1
    int_err(3)=lldx
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (lldy.lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=8
    int_err(2)=1
    int_err(3)=lldy
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  endif

  if (alpha.eq.zzero) then
    if (beta.eq.zzero) then
      do j=1, n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = zzero
        enddo
      enddo
    else if (beta.eq.zone) then
      !
      !        Do nothing!
      !

    else if (beta.eq.-zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = - y(i,j)
        enddo
      enddo
    else
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) =  beta*y(i,j)
        enddo
      enddo
    endif

  else if (alpha.eq.zone) then

    if (beta.eq.zzero) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = x(i,j) - y(i,j)
        enddo
      enddo
    else
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = x(i,j) + beta*y(i,j)
        enddo
      enddo
    endif

  else if (alpha.eq.-zone) then

    if (beta.eq.zzero) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = -x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = -x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = -x(i,j) - y(i,j)
        enddo
      enddo
    else
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = -x(i,j) + beta*y(i,j)
        enddo
      enddo
    endif

  else

    if (beta.eq.zzero) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = alpha*x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = alpha*x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = alpha*x(i,j) - y(i,j)
        enddo
      enddo
    else
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          y(i,j) = alpha*x(i,j) + beta*y(i,j)
        enddo
      enddo
    endif

  endif

  return

9999 continue
  call fcpsb_serror()
  return

end subroutine zaxpby

subroutine  zaxpbyv2(m, n, alpha, X, lldx, beta, Y, lldy, Z, lldz, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_) :: n, m, lldx, lldy, lldz, info
  complex(psb_dpk_) X(lldx,*), Y(lldy,*), Z(lldy,*)
  complex(psb_dpk_) alpha, beta
  integer(psb_ipk_) :: i, j
  integer(psb_ipk_) :: int_err(5)
  character  name*20
  name='zaxpby'


  !
  !     Error handling
  !
  info = psb_success_
  if (m.lt.0) then
    info=psb_err_iarg_neg_
    int_err(1)=1
    int_err(2)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (n.lt.0) then
    info=psb_err_iarg_neg_
    int_err(1)=1
    int_err(2)=n
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (lldx.lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=5
    int_err(2)=1
    int_err(3)=lldx
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (lldy.lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=8
    int_err(2)=1
    int_err(3)=lldy
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (lldz.lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=8
    int_err(2)=1
    int_err(3)=lldz
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  endif

  if (alpha.eq.zzero) then
    if (beta.eq.zzero) then
      do j=1, n
        do i=1,m
          Z(i,j) = zzero
        enddo
      enddo
    else if (beta.eq.zone) then
      !
      !        Do nothing!
      !

    else if (beta.eq.-zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = - y(i,j)
        enddo
      enddo
    else
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) =  beta*y(i,j)
        enddo
      enddo
    endif

  else if (alpha.eq.zone) then

    if (beta.eq.zzero) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = x(i,j) - y(i,j)
        enddo
      enddo
    else
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = x(i,j) + beta*y(i,j)
        enddo
      enddo
    endif

  else if (alpha.eq.-zone) then

    if (beta.eq.zzero) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = -x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = -x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = -x(i,j) - y(i,j)
        enddo
      enddo
    else
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = -x(i,j) + beta*y(i,j)
        enddo
      enddo
    endif

  else

    if (beta.eq.zzero) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = alpha*x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = alpha*x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = alpha*x(i,j) - y(i,j)
        enddo
      enddo
    else
      do j=1,n
        !$omp parallel do private(i)
        do i=1,m
          Z(i,j) = alpha*x(i,j) + beta*y(i,j)
        enddo
      enddo
    endif

  endif

  return

9999 continue
  call fcpsb_serror()
  return

end subroutine zaxpbyv2

subroutine psi_z_upd_xyz(m,alpha, beta, gamma,delta,x, y, z, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)      :: m
  complex(psb_dpk_), intent (in)       ::  x(:)
  complex(psb_dpk_), intent (inout)    ::  y(:)
  complex(psb_dpk_), intent (inout)    ::  z(:)
  complex(psb_dpk_), intent (in)       :: alpha, beta, gamma, delta
  integer(psb_ipk_), intent(out)     :: info

  integer(psb_ipk_) :: i
  integer(psb_ipk_) :: int_err(5)
  character  name*20
  name='z_upd_xyz'

  info = psb_success_
  if (m.lt.0) then
    info=psb_err_iarg_neg_
    int_err(1)=1
    int_err(2)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (size(x).lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=6
    int_err(2)=1
    int_err(3)=size(x)
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (size(y).lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=7
    int_err(2)=1
    int_err(3)=size(y)
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (size(z).lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=8
    int_err(2)=1
    int_err(3)=size(z)
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  endif
 
  if (beta == zzero) then
    if (gamma == zzero) then
      if (alpha == zzero) then
        if (delta == zzero) then
          !  a 0   b 0 g 0 d 0 
          !$omp parallel do private(i)
          do i=1,m
            y(i) = zzero
            z(i) = zzero
          end do
        else if (delta /= zzero) then
          !  a 0   b 0 g 0 d n 
          !$omp parallel do private(i)
          do i=1,m
            y(i) = zzero
            z(i) = delta*z(i)
          end do
        end if
      else if (alpha /= zzero) then
        if (delta == zzero) then
          !  a n   b 0 g 0 d 0 
          !$omp parallel do private(i)
          do i=1,m
            y(i) = alpha*x(i)
            z(i) = zzero
          end do
        else if (delta /= zzero) then
          !  a n   b 0 g 0 d n 
          !$omp parallel do private(i)
          do i=1,m
            y(i) = alpha*x(i)
            z(i) = delta*z(i)
          end do
          
        end if

      end if

    else  if (gamma /= zzero) then

      if (alpha == zzero) then
      
        if (delta == zzero) then
          !  a 0   b 0 g n d 0
          !$omp parallel do private(i)
          do i=1,m
            y(i) = zzero
            z(i) = zzero  ! gamma*y(i)
          end do
          
        else if (delta /= zzero) then
          !  a 0   b 0 g n d n
          !$omp parallel do private(i)
          do i=1,m
            y(i) = zzero
            z(i) = delta*z(i)
          end do
          
        end if

      else if (alpha /= zzero) then
        
        if (delta == zzero) then
          !  a n   b 0 g n d 0
          !$omp parallel do private(i)
          do i=1,m
            y(i) = alpha*x(i)
            z(i) = gamma*y(i)
          end do
          
        else if (delta /= zzero) then
          !  a n   b 0 g n d n
          !$omp parallel do private(i)
          do i=1,m
            y(i) = alpha*x(i)
            z(i) = gamma*y(i)+delta*z(i)
          end do
          
        end if

      end if

    end if

  else  if (beta /= zzero) then
    
    if (gamma == zzero) then
      if (alpha == zzero) then
        if (delta == zzero) then
          !  a 0   b n g 0 d 0
          !$omp parallel do private(i)
          do i=1,m
            y(i) = beta*y(i)
            z(i) = zzero
          end do
          
        else  if (delta /= zzero) then
          !  a 0   b n g 0 d n
          !$omp parallel do private(i)
          do i=1,m
            y(i) = beta*y(i)
            z(i) = delta*z(i)
          end do
          
        end if

      else  if (alpha /= zzero) then
        if (delta == zzero) then
          !  a n  b n g 0 d 0
          !$omp parallel do private(i)
          do i=1,m
            y(i) = alpha*x(i)+beta*y(i)
            z(i) = zzero
          end do
          
        else if (delta /= zzero) then
          !  a n  b n g 0 d n
          !$omp parallel do private(i)
          do i=1,m
            y(i) = alpha*x(i)+beta*y(i)
            z(i) = delta*z(i)
          end do
          
        end if

      end if
    else  if (gamma /= zzero) then
      if (alpha == zzero) then
        if (delta == zzero) then
          !  a 0  b n g n d 0
          !$omp parallel do private(i)
          do i=1,m
            y(i) = beta*y(i)
            z(i) = gamma*y(i)
          end do
          
        else if (delta /= zzero) then
          !  a 0  b n g n d n
          !$omp parallel do private(i)
          do i=1,m
            y(i) = beta*y(i)
            z(i) = gamma*y(i)+delta*z(i)
          end do

        end if

      else if (alpha /= zzero) then
        if (delta == zzero) then
          !  a n b n g n d 0
          !$omp parallel do private(i)
          do i=1,m
            y(i) = alpha*x(i)+beta*y(i)
            z(i) = gamma*y(i)
          end do
          
        else if (delta /= zzero) then
          !  a n b n g n d n
          !$omp parallel do private(i)
          do i=1,m
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

end subroutine psi_z_upd_xyz

subroutine psi_zxyzw(m,a,b,c,d,e,f,x, y, z,w, info)
  use psb_const_mod
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in)      :: m
  complex(psb_dpk_), intent (in)       :: x(:)
  complex(psb_dpk_), intent (inout)    :: y(:)
  complex(psb_dpk_), intent (inout)    :: z(:)
  complex(psb_dpk_), intent (inout)    :: w(:)
  complex(psb_dpk_), intent (in)       :: a,b,c,d,e,f
  integer(psb_ipk_), intent(out)     :: info

  integer(psb_ipk_) :: i
  integer(psb_ipk_) :: int_err(5)
  character  name*20
  name='z_xyzw'

  info = psb_success_
  if (m.lt.0) then
    info=psb_err_iarg_neg_
    int_err(1)=1
    int_err(2)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (size(x).lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=6
    int_err(2)=1
    int_err(3)=size(x)
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (size(y).lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=7
    int_err(2)=1
    int_err(3)=size(y)
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  else if (size(z).lt.max(1,m)) then
    info=psb_err_iarg_not_gtia_ii_
    int_err(1)=8
    int_err(2)=1
    int_err(3)=size(z)
    int_err(4)=m
    call fcpsb_errpush(info,name,int_err)
    goto 9999
  endif

  if ((a==zzero).or.(b==zzero).or. &
       & (c==zzero).or.(d==zzero).or.&
       & (e==zzero).or.(f==zzero)) then
    write(0,*) 'XYZW assumes  a,b,c,d,e,f are all nonzero'
  else
    !$omp parallel do private(i)
    do i=1,m
      y(i) = a*x(i)+b*y(i)
      z(i) = c*y(i)+d*z(i)
      w(i) = e*z(i)+f*w(i)
    end do
    
  end if

  return

9999 continue
  call fcpsb_serror()
  return

end subroutine psi_zxyzw
