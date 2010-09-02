!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2010
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
! File:  imsr.f90 
 ! Subroutine: 
 ! Parameters:
subroutine imsr(n,x,idir)
  use psb_serial_mod
  use psb_ip_reord_mod
  implicit none

  integer :: n, idir
  integer :: x(n)
  
  
  integer, allocatable :: iaux(:)
  
  integer :: iswap, iret, info, lp, k
  integer :: lswap

  if (n<0) then 
    return
  endif
  
  if (n<=1) return
  
  allocate(iaux(0:n+1),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_alloc_dealloc_,r_name='imsr')
    call psb_error()
  endif
  
  if (idir == psb_sort_up_) then 
    call msort_up(n,x,iaux,iret)
  else
    call msort_dw(n,x,iaux,iret)
  end if
  
  if (iret == 0) call psb_ip_reord(n,x,iaux)

  deallocate(iaux,stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_alloc_dealloc_,r_name='imsr')
    call psb_error()
  endif
  return
end subroutine imsr
