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
! File:  imsr.f90 
 ! Subroutine: 
 ! Parameters:
subroutine imsr(n,x,idir)
  use psb_serial_mod
  integer :: n, idir
  integer :: x(n)
  
  
  integer, allocatable :: iaux(:)
  
  integer :: iswap, iret, info, lp, k
  integer :: lswap

  if (n<0) then 
!!$    write(0,*) 'Error: IMSR: N<0'
    return
  endif
  
  if (n<=1) return
  
  allocate(iaux(0:n+1),stat=info)
  if (info/=0) then 
    write(0,*) 'IMSR: memory allocation failed',info
    return
  endif
  
  if (idir==psb_sort_up_) then 
    call mrgsrt(n,x,iaux,iret)
  else
    call mrgsrtd(n,x,iaux,iret)
  end if
  
  if (iret == 0) then 
    lp = iaux(0)
    k  = 1
    do 
      if ((lp==0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      iswap    = x(lp)
      x(lp)    = x(k)
      x(k)     = iswap
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lp
      lp = lswap 
      k  = k + 1
    enddo
  end if

  deallocate(iaux,stat=info)
  if (info/=0) then 
    write(0,*) 'IMSR: memory deallocation failed',info
  endif
  return
end subroutine imsr
