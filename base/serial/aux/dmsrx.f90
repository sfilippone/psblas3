!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
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
! File:  imsrx.f90 
 ! Subroutine: 
 ! Parameters:
subroutine dmsrx(n,x,indx,idir,flag)
  use psb_serial_mod
  implicit none
  integer :: n,idir,flag
  real(psb_dpk_) :: x(n)
  integer :: indx(n)

  integer, allocatable :: iaux(:)

  integer :: iret, info, lp, k,lswap, ixswap
  real(psb_dpk_) :: swap

  if (n<0) then 
!!$    write(0,*) 'Error: DMSRX: N<0'
    return
  endif

  if (n==0) return

  if (flag == psb_sort_ovw_idx_) then 
    do k=1,n
      indx(k) = k
    enddo
  end if

  if (n==1) return

  allocate(iaux(0:n+1),stat=info)
  if (info/=0) then 
    call psb_errpush(4000,r_name='dmsrx')
    call psb_error()
  endif

  if (idir == psb_sort_up_) then 
    call dmsort_up(n,x,iaux,iret)
  else
    call dmsort_dw(n,x,iaux,iret)
  end if

  if (iret /= 1) then 
    lp = iaux(0)
    k  = 1
    do 
      if ((lp==0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      swap     = x(lp)
      x(lp)    = x(k)
      x(k)     = swap
      ixswap   = indx(lp)
      indx(lp) = indx(k)
      indx(k)  = ixswap
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lp
      lp = lswap 
      k  = k + 1
    enddo
  end if

  deallocate(iaux,stat=info)
  if (info/=0) then 
    call psb_errpush(4000,r_name='dmsrx')
    call psb_error()
  endif
  return
end subroutine dmsrx
