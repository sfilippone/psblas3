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

subroutine dsr(n,x,dir)
  use psb_serial_mod
  implicit none
  !
  !  Quicksort.
  !  Adapted from a number of sources, including Don Knuth's TAOCP.
  !
  !     .. Scalar Arguments ..
  integer, intent(in) :: n, dir 
  real(psb_dpk_) ::  x(n)
  !     ..
  !     .. Local Scalars ..
  real(psb_dpk_) :: piv, xt, xk
  integer i, j, ilx, iux, istp, lpiv
  integer n1, n2

  integer, parameter :: maxstack=64,nparms=3,ithrs=16
  integer            :: istack(nparms,maxstack)
  !     ..

  !
  !     small inputs will only get through insertion sort. 
  !
  select case(dir) 

  case(psb_sort_up_)

    if (n > ithrs) then          
      !
      !     Init stack pointer
      !
      istp = 1
      istack(1,istp) = 1
      istack(2,istp) = n

      do 
        if (istp <= 0) exit
        ilx  = istack(1,istp)
        iux  = istack(2,istp)
        istp = istp - 1
        !
        !       Choose a pivot with median-of-three heuristics, leave it 
        !       in the LPIV location
        !            
        i = ilx
        j = iux 
        lpiv = (i+j)/2
        piv  = x(lpiv)
        if (piv < x(i)) then
          xt = x(i)
          x(i) = x(lpiv)
          x(lpiv) = xt
          piv = x(lpiv)
        endif
        if (piv > x(j)) then
          xt = x(j)
          x(j) = x(lpiv)
          x(lpiv) = xt
          piv = x(lpiv)
        endif
        if (piv < x(i)) then
          xt = x(i)
          x(i) = x(lpiv)
          x(lpiv) = xt
          piv = x(lpiv)
        endif
        !
        !     now piv is correct;  place it into first location

        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt

        i = ilx - 1 
        j = iux + 1 

        outer_up: do
          in_up1: do
            i = i + 1
            xk = x(i)
            if (xk >= piv) exit in_up1
          end do in_up1
          !
          !     Ensure finite termination for next loop
          !
          xt  = xk
          x(i) = piv
          in_up2:do 
            j = j - 1
            xk = x(j)
            if (xk <= piv) exit in_up2
          end do in_up2
          x(i) = xt  

          if (j > i) then
            xt  = x(i)
            x(i) = x(j)
            x(j) = xt 
          else
            exit outer_up
          end if
        end do outer_up
        if (i == ilx) then 
          if (x(i) /= piv) then
            write(0,*) 'Should never ever get here????!!!!'
            stop
          endif
          i = i + 1 
        endif

        n1 = (i-1)-ilx+1
        n2 = iux-(i)+1
        if (n1 > n2) then
          if (n1 > ithrs) then 
            istp = istp + 1
            istack(1,istp) = ilx
            istack(2,istp) = i-1
          else
            call disr_up(n1,x(ilx:i-1))
          endif
          if (n2 > ithrs) then
            istp = istp + 1
            istack(1,istp) = i
            istack(2,istp) = iux
          else
            call disr_up(n2,x(i:iux))
          endif
        else
          if (n2 > ithrs) then
            istp = istp + 1
            istack(1,istp) = i
            istack(2,istp) = iux
          else
            call disr_up(n2,x(i:iux))
          endif
          if (n1 > ithrs) then 
            istp = istp + 1
            istack(1,istp) = ilx
            istack(2,istp) = i-1
          else
            call disr_up(n1,x(ilx:i-1))
          endif
        endif
      enddo
    else
      call disr_up(n,x)
    endif

  case(psb_sort_down_) 
   

    if (n > ithrs) then          
      !
      !     Init stack pointer
      !
      istp = 1
      istack(1,istp) = 1
      istack(2,istp) = n

      do 
        if (istp <= 0) exit
        ilx  = istack(1,istp)
        iux  = istack(2,istp)
        istp = istp - 1
        !
        !       Choose a pivot with median-of-three heuristics, leave it 
        !       in the LPIV location
        !            
        i = ilx
        j = iux 
        lpiv = (i+j)/2
        piv  = x(lpiv)
        if (piv > x(i)) then
          xt = x(i)
          x(i) = x(lpiv)
          x(lpiv) = xt
          piv = x(lpiv)
        endif
        if (piv < x(j)) then
          xt = x(j)
          x(j) = x(lpiv)
          x(lpiv) = xt
          piv = x(lpiv)
        endif
        if (piv > x(i)) then
          xt = x(i)
          x(i) = x(lpiv)
          x(lpiv) = xt
          piv = x(lpiv)
        endif
        !
        !     now piv is correct;  place it into first location

        xt = x(i)
        x(i) = x(lpiv)
        x(lpiv) = xt

        i = ilx - 1 
        j = iux + 1 

        outer_dw: do
          in_dw1: do
            i = i + 1
            xk = x(i)
            if (xk <= piv) exit in_dw1
          end do in_dw1
          !
          !     Ensure finite termination for next loop
          !
          xt  = xk
          x(i) = piv
          in_dw2:do 
            j = j - 1
            xk = x(j)
            if (xk >= piv) exit in_dw2
          end do in_dw2
          x(i) = xt  

          if (j > i) then
            xt  = x(i)
            x(i) = x(j)
            x(j) = xt  
          else
            exit outer_dw
          end if
        end do outer_dw
        if (i == ilx) then 
          if (x(i) /= piv) then
            write(0,*) 'Should never ever get here????!!!!'
            stop
          endif
          i = i + 1 
        endif

        n1 = (i-1)-ilx+1
        n2 = iux-(i)+1
        if (n1 > n2) then
          if (n1 > ithrs) then 
            istp = istp + 1
            istack(1,istp) = ilx
            istack(2,istp) = i-1
          else
            call disr_dw(n1,x(ilx:i-1))
          endif
          if (n2 > ithrs) then
            istp = istp + 1
            istack(1,istp) = i
            istack(2,istp) = iux
          else
            call disr_dw(n2,x(i:iux))
          endif
        else
          if (n2 > ithrs) then
            istp = istp + 1
            istack(1,istp) = i
            istack(2,istp) = iux
          else
            call disr_dw(n2,x(i:iux))
          endif
          if (n1 > ithrs) then 
            istp = istp + 1
            istack(1,istp) = ilx
            istack(2,istp) = i-1
          else
            call disr_dw(n1,x(ilx:i-1))
          endif
        endif
      enddo
    else
      call disr_dw(n,x)
    endif

  case default
    write(0,*) 'isr error !',dir
  end select


  return

contains

  subroutine disr_up(n,x)
    implicit none
    integer :: n
    real(psb_dpk_) :: x(n)
    integer :: i,j
    real(psb_dpk_) :: xx

    do j=n-1,1,-1
      if (x(j+1) < x(j)) then
        xx = x(j)
        i=j+1
        do 
          x(i-1) = x(i)
          i = i+1
          if (i>n) exit          
          if (x(i) >= xx) exit
        end do
        x(i-1) = xx
      endif
    enddo
  end subroutine disr_up

  subroutine disr_dw(n,x)
    implicit none
    integer :: n
    real(psb_dpk_) :: x(n)
    integer :: i,j
    real(psb_dpk_) :: xx

    do j=n-1,1,-1
      if (x(j+1) > x(j)) then
        xx = x(j)
        i=j+1
        do 
          x(i-1) = x(i)
          i = i+1
          if (i>n) exit          
          if (x(i) <= xx) exit
        end do
        x(i-1) = xx
      endif
    enddo
  end subroutine disr_dw

end subroutine dsr
