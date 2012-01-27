!
!             Parallel Sparse BLAS  version 3.0
!   (C) Copyright 2006, 2007, 2008, 2009, 2010
!                      Salvatore Filippone    University of Rome Tor Vergata
!                      Alfredo Buttari        CNRS-IRIT, Toulouse
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!   1. Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!   2. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions, and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!   3. The name of the PSBLAS group or the names of its contributors may
!      not be used to endorse or promote products derived from this
!      software without specific written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! File:  msort_dw.f90
!
! Subroutine: msort_dw
!   This subroutine sorts an integer array into ascending order.
!
! Arguments:
!   n         -  integer                   Input: size of the array 
!   k         -  real(*)                  input: array of keys to be sorted
!   l         -  integer(0:n+1)           output: link list 
!   iret      -  integer                  output: 0 Normal termination
!                                                 1 the array was already sorted 
!                                                                     *
! REFERENCES  = (1) D. E. Knuth                                       *
!                   The Art of Computer Programming,                  *
!                     vol.3: Sorting and Searching                    *
!                   Addison-Wesley, 1973                              *
!                                                                     *
!  call msort_dw(n,x,iaux,iret)
!  
!  if (iret == 0) then 
!    lp = iaux(0)
!    k  = 1
!    do 
!      if ((lp == 0).or.(k>n)) exit
!      do 
!        if (lp >= k) exit
!        lp = iaux(lp)
!      end do
!      iswap    = x(lp)
!      x(lp)    = x(k)
!      x(k)     = iswap
!      lswap    = iaux(lp)
!      iaux(lp) = iaux(k)
!      iaux(k)  = lp
!      lp = lswap 
!      k  = k + 1
!    enddo
!  end if
!
!
subroutine smsort_dw(n,k,l,iret)
  use psb_const_mod
  implicit none
  integer(psb_ipk_) :: n, iret
  real(psb_spk_) :: k(n)
  integer(psb_ipk_) :: l(0:n+1)
  !
  integer(psb_ipk_) :: p,q,s,t
  intrinsic iabs,isign
  !     ..
  iret = 0
  !  first step: we are preparing ordered sublists, exploiting
  !  what order was already in the input data; negative links
  !  mark the end of the sublists
  l(0) = 1
  t = n + 1
  do  p = 1,n - 1
    if (k(p) >= k(p+1)) then
      l(p) = p + 1
    else
      l(t) = - (p+1)
      t = p
    end if
  end do
  l(t) = 0
  l(n) = 0
  ! see if the input was already sorted
  if (l(n+1) == 0) then
    iret = 1
    return 
  else
    l(n+1) = iabs(l(n+1))
  end if

  mergepass: do 
    ! otherwise, begin a pass through the list.
    ! throughout all the subroutine we have:
    !  p, q: pointing to the sublists being merged
    !  s: pointing to the most recently processed record
    !  t: pointing to the end of previously completed sublist
    s = 0
    t = n + 1
    p = l(s)
    q = l(t)
    if (q == 0) exit mergepass

    outer: do 

      if (k(p) < k(q)) then 

        l(s) = isign(q,l(s))
        s = q
        q = l(q)
        if (q > 0) then 
          do 
            if (k(p) >= k(q)) cycle outer
            s = q
            q = l(q)
            if (q <= 0) exit
          end do
        end if
        l(s) = p
        s = t
        do 
          t = p
          p = l(p)
          if (p <= 0) exit
        end do

      else 

        l(s) = isign(p,l(s))
        s = p
        p = l(p)
        if (p>0) then 
          do 
            if (k(p) < k(q)) cycle outer 
            s = p
            p = l(p)
            if (p <= 0) exit
          end do
        end if
        !  otherwise, one sublist ended, and we append to it the rest
        !  of the other one.
        l(s) = q
        s = t
        do 
          t = q
          q = l(q)
          if (q <= 0) exit
        end do
      end if

      p = -p
      q = -q
      if (q == 0) then
        l(s) = isign(p,l(s))
        l(t) = 0
        exit outer 
      end if
    end do outer
  end do mergepass

end subroutine smsort_dw
