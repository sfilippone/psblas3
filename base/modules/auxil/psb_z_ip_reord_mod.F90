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
!  Reorder (an) input vector(s) based on a list sort output.
!  Based on: D. E. Knuth: The Art of Computer Programming
!            vol. 3: Sorting and Searching, Addison Wesley, 1973
!            ex. 5.2.12
!
!
module psb_z_ip_reord_mod 
  use psb_const_mod
  
  interface psb_ip_reord
    module procedure psb_ip_reord_z1m,&
         & psb_ip_reord_z1m1, psb_ip_reord_z1m2,&
         & psb_ip_reord_z1m3
    module procedure psb_ip_reord_z1e,&
         & psb_ip_reord_z1e1, psb_ip_reord_z1e2,&
         & psb_ip_reord_z1e3


  end interface

contains
  
  subroutine psb_ip_reord_z1m(n,x,iaux)
    integer(psb_ipk_), intent(in) :: n
    integer(psb_mpk_) :: iaux(0:*) 
    complex(psb_dpk_)   :: x(*)
    integer(psb_mpk_) :: lswap, lp, k
    complex(psb_dpk_) :: swap
    
    lp = iaux(0)
    k  = 1
    do 
      if ((lp == 0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      swap     = x(lp)
      x(lp)    = x(k)
      x(k)     = swap
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lp
      lp = lswap 
      k  = k + 1
    enddo
    return
  end subroutine psb_ip_reord_z1m
  
  subroutine psb_ip_reord_z1m1(n,x,indx,iaux)
    integer(psb_ipk_), intent(in) :: n
    integer(psb_mpk_) :: iaux(0:*) 
    complex(psb_dpk_)   :: x(*)
    integer(psb_mpk_) :: indx(*) 
    integer(psb_mpk_) :: lswap, lp, k, ixswap
    complex(psb_dpk_) :: swap

    lp = iaux(0)
    k  = 1
    do 
      if ((lp == 0).or.(k>n)) exit
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
    return
  end subroutine psb_ip_reord_z1m1
  
  subroutine psb_ip_reord_z1m2(n,x,i1,i2,iaux)
    integer(psb_ipk_), intent(in) :: n
    integer(psb_mpk_) :: iaux(0:*) 
    complex(psb_dpk_)   :: x(*)
    integer(psb_mpk_) :: i1(*), i2(*) 

    
    integer(psb_mpk_) :: lswap, lp, k, isw1, isw2
    complex(psb_dpk_)  :: swap

    lp = iaux(0)
    k  = 1
    do 
      if ((lp == 0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      swap     = x(lp)
      x(lp)    = x(k)
      x(k)     = swap
      isw1     = i1(lp)
      i1(lp)   = i1(k)
      i1(k)    = isw1
      isw2     = i2(lp)
      i2(lp)   = i2(k)
      i2(k)    = isw2
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lp
      lp = lswap 
      k  = k + 1
    enddo
    return
  end subroutine psb_ip_reord_z1m2
  
  subroutine psb_ip_reord_z1m3(n,x,i1,i2,i3,iaux)
    integer(psb_ipk_), intent(in) :: n
    integer(psb_mpk_) :: iaux(0:*) 
    complex(psb_dpk_)   :: x(*)
    integer(psb_mpk_) :: i1(*), i2(*), i3(*) 
    
    integer(psb_mpk_) :: lswap, lp, k, isw1, isw2, isw3
    complex(psb_dpk_)  :: swap

    lp = iaux(0)
    k  = 1
    do 
      if ((lp == 0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      swap     = x(lp)
      x(lp)    = x(k)
      x(k)     = swap
      isw1     = i1(lp)
      i1(lp)   = i1(k)
      i1(k)    = isw1
      isw2     = i2(lp)
      i2(lp)   = i2(k)
      i2(k)    = isw2
      isw3     = i3(lp)
      i3(lp)   = i3(k)
      i3(k)    = isw3
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lp
      lp = lswap 
      k  = k + 1
    enddo
    return
  end subroutine psb_ip_reord_z1m3

    
  subroutine psb_ip_reord_z1e(n,x,iaux)
    integer(psb_ipk_), intent(in) :: n
    integer(psb_epk_) :: iaux(0:*) 
    complex(psb_dpk_)   :: x(*)
    integer(psb_epk_) :: lswap, lp, k
    complex(psb_dpk_) :: swap
    
    lp = iaux(0)
    k  = 1
    do 
      if ((lp == 0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      swap     = x(lp)
      x(lp)    = x(k)
      x(k)     = swap
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lp
      lp = lswap 
      k  = k + 1
    enddo
    return
  end subroutine psb_ip_reord_z1e
  
  subroutine psb_ip_reord_z1e1(n,x,indx,iaux)
    integer(psb_ipk_), intent(in) :: n
    integer(psb_epk_) :: iaux(0:*) 
    complex(psb_dpk_)   :: x(*)
    integer(psb_epk_) :: indx(*) 
    integer(psb_epk_) :: lswap, lp, k, ixswap
    complex(psb_dpk_) :: swap

    lp = iaux(0)
    k  = 1
    do 
      if ((lp == 0).or.(k>n)) exit
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
    return
  end subroutine psb_ip_reord_z1e1
  
  subroutine psb_ip_reord_z1e2(n,x,i1,i2,iaux)
    integer(psb_ipk_), intent(in) :: n
    integer(psb_epk_) :: iaux(0:*) 
    complex(psb_dpk_)   :: x(*)
    integer(psb_epk_) :: i1(*), i2(*) 

    
    integer(psb_epk_) :: lswap, lp, k, isw1, isw2
    complex(psb_dpk_)  :: swap

    lp = iaux(0)
    k  = 1
    do 
      if ((lp == 0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      swap     = x(lp)
      x(lp)    = x(k)
      x(k)     = swap
      isw1     = i1(lp)
      i1(lp)   = i1(k)
      i1(k)    = isw1
      isw2     = i2(lp)
      i2(lp)   = i2(k)
      i2(k)    = isw2
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lp
      lp = lswap 
      k  = k + 1
    enddo
    return
  end subroutine psb_ip_reord_z1e2
  
  subroutine psb_ip_reord_z1e3(n,x,i1,i2,i3,iaux)
    integer(psb_ipk_), intent(in) :: n
    integer(psb_epk_) :: iaux(0:*) 
    complex(psb_dpk_)   :: x(*)
    integer(psb_epk_) :: i1(*), i2(*), i3(*) 
    
    integer(psb_epk_) :: lswap, lp, k, isw1, isw2, isw3
    complex(psb_dpk_)  :: swap

    lp = iaux(0)
    k  = 1
    do 
      if ((lp == 0).or.(k>n)) exit
      do 
        if (lp >= k) exit
        lp = iaux(lp)
      end do
      swap     = x(lp)
      x(lp)    = x(k)
      x(k)     = swap
      isw1     = i1(lp)
      i1(lp)   = i1(k)
      i1(k)    = isw1
      isw2     = i2(lp)
      i2(lp)   = i2(k)
      i2(k)    = isw2
      isw3     = i3(lp)
      i3(lp)   = i3(k)
      i3(k)    = isw3
      lswap    = iaux(lp)
      iaux(lp) = iaux(k)
      iaux(k)  = lp
      lp = lswap 
      k  = k + 1
    enddo
    return
  end subroutine psb_ip_reord_z1e3

end module psb_z_ip_reord_mod
