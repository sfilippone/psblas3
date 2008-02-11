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
module psi_serial_mod
  
  
  interface psi_gth
    module procedure &
         & psi_igthv,  psi_dgthv,  psi_zgthv,&
         & psi_igthzv, psi_dgthzv, psi_zgthzv,&
         & psi_igthzmv, psi_dgthzmv, psi_zgthzmv 
  end interface

  interface psi_sct
    module procedure psi_isctmv, psi_isctv,&
         & psi_dsctmv, psi_dsctv,&
         & psi_zsctmv, psi_zsctv
  end interface


contains


  subroutine psi_igthv(n,idx,alpha,x,beta,y)

    use psb_const_mod
    implicit none

    integer :: n, idx(:)
    integer :: x(:), y(:), alpha, beta

    ! Locals
    integer :: i
    if (beta == izero) then 
      if (alpha == ione) then 
        do i=1,n
          y(i) = x(idx(i))
        end do
      else if (alpha == -ione) then 
        do i=1,n
          y(i) = -x(idx(i))
        end do
      else  
        do i=1,n
          y(i) = alpha*x(idx(i))
        end do
      end if
    else 
      if (beta == ione) then 
        ! Do nothing
      else if (beta == -ione) then 
        y(1:n) = -y(1:n) 
      else
        y(1:n) = beta*y(1:n) 
      end if
        
      if (alpha == ione) then 
        do i=1,n
          y(i) = y(i) + x(idx(i))
        end do
      else if (alpha == -ione) then 
        do i=1,n
          y(i) = y(i) - x(idx(i))
        end do
      else  
        do i=1,n
          y(i) = y(i) + alpha*x(idx(i))
        end do
      end if
    end if

  end subroutine psi_igthv

  subroutine psi_dgthv(n,idx,alpha,x,beta,y)

    use psb_const_mod
    implicit none

    integer :: n, idx(:)
    real(kind(1.d0)) :: x(:), y(:), alpha, beta

    ! Locals
    integer :: i
    if (beta == dzero) then 
      if (alpha == done) then 
        do i=1,n
          y(i) = x(idx(i))
        end do
      else if (alpha == -done) then 
        do i=1,n
          y(i) = -x(idx(i))
        end do
      else 
        do i=1,n
          y(i) = alpha*x(idx(i))
        end do
      end if
    else 
      if (beta == done) then 
        ! Do nothing
      else if (beta == -done) then 
        y(1:n) = -y(1:n) 
      else
        y(1:n) = beta*y(1:n) 
      end if
        
      if (alpha == done) then 
        do i=1,n
          y(i) = y(i) + x(idx(i))
        end do
      else if (alpha == -done) then 
        do i=1,n
          y(i) = y(i) - x(idx(i))
        end do
      else  
        do i=1,n
          y(i) = y(i) + alpha*x(idx(i))
        end do
      end if
    end if

  end subroutine psi_dgthv

  subroutine psi_zgthv(n,idx,alpha,x,beta,y)

    use psb_const_mod
    implicit none

    integer :: n, idx(:)
    complex(kind(1.d0)) :: x(:), y(:),alpha,beta

    ! Locals
    integer :: i
    if (beta == zzero) then 
      if (alpha == zone) then 
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
        
      if (alpha == zone) then 
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




  subroutine psi_dgthzmv(n,k,idx,x,y)

    use psb_const_mod
    implicit none

    integer :: n, k, idx(:)
    real(kind(1.d0)) :: x(:,:), y(:)

    ! Locals
    integer :: i, j, pt

    pt=0
    do j=1,k
      do i=1,n
        pt=pt+1
        y(pt)=x(idx(i),j)
      end do
    end do

  end subroutine psi_dgthzmv


  subroutine psi_igthzmv(n,k,idx,x,y)

    use psb_const_mod
    implicit none

    integer :: n, k, idx(:)
    integer :: x(:,:), y(:)

    ! Locals
    integer :: i, j, pt

    pt=0
    do j=1,k
      do i=1,n
        pt=pt+1
        y(pt)=x(idx(i),j)
      end do
    end do

  end subroutine psi_igthzmv


  subroutine psi_zgthzmv(n,k,idx,x,y)

    use psb_const_mod
    implicit none

    integer :: n, k, idx(:)
    complex(kind(1.d0)) :: x(:,:), y(:)

    ! Locals
    integer :: i, j, pt

    pt=0
    do j=1,k
      do i=1,n
        pt=pt+1
        y(pt)=x(idx(i),j)
      end do
    end do

  end subroutine psi_zgthzmv

  subroutine psi_dgthzv(n,idx,x,y)

    use psb_const_mod
    implicit none

    integer :: n, idx(:)
    real(kind(1.d0)) :: x(:), y(:)

    ! Locals
    integer :: i

    do i=1,n
      y(i)=x(idx(i))
    end do

  end subroutine psi_dgthzv

  subroutine psi_igthzv(n,idx,x,y)

    use psb_const_mod
    implicit none

    integer :: n, idx(:)
    integer :: x(:), y(:)

    ! Locals
    integer :: i

    do i=1,n
      y(i)=x(idx(i))
    end do

  end subroutine psi_igthzv

  subroutine psi_zgthzv(n,idx,x,y)

    use psb_const_mod
    implicit none

    integer :: n, idx(:)
    complex(kind(1.d0)) :: x(:), y(:)

    ! Locals
    integer :: i

    do i=1,n
      y(i)=x(idx(i))
    end do

  end subroutine psi_zgthzv


  subroutine psi_dsctmv(n,k,idx,x,beta,y)

    use psb_const_mod
    implicit none

    integer :: n, k, idx(:)
    real(kind(1.d0)) :: beta, x(:), y(:,:)

    ! Locals
    integer :: i, j, pt

    if (beta == dzero) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(idx(i),j) = x(pt)
        end do
      end do
    else if (beta == done) then
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
  end subroutine psi_dsctmv

  subroutine psi_dsctv(n,idx,x,beta,y)

    use psb_const_mod
    implicit none

    integer :: n, idx(:)
    real(kind(1.d0)) :: beta, x(:), y(:)

    ! Locals
    integer :: i

    if (beta == dzero) then
      do i=1,n
        y(idx(i)) = x(i)
      end do
    else if (beta == done) then
      do i=1,n
        y(idx(i)) = y(idx(i))+x(i)
      end do
    else
      do i=1,n
        y(idx(i)) = beta*y(idx(i))
      end do
      do i=1,n
        y(idx(i)) = y(idx(i))+x(i)
      end do
    end if
  end subroutine psi_dsctv

  subroutine psi_isctmv(n,k,idx,x,beta,y)

    use psb_const_mod
    implicit none

    integer :: n, k, idx(:)
    integer :: beta, x(:), y(:,:)

    ! Locals
    integer :: i, j, pt

    if (beta == izero) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(idx(i),j) = x(pt)
        end do
      end do
    else if (beta == ione) then
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
  end subroutine psi_isctmv

  subroutine psi_isctv(n,idx,x,beta,y)

    use psb_const_mod
    implicit none

    integer :: n, idx(:)
    integer :: beta, x(:), y(:)

    ! Locals
    integer :: i

    if (beta == izero) then
      do i=1,n
        y(idx(i)) = x(i)
      end do
    else if (beta == ione) then
      do i=1,n
        y(idx(i)) = y(idx(i))+x(i)
      end do
    else
      do i=1,n
        y(idx(i)) = beta*y(idx(i))+x(i)
      end do
    end if
  end subroutine psi_isctv

  subroutine psi_zsctmv(n,k,idx,x,beta,y)

    use psb_const_mod
    implicit none

    integer :: n, k, idx(:)
    complex(kind(1.d0)) :: beta, x(:), y(:,:)

    ! Locals
    integer :: i, j, pt

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

    integer :: n, idx(:)
    complex(kind(1.d0)) :: beta, x(:), y(:)

    ! Locals
    integer :: i

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
  

end module psi_serial_mod
