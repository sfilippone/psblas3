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
module psi_gthsct_mod

  interface psi_gth
    module procedure psi_igthm, psi_igthv,&
         & psi_dgthm, psi_dgthv,&
         & psi_zgthm, psi_zgthv
  end interface

  interface psi_sct
    module procedure psi_isctm, psi_isctv,&
         & psi_dsctm, psi_dsctv,&
         & psi_zsctm, psi_zsctv
  end interface

contains

  subroutine psi_dgthm(n,k,idx,x,y)

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

  end subroutine psi_dgthm

  subroutine psi_dgthv(n,idx,x,y)

    implicit none

    integer :: n, idx(:)
    real(kind(1.d0)) :: x(:), y(:)

    ! Locals
    integer :: i, j

    do i=1,n
      y(i)=x(idx(i))
    end do

  end subroutine psi_dgthv


  subroutine psi_dsctm(n,k,idx,x,beta,y)

    implicit none

    integer :: n, k, idx(:)
    real(kind(1.d0)) :: beta, x(:), y(:,:)

    ! Locals
    integer :: i, j, pt

    if (beta.eq.0.d0) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(idx(i),j) = x(pt)
        end do
      end do
    else if (beta.eq.1.d0) then
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
  end subroutine psi_dsctm

  subroutine psi_dsctv(n,idx,x,beta,y)

    implicit none

    integer :: n, k, idx(:)
    real(kind(1.d0)) :: beta, x(:), y(:)

    ! Locals
    integer :: i, j, pt

    if (beta.eq.0.d0) then
      do i=1,n
        y(idx(i)) = x(i)
      end do
    else if (beta.eq.1.d0) then
      do i=1,n
        y(idx(i)) = y(idx(i))+x(i)
      end do
    else
      do i=1,n
        y(idx(i)) = beta*y(idx(i))+x(i)
      end do
    end if
  end subroutine psi_dsctv


  subroutine psi_igthm(n,k,idx,x,y)

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

  end subroutine psi_igthm


  subroutine psi_igthv(n,idx,x,y)

    implicit none

    integer :: n, idx(:)
    integer :: x(:), y(:)

    ! Locals
    integer :: i, j

    do i=1,n
      y(i)=x(idx(i))
    end do

  end subroutine psi_igthv



  subroutine psi_isctm(n,k,idx,x,beta,y)

    implicit none

    integer :: n, k, idx(:)
    integer :: beta, x(:), y(:,:)

    ! Locals
    integer :: i, j, pt

    if (beta.eq.0.d0) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(idx(i),j) = x(pt)
        end do
      end do
    else if (beta.eq.1.d0) then
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
  end subroutine psi_isctm

  subroutine psi_isctv(n,idx,x,beta,y)

    implicit none

    integer :: n, k, idx(:)
    integer :: beta, x(:), y(:)

    ! Locals
    integer :: i, j, pt

    if (beta.eq.0.d0) then
      do i=1,n
        y(idx(i)) = x(i)
      end do
    else if (beta.eq.1.d0) then
      do i=1,n
        y(idx(i)) = y(idx(i))+x(i)
      end do
    else
      do i=1,n
        y(idx(i)) = beta*y(idx(i))+x(i)
      end do
    end if
  end subroutine psi_isctv


  subroutine psi_zgthm(n,k,idx,x,y)

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

  end subroutine psi_zgthm


  subroutine psi_zgthv(n,idx,x,y)

    implicit none

    integer :: n, idx(:)
    complex(kind(1.d0)) :: x(:), y(:)

    ! Locals
    integer :: i, j

    do i=1,n
      y(i)=x(idx(i))
    end do

  end subroutine psi_zgthv

  subroutine psi_zsctm(n,k,idx,x,beta,y)

    implicit none

    integer :: n, k, idx(:)
    complex(kind(1.d0)) :: beta, x(:), y(:,:)

    ! Locals
    integer :: i, j, pt

    if (beta.eq.0.d0) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(idx(i),j) = x(pt)
        end do
      end do
    else if (beta.eq.1.d0) then
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
  end subroutine psi_zsctm


  subroutine psi_zsctv(n,idx,x,beta,y)

    implicit none

    integer :: n, k, idx(:)
    complex(kind(1.d0)) :: beta, x(:), y(:)

    ! Locals
    integer :: i, j, pt

    if (beta.eq.0.d0) then
      do i=1,n
        y(idx(i)) = x(i)
      end do
    else if (beta.eq.1.d0) then
      do i=1,n
        y(idx(i)) = y(idx(i))+x(i)
      end do
    else
      do i=1,n
        y(idx(i)) = beta*y(idx(i))+x(i)
      end do
    end if
  end subroutine psi_zsctv

end module psi_gthsct_mod
