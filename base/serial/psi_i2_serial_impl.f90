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
subroutine psi_i2axpby(m,n,alpha, x, beta, y, info)
  
  use psb_const_mod
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(in)      :: m, n
  integer(psb_i2pk_), intent (in)       ::  x(:,:)
  integer(psb_i2pk_), intent (inout)    ::  y(:,:)
  integer(psb_i2pk_), intent (in)       ::  alpha, beta
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly
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

  if ((m>0).and.(n>0)) call i2axpby(m,n,alpha,x,lx,beta,y,ly,info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine psi_i2axpby

subroutine psi_i2axpbyv(m,alpha, x, beta, y, info)
  
  use psb_const_mod
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(in)      :: m
  integer(psb_i2pk_), intent (in)       ::  x(:)
  integer(psb_i2pk_), intent (inout)    ::  y(:)
  integer(psb_i2pk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly
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

  if (m>0) call i2axpby(m,ione,alpha,x,lx,beta,y,ly,info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psi_i2axpbyv


subroutine psi_i2gthmv(n,k,idx,alpha,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  integer(psb_i2pk_) :: x(:,:), y(:),alpha,beta

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  if (beta == i2zero) then 
    if (alpha == i2zero) then 
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = i2zero
        end do
      end do
    else if (alpha == i2one) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = x(idx(i),j)
        end do
      end do
    else if (alpha == -i2one) then 
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
    if (beta == i2one) then 
      ! Do nothing
    else if (beta == -i2one) then 
      y(1:n*k) = -y(1:n*k) 
    else
      y(1:n*k) = beta*y(1:n*k) 
    end if

    if (alpha == i2zero) then 
      ! do nothing
    else if (alpha == i2one) then 
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = y(pt) + x(idx(i),j)
        end do
      end do
    else if (alpha == -i2one) then
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

end subroutine psi_i2gthmv

subroutine psi_i2gthv(n,idx,alpha,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, idx(:)
  integer(psb_i2pk_) :: x(:), y(:),alpha,beta

  ! Locals
  integer(psb_ipk_) :: i
  if (beta == i2zero) then 
    if (alpha == i2zero) then 
      do i=1,n
        y(i) = i2zero
      end do
    else if (alpha == i2one) then 
      do i=1,n
        y(i) = x(idx(i))
      end do
    else if (alpha == -i2one) then 
      do i=1,n
        y(i) = -x(idx(i))
      end do
    else  
      do i=1,n
        y(i) = alpha*x(idx(i))
      end do
    end if
  else 
    if (beta == i2one) then 
      ! Do nothing
    else if (beta == -i2one) then 
      y(1:n) = -y(1:n) 
    else
      y(1:n) = beta*y(1:n) 
    end if

    if (alpha == i2zero) then 
      ! do nothing
    else if (alpha == i2one) then 
      do i=1,n
        y(i) = y(i) + x(idx(i))
      end do
    else if (alpha == -i2one) then 
      do i=1,n
        y(i) = y(i) - x(idx(i))
      end do
    else  
      do i=1,n
        y(i) = y(i) + alpha*x(idx(i))
      end do
    end if
  end if

end subroutine psi_i2gthv

subroutine psi_i2gthzmm(n,k,idx,x,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  integer(psb_i2pk_) :: x(:,:), y(:,:)

  ! Locals
  integer(psb_ipk_) :: i

  
  do i=1,n
    y(i,1:k)=x(idx(i),1:k)
  end do

end subroutine psi_i2gthzmm

subroutine psi_i2gthzmv(n,k,idx,x,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  integer(psb_i2pk_) :: x(:,:), y(:)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  pt=0
  do j=1,k
    do i=1,n
      pt=pt+1
      y(pt)=x(idx(i),j)
    end do
  end do

end subroutine psi_i2gthzmv

subroutine psi_i2gthzv(n,idx,x,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, idx(:)
  integer(psb_i2pk_) :: x(:), y(:)

  ! Locals
  integer(psb_ipk_) :: i

  do i=1,n
    y(i)=x(idx(i))
  end do

end subroutine psi_i2gthzv

subroutine psi_i2sctmm(n,k,idx,x,beta,y)
  
  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  integer(psb_i2pk_) :: beta, x(:,:), y(:,:)

  ! Locals
  integer(psb_ipk_) :: i, j

  if (beta == i2zero) then
    do i=1,n
      y(idx(i),1:k) = x(i,1:k)
    end do
  else if (beta == i2one) then
    do i=1,n
      y(idx(i),1:k) = y(idx(i),1:k)+x(i,1:k)
    end do
  else
    do i=1,n
      y(idx(i),1:k) = beta*y(idx(i),1:k)+x(i,1:k)
    end do
  end if
end subroutine psi_i2sctmm

subroutine psi_i2sctmv(n,k,idx,x,beta,y)
  
  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  integer(psb_i2pk_) :: beta, x(:), y(:,:)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  if (beta == i2zero) then
    pt=0
    do j=1,k
      do i=1,n
        pt=pt+1
        y(idx(i),j) = x(pt)
      end do
    end do
  else if (beta == i2one) then
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
end subroutine psi_i2sctmv

subroutine psi_i2sctv(n,idx,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, idx(:)
  integer(psb_i2pk_) :: beta, x(:), y(:)

  ! Locals
  integer(psb_ipk_) :: i

  if (beta == i2zero) then
    do i=1,n
      y(idx(i)) = x(i)
    end do
  else if (beta == i2one) then
    do i=1,n
      y(idx(i)) = y(idx(i))+x(i)
    end do
  else
    do i=1,n
      y(idx(i)) = beta*y(idx(i))+x(i)
    end do
  end if
end subroutine psi_i2sctv

subroutine  i2axpby(m, n, alpha, X, lldx, beta, Y, lldy, info)
  use psb_const_mod
  use psb_error_mod
  implicit none 
  integer(psb_ipk_) :: n, m, lldx, lldy, info
  integer(psb_i2pk_) X(lldx,*), Y(lldy,*)
  integer(psb_i2pk_) alpha, beta
  integer(psb_ipk_) :: i, j
  integer(psb_ipk_) :: int_err(5)
  character  name*20
  name='i2axpby'


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

  if (alpha.eq.i2zero) then 
    if (beta.eq.i2zero) then 
      do j=1, n 
        do i=1,m 
          y(i,j) = i2zero
        enddo
      enddo
    else if (beta.eq.i2one) then
      !   
      !        Do nothing! 
      !               

    else if (beta.eq.-i2one) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = - y(i,j)
        enddo
      enddo
    else  
      do j=1,n 
        do i=1,m 
          y(i,j) =  beta*y(i,j)
        enddo
      enddo
    endif

  else if (alpha.eq.i2one) then

    if (beta.eq.i2zero) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = x(i,j)
        enddo
      enddo
    else if (beta.eq.i2one) then
      do j=1,n 
        do i=1,m 
          y(i,j) = x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-i2one) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = x(i,j) - y(i,j)
        enddo
      enddo
    else  
      do j=1,n 
        do i=1,m 
          y(i,j) = x(i,j) + beta*y(i,j)
        enddo
      enddo
    endif

  else if (alpha.eq.-i2one) then 

    if (beta.eq.i2zero) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = -x(i,j)
        enddo
      enddo
    else if (beta.eq.i2one) then
      do j=1,n 
        do i=1,m 
          y(i,j) = -x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-i2one) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = -x(i,j) - y(i,j)
        enddo
      enddo
    else  
      do j=1,n 
        do i=1,m 
          y(i,j) = -x(i,j) + beta*y(i,j)
        enddo
      enddo
    endif

  else  

    if (beta.eq.i2zero) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = alpha*x(i,j)
        enddo
      enddo
    else if (beta.eq.i2one) then
      do j=1,n 
        do i=1,m 
          y(i,j) = alpha*x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-i2one) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = alpha*x(i,j) - y(i,j)
        enddo
      enddo
    else  
      do j=1,n 
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

end subroutine i2axpby
