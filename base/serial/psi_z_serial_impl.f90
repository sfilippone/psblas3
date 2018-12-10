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

  if (m>0) call zaxpby(m,ione,alpha,x,lx,beta,y,ly,info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psi_zaxpbyv


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

  else if (alpha.eq.zone) then

    if (beta.eq.zzero) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n 
        do i=1,m 
          y(i,j) = x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then 
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

  else if (alpha.eq.-zone) then 

    if (beta.eq.zzero) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = -x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n 
        do i=1,m 
          y(i,j) = -x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then 
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

    if (beta.eq.zzero) then 
      do j=1,n 
        do i=1,m 
          y(i,j) = alpha*x(i,j)
        enddo
      enddo
    else if (beta.eq.zone) then
      do j=1,n 
        do i=1,m 
          y(i,j) = alpha*x(i,j) + y(i,j)
        enddo
      enddo

    else if (beta.eq.-zone) then 
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

end subroutine zaxpby
