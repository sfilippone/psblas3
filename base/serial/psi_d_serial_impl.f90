subroutine psi_daxpby(m,n,alpha, x, beta, y, info)
  
  use psb_const_mod
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(in)      :: m, n
  real(psb_dpk_), intent (in)       ::  x(:,:)
  real(psb_dpk_), intent (inout)    ::  y(:,:)
  real(psb_dpk_), intent (in)       ::  alpha, beta
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly
  integer(psb_ipk_) :: ierr(5)
  character(len=20)        :: name, ch_err

  name='psb_geaxpby'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

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
end subroutine psi_daxpby

subroutine psi_daxpbyv(m,alpha, x, beta, y, info)
  
  use psb_const_mod
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(in)      :: m
  real(psb_dpk_), intent (in)       ::  x(:)
  real(psb_dpk_), intent (inout)    ::  y(:)
  real(psb_dpk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: lx, ly
  integer(psb_ipk_) :: ierr(5)
  character(len=20)        :: name, ch_err

  name='psb_geaxpby'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

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

end subroutine psi_daxpbyv


subroutine psi_dgthmv(n,k,idx,alpha,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  real(psb_dpk_) :: x(:,:), y(:),alpha,beta

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  if (beta == dzero) then 
    if (alpha == dzero) then 
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = dzero
        end do
      end do
    else if (alpha == done) then
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = x(idx(i),j)
        end do
      end do
    else if (alpha == -done) then 
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
    if (beta == done) then 
      ! Do nothing
    else if (beta == -done) then 
      y(1:n*k) = -y(1:n*k) 
    else
      y(1:n*k) = beta*y(1:n*k) 
    end if

    if (alpha == dzero) then 
      ! do nothing
    else if (alpha == done) then 
      pt=0
      do j=1,k
        do i=1,n
          pt=pt+1
          y(pt) = y(pt) + x(idx(i),j)
        end do
      end do
    else if (alpha == -done) then
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

end subroutine psi_dgthmv

subroutine psi_dgthv(n,idx,alpha,x,beta,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, idx(:)
  real(psb_dpk_) :: x(:), y(:),alpha,beta

  ! Locals
  integer(psb_ipk_) :: i
  if (beta == dzero) then 
    if (alpha == dzero) then 
      do i=1,n
        y(i) = dzero
      end do
    else if (alpha == done) then 
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

    if (alpha == dzero) then 
      ! do nothing
    else if (alpha == done) then 
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

subroutine psi_dgthzmv(n,k,idx,x,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  real(psb_dpk_) :: x(:,:), y(:)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

  pt=0
  do j=1,k
    do i=1,n
      pt=pt+1
      y(pt)=x(idx(i),j)
    end do
  end do

end subroutine psi_dgthzmv

subroutine psi_dgthzv(n,idx,x,y)

  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, idx(:)
  real(psb_dpk_) :: x(:), y(:)

  ! Locals
  integer(psb_ipk_) :: i

  do i=1,n
    y(i)=x(idx(i))
  end do

end subroutine psi_dgthzv

subroutine psi_dsctmv(n,k,idx,x,beta,y)
  
  use psb_const_mod
  implicit none

  integer(psb_ipk_) :: n, k, idx(:)
  real(psb_dpk_) :: beta, x(:), y(:,:)

  ! Locals
  integer(psb_ipk_) :: i, j, pt

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

  integer(psb_ipk_) :: n, idx(:)
  real(psb_dpk_) :: beta, x(:), y(:)

  ! Locals
  integer(psb_ipk_) :: i

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
      y(idx(i)) = beta*y(idx(i))+x(i)
    end do
  end if
end subroutine psi_dsctv
