subroutine psi_dgthmm(n,k,idx,x,y,myrow,icontxt)

  implicit none

  integer :: n, k, idx(:),myrow,icontxt
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

end subroutine psi_dgthmm

subroutine psi_dgthm(n,k,idx,x,y)

  implicit none

  integer :: n, k, idx(:)
  real(kind(1.d0)) :: x(:,:), y(:)

  ! Locals
  integer :: i, j, pt

  write(0,'("Inside gth ",5(i6,2x))')n,k,size(idx),size(x),size(y)
  pt=0
  do j=1,k
     do i=1,n
        pt=pt+1
        y(pt)=x(idx(i),j)
     end do
  end do
  write(0,'("Leaving gth")')

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
