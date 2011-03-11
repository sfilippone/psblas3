subroutine psb_s_aspxpby(alpha, nx, ix, x, beta, y, info)
  use psb_const_mod
  integer, intent(in)               :: nx
  integer, intent(in)               :: ix(:)
  real(psb_spk_), intent (in)       :: x(:)
  real(psb_spk_), intent (inout)    :: y(:)
  real(psb_spk_), intent (in)       :: alpha, beta
  integer, intent(out)              :: info
  integer :: i, ip 

  info=psb_success_

  if (nx > max(size(ix),size(x))) then 
    info = -2
    return
  end if

  if (beta /= sone) then 
    if (beta == -sone) then 
      y(:) = -y(:)
    else if (beta == szero) then 
      y(:) = szero
    else
      y(:) = beta * y(:)
    end if
  end if

  if (alpha == szero) return 

  if (alpha == sone) then 
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) + x(i)
    end do
  else   if (alpha == -sone) then 
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) - x(i)
    end do
  else
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) + alpha*x(i)
    end do
  end if


end subroutine psb_s_aspxpby

subroutine psb_d_aspxpby(alpha, nx, ix, x, beta, y, info)
  use psb_const_mod
  integer, intent(in)               :: nx
  integer, intent(in)               :: ix(:)
  real(psb_dpk_), intent (in)       :: x(:)
  real(psb_dpk_), intent (inout)    :: y(:)
  real(psb_dpk_), intent (in)       :: alpha, beta
  integer, intent(out)              :: info
  integer :: i, ip 

  info=psb_success_
  
  if (nx > max(size(ix),size(x))) then 
    info = -2
    return
  end if

  if (beta /= done) then 
    if (beta == -done) then 
      y(:) = -y(:)
    else if (beta == dzero) then 
      y(:) = dzero
    else
      y(:) = beta * y(:)
    end if
  end if

  if (alpha == dzero) return 

  if (alpha == done) then 
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) + x(i)
    end do
  else   if (alpha == -done) then 
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) - x(i)
    end do
  else
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) + alpha*x(i)
    end do
  end if


end subroutine psb_d_aspxpby
subroutine psb_c_aspxpby(alpha, nx, ix, x, beta, y, info)
  use psb_const_mod
  integer, intent(in)               :: nx
  integer, intent(in)               :: ix(:)
  complex(psb_spk_), intent (in)    :: x(:)
  complex(psb_spk_), intent (inout) :: y(:)
  complex(psb_spk_), intent (in)    :: alpha, beta
  integer, intent(out)              :: info
  integer :: i, ip 

  info=psb_success_
  
  if (nx > max(size(ix),size(x))) then 
    info = -2
    return
  end if
  
  if (beta /= cone) then 
    if (beta == -cone) then 
      y(:) = -y(:)
    else if (beta == czero) then 
      y(:) = czero
    else
      y(:) = beta * y(:)
    end if
  end if

  if (alpha == czero) return 

  if (alpha == cone) then 
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) + x(i)
    end do
  else   if (alpha == -cone) then 
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) - x(i)
    end do
  else
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) + alpha*x(i)
    end do
  end if


end subroutine psb_c_aspxpby

subroutine psb_z_aspxpby(alpha, nx, ix, x, beta, y, info)
  use psb_const_mod
  integer, intent(in)               :: nx
  integer, intent(in)               :: ix(:)
  complex(psb_dpk_), intent (in)    :: x(:)
  complex(psb_dpk_), intent (inout) :: y(:)
  complex(psb_dpk_), intent (in)    :: alpha, beta
  integer, intent(out)              :: info
  integer :: i, ip 

  info=psb_success_
  
  if (nx > max(size(ix),size(x))) then 
    info = -2
    return
  end if
  
  if (beta /= zone) then 
    if (beta == -zone) then 
      y(:) = -y(:)
    else if (beta == zzero) then 
      y(:) = zzero
    else
      y(:) = beta * y(:)
    end if
  end if

  if (alpha == zzero) return 

  if (alpha == zone) then 
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) + x(i)
    end do
  else if (alpha == -zone) then 
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) - x(i)
    end do
  else
    do i=1, nx
      ip = ix(i) 
      y(ip) = y(ip) + alpha*x(i)
    end do
  end if

end subroutine psb_z_aspxpby
