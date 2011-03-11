subroutine psb_s_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
  use psb_const_mod
  integer, intent(out)              :: nz
  integer, intent(out)              :: iz(:)
  real(psb_spk_), intent (out)      :: z(:)
  integer, intent(in)               :: nx, ny
  integer, intent(in)               :: ix(:), iy(:)
  real(psb_spk_), intent (in)       :: x(:), y(:)
  real(psb_spk_), intent (in)       :: alpha, beta
  integer, intent(out)              :: info

  integer        :: i,j,k, ipx, ipy, isz
  real(psb_spk_) :: acc

  info=psb_success_
  nz   = 0
  ipx  = 1
  ipy  = 1 
  isz  = min(size(iz),size(z))

  if (beta == szero) then 
    if (alpha == szero) return 
    nz = nx
    if (nz > isz) then 
      info = -1 
      return
    endif
    iz(1:nz) = ix(1:nx)
    z(1:nz)  = alpha*x(1:nx) 

  else if (alpha == szero) then

    nz = ny
    if (nz > isz) then 
      info = -1 
      return
    endif
    iz(1:nz) = iy(1:ny)
    z(1:nz)  = beta*y(1:ny) 

  else
    
    do 
      if (ipx > nx) exit
      if (ipy > ny) exit
      if (ix(ipx) == iy(ipy)) then 
        acc = beta*y(ipy) + alpha*x(ipx) 
        if (acc /= szero) then 
          nz = nz + 1 
          if (nz > isz) then 
            info = -1 
            return
          endif
          iz(nz) = ix(ipx) 
          z(nz)  = acc
          ipx    = ipx + 1
          ipy    = ipy + 1
        end if
      else 
        nz = nz + 1 
        if (nz > isz) then 
          info = -1 
          return
        endif
        if (ix(ipx) < iy(ipy)) then 
          iz(nz) = ix(ipx) 
          z(nz)  = alpha*x(ipx) 
          ipx    = ipx + 1
        else
          iz(nz) = iy(ipy) 
          z(nz)  = beta*y(ipy) 
          ipy    = ipy + 1
        end if
      end if
    end do
    do 
      if (ipx > nx) exit
      nz = nz + 1 
      if (nz > isz) then 
        info = -1 
        return
      endif
      iz(nz) = ix(ipx) 
      z(nz)  = alpha*x(ipx) 
      ipx    = ipx + 1
    end do
    do 
      if (ipy > ny) exit
      nz = nz + 1 
      if (nz > isz) then 
        info = -1 
        return
      endif
      iz(nz) = iy(ipy) 
      z(nz)  = beta*y(ipy) 
      ipy    = ipy + 1
    end do
  end if

end subroutine psb_s_nspaxpby

function psb_s_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
  use psb_const_mod
  integer, intent(in) :: nv1,nv2
  integer, intent(in) :: iv1(*), iv2(*)
  real(psb_spk_), intent(in) :: v1(*),v2(*)
  real(psb_spk_)      :: dot

  integer :: i,j,k, ip1, ip2

  dot = szero 
  ip1 = 1
  ip2 = 1

  do 
    if (ip1 > nv1) exit
    if (ip2 > nv2) exit
    if (iv1(ip1) == iv2(ip2)) then 
      dot = dot + v1(ip1)*v2(ip2)
      ip1 = ip1 + 1
      ip2 = ip2 + 1
    else if (iv1(ip1) < iv2(ip2)) then 
      ip1 = ip1 + 1 
    else
      ip2 = ip2 + 1 
    end if
  end do

end function psb_s_spdot_srtd

subroutine psb_d_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
  use psb_const_mod
  integer, intent(out)              :: nz
  integer, intent(out)              :: iz(:)
  real(psb_dpk_), intent (out)      :: z(:)
  integer, intent(in)               :: nx, ny
  integer, intent(in)               :: ix(:), iy(:)
  real(psb_dpk_), intent (in)       :: x(:), y(:)
  real(psb_dpk_), intent (in)       :: alpha, beta
  integer, intent(out)              :: info

  integer        :: i,j,k, ipx, ipy, isz
  real(psb_dpk_) :: acc

  info=psb_success_
  nz   = 0
  ipx  = 1
  ipy  = 1 
  isz  = min(size(iz),size(z))

  if (beta == dzero) then 
    if (alpha == dzero) return 
    nz = nx
    if (nz > isz) then 
      info = -1 
      return
    endif
    iz(1:nz) = ix(1:nx)
    z(1:nz)  = alpha*x(1:nx) 

  else if (alpha == dzero) then

    nz = ny
    if (nz > isz) then 
      info = -1 
      return
    endif
    iz(1:nz) = iy(1:ny)
    z(1:nz)  = beta*y(1:ny) 

  else
    
    do 
      if (ipx > nx) exit
      if (ipy > ny) exit
      if (ix(ipx) == iy(ipy)) then 
        acc = beta*y(ipy) + alpha*x(ipx) 
        if (acc /= dzero) then 
          nz = nz + 1 
          if (nz > isz) then 
            info = -1 
            return
          endif
          iz(nz) = ix(ipx) 
          z(nz)  = acc
          ipx    = ipx + 1
          ipy    = ipy + 1
        end if
      else 
        nz = nz + 1 
        if (nz > isz) then 
          info = -1 
          return
        endif
        if (ix(ipx) < iy(ipy)) then 
          iz(nz) = ix(ipx) 
          z(nz)  = alpha*x(ipx) 
          ipx    = ipx + 1
        else
          iz(nz) = iy(ipy) 
          z(nz)  = beta*y(ipy) 
          ipy    = ipy + 1
        end if
      end if
    end do
    do 
      if (ipx > nx) exit
      nz = nz + 1 
      if (nz > isz) then 
        info = -1 
        return
      endif
      iz(nz) = ix(ipx) 
      z(nz)  = alpha*x(ipx) 
      ipx    = ipx + 1
    end do
    do 
      if (ipy > ny) exit
      nz = nz + 1 
      if (nz > isz) then 
        info = -1 
        return
      endif
      iz(nz) = iy(ipy) 
      z(nz)  = beta*y(ipy) 
      ipy    = ipy + 1
    end do
  end if

end subroutine psb_d_nspaxpby



function psb_d_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
  use psb_const_mod
  integer, intent(in) :: nv1,nv2
  integer, intent(in) :: iv1(*), iv2(*)
  real(psb_dpk_), intent(in) :: v1(*),v2(*)
  real(psb_dpk_)      :: dot

  integer :: i,j,k, ip1, ip2

  dot = dzero 
  ip1 = 1
  ip2 = 1
  if (nv1 == 0) return
  if (nv2 == 0) return
  do 
    if (iv1(ip1) == iv2(ip2)) then 
      dot = dot + v1(ip1)*v2(ip2)
      ip1 = ip1 + 1
      ip2 = ip2 + 1
    else if (iv1(ip1) < iv2(ip2)) then 
      ip1 = ip1 + 1 
    else
      ip2 = ip2 + 1 
    end if
    if ((ip1 > nv1) .or. (ip2 > nv2)) exit
  end do

end function psb_d_spdot_srtd

subroutine psb_c_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
  use psb_const_mod
  integer, intent(out)              :: nz
  integer, intent(out)              :: iz(:)
  complex(psb_spk_), intent (out)      :: z(:)
  integer, intent(in)               :: nx, ny
  integer, intent(in)               :: ix(:), iy(:)
  complex(psb_spk_), intent (in)       :: x(:), y(:)
  complex(psb_spk_), intent (in)       :: alpha, beta
  integer, intent(out)              :: info

  integer        :: i,j,k, ipx, ipy, isz
  complex(psb_spk_) :: acc

  info=psb_success_
  nz   = 0
  ipx  = 1
  ipy  = 1 
  isz  = min(size(iz),size(z))

  if (beta == czero) then 
    if (alpha == czero) return 
    nz = nx
    if (nz > isz) then 
      info = -1 
      return
    endif
    iz(1:nz) = ix(1:nx)
    z(1:nz)  = alpha*x(1:nx) 

  else if (alpha == czero) then

    nz = ny
    if (nz > isz) then 
      info = -1 
      return
    endif
    iz(1:nz) = iy(1:ny)
    z(1:nz)  = beta*y(1:ny) 

  else
    
    do 
      if (ipx > nx) exit
      if (ipy > ny) exit
      if (ix(ipx) == iy(ipy)) then 
        acc = beta*y(ipy) + alpha*x(ipx) 
        if (acc /= czero) then 
          nz = nz + 1 
          if (nz > isz) then 
            info = -1 
            return
          endif
          iz(nz) = ix(ipx) 
          z(nz)  = acc
          ipx    = ipx + 1
          ipy    = ipy + 1
        end if
      else 
        nz = nz + 1 
        if (nz > isz) then 
          info = -1 
          return
        endif
        if (ix(ipx) < iy(ipy)) then 
          iz(nz) = ix(ipx) 
          z(nz)  = alpha*x(ipx) 
          ipx    = ipx + 1
        else
          iz(nz) = iy(ipy) 
          z(nz)  = beta*y(ipy) 
          ipy    = ipy + 1
        end if
      end if
    end do
    do 
      if (ipx > nx) exit
      nz = nz + 1 
      if (nz > isz) then 
        info = -1 
        return
      endif
      iz(nz) = ix(ipx) 
      z(nz)  = alpha*x(ipx) 
      ipx    = ipx + 1
    end do
    do 
      if (ipy > ny) exit
      nz = nz + 1 
      if (nz > isz) then 
        info = -1 
        return
      endif
      iz(nz) = iy(ipy) 
      z(nz)  = beta*y(ipy) 
      ipy    = ipy + 1
    end do
  end if

end subroutine psb_c_nspaxpby

function psb_c_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
  use psb_const_mod
  integer, intent(in) :: nv1,nv2
  integer, intent(in) :: iv1(*), iv2(*)
  complex(psb_spk_), intent(in) :: v1(*),v2(*)
  complex(psb_spk_)      :: dot

  integer :: i,j,k, ip1, ip2

  dot = czero 
  ip1 = 1
  ip2 = 1

  do 
    if (ip1 > nv1) exit
    if (ip2 > nv2) exit
    if (iv1(ip1) == iv2(ip2)) then 
      dot = dot + conjg(v1(ip1))*v2(ip2)
      ip1 = ip1 + 1
      ip2 = ip2 + 1
    else if (iv1(ip1) < iv2(ip2)) then 
      ip1 = ip1 + 1 
    else
      ip2 = ip2 + 1 
    end if
  end do

end function psb_c_spdot_srtd

subroutine psb_z_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
  use psb_const_mod
  integer, intent(out)              :: nz
  integer, intent(out)              :: iz(:)
  complex(psb_dpk_), intent (out)      :: z(:)
  integer, intent(in)               :: nx, ny
  integer, intent(in)               :: ix(:), iy(:)
  complex(psb_dpk_), intent (in)       :: x(:), y(:)
  complex(psb_dpk_), intent (in)       :: alpha, beta
  integer, intent(out)              :: info

  integer        :: i,j,k, ipx, ipy, isz
  complex(psb_dpk_) :: acc

  info=psb_success_
  nz   = 0
  ipx  = 1
  ipy  = 1 
  isz  = min(size(iz),size(z))

  if (beta == zzero) then 
    if (alpha == zzero) return 
    nz = nx
    if (nz > isz) then 
      info = -1 
      return
    endif
    iz(1:nz) = ix(1:nx)
    z(1:nz)  = alpha*x(1:nx) 

  else if (alpha == zzero) then

    nz = ny
    if (nz > isz) then 
      info = -1 
      return
    endif
    iz(1:nz) = iy(1:ny)
    z(1:nz)  = beta*y(1:ny) 

  else
    
    do 
      if (ipx > nx) exit
      if (ipy > ny) exit
      if (ix(ipx) == iy(ipy)) then 
        acc = beta*y(ipy) + alpha*x(ipx) 
        if (acc /= zzero) then 
          nz = nz + 1 
          if (nz > isz) then 
            info = -1 
            return
          endif
          iz(nz) = ix(ipx) 
          z(nz)  = acc
          ipx    = ipx + 1
          ipy    = ipy + 1
        end if
      else 
        nz = nz + 1 
        if (nz > isz) then 
          info = -1 
          return
        endif
        if (ix(ipx) < iy(ipy)) then 
          iz(nz) = ix(ipx) 
          z(nz)  = alpha*x(ipx) 
          ipx    = ipx + 1
        else
          iz(nz) = iy(ipy) 
          z(nz)  = beta*y(ipy) 
          ipy    = ipy + 1
        end if
      end if
    end do
    do 
      if (ipx > nx) exit
      nz = nz + 1 
      if (nz > isz) then 
        info = -1 
        return
      endif
      iz(nz) = ix(ipx) 
      z(nz)  = alpha*x(ipx) 
      ipx    = ipx + 1
    end do
    do 
      if (ipy > ny) exit
      nz = nz + 1 
      if (nz > isz) then 
        info = -1 
        return
      endif
      iz(nz) = iy(ipy) 
      z(nz)  = beta*y(ipy) 
      ipy    = ipy + 1
    end do
  end if

end subroutine psb_z_nspaxpby

function psb_z_spdot_srtd(nv1,iv1,v1,nv2,iv2,v2) result(dot) 
  use psb_const_mod
  integer, intent(in) :: nv1,nv2
  integer, intent(in) :: iv1(*), iv2(*)
  complex(psb_dpk_), intent(in) :: v1(*),v2(*)
  complex(psb_dpk_)      :: dot

  integer :: i,j,k, ip1, ip2

  dot = zzero 
  ip1 = 1
  ip2 = 1

  do 
    if (ip1 > nv1) exit
    if (ip2 > nv2) exit
    if (iv1(ip1) == iv2(ip2)) then 
      dot = dot + conjg(v1(ip1))*v2(ip2)
      ip1 = ip1 + 1
      ip2 = ip2 + 1
    else if (iv1(ip1) < iv2(ip2)) then 
      ip1 = ip1 + 1 
    else
      ip2 = ip2 + 1 
    end if
  end do

end function psb_z_spdot_srtd
