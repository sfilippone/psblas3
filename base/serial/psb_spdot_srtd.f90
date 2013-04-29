!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File:  psb_spdot_srtd 
! Subroutine: X_nspaxpby
!   Computes Z = alpha*X + beta*Y
!    where X,Y and Z are all sparse vectors
!    
! Function: X_spdot_srtd
!   Computes dot product (V1,V2) 
!    where both V1 and V2 are  sparse vectors; assumes the
!    vectors' entries are sorted in increasing order. 
!
!
!
subroutine psb_s_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
  use psb_const_mod
  integer(psb_ipk_), intent(out)              :: nz
  integer(psb_ipk_), intent(out)              :: iz(:)
  real(psb_spk_), intent (out)      :: z(:)
  integer(psb_ipk_), intent(in)               :: nx, ny
  integer(psb_ipk_), intent(in)               :: ix(:), iy(:)
  real(psb_spk_), intent (in)       :: x(:), y(:)
  real(psb_spk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)              :: info

  integer(psb_ipk_) :: i,j,k, ipx, ipy, isz
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
        nz = nz + 1 
        if (nz > isz) then 
          info = -1 
          return
        endif
        iz(nz) = ix(ipx) 
        z(nz)  = acc
        ipx    = ipx + 1
        ipy    = ipy + 1
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
  integer(psb_ipk_), intent(in) :: nv1,nv2
  integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
  real(psb_spk_), intent(in) :: v1(*),v2(*)
  real(psb_spk_)      :: dot

  integer(psb_ipk_) :: i,j,k, ip1, ip2

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
  integer(psb_ipk_), intent(out)              :: nz
  integer(psb_ipk_), intent(out)              :: iz(:)
  real(psb_dpk_), intent (out)      :: z(:)
  integer(psb_ipk_), intent(in)               :: nx, ny
  integer(psb_ipk_), intent(in)               :: ix(:), iy(:)
  real(psb_dpk_), intent (in)       :: x(:), y(:)
  real(psb_dpk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)              :: info

  integer(psb_ipk_) :: i,j,k, ipx, ipy, isz
  real(psb_dpk_) :: acc

  info = psb_success_
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
        nz = nz + 1 
        if (nz > isz) then 
          info = -1 
          return
        endif
        iz(nz) = ix(ipx) 
        z(nz)  = acc
        ipx    = ipx + 1
        ipy    = ipy + 1
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
  integer(psb_ipk_), intent(in) :: nv1,nv2
  integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
  real(psb_dpk_), intent(in) :: v1(*), v2(*)
  real(psb_dpk_)      :: dot

  integer(psb_ipk_) :: i,j,k, ip1, ip2, im1, im2, ix1, ix2

  dot = dzero 
  ip1 = 1
  ip2 = 1
  if (nv1 == 0) return
  if (nv2 == 0) return
  im1 = iv1(nv1)
  im2 = iv2(nv2)
  ix1 = iv1(ip1)
  ix2 = iv2(ip2)
  do 
    if (ix1>im2) exit
    if (ix2>im1) exit

    if (ix1 == ix2) then 
      dot = dot + v1(ip1)*v2(ip2)
      ip1 = ip1 + 1
      if (ip1 > nv1) exit
      ix1 = iv1(ip1)
      ip2 = ip2 + 1
      if (ip2 > nv2) exit
      ix2 = iv2(ip2)
    else if (ix1 < ix2) then 
      ip1 = ip1 + 1 
      if (ip1 > nv1) exit
      ix1 = iv1(ip1)
    else
      ip2 = ip2 + 1 
      if (ip2 > nv2) exit
      ix2 = iv2(ip2)
    end if
  end do

end function psb_d_spdot_srtd

subroutine psb_c_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
  use psb_const_mod
  integer(psb_ipk_), intent(out)              :: nz
  integer(psb_ipk_), intent(out)              :: iz(:)
  complex(psb_spk_), intent (out)      :: z(:)
  integer(psb_ipk_), intent(in)               :: nx, ny
  integer(psb_ipk_), intent(in)               :: ix(:), iy(:)
  complex(psb_spk_), intent (in)       :: x(:), y(:)
  complex(psb_spk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)              :: info

  integer(psb_ipk_) :: i,j,k, ipx, ipy, isz
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
        nz = nz + 1 
        if (nz > isz) then 
          info = -1 
          return
        endif
        iz(nz) = ix(ipx) 
        z(nz)  = acc
        ipx    = ipx + 1
        ipy    = ipy + 1
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
  integer(psb_ipk_), intent(in) :: nv1,nv2
  integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
  complex(psb_spk_), intent(in) :: v1(*),v2(*)
  complex(psb_spk_)      :: dot

  integer(psb_ipk_) :: i,j,k, ip1, ip2

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
  integer(psb_ipk_), intent(out)              :: nz
  integer(psb_ipk_), intent(out)              :: iz(:)
  complex(psb_dpk_), intent (out)      :: z(:)
  integer(psb_ipk_), intent(in)               :: nx, ny
  integer(psb_ipk_), intent(in)               :: ix(:), iy(:)
  complex(psb_dpk_), intent (in)       :: x(:), y(:)
  complex(psb_dpk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)              :: info

  integer(psb_ipk_) :: i,j,k, ipx, ipy, isz
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
        nz = nz + 1 
        if (nz > isz) then 
          info = -1 
          return
        endif
        iz(nz) = ix(ipx) 
        z(nz)  = acc
        ipx    = ipx + 1
        ipy    = ipy + 1
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
  integer(psb_ipk_), intent(in) :: nv1,nv2
  integer(psb_ipk_), intent(in) :: iv1(*), iv2(*)
  complex(psb_dpk_), intent(in) :: v1(*),v2(*)
  complex(psb_dpk_)      :: dot

  integer(psb_ipk_) :: i,j,k, ip1, ip2

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
