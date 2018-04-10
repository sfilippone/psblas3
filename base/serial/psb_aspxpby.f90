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
subroutine psb_s_aspxpby(alpha, nx, ix, x, beta, y, info)
  use psb_const_mod
  integer(psb_ipk_), intent(in)               :: nx
  integer(psb_ipk_), intent(in)               :: ix(:)
  real(psb_spk_), intent (in)       :: x(:)
  real(psb_spk_), intent (inout)    :: y(:)
  real(psb_spk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_) :: i, ip 

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
  integer(psb_ipk_), intent(in)               :: nx
  integer(psb_ipk_), intent(in)               :: ix(:)
  real(psb_dpk_), intent (in)       :: x(:)
  real(psb_dpk_), intent (inout)    :: y(:)
  real(psb_dpk_), intent (in)       :: alpha, beta
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_) :: i, ip 

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
  integer(psb_ipk_), intent(in)               :: nx
  integer(psb_ipk_), intent(in)               :: ix(:)
  complex(psb_spk_), intent (in)    :: x(:)
  complex(psb_spk_), intent (inout) :: y(:)
  complex(psb_spk_), intent (in)    :: alpha, beta
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_) :: i, ip 

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
  integer(psb_ipk_), intent(in)               :: nx
  integer(psb_ipk_), intent(in)               :: ix(:)
  complex(psb_dpk_), intent (in)    :: x(:)
  complex(psb_dpk_), intent (inout) :: y(:)
  complex(psb_dpk_), intent (in)    :: alpha, beta
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_) :: i, ip 

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
