!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
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
! File:  psb_spge_dot.f90 
! Function: X_spge_dot
!
! Compute the dot product between a sparse vector and a dense vector. 
!
function psb_s_spge_dot(nv1,iv1,v1,v2) result(dot) 
  use psb_const_mod
  integer(psb_ipk_), intent(in) :: nv1
  integer(psb_ipk_), intent(in) :: iv1(*)
  real(psb_spk_), intent(in) :: v1(*),v2(*)
  real(psb_spk_)      :: dot

  integer(psb_ipk_) :: i,j,k, ip1, ip2

  dot = szero 
  ip1 = 1

  do i=1, nv1
    dot = dot + v1(i)*v2(iv1(i))
  end do

end function psb_s_spge_dot

function psb_d_spge_dot(nv1,iv1,v1,v2) result(dot) 
  use psb_const_mod
  integer(psb_ipk_), intent(in) :: nv1
  integer(psb_ipk_), intent(in) :: iv1(*)
  real(psb_dpk_), intent(in) :: v1(*),v2(*)
  real(psb_dpk_)      :: dot

  integer(psb_ipk_) :: i,j,k, ip1, ip2

  dot = dzero 

  do i=1, nv1
    dot = dot + v1(i)*v2(iv1(i))
  end do

end function psb_d_spge_dot

function psb_c_spge_dot(nv1,iv1,v1,v2) result(dot) 
  use psb_const_mod
  integer(psb_ipk_), intent(in) :: nv1
  integer(psb_ipk_), intent(in) :: iv1(*)
  complex(psb_spk_), intent(in) :: v1(*),v2(*)
  complex(psb_spk_)      :: dot

  integer(psb_ipk_) :: i,j,k, ip1, ip2

  dot = czero 

  do i=1, nv1
    dot = dot + v1(i)*v2(iv1(i))
  end do

end function psb_c_spge_dot

function psb_z_spge_dot(nv1,iv1,v1,v2) result(dot) 
  use psb_const_mod
  integer(psb_ipk_), intent(in) :: nv1
  integer(psb_ipk_), intent(in) :: iv1(*)
  complex(psb_dpk_), intent(in) :: v1(*),v2(*)
  complex(psb_dpk_)      :: dot

  integer(psb_ipk_) :: i,j,k, ip1, ip2

  dot = zzero 

  do i=1, nv1
    dot = dot + v1(i)*v2(iv1(i))
  end do

end function psb_z_spge_dot
