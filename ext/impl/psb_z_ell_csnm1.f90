!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
  

function psb_z_ell_csnm1(a) result(res)
  
  use psb_base_mod
  use psb_z_ell_mat_mod, psb_protect_name => psb_z_ell_csnm1

  implicit none 
  class(psb_z_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_)                          :: res

  integer(psb_ipk_)  :: i,j,k,m,n, nnz, ir, jc, nc, info
  real(psb_dpk_), allocatable :: vt(:)
  logical            :: tra
  Integer(Psb_ipk_)  :: err_act
  character(len=20)  :: name='z_ell_csnm1'
  logical, parameter :: debug=.false.


  if (a%is_dev()) call a%sync()
  res = dzero 
  nnz = a%get_nzeros()
  m   = a%get_nrows()
  n   = a%get_ncols()
  allocate(vt(n),stat=info)
  if (info /= 0) return
  if (a%is_unit()) then 
    vt(:) = done
  else
    vt(:) = dzero
  end if
  do i=1, m
    do j=1,a%irn(i)
      k = a%ja(i,j)
      vt(k) = vt(k) + abs(a%val(i,j))
    end do
  end do
  res = maxval(vt(1:n))
  deallocate(vt,stat=info)

  return

end function psb_z_ell_csnm1
