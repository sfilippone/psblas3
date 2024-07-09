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
  

subroutine psb_d_cp_ell_to_coo(a,b,info) 
  
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_cp_ell_to_coo
  implicit none 

  class(psb_d_ell_sparse_mat), intent(in)    :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)             :: info

  !locals
  Integer(Psb_ipk_)   :: i, j, k, nr, nc, nza

  info = psb_success_

  if (a%is_dev()) call a%sync()
  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  b%psb_d_base_sparse_mat = a%psb_d_base_sparse_mat

  k=0
  do i=1, nr
    do j=1,a%irn(i)
      k = k + 1
      b%ia(k)  = i
      b%ja(k)  = a%ja(i,j)
      b%val(k) = a%val(i,j)
    end do
  end do
  call b%set_nzeros(a%get_nzeros())
  call b%fix(info)
  call b%set_host()

end subroutine psb_d_cp_ell_to_coo
