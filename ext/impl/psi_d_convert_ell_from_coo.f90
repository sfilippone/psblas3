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
subroutine psi_d_convert_ell_from_coo(a,tmp,info,hacksize) 

  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psi_d_convert_ell_from_coo
  use psi_ext_util_mod
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: tmp
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_), intent(in), optional    :: hacksize

  !locals
  Integer(Psb_ipk_) :: nza, nr, i,j,k, idl,err_act, nc, nzm, &
       & ir, ic, hsz_, ldv

  info = psb_success_

  nr  = tmp%get_nrows()
  nc  = tmp%get_ncols()
  nza = tmp%get_nzeros()

  hsz_ = 1
  if (present(hacksize)) then
    if (hacksize> 1) hsz_ = hacksize
  end if
  ! Make ldv a multiple of hacksize 
  ldv = ((nr+hsz_-1)/hsz_)*hsz_

  ! If it is sorted then we can lessen memory impact 
  a%psb_d_base_sparse_mat = tmp%psb_d_base_sparse_mat

  ! First compute the number of nonzeros in each row.
  call psb_realloc(nr,a%irn,info) 
  if (info /= psb_success_) return
  a%irn = 0
  do i=1, nza
    ir = tmp%ia(i)
    a%irn(ir) = a%irn(ir) + 1
  end do
  nzm = 0
  a%nzt = 0 
  do i=1,nr
    nzm = max(nzm,a%irn(i))
    a%nzt = a%nzt + a%irn(i)
  end do
  ! Allocate and extract.
  call psb_realloc(nr,a%idiag,info) 
  if (info == psb_success_) call psb_realloc(ldv,nzm,a%ja,info) 
  if (info == psb_success_) call psb_realloc(ldv,nzm,a%val,info)
  if (info /= psb_success_) return

  call psi_d_xtr_ell_from_coo(ione,nr,nzm,tmp%ia,tmp%ja,tmp%val,&
       & a%ja,a%val,a%irn,a%idiag,ldv)

end subroutine psi_d_convert_ell_from_coo

