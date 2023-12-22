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
subroutine psi_z_xtr_ell_from_coo(i,nr,mxrwl,iac,jac,valc, &
     & ja,val,irn,diag,ld)
  use psb_base_mod, only : psb_ipk_, psb_success_, psb_dpk_, zzero
  
  implicit none 
  integer(psb_ipk_) :: i,nr,mxrwl,ld
  integer(psb_ipk_) :: iac(*),jac(*),ja(ld,*),irn(*),diag(*)
  complex(psb_dpk_)    :: valc(*), val(ld,*)
  
  integer(psb_ipk_) :: ii,jj,kk, kc,nc, ir, ic 
  kc = 1
  do ii = 1, nr
    nc = irn(ii)
    do jj=1,nc
      !if (iac(kc) /= i+ii-1) write(0,*) 'Copy mismatch',iac(kc),i,ii,i+ii-1
      ir = iac(kc)
      ic = jac(kc)
      if (ir == ic) diag(ii) = jj
      ja(ii,jj)  = ic
      val(ii,jj) = valc(kc) 
      kc = kc + 1
    end do
    ! We are assuming that jac contains at least one valid entry
    ! If the previous loop did not have any entries, pick one valid
    ! value.
    if (nc == 0) ic = jac(1)
    do jj = nc+1,mxrwl
      ja(ii,jj)  = ic
      val(ii,jj) = zzero
    end do
  end do
end subroutine psi_z_xtr_ell_from_coo

