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
  

subroutine psb_c_cp_hll_to_coo(a,b,info) 
  
  use psb_base_mod
  use psb_c_hll_mat_mod, psb_protect_name => psb_c_cp_hll_to_coo
  implicit none 

  class(psb_c_hll_sparse_mat), intent(in)    :: a
  class(psb_c_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)             :: info

  !locals
  Integer(Psb_ipk_)   :: nza, nr, nc,i,j, jj,k,ir, isz,err_act,  hksz, hk, mxrwl,&
       & irs, nzblk, kc
  integer(psb_ipk_)   :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  if (a%is_dev()) call a%sync()
  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%alloc(nr,nc,nza)
  b%psb_c_base_sparse_mat = a%psb_c_base_sparse_mat
  
  j    = 1
  kc   = 1 
  k    = 1
  hksz = a%hksz
  do i=1, nr,hksz
    ir    = min(hksz,nr-i+1) 
    irs   = (i-1)/hksz
    hk    = irs + 1
    isz   = (a%hkoffs(hk+1)-a%hkoffs(hk))
    nzblk = sum(a%irn(i:i+ir-1))
    call inner_copy(i,ir,b%ia(kc:kc+nzblk-1),&
         & b%ja(kc:kc+nzblk-1),b%val(kc:kc+nzblk-1),&
         & a%ja(k:k+isz-1),a%val(k:k+isz-1),a%irn(i:i+ir-1),&
         & hksz)
    k  = k + isz
    kc = kc + nzblk
    
  enddo

  call b%set_nzeros(nza)
  call b%set_host()
  call b%fix(info)

contains

  subroutine  inner_copy(i,ir,iac,&
       & jac,valc,ja,val,irn,ld)
    integer(psb_ipk_) :: i,ir,ld
    integer(psb_ipk_) :: iac(*),jac(*),ja(ld,*),irn(*)
    complex(psb_spk_)    :: valc(*), val(ld,*)
    
    integer(psb_ipk_) :: ii,jj,kk, kc,nc
    kc = 1
    do ii = 1, ir
      nc = irn(ii)
      do jj=1,nc
        iac(kc)  = i+ii-1
        jac(kc)  = ja(ii,jj) 
        valc(kc) = val(ii,jj) 
        kc = kc + 1
      end do
    end do
        
  end subroutine inner_copy

end subroutine psb_c_cp_hll_to_coo
