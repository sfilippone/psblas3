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
subroutine psi_c_convert_hll_from_coo(a,hksz,tmp,info)
  use psb_base_mod
  use psb_c_hll_mat_mod, psb_protect_name => psi_c_convert_hll_from_coo
  use psi_ext_util_mod
  implicit none 
  class(psb_c_hll_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(in)    :: tmp
  integer(psb_ipk_), intent(in)               :: hksz
  integer(psb_ipk_), intent(out)             :: info

  !locals
  Integer(Psb_ipk_)   :: nza, nr, i,j,irw, idl,err_act, nc, isz,irs
  integer(psb_ipk_)   :: nzm, ir, ic, k, hk, mxrwl, noffs, kc


  if (.not.tmp%is_by_rows()) then 
    info = -98765
    return
  end if


  nr  = tmp%get_nrows()
  nc  = tmp%get_ncols()
  nza = tmp%get_nzeros()
  ! If it is sorted then we can lessen memory impact 
  a%psb_c_base_sparse_mat = tmp%psb_c_base_sparse_mat

  ! First compute the number of nonzeros in each row.
  call psb_realloc(nr,a%irn,info) 
  if (info /= 0) return
  a%irn = 0
  do i=1, nza
    a%irn(tmp%ia(i)) = a%irn(tmp%ia(i)) + 1
  end do

  a%nzt = nza
  ! Second. Figure out the block offsets. 
  call a%set_hksz(hksz)
  noffs = (nr+hksz-1)/hksz
  call psb_realloc(noffs+1,a%hkoffs,info) 
  if (info /= 0) return
  a%hkoffs(1) = 0
  j=1
  do i=1,nr,hksz
    ir    = min(hksz,nr-i+1) 
    mxrwl = a%irn(i)
    do k=1,ir-1
      mxrwl = max(mxrwl,a%irn(i+k))
    end do
    a%hkoffs(j+1) = a%hkoffs(j) + mxrwl*hksz
    j = j + 1 
  end do

  !
  ! At this point a%hkoffs(noffs+1) contains the allocation
  ! size a%ja a%val. 
  ! 
  isz = a%hkoffs(noffs+1)
  call psb_realloc(nr,a%idiag,info) 
  if (info == 0) call psb_realloc(isz,a%ja,info) 
  if (info == 0) call psb_realloc(isz,a%val,info) 
  if (info /= 0) return
  ! Init last chunk of data
  nzm = a%hkoffs(noffs+1)-a%hkoffs(noffs)
  a%val(isz-(nzm-1):isz) = czero
  a%ja(isz-(nzm-1):isz)  = nr
  !
  ! Now copy everything, noting the position of the diagonal. 
  !
  kc = 1 
  k  = 1
  do i=1, nr,hksz
    ir    = min(hksz,nr-i+1) 
    irs   = (i-1)/hksz
    hk    = irs + 1
    isz   = (a%hkoffs(hk+1)-a%hkoffs(hk))
    mxrwl = isz/hksz
    nza   = sum(a%irn(i:i+ir-1))
    call psi_c_xtr_ell_from_coo(i,ir,mxrwl,tmp%ia(kc:kc+nza-1),&
         & tmp%ja(kc:kc+nza-1),tmp%val(kc:kc+nza-1),&
         & a%ja(k:k+isz-1),a%val(k:k+isz-1),a%irn(i:i+ir-1),&
         & a%idiag(i:i+ir-1),hksz)
    k  = k + isz
    kc = kc + nza

  enddo

  ! Third copy the other stuff
  if (info /= 0) return
  call a%set_sorted(.true.)

end subroutine psi_c_convert_hll_from_coo
