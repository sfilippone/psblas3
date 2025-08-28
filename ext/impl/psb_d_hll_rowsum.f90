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
  

subroutine psb_d_hll_rowsum(d,a) 
  
  use psb_base_mod
  use psb_d_hll_mat_mod, psb_protect_name => psb_d_hll_rowsum
  implicit none 
  class(psb_d_hll_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)

  integer(psb_ipk_)  :: i,j,k,m,n, nnz, ir, jc, nc, hksz, mxrwl
  logical            :: tra
  Integer(Psb_ipk_)  :: err_act, info, int_err(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = 0
  if (a%is_dev()) call a%sync()

  m = a%get_nrows()
  n = a%get_ncols()
  if (size(d) < m) then 
    info=psb_err_input_asize_small_i_
    int_err(1) = 1
    int_err(2) = size(d)
    int_err(3) = m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if


  if (a%is_unit()) then 
    d =  done
  else
    d = dzero
  end if
  hksz = a%get_hksz()
  j    = 1
  do i=1,m,hksz
    ir    = min(hksz,m-i+1) 
    mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
    k     = a%hkoffs(j) + 1
    call d_hll_rowsum(i,ir,mxrwl,a%irn(i),&
         & a%ja(k),hksz,a%val(k),hksz, &
         & d,info) 
    if (info /= psb_success_) goto 9999
    j = j + 1 
  end do
  
  

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(err_act)
  return

contains

  subroutine d_hll_rowsum(ir,m,n,irn,ja,ldj,val,ldv,&
       & d,info) 
    integer(psb_ipk_), intent(in)   :: ir,m,n,ldj,ldv,ja(ldj,*),irn(*)
    real(psb_dpk_), intent(in)     :: val(ldv,*)
    real(psb_dpk_), intent(inout) :: d(*)
    integer(psb_ipk_), intent(out)  :: info

    integer(psb_ipk_)   :: i,j,k, m4, jc
    real(psb_dpk_) :: acc(4), tmp

    info = psb_success_
    do i=1,m
      do j=1, irn(i)
        d(ir+i-1) = d(ir+i-1) + (val(i,j))
      end do
    end do

  end subroutine d_hll_rowsum

end subroutine psb_d_hll_rowsum
