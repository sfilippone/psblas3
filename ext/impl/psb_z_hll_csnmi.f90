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
  

function psb_z_hll_csnmi(a) result(res)
  
  use psb_base_mod
  use psb_z_hll_mat_mod, psb_protect_name => psb_z_hll_csnmi
  implicit none 
  class(psb_z_hll_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_)  :: i,j,k,m,n, nr, ir, jc, nc, hksz, mxrwl, info 
  Integer(Psb_ipk_)  :: err_act
  logical            :: is_unit
  character(len=20)  :: name='z_csnmi'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = 0
  res = dzero 
  if (a%is_dev()) call a%sync()
  
  n = a%get_ncols()
  m = a%get_nrows()
  is_unit = a%is_unit()
  hksz = a%get_hksz()
  j=1
  do i=1,m,hksz
    ir    = min(hksz,m-i+1) 
    mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
    k     = a%hkoffs(j) + 1
    call psb_z_hll_csnmi_inner(i,ir,mxrwl,a%irn(i),&
         & a%ja(k),hksz,a%val(k),hksz,&
         & res,is_unit,info) 
    if (info /= psb_success_) goto 9999
    j = j + 1 
  end do
  
  call psb_erractionrestore(err_act)
  return
  

9999 call psb_error_handler(err_act)
  return

contains

  subroutine psb_z_hll_csnmi_inner(ir,m,n,irn,ja,ldj,val,ldv,&
       & res,is_unit,info) 
    integer(psb_ipk_), intent(in)    :: ir,m,n,ldj,ldv,ja(ldj,*),irn(*)
    complex(psb_dpk_), intent(in)      :: val(ldv,*)
    real(psb_dpk_), intent(inout)   :: res
    logical                          :: is_unit
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_) :: i,j,k, m4, jc
    real(psb_dpk_)   :: tmp, acc

    info = psb_success_
    if (is_unit) then 
      tmp = done
    else
      tmp = dzero
    end if
    do i=1,m
      acc = tmp
      do j=1, irn(i)
        acc = acc + abs(val(i,j))
      end do
      res = max(acc,res)
    end do
  end subroutine psb_z_hll_csnmi_inner

end function psb_z_hll_csnmi
