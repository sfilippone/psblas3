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
  

subroutine psb_c_hll_get_diag(a,d,info) 
  
  use psb_base_mod
  use psb_c_hll_mat_mod, psb_protect_name => psb_c_hll_get_diag
  implicit none 
  class(psb_c_hll_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(out)             :: d(:)
  integer(psb_ipk_), intent(out)          :: info

  Integer(Psb_ipk_)  :: err_act, mnm, i, j, k, ke, hksz, ld,ir, mxrwl
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_dev()) call a%sync()

  mnm = min(a%get_nrows(),a%get_ncols())
  ld = size(d)
  if (ld< mnm) then 
    info=psb_err_input_asize_invalid_i_
    call psb_errpush(info,name,i_err=(/2*ione,ld/))
    goto 9999
  end if

  if (a%is_triangle().and.a%is_unit()) then 
    d(1:mnm) = cone
  else

    hksz = a%get_hksz()
    j=1
    do i=1,mnm,hksz
      ir    = min(hksz,mnm-i+1) 
      mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
      k     = a%hkoffs(j) + 1
      ke    = a%hkoffs(j+1)
      call psb_c_hll_get_diag_inner(ir,a%irn(i:i+ir-1),&
           & a%ja(k:ke),hksz,a%val(k:ke),hksz,&
           & a%idiag(i:i+ir-1),d(i:i+ir-1),info) 
      if (info /= psb_success_) goto 9999
      j = j + 1 
    end do

  end if

  do i=mnm+1,size(d) 
    d(i) = czero
  end do
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains

  subroutine psb_c_hll_get_diag_inner(m,irn,ja,ldj,val,ldv,&
       & idiag,d,info) 
    integer(psb_ipk_), intent(in)    :: m,ldj,ldv,ja(ldj,*),irn(*), idiag(*)
    complex(psb_spk_), intent(in)      :: val(ldv,*)
    complex(psb_spk_), intent(inout)   :: d(*)
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_) :: i,j,k, m4, jc

    info = psb_success_
    
    do i=1,m
      if (idiag(i) /= 0) then 
        d(i) = val(i,idiag(i))
      else
        d(i) = czero
      end if
    end do

  end subroutine psb_c_hll_get_diag_inner

end subroutine psb_c_hll_get_diag
