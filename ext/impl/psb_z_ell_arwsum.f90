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
  

subroutine psb_z_ell_arwsum(d,a) 
  
  use psb_base_mod
  use psb_z_ell_mat_mod, psb_protect_name => psb_z_ell_arwsum
  implicit none 
  class(psb_z_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)

  integer(psb_ipk_)  :: i,j,k,m,n, nnz, ir, jc, nc
  logical            :: tra, is_unit
  Integer(Psb_ipk_)  :: err_act, info, int_err(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (a%is_dev()) call a%sync()

  m = a%get_nrows()
  if (size(d) < m) then 
    info=psb_err_input_asize_small_i_
    int_err(1) = 1
    int_err(2) = size(d)
    int_err(3) = m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  is_unit = a%is_unit()

  do i = 1, a%get_nrows()
    if (is_unit) then 
      d(i) = done
    else
      d(i) = dzero
    end if
    do j=1,a%irn(i)
      d(i) = d(i) + abs(a%val(i,j))
    end do
  end do

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(err_act)
  return

end subroutine psb_z_ell_arwsum
