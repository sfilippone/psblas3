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
subroutine psb_d_dia_colsum(d,a) 
  
  use psb_base_mod
  use psb_d_dia_mat_mod, psb_protect_name => psb_d_dia_colsum
  implicit none 
  class(psb_d_dia_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc, ir1,ir2, nr
  logical           :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
  character(len=20) :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (a%is_dev()) call a%sync()

  m = a%get_nrows()
  n = a%get_ncols()
  if (size(d) < n) then 
    info=psb_err_input_asize_small_i_
    int_err(1) = 1
    int_err(2) = size(d)
    int_err(3) = n
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (a%is_unit()) then 
    d = done
  else
    d = dzero
  end if

  nr = size(a%data,1)
  nc = size(a%data,2)
  do j=1,nc
    jc = a%offset(j) 
    if (jc > 0) then 
      ir1 = 1
      ir2 = nr - jc
    else
      ir1 = 1 - jc
      ir2 = nr
    end if
    do i=ir1, ir2
      d(i+jc) = d(i+jc) + a%data(i,j)
    enddo
  enddo

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(err_act)
  return

end subroutine psb_d_dia_colsum
