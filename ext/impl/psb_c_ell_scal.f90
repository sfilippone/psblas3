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
  

subroutine psb_c_ell_scal(d,a,info,side) 
  
  use psb_base_mod
  use psb_c_ell_mat_mod, psb_protect_name => psb_c_ell_scal
  implicit none 
  class(psb_c_ell_sparse_mat), intent(inout) :: a
  complex(psb_spk_), intent(in)     :: d(:)
  integer(psb_ipk_), intent(out)  :: info
  character, intent(in), optional :: side

  Integer(Psb_ipk_)  :: err_act,mnm, i, j, m, n, ierr(5)
  character(len=20)  :: name='scal'
  character          :: side_
  logical            :: left 
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_dev()) call a%sync()
  if (a%is_unit()) then 
    call a%make_nonunit()
  end if

  side_ = 'L'
  if (present(side)) then 
    side_ = psb_toupper(side)
  end if

  left = (side_ == 'L')

  if (left) then 
    m = a%get_nrows()
    if (size(d) < m) then 
      info=psb_err_input_asize_invalid_i_
      call psb_errpush(info,name,i_err=(/2*ione,size(d,kind=psb_ipk_)/))
      goto 9999
    end if
    
    do i=1, m 
      a%val(i,:) = a%val(i,:) * d(i)
    enddo
  else
    n = a%get_ncols()
    if (size(d) < n) then 
      info=psb_err_input_asize_invalid_i_
      ierr(1) = 2; ierr(2) = size(d); 
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end if

    do i=1, m 
      do j=1, a%irn(i) 
        a%val(i,j) = a%val(i,j) * d(a%ja(i,j))
      end do
    enddo
    
  end if

  call a%set_host()
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_c_ell_scal
