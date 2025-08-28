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
  

subroutine psb_z_hll_scal(d,a,info,side) 

  use psb_base_mod
  use psb_z_hll_mat_mod, psb_protect_name => psb_z_hll_scal
  implicit none 
  class(psb_z_hll_sparse_mat), intent(inout) :: a
  complex(psb_dpk_), intent(in)                :: d(:)
  integer(psb_ipk_), intent(out)            :: info
  character, intent(in), optional           :: side

  Integer(Psb_ipk_) :: err_act,mnm, i, j, m, n, ierr(5), ld, k, mxrwl, hksz, ir
  character(len=20) :: name='scal'
  character         :: side_
  logical           :: left 
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_dev()) call a%sync()

  info = psb_err_missing_override_method_
  call psb_errpush(info,name,i_err=ierr)
  goto 9999

  side_ = 'L'
  if (present(side)) then 
    side_ = psb_toupper(side)
  end if

  left = (side_ == 'L')

  ld = size(d)
  if (left) then 
    m = a%get_nrows()
    if (ld < m) then 
      ierr(1) = 2; ierr(2) = ld; 
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end if
  else
    n = a%get_ncols()
    if (ld < n) then 
      info=psb_err_input_asize_invalid_i_
      ierr(1) = 2; ierr(2) = ld; 
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end if
  end if

  hksz = a%get_hksz()
  j    = 1
  do i=1,m,hksz
    ir    = min(hksz,m-i+1) 
    mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
    k     = a%hkoffs(j) + 1
    call psb_z_hll_scal_inner(i,ir,mxrwl,a%irn(i),&
         & a%ja(k),hksz,a%val(k),hksz,&
         & left,d,info) 
    if (info /= psb_success_) goto 9999
    j = j + 1 
  end do


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains

  subroutine psb_z_hll_scal_inner(ir,m,n,irn,ja,ldj,val,ldv,left,d,info) 
    integer(psb_ipk_), intent(in)    :: ir,m,n,ldj,ldv,ja(ldj,*),irn(*)
    complex(psb_dpk_), intent(in)      :: d(*)
    complex(psb_dpk_), intent(inout)   :: val(ldv,*)
    logical, intent(in)              :: left
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_) :: i,j,k, m4, jc

    info = psb_success_

    if (left) then 
      do i=1,m
        do j=1, irn(i)
          val(i,j) = val(i,j)*d(ir+i-1)
        end do
      end do
    else
      do i=1,m
        do j=1, irn(i)
          jc       = ja(i,j)
          val(i,j) = val(i,j)*d(jc)
        end do
      end do

    end if

  end subroutine psb_z_hll_scal_inner


end subroutine psb_z_hll_scal
