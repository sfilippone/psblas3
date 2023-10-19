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
  

subroutine psb_c_diag_csmv(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use diagdev_mod
  use psb_vectordev_mod
  use psb_c_diag_mat_mod, psb_protect_name => psb_c_diag_csmv
#else 
  use psb_c_diag_mat_mod
#endif
  implicit none 
  class(psb_c_diag_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  complex(psb_spk_)  :: acc
  type(c_ptr)    :: gpX, gpY
  logical        :: tra
  Integer        :: err_act
  character(len=20)  :: name='c_diag_csmv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')

  if (tra) then 
    m = a%get_ncols()
    n = a%get_nrows()
  else
    n = a%get_ncols()
    m = a%get_nrows()
  end if

  if (size(x,1)<n) then 
    info = 36
    call psb_errpush(info,name,i_err=(/3*ione,n,izero,izero,izero/))
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/5*ione,m,izero,izero,izero/))
    goto 9999
  end if

#ifdef HAVE_SPGPU
  if (tra) then 
    call a%psb_c_dia_sparse_mat%spmm(alpha,x,beta,y,info,trans) 
  else
    !
    ! Just to test, move X/Y to/from the GPU.
    ! 
    if (info == 0) &
         & info = FallocMultiVecDevice(gpX,1,size(x,1),spgpu_type_double)
    if (alpha /= dzero) then 
      if (info == 0) &
           & info = writeMultiVecDevice(gpX,x)
    end if
    if (info == 0) &
         & info = FallocMultiVecDevice(gpY,1,size(y,1),spgpu_type_double)
    if (beta /= dzero) then 
      if (info == 0) &
           & info = writeMultiVecDevice(gpY,y)
    end if
    if (info == 0)  &
         & info = spmvDiagDevice(a%deviceMat,alpha,gpX,beta,gpY)
    if (info == 0) &
         & info = readMultiVecDevice(gpY,y)
    if (info /= 0) goto 9999
    call freeMultiVecDevice(gpX)
    call freeMultiVecDevice(gpY)
  endif
#else
  call a%psb_c_dia_sparse_mat%spmm(alpha,x,beta,y,info,trans) 
#endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return


end subroutine psb_c_diag_csmv
