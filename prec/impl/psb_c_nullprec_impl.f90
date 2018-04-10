!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
!    
subroutine psb_c_null_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_c_nullprec, psb_protect_name => psb_c_null_apply_vect
  implicit none 
  type(psb_desc_type),intent(in)       :: desc_data
  class(psb_c_null_prec_type), intent(inout)  :: prec
  type(psb_c_vect_type),intent(inout)  :: x
  complex(psb_spk_),intent(in)         :: alpha, beta
  type(psb_c_vect_type),intent(inout)  :: y
  integer(psb_ipk_), intent(out)                 :: info
  character(len=1), optional           :: trans
  complex(psb_spk_),intent(inout), optional, target :: work(:)
  integer(psb_ipk_) :: err_act, nrow, ierr(5)
  character(len=20)  :: name='c_null_prec_apply'

  call psb_erractionsave(err_act)

  !
  ! This is the base version and we should throw an error. 
  ! Or should it be the NULL preonditioner???
  !
  info = psb_success_

  nrow = desc_data%get_local_rows()
  if (x%get_nrows() < nrow) then 
    info = 36; ierr(1) = 2; ierr(2) = nrow;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (y%get_nrows() < nrow) then 
    info = 36; ierr(1) = 3; ierr(2) = nrow;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  call psb_geaxpby(alpha,x,beta,y,desc_data,info)
  if (info /= psb_success_ ) then 
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="psb_geaxpby")
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_c_null_apply_vect

subroutine psb_c_null_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_c_nullprec, psb_protect_name => psb_c_null_apply
  implicit none 
  type(psb_desc_type),intent(in)       :: desc_data
  class(psb_c_null_prec_type), intent(inout)  :: prec
  complex(psb_spk_),intent(inout)      :: x(:)
  complex(psb_spk_),intent(in)         :: alpha, beta
  complex(psb_spk_),intent(inout)      :: y(:)
  integer(psb_ipk_), intent(out)                 :: info
  character(len=1), optional           :: trans
  complex(psb_spk_),intent(inout), optional, target :: work(:)
  integer(psb_ipk_) :: err_act, nrow, ierr(5)
  character(len=20)  :: name='c_null_prec_apply'

  call psb_erractionsave(err_act)

  !
  !
  info = psb_success_

  nrow = desc_data%get_local_rows()
  if (size(x) < nrow) then 
    info = 36; ierr(1) = 2; ierr(2) = nrow;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (size(y) < nrow) then 
    info = 36; ierr(1) = 3; ierr(2) = nrow;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  call psb_geaxpby(alpha,x,beta,y,desc_data,info)
  if (info /= psb_success_ ) then 
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="psb_geaxpby")
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_c_null_apply
