subroutine psb_d_null_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_d_nullprec, psb_protect_name => psb_d_null_apply_vect
  implicit none 
  type(psb_desc_type),intent(in)       :: desc_data
  class(psb_d_null_prec_type), intent(inout)  :: prec
  type(psb_d_vect_type),intent(inout)  :: x
  real(psb_dpk_),intent(in)         :: alpha, beta
  type(psb_d_vect_type),intent(inout)  :: y
  integer(psb_ipk_), intent(out)                 :: info
  character(len=1), optional           :: trans
  real(psb_dpk_),intent(inout), optional, target :: work(:)
  integer(psb_ipk_) :: err_act, nrow, ierr(5)
  character(len=20)  :: name='c_null_prec_apply'

  call psb_erractionsave(err_act)

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

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_null_apply_vect

subroutine psb_d_null_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_d_nullprec, psb_protect_name => psb_d_null_apply
  implicit none 
  type(psb_desc_type),intent(in)       :: desc_data
  class(psb_d_null_prec_type), intent(in)  :: prec
  real(psb_dpk_),intent(inout)      :: x(:)
  real(psb_dpk_),intent(in)         :: alpha, beta
  real(psb_dpk_),intent(inout)      :: y(:)
  integer(psb_ipk_), intent(out)                 :: info
  character(len=1), optional           :: trans
  real(psb_dpk_),intent(inout), optional, target :: work(:)
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

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_null_apply
