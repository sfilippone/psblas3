
subroutine psb_z_diag_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_z_diagprec, psb_protect_name =>  psb_z_diag_apply_vect
  implicit none 
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_z_diag_prec_type), intent(inout)  :: prec
  type(psb_z_vect_type),intent(inout)   :: x
  complex(psb_dpk_),intent(in)         :: alpha, beta
  type(psb_z_vect_type),intent(inout)   :: y
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans
  complex(psb_dpk_),intent(inout), optional, target :: work(:)
  integer(psb_ipk_) :: err_act, nrow, ierr(5)
  character(len=20)  :: name='z_diag_prec_apply'
  complex(psb_dpk_), pointer :: ww(:)
  class(psb_z_base_vect_type), allocatable :: dw

  call psb_erractionsave(err_act)

  !
  ! This is the base version and we should throw an error. 
  ! Or should it be the DIAG preonditioner???
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
  if (.not.allocated(prec%d)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if
  if (size(prec%d) < nrow) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if

  if (size(work) >= x%get_nrows()) then 
    ww => work
  else
    allocate(ww(x%get_nrows()),stat=info)
    if (info /= psb_success_) then 
      ierr(1) = x%get_nrows()
      call psb_errpush(psb_err_alloc_request_,name,&
           & i_err=ierr,a_err='complex(psb_dpk_)')
      goto 9999      
    end if
  end if


  call y%mlt(alpha,prec%dv,x,beta,info,conjgx=trans)

  if (size(work) < x%get_nrows()) then 
    deallocate(ww,stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_, &
           & name,a_err='Deallocate')
      goto 9999      
    end if
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

end subroutine psb_z_diag_apply_vect


subroutine psb_z_diag_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_z_diagprec, psb_protect_name =>  psb_z_diag_apply
  implicit none 

  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_z_diag_prec_type), intent(in)  :: prec
  complex(psb_dpk_),intent(inout)      :: x(:)
  complex(psb_dpk_),intent(in)         :: alpha, beta
  complex(psb_dpk_),intent(inout)      :: y(:)
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans
  complex(psb_dpk_),intent(inout), optional, target :: work(:)
  integer(psb_ipk_) :: err_act, nrow, ierr(5)
  character :: trans_
  character(len=20)  :: name='z_diag_prec_apply'
  complex(psb_dpk_), pointer :: ww(:)

  call psb_erractionsave(err_act)

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
  if (.not.allocated(prec%d)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if
  if (size(prec%d) < nrow) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if
  if (present(trans)) then 
    trans_ = psb_toupper(trans)
  else
    trans_='N'
  end if

  select case(trans_)
  case('N')
  case('T','C')
  case default
    info=psb_err_iarg_invalid_i_
    ierr(1) = 6
    call psb_errpush(info,name,i_err=ierr,a_err=trans_)
    goto 9999
  end select

  if (size(work) >= size(x)) then 
    ww => work
  else
    allocate(ww(size(x)),stat=info)
    if (info /= psb_success_) then 
      ierr(1) = size(x)
      call psb_errpush(psb_err_alloc_request_,name,&
           & i_err=ierr,a_err='complex(psb_dpk_)')
      goto 9999      
    end if
  end if


  if (trans_ == 'C') then 
    ww(1:nrow) = x(1:nrow)*conjg(prec%d(1:nrow))
  else
    ww(1:nrow) = x(1:nrow)*prec%d(1:nrow)
  endif
  call psb_geaxpby(alpha,ww,beta,y,desc_data,info)

  if (size(work) < size(x)) then 
    deallocate(ww,stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Deallocate')
      goto 9999      
    end if
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

end subroutine psb_z_diag_apply


subroutine psb_z_diag_precbld(a,desc_a,prec,info,upd,amold,afmt,vmold)
  use psb_base_mod
  use psb_z_diagprec, psb_protect_name =>  psb_z_diag_precbld

  Implicit None

  type(psb_zspmat_type), intent(in), target :: a
  type(psb_desc_type), intent(in), target   :: desc_a
  class(psb_z_diag_prec_type),intent(inout) :: prec
  integer(psb_ipk_), intent(out)                      :: info
  character, intent(in), optional           :: upd
  character(len=*), intent(in), optional    :: afmt
  class(psb_z_base_sparse_mat), intent(in), optional :: amold
  class(psb_z_base_vect_type), intent(in), optional  :: vmold
  integer(psb_ipk_) :: err_act, nrow,i
  character(len=20)  :: name='z_diag_precbld'

  call psb_erractionsave(err_act)

  info = psb_success_
  nrow = desc_a%get_local_cols()
  if (allocated(prec%d)) then 
    if (size(prec%d) < nrow) then 
      deallocate(prec%d,stat=info)
    end if
  end if
  if ((info == psb_success_).and.(.not.allocated(prec%d))) then 
    allocate(prec%d(nrow), stat=info)
  end if
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  call a%get_diag(prec%d,info) 
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='get_diag')
    goto 9999
  end if

  do i=1,nrow
    if (prec%d(i) == dzero) then
      prec%d(i) = done
    else
      prec%d(i) = done/prec%d(i)
    endif
  end do
  allocate(prec%dv,stat=info) 
  if (info == 0) then 
    if (present(vmold)) then 
      allocate(prec%dv%v,mold=vmold,stat=info) 
    else
      allocate(psb_z_base_vect_type :: prec%dv%v,stat=info) 
    end if
  end if
  if (info == 0) then 
    call prec%dv%bld(prec%d)
  else 
    write(0,*) 'Error on precbld ',info
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
end subroutine psb_z_diag_precbld

