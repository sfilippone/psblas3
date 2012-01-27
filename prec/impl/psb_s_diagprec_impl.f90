
subroutine psb_s_diag_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_s_diagprec, psb_protect_name =>  psb_s_diag_apply_vect
  implicit none 
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_s_diag_prec_type), intent(inout)  :: prec
  type(psb_s_vect_type),intent(inout)   :: x
  real(psb_spk_),intent(in)         :: alpha, beta
  type(psb_s_vect_type),intent(inout)   :: y
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans
  real(psb_spk_),intent(inout), optional, target :: work(:)
  integer(psb_ipk_) :: err_act, nrow
  character(len=20)  :: name='s_diag_prec_apply'
  real(psb_spk_), pointer :: ww(:)
  class(psb_s_base_vect_type), allocatable :: dw

  call psb_erractionsave(err_act)

  !
  ! This is the base version and we should throw an error. 
  ! Or should it be the DIAG preonditioner???
  !
  info = psb_success_

  nrow = desc_data%get_local_rows()
  if (x%get_nrows() < nrow) then 
    info = 36
    call psb_errpush(info,name,i_err=(/2,nrow,0,0,0/))
    goto 9999
  end if
  if (y%get_nrows() < nrow) then 
    info = 36
    call psb_errpush(info,name,i_err=(/3,nrow,0,0,0/))
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
      call psb_errpush(psb_err_alloc_request_,name,&
           & i_err=(/x%get_nrows(),0,0,0,0/),a_err='real(psb_spk_)')
      goto 9999      
    end if
  end if


  call y%mlt(alpha,prec%dv,x,beta,info,conjgx=trans)

  if (size(work) < x%get_nrows()) then 
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

end subroutine psb_s_diag_apply_vect


subroutine psb_s_diag_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_s_diagprec, psb_protect_name =>  psb_s_diag_apply
  implicit none 

  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_s_diag_prec_type), intent(in)  :: prec
  real(psb_spk_),intent(inout)      :: x(:)
  real(psb_spk_),intent(in)         :: alpha, beta
  real(psb_spk_),intent(inout)      :: y(:)
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans
  real(psb_spk_),intent(inout), optional, target :: work(:)
  integer(psb_ipk_) :: err_act, nrow
  character :: trans_
  character(len=20)  :: name='s_diag_prec_apply'
  real(psb_spk_), pointer :: ww(:)

  call psb_erractionsave(err_act)

  info = psb_success_

  nrow = desc_data%get_local_rows()
  if (size(x) < nrow) then 
    info = 36
    call psb_errpush(info,name,i_err=(/2,nrow,0,0,0/))
    goto 9999
  end if
  if (size(y) < nrow) then 
    info = 36
    call psb_errpush(info,name,i_err=(/3,nrow,0,0,0/))
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
    call psb_errpush(info,name,&
         & i_err=(/6,0,0,0,0/),a_err=trans_)
    goto 9999
  end select

  if (size(work) >= size(x)) then 
    ww => work
  else
    allocate(ww(size(x)),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_request_,name,&
           & i_err=(/size(x),0,0,0,0/),a_err='real(psb_spk_)')
      goto 9999      
    end if
  end if


  if (trans_ == 'C') then 
    ww(1:nrow) = x(1:nrow)*(prec%d(1:nrow))
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

end subroutine psb_s_diag_apply


subroutine psb_s_diag_precbld(a,desc_a,prec,info,upd,amold,afmt,vmold)
  use psb_base_mod
  use psb_s_diagprec, psb_protect_name =>  psb_s_diag_precbld

  Implicit None

  type(psb_sspmat_type), intent(in), target :: a
  type(psb_desc_type), intent(in), target   :: desc_a
  class(psb_s_diag_prec_type),intent(inout) :: prec
  integer(psb_ipk_), intent(out)                      :: info
  character, intent(in), optional           :: upd
  character(len=*), intent(in), optional    :: afmt
  class(psb_s_base_sparse_mat), intent(in), optional :: amold
  class(psb_s_base_vect_type), intent(in), optional  :: vmold
  integer(psb_ipk_) :: err_act, nrow,i
  character(len=20)  :: name='s_diag_precbld'

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
      allocate(psb_s_base_vect_type :: prec%dv%v,stat=info) 
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
end subroutine psb_s_diag_precbld

