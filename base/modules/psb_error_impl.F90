! checks wether an error has occurred on one of the porecesses in the execution pool
subroutine psb_errcomm_i(ctxt, err)
  use psb_error_mod, psb_protect_name => psb_errcomm
  use psb_penv_mod
  type(psb_ctxt_type), intent(in)  :: ctxt
  integer(psb_ipk_), intent(inout) :: err

  if (psb_get_global_checks()) call psb_amx(ctxt, err)

end subroutine psb_errcomm_i

#if defined(IPK8)

subroutine psb_errcomm_m(ctxt, err)
  use psb_error_mod, psb_protect_name => psb_errcomm
  use psb_penv_mod
  type(psb_ctxt_type), intent(in)  :: ctxt
  integer(psb_ipk_), intent(inout) :: err
  
  if (psb_get_global_checks()) call psb_amx(ctxt, err)

end subroutine psb_errcomm_m
#endif

subroutine psb_ser_error_handler(err_act)
  use psb_error_mod, psb_protect_name => psb_ser_error_handler
  use psb_penv_mod
  implicit none 
  integer(psb_ipk_), intent(inout) ::  err_act
  
  call psb_erractionrestore(err_act)
  
  if (err_act /= psb_act_ret_)     &
       &  call psb_error()
  if (err_act == psb_act_abort_) stop

  return 
end subroutine psb_ser_error_handler

subroutine psb_par_error_handler(ctxt,err_act)
  use psb_error_mod, psb_protect_name => psb_par_error_handler
  use psb_penv_mod
  implicit none 
  type(psb_ctxt_type), intent(in) ::  ctxt
  integer(psb_ipk_), intent(in)   ::  err_act

  call psb_erractionrestore(err_act)

  if (err_act == psb_act_print_)    then 
    call psb_error(ctxt, abrt=.false.)
  else if (err_act == psb_act_ret_)     then
    ! do nothing 
  else 
    call psb_error(ctxt, abrt=.true.)
  end if
  
  return 

end subroutine psb_par_error_handler

subroutine psb_par_error_print_stack(ctxt)
  use psb_error_mod, psb_protect_name => psb_par_error_print_stack
  use psb_penv_mod
  type(psb_ctxt_type), intent(in) ::  ctxt

  call psb_error(ctxt, abrt=.false.)
  
end subroutine psb_par_error_print_stack

subroutine psb_ser_error_print_stack()
  use psb_error_mod, psb_protect_name => psb_ser_error_print_stack

  call psb_error()
end subroutine psb_ser_error_print_stack




! handles the occurence of an error in a serial routine
subroutine psb_serror()
  use psb_const_mod
  use psb_error_mod
  implicit none 
  integer(psb_ipk_) ::  err_c
  character(len=20) ::  r_name
  character(len=40) ::  a_e_d
  integer(psb_epk_) ::  e_e_d(5)

  if (psb_errstatus_fatal()) then
    if(psb_get_errverbosity() > 1) then

      do while (psb_get_numerr() > izero)
        write(psb_err_unit,'(50("="))')
        call psb_errpop(err_c, r_name, e_e_d, a_e_d)
        call psb_errmsg(psb_err_unit,err_c, r_name, e_e_d, a_e_d)
        !            write(psb_err_unit,'(50("="))')
      end do

    else

      call psb_errpop(err_c, r_name, e_e_d, a_e_d)
      call psb_errmsg(psb_err_unit,err_c, r_name, e_e_d, a_e_d)
      do while (psb_get_numerr() > 0)
        call psb_errpop(err_c, r_name, e_e_d, a_e_d)
      end do
    end if
  end if
#if defined(HAVE_FLUSH_STMT)
  flush(psb_err_unit) 
#endif


end subroutine psb_serror


! handles the occurence of an error in a parallel routine
subroutine psb_perror(ctxt,abrt)
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none 
  type(psb_ctxt_type), intent(in) :: ctxt
  logical, intent(in), optional   :: abrt

  integer(psb_ipk_) :: err_c
  character(len=20) :: r_name
  character(len=40) :: a_e_d
  integer(psb_epk_) :: e_e_d(5)
  integer(psb_mpk_) :: iam, np
  logical :: abrt_

  abrt_=.true.
  if (present(abrt)) abrt_=abrt
  call psb_info(ctxt,iam,np)
  
  if (psb_errstatus_fatal()) then
    if (psb_get_errverbosity() > 1) then

      do while (psb_get_numerr() > izero)
        write(psb_err_unit,'(50("="))')
        call psb_errpop(err_c, r_name, e_e_d, a_e_d)
        call psb_errmsg(psb_err_unit,err_c, r_name, e_e_d, a_e_d,iam)
        !            write(psb_err_unit,'(50("="))')
      end do
#if defined(HAVE_FLUSH_STMT)
      flush(psb_err_unit) 
#endif
      
      if (abrt_) call psb_abort(ctxt,-1)
      
    else

      call psb_errpop(err_c, r_name, e_e_d, a_e_d)
      call psb_errmsg(psb_err_unit,err_c, r_name, e_e_d, a_e_d,iam)
      do while (psb_get_numerr() > 0)
        call psb_errpop(err_c, r_name, e_e_d, a_e_d)
      end do
#if defined(HAVE_FLUSH_STMT)
      flush(psb_err_unit) 
#endif

      if (abrt_) call psb_abort(ctxt,-1)

    end if
  end if

end subroutine psb_perror

