! checks wether an error has occurred on one of the porecesses in the execution pool
subroutine psb_errcomm(ictxt, err)
  use psb_error_mod, psb_protect_name => psb_errcomm
  use psb_penv_mod
  integer(psb_mpik_), intent(in)   :: ictxt
  integer(psb_ipk_), intent(inout):: err
  integer(psb_ipk_) :: temp(2)
  
  call psb_amx(ictxt, err)

end subroutine psb_errcomm

! handles the occurence of an error in a serial routine
subroutine psb_serror()
  use psb_const_mod
  use psb_error_mod
  implicit none 
  integer(psb_ipk_) ::  err_c
  character(len=20)       ::  r_name
  character(len=40)       ::  a_e_d
  integer(psb_ipk_) ::  i_e_d(5)

  if (psb_errstatus_fatal()) then
    if(psb_get_errverbosity() > 1) then

      do while (psb_get_numerr() > izero)
        write(psb_err_unit,'(50("="))')
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        call psb_errmsg(err_c, r_name, i_e_d, a_e_d)
        !            write(psb_err_unit,'(50("="))')
      end do

    else

      call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      call psb_errmsg(err_c, r_name, i_e_d, a_e_d)
      do while (psb_get_numerr() > 0)
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      end do
    end if
  end if
#if defined(HAVE_FLUSH_STMT)
  flush(0) 
#endif


end subroutine psb_serror


! handles the occurence of an error in a parallel routine
subroutine psb_perror(ictxt)
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none 
  integer(psb_mpik_), intent(in)     :: ictxt
  integer(psb_ipk_) :: err_c
  character(len=20)       :: r_name
  character(len=40)       :: a_e_d
  integer(psb_ipk_) :: i_e_d(5)
  integer(psb_mpik_) :: iam, np

  call psb_info(ictxt,iam,np)
  
  if (psb_errstatus_fatal()) then
    if (psb_get_errverbosity() > 1) then

      do while (psb_get_numerr() > izero)
        write(psb_err_unit,'(50("="))')
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        call psb_errmsg(err_c, r_name, i_e_d, a_e_d,iam)
        !            write(psb_err_unit,'(50("="))')
      end do
#if defined(HAVE_FLUSH_STMT)
      flush(0) 
#endif

      call psb_abort(ictxt,-1)

    else

      call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      call psb_errmsg(err_c, r_name, i_e_d, a_e_d,iam)
      do while (psb_get_numerr() > 0)
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      end do
#if defined(HAVE_FLUSH_STMT)
      flush(0) 
#endif

      call psb_abort(ictxt,-1)

    end if
  end if

end subroutine psb_perror

