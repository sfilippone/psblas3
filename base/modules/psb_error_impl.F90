! checks wether an error has occurred on one of the porecesses in the execution pool
subroutine psb_errcomm(ictxt, err)
  use psb_error_mod, psb_protect_name => psb_errcomm
  use psb_penv_mod
  integer, intent(in)   :: ictxt
  integer, intent(inout):: err
  integer :: temp(2)
  ! Cannot use psb_amx or otherwise we have a recursion in module usage
#if !defined(SERIAL_MPI)
  call psb_amx(ictxt, err)
#endif    
end subroutine psb_errcomm
! handles the occurence of an error in a serial routine
subroutine psb_serror()
  use psb_error_mod!, psb_protect_name => psb_serror

  integer                 ::  err_c
  character(len=20)       ::  r_name
  character(len=40)       ::  a_e_d
  integer                 ::  i_e_d(5)

  if(error_status > 0) then
    if(verbosity_level > 1) then

      do while (psb_get_numerr() > izero)
        write(0,'(50("="))')
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        call psb_errmsg(err_c, r_name, i_e_d, a_e_d)
        !            write(0,'(50("="))')
      end do

    else

      call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      call psb_errmsg(err_c, r_name, i_e_d, a_e_d)
      do while (psb_get_numerr() > 0)
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      end do
    end if
  end if

end subroutine psb_serror


! handles the occurence of an error in a parallel routine
subroutine psb_perror(ictxt)
  use psb_error_mod!, psb_protect_name => psb_perror
  use psb_penv_mod

  integer, intent(in)     :: ictxt
  integer                 :: err_c
  character(len=20)       :: r_name
  character(len=40)       :: a_e_d
  integer                 :: i_e_d(5)
  integer                 :: iam, np

#if defined(SERIAL_MPI)
  me = -1
#else        
  call psb_info(ictxt,iam,np)
#endif


  if(error_status > 0) then
    if(verbosity_level > 1) then

      do while (psb_get_numerr() > izero)
        write(0,'(50("="))')
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        call psb_errmsg(err_c, r_name, i_e_d, a_e_d,me)
        !            write(0,'(50("="))')
      end do
#if defined(SERIAL_MPI)
      stop 
#else        
      call psb_abort(ictxt,-1)
#endif
    else

      call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      call psb_errmsg(err_c, r_name, i_e_d, a_e_d,me)
      do while (psb_get_numerr() > 0)
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      end do
#if defined(SERIAL_MPI)
      stop 
#else        
      call psb_abort(ictxt,-1)
#endif
    end if
  end if

  if(error_status > izero) then
#if defined(SERIAL_MPI)
    stop 
#else        
    call psb_abort(ictxt,err_c)
#endif
  end if


end subroutine psb_perror

