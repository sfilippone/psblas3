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
  use psb_const_mod
  use psb_error_mod
  implicit none 
  integer                 ::  err_c
  character(len=20)       ::  r_name
  character(len=40)       ::  a_e_d
  integer                 ::  i_e_d(5)

  if(psb_get_errstatus() > 0) then
    if(psb_get_errverbosity() > 1) then

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
  integer, intent(in)     :: ictxt
  integer                 :: err_c
  character(len=20)       :: r_name
  character(len=40)       :: a_e_d
  integer                 :: i_e_d(5)
  integer                 :: iam, np

#if defined(SERIAL_MPI)
  iam = -1
#else        
  call psb_info(ictxt,iam,np)
#endif
  
  
  if(psb_get_errstatus() > 0) then
    if(psb_get_errverbosity() > 1) then

      do while (psb_get_numerr() > izero)
        write(0,'(50("="))')
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        call psb_errmsg(err_c, r_name, i_e_d, a_e_d,iam)
        !            write(0,'(50("="))')
      end do
#if defined(HAVE_FLUSH_STMT)
      flush(0) 
#endif
#if defined(SERIAL_MPI)
      stop 
#else        
      call psb_abort(ictxt,-1)
#endif
    else

      call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      call psb_errmsg(err_c, r_name, i_e_d, a_e_d,iam)
      do while (psb_get_numerr() > 0)
        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
      end do
#if defined(HAVE_FLUSH_STMT)
      flush(0) 
#endif
#if defined(SERIAL_MPI)
      stop 
#else        
      call psb_abort(ictxt,-1)
#endif
    end if
  end if

  if(psb_get_errstatus() > izero) then
#if defined(SERIAL_MPI)
    stop 
#else        
    call psb_abort(ictxt,err_c)
#endif
  end if


end subroutine psb_perror

