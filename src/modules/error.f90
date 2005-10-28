!
!  Wrapper subroutines to provide error tools to F77 and C code
!

subroutine FCpsb_errcomm(icontxt, err)
  use psb_error_mod
  integer, intent(in)   :: icontxt
  integer, intent(inout):: err

  call psb_errcomm(icontxt, err)

end subroutine FCpsb_errcomm

subroutine FCpsb_errpush(err_c, r_name, i_err)
  use psb_error_mod
  implicit none
  
  integer, intent(in)              ::  err_c
  character(len=20), intent(in)    ::  r_name
  integer                          ::  i_err(5)

  call psb_errpush(err_c, r_name, i_err)
  
end subroutine FCpsb_errpush



subroutine FCpsb_serror()
  use psb_error_mod
  implicit none

  call psb_error()

end subroutine FCpsb_serror





subroutine FCpsb_perror(icontxt)
  use psb_error_mod
  implicit none

  integer, intent(in)   :: icontxt

  call psb_error(icontxt)

end subroutine FCpsb_perror





function FCpsb_get_errstatus()
  use psb_error_mod
  implicit none

  integer :: FCpsb_get_errstatus

  FCpsb_get_errstatus = psb_get_errstatus()

end function FCpsb_get_errstatus





subroutine FCpsb_get_errverbosity(v)
  use psb_error_mod
  implicit none

  integer, intent(out)   :: v

  call psb_get_errverbosity(v)

end subroutine FCpsb_get_errverbosity




subroutine FCpsb_set_errverbosity(v)
  use psb_error_mod
  implicit none

  integer, intent(inout)   :: v

  call psb_set_errverbosity(v)

end subroutine FCpsb_set_errverbosity





subroutine FCpsb_erractionsave(err_act)
  use psb_error_mod
  implicit none

  integer, intent(out) :: err_act

  call psb_erractionsave(err_act)

end subroutine FCpsb_erractionsave


subroutine FCpsb_get_erraction(err_act)
  use psb_error_mod
  implicit none
  integer, intent(out) :: err_act 

  call psb_get_erraction(err_act)
end subroutine FCpsb_get_erraction



subroutine FCpsb_erractionrestore(err_act)
  use psb_error_mod
  implicit none

  integer, intent(in) :: err_act

  call psb_erractionrestore(err_act)

end subroutine FCpsb_erractionrestore






