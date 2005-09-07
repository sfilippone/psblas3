module psb_string_mod
  interface tolower
    function  tolowerc(string)
      character(len=*), intent(in)  :: string
      character(len=len(string))    :: tolowerc
    end function tolowerc
  end interface
  interface toupper
    function  toupperc(string)
      character(len=*), intent(in)  :: string
      character(len=len(string))    :: toupperc
    end function toupperc
  end interface
  interface touppers
    subroutine   sub_toupperc(string,strout)
      character(len=*), intent(in)  :: string
      character(len=*), intent(out)  :: strout
    end subroutine sub_toupperc
  end interface
end module psb_string_mod
