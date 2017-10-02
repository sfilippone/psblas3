module psb_base_string_cbind_mod
  use iso_c_binding

contains

  subroutine stringc2f(cstring,fstring) 
    character(c_char)        :: cstring(*)
    character(len=*)         :: fstring
    integer :: i
    
    i = 1
    do 
      if (cstring(i) == c_null_char) exit
      if (i > len(fstring)) exit
      fstring(i:i) = cstring(i)
      i = i + 1 
    end do
    do 
      if (i > len(fstring)) exit
      fstring(i:i) = " "
      i = i + 1 
    end do
    return
  end subroutine stringc2f

  subroutine stringf2c(fstring,cstring)
    character(c_char)        :: cstring(*)
    character(len=*)         :: fstring
    integer :: i
    
    do i=1, len(fstring)
      cstring(i) = fstring(i:i)
    end do
    cstring(len(fstring)+1) = c_null_char
    return
  end subroutine stringf2c

end module psb_base_string_cbind_mod
