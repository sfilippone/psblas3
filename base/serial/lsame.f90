function lsame(a,b)
  use psb_string_mod
  logical :: lsame
  character(len=1) :: a, b 
  
  lsame = (tolower(a) == tolower(b))
end function lsame
